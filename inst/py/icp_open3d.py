import open3d as o3d
import numpy as np
import laspy
import os
import pyproj


class Open3DICP:
    """
    ICP alignment using Open3D with LAS/LAZ IO via laspy, preserving key metadata.

    Improvements vs. your original:
      - multi-scale ICP (coarse -> fine)
      - threshold scales with voxel size
      - point-to-plane attempts first (if requested) and falls back to point-to-point on singular/failed solve
      - normals radius tied to voxel size (more stable)
      - optional centering to improve numeric stability
    """

    def __init__(
        self,
        source_path,
        target_path,
        voxel_size=0.05,
        icp_method="point-to-plane",
        max_iteration=50,
        threshold_factor=15.0,
        multiscale_factors=(4.0, 2.0, 1.0),
        multiscale_iters=(50, 30, 20),
        normal_radius_factor=6.0,
        normal_max_nn=50,
        min_correspondences=50,
        center=True,
    ):
        """
        :param source_path: Path to source point cloud (.las or .laz)
        :param target_path: Path to target point cloud (.las or .laz)
        :param voxel_size: Base voxel downsample size (meters)
        :param icp_method: "point-to-plane" or "point-to-point"
        :param max_iteration: (legacy single-scale) max iters; used if multiscale is disabled
        :param threshold_factor: threshold = voxel_size * threshold_factor
        :param multiscale_factors: tuple of voxel multipliers for coarse->fine
        :param multiscale_iters: tuple of iters corresponding to multiscale_factors
        :param normal_radius_factor: radius = voxel * normal_radius_factor
        :param normal_max_nn: max neighbors for normals
        :param min_correspondences: minimum correspondences required to proceed at a scale
        :param center: if True, center clouds before ICP for numeric stability
        """
        self.source_path = source_path
        self.target_path = target_path
        self.aligned_path = self._generate_aligned_path(source_path)

        self.voxel_size = float(voxel_size)
        self.icp_method = str(icp_method)
        self.max_iteration = int(max_iteration)
        self.threshold_factor = float(threshold_factor)

        self.multiscale_factors = tuple(float(x) for x in multiscale_factors)
        self.multiscale_iters = tuple(int(x) for x in multiscale_iters)
        if len(self.multiscale_factors) != len(self.multiscale_iters):
            raise ValueError("multiscale_factors and multiscale_iters must have the same length.")

        self.normal_radius_factor = float(normal_radius_factor)
        self.normal_max_nn = int(normal_max_nn)
        self.min_correspondences = int(min_correspondences)
        self.center = bool(center)

        self.transformation = None
        self.rmse = None
        self.original_crs = None
        self.original_header = None

    def align(self):
        """
        Performs ICP alignment, saves aligned LAS/LAZ, returns (aligned_path, message).
        """
        source_pcd, source_metadata = self._load_point_cloud(self.source_path)
        target_pcd, _ = self._load_point_cloud(self.target_path)

        if source_pcd is None or target_pcd is None:
            return None, "Error: Failed to load point clouds."

        # Optional centering for numeric stability, compose back after ICP
        T_center = np.eye(4, dtype=np.float64)
        if self.center:
            src_pts = np.asarray(source_pcd.points)
            tgt_pts = np.asarray(target_pcd.points)
            if src_pts.size == 0 or tgt_pts.size == 0:
                return None, "Error: One of the point clouds is empty."
            src_c = src_pts.mean(axis=0)
            tgt_c = tgt_pts.mean(axis=0)
            source_pcd.translate(-src_c)
            target_pcd.translate(-tgt_c)
            # To map source->target in original coords after ICP on centered:
            # x' = (x - src_c) -> ICP -> + tgt_c
            # Equivalent transform: T = [R, t + tgt_c - R*src_c]
            # We'll compute final at the end.
            T_center[:3, 3] = (tgt_c - src_c)

        # Run multi-scale ICP to get robust T
        try:
            T, rmse = self._multiscale_icp(source_pcd, target_pcd)
        except Exception as e:
            return None, f"Error during ICP alignment: {str(e)}"

        # Compose centering correction if used
        if self.center:
            # T maps centered source -> centered target
            # Final transform in original coords: x_tgt = R * x_src + (t + tgt_c - R*src_c)
            # Our T_center encodes (tgt_c - src_c) in translation; but must account for rotation on src_c.
            # We'll recompute properly:
            R = T[:3, :3]
            t = T[:3, 3]
            # Recover centroids used above
            # (We can re-read means from centered state by reversing translate, but easier:
            # just re-compute from original files would be expensive; instead, store them.)
            # So we store them by capturing just before translating:
            # -> simplest: re-load just centroids from metadata points we already had at load time.
            # We'll store them in self during centering:
            # (Handled below in _center_info)
            src_c = self._center_info["src_c"]
            tgt_c = self._center_info["tgt_c"]
            t_final = t + tgt_c - R @ src_c

            T_final = np.eye(4, dtype=np.float64)
            T_final[:3, :3] = R
            T_final[:3, 3] = t_final
        else:
            T_final = T

        self.transformation = T_final
        self.rmse = float(rmse)

        # Apply transformation to full-resolution source
        source_pcd_full, source_metadata_full = self._load_point_cloud(self.source_path)
        if source_pcd_full is None:
            return None, "Error: Failed to reload source point cloud for full-res transform."
        source_pcd_full.transform(self.transformation)

        # Merge metadata back
        aligned_data = self._merge_metadata(source_pcd_full, source_metadata_full)

        # Save aligned LAS/LAZ
        self._save_as_laz(aligned_data)

        output_str = (
            "ICP Alignment Completed.\n"
            f"Method requested: {self.icp_method}\n"
            f"RMSE: {self.rmse:.6f}\n"
            "Transformation Matrix:\n"
            f"{np.array_str(self.transformation, precision=6, suppress_small=True)}"
        )

        return self.aligned_path, output_str

    # -------------------- ICP internals --------------------

    def _multiscale_icp(self, source_pcd, target_pcd):
        """
        Coarse->fine ICP.
        Attempts point-to-plane first if requested; falls back to point-to-point if singular/fails.
        Returns (T, rmse).
        """
        T = np.eye(4, dtype=np.float64)
        last_rmse = None

        # Keep a record of centroids if centering is enabled
        if self.center:
            self._center_info = {
                "src_c": np.asarray(source_pcd.points).mean(axis=0),
                "tgt_c": np.asarray(target_pcd.points).mean(axis=0),
            }

        for factor, iters in zip(self.multiscale_factors, self.multiscale_iters):
            vs = self.voxel_size * factor
            if vs <= 0:
                raise ValueError("voxel_size must be > 0.")

            src = source_pcd.voxel_down_sample(vs)
            tgt = target_pcd.voxel_down_sample(vs)

            threshold = vs * self.threshold_factor

            # Sanity check: enough points?
            if len(src.points) < 100 or len(tgt.points) < 100:
                # Too few points at this scale; skip to next (finer) scale
                continue

            # Choose estimation, attempt p2l then fallback p2p
            use_p2l = (self.icp_method.lower() == "point-to-plane")

            def estimate_normals(pcd, voxel):
                rad = voxel * self.normal_radius_factor
                pcd.estimate_normals(
                    o3d.geometry.KDTreeSearchParamHybrid(radius=rad, max_nn=self.normal_max_nn)
                )
                pcd.normalize_normals()

            if use_p2l:
                estimate_normals(src, vs)
                estimate_normals(tgt, vs)

            # Ensure we have correspondences before solving
            eval0 = o3d.pipelines.registration.evaluate_registration(src, tgt, threshold, T)
            if eval0.correspondence_set is None or len(eval0.correspondence_set) < self.min_correspondences:
                # If too few correspondences, try increasing threshold slightly at this scale
                threshold2 = threshold * 2.0
                eval1 = o3d.pipelines.registration.evaluate_registration(src, tgt, threshold2, T)
                if eval1.correspondence_set is None or len(eval1.correspondence_set) < self.min_correspondences:
                    # Still too few; skip this scale
                    continue
                threshold = threshold2

            # Run ICP
            crit = o3d.pipelines.registration.ICPConvergenceCriteria(max_iteration=int(iters))

            if use_p2l:
                try:
                    result = o3d.pipelines.registration.registration_icp(
                        src,
                        tgt,
                        threshold,
                        T,
                        o3d.pipelines.registration.TransformationEstimationPointToPlane(),
                        crit,
                    )
                except RuntimeError:
                    # Fallback to point-to-point at this scale
                    result = o3d.pipelines.registration.registration_icp(
                        src,
                        tgt,
                        threshold,
                        T,
                        o3d.pipelines.registration.TransformationEstimationPointToPoint(),
                        crit,
                    )
            else:
                result = o3d.pipelines.registration.registration_icp(
                    src,
                    tgt,
                    threshold,
                    T,
                    o3d.pipelines.registration.TransformationEstimationPointToPoint(),
                    crit,
                )

            T = result.transformation
            last_rmse = float(result.inlier_rmse)

        if last_rmse is None:
            raise RuntimeError(
                "ICP failed: no scale produced enough correspondences. "
                "Try increasing voxel_size and/or threshold_factor."
            )

        return T, last_rmse

    # -------------------- IO + metadata --------------------

    def _load_point_cloud(self, file_path):
        """
        Loads LAS/LAZ into Open3D, capturing metadata arrays.
        Returns (pcd, metadata_dict).
        """
        if file_path.endswith(".las") or file_path.endswith(".laz"):
            with laspy.open(file_path) as las_file:
                las = las_file.read()

                points = np.vstack((las.x, las.y, las.z)).T.astype(np.float64)
                pcd = o3d.geometry.PointCloud()
                pcd.points = o3d.utility.Vector3dVector(points)

                # Store original LAS header and CRS
                self.original_header = las.header
                self.original_crs = self._parse_crs_safe(las.header)

                metadata = {
                    "classification": las.classification.copy(),
                    "intensity": las.intensity.copy(),
                    "return_number": las.return_number.copy(),
                    "num_returns": las.num_returns.copy(),
                }
                return pcd, metadata

        if file_path.endswith(".ply"):
            pcd = o3d.io.read_point_cloud(file_path)
            return pcd, None

        return None, None

    def _parse_crs_safe(self, header):
        try:
            crs = header.parse_crs()
            if crs is None:
                return pyproj.CRS.from_epsg(26917)
            if isinstance(crs, pyproj.CRS):
                if crs.is_compound:
                    return crs.sub_crs_list[0]
                return crs
            return pyproj.CRS.from_user_input(crs)
        except Exception:
            return pyproj.CRS.from_epsg(26917)

    def _generate_aligned_path(self, source_path):
        base, _ = os.path.splitext(source_path)
        return f"{base}_aligned.laz"

    def _merge_metadata(self, transformed_pcd, metadata):
        """
        Creates structured array for LAS output with transformed XYZ + original metadata.
        Assumes point count unchanged (true for rigid transform).
        """
        pts = np.asarray(transformed_pcd.points)
        num_points = pts.shape[0]

        if metadata is None:
            # No metadata available, still output XYZ
            new_data = np.zeros(num_points, dtype=[("x", "f8"), ("y", "f8"), ("z", "f8")])
            new_data["x"] = pts[:, 0]
            new_data["y"] = pts[:, 1]
            new_data["z"] = pts[:, 2]
            return new_data

        # Guard: if mismatch, truncate
        for k in list(metadata.keys()):
            if len(metadata[k]) != num_points:
                min_len = min(len(metadata[k]), num_points)
                pts = pts[:min_len, :]
                num_points = min_len
                for kk in metadata:
                    metadata[kk] = metadata[kk][:min_len]
                break

        new_data = np.zeros(
            num_points,
            dtype=[
                ("x", "f8"),
                ("y", "f8"),
                ("z", "f8"),
                ("classification", "u1"),
                ("intensity", "u2"),
                ("return_number", "u1"),
                ("num_returns", "u1"),
            ],
        )

        new_data["x"] = pts[:, 0]
        new_data["y"] = pts[:, 1]
        new_data["z"] = pts[:, 2]
        new_data["classification"] = metadata["classification"]
        new_data["intensity"] = metadata["intensity"]
        new_data["return_number"] = metadata["return_number"]
        new_data["num_returns"] = metadata["num_returns"]

        return new_data

    def _save_as_laz(self, aligned_data):
        """
        Saves transformed LAS/LAZ with CRS applied.
        Uses point_format=3, LAS 1.4; keeps it simple and stable.
        """
        has_meta = ("classification" in aligned_data.dtype.names)

        header = laspy.LasHeader(point_format=3, version="1.4")
        # Use reasonable scales; adjust if you want finer than 1 cm
        header.scales = np.array([0.01, 0.01, 0.01], dtype=np.float64)
        header.offsets = np.min(
            np.vstack((aligned_data["x"], aligned_data["y"], aligned_data["z"])).T, axis=0
        )

        # Apply CRS
        try:
            if self.original_crs is not None:
                header.add_crs(self.original_crs)
            else:
                header.add_crs(pyproj.CRS.from_epsg(26917))
        except Exception:
            header.add_crs(pyproj.CRS.from_epsg(26917))

        las = laspy.LasData(header)
        las.x = aligned_data["x"]
        las.y = aligned_data["y"]
        las.z = aligned_data["z"]

        if has_meta:
            las.classification = aligned_data["classification"]
            las.intensity = aligned_data["intensity"]
            las.return_number = aligned_data["return_number"]
            las.num_returns = aligned_data["num_returns"]

        las.write(self.aligned_path)
        print(f"Saved aligned LAS file to {self.aligned_path} with CRS applied.")
