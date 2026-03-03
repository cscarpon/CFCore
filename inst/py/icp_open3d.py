import open3d as o3d
import numpy as np
import laspy
import os
import pyproj

class Open3DICP:
    """
    ICP alignment using Open3D with LAS/LAZ IO via laspy.
    Corrected to align TARGET to SOURCE.
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
        use_gpu=True, 
    ):
        self.source_path = source_path
        self.target_path = target_path
        
        # Output should be based on the target filename, since we are moving the target
        self.aligned_path = self._generate_aligned_path(target_path)

        self.voxel_size = float(voxel_size)
        self.icp_method = str(icp_method)
        self.max_iteration = int(max_iteration)
        self.threshold_factor = float(threshold_factor)

        self.multiscale_factors = tuple(float(x) for x in multiscale_factors)
        self.multiscale_iters = tuple(int(x) for x in multiscale_iters)

        self.normal_radius_factor = float(normal_radius_factor)
        self.normal_max_nn = int(normal_max_nn)
        self.min_correspondences = int(min_correspondences)
        self.center = bool(center)

        self.transformation = None
        self.rmse = None
        self.original_crs = None
        self.original_header = None
        
        self.use_gpu = use_gpu
        self.is_cuda = False
        
        if self.use_gpu and hasattr(o3d.core, "cuda") and o3d.core.cuda.is_available():
            self.device = o3d.core.Device("CUDA:0")
            self.is_cuda = True
            print("Open3D ICP: CUDA is available. Using GPU Tensor API.")
        else:
            self.device = o3d.core.Device("CPU:0")
            if self.use_gpu:
                print("Open3D ICP: CUDA not available in this Python build. Falling back to CPU.")

    def align(self):
        source_pcd, _ = self._load_point_cloud(self.source_path)
        target_pcd, _ = self._load_point_cloud(self.target_path)

        if source_pcd is None or target_pcd is None:
            return None, "Error: Failed to load point clouds."

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
            
            # Store for un-centering later
            self._center_info = {"src_c": src_c, "tgt_c": tgt_c}

        try:
            # Get transformation mapping TARGET to SOURCE
            T, rmse = self._multiscale_icp(source_pcd, target_pcd)
        except Exception as e:
            return None, f"Error during ICP alignment: {str(e)}"

        if self.center:
            R = T[:3, :3]
            t = T[:3, 3]
            
            # Math to convert centered Target->Source mapping back to original coordinates
            src_c = self._center_info["src_c"]
            tgt_c = self._center_info["tgt_c"]
            t_final = t + src_c - R @ tgt_c

            T_final = np.eye(4, dtype=np.float64)
            T_final[:3, :3] = R
            T_final[:3, 3] = t_final
        else:
            T_final = T

        self.transformation = T_final
        self.rmse = float(rmse)

        # Reload the FULL TARGET point cloud (Epoch 2) and apply the transformation
        target_pcd_full, target_metadata_full = self._load_point_cloud(self.target_path)
        if target_pcd_full is None:
            return None, "Error: Failed to reload target point cloud."
        
        target_pcd_full.transform(self.transformation)

        # Save the transformed Target data
        aligned_data = self._merge_metadata(target_pcd_full, target_metadata_full)
        self._save_as_laz(aligned_data)

        output_str = (
            "ICP Alignment Completed (Target aligned to Source).\n"
            f"Method requested: {self.icp_method}\n"
            f"RMSE: {self.rmse:.6f}\n"
            "Transformation Matrix:\n"
            f"{np.array_str(self.transformation, precision=6, suppress_small=True)}"
        )

        return self.aligned_path, output_str

    def _multiscale_icp(self, source_pcd, target_pcd):
        T = np.eye(4, dtype=np.float64)
        last_rmse = None
        use_p2l = (self.icp_method.lower() == "point-to-plane")
        
        if self.is_cuda:
            t_source = o3d.t.geometry.PointCloud.from_legacy(source_pcd, device=self.device)
            t_target = o3d.t.geometry.PointCloud.from_legacy(target_pcd, device=self.device)
            T_tensor = o3d.core.Tensor(T, dtype=o3d.core.Dtype.Float64, device=self.device)

            for factor, iters in zip(self.multiscale_factors, self.multiscale_iters):
                vs = self.voxel_size * factor
                
                t_src = t_source.voxel_down_sample(vs)
                t_tgt = t_target.voxel_down_sample(vs)
                
                if t_src.point.positions.shape[0] < 100 or t_tgt.point.positions.shape[0] < 100:
                    continue
                
                if use_p2l:
                    rad = vs * self.normal_radius_factor
                    t_src.estimate_normals(max_nn=self.normal_max_nn, radius=rad)
                    t_tgt.estimate_normals(max_nn=self.normal_max_nn, radius=rad)

                threshold = vs * self.threshold_factor
                crit = o3d.t.pipelines.registration.ICPConvergenceCriteria(max_iteration=int(iters))
                est = (o3d.t.pipelines.registration.TransformationEstimationPointToPlane() 
                       if use_p2l else o3d.t.pipelines.registration.TransformationEstimationPointToPoint())

                # IMPORTANT: Align Target TO Source (t_tgt, t_src)
                result = o3d.t.pipelines.registration.icp(
                    t_tgt, t_src, threshold, T_tensor, est, crit
                )

                T_tensor = result.transformation
                last_rmse = float(result.inlier_rmse)

            T = T_tensor.cpu().numpy()
            
        else:
            for factor, iters in zip(self.multiscale_factors, self.multiscale_iters):
                vs = self.voxel_size * factor
                src = source_pcd.voxel_down_sample(vs)
                tgt = target_pcd.voxel_down_sample(vs)
                
                if len(src.points) < 100 or len(tgt.points) < 100:
                    continue

                def estimate_normals(pcd, voxel):
                    rad = voxel * self.normal_radius_factor
                    pcd.estimate_normals(
                        o3d.geometry.KDTreeSearchParamHybrid(radius=rad, max_nn=self.normal_max_nn)
                    )
                    pcd.normalize_normals()

                if use_p2l:
                    estimate_normals(src, vs)
                    estimate_normals(tgt, vs)

                threshold = vs * self.threshold_factor
                
                # IMPORTANT: Evaluate Target TO Source (tgt, src)
                eval0 = o3d.pipelines.registration.evaluate_registration(tgt, src, threshold, T)
                if eval0.correspondence_set is None or len(eval0.correspondence_set) < self.min_correspondences:
                    threshold2 = threshold * 2.0
                    eval1 = o3d.pipelines.registration.evaluate_registration(tgt, src, threshold2, T)
                    if eval1.correspondence_set is None or len(eval1.correspondence_set) < self.min_correspondences:
                        continue
                    threshold = threshold2

                crit = o3d.pipelines.registration.ICPConvergenceCriteria(max_iteration=int(iters))

                if use_p2l:
                    try:
                        # IMPORTANT: Register Target TO Source (tgt, src)
                        result = o3d.pipelines.registration.registration_icp(
                            tgt, src, threshold, T,
                            o3d.pipelines.registration.TransformationEstimationPointToPlane(), crit
                        )
                    except RuntimeError:
                        result = o3d.pipelines.registration.registration_icp(
                            tgt, src, threshold, T,
                            o3d.pipelines.registration.TransformationEstimationPointToPoint(), crit
                        )
                else:
                    result = o3d.pipelines.registration.registration_icp(
                        tgt, src, threshold, T,
                        o3d.pipelines.registration.TransformationEstimationPointToPoint(), crit
                    )

                T = result.transformation
                last_rmse = float(result.inlier_rmse)

        if last_rmse is None:
            raise RuntimeError("ICP failed: no scale produced enough correspondences.")

        return T, last_rmse

    def _load_point_cloud(self, file_path):
        if file_path.endswith(".las") or file_path.endswith(".laz"):
            with laspy.open(file_path) as las_file:
                las = las_file.read()

                points = np.vstack((las.x, las.y, las.z)).T.astype(np.float64)
                pcd = o3d.geometry.PointCloud()
                pcd.points = o3d.utility.Vector3dVector(points)

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

    def _generate_aligned_path(self, target_path):
        base, _ = os.path.splitext(target_path)
        return f"{base}_aligned.laz"

    def _merge_metadata(self, transformed_pcd, metadata):
        pts = np.asarray(transformed_pcd.points)
        num_points = pts.shape[0]

        if metadata is None:
            new_data = np.zeros(num_points, dtype=[("x", "f8"), ("y", "f8"), ("z", "f8")])
            new_data["x"] = pts[:, 0]
            new_data["y"] = pts[:, 1]
            new_data["z"] = pts[:, 2]
            return new_data

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
        has_meta = ("classification" in aligned_data.dtype.names)
        header = laspy.LasHeader(point_format=3, version="1.4")
        header.scales = np.array([0.01, 0.01, 0.01], dtype=np.float64)
        header.offsets = np.min(
            np.vstack((aligned_data["x"], aligned_data["y"], aligned_data["z"])).T, axis=0
        )

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
