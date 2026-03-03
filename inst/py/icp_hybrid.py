import numpy as np
import laspy
import os
import pyproj

# --- HYBRID IMPORT LOGIC ---
# Try to load GPU-accelerated Cupoch first
try:
    import cupoch as cph
    HAS_CUPOCH = True
except ImportError:
    HAS_CUPOCH = False

# Always load Open3D as the fallback
import open3d as o3d

class HybridICP:
    """
    Hybrid ICP alignment for CFCore.
    Attempts to use Cupoch (CUDA GPU) for maximum speed. 
    Gracefully falls back to Open3D (CPU) if a GPU or Cupoch is unavailable.
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
        use_gpu=True, 
    ):
        self.source_path = source_path
        self.target_path = target_path
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
        
        self.transformation = None
        self.rmse = None
        self.original_crs = None
        self.original_header = None
        
        # Determine Execution Engine
        if use_gpu and HAS_CUPOCH:
            self.engine = "CUPOCH_GPU"
            print("CFCore ICP: Cupoch detected. Using CUDA GPU Acceleration.")
        else:
            self.engine = "OPEN3D_CPU"
            if use_gpu and not HAS_CUPOCH:
                print("CFCore ICP: GPU requested, but Cupoch is not installed. Falling back to Open3D CPU.")
            else:
                print("CFCore ICP: Using Open3D CPU processing.")

    def align(self):
        # 1. Load raw data
        src_pts, src_metadata = self._load_las_data(self.source_path)
        tgt_pts, tgt_metadata = self._load_las_data(self.target_path)

        if src_pts is None or tgt_pts is None:
            return None, "Error: Failed to load point clouds."

        # 2. Universal Centering (Required for both engines to maintain precision)
        src_c = src_pts.mean(axis=0)
        tgt_c = tgt_pts.mean(axis=0)
        
        src_pts_centered = src_pts - src_c
        tgt_pts_centered = tgt_pts - tgt_c

        # 3. Route to the correct execution engine
        try:
            if self.engine == "CUPOCH_GPU":
                T_centered, rmse = self._run_cupoch_gpu(src_pts_centered, tgt_pts_centered)
            else:
                T_centered, rmse = self._run_open3d_cpu(src_pts_centered, tgt_pts_centered)
        except Exception as e:
            return None, f"Error during ICP alignment ({self.engine}): {str(e)}"

        # 4. Uncenter the transformation matrix
        R = T_centered[:3, :3]
        t = T_centered[:3, 3]
        t_final = t + src_c - R @ tgt_c

        T_final = np.eye(4, dtype=np.float64)
        T_final[:3, :3] = R
        T_final[:3, 3] = t_final

        self.transformation = T_final
        self.rmse = float(rmse)

        # 5. Apply transformation to the ORIGINAL target points
        tgt_pts_transformed = (tgt_pts @ R.T) + t_final
        
        # 6. Save outputs
        aligned_data = self._merge_metadata(tgt_pts_transformed, tgt_metadata)
        self._save_as_laz(aligned_data)

        output_str = (
            f"ICP Alignment Completed via {self.engine}.\n"
            f"Method: {self.icp_method} | RMSE: {self.rmse:.6f}\n"
        )
        return self.aligned_path, output_str

    # --- CUPOCH GPU ENGINE ---
    def _run_cupoch_gpu(self, src_pts_centered, tgt_pts_centered):
        source_pcd = cph.geometry.PointCloud()
        source_pcd.points = cph.utility.Vector3fVector(src_pts_centered.astype(np.float32))
        
        target_pcd = cph.geometry.PointCloud()
        target_pcd.points = cph.utility.Vector3fVector(tgt_pts_centered.astype(np.float32))

        T = np.eye(4, dtype=np.float32)
        last_rmse = None
        use_p2l = (self.icp_method.lower() == "point-to-plane")
        
        for factor, iters in zip(self.multiscale_factors, self.multiscale_iters):
            vs = self.voxel_size * factor
            src = source_pcd.voxel_down_sample(vs)
            tgt = target_pcd.voxel_down_sample(vs)
            
            src_points = np.asarray(src.points.cpu()) if hasattr(src.points, 'cpu') else np.asarray(src.points)
            tgt_points = np.asarray(tgt.points.cpu()) if hasattr(tgt.points, 'cpu') else np.asarray(tgt.points)
            
            if len(src_points) < 100 or len(tgt_points) < 100:
                continue
            
            if use_p2l:
                search_param = cph.geometry.KDTreeSearchParamKNN(knn=int(self.normal_max_nn))
                src.estimate_normals(search_param)
                tgt.estimate_normals(search_param)

            threshold = vs * self.threshold_factor
            crit = cph.registration.ICPConvergenceCriteria(max_iteration=int(iters))

            if use_p2l:
                try:
                    result = cph.registration.registration_icp(
                        tgt, src, threshold, T,
                        cph.registration.TransformationEstimationPointToPlane(), crit
                    )
                except RuntimeError:
                    result = cph.registration.registration_icp(
                        tgt, src, threshold, T, cph.registration.TransformationEstimationPointToPoint(), crit
                    )
            else:
                result = cph.registration.registration_icp(
                    tgt, src, threshold, T, cph.registration.TransformationEstimationPointToPoint(), crit
                )

            T = result.transformation
            last_rmse = float(result.inlier_rmse)

        if last_rmse is None:
            raise RuntimeError("ICP failed.")
            
        if hasattr(T, 'cpu'):
            T_out = np.asarray(T.cpu()).astype(np.float64)
        else:
            T_out = np.asarray(T).astype(np.float64)
            
        return T_out, last_rmse

    # --- OPEN3D CPU ENGINE ---
    def _run_open3d_cpu(self, src_pts_centered, tgt_pts_centered):
        source_pcd = o3d.geometry.PointCloud()
        source_pcd.points = o3d.utility.Vector3dVector(src_pts_centered.astype(np.float64))
        
        target_pcd = o3d.geometry.PointCloud()
        target_pcd.points = o3d.utility.Vector3dVector(tgt_pts_centered.astype(np.float64))

        T = np.eye(4, dtype=np.float64)
        last_rmse = None
        use_p2l = (self.icp_method.lower() == "point-to-plane")

        for factor, iters in zip(self.multiscale_factors, self.multiscale_iters):
            vs = self.voxel_size * factor
            src = source_pcd.voxel_down_sample(vs)
            tgt = target_pcd.voxel_down_sample(vs)

            if len(src.points) < 100 or len(tgt.points) < 100:
                continue

            if use_p2l:
                src.estimate_normals(o3d.geometry.KDTreeSearchParamHybrid(radius=vs*self.normal_radius_factor, max_nn=self.normal_max_nn))
                tgt.estimate_normals(o3d.geometry.KDTreeSearchParamHybrid(radius=vs*self.normal_radius_factor, max_nn=self.normal_max_nn))

            threshold = vs * self.threshold_factor
            crit = o3d.pipelines.registration.ICPConvergenceCriteria(max_iteration=int(iters))

            if use_p2l:
                try:
                    result = o3d.pipelines.registration.registration_icp(
                        tgt, src, threshold, T, o3d.pipelines.registration.TransformationEstimationPointToPlane(), crit
                    )
                except RuntimeError:
                    result = o3d.pipelines.registration.registration_icp(
                        tgt, src, threshold, T, o3d.pipelines.registration.TransformationEstimationPointToPoint(), crit
                    )
            else:
                result = o3d.pipelines.registration.registration_icp(
                    tgt, src, threshold, T, o3d.pipelines.registration.TransformationEstimationPointToPoint(), crit
                )

            T = result.transformation
            last_rmse = float(result.inlier_rmse)

        if last_rmse is None:
            raise RuntimeError("ICP failed.")

        return T, last_rmse

    # --- IO METHODS (Same for both) ---
    def _load_las_data(self, file_path):
        if file_path.endswith(".las") or file_path.endswith(".laz"):
            with laspy.open(file_path) as las_file:
                las = las_file.read()
                points = np.vstack((las.x, las.y, las.z)).T.astype(np.float64)

                self.original_header = las.header
                self.original_crs = self._parse_crs_safe(las.header)

                metadata = {
                    "classification": las.classification.copy(),
                    "intensity": las.intensity.copy(),
                    "return_number": las.return_number.copy(),
                    "num_returns": las.num_returns.copy(),
                }
                return points, metadata
        return None, None

    def _parse_crs_safe(self, header):
        try:
            crs = header.parse_crs()
            if crs is None: return pyproj.CRS.from_epsg(26917)
            if isinstance(crs, pyproj.CRS):
                if crs.is_compound: return crs.sub_crs_list[0]
                return crs
            return pyproj.CRS.from_user_input(crs)
        except Exception:
            return pyproj.CRS.from_epsg(26917)

    def _generate_aligned_path(self, target_path):
        base, _ = os.path.splitext(target_path)
        return f"{base}_aligned.laz"

    def _merge_metadata(self, transformed_pts, metadata):
        num_points = transformed_pts.shape[0]

        if metadata is None:
            new_data = np.zeros(num_points, dtype=[("x", "f8"), ("y", "f8"), ("z", "f8")])
            new_data["x"] = transformed_pts[:, 0]
            new_data["y"] = transformed_pts[:, 1]
            new_data["z"] = transformed_pts[:, 2]
            return new_data

        for k in list(metadata.keys()):
            if len(metadata[k]) != num_points:
                min_len = min(len(metadata[k]), num_points)
                transformed_pts = transformed_pts[:min_len, :]
                num_points = min_len
                for kk in metadata:
                    metadata[kk] = metadata[kk][:min_len]
                break

        new_data = np.zeros(
            num_points,
            dtype=[("x", "f8"), ("y", "f8"), ("z", "f8"),
                   ("classification", "u1"), ("intensity", "u2"),
                   ("return_number", "u1"), ("num_returns", "u1")],
        )
        new_data["x"] = transformed_pts[:, 0]
        new_data["y"] = transformed_pts[:, 1]
        new_data["z"] = transformed_pts[:, 2]
        new_data["classification"] = metadata["classification"]
        new_data["intensity"] = metadata["intensity"]
        new_data["return_number"] = metadata["return_number"]
        new_data["num_returns"] = metadata["num_returns"]

        return new_data

    def _save_as_laz(self, aligned_data):
        has_meta = ("classification" in aligned_data.dtype.names)
        
        # --- THE FIX IS HERE ---
        # Dynamically match the original file's format to support >7 returns
        if self.original_header is not None:
            header = laspy.LasHeader(
                point_format=self.original_header.point_format, 
                version=self.original_header.version
            )
        else:
            header = laspy.LasHeader(point_format=6, version="1.4") # Format 6 supports up to 15 returns
        # -----------------------

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
