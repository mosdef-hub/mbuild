"""GPU (CUDA) equivalents of path_utils numba functions.

This module provides drop-in replacements for the @njit functions in path_utils.py,
using numba's cuda.jit so that work runs on the GPU.

Notes:
- Device functions (norm, rotate_vector) are compiled with device=True and run
  only on the device; they are used by kernels or other device code.
- Kernels write results into device arrays. The
  public API functions allocate device arrays, launch the kernel, and copy
  results back to the host so the interface matches path_utils.
"""

import math

import numpy as np
from numba import cuda


@cuda.jit(device=True)
def _norm(vec):
    """Vector 2-norm. Replacement for path_utils.norm."""
    s = 0.0
    for i in range(vec.shape[0]):
        s += vec[i] * vec[i]
    return math.sqrt(s)


@cuda.jit
def _check_path_split_kernel(points, candidates, min_sq_dist, valid):
    """Check candidates against points. Sets valid[i]=0 if candidate i contains overlaps."""
    cand_i = cuda.blockIdx.x
    point_i = cuda.threadIdx.x + cuda.blockIdx.y * cuda.blockDim.x

    if cand_i >= candidates.shape[0]:
        return
    if point_i >= points.shape[0]:
        return

    # Skip if already marked invalid
    if valid[cand_i] == 0:
        return

    dx = points[point_i, 0] - candidates[cand_i, 0]
    dy = points[point_i, 1] - candidates[cand_i, 1]
    dz = points[point_i, 2] - candidates[cand_i, 2]
    dist_sq = dx * dx + dy * dy + dz * dz

    if dist_sq < min_sq_dist:
        cuda.atomic.min(valid, cand_i, 0)


@cuda.jit
def _check_path_kernel(existing_points, new_point, min_sq_dist, collision_out):
    """One thread per existing point. If any dist_sq < min_sq_dist, set collision_out[0] = 1."""
    i = cuda.grid(1)
    n = existing_points.shape[0]
    if i >= n:
        return
    dx = existing_points[i, 0] - new_point[0]
    dy = existing_points[i, 1] - new_point[1]
    dz = existing_points[i, 2] - new_point[2]
    dist_sq = dx * dx + dy * dy + dz * dz
    if dist_sq < min_sq_dist:
        cuda.atomic.max(collision_out, 0, 1)


@cuda.jit
def _target_sq_distances_kernel(
    target_coordinate, new_points, pbc, box_lengths, sq_distances
):
    """One thread per new_point: compute squared distance with PBC and write to sq_distances."""
    i = cuda.grid(1)
    n = new_points.shape[0]
    if i >= n:
        return
    dx = target_coordinate[0] - new_points[i, 0]
    dy = target_coordinate[1] - new_points[i, 1]
    dz = target_coordinate[2] - new_points[i, 2]
    if pbc[0]:
        dx -= np.round(dx / box_lengths[0]) * box_lengths[0]
    if pbc[1]:
        dy -= np.round(dy / box_lengths[1]) * box_lengths[1]
    if pbc[2]:
        dz -= np.round(dz / box_lengths[2]) * box_lengths[2]
    sq_distances[i] = dx * dx + dy * dy + dz * dz


@cuda.jit
def _local_density_kernel(candidate, target_coords, r2_cut, density_out):
    """One thread per target. Each thread adds 1 to density_out[0] if within r_cut (atomic add)."""
    i = cuda.grid(1)
    n = target_coords.shape[0]
    if i >= n:
        return
    dx = candidate[0] - target_coords[i, 0]
    dy = candidate[1] - target_coords[i, 1]
    dz = candidate[2] - target_coords[i, 2]
    dist2 = dx * dx + dy * dy + dz * dz
    if dist2 < r2_cut:
        cuda.atomic.add(density_out, 0, 1)


@cuda.jit
def _target_density_kernel(candidates, target_coords, r2_cut, out):
    """One block per candidate. Each block runs _local_density logic and writes to out[blockIdx.x]."""
    cand_i = cuda.blockIdx.x
    t = cuda.threadIdx.x
    n_targets = target_coords.shape[0]
    n_threads = cuda.blockDim.x
    count = 0
    while t < n_targets:
        dx = candidates[cand_i, 0] - target_coords[t, 0]
        dy = candidates[cand_i, 1] - target_coords[t, 1]
        dz = candidates[cand_i, 2] - target_coords[t, 2]
        if dx * dx + dy * dy + dz * dz < r2_cut:
            count += 1
        t += n_threads
    shared = cuda.shared.array(256, dtype=np.int32)
    tid = cuda.threadIdx.x
    shared[tid] = count
    cuda.syncthreads()
    # Single-thread reduction over shared[0 : n_threads]
    if tid == 0:
        total = 0
        for k in range(n_threads):
            total += shared[k]
        out[cand_i] = np.float32(total)


def check_path_split(d_static_points, dynamic_points, candidates, radius, tolerance):
    """
    Check candidates against both static (already on GPU) and dynamic points.

    Parameters
    ----------
    d_static_points : cuda device array
        Static points already on GPU (e.g., previous chains). Do NOT transfer.
    dynamic_points : numpy array
        Dynamic points on CPU (e.g., current chain so far). Will be transferred.
    candidates : numpy array
        Candidate points to check (Nx3).
    radius : float
        Collision radius.
    tolerance : float
        Tolerance for collision detection.

    Returns
    -------
    numpy array of bool
        Mask where True = candidate is valid (no collision).
    """
    candidates = np.asarray(candidates, dtype=np.float32)
    dynamic_points = np.asarray(dynamic_points, dtype=np.float32)
    min_sq_dist = np.float32(radius - tolerance) ** np.float32(2.0)

    n_candidates = candidates.shape[0]
    if n_candidates == 0:
        return np.zeros(0, dtype=bool)

    n_static = d_static_points.shape[0] if d_static_points is not None else 0
    n_dynamic = dynamic_points.shape[0]

    # Initialize all as valid
    valid = np.ones(n_candidates, dtype=np.int32)
    d_valid = cuda.to_device(valid)
    d_candidates = cuda.to_device(candidates)

    threads_per_block = 256

    # Check against static points (already on GPU)
    if n_static > 0:
        blocks_y = (n_static + threads_per_block - 1) // threads_per_block
        blocks = (n_candidates, blocks_y)
        _check_path_split_kernel[blocks, threads_per_block](
            d_static_points, d_candidates, min_sq_dist, d_valid
        )

    # Check against dynamic points (transfer to GPU)
    if n_dynamic > 0:
        d_dynamic = cuda.to_device(dynamic_points)
        blocks_y = (n_dynamic + threads_per_block - 1) // threads_per_block
        blocks = (n_candidates, blocks_y)
        _check_path_split_kernel[blocks, threads_per_block](
            d_dynamic, d_candidates, min_sq_dist, d_valid
        )

    # Copy result back
    d_valid.copy_to_host(valid)
    return valid.astype(bool)


def check_path(existing_points, new_point, radius, tolerance):
    """GPU version of path_utils.check_path. Returns True if path is valid (no overlaps)."""
    existing_points = np.asarray(existing_points, dtype=np.float32)
    new_point = np.asarray(new_point, dtype=np.float32)
    min_sq_dist = np.float32(radius - tolerance) ** np.float32(2.0)

    n = existing_points.shape[0]
    collision_out = np.zeros(1, dtype=np.int32)

    d_existing = cuda.to_device(existing_points)
    d_new = cuda.to_device(new_point)
    d_collision = cuda.to_device(collision_out)

    threads_per_block = 256
    blocks = (n + threads_per_block - 1) // threads_per_block
    _check_path_kernel[blocks, threads_per_block](
        d_existing, d_new, min_sq_dist, d_collision
    )
    d_collision.copy_to_host(collision_out)
    return collision_out[0] == 0


def target_sq_distances(
    target_coordinate,
    new_points,
    pbc=None,
    box_lengths=None,
):
    """GPU version of path_utils.target_sq_distances. Same default pbc/box_lengths when omitted."""
    if pbc is None:
        pbc = np.array([False, False, False], dtype=bool)
    if box_lengths is None:
        box_lengths = np.array([np.inf, np.inf, np.inf], dtype=np.float32)

    target_coordinate = np.asarray(target_coordinate, dtype=np.float32)
    new_points = np.asarray(new_points, dtype=np.float32)
    pbc = np.asarray(pbc, dtype=bool)
    box_lengths = np.asarray(box_lengths, dtype=np.float32)

    n = new_points.shape[0]
    sq_distances = np.empty(n, dtype=np.float32)

    d_target = cuda.to_device(target_coordinate)
    d_new = cuda.to_device(new_points)
    d_pbc = cuda.to_device(pbc)
    d_box = cuda.to_device(box_lengths)
    d_out = cuda.device_array(n, dtype=np.float32)

    threads_per_block = 256
    blocks = (n + threads_per_block - 1) // threads_per_block
    _target_sq_distances_kernel[blocks, threads_per_block](
        d_target, d_new, d_pbc, d_box, d_out
    )
    d_out.copy_to_host(sq_distances)
    return sq_distances


def local_density(candidate, target_coords, r_cut):
    """GPU version of path_utils.local_density. Returns count of targets within r_cut."""
    candidate = np.asarray(candidate, dtype=np.float32)
    target_coords = np.asarray(target_coords, dtype=np.float32)
    r2_cut = np.float32(r_cut) * np.float32(r_cut)

    density_out = np.zeros(1, dtype=np.int32)
    n = target_coords.shape[0]

    d_candidate = cuda.to_device(candidate)
    d_targets = cuda.to_device(target_coords)
    d_out = cuda.to_device(density_out)

    threads_per_block = 256
    blocks = (n + threads_per_block - 1) // threads_per_block
    _local_density_kernel[blocks, threads_per_block](
        d_candidate, d_targets, r2_cut, d_out
    )
    d_out.copy_to_host(density_out)
    return density_out[0]


def target_density(candidates, target_coords, r_cut):
    """GPU version of path_utils.target_density. One block per candidate."""
    candidates = np.asarray(candidates, dtype=np.float32)
    target_coords = np.asarray(target_coords, dtype=np.float32)
    r2_cut = np.float32(r_cut) * np.float32(r_cut)

    n = candidates.shape[0]
    out = np.empty(n, dtype=np.float32)

    d_candidates = cuda.to_device(candidates)
    d_targets = cuda.to_device(target_coords)
    d_out = cuda.device_array(n, dtype=np.float32)

    threads_per_block = 256
    _target_density_kernel[n, threads_per_block](d_candidates, d_targets, r2_cut, d_out)
    d_out.copy_to_host(out)
    return out


def norm(vec):
    """Not launched from host; provided for API parity. Use on CPU or call from device code.
    For host-side use this runs on CPU (numpy). For device-side use, use _norm in kernels.
    All math is done in float32.
    """
    vec = np.asarray(vec, dtype=np.float32)
    return float(np.linalg.norm(vec).astype(np.float32))


def rotate_vector(v, axis, theta):
    """GPU-capable only when called from device. Host version using numpy for API parity."""
    v = np.asarray(v, dtype=np.float32)
    axis = np.asarray(axis, dtype=np.float32)
    theta = np.float32(theta)
    c = np.cos(theta)
    s = np.sin(theta)
    k_dot_v = np.dot(axis, v)
    cross = np.cross(axis, v)
    return (v * c + cross * s + axis * k_dot_v * (np.float32(1.0) - c)).astype(
        np.float32
    )


def _prepare_gpu_static_points(self):
    """Transfer static points to GPU once at the start of generation."""
    if not self.run_on_gpu:
        return

    from numba import cuda

    static_parts = []

    # Previous path coordinates (from start_from_path)
    if self._init_count > 0:
        static_parts.append(self.coordinates[: self._init_count])

    # Include compound coordinates
    if self.include_compound is not None:
        static_parts.append(self.include_compound.xyz)

    # Combine and transfer to GPU
    if static_parts:
        static_points = np.concatenate(static_parts).astype(np.float32)
        self._gpu_static_points = cuda.to_device(static_points)
    else:
        self._gpu_static_points = None
