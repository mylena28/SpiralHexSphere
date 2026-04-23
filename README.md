# SpiralHexSphere

Fortran 90 mesh generator for a helical tube capped with two hemispherical ends, producing a watertight triangulated surface mesh.

---

## Repository Structure

```
SpiralHexSphere/
├── params.f90        # Global parameters (geometry and refinement settings)
├── geometry.f90      # Mesh generation: Frenet frame, subdivision, projection
├── main.f90          # Entry point: orchestrates generation and output
├── io.f90            # Writes mesh to .dat file
├── plot_mesh.py      # Python visualization (matplotlib)
└── Makefile
```

---

## Quick Start

```bash
make            # compile
make run        # generate mesh  →  spiral_hex_sphere.dat
make plot       # visualize with Python
make all_run    # compile + run + plot in one step
make clean
```

---

## Parameters (`params.f90`)

All parameters require recompilation after changes.

| Parameter       | Default | Description |
|-----------------|---------|-------------|
| `helix_R`       | `1.0`   | Helix radius (distance from axis to centerline) |
| `helix_h`       | `6.0`   | Total axial height |
| `n_turns`       | `4.0`   | Number of full turns |
| `tube_r`        | `0.3`   | Tube cross-section radius |
| `num_sides`     | `6`     | Vertices per cross-section ring |
| `refine_level`  | `1`     | Subdivision levels (each level splits every triangle into 4) |
| `outfile`       | `'spiral_hex_sphere.dat'` | Output filename |

`n_s` (number of axial rings) is computed at runtime from the equilateral-triangle condition:
```
n_s = round(helix_R * n_turns * num_sides / tube_r)
```

---

## Output File Format

```
nv  nf  nv_spiral  pts_per_ring
x1  y1  z1
...                              ← nv vertex coordinates (double precision)
i1  j1  k1
...                              ← nf triangle faces (1-based indices)
```

- `nv_spiral` — number of tube-body vertices (remaining vertices belong to the two caps)
- Face indices give outward-pointing normals

---

## Algorithms

### Frenet-Serret frame
At each axial station along the helix, a local orthonormal frame **{T, N, B}** is computed. Ring vertices are placed at:
```
v = P + tube_r * (N·cos θ + B·sin θ)
```

### 4-1 midpoint subdivision
Each triangle is split into 4 children by inserting the midpoint of every edge. After subdivision, every new midpoint is **projected back onto the exact surface**:
- *Tube edges* — the chord midpoint is projected onto the tube's circular cross-section at the interpolated arc parameter, by removing the tangential component and renormalising in the N-B plane.
- *Cap edges* — the chord midpoint is projected onto the hemisphere sphere of radius `tube_r`.

This keeps cross-sections exactly circular and caps exactly spherical at every refinement level.

### Hemisphere caps
Each cap uses 3 interior vertices placed at an optimal elevation (`sin φ = (√7 − 1)/3`) to minimise edge-length variation, giving 10 triangles per cap. These interior points lie on the hemisphere sphere and are projected correctly during subdivision.
