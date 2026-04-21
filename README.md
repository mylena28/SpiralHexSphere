# SpiralHexSphere

A Fortran 90 mesh generator for a helical tubular surface capped with two hemispherical ends. The geometry follows the Zinchenko 1997 construction and produces a triangulated surface mesh forming a spiral with caps at the end.

---

## Overview

The code builds a closed triangular surface mesh in three parts:

1. **Spiral tube** — a helical tube swept along a 3D helix centerline using the Frenet-Serret frame. The cross-section polygon is swept and subdivided using a Zinchenko Level-1 refinement, producing a fine hexagonal/triangular tiling.
2. **Start hemisphere cap** — a hemispherical endcap closing the first end of the tube.
3. **End hemisphere cap** — a hemispherical endcap closing the second end.

The result is a watertight, consistently-oriented surface mesh written to a plain-text `.dat` file.

---

## Repository Structure

```
SpiralHexSphere/
├── main.f90          # Entry point; orchestrates mesh generation and output
├── params.f90        # Global parameters module (geometry and mesh settings)
├── geometry.f90      # Core algorithms: Frenet-Serret frame, vertex generation, triangulation
├── io.f90            # Output routines; writes vertex/face data to file
├── plot_mesh.py      # Python visualization script (matplotlib 3D)
└── Makefile          # Build, run, and plot targets
```

---

## Requirements

### Fortran (mesh generation)
- `gfortran` (GCC Fortran compiler, version 7 or later)

### Python (visualization, optional)
- Python 3.x
- NumPy
- Matplotlib

---

## Building and Running

### Compile
```bash
make
```

### Generate the mesh
```bash
make run
# or directly:
./spiral_hex_sphere
```

This writes the mesh to `spiral_hex_sphere.dat`.

### Visualize the mesh
```bash
make plot
# or directly:
python3 plot_mesh.py
```

Spiral vertices are shown in blue; hemisphere cap vertices in red.

### Compile, run, and plot in one step
```bash
make all_run
```

### Clean build artifacts
```bash
make clean
```

---

## Configuration

All parameters are set in `params.f90`. Recompile after any changes.

| Parameter  | Default | Description |
|------------|---------|-------------|
| `helix_R`  | `1.0`   | Major radius of the helix (distance from the helix axis) |
| `helix_h`  | `6.0`   | Total axial height of the helix |
| `n_turns`  | `4.0`   | Number of full 360° turns |
| `tube_r`   | `0.3`   | Cross-sectional radius of the tube |
| `num_sides`| `3`     | Number of polygon sides per cross-section ring (3 → hexagonal tiling after subdivision) |
| `n_s`      | `144`   | Number of axial sections along the helix (controls spiral resolution) |
| `n_hemi`   | `6`     | Number of latitude rings per hemisphere cap |
| `outfile`  | `'spiral_hex_sphere.dat'` | Output filename |

### Derived mesh sizes (default parameters)

| Quantity | Value |
|----------|-------|
| Rings along spiral | 287 |
| Vertices per ring | 6 |
| Spiral vertices | 1,722 |
| Total vertices | ~2,500 |
| Total triangular faces | ~5,000 |

---

## Output File Format

The file `spiral_hex_sphere.dat` uses the following plain-text layout:

```
nv  nf  nv_spiral
x1  y1  z1
x2  y2  z2
...                   ← nv vertex lines (double precision)
i1  j1  k1
i2  j2  k2
...                   ← nf face lines (1-based triangle indices)
```

- `nv` — total number of vertices
- `nf` — total number of triangular faces
- `nv_spiral` — number of vertices belonging to the spiral portion (the remaining vertices are the two hemisphere caps)
- Face indices follow Fortran 1-based convention and are ordered for outward-pointing normals

---

## Algorithms

### Frenet-Serret Frame
At each axial position along the helix, the code computes a local orthonormal frame:
- **T** — unit tangent to the centerline
- **N** — principal normal (points toward the helix axis)
- **B = T × N** — binormal

Tube vertices are then placed at:
```
v = P + tube_r * (N * cos(θ) + B * sin(θ))
```
where `θ` varies uniformly around the cross-section.

### Zinchenko Level-1 Subdivision
Each 2×2 block of neighboring rings is subdivided into 8 triangles by inserting edge midpoints. This produces the characteristic hexagonal/triangular tiling associated with the Zinchenko mesh construction.

### Hemisphere Caps
Each cap is built as a series of latitude rings whose radius shrinks from `tube_r` at the equator to zero at the pole. Rings are triangulated by connecting adjacent latitude rings with two triangles per quad. A fan of triangles closes the pole.

---

## License

Please contact the author before reuse or redistribution.
