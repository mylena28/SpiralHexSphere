import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

DAT_FILE = 'spiral.dat'

# ── Read file ─────────────────────────────────────────────────────────────────
try:
    with open(DAT_FILE, 'r') as fh:
        lines = fh.readlines()

    header = lines[0].split()
    nv, nf, nv_spiral, pts_per_ring = (int(header[0]), int(header[1]),
                                        int(header[2]), int(header[3]))

    # Locate the 'v' and 'f' marker lines
    v_line = next(i for i, l in enumerate(lines) if l.strip() == 'v')
    f_line = next(i for i, l in enumerate(lines) if l.strip() == 'f')

    # Vertices: columns 0-2 = position, columns 3-5 = normal
    vert_data = np.array(
        [lines[v_line + 1 + i].split() for i in range(nv)], dtype=float)
    vertices = vert_data[:, :3]
    normals  = vert_data[:, 3:6]

    # Faces (1-based → 0-based)
    faces = np.array(
        [lines[f_line + 1 + i].split() for i in range(nf)], dtype=int) - 1

except FileNotFoundError:
    print(f"ERROR: '{DAT_FILE}' not found.  Run 'make run' first.")
    sys.exit(1)
except Exception as exc:
    print(f"ERROR reading file: {exc}")
    sys.exit(1)

print(f"Loaded  {nv} vertices,  {nf} faces")
print(f"  nv_spiral = {nv_spiral}")

# ── Split faces: spiral vs caps ───────────────────────────────────────────────
is_cap   = np.any(faces >= nv_spiral, axis=1)
spiral_f = faces[~is_cap]
cap_f    = faces[ is_cap]
print(f"  Spiral faces : {len(spiral_f)}")
print(f"  Cap faces    : {len(cap_f)}")

# ── Helpers ───────────────────────────────────────────────────────────────────
def make_tri_collection(f_idx, facecolor, alpha):
    tris = vertices[f_idx]          # (n, 3, 3)
    col  = Poly3DCollection(tris)
    col.set_facecolor(facecolor)
    col.set_edgecolor('black')
    col.set_linewidth(0.2)
    col.set_alpha(alpha)
    return col

def set_equal_axes(ax, pts):
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    mid    = 0.5 * (lo + hi)
    half   = 0.5 * (hi - lo).max()
    ax.set_xlim(mid[0] - half, mid[0] + half)
    ax.set_ylim(mid[1] - half, mid[1] + half)
    ax.set_zlim(mid[2] - half, mid[2] + half)

# ── Figure 1: 3-D mesh ────────────────────────────────────────────────────────
fig = plt.figure(figsize=(11, 9))
ax  = fig.add_subplot(111, projection='3d')
ax.set_title('SpiralHexSphere mesh')

if len(spiral_f) > 0:
    ax.add_collection3d(make_tri_collection(spiral_f, 'steelblue', 0.75))
if len(cap_f) > 0:
    ax.add_collection3d(make_tri_collection(cap_f, 'tomato', 0.85))

set_equal_axes(ax, vertices)
ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
plt.tight_layout()

# ── Figure 2: vertex normals (quiver, sub-sampled) ───────────────────────────
# fig2 = plt.figure(figsize=(11, 9))
# ax2  = fig2.add_subplot(111, projection='3d')
# ax2.set_title('Vertex normals (sub-sampled)')
#
# # Sub-sample to keep the plot readable
# step = max(1, nv // 800)
# idx  = np.arange(0, nv, step)
# scale = 0.08   # arrow length relative to geometry
# ax2.quiver(vertices[idx, 0], vertices[idx, 1], vertices[idx, 2],
#            normals[idx, 0],  normals[idx, 1],  normals[idx, 2],
#            length=scale, normalize=False, color='navy', linewidth=0.6, alpha=0.7)
#
# if len(spiral_f) > 0:
#     ax2.add_collection3d(make_tri_collection(spiral_f, 'steelblue', 0.25))
# if len(cap_f) > 0:
#     ax2.add_collection3d(make_tri_collection(cap_f, 'tomato', 0.35))
#
# set_equal_axes(ax2, vertices)
# ax2.set_xlabel('X'); ax2.set_ylabel('Y'); ax2.set_zlabel('Z')
# plt.tight_layout()

plt.show()
