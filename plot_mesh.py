import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Patch

DAT_FILE = 'spiral_hex_sphere.dat'

# Read file
try:
    with open(DAT_FILE, 'r') as fh:
        header = fh.readline().split()
        nv, nf, nv_spiral = int(header[0]), int(header[1]), int(header[2])

        vertices = np.array(
            [fh.readline().split() for _ in range(nv)], dtype=float)

        faces = np.array(
            [fh.readline().split() for _ in range(nf)], dtype=int) - 1  # 0-indexed

except FileNotFoundError:
    print(f"ERROR: '{DAT_FILE}' not found.  Run 'make run' first.")
    sys.exit(1)
except Exception as exc:
    print(f"ERROR reading file: {exc}")
    sys.exit(1)

print(f"Loaded  {nv} vertices,  {nf} faces")
print(f"  nv_spiral = {nv_spiral}  (caps use vertices >= {nv_spiral})")

# Split faces into spiral / hemisphere caps
is_cap   = np.any(faces >= nv_spiral, axis=1)
spiral_f = faces[~is_cap]
cap_f    = faces[ is_cap]

print(f"  Spiral faces : {len(spiral_f)}")
print(f"  Cap faces    : {len(cap_f)}")

# Build triangle collections
def make_tri_collection(f_idx, facecolor, alpha):
    tris = vertices[f_idx]          # (n, 3, 3) array of triangle vertex coords
    col  = Poly3DCollection(tris)
    col.set_facecolor(facecolor)
    col.set_edgecolor('black')
    col.set_linewidth(0.2)
    col.set_alpha(alpha)
    return col

# Plot
fig = plt.figure(figsize=(11, 9))
ax  = fig.add_subplot(111, projection='3d')

if len(spiral_f) > 0:
    ax.add_collection3d(make_tri_collection(spiral_f, 'steelblue', 0.75))
if len(cap_f) > 0:
    ax.add_collection3d(make_tri_collection(cap_f, 'tomato', 0.85))

# Equal-aspect bounding box
lo, hi = vertices.min(axis=0), vertices.max(axis=0)
mid    = 0.5 * (lo + hi)
half   = 0.5 * (hi - lo).max()

ax.set_xlim(mid[0] - half, mid[0] + half)
ax.set_ylim(mid[1] - half, mid[1] + half)
ax.set_zlim(mid[2] - half, mid[2] + half)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')


plt.tight_layout()
plt.show()
