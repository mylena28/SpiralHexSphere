import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

OBJ_FILE = 'spiral.obj'

try:
    vertices = []
    normals  = []
    faces    = []

    with open(OBJ_FILE, 'r') as fh:
        for line in fh:
            parts = line.split()
            if not parts:
                continue
            if parts[0] == 'v':
                vertices.append([float(x) for x in parts[1:4]])
            elif parts[0] == 'vn':
                normals.append([float(x) for x in parts[1:4]])
            elif parts[0] == 'f':
                # face entries are v//vn; take the vertex index (1-based → 0-based)
                faces.append([int(p.split('//')[0]) - 1 for p in parts[1:4]])

    vertices = np.array(vertices, dtype=float)
    normals  = np.array(normals,  dtype=float)
    faces    = np.array(faces,    dtype=int)

except FileNotFoundError:
    print(f"ERROR: '{OBJ_FILE}' not found.  Run 'make run' first.")
    sys.exit(1)
except Exception as exc:
    print(f"ERROR reading file: {exc}")
    sys.exit(1)

nv = len(vertices)
nf = len(faces)
print(f"Loaded  {nv} vertices,  {nf} faces")

def set_equal_axes(ax, pts):
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    mid    = 0.5 * (lo + hi)
    half   = 0.5 * (hi - lo).max()
    ax.set_xlim(mid[0] - half, mid[0] + half)
    ax.set_ylim(mid[1] - half, mid[1] + half)
    ax.set_zlim(mid[2] - half, mid[2] + half)

fig = plt.figure(figsize=(11, 9))
ax  = fig.add_subplot(111, projection='3d')
ax.set_title('SpiralHexSphere mesh')

tris = vertices[faces]
col  = Poly3DCollection(tris)
col.set_facecolor('steelblue')
col.set_edgecolor('black')
col.set_linewidth(0.2)
col.set_alpha(0.75)
ax.add_collection3d(col)

set_equal_axes(ax, vertices)
ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
plt.tight_layout()
plt.show()
