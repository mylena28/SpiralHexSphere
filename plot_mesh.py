import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import Patch
from matplotlib.collections import PolyCollection

DAT_FILE = 'spiral_hex_sphere.dat'

# Read file
try:
    with open(DAT_FILE, 'r') as fh:
        header = fh.readline().split()
        nv, nf, nv_spiral, pts_per_ring = int(header[0]), int(header[1]), int(header[2]), int(header[3])

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

#  Figure 1: 3-D mesh
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

#  Figure 2: tube cross-section point distribution
# total_rings = nv_spiral // pts_per_ring
# print(f"  pts_per_ring = {pts_per_ring},  total_rings = {total_rings}")
#
# # Use the middle ring so the Frenet frame is well away from end effects
# mid_ring = total_rings // 2
# ring_pts = vertices[mid_ring * pts_per_ring : (mid_ring + 1) * pts_per_ring]
#
# # Project onto the local plane via SVD
# center   = ring_pts.mean(axis=0)
# centered = ring_pts - center
# _, _, Vt = np.linalg.svd(centered)
# e1, e2   = Vt[0], Vt[1]
# u = centered @ e1
# v = centered @ e2
#
# fig2, ax2 = plt.subplots(figsize=(6, 6))
#
# ax2.scatter(u, v, s=100, c='steelblue', zorder=5)
#
# # Reference circle (mean radius)
# r_ref   = np.hypot(u, v).mean()
# theta_r = np.linspace(0, 2 * np.pi, 300)
# ax2.plot(r_ref * np.cos(theta_r), r_ref * np.sin(theta_r),
#          'k--', lw=0.8, alpha=0.5, label=f'Mean radius = {r_ref:.4f}')
#
# for i, (ui, vi) in enumerate(zip(u, v)):
#     ax2.annotate(str(i + 1), (ui, vi),
#                  textcoords='offset points', xytext=(5, 4), fontsize=8)
#
# ax2.set_aspect('equal')
# ax2.set_xlabel('Local N  (normal)')
# ax2.set_ylabel('Local B  (binormal)')
# ax2.set_title(f'Tube cross-section  —  ring {mid_ring + 1}/{total_rings}'
#               f'\n{pts_per_ring} pts per ring')
# ax2.legend()
# ax2.grid(True, alpha=0.3)
# plt.tight_layout()
#
# #  Figure 3: start cap 2D verification
# # Start-cap faces reference only interior vertices nv_spiral..nv_spiral+2 (0-based),
# # so max vertex index < nv_spiral + 3.
# start_cap_f = cap_f[np.max(cap_f, axis=1) < nv_spiral + 3]
#
# # Local frame from ring 1 (equatorial plane of start cap)
# ring1    = vertices[:pts_per_ring]
# c0       = ring1.mean(axis=0)
# _, _, Vt = np.linalg.svd(ring1 - c0)
# e1, e2   = Vt[0], Vt[1]
#
# # Project every vertex that appears in these faces
# cap_vidx = np.unique(start_cap_f)
# uv = {i: ((vertices[i] - c0) @ e1, (vertices[i] - c0) @ e2) for i in cap_vidx}
#
# # Build polygon list for PolyCollection
# polys = [np.array([uv[v] for v in face]) for face in start_cap_f]
#
# fig3, ax3 = plt.subplots(figsize=(6, 6))
# col = PolyCollection(polys, facecolor='lightsteelblue', edgecolor='black',
#                      linewidth=1.2, alpha=0.75)
# ax3.add_collection(col)
#
# # Ring vs interior vertices
# ring_vidx = cap_vidx[cap_vidx <  nv_spiral]
# int_vidx  = cap_vidx[cap_vidx >= nv_spiral]
#
# ax3.scatter([uv[i][0] for i in ring_vidx], [uv[i][1] for i in ring_vidx],
#             s=80, c='steelblue', zorder=5, label='Ring vertices')
# ax3.scatter([uv[i][0] for i in int_vidx],  [uv[i][1] for i in int_vidx],
#             s=100, c='tomato', marker='^', zorder=5, label='Cap interior pts')
#
# for i in cap_vidx:
#     ax3.annotate(str(i + 1), uv[i], textcoords='offset points',
#                  xytext=(4, 4), fontsize=9)
#
# ax3.set_aspect('equal')
# ax3.autoscale_view()
# ax3.set_xlabel('Local N')
# ax3.set_ylabel('Local B')
# ax3.set_title(f'Start cap  —  2D projection\n{len(start_cap_f)} triangles')
# ax3.legend()
# ax3.grid(True, alpha=0.3)
# plt.tight_layout()
#
# #  Figure 4: start cap 3D inspection
# # Reutilizamos as faces identificadas na Figure 3
# if len(start_cap_f) > 0:
#     fig4 = plt.figure(figsize=(9, 7))
#     ax4 = fig4.add_subplot(111, projection='3d')
#
#     # Adiciona a coleção de triângulos do cap inicial
#     ax4.add_collection3d(make_tri_collection(start_cap_f, 'tomato', 0.9))
#
#     # Identifica apenas os vértices que compõem o cap para ajustar o zoom
#     cap_vidx = np.unique(start_cap_f)
#     cap_coords = vertices[cap_vidx]
#
#     # Ajuste do Bounding Box focado apenas no Cap
#     c_lo, c_hi = cap_coords.min(axis=0), cap_coords.max(axis=0)
#     c_mid = 0.5 * (c_lo + c_hi)
#     c_half = 0.6 * (c_hi - c_lo).max()  # 0.6 para dar uma pequena margem
#
#     ax4.set_xlim(c_mid[0] - c_half, c_mid[0] + c_half)
#     ax4.set_ylim(c_mid[1] - c_half, c_mid[1] + c_half)
#     ax4.set_zlim(c_mid[2] - c_half, c_mid[2] + c_half)
#
#     # Adiciona os números dos vértices para facilitar a conferência com o Fortran
#     for idx in cap_vidx:
#         ax4.text(vertices[idx, 0], vertices[idx, 1], vertices[idx, 2],
#                  str(idx + 1), fontsize=8, color='black')
#
#     ax4.set_title(f'Start Cap 3D Detail\nVertices {nv_spiral+1} a {nv_spiral+3} (interior)')
#     ax4.set_xlabel('X')
#     ax4.set_ylabel('Y')
#     ax4.set_zlabel('Z')
#     plt.tight_layout()

plt.show()
