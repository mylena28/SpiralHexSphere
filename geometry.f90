module geometry
    use params
    implicit none

contains

    ! Project vertex v onto the sphere of radius r centred at c.
    subroutine normalize_to_sphere(v, c, r)
        real(dp), intent(inout) :: v(3)
        real(dp), intent(in)    :: c(3), r
        real(dp) :: d(3), len
        d   = v - c
        len = sqrt(d(1)**2 + d(2)**2 + d(3)**2)
        v   = c + r * d / len
    end subroutine normalize_to_sphere

    subroutine generate_full_mesh(total_nv, total_nf, nv_spiral, verts, faces)
        integer,  intent(out) :: total_nv, total_nf, nv_spiral
        real(dp), allocatable, intent(out) :: verts(:,:)   ! (3, total_nv)
        integer,  allocatable, intent(out) :: faces(:,:)   ! (3, total_nf)

        integer  :: total_rings, pts_per_ring
        integer  :: nv_hemi, nf_spiral, nf_hemi
        integer  :: i, j, k, idx, iface, next_k
        integer  :: p1, p2, p3, p4
        integer  :: pm12, pm13, pm34, pm24, pm_diag
        integer  :: pole_start, pole_end
        integer  :: A, B, mAB, mAC, mBC, C_idx

        ! Frenet frame along helix
        real(dp) :: omega, norm_v
        real(dp) :: s_frac, t_param, rot_off, angle
        real(dp) :: vec_P(3), vec_T(3), vec_N(3), vec_B(3)

        ! saved frames at start and end rings (only P and T needed for cap poles)
        real(dp) :: vec_P_first(3) = 0.0_dp
        real(dp) :: vec_T_first(3) = 0.0_dp
        real(dp) :: vec_P_last(3)  = 0.0_dp
        real(dp) :: vec_T_last(3)  = 0.0_dp

        ! ── Derived dimensions ───────────────────────────────────────────────
        total_rings  = (n_s - 1) * 2 + 1
        pts_per_ring = num_sides * 2

        nv_spiral = total_rings * pts_per_ring
        nf_spiral = (n_s - 1) * num_sides * 8

        ! Hemisphere cap (Zinchenko Level-1 icosahedral subdivision):
        !   new vertices per cap = num_sides pole-edge midpoints + 1 pole
        !   new faces    per cap = 4 sub-triangles × num_sides base triangles
        nv_hemi = num_sides + 1
        nf_hemi = 4 * num_sides

        total_nv = nv_spiral + 2 * nv_hemi
        total_nf = nf_spiral + 2 * nf_hemi

        allocate(verts(3, total_nv))
        allocate(faces(3, total_nf))

        omega  = 2.0_dp * PI * n_turns / helix_h
        norm_v = sqrt((helix_R * omega)**2 + 1.0_dp)

        ! ════════════════════════════════════════════════════════════════════
        ! PART 1 — SPIRAL VERTICES
        ! ════════════════════════════════════════════════════════════════════
        do i = 1, total_rings
            s_frac  = real(i - 1, dp) / real(total_rings - 1, dp)
            t_param = 2.0_dp * PI * n_turns * s_frac
            rot_off = real(i - 1, dp) * (PI / real(pts_per_ring, dp))

            vec_P = [ helix_R * cos(t_param),  &
                      helix_R * sin(t_param),  &
                      helix_h * s_frac          ]

            vec_N = [ -cos(t_param), -sin(t_param), 0.0_dp ]

            vec_T(1) = -helix_R * omega * sin(t_param) / norm_v
            vec_T(2) =  helix_R * omega * cos(t_param) / norm_v
            vec_T(3) =  1.0_dp / norm_v

            ! Binormal  B = T × N
            vec_B(1) = vec_T(2)*vec_N(3) - vec_T(3)*vec_N(2)
            vec_B(2) = vec_T(3)*vec_N(1) - vec_T(1)*vec_N(3)
            vec_B(3) = vec_T(1)*vec_N(2) - vec_T(2)*vec_N(1)

            do j = 1, pts_per_ring
                angle = 2.0_dp * PI * real(j - 1, dp) / real(pts_per_ring, dp) + rot_off
                idx = (i - 1) * pts_per_ring + j
                verts(1, idx) = vec_P(1) + tube_r*(vec_N(1)*cos(angle) + vec_B(1)*sin(angle))
                verts(2, idx) = vec_P(2) + tube_r*(vec_N(2)*cos(angle) + vec_B(2)*sin(angle))
                verts(3, idx) = vec_P(3) + tube_r*(vec_N(3)*cos(angle) + vec_B(3)*sin(angle))
            end do

            if (i == 1) then
                vec_P_first = vec_P;  vec_T_first = vec_T
            end if
            if (i == total_rings) then
                vec_P_last = vec_P;  vec_T_last = vec_T
            end if
        end do

        ! ════════════════════════════════════════════════════════════════════
        ! PART 2 — HEMISPHERE CAP VERTICES (Zinchenko Level-1 subdivision)
        !
        ! Coarse base mesh for each cap: a fan pyramid
        !   - num_sides "original" equatorial vertices (odd-j positions in end ring)
        !   - 1 pole vertex
        !
        ! After Level-1 midpoint subdivision each base triangle → 4:
        !   mAB  equatorial edge midpoint  → already stored at even-j in end ring
        !   mAC  pole-edge midpoint        → NEW, normalized to sphere
        !   mBC  pole-edge midpoint        → NEW, normalized to sphere
        !
        ! In start ring (i=1):
        !   eq_orig(k) → global index  2*k-1         (k = 1..num_sides)
        !   eq_mid(k)  → global index  2*k
        ! In end ring (i=total_rings):
        !   eq_orig(k) → (total_rings-1)*pts_per_ring + 2*k-1
        !   eq_mid(k)  → (total_rings-1)*pts_per_ring + 2*k
        !
        ! New vertex layout:
        !   nv_spiral + k               → start-cap pole-edge mid k  (k=1..num_sides)
        !   nv_spiral + num_sides + 1   → start-cap pole
        !   nv_spiral + nv_hemi + k     → end-cap   pole-edge mid k
        !   nv_spiral + 2*nv_hemi       → end-cap   pole
        ! ────────────────────────────────────────────────────────────────────

        ! ── START CAP ──
        pole_start = nv_spiral + num_sides + 1
        verts(:, pole_start) = vec_P_first - tube_r * vec_T_first

        do k = 1, num_sides
            idx = nv_spiral + k
            j   = 2*k - 1                          ! eq_orig(k) in ring 1
            verts(:, idx) = 0.5_dp * (verts(:, j) + verts(:, pole_start))
            call normalize_to_sphere(verts(:, idx), vec_P_first, tube_r)
        end do

        ! ── END CAP ──
        pole_end = nv_spiral + 2 * nv_hemi
        verts(:, pole_end) = vec_P_last + tube_r * vec_T_last

        do k = 1, num_sides
            idx = nv_spiral + nv_hemi + k
            j   = (total_rings - 1)*pts_per_ring + 2*k - 1   ! eq_orig(k) in last ring
            verts(:, idx) = 0.5_dp * (verts(:, j) + verts(:, pole_end))
            call normalize_to_sphere(verts(:, idx), vec_P_last, tube_r)
        end do

        ! ════════════════════════════════════════════════════════════════════
        ! PART 3 — SPIRAL FACES (Zinchenko Level-1 subdivision, unchanged)
        ! ════════════════════════════════════════════════════════════════════
        !
        !   Layout (left=col j, centre=col j+1, right=col j+2):
        !
        !     p3 ──pm34── p4
        !     │    │    │
        !    pm13 pm_diag pm24
        !     │    │    │
        !     p1 ──pm12── p2
        !
        iface = 0
        do i = 1, total_rings - 2, 2
            do j = 1, pts_per_ring, 2

                p1   = (i - 1)*pts_per_ring + j
                pm12 = (i - 1)*pts_per_ring + j + 1
                p2   = (i - 1)*pts_per_ring + mod(j + 1, pts_per_ring) + 1

                p3   = (i + 1)*pts_per_ring + j
                pm34 = (i + 1)*pts_per_ring + j + 1
                p4   = (i + 1)*pts_per_ring + mod(j + 1, pts_per_ring) + 1

                pm13    = i*pts_per_ring + j
                pm_diag = i*pts_per_ring + j + 1
                pm24    = i*pts_per_ring + mod(j + 1, pts_per_ring) + 1

                ! Left half of cell
                iface = iface + 1;  faces(:, iface) = [p1,      pm12,    pm13   ]
                iface = iface + 1;  faces(:, iface) = [pm12,    p2,      pm_diag]
                iface = iface + 1;  faces(:, iface) = [pm13,    pm_diag, p3     ]
                iface = iface + 1;  faces(:, iface) = [pm12,    pm_diag, pm13   ]

                ! Right half of cell
                iface = iface + 1;  faces(:, iface) = [p2,      pm24,    pm_diag]
                iface = iface + 1;  faces(:, iface) = [pm24,    p4,      pm34   ]
                iface = iface + 1;  faces(:, iface) = [pm_diag, pm34,    p3     ]
                iface = iface + 1;  faces(:, iface) = [pm24,    pm34,    pm_diag]

            end do
        end do

        ! ════════════════════════════════════════════════════════════════════
        ! PART 4 — HEMISPHERE CAP FACES (Zinchenko Level-1 subdivision)
        !
        ! For the k-th base triangle the vertices are:
        !   A, B  — two adjacent original equatorial vertices
        !   C     — pole
        ! and midpoints:
        !   mAB   — equatorial edge midpoint  (already in spiral end ring)
        !   mAC   — pole-edge midpoint from A (new vertex)
        !   mBC   — pole-edge midpoint from B (new vertex)
        !
        ! Standard CCW midpoint subdivision of triangle (A, B, C) → 4 sub-triangles:
        !   (A,   mAB, mAC)   corner A
        !   (mAB, B,   mBC)   corner B
        !   (mAC, mBC, C  )   corner C
        !   (mAB, mBC, mAC)   centre
        !
        ! START CAP winding:
        !   Equatorial ring viewed from outside (-T) appears CW.
        !   Base triangle is (A=eq_orig(k+1), B=eq_orig(k), C=pole) — reversed —
        !   so it is CCW from outside, giving outward normals in -T direction.
        !   Then: mAB=eq_mid(k), mAC=pole_edge_mid(k+1), mBC=pole_edge_mid(k).
        !
        ! END CAP winding:
        !   Equatorial ring viewed from outside (+T) appears CCW.
        !   Base triangle is (A=eq_orig(k), B=eq_orig(k+1), C=pole) — CCW from outside.
        !   Then: mAB=eq_mid(k), mAC=pole_edge_mid(k), mBC=pole_edge_mid(k+1).
        ! ════════════════════════════════════════════════════════════════════

        ! ── START CAP FACES ──
        do k = 1, num_sides
            next_k = mod(k, num_sides) + 1
            A     = 2*next_k - 1                    ! eq_orig(k+1) in ring 1
            B     = 2*k - 1                         ! eq_orig(k)   in ring 1
            mAB   = 2*k                             ! eq_mid(k)    in ring 1
            C_idx = pole_start
            mAC   = nv_spiral + next_k              ! pole-edge mid for eq_orig(k+1)
            mBC   = nv_spiral + k                   ! pole-edge mid for eq_orig(k)

            iface = iface + 1;  faces(:, iface) = [A,   mAB, mAC  ]
            iface = iface + 1;  faces(:, iface) = [mAB, B,   mBC  ]
            iface = iface + 1;  faces(:, iface) = [mAC, mBC, C_idx]
            iface = iface + 1;  faces(:, iface) = [mAB, mBC, mAC  ]
        end do

        ! ── END CAP FACES ──
        do k = 1, num_sides
            next_k = mod(k, num_sides) + 1
            A     = (total_rings-1)*pts_per_ring + 2*k - 1       ! eq_orig(k)   in last ring
            B     = (total_rings-1)*pts_per_ring + 2*next_k - 1  ! eq_orig(k+1) in last ring
            mAB   = (total_rings-1)*pts_per_ring + 2*k           ! eq_mid(k)    in last ring
            C_idx = pole_end
            mAC   = nv_spiral + nv_hemi + k
            mBC   = nv_spiral + nv_hemi + next_k

            iface = iface + 1;  faces(:, iface) = [A,   mAB, mAC  ]
            iface = iface + 1;  faces(:, iface) = [mAB, B,   mBC  ]
            iface = iface + 1;  faces(:, iface) = [mAC, mBC, C_idx]
            iface = iface + 1;  faces(:, iface) = [mAB, mBC, mAC  ]
        end do

    end subroutine generate_full_mesh

end module geometry
