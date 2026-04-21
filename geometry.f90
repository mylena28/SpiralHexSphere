module geometry
    use params
    implicit none

contains

    subroutine generate_full_mesh(total_nv, total_nf, nv_spiral, verts, faces)
        integer,  intent(out) :: total_nv, total_nf, nv_spiral
        real(dp), allocatable, intent(out) :: verts(:,:)   ! (3, total_nv)
        integer,  allocatable, intent(out) :: faces(:,:)   ! (3, total_nf)

        integer  :: total_rings, pts_per_ring
        integer  :: nv_hemi, nf_spiral, nf_hemi
        integer  :: i, j, k, idx, iface
        integer  :: p1, p2, p3, p4
        integer  :: cs1, cs2, cs3   ! start-cap interior pts (indices)
        integer  :: ce1, ce2, ce3   ! end-cap   interior pts (indices)
        integer  :: Sb              ! base index of last ring

        ! Frenet frame
        real(dp) :: omega, norm_v
        real(dp) :: s_frac, t_param, rot_off, angle
        real(dp) :: vec_P(3), vec_T(3), vec_N(3), vec_B(3)
        real(dp) :: d_vec(3)

        real(dp) :: vec_P_first(3) = 0.0_dp, vec_T_first(3) = 0.0_dp
        real(dp) :: vec_N_first(3) = 0.0_dp, vec_B_first(3) = 0.0_dp
        real(dp) :: vec_P_last(3)  = 0.0_dp, vec_T_last(3)  = 0.0_dp
        real(dp) :: vec_N_last(3)  = 0.0_dp, vec_B_last(3)  = 0.0_dp
        real(dp) :: sinp = 0.0_dp, cosp = 0.0_dp

        ! Level-0: one ring per original section, one vertex per polygon side
        total_rings  = n_s
        pts_per_ring = num_sides

        nv_spiral = total_rings * pts_per_ring
        nf_spiral = (n_s - 1) * num_sides * 2

        ! Cap: 3 interior pts + 0 pole (works for num_sides = 6)
        !   outer triangles : num_sides      = 6
        !   bridge triangles: num_sides/2    = 3
        !   center triangle :               = 1
        !   total per cap                   = 10
        nv_hemi = 3
        nf_hemi = 10

        total_nv = nv_spiral + 2 * nv_hemi
        total_nf = nf_spiral + 2 * nf_hemi

        allocate(verts(3, total_nv))
        allocate(faces(3, total_nf))

        ! Optimal elevation for minimum edge-length range on the cap hemisphere.
        ! Derived from L2 = L4  =>  sin(phi) = (sqrt(7)-1)/3  ≈ 0.5486
        ! This gives range(L2,L3,L4) ≈ 0.255  vs 0.460 at 45°.
        sinp = (sqrt(7.0_dp) - 1.0_dp) / 3.0_dp
        cosp = sqrt(1.0_dp - sinp**2)
        omega  = 2.0_dp * PI * n_turns / helix_h
        norm_v = sqrt((helix_R * omega)**2 + 1.0_dp)

        ! ════════════════════════════════════════════════════════════════════
        ! PART 1 — SPIRAL VERTICES
        ! ════════════════════════════════════════════════════════════════════
        do i = 1, total_rings
            s_frac  = real(i - 1, dp) / real(total_rings - 1, dp)
            t_param = 2.0_dp * PI * n_turns * s_frac
            rot_off = real(i - 1, dp) * (PI / real(pts_per_ring, dp))

            vec_P = [ helix_R * cos(t_param), helix_R * sin(t_param), helix_h * s_frac ]
            vec_N = [ -cos(t_param), -sin(t_param), 0.0_dp ]
            vec_T(1) = -helix_R * omega * sin(t_param) / norm_v
            vec_T(2) =  helix_R * omega * cos(t_param) / norm_v
            vec_T(3) =  1.0_dp / norm_v
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
                vec_N_first = vec_N;  vec_B_first = vec_B
            end if
            if (i == total_rings) then
                vec_P_last = vec_P;  vec_T_last = vec_T
                vec_N_last = vec_N;  vec_B_last = vec_B
            end if
        end do

        ! ════════════════════════════════════════════════════════════════════
        ! PART 2 — CAP INTERIOR VERTICES (Optimized Spherical Projection)
        ! ════════════════════════════════════════════════════════════════════

        ! --- START CAP (Elevated toward -T_first) ---
        cs1 = nv_spiral + 1
        cs2 = nv_spiral + 2
        cs3 = nv_spiral + 3
        do k = 0, 2
            j     = 2*k + 1                                 ! Azimuths: 1, 3, 5
            angle = 2.0_dp * PI * real(j-1, dp) / real(num_sides, dp)
            d_vec = vec_N_first * cos(angle) + vec_B_first * sin(angle)
            idx   = nv_spiral + k + 1
            ! Projeção: Componente Radial (sinp) e Componente Axial (cosp)
            verts(:, idx) = vec_P_first + tube_r * (sinp * d_vec - cosp * vec_T_first)
        end do

        ! --- END CAP (Elevated toward +T_last) ---
        ! Sincroniza o rot_off para que os 3 pontos internos acompanhem o giro do último anel
        rot_off = real(total_rings - 1, dp) * (PI / real(pts_per_ring, dp))
        ce1 = nv_spiral + 4
        ce2 = nv_spiral + 5
        ce3 = nv_spiral + 6
        do k = 0, 2
            j     = 2*k + 1
            angle = 2.0_dp * PI * real(j-1, dp) / real(num_sides, dp) + rot_off
            d_vec = vec_N_last * cos(angle) + vec_B_last * sin(angle)
            idx   = nv_spiral + 3 + k + 1
            ! Agora usa sinp e cosp corretamente para evitar o achatamento
            verts(:, idx) = vec_P_last + tube_r * (sinp * d_vec + cosp * vec_T_last)
        end do

        ! ════════════════════════════════════════════════════════════════════
        ! PART 3 — SPIRAL FACES  (quad → 2 triangles per column)
        ! ════════════════════════════════════════════════════════════════════
        iface = 0
        do i = 1, total_rings - 1
            do j = 1, pts_per_ring
                p1 = (i - 1)*pts_per_ring + j
                p2 = (i - 1)*pts_per_ring + mod(j, pts_per_ring) + 1
                p3 = i*pts_per_ring + j
                p4 = i*pts_per_ring + mod(j, pts_per_ring) + 1
                iface = iface + 1;  faces(:, iface) = [p1, p2, p3]
                iface = iface + 1;  faces(:, iface) = [p2, p4, p3]
            end do
        end do

        ! ════════════════════════════════════════════════════════════════════
        ! PART 4 — CAP FACES  (10 triangles each, for num_sides = 6)
        !
        ! Start cap ring  : R1..R6  (vertices 1..6  of ring 1)
        !   Q1 = cs1 (azimuths of R1), Q2 = cs2 (R3), Q3 = cs3 (R5)
        !
        ! Connectivity (CCW from outside, i.e., from  -T direction):
        !   outer  : [R6,R1,Q1] [R1,R2,Q1] [R2,R3,Q2] [R3,R4,Q2] [R4,R5,Q3] [R5,R6,Q3]
        !   bridge : [R2,Q2,Q1] [R4,Q3,Q2] [R6,Q1,Q3]
        !   center : [Q1,Q2,Q3]
        ! ════════════════════════════════════════════════════════════════════

        ! ── START CAP ──
        iface = iface + 1;  faces(:, iface) = [6,   1,   cs1]
        iface = iface + 1;  faces(:, iface) = [1,   2,   cs1]
        iface = iface + 1;  faces(:, iface) = [2,   3,   cs2]
        iface = iface + 1;  faces(:, iface) = [3,   4,   cs2]
        iface = iface + 1;  faces(:, iface) = [4,   5,   cs3]
        iface = iface + 1;  faces(:, iface) = [5,   6,   cs3]
        iface = iface + 1;  faces(:, iface) = [2,   cs2, cs1]
        iface = iface + 1;  faces(:, iface) = [4,   cs3, cs2]
        iface = iface + 1;  faces(:, iface) = [6,   cs1, cs3]
        iface = iface + 1;  faces(:, iface) = [cs1, cs2, cs3]

        ! ── END CAP  (same pattern, reversed vertex order for outward +T normal) ──
        Sb = (total_rings - 1) * num_sides   ! offset of last ring
        iface = iface + 1;  faces(:, iface) = [ce1, Sb+1, Sb+6]
        iface = iface + 1;  faces(:, iface) = [ce1, Sb+2, Sb+1]
        iface = iface + 1;  faces(:, iface) = [ce2, Sb+3, Sb+2]
        iface = iface + 1;  faces(:, iface) = [ce2, Sb+4, Sb+3]
        iface = iface + 1;  faces(:, iface) = [ce3, Sb+5, Sb+4]
        iface = iface + 1;  faces(:, iface) = [ce3, Sb+6, Sb+5]
        iface = iface + 1;  faces(:, iface) = [ce1, ce2,  Sb+2]
        iface = iface + 1;  faces(:, iface) = [ce2, ce3,  Sb+4]
        iface = iface + 1;  faces(:, iface) = [ce3, ce1,  Sb+6]
        iface = iface + 1;  faces(:, iface) = [ce3, ce2,  ce1 ]

    end subroutine generate_full_mesh

end module geometry
