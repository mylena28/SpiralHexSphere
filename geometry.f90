module geometry
    use params
    implicit none

contains

    subroutine generate_full_mesh(total_nv, total_nf, nv_spiral, verts, faces, s_frac, cap_flag)
        integer,  intent(out) :: total_nv, total_nf, nv_spiral
        real(dp), allocatable, intent(out) :: verts(:,:)   ! (3, total_nv)
        integer,  allocatable, intent(out) :: faces(:,:)   ! (3, total_nf)
        real(dp), allocatable, intent(out) :: s_frac(:)    ! helix arc parameter per vertex (0..1)
        integer,  allocatable, intent(out) :: cap_flag(:)  ! 0=tube, -1=start cap, +1=end cap

        integer  :: total_rings, pts_per_ring
        integer  :: nv_hemi, nf_spiral, nf_hemi
        integer  :: i, j, k, idx, iface
        integer  :: p1, p2, p3, p4
        integer  :: cs1, cs2, cs3   ! start-cap interior pts (indices)
        integer  :: ce1, ce2, ce3   ! end-cap   interior pts (indices)
        integer  :: Sb              ! base index of last ring

        ! Frenet frame
        real(dp) :: omega, norm_v
        real(dp) :: s_param, t_param, rot_off, angle
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
        allocate(s_frac(total_nv))
        allocate(cap_flag(total_nv))

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
            s_param = real(i - 1, dp) / real(total_rings - 1, dp)
            t_param = 2.0_dp * PI * n_turns * s_param
            rot_off = real(i - 1, dp) * (PI / real(pts_per_ring, dp))

            vec_P = [ helix_R * cos(t_param), helix_R * sin(t_param), helix_h * s_param ]
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
                s_frac(idx)   = s_param
                cap_flag(idx) = 0
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
            verts(:, idx)    = vec_P_first + tube_r * (sinp * d_vec - cosp * vec_T_first)
            s_frac(idx)      = 0.0_dp
            cap_flag(idx)    = -1
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
            verts(:, idx)    = vec_P_last + tube_r * (sinp * d_vec + cosp * vec_T_last)
            s_frac(idx)      = 1.0_dp
            cap_flag(idx)    = 1
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

    ! ════════════════════════════════════════════════════════════════════════
    ! SUBDIVIDE_MESH — 4-1 midpoint subdivision (one level)
    !
    ! Each triangle is split into 4 by inserting the midpoint of every edge:
    !
    !        v1                v1
    !       /  \              / \
    !      /    \    →    m31 --- m12
    !     /      \          / \ / \
    !    v3------v2        v3-m23--v2
    !
    ! Intra-spiral edges are inserted first so new spiral midpoints keep
    ! consecutive indices, preserving the nv_spiral boundary for plot_mesh.py.
    ! ════════════════════════════════════════════════════════════════════════
    subroutine subdivide_mesh(nv_in, nf_in, nv_sp_in, verts_in, faces_in, &
                              s_frac_in, cap_flag_in, &
                              nv_out, nf_out, nv_sp_out, verts_out, faces_out, &
                              s_frac_out, cap_flag_out)
        integer,  intent(in)  :: nv_in, nf_in, nv_sp_in
        real(dp), intent(in)  :: verts_in(3, nv_in)
        integer,  intent(in)  :: faces_in(3, nf_in)
        real(dp), intent(in)  :: s_frac_in(nv_in)
        integer,  intent(in)  :: cap_flag_in(nv_in)
        integer,  intent(out) :: nv_out, nf_out, nv_sp_out
        real(dp), allocatable, intent(out) :: verts_out(:,:)
        integer,  allocatable, intent(out) :: faces_out(:,:)
        real(dp), allocatable, intent(out) :: s_frac_out(:)
        integer,  allocatable, intent(out) :: cap_flag_out(:)

        integer :: max_edges, n_edges, n_sp_edges, e, f, iout
        integer :: v1, v2, v3, m12, m23, m31
        integer, allocatable :: ev1(:), ev2(:), emid(:)

        ! Surface projection variables
        real(dp) :: omega, norm_v
        real(dp) :: chord_mid(3), s_mid, d_vec(3), d_NB(3), d_norm
        real(dp) :: t_param, vec_P(3), vec_T(3), vec_N(3), vec_B(3)
        real(dp) :: vec_P_first(3), vec_P_last(3)
        integer  :: v1e, v2e, flag1, flag2

        max_edges = 3 * nf_in
        allocate(ev1(max_edges), ev2(max_edges), emid(max_edges))
        n_edges = 0

        ! Phase 1: intra-spiral edges (both endpoints inside spiral region)
        do f = 1, nf_in
            v1 = faces_in(1,f);  v2 = faces_in(2,f);  v3 = faces_in(3,f)
            if (v1 <= nv_sp_in .and. v2 <= nv_sp_in) &
                call add_edge(ev1, ev2, emid, n_edges, nv_in, v1, v2)
            if (v2 <= nv_sp_in .and. v3 <= nv_sp_in) &
                call add_edge(ev1, ev2, emid, n_edges, nv_in, v2, v3)
            if (v3 <= nv_sp_in .and. v1 <= nv_sp_in) &
                call add_edge(ev1, ev2, emid, n_edges, nv_in, v3, v1)
        end do
        n_sp_edges = n_edges   ! midpoints of these edges stay in the spiral region

        ! Phase 2: boundary and intra-cap edges
        do f = 1, nf_in
            v1 = faces_in(1,f);  v2 = faces_in(2,f);  v3 = faces_in(3,f)
            if (.not. (v1 <= nv_sp_in .and. v2 <= nv_sp_in)) &
                call add_edge(ev1, ev2, emid, n_edges, nv_in, v1, v2)
            if (.not. (v2 <= nv_sp_in .and. v3 <= nv_sp_in)) &
                call add_edge(ev1, ev2, emid, n_edges, nv_in, v2, v3)
            if (.not. (v3 <= nv_sp_in .and. v1 <= nv_sp_in)) &
                call add_edge(ev1, ev2, emid, n_edges, nv_in, v3, v1)
        end do

        nv_sp_out = nv_sp_in + n_sp_edges
        nv_out    = nv_in    + n_edges
        nf_out    = 4 * nf_in

        allocate(verts_out(3, nv_out))
        allocate(faces_out(3, nf_out))
        allocate(s_frac_out(nv_out))
        allocate(cap_flag_out(nv_out))

        ! Keep all original vertices in place
        verts_out(:, 1:nv_in) = verts_in
        s_frac_out(1:nv_in)   = s_frac_in
        cap_flag_out(1:nv_in) = cap_flag_in

        ! Helix constants for tube circle projection
        omega  = 2.0_dp * PI * n_turns / helix_h
        norm_v = sqrt((helix_R * omega)**2 + 1.0_dp)
        vec_P_first = [helix_R, 0.0_dp, 0.0_dp]
        vec_P_last  = [helix_R*cos(2.0_dp*PI*n_turns), &
                       helix_R*sin(2.0_dp*PI*n_turns), helix_h]

        ! New midpoint vertices: projected onto the correct curved surface
        do e = 1, n_edges
            v1e  = ev1(e);  v2e  = ev2(e)
            s_mid     = 0.5_dp * (s_frac_in(v1e) + s_frac_in(v2e))
            chord_mid = 0.5_dp * (verts_in(:, v1e) + verts_in(:, v2e))
            flag1 = cap_flag_in(v1e);  flag2 = cap_flag_in(v2e)
            s_frac_out(nv_in + e) = s_mid

            if (flag1 == 1 .or. flag2 == 1) then
                ! End cap: project onto sphere of radius tube_r centred at vec_P_last
                cap_flag_out(nv_in + e) = 1
                d_vec  = chord_mid - vec_P_last
                d_norm = sqrt(d_vec(1)**2 + d_vec(2)**2 + d_vec(3)**2)
                verts_out(:, nv_in+e) = vec_P_last + tube_r * d_vec / d_norm

            else if (flag1 == -1 .or. flag2 == -1) then
                ! Start cap: project onto sphere of radius tube_r centred at vec_P_first
                cap_flag_out(nv_in + e) = -1
                d_vec  = chord_mid - vec_P_first
                d_norm = sqrt(d_vec(1)**2 + d_vec(2)**2 + d_vec(3)**2)
                verts_out(:, nv_in+e) = vec_P_first + tube_r * d_vec / d_norm

            else
                ! Tube body: project onto the cross-section circle at s_mid.
                ! Remove the tangent component of (chord_mid - vec_P), then
                ! normalise in the N-B plane and scale by tube_r.
                cap_flag_out(nv_in + e) = 0
                t_param = 2.0_dp * PI * n_turns * s_mid
                vec_P = [helix_R*cos(t_param), helix_R*sin(t_param), helix_h*s_mid]
                vec_N = [-cos(t_param), -sin(t_param), 0.0_dp]
                vec_T(1) = -helix_R * omega * sin(t_param) / norm_v
                vec_T(2) =  helix_R * omega * cos(t_param) / norm_v
                vec_T(3) =  1.0_dp / norm_v
                vec_B(1) = vec_T(2)*vec_N(3) - vec_T(3)*vec_N(2)
                vec_B(2) = vec_T(3)*vec_N(1) - vec_T(1)*vec_N(3)
                vec_B(3) = vec_T(1)*vec_N(2) - vec_T(2)*vec_N(1)
                d_vec = chord_mid - vec_P
                d_norm = d_vec(1)*vec_T(1) + d_vec(2)*vec_T(2) + d_vec(3)*vec_T(3)
                d_NB(1) = d_vec(1) - d_norm * vec_T(1)
                d_NB(2) = d_vec(2) - d_norm * vec_T(2)
                d_NB(3) = d_vec(3) - d_norm * vec_T(3)
                d_norm = sqrt(d_NB(1)**2 + d_NB(2)**2 + d_NB(3)**2)
                verts_out(:, nv_in+e) = vec_P + tube_r * d_NB / d_norm
            end if
        end do

        ! Build 4 child triangles per parent
        iout = 0
        do f = 1, nf_in
            v1 = faces_in(1,f);  v2 = faces_in(2,f);  v3 = faces_in(3,f)
            call find_edge(ev1, ev2, emid, n_edges, v1, v2, m12)
            call find_edge(ev1, ev2, emid, n_edges, v2, v3, m23)
            call find_edge(ev1, ev2, emid, n_edges, v3, v1, m31)
            iout = iout+1;  faces_out(:,iout) = [v1,  m12, m31]
            iout = iout+1;  faces_out(:,iout) = [m12, v2,  m23]
            iout = iout+1;  faces_out(:,iout) = [m31, m23, v3 ]
            iout = iout+1;  faces_out(:,iout) = [m12, m23, m31]
        end do

        deallocate(ev1, ev2, emid)

    end subroutine subdivide_mesh

    ! ════════════════════════════════════════════════════════════════════════
    ! CHECK_PROJECTION_ERROR
    !
    ! For every vertex, computes the distance from the vertex to its target
    ! surface and returns the maximum |distance - tube_r|:
    !   cap_flag =  0  →  NB-plane radius at s_frac should equal tube_r
    !   cap_flag = -1  →  3-D distance to vec_P_first should equal tube_r
    !   cap_flag = +1  →  3-D distance to vec_P_last  should equal tube_r
    ! ════════════════════════════════════════════════════════════════════════
    subroutine check_projection_error(nv, verts, s_frac, cap_flag, max_err, max_idx)
        integer,  intent(in)  :: nv
        real(dp), intent(in)  :: verts(3, nv)
        real(dp), intent(in)  :: s_frac(nv)
        integer,  intent(in)  :: cap_flag(nv)
        real(dp), intent(out) :: max_err
        integer,  intent(out) :: max_idx

        integer  :: i
        real(dp) :: err, t_param, d_vec(3), d_dot
        real(dp) :: vec_P(3), vec_T(3), vec_P_first(3), vec_P_last(3)
        real(dp) :: omega, norm_v

        omega  = 2.0_dp * PI * n_turns / helix_h
        norm_v = sqrt((helix_R * omega)**2 + 1.0_dp)
        vec_P_first = [helix_R, 0.0_dp, 0.0_dp]
        vec_P_last  = [helix_R*cos(2.0_dp*PI*n_turns), &
                       helix_R*sin(2.0_dp*PI*n_turns), helix_h]

        max_err = 0.0_dp
        max_idx = 0

        do i = 1, nv
            select case (cap_flag(i))
            case (-1)
                d_vec = verts(:,i) - vec_P_first
                err   = abs(sqrt(d_vec(1)**2 + d_vec(2)**2 + d_vec(3)**2) - tube_r)
            case (1)
                d_vec = verts(:,i) - vec_P_last
                err   = abs(sqrt(d_vec(1)**2 + d_vec(2)**2 + d_vec(3)**2) - tube_r)
            case default
                t_param  = 2.0_dp * PI * n_turns * s_frac(i)
                vec_P    = [helix_R*cos(t_param), helix_R*sin(t_param), helix_h*s_frac(i)]
                vec_T(1) = -helix_R * omega * sin(t_param) / norm_v
                vec_T(2) =  helix_R * omega * cos(t_param) / norm_v
                vec_T(3) =  1.0_dp / norm_v
                d_vec    = verts(:,i) - vec_P
                d_dot    = d_vec(1)*vec_T(1) + d_vec(2)*vec_T(2) + d_vec(3)*vec_T(3)
                err = abs(sqrt(max(0.0_dp, &
                    d_vec(1)**2+d_vec(2)**2+d_vec(3)**2 - d_dot**2)) - tube_r)
            end select
            if (err > max_err) then
                max_err = err
                max_idx = i
            end if
        end do

    end subroutine check_projection_error

    ! ── Private helpers ──────────────────────────────────────────────────────

    subroutine add_edge(ev1, ev2, emid, n, nv_base, a, b)
        integer, intent(inout) :: ev1(:), ev2(:), emid(:), n
        integer, intent(in)    :: nv_base, a, b
        integer :: lo, hi, e
        lo = min(a,b);  hi = max(a,b)
        do e = 1, n
            if (ev1(e)==lo .and. ev2(e)==hi) return
        end do
        n = n + 1
        ev1(n) = lo;  ev2(n) = hi;  emid(n) = nv_base + n
    end subroutine add_edge

    subroutine find_edge(ev1, ev2, emid, n, a, b, mid_idx)
        integer, intent(in)  :: ev1(:), ev2(:), emid(:), n, a, b
        integer, intent(out) :: mid_idx
        integer :: lo, hi, e
        lo = min(a,b);  hi = max(a,b)
        do e = 1, n
            if (ev1(e)==lo .and. ev2(e)==hi) then
                mid_idx = emid(e);  return
            end if
        end do
        mid_idx = -1   ! unreachable if mesh is consistent
    end subroutine find_edge

end module geometry
