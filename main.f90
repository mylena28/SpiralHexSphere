program main
    use params
    use geometry
    use io
    implicit none

    integer  :: total_nv, total_nf, nv_spiral
    integer  :: total_rings, pts_per_ring
    integer  :: lev, nv_new, nf_new, nv_sp_new
    real(dp), allocatable :: verts(:,:), verts_new(:,:)
    integer,  allocatable :: faces(:,:), faces_new(:,:)
    real(dp), allocatable :: s_frac(:), s_frac_new(:)
    integer,  allocatable :: cap_flag(:), cap_flag_new(:)
    real(dp), allocatable :: normals(:,:)

    n_s = nint(helix_R * n_turns * real(num_sides, dp) / tube_r)

    total_rings  = n_s
    pts_per_ring = num_sides

    call generate_full_mesh(total_nv, total_nf, nv_spiral, verts, faces, s_frac, cap_flag)

    ! Apply 4-1 midpoint subdivision refine_level times
    do lev = 1, refine_level
        call subdivide_mesh(total_nv, total_nf, nv_spiral, verts, faces, &
                            s_frac, cap_flag, &
                            nv_new, nf_new, nv_sp_new, verts_new, faces_new, &
                            s_frac_new, cap_flag_new)
        call move_alloc(verts_new, verts)
        call move_alloc(faces_new, faces)
        call move_alloc(s_frac_new, s_frac)
        call move_alloc(cap_flag_new, cap_flag)
        total_nv  = nv_new
        total_nf  = nf_new
        nv_spiral = nv_sp_new
    end do

    allocate(normals(3, total_nv))
    call compute_vertex_normals(total_nv, total_nf, verts, faces, normals)
    call write_mesh(trim(outfile), total_nv, total_nf,verts, faces, normals)

    deallocate(verts, faces, s_frac, cap_flag, normals)

end program main
