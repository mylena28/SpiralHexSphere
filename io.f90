! io.f90 — I/O routines for SpiralHexSphere
module io
    use params
    implicit none

contains

    ! ──────────────────────────────────────────────────────────────────────────
    ! Write the combined mesh in the simple binary-free text format:
    !
    !   Line 1 :  nv  nf
    !   Lines 2…nv+1 :  x  y  z   (one vertex per line)
    !   Lines nv+2… :  i  j  k   (1-based face indices, one face per line)
    !
    ! This format is identical to what Spiral/io.f90 and Sphere/write_mesh.f
    ! produce, so plot_mesh.py can read it directly.
    ! ──────────────────────────────────────────────────────────────────────────
    subroutine write_mesh(filename, total_nv, total_nf, nv_spiral, pts_per_ring, verts, faces)
        character(len=*), intent(in) :: filename
        integer,  intent(in) :: total_nv, total_nf, nv_spiral, pts_per_ring
        real(dp), intent(in) :: verts(3, total_nv)
        integer,  intent(in) :: faces(3, total_nf)

        integer :: iunit, i

        open(newunit=iunit, file=filename, status='replace', action='write')

        ! Header: vertex count, face count, spiral vertex count, pts per ring
        write(iunit, '(4(I10,1X))') total_nv, total_nf, nv_spiral, pts_per_ring

        ! Vertices
        do i = 1, total_nv
            write(iunit, '(3(F14.8,1X))') verts(1,i), verts(2,i), verts(3,i)
        end do

        ! Faces (1-based indices)
        do i = 1, total_nf
            write(iunit, '(3(I10,1X))') faces(1,i), faces(2,i), faces(3,i)
        end do

        close(iunit)
        write(*,'(A,A)') '  Mesh written to: ', trim(filename)

    end subroutine write_mesh

end module io
