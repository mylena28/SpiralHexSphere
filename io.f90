module io
    use params
    implicit none

contains

    subroutine write_mesh(filename, total_nv, total_nf,verts, faces, normals)
        character(len=*), intent(in) :: filename
        integer,  intent(in) :: total_nv, total_nf
        real(dp), intent(in) :: verts(3, total_nv)
        integer,  intent(in) :: faces(3, total_nf)
        real(dp), intent(in) :: normals(3, total_nv)

        integer :: iunit, i

        open(newunit=iunit, file=filename, status='replace', action='write')

        do i = 1, total_nv
            write(iunit, '(A,3(1X,ES18.9))') 'v', verts(1,i), verts(2,i), verts(3,i)
        end do

        do i = 1, total_nv
            write(iunit, '(A,3(1X,ES18.9))') 'vn', normals(1,i), normals(2,i), normals(3,i)
        end do

        do i = 1, total_nf
            write(iunit, '(A,3(1X,I0,A,I0))') 'f', &
                faces(1,i), '//', faces(1,i), &
                faces(2,i), '//', faces(2,i), &
                faces(3,i), '//', faces(3,i)
        end do

        close(iunit)

    end subroutine write_mesh

end module io
