module pgm

! Module for reading and saving pgm images.

use basic, only: stop_error
implicit none
private
public loadpgm, savepgm

contains

    subroutine loadpgm(filename, img)
        ! Loads the image 'img' from 'filename' in pgm format.
        character(len=*), intent(in) :: filename
        ! img(c, x, y) is the c-component (c=1,2,3 is R, G, B colors) of the
        ! (x, y) pixel. The image is in the 1st quadrant, so (1, 1) is the lower left
        ! corner, (w, h) is the upper right corner.
        integer, intent(out), allocatable :: img(:, :)

        character(len=2) :: signature
        character :: ccode
        integer :: w, h, i, j, ncol
        integer :: u, offset

        open(newunit=u, file=filename, access="stream", form="formatted", status="old")
        read(u, '(a2)') signature
        if (signature /= "P5") call stop_error("Invalid format. expected 'P5'")
        read(u, *) w, h
        read(u, *) ncol

        !if (ncol /= 255) call stop_error("Unsupported color range.")
        inquire(u, pos=offset)
        close(u)

        open(newunit=u, file=filename, access="stream", status="old")
        ! Move the corrent position to "offset":
        read(u, pos=offset-1) ccode
        allocate(img(w, h))

        do j = 1, h
            do i = 1, w
                read(u) ccode
                img(i, j) = iachar(ccode)
            end do
        end do

        close(u)

    end subroutine

    subroutine savepgm(filename, img)
        ! Saves the image 'img' into 'filename' in pgm format.
        character(len=*), intent(in) :: filename
        ! img(c, x, y) is the c-component (c=1,2,3 is R, G, B colors) of the
        ! (x, y) pixel. The image is in the 1st quadrant, so (1, 1) is the lower left
        ! corner, (w, h) is the upper right corner.
        integer, intent(in) :: img(:, :)

        integer :: h, w, i, j
        integer :: u

        w = size(img, 1)
        h = size(img, 2)

        open(newunit=u, file=filename, status="replace")
        write(u, '(a2)') "P5"
        write(u, '(i0," ",i0)') w, h
        write(u, '(i0)') 255

        do j = 1, h
            do i = 1, w
                write(u, '(a1)', advance='no') achar(img(i, j))
            end do
        end do

        close(u)

    end subroutine

end module
