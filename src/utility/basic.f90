module basic
implicit none


contains

    ! Generates evenly spaced numbers from `from` to `to` (inclusive).
    ! 
    ! Inputs:
    ! -------
    ! 
    ! from, to : the lower and upper boundaries of the numbers to generate
    ! 
    ! Outputs:
    ! -------
    ! 
    ! array : Array of evenly spaced numbers
    ! 
    subroutine linspace(from, to, array)
    
        real, intent(in) :: from, to
        real, intent(out) :: array(:)
        real :: range
        integer :: n, i
        n = size(array)
        range = to - from

        if (n == 0) return

        if (n == 1) then
            array(1) = from
            return
        end if


        do i=1, n
            array(i) = from + range * (i - 1) / (n - 1)
        end do
    end subroutine  

    subroutine stop_error(msg)
    ! Aborts the program with nonzero exit code
    !
    ! The statement "stop msg" will return 0 exit code when compiled using
    ! gfortran. stop_error() uses the statement "stop 1" which returns an exit code
    ! 1 and a print statement to print the message.
    !
    ! Example
    ! -------
    !
    ! call stop_error("Invalid argument")

        character(len=*) :: msg ! Message to print on stdout
        print *, msg
        stop 1
    end subroutine


    subroutine assert(condition)
    ! If condition == .false., it aborts the program.
    !
    ! Arguments
    ! ---------
    !
        logical, intent(in) :: condition
    !
    ! Example
    ! -------
    !
    ! call assert(a == 5)

        if (.not. condition) call stop_error("Assert failed.")
    end subroutine



    function conv_1d(d_sgl,kernel, direction) result(cv)
        implicit none

        real, intent(in), dimension(:) ::d_sgl ,kernel
        integer , intent(in), optional :: direction!0= LTR, 1=RTL
        real, dimension(:), allocatable ::d_signal, padded_signal, cv
        integer :: dir, n, k, i, q

        dir=0
        if ( present(direction) ) dir = direction

        n = size(d_sgl)
        k = size(kernel)

        d_signal = d_sgl
        if(dir == 1) d_signal = d_sgl(n:1:-1)

        !write(*,*) n,k

        allocate(padded_signal(n + k))

        padded_signal(:n)=d_signal
        
        q = int((n+k)/n)

        if(n>=k)then
            padded_signal((n+1):(n+k)) = d_signal(1:k)
        else
            do i = 0,q-1,1
                padded_signal((i*n +1 ): ((i+1)*n)) = d_signal
            end do
            padded_signal((q*n+1):) = d_signal(1: (n+k - q*n ))
        end if
        
        !padded_signal=0
        !padded_signal((int(k/2)+1):(n+int(k/2)))=d_signal
        
        if(dir == 1) padded_signal = padded_signal((n+k):1:-1)

       ! write(*,*) d_sgl
        !write(*,*) d_sgl((n-k)+1:)
        !write(*,*) padded_signal

        !write(*,*) ''

        cv = (/( dot_product(padded_signal(i:(i+k-1)) , kernel(k:1:-1)), i = 1, n+1, 1)/)


        !write(*,*) size(cv)

        if(dir == 1) then
            cv = cv(1:(n))
        else
            cv = cv(1:n)
        end if        

        !write(*,*) size(cv)

    end function
end module