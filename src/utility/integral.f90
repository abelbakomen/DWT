module integrals
    implicit none 

    interface trapz
        module procedure trapz1D
        module procedure trapz2D
    end interface
    
contains

function milieu1D(c_1D_func, a, b, step) result (I)
    implicit none
    
    real, intent(in) :: a,b
    real, optional, intent(in) :: step
    real :: I,h
    integer :: j, n

    interface
        function c_1D_func (t) result (y)
            real, intent(in) :: t
            real :: y 
        end function
    end interface

    
    h = b-a
    if(present(step)) h = step
    
    n = int((b-a)/h)

    I=0.
    do j = 1,n
        I = I + h * c_1D_func(a + j*h/2)
    end do
end function

function trapz1D(f,step) result (I)
    implicit none
    
    real, dimension(:), intent(in) :: f
    real, optional, intent(in) :: step
    real :: I

    real ::  h
    integer :: n
    
    h = 1
    if(present(step)) h = step
    
    n = size(f)
    I = h*((f(1) + f(n))/2 + sum( f(2: (n-1)) ) )

end function

function trapz2D(f,step1,step2) result (I)
    implicit none
    
    real, dimension(:,:), intent(in) :: f
    real, optional, intent(in) :: step1, step2
    real :: I

    !real, dimension (:), allocatable :: Iy 
    real:: h1,h2
    integer :: j
    
    h1 = 1
    h2 = 1 
    if(present(step1)) h1=step1
    if(present(step2)) h2=step2
    
    !allocate(Iy(size(f,2)))

    !do j = 1, size(f,2)
    !  Iy(j) = trapz( f(:, j), step1)
    !end do        
 

    !I=trapz(Iy,step2)

    I=sum(f)

end function

function h_phi_(father_scaling, k, a, b, step) result (y)
    implicit none
    

    integer, intent(in) :: k
    real, intent(in) :: a, b
    real, optional, intent(in) :: step
    real :: h,y
    integer :: j, n

    interface
        function father_scaling (t) result (y)
            real, intent(in) :: t
            real :: y 
        end function
    end interface

    
    h = b-a
    if(present(step)) h = step
    
    n = int((b-a)/h)

    y=0.
    do j = 1,n
        y = y + h * (father_scaling(a + j*h/2) * sqrt(2.) * father_scaling(2*(a + j*h/2)-k))
    end do
end function

end module
