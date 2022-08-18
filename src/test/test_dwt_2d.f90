program main
    !test 2d wavelet transform
    CHARACTER(100) :: num1char
    integer :: num1,num2

    CALL GET_COMMAND_ARGUMENT(1,num1char)
    READ(num1char,*) num1 
    CALL GET_COMMAND_ARGUMENT(2,num1char)
    READ(num1char,*) num2 
    
    call test_fdwt_2d(num1,num2)
    
end program main

subroutine test_fdwt_2d(K, depth)
    use wavelet_transform
    use basic
    use wavelet_bank
    use pgm 

    
    integer, intent(in) ::K, depth

    real , dimension(:), allocatable:: h_phi, h_psi, g_phi, g_psi
    character (len=*), parameter :: img_in_file='data/ruler.pgm'
    integer, dimension(:,:), allocatable :: img
    real, dimension(:,:), allocatable :: img_real, w_coefs, inv

    write(*,*) 'loading data' 
    
    call loadpgm(img_in_file, img)
   
    h_phi= daubechies_father_scaling_coefs(K)
    !h_phi = haar_father_scaling_coefs()
    h_psi = wavelet_function_coefs(h_phi)

    g_phi = h_phi(size(h_phi): 1: -1)
    g_psi = h_psi(size(h_psi): 1: -1)
    

    img_real = real(img)

    write(*,*) 'computing wlet transform' 

    !write(*,*) 'img', shape(img_real)
    w_coefs = fdwt_2d(img_real, h_phi, h_psi, depth)

    write(*,*) 'saving'

    call savepgm('output/ruler_wt.pgm', int(w_coefs))

    write(*,*) 'saving'
    write(*,*) 'computing inverse wlet transform' 
    

    N=512
    M=512
    w_coefs(1:int(N/2):1, 1:int(M/2):1)=0

    inv = fidwt_2d(w_coefs, h_phi, h_psi, depth)

    call savepgm('output/ruler_iwt.pgm', int(inv))

    write(*,*) 'all done!'

end subroutine test_fdwt_2d