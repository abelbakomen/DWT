module wavelet_transform
implicit none

contains
    recursive function fdwt_1d(d_signal, h_phi, h_psi) result (w_coefs)
        use basic
        
        implicit none

        real, intent(in), dimension(:) :: d_signal, h_phi, h_psi
        real, dimension(:), allocatable :: A, D, w_coefs
        integer :: N,J,i,k

        N = size(d_signal)
        J = int(log(0.+N)/log(2.))
        k = size(h_phi)

        A = conv_1d(d_signal, h_phi(k:1:-1))
        A = A(2: N: 2)
        D = conv_1d(d_signal, h_psi(k:1:-1))
        D = D(2: N: 2)

        allocate(w_coefs(N))

        if(N <=  2) then
            w_coefs(1 : (int(N/2)) ) = A 
        else
            A = fdwt_1d(A, h_phi, h_psi)
            w_coefs(1 : (int(N/2)) ) = A(:int(N/2))
        end if

        w_coefs((int(N/2)+1):N ) = D
        
    end function
    
    recursive function fidwt_1d(w_coefs, g_phi, g_psi) result (d_signal)
        use basic
        
        implicit none

        real, intent(in), dimension(:) :: w_coefs, g_phi, g_psi
        real, dimension(:), allocatable :: A, D, d_signal
        integer :: N,i


        N = size(w_coefs)
        
        allocate(D(N))

        allocate(A(N))
        allocate(d_signal(N))

        D=0
        A=0

        D(1: N: 2) = w_coefs( (int(N/2)+1): N: 1) 

        if(N<=2) then
            A(1: N: 2) = w_coefs(1: int(N/2): 1)
        else
            A( 1: N: 2) = fidwt_1d( w_coefs( 1: int(N/2): 1), g_phi, g_psi) 
        end if 

        D = conv_1d(A, g_phi, 1) + conv_1d( D, g_psi, 1 ) 
        
        d_signal = D(:N)
        
    end function

    recursive function fdwt_2d(d_signal, h_phi, h_psi, dpth) result (w_coefs)
        use basic
        
        implicit none

        real, intent(in), dimension(:,:) :: d_signal
        real, intent(in), dimension(:) ::  h_phi, h_psi
        integer, intent(in) :: dpth
        real, dimension(:), allocatable ::  Temp
        real, dimension(:,:), allocatable :: A0, A, D0, D, H, V, w_coefs
        integer :: N, M, J, i, depth

        !write(*,*) 'd_signal', shape(d_signal)

        N = size(d_signal, 1)
        M = size(d_signal, 2)

        depth = max(dpth, 1)
        depth = min(depth, int(log(N+0.)/log(2.)))

        allocate(A0(int(N/2), M))
        allocate(D0(int(N/2), M))

        do  i = 1, M
            Temp = conv_1d(d_signal(:,i), h_phi(size(h_phi): 1: -1))
            
            A0(:,i) = Temp(2:N:2)

            Temp = conv_1d(d_signal(:,i), h_psi(size(h_psi): 1: -1))
            D0(:,i) = Temp( 2:N:2 )
        end do
!--------
        allocate(A(int(N/2), int(M/2)))
        allocate(D(int(N/2), int(M/2)))
        
        allocate(V(int(N/2), int(M/2)))
        allocate(H(int(N/2), int(M/2)))
        
        do  i = 1, int(N/2)
            Temp = conv_1d(A0(i, :), h_phi(size(h_phi): 1: -1))
            A(i, :) = Temp(2:M:2 )

            Temp =  conv_1d(A0(i, :), h_psi(size(h_psi): 1: -1))
            H(i, :) = Temp( 2:M:2 )

            Temp =  conv_1d(D0(i, :), h_phi(size(h_phi): 1: -1))
            V(i, :) = Temp( 2:M:2 )

            Temp =  conv_1d(D0(i, :), h_psi(size(h_psi): 1: -1))
            D(i, :) = Temp( 2:M:2 )
        end do

        if(depth > 1) A = fdwt_2d(A, h_phi, h_psi, depth-1)

        allocate(w_coefs(N, M))
            
        w_coefs(1 : (int(N/2)), 1 : (int(M/2)) ) = A
        w_coefs(1 : (int(N/2)), (int(M/2)+ 1) : M ) = V
        w_coefs((int(N/2)+ 1) : N, 1 : (int(M/2))) = H
        w_coefs((int(N/2)+ 1) : N, (int(M/2)+ 1) : M ) = D

        deallocate(A0)
        deallocate(D0)
        deallocate(A)
        deallocate(H)
        deallocate(V)
        deallocate(D)
    end function


    recursive function fidwt_2d(w_coefs, g_phi, g_psi, dpth) result (d_signal)
        use basic
        
        implicit none

        real, intent(in), dimension(:,:) :: w_coefs
        real, intent(in), dimension(:) ::  g_phi, g_psi
        integer, intent(in):: dpth
        real, dimension(:), allocatable :: Temp
        real, dimension(:,:), allocatable :: A0, A, D0, D, H, V, d_signal
        integer :: N,M,i,depth

    
        N = size(w_coefs, 1)
        M = size(w_coefs, 2)

 
        depth = max(dpth, 1)
        depth = min(depth, int(log(N+0.)/log(2.)))

        write(*,*) 'ok there', N, M
        allocate(A(int(N/2), M))
        allocate(D(int(N/2), M))
        allocate(V(int(N/2), M))
        allocate(H(int(N/2), M))
        
        A = 0.
        V = 0.
        H = 0.
        D = 0.

        A(:, 1: M: 2) = w_coefs(1 : (int(N/2)), 1 : (int(M/2)))     
        
        if(depth > 1) A(:, 1: M: 2) = fidwt_2d(w_coefs(1 : (int(N/2)), 1 : (int(M/2)) ), g_phi, g_psi, depth-1)
        !if(N==32 .or. N==128) A=10 
        V(:, 1: M: 2) = w_coefs(1 : (int(N/2)), (int(M/2)+ 1) : M ) 
        H(:, 1: M: 2) =  w_coefs((int(N/2)+ 1) : N, 1 : (int(M/2)) )  
        D(:, 1: M: 2) = w_coefs((int(N/2)+ 1) : N, (int(M/2)+ 1) : M )
        

        allocate(A0(N, M))
        allocate(D0(N, M))

        A0 = 0.
        D0 = 0.

        do  i = 1, int(N/2)
            Temp = conv_1d(A(i, :), g_phi, 1) + conv_1d(H(i, :), g_psi, 1)
            A0(2*i-1, :) = Temp(1:M)

            Temp =conv_1d(V(i, :), g_phi, 1) + conv_1d(D(i,:), g_psi, 1)
            D0(2*i-1 ,:) = Temp(1:M)
        end do

        allocate(d_signal(N, M))
        do  i = 1, M
            Temp = conv_1d(A0(:,i), g_phi, 1) + conv_1d(D0(:,i), g_psi, 1)
            d_signal(:, i) = Temp(1:N)
        end do
        
        deallocate(A0)
        deallocate(D0)
        deallocate(A)
        deallocate(H)
        deallocate(V)
        deallocate(D)
    end function


end module 