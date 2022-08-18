module wavelet_bank
    implicit none

contains
  
    function wavelet_function_coefs(fthr_scalig_coef) result(h1)
        real, dimension(:), intent(in) :: fthr_scalig_coef 
        real, dimension(:), allocatable :: h1 
        integer :: K,n

        K = size(fthr_scalig_coef)

        h1 = (/(((-1)**n) * fthr_scalig_coef(K-n), n=0,K-1, 1)/)

    end function
    

    function haar_father_scaling_coefs() result(h)
        real, dimension(2) :: h 

        h(1) = 1/sqrt(2.)
        h(2) = 1/sqrt(2.)

    end function

    function daubechies_father_scaling_coefs(n) result(h)
        use basic
        integer, intent(in) :: n 
        real, dimension(:), allocatable :: h 

        allocate(h(n))
        select case (n)
            case (2)
                h = (1/sqrt(2.))*[1., 1.]
            case (4)
                h = (1/sqrt(2.))*[0.6830127, 1.1830127, 0.3169873, -0.1830127]
            case (6)
                h = (1/sqrt(2.))*[0.47046721, 1.14111692, 0.650365, -0.19093442, -0.12083221, 0.0498175]
            case (8)
                h = (1/sqrt(2.))*[0.32580343, 1.01094572, 0.8922014, -0.03957503, -0.26450717, 0.0436163, 0.0465036, -0.01498699]
            case (10)
                h = (1/sqrt(2.))*[0.22641898, 0.85394354, 1.02432694, 0.19576696, -0.34265671,&
                 -0.04560113, 0.10970265, -0.0088268, -0.01779187, 0.00471742793]
            case (12)
                h = (1/sqrt(2.))*[0.15774243, 0.69950381, 1.06226376, 0.44583132, -0.3199866,&
                -0.18351806, 0.13788809, 0.03892321, -0.04466375, 0.000783251152,&
                 0.00675606236, -0.00152353381]
            case (14)
                h = (1/sqrt(2.))*[0.11009943, 0.56079128, 1.03114849, 0.66437248, -0.20351382,&
                 -0.31683501, 0.1008467, 0.11400345, -0.05378245, -0.02343994,&
                  0.01774979, 0.000607514995, -0.00254790472, 0.000500226853]
            case (16)
                h = (1/sqrt(2.))*[0.07695562, 0.44246725, 0.95548615, 0.82781653, -0.02238574,&
                 -0.40165863, 0.000668194092, 0.18207636, -0.0245639, -0.06235021,&
                  0.01977216, 0.01236884, -0.00688771926, -0.000554004549,&
                   0.000955229711, -0.000166137261]
            case (18)
                h = (1/sqrt(2.))*[0.05385035, 0.3448343, 0.85534906, 0.92954571, 0.18836955,&
                 -0.41475176, -0.13695355, 0.21006834, 0.043452675, -0.09564726,&
                  0.000354892813, 0.03162417, -0.00667962023, -0.00605496058,&
                  0.00261296728, 0.000325814671, -0.000356329759, 5.5645514e-05]
            case (20)
                h = (1/sqrt(2.))*[0.03771716, 0.26612218, 0.74557507, 0.97362811, 0.39763774,&
                 -0.3533362, -0.27710988, 0.18012745, 0.13160299, -0.10096657,&
                  -0.04165925, 0.04696981, 0.00510043697, -0.015179, 0.00197332536,&
                   0.00281768659, -0.00096994784, -0.000164709006, 0.000132354367, -1.875841e-05]

           case default
                call stop_error('Unsuported Length')

        end select


    end function

end module

