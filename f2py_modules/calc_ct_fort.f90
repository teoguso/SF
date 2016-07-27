!!! Written by Matteo Guzzo ###
!!! A.D. MMXVI (2016)       ###
      subroutine calc_ct_fort(ct,im,en,nen,xt,nt)
      implicit none
      !integer nen,npoles
      !real*8 enexp(0:nen-1),akb(0:npoles-1),omegap(0:npoles-1)
      integer :: nen,nt,it,ien
      double precision, dimension (0:nen-1), intent(in) :: en, im
      complex(8), dimension (0:nt-1), intent(out) :: ct
      double precision, dimension (0:nen-1) :: den
      double precision, dimension (0:nt-1) :: xt
      complex(8), dimension (0:nen-1) :: integrand
      complex(8), PARAMETER :: xj=(0,1d0)
      double precision, PARAMETER :: pi=3.141592653589793
      double precision, parameter :: half=1d0/2d0
      double precision, parameter :: fourth=1d0/4d0
      double precision, parameter :: invpi=1d0/pi

!!!$OMP PARALLEL DEFAULT(SHARED) REDUCTION (+:ct) !!! PRIVATE(integrand)
      den = 1d0/en/en
      ct = 0d0
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(integrand,it,ien)
!$OMP DO 
      do it=0,nt-1
          do ien=0,nen-1
              integrand(ien) = im(ien)*(cos(en(ien)*xt(it))+xj*sin(en(ien)*xt(it))- xj*en(ien)*xt(it) - 1d0)*den(ien)
!         integrand(it) = im*(exp(xj*en*xt(it)) - xj*en*xt(it) - 1)*den
          end do
         
!         ct(it) = sum((integrand(1+1:nen-0) + integrand(1+0:nen-1))*(en(1+1:nen-0) - en(1+0:nen-1)))
          do ien=1,nen-2
            ct(it) = ct(it) + integrand(ien) 
          end do
!         ct(it) = sum(integrand) 
          ct(it) = ct(it) + 0.5d0*(integrand(1)+integrand(nen-1))
          ct(it) = ct(it)*invpi*half*(en(1)-en(0))
         
      end do
!     r = sum((y(1+1:n-0) + y(1+0:n-1))*(x(1+1:n-0) - x(1+0:n-1)))/2
!$OMP END DO 
!$OMP END PARALLEL         
!      ct = ct*invpi*half*(en(1)-en(0))
    
      !return ct
      end subroutine calc_ct_fort
