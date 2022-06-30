module quadratures
USE GlobalParam
  USE precisions
  IMPLICIT NONE
CONTAINS
subroutine trapiter(T,n,a,b,f)
  INTEGER::  nbp
  REAL (KIND=DP):: h, somme
  real, intent(in) :: a,b
  real, intent(inout) :: T
  integer, intent(in) :: n
if (n==1) then
T=0.5d0*(b-a)*(f(a)+f(b))
else
nbp=2**(n-2)
h=(b-a)/nbp
somme=0.d0
do j=1,nbp
somme=somme+f(a+(j-0.5d0)*h)
end do
T=0.5d0*(T+somme*h)
end if
end subroutine trapiter

function quasimp(a,b,f,er,nmax) result(S)
real :: S,Sa,T,Ta
integer :: j
 call trapiter(T,1,a,b,f)
Ta=T
 call trapiter(T,2,a,b,f)
Sa=(4.*T-Ta)/3.
j=3
Ta=T
 call trapiter(T,3,a,b,f)
S=(4.*T-Ta)/3.
do while (abs(S-Sa)>er*abs(Sa) .and. j<nmax)
  j=j+1
Ta=T
 call trapiter(T,j,a,b,f)
Sa=S
S=(4.*T-Ta)/3.
end do
er=abs(S/Sa-1)
nmax=j
end function quasimp

Subroutine simpson(f,a,b,integral,n)
implicit none
REAL (KIND=DP)::f, a, b, integral,s
REAL (KIND=DP):: h, x
integer:: nint
integer:: n, i
!==========================================================
! Integration of f(x) on [a,b]
! Method: Simpson rule for n intervals  
!==========================================================

! if n is odd we add +1 to make it even
if((n/2)*2.ne.n) n=n+1

! loop over n (number of intervals)
s = 0.0
h = (b-a)/dfloat(n)
do i=2, n-2, 2
   x   = a+dfloat(i)*h
   s = s + 2.0*f(x) + 4.0*f(x+h)
end do
integral = (s + f(a) + f(b) + 4.0*f(a+h))*h/3.0
return
end subroutine simpson

 SUBROUTINE quadrature_element(e)  ! quadrature points and weights
    CLASS(element), INTENT(inout):: e
    REAL(8):: w,xo,yo,zo,s
    INTEGER:: nquad, nquad_1

    SELECT CASE(e%itype)
    CASE(1)
       e%nquad=2
       nquad=e%nquad
       ALLOCATE(e%quad(2,e%nquad),e%weight(e%nquad))

       ! quadrature for edges (Gauss formula, 3 points)

       e%quad(1,1)=1.0
       e%quad(2,1)=0.
       e%weight(1) = 0.5

       e%quad(1,2)=0.0
       e%quad(2,2)=1.0
       e%weight(2) = 0.5


    CASE(2,3)!! exact for degree 5
       e%nquad=3
       nquad=e%nquad
       ALLOCATE(e%quad(2,e%nquad),e%weight(e%nquad))



       s=0.7745966692414834  !SQRT(0.6d0)
       e%quad(1,1)=0.5d0*(1.0d0 - s)
       e%quad(2,1)=0.5d0*(1.0d0 + s)
       e%weight(1) = 0.2777777777777778 !5.0d0/18.0d0

       e%quad(1,2)=0.5d0*(1.0d0 + s)
       e%quad(2,2)=0.5d0*(1.0d0 - s)
       e%weight(2) = 0.2777777777777778 !5.0d0/18.0d0

       e%quad(1,3)=0.5d0
       e%quad(2,3)=0.5d0
       e%weight(3)= 0.4444444444444444 ! 8.0d0/18.0d0

    CASE(4,5) !! exact for degree 7
       e%nquad=4 ! ordre 7
       nquad=e%nquad
       ALLOCATE(e%quad(2,e%nquad),e%weight(e%nquad) )
       s=0.1739274225687269  !(18.- SQRT(30.))/72. !SQRT(30.d0)*(-5.0d0+3.d0*SQRT(30.d0))/360.d0
       e%weight(1:2)=s
       s=0.8611363115940526!SQRT(525.d0+70.d0*SQRT(30.d0))/35.d0
       e%quad(1,1)=0.5d0*(1.0d0-s)
       e%quad(1,2)=0.5d0*(1.0d0+s)
       e%quad(2,1)=e%quad(1,2) !0.5d0*(1.0d0+s)
       e%quad(2,2)=e%quad(1,1) !0.5d0*(1.0d0-s)

       s=0.3260725774312731 !(18.+  SQRT(30.))/72.! SQRT(30.d0)*(5.0d0+3.d0*SQRT(30.d0))/360.d0
       e%weight(3:4)=s
       s=0.3399810435848563 !SQRT(525.d0-70.d0*SQRT(30.d0))/35.d0
       e%quad(1,3)=0.5d0*(1.0d0-s)
       e%quad(1,4)=0.5d0*(1.0d0+s)
       e%quad(2,3)=e%quad(1,4) !0.5d0*(1.0d0+s)
       e%quad(2,4)=e%quad(1,3) !0.5d0*(1.0d0-s)


    CASE(6) ! exact for degree 9

       e%nquad=5
       nquad=e%nquad
       ALLOCATE(e%quad(2,e%nquad),e%weight(e%nquad))

       s=1.d0/3.d0*SQRT(5.d0-2.d0*SQRT(10.d0/7.d0))
       e%quad(1,1)=0.5d0*(1.0d0 - s)
       e%quad(2,1)=0.5d0*(1.0d0 + s)
       e%weight(1) = (322.d0+13.d0*SQRT(70.0d0))/1800.0d0
       e%quad(1,2)=0.5d0*(1.0d0 + s)
       e%quad(2,2)=0.5d0*(1.0d0 - s)
       e%weight(2) = (322.d0+13.d0*SQRT(70.0d0))/1800.0d0

       s=1.d0/3.d0*SQRT(5.d0+2.d0*SQRT(10.0d0/7.0d0))
       e%quad(1,3)=0.5d0*(1.0 - s)
       e%quad(2,3)=0.5d0*(1.0 + s)
       e%weight(3) = (322.0d0-13.0d0*SQRT(70.0d0))/1800.0d0
       e%quad(1,4)=0.5d0*(1.0d0 + s)
       e%quad(2,4)=0.5d0*(1.0d0 - s)
       e%weight(4) = (322.0d0-13.0d0*SQRT(70.0d0))/1800.0d0

       e%quad(1,5)=0.5d0
       e%quad(2,5)=0.5d0
       e%weight(5)=64.0d0/225.0d0
    END SELECT
  END SUBROUTINE quadrature_element

end module quadratures
! ----------------------------------------------------------------------
