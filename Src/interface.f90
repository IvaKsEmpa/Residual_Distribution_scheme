MODULE ModuleInterface
  USE GlobalParam
  USE precisions
  IMPLICIT NONE
CONTAINS

  FUNCTION fex(p4, p1, p5, rho1, rho5) RESULT(f)
    IMPLICIT NONE
    !shock tube equation
    !declare the pass
    REAL(kind=dp):: p4,p1,p5,rho1,rho5, f
    !local variables
    REAL(kind=dp):: z,c1,c5,gm1,gp1,g2,fact
    z = (p4/p5 - 1.d0)
    c1 = dsqrt(gamma*p1/rho1)
    c5 = dsqrt(gamma*p5/rho5)

    gm1 = gamma - 1.d0
    gp1 = gamma + 1.d0
    g2  = 2.d0*gamma

    fact = gm1/g2*(c5/c1)*z/dsqrt(1.d0 + gp1/g2*z)
    fact = (1.d0 - fact)**(g2/gm1)
    f = p1*fact - p4
    RETURN
  END FUNCTION fex
  !--------------------------------------------------------
  FUNCTION InternalEn(ro, press) RESULT(InternalE)
    IMPLICIT NONE
    REAL (KIND = DP):: ro, press, InternalE
    InternalE =  press/((gamma - 1.d0)) ! *ro

  END FUNCTION InternalEn
  !--------------------------------------------------------
  FUNCTION Pressure(ro, ein) RESULT( press)
    IMPLICIT NONE
    REAL (KIND = DP):: ro, ein,  press
    press =  ein*(gamma - 1.d0)!*ro
  END FUNCTION Pressure
  !--------------------------------------------------------
  FUNCTION TotalEn(ro, u, InternalE) RESULT(TotalE)
    IMPLICIT NONE
    REAL (KIND = DP):: ro, u, InternalE, TotalE
    TotalE =  InternalE + ro*U**2.D0/2.D0
  END FUNCTION TotalEn
  !--------------------------------------------------------
  FUNCTION SoundPolGas(ro, press) RESULT(Sound)
    IMPLICIT NONE
    REAL (KIND = DP):: ro, press, Sound
    Sound =  dsqrt(gamma*press/ro)
  END FUNCTION SoundPolGas
  !--------------------------------------------------------

Subroutine simpson(f,a,b,integral,n)
IMPLICIT NONE
REAL (KIND=DP)::f, a, b, integral,s
REAL (KIND=DP):: h, x
integer:: nint
integer:: n, i
!==========================================================
! Integration of f(x) on [a,b]
! Method: Simpson rule for n intervals  
!----------------------------------------------------------
! IN:
! f   - Function to integrate (supplied by a user)
! a   - Lower limit of integration
! b   - Upper limit of integration
! n   - number of intervals
! OUT:
! integral - Result of subintegration
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
! ----------------------------------------------------------------------
END MODULE ModuleInterface
!----------------------------------------------------------------------------------------------------
