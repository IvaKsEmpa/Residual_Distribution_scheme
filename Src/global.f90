!----------------------------------------------------------------------------------------------------	
MODULE GlobalParam
  USE precisions
  IMPLICIT NONE
  INTEGER:: argunit = 6,  MyUnit = 30
  INTEGER::  Nx, iterfinal, cond_lim, solver, ImpreE, ImpreF
  INTEGER:: nt, ns, order
  REAL (KIND = DP)::  CFL, TIMEOUT
  REAL (KIND = DP)::  Lx, period_time 
  REAL (KIND = DP), PARAMETER:: EPS = 1.d-8
  REAL (KIND = DP), PARAMETER:: pi = 4.0d0*ATAN(1.0d0)
  REAL (KIND = DP), PARAMETER:: gamma = 1.4d0
  INTEGER, PUBLIC, PARAMETER:: n_dim  = 1   ! number of physical dimensions

  INTEGER , PARAMETER  :: Nv_Prim = 4, ro_pv= 1
  INTEGER , PARAMETER  ::  ein_pv = 2, press_pv= 3
  INTEGER , PARAMETER  :: sound_pv= 4
  INTEGER , PARAMETER  :: Lim_MinMod=1, Lim_SuperBee=2, Lim_VanLeer=3

 !_____________________________________________________________________________
  TYPE, PUBLIC:: variable_ther
     REAL(dp):: v(Nv_prim)!ro,e,p,c
  END TYPE variable_ther
 !_____________________________________________________________________________
  TYPE, PUBLIC:: variable
     TYPE(variable_ther), DIMENSION(:), ALLOCATABLE:: prim, prim_n, res_p
     REAL(dp), DIMENSION(:), ALLOCATABLE:: v,v_n,  res_v
  END TYPE variable
   !_____________________________________________________________________________
  TYPE, PUBLIC:: maillage
     INTEGER:: ns, nt, nseg
     REAL(dp), DIMENSION(:), ALLOCATABLE:: x
     INTEGER, DIMENSION(:,:), ALLOCATABLE:: elem
     INTEGER, DIMENSION(:,:), ALLOCATABLE:: nutv, edge
  END TYPE maillage
 !_____________________________________________________________________________

   TYPE, PUBLIC:: element
     INTEGER:: type_flux=-10
     INTEGER:: diag=-1
     INTEGER:: nsommets, itype, nvertex ! nbre de dofs, type element:
     !1->P1
     !2->B2
     !3->P2
     !4->P3
     REAL(8),  DIMENSION(:), POINTER :: coor =>NULL() ! 
     !   *-----*-----*
     !   1     3     2
     REAL(8), DIMENSION(:,:), POINTER:: x=> NULL()
     INTEGER, DIMENSION(:), POINTER :: nu =>NULL() ! local connectivity table, nu(1:3) are the 
     ! indices of the 3 vertices, nu(4) corresponds to the dof #4 and so on
     REAL(8)                           :: volume =0.0  ! volume
     REAL(8),  DIMENSION(:), POINTER :: n =>NULL()     ! external normals 
     INTEGER                        :: log    ! logique element
     !  : this for boundary conditions
!!!!   quadrature de surface
     REAL(8),   DIMENSION(:,:),POINTER :: quad =>NULL()  ! point de quadrature 
     REAL(8),     DIMENSION(:),POINTER :: weight=>NULL() ! poids 
     INTEGER                        :: nquad=0  ! nbre de points de quadrature
     REAL(8),DIMENSION(:,:),POINTER:: base0=>NULL(),base1=>NULL()
     INTEGER, DIMENSION(:), POINTER :: dof2ind
     CONTAINS
     PRIVATE

     PROCEDURE, PUBLIC:: base=>base_element

     FINAL:: clean
   
  END TYPE element

  CONTAINS
 !_____________________________________________________________________________
  REAL(8) FUNCTION aire_element(e) ! area of a element
    CLASS(element), INTENT(in):: e
    REAL(8), DIMENSION(1):: a,b
    a= e%coor(2)-e%coor(1)
    aire_element=ABS ( a(1))
  END FUNCTION aire_element
 !_____________________________________________________________________________
  FUNCTION normale_element(e) RESULT(n) ! inward normals
    CLASS(element), INTENT(in)::e
    REAL(8), DIMENSION(2):: n
    INTEGER:: l
    DO l=1,2
       n(1)= 1.0 ! at e%nu(1), i.e. x_{i+1}
       n(2)=-1.0 ! at e%nu(2), i.e. x_{i}
    ENDDO
  END FUNCTION normale_element
 !_____________________________________________________________________________
  REAL(8) FUNCTION base_element(e,k,x) ! basis functions  
    CLASS(element), INTENT(in):: e
    INTEGER, INTENT(in):: k ! index of basis function
    REAL(8), DIMENSION(2), INTENT(in):: x ! barycentric coordinate

    SELECT CASE(e%itype)
    CASE(1) ! P1
       SELECT CASE(k)
       CASE(2)
          base_element=x(1)
       CASE(1)
          base_element=x(2)
       CASE default
          PRINT*, "P1, numero base ", k
          STOP
       END SELECT

    CASE(2) ! B2 Bezier
       SELECT CASE(k)
       CASE(1)
          base_element=(1.0-x(1))*(1.0-x(1))
       CASE(2)
          base_element=x(1)*x(1)
       CASE(3)
          base_element=2.0*x(1)*(1.0-x(1))
       CASE default
          PRINT*, "B2, numero base ", k
          STOP
       END SELECT

    CASE(3)! P2
       SELECT CASE(k)
       CASE(1)
          base_element=(1.0-x(1))*(1.00-2.0*x(1))
       CASE(2)
          base_element=x(1)*(2.0*x(1)-1.00)
       CASE(3)
          base_element=4.00*x(1)*(1.00-x(1))
       CASE default
          PRINT*, "P2, numero base ", k
          STOP
       END SELECT

    CASE(4) ! P3
       SELECT CASE(k)
       CASE(1)
          base_element = -0.50*(3.0*x(1)-1.00)*(3.*x(1)-2.0)*(x(1)-1.0)
       CASE(2)
          base_element = 0.5*x(1)*(3.*x(1)-1.0)*(3.0*x(1)-2.0)
       CASE(3)
          base_element = 1.5*x(1)*(3.*x(1)-2.0)*(3.0*x(1)-3.0)
       CASE(4)
          base_element = -1.5*x(1)*(3.*x(1)-1.0)*(3.*x(1)-3.0)
       CASE default
          PRINT*, "P3, numero base ", k
          STOP
       END SELECT

    CASE(5) ! B3
       SELECT CASE(k)
       CASE(1)
          base_element =x(2)**3.0! (1.-x(1))**3
       CASE(2)
          base_element = x(1)**3.0
       CASE(3)
          base_element = 3.0*x(1)*x(2)*x(2)!3.*x(1)*( (1-x(1))**2)
       CASE(4)
          base_element = 3.0*x(1)*x(1)*x(2)!3.*x(1)*x(1)*((1-x(1)))
       CASE default
          PRINT*, "P3, numero base ", k
          STOP
       END SELECT
    CASE(6) ! B4
       SELECT CASE(k)
       CASE(1)
          base_element = x(2)**4.0! (1.-x(1))**3
       CASE(2)
          base_element = x(1)**4.0
       CASE(3)
          base_element = 4.0*x(1)*x(2)*x(2)*x(2)!3.*x(1)*( (1-x(1))**2)
       CASE(4)
          base_element = 6.0*x(1)*x(1)*x(2)*x(2)!3.*x(1)*x(1)*((1-x(1)))
       CASE(5)
          base_element = 4.0*x(1)*x(1)*x(1)*x(2)!3.*x(1)*x(1)*((1-x(1)))
       CASE default
          PRINT*, "B4, numero base ", k
          STOP
       END SELECT
    CASE default
       PRINT*, "Type non existant", e%itype
       STOP
    END SELECT

  END FUNCTION base_element
  !_____________________________________________________________________________

  FUNCTION gradient_element(e,k,x) RESULT (grad) ! gradient in reference element
    CLASS(element), INTENT(in):: e
    INTEGER, INTENT(in):: k ! numero de la fonction de base
    REAL(8), DIMENSION(2), INTENT(in):: x ! coordonnees barycentriques
    REAL(8),DIMENSION(n_dim):: grad
    REAL(8):: fx,fy,fz

    SELECT CASE(e%itype)
    CASE(1)! P1
       SELECT CASE(k)
       CASE(1)
          fx=-1.0
       CASE(2)
          fx=1.0
       END SELECT
    CASE(3) !P2
       SELECT CASE(k)
       CASE(1)
          fx=-3.0+4.0*x(1)
       CASE(2)
          fx=4.0*x(1)-1.0
       CASE(3)
          fx=4.0-8.0*x(1)
       END SELECT
    CASE(2) ! B2
       SELECT CASE(k)
       CASE(1)
          fx=-2.0*(1.0-x(1))
       CASE(2)
          fx=2.0*x(1)
       CASE(3)
          fx=2.0-4.0*x(1)
       END SELECT
    CASE(4) ! P3
       SELECT CASE(k)
       CASE(1)
          fx=-11./2.+(-27./2.*x(1)+18.)*x(1)
       CASE(2)
          fx=1.0+x(1)*(27./2.*x(1)-9.)
       CASE(3)
          fx=9.+x(1)*(81./2.*x(1)-45.)
       CASE(4)
          fx=-9./2.+x(1)*(-81./2.*x(1)+36.)
       END SELECT
    CASE(5)
       SELECT CASE(k)
       CASE(1)
          fx=-3.0*x(1)*x(1)+6.0*x(1)-3.0

       CASE(2)
          fx=3.0*x(1)*x(1)
       
       CASE(3)
          fx=9.0*x(1)*x(1)-12.0*x(1)+3.0
      
       CASE(4)
          fx= -9.0*x(1)*x(1)+6.0*x(1)

       END SELECT
    CASE(6)
       SELECT CASE(k)
       CASE(1)
      
          fx=-4.0*x(2)*x(2)*x(2)
       CASE(2)
     
          fx= 4.0* x(1)*x(1)*x(1)
       CASE(3)
     
          fx= 4.0*( x(2)-3.0*x(1) )*x(2)*x(2)
       CASE(4) 
     
          fx=12.0 * x(1)*x(2) * ( x(2)-x(1) )
       CASE(5)
          
          fx=4.0*( 3.*x(2)-x(1) )*x(1)*x(1)
       END SELECT
    CASE default
       PRINT*, "Type non existant", e%itype
       STOP
    END SELECT
    grad=fx/e%volume
  END FUNCTION gradient_element
!______________________________________________________________________________

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

!       e%quad(2,:)=1.0d0-e%quad(1,:)

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

  !  e%weight = e%weight/SUM(e%weight)
  END SUBROUTINE quadrature_element

 !_____________________________________
   REAL(8) FUNCTION eval_func_element(e, u, y)
    CLASS(element),                       INTENT(in):: e
    REAL(8),    DIMENSION(e%nsommets), INTENT(in):: u
    REAL(8),       DIMENSION(2),          INTENT(in):: y
    REAL(8),       DIMENSION(2)                     :: x
    REAL(8),       DIMENSION(e%nsommets)            :: alpha, beta, base
    REAL(8)                                     :: a,b,c, aa, cc
    INTEGER                                         :: l
    LOGICAL                                         :: flag

       DO l=1, e%nsommets
          base(l)=e%base(l,y)
       ENDDO
       DO l=1,e%nsommets
          eval_func_element=SUM(base(:)*u(:))
       END DO
  END FUNCTION eval_func_element

 !_____________________________________________________________________________

  SUBROUTINE clean(e)
    TYPE(element), INTENT(inout)::e
    IF (ASSOCIATED(e%coor))  NULLIFY(e%coor)
    IF (ASSOCIATED(e%n))  NULLIFY(e%nu)
    IF (ASSOCIATED(e%n))  NULLIFY(e%n)
    IF (ASSOCIATED(e%quad))  NULLIFY(e%quad)
    IF (ASSOCIATED(e%weight))  NULLIFY(e%weight)
    IF (ASSOCIATED(e%base0))  NULLIFY(e%base0)
    IF (ASSOCIATED(e%base1))  NULLIFY(e%base1)
  END SUBROUTINE clean

 !_____________________________________________________________________________

END MODULE GlobalParam
!-------------------------------------------------------------------------------------------------------------

