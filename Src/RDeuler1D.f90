
PROGRAM main
  USE precisions
  USE GlobalParam
  USE InOutputs
  USE ModuleInterface
  USE scheme
  IMPLICIT NONE  
  INTEGER:: iv,  IT !I,ix,, n_thermo
  TYPE(maillage):: mesh
  TYPE(variable):: var

  REAL (KIND=DP), ALLOCATABLE, DIMENSION(:):: MaxVP, MinVp, einstar, alpha
  REAL (KIND=DP), ALLOCATABLE, DIMENSION(:):: rostar, pstar!, ronplus1, ubar
  !  REAL (KIND=DP), ALLOCATABLE, DIMENSION(:):: utildatilda, rotildatilda, deltaro
  REAL (KIND=DP), ALLOCATABLE, DIMENSION(:):: pente!, rotilde, un, ron
  REAL(dp), DIMENSION(:,:), ALLOCATABLE:: flux
  REAL (KIND=DP) :: T1_CPU, T2_CPU, TIME,  TIME2
  REAL (KIND=DP):: Dx,  DT, dt2
  REAL (KIND=DP):: Cmax, lambda
  REAL(KIND=DP):: Xmax, XMIN,  UMAX, v0, vmax
  CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:):: NamesOfPrimVar

  !----------------------------------------------------------------------------------------------------
  WRITE(6,*) " Begin of the Program "
  CALL CPU_TIME(T1_CPU)
  !----------------------------------------------------------------------------------------------------
  CALL LECTURE_DONNEES()

  Mesh%ns=nx+1
  mesh%nt=nx
  mesh%nseg=mesh%ns

  ALLOCATE( var%prim(mesh%nt), var%res_p(mesh%nt), var%v(mesh%ns),var%res_v(mesh%ns))
  ALLOCATE(var%prim_n(mesh%nt), var%v_n(mesh%ns))
  ALLOCATE( NamesOfPrimVar(Nv_Prim), PENTE(1:mesh%nseg), FLUX(3,1:mesh%nseg))
  ALLOCATE( MaxVP(Nv_Prim), MinVp(Nv_Prim))
  ALLOCATE( rostar(Mesh%nseg), pstar(Mesh%nseg), einstar(mesh%nseg), alpha(Mesh%nseg) )

  Dx = Lx/REAL(mesh%nt,dp)

  PRINT*, 'Nx =', Nx, 'Dx =', Dx
  !----------------------------------------------------------------------------------------------------  
!!$  NamesOfPrimVar(ro_pv)    = "Ro"
!!$  NamesOfPrimVar(U_pv)     = "Velocity x"
!!$  NamesOfPrimVar(ein_pv)   = "int energy"
!!$  NamesOfPrimVar(press_pv) = "pressure"
!!$  NamesOfPrimVar(sound_pv) = "sound speed"
  !----------------------------------------------------------------------------------------------------
  DO iv = 1, Nv_Prim
     WRITE(6,*) NamesOfPrimVar(iv) , " <<<< position in 'Prim' array == ", iv, Nx
  ENDDO
  WRITE(6,*) " >>> End of LECTURE_DONNEES"

  !----------------------------------------------------------------------------------------------------
  !INITIALISATION  
  CALL INITIALISATION(var, mesh, dx)
  CALL lim_cond(var,mesh)
  CALL state(Var%prim)
  TIME = 0.D0
  TIME2  = period_time ;
  IT = 1
  CALL PutonScreen()
  !----------------------------------------------------------------------------------------------------
  XMIN = MINVAL(mesh%X(:))
  XMAX = MAXVAL(mesh%X(:))
  WRITE(6,'( 2(A10,E15.6))')  " XMIN = ", XMIN, &
       &                      " Xmax = ", XMAX  
  !---------------------------------------------------------------------------------------------------- 
  CALL Ecriture_donnees(var,mesh)
  Umax= -1.d0
  ! BOUCLE SUR LE TEMPS
  !----------------------------------------------------------------------------------------------------
  Time_loop: DO WHILE ((IT .LE. iterfinal) .AND. (TIME .LE. TIMEOUT))
     !----------------------------------------------------------------------------------------------------
  
     v0=MAXVAL( ABS(var%v) )
     Umax=MAX(Umax,  v0) ! MAXVAL(ABS(var%v-var%v(sound_pv))) !
     Cmax=MAXVAL( ABS( var%Prim(:)%v(sound_pv) ) ) ! MAXVAL(ABS(var%v+var%v(sound_pv))) !
     vmax=umax+cmax!MAX(Umax,Cmax)
    !vmax=maxval(abs(var%v)+var%Prim(:)%v(sound_pv))

     
     !----------------------------------------------------------------------------------------------------
     CALL PENTE_U( DX, Var,mesh, PENTE)

    if (solver == 1) then 
     CALL HLL(dx, var, mesh, FLUX,  pente, rostar, pstar, einstar, alpha)
    elseif (solver == 2) then 
     CALL HLLC(dx, var, mesh,  FLUX, pente, rostar, pstar, einstar, alpha)
   endif

    vmax= max( maxval(alpha), umax+cmax)
     !print*, maxval(alpha), umax+vmax
     dt = Dx/vmax!(umax + cmax)! vmax !
     DT = CFL*dt
     DT=MIN(DT,TIMEOUT-TIME)
     IF (dt.LE.1.e-10) PRINT*, 'dt.LE.1.e-10, program stop'
     IF (dt.LE.1.e-10) EXIT
     dt2 = 0.5d0*dt
     lambda = dt/dx
     IF (MOD(it,1000)==0)  WRITE(6,*) " Dt = ", Dt, 'time=', time, 'lambda =', lambda
     TIME = TIME + DT

     ! vmax= max( maxval(alpha), umax+cmax)
     ! !print*, maxval(alpha), umax+vmax
     ! dt = Dx/vmax!(umax + cmax)! vmax !
     ! DT = CFL*dt
     ! DT=MIN(DT,TIMEOUT-TIME)
     ! IF (dt.LE.1.e-10) PRINT*, 'dt.LE.1.e-10, program stop'
     ! IF (dt.LE.1.e-10) EXIT
     ! dt2 = 0.5d0*dt
     ! lambda = dt/dx
     ! IF (MOD(it,1000)==0)  WRITE(6,*) " Dt = ", Dt, 'time=', time, 'lambda =', lambda
     ! TIME = TIME + DT

     CALL shema(var,mesh, FLUX, lambda, rostar, pstar, einstar, alpha)
     CALL lim_cond(var, mesh)
     CALL state(var%prim)
     IT = IT + 1
     !----------------------------------------------------------------------------------------------------
     IF (TIME2.LE.TIME) THEN
        CALL Ecriture_donnees(var,mesh)
        CALL PutonScreen()
     END IF
     !----------------------------------------------------------------------------------------------------
     IF (TIME2.LE.TIME) THEN
        PRINT*, 'EN', IT, 'ITERATIONS, ', ' TIME:', TIME
        TIME2 = TIME + period_time
     END IF

     !----------------------------------------------------------------------------------------------------
  ENDDO TIME_LOOP
  !! FIN BOUCLE SUR LE TEMPS


  CALL Ecriture_donnees(var,mesh)

  
  !CALL Exact_sol_sod(mesh, TIME)

   !CALL Exact_Riemann_solver(mesh, time)



  !----------------------------------------------------------------------------------------------------
  CALL CPU_TIME(T2_CPU)
  !----------------------------------------------------------------------------------------------------
  PRINT*, 'L EXECUTION DU PROGRAMME A PRIS', T2_CPU - T1_CPU
  PRINT*, 'EN', IT-1, 'ITERATIONS, ', ' TIME:', TIME
CONTAINS

SUBROUTINE PutonScreen()
    IMPLICIT NONE
    REAL (KIND = DP):: MinVp(Nv_Prim), MaxVp(Nv_Prim)
    INTEGER:: iv

    DO iv = 1, Nv_Prim
       MinVp(iv)     = MINVAL(var%Prim(:)%v(iv))
       MaxVp(iv)     = MAXVAL(var%Prim(:)%v(iv))
    ENDDO

    WRITE(argunit,*)
    WRITE(argunit,'(65(":"))') 
    WRITE(argunit,'("::  ",47("-"),12(" "),"::")')
    WRITE(argunit,'(":: | Time  = ", E10.3, 1x, &
         &   "  ||     kt =  ",I9,3x,"| ",10x,"::")' ) TIME, it
    WRITE(argunit,'(":: |    Dt = ", E10.3,1x, &
         &   "  ||    Cfl =  ",E10.3,"  | ",10x,"::")' ) Dt , CFL
    WRITE(argunit,'("::  ",47("-"),12(" "),"::")')
    WRITE(argunit,'("::",61(" "),"::")') 
    WRITE(argunit,'(":: ",12x, " ",2(3x,A12), 11x," ::")') "Minimum","Maximum"

    DO iv = 1, Nv_Prim
       WRITE(argunit,'("::   ",A12, " => ",  2(3x,E12.5), 11x," ::" )') &
            & TRIM(NamesOfPrimVar(iv)), MinVp(iv) , MaxVp(iv)
    END DO

    WRITE(argunit,'("::   ",A12, " => ",  2(3x,E12.5), 11x," ::" )')  
    WRITE(argunit,'(65(":"))') 

  END SUBROUTINE PutonScreen
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE Exact_sol_sod(mesh, time)
    IMPLICIT NONE
    !.solves the exact sod problem
    INTEGER :: iter,i
    TYPE(maillage), INTENT(in):: mesh
    REAL(kind=dp), intent(in):: time
    REAL(kind=dp)::rho(1:nx),u(1:nx),p(1:nx)
    REAL(kind=dp)::rhol,pl,ul,rhor,pr,ur,xi,t, xl,xr
    REAL(kind=dp)::rho1,p1,u1,rho5,p5,u5,p40,p41,f0 !,f
    REAL(kind=dp)::f1,p4,error,z,c5,gm1,gp1,gmfac1,gmfac2,fact
    REAL(kind=dp)::u4,rho4,w,p3,u3,rho3,c1,c3,xsh,xcd,xft,xhd
    !dx

    !.define initial conditions
    !state at left of discontinuity
    rhol = 1.d0
    pl   = 1.0d0
    ul   = 0.0d0

    !.state at right of discontinuity
    rhor = 0.125d0
    pr   = 0.1d0
    ur   = 0.d0        
    !.location of discontinuity at t = 0
    xi = 0.5d0*Lx

    !.time at which solution is desired
    t = time

    !.spatial interval over which to compute solution
    xl = 0.0d0
    xr = 1.0d0
    IF (xr .LT. xl) THEN
       WRITE (6,*) 'xr must be greater than xl'
       STOP
    ENDIF

    !..begin solution
    IF (pl .GT. pr) THEN
       rho1 = rhol
       p1   = pl
       u1   = ul
       rho5 = rhor
       p5   = pr
       u5   = ur
    ELSE
       rho1 = rhor
       p1   = pr
       u1   = ur
       rho5 = rhol
       p5   = pl
       u5   = ul
    ENDIF

    !.solve for post-shock pressure by secant method
    !initial guesses

    p40 = p1
    p41 = p5
    f0  = fex(p40, p1, p5, rho1, rho5)

    DO iter = 1, iterfinal
       f1 = fex(p41, p1, p5, rho1, rho5)
       IF (f1 .EQ. f0) go to 10

       p4 = p41 - (p41 - p40)*f1/(f1 - f0)

       error = dabs(p4 - p41)/p41

       IF (error .LT. eps) go to 10

       p40 = p41
       p41 = p4
       f0  = f1
    ENDDO
    WRITE (6,*) 'iteration failed to converge'
    STOP 'abnormal termination'

10  CONTINUE


    !.compute post-shock density and velocity
    z  = (p4/p5 - 1.d0)
    c5 = dsqrt(gamma*p5/rho5)

    gm1 = gamma - 1.d0
    gp1 = gamma + 1.d0
    gmfac1 = 0.5d0*gm1/gamma
    gmfac2 = 0.5d0*gp1/gamma

    fact = dsqrt(1.d0 + gmfac2*z)

    u4 = c5*z/(gamma*fact)
    rho4 = rho5*(1.d0+ gmfac2*z)/(1.d0 + gmfac1*z)

    !.shock speed
    w = c5*fact
    !compute values at foot of rarefaction
    p3 = p4
    u3 = u4
    rho3 = rho1*(p3/p1)**(1.d0/gamma)

    !.compute positions of waves
    IF (pl .GT. pr) THEN
       c1 = dsqrt(gamma * p1/rho1)
       c3 = dsqrt(gamma * p3/rho3)

       xsh = xi + w*t
       xcd = xi + u3*t
       xft = xi + (u3 - c3)*t
       xhd = xi - c1*t

       !.and do say what we found
       WRITE (6, 500)
       WRITE (6, 501) rho1, p1, u1
       WRITE (6, 502) 
       WRITE (6, 503) rho3, p3, u3
       WRITE (6, 504) rho4, p4, u4
       WRITE (6, 505) rho5, p5, u5

       WRITE (6, 506) xhd
       WRITE (6, 507) xft
       WRITE (6, 508) xcd
       WRITE (6, 509) xsh

500    FORMAT (// 2x, 'Region', 4x, 'Density', 8x, 'Pressure', 8x, 'Velocity')
501    FORMAT (5x, '1', 3(2x,1pe14.7))
502    FORMAT (5x, '2' ,20x, 'RAREFACTION')
503    FORMAT (5x, '3', 3(2x,1pe14.7))
504    FORMAT (5x, '4', 3(2x,1pe14.7))
505    FORMAT (5x, '5', 3(2x,1pe14.7)//)

506    FORMAT (2x, 'Head Of Rarefaction    x = ', 1pe14.7)
507    FORMAT (2x, 'Foot Of Rarefaction    x = ', 1pe14.7)
508    FORMAT (2x, 'Contact Discontinuity  x = ', 1pe14.7)
509    FORMAT (2x, 'Shock                  x = ', 1pe14.7//)

       !compute solution as a function of position

       DO i = 1, nx
          IF (mesh%x(i) .LT. xhd) THEN
             rho(i) = rho1
             p(i)   = p1
             u(i)   = u1
          ELSE IF (mesh%x(i) .LT. xft) THEN
             u(i)   = 2.d0/ gp1 * (c1 + (mesh%x(i) - xi)/t)
             fact   = 1.d0- 0.5d0 * gm1 * u(i)/c1
             rho(i) = rho1 * fact ** (2.d0/gm1)
             p(i)   = p1 * fact ** (2.d0*gamma/gm1)
          ELSE IF (mesh%x(i) .LT. xcd) THEN
             rho(i) = rho3
             p(i)   = p3
             u(i)   = u3
          ELSE IF (mesh%x(i) .LT. xsh) THEN
             rho(i) = rho4
             p(i)   = p4
             u(i)   = u4
          ELSE
             rho(i) = rho5
             p(i)   = p5
             u(i)   = u5
          ENDIF
       ENDDO
    ENDIF

    !.if pr > pl, reverse solution
    IF (pr .GT. pl) THEN
       c1 = dsqrt(gamma*p1/rho1)
       c3 = dsqrt(gamma*p3/rho3)

       xsh = xi - w*t
       xcd = xi - u3*t
       xft = xi - (u3 - c3)*t
       xhd = xi + c1*t

       !.and do say what we found
       WRITE (6, 500)
       WRITE (6, 501) rho5, p5, u5
       WRITE (6, 602) rho4, p4, u4
       WRITE (6, 503) rho3, p3, u3
       WRITE (6, 604) 
       WRITE (6, 505) rho1, p1, u1

       WRITE (6, 609) xsh
       WRITE (6, 508) xcd
       WRITE (6, 507) xft
       WRITE (6, 606) xhd

602    FORMAT (5x, '2', 3(2x,1pe14.7))
604    FORMAT (5x, '4' ,20x, 'RAREFACTION')
606    FORMAT (2x, 'Head Of Rarefaction    x = ', 1pe14.7//)
609    FORMAT (2x, 'Shock                  x = ', 1pe14.7)

       DO i = 1, nx
          IF (mesh%x(i) .LT. xsh) THEN
             rho(i) = rho5
             p(i)   = p5
             u(i)   = -u5
          ELSE IF (mesh%x(i) .LT. xcd) THEN
             rho(i) = rho4
             p(i)   = p4
             u(i)   = -u4
          ELSE IF (mesh%x(i) .LT. xft) THEN
             rho(i) = rho3
             p(i)   = p3
             u(i)   = -u3
          ELSE IF (mesh%x(i) .LT. xhd) THEN
             u(i)   = -2.d0 /gp1*(c1 + (xi - mesh%x(i))/t)
             fact   = 1.d0 + 0.5d0*gm1*u(i)/c1
             rho(i) = rho1 * fact ** (2.d0/gm1)
             p(i)   = p1*fact**(2.d0*gamma/gm1)
          ELSE
             rho(i) = rho1
             p(i)   = p1
             u(i)   = -u1
          ENDIF
       ENDDO
    ENDIF
DO i = 1, nx
  if (rho(i) .le. eps) then 
  rho(i) = 0.d0
  endif
  if (u(i) .le. eps) then 
  u(i) = 0.d0
  endif
  if (p(i) .le. eps) then 
  p(i) = 0.d0
  endif
  
enddo
  

    OPEN (unit=1, file = './resu/TESTsod.out',status='unknown')
    WRITE (1, 1000)
1000 FORMAT (10x, 'x', 12x, 'density', 8x, 'pressure', 8x, 'velocity', 8x, 'internal energy'/)
1001 FORMAT (5(2x, 1pe14.7))
    DO i = 1, nx-1
       WRITE (1,1001) 0.5d0*(mesh%x(i) + mesh%x(i+1)), rho(i), 0.5d0*(u(i)+ u(i+1)),p(i),  p(i)/(rho(i)*(gamma - 1.d0))
    ENDDO
    CLOSE (1)
  END SUBROUTINE Exact_sol_sod
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  SUBROUTINE Exact_Riemann_solver(mesh, time)
    IMPLICIT NONE
    !.solves the exact sod problem
    INTEGER :: ix
    type(maillage), intent(in)::mesh 
    REAL(kind=dp), intent(in):: time
    REAL(kind=dp)::rhos,us,ps
    REAL(kind=dp)::rhol,pl,ul,rhor,pr,ur,diaph
    REAL(kind=dp):: G1, G2, G3,G4,G5,G6,G7,G8
    REAL(kind=dp):: cl, cr, pstar, ustar, xpos, s
    !.define initial conditions
    !state at left of discontinuity
  
    rhol = 1.d0
    ul   = -2.0d0
    pl   = 0.4d0

    !.state at right of discontinuity
    rhor = 1.0d0
    ur   = 2.d0 
    pr   = 0.4d0       
    !.location of discontinuity at t = 0
    diaph = 0.5d0*Lx

    !..begin solution
    ! gamma related constants

    G1 = (gamma - 1.d0)/(2.d0*gamma)
    G2 = (gamma + 1.d0)/(2.d0*gamma)
    G3 = (2.d0*gamma)/(gamma - 1.d0)
    G4 = 2.d0/(gamma - 1.d0)
    G5 = 2.d0/(gamma + 1.d0)
    G6 = (gamma - 1.d0)/(gamma + 1.d0)
    G7 = (gamma - 1.d0)/2.d0
    G8 = gamma - 1.d0

    ! compute sound speeds
    Cl = dsqrt(gamma*pl/rhol)
    Cr = dsqrt(gamma*pr/rhor)

    ! pressure positivity condition 
     if (g4*(cl+cr) .le. (ur - ul)) then 
      print*, '***************vacuum is generated by data***************'
      print*, '***************program stop***************'
      STOP
     endif
! exact solution in star region for pressure and velocity 

  CALL PUSTAR(rhol,pl,ul,rhor,pr,ur,cl,cr,pstar,ustar)
  OPEN (unit=1, file = './resu/Strong_exact.out',status='unknown')
   ! WRITE (1, 1000)
!1000 FORMAT (10x, 'x', 12x, 'density', 8x, 'velocity',  8x, 'internal energy,' 8x , 'pressure'/)
1001 FORMAT (5(2x, 1pe14.7))
  ! complerte solution at time = timeout 
  DO  ix =  1, nx
    xpos = (mesh%x(ix) - 0.5d0)*dx
    s = (xpos - diaph)/time
    ! solution at point (x,t) = (xpos - diaph, timeout)

    call sample(pstar,ustar,s,rhol,pl,ul,rhor,pr,ur,cl,cr,rhos,us,ps)

    WRITE (1,1001) xpos, rhos, us,  ps/(rhos*g8), ps
  ENDDO 
    CLOSE (1)
   RETURN
  END SUBROUTINE Exact_Riemann_solver
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE PUSTAR(rhol,pl,ul,rhor,pr,ur,cl,cr,pstar,ustar)
IMPLICIT NONE
! purpose to compute the solution for pressure and velocity 
! in the star region 
    INTEGER :: iter
    REAL(kind=dp), intent(in):: rhol,pl,ul,rhor,pr,ur,cl,cr
    REAL(kind=dp), intent(out):: pstar, ustar
    REAL(kind=dp):: pold, udiff, fl, fr, fld, frd 
    REAL(kind=dp):: p, change

    ! guessed value pstar is computed

    CALL GUESSP(pstar, rhol,pl,ul,rhor,pr,ur,cl,cr)
    pold = pstar 
    udiff = ur - ul 
     print*, 'iteration number    change'
  do iter = 1, iterfinal
     call prefun(fl, fld, pold, rhol, pl, cl)
     call prefun(fr, frd, pold, rhor, pr, cr)

     p = pold - (fl + fr + udiff)/(fld+ frd)
     change = 2.d0*dabs((p - pold)/(p + pold))
     print*, iter, change 
     if (change .le. eps) goto 20
     if (p .lt. 0.d0) p = eps
     pold = p
   ! print*, 'divergence in Newton-Raphson method'
    
  enddo 

    20 CONTINUE

    ! compute velocity in star region 
    ustar = 0.5d0*(ul + ur + fr - fl)
    pstar = p 
    print*, 'pressure   velocity in star region'
    print*, pstar, ustar 

    RETURN
END SUBROUTINE PUSTAR
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE GUESSP(pg, rhol,pl,ul,rhor,pr,ur,cl,cr)
IMPLICIT NONE
REAL(kind=dp), intent(in):: rhol,pl,ul,rhor,pr,ur, cl, cr
REAL(kind=dp), intent(out):: pg
REAL(kind=dp):: G1, G2, G3,G4,G5,G6,G7,G8
REAL(kind=dp):: pmax, pmin, quser, cup , ppv, qmax, um 
REAL(kind=dp):: PQ, PTL, PTR, gel, ger 
! purpose : to provide a guess value for pressure pg in the star region.

    G1 = (gamma - 1.d0)/(2.d0*gamma)
    G2 = (gamma + 1.d0)/(2.d0*gamma)
    G3 = (2.d0*gamma)/(gamma - 1.d0)
    G4 = 2.d0/(gamma - 1.d0)
    G5 = 2.d0/(gamma + 1.d0)
    G6 = (gamma - 1.d0)/(gamma + 1.d0)
    G7 = (gamma - 1.d0)/2.d0
    G8 = gamma - 1.d0
    quser = 2.d0
    ! compute guess pressure from pvrs riemann solver

    cup = 0.25d0*(rhol + rhor)*(cl + cr)
    ppv = 0.5d0*(pl + pr) + 0.5d0*(ul - ur)*cup
    ppv = Max(0.d0, ppv) 
    pmin= min(pl, pr)
    pmax= max(pl, pr)
    qmax= pmax/pmin
    
    if (qmax .le. quser .and. (pmin .le. ppv .and. ppv .le. pmax)) then 
      !select pvrs riemann solver 
      pg = ppv 
    else
      if (ppv .lt. pmin) then 
       !select two-rarefaction riemann solver 

       PQ = (pl/pr)**G1
       UM = (PQ*UL/cl + ur/cr +G4*(PQ - 1.d0))/(PQ/cl + 1.d0/cr)
       PTL = 1.d0 + G7*(ul - um)/cl
       PTR = 1.d0 + G7*(um - ur)/cr
       PG = 0.5d0*(pl*ptl**g3 + pr*ptr**g3)
      else
      ! select two-shock riemann solver with pvrs as estimate 
       gel = sqrt((g5/rhol)/(g6*pl +ppv))
       ger = sqrt((g5/rhor)/(g6*pr +ppv))
       PG = (gel*pl + ger*pr - (ur - ul))/(gel + ger)
      endif
    endif
RETURN
END SUBROUTINE GUESSP
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE  prefun(f, fd, p, rhok, pk, ck)
IMPLICIT NONE
REAL(kind=dp), INTENT(out)::f, fd 
REAL(kind=dp), INTENT(in)::p, rhok, pk, ck
REAL(kind=dp):: G1, G2, G3,G4,G5,G6,G7,G8
REAL(kind=dp):: prat , ak, bk , qrt  
! purpose: to evaluate the pressure functions FL and FR in exact riemann solver 
    G1 = (gamma - 1.d0)/(2.d0*gamma)
    G2 = (gamma + 1.d0)/(2.d0*gamma)
    G3 = (2.d0*gamma)/(gamma - 1.d0)
    G4 = 2.d0/(gamma - 1.d0)
    G5 = 2.d0/(gamma + 1.d0)
    G6 = (gamma - 1.d0)/(gamma + 1.d0)
    G7 = (gamma - 1.d0)/2.d0
    G8 = gamma - 1.d0

    if (p .le. pk) then 
     ! rarefaction wave 
       prat = p/pk 
       f = g4*ck*(prat**g1 - 1.d0)
       fd = (1.d0/(rhok*ck))*prat**(-g2) 
    else
      ! shock wave 
      ak = g5/rhok
      bk = g6*pk 
      qrt = sqrt(ak/(bk + p))
      f = (p - pk)*qrt 
      fd = (1.d0 - 0.5d0*(p - pk)/(bk + p))*qrt
    endif

RETURN
END SUBROUTINE prefun
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE sample(pm,um,s,rhol,pl,ul,rhor,pr,ur,cl,cr,rho,u,p)
IMPLICIT NONE
REAL(kind=dp),INTENT(out):: rho, u, p
REAL(kind=dp),INTENT(in)::  pm, um, s
REAL(kind=dp),INTENT(in)::  rhol,pl,ul,rhor,pr,ur,cl,cr
REAL(kind=dp):: G1, G2, G3,G4,G5,G6,G7,G8
REAL(kind=dp):: shl, cml, stl , c, pml , sl
REAL(kind=dp):: pmr, sr, shr , cmr, str 

! purpose: to sample the solution throughout the wave pattern. Pressure PM 
! and velocity UM in the star region are known. Sampling is performed in terms of the speed s = x/t. sampled 
! values are rho, u, p 
    G1 = (gamma - 1.d0)/(2.d0*gamma)
    G2 = (gamma + 1.d0)/(2.d0*gamma)
    G3 = (2.d0*gamma)/(gamma - 1.d0)
    G4 = 2.d0/(gamma - 1.d0)
    G5 = 2.d0/(gamma + 1.d0)
    G6 = (gamma - 1.d0)/(gamma + 1.d0)
    G7 = (gamma - 1.d0)/2.d0
    G8 = gamma - 1.d0

if (s .le. um) then
! sampling point lies to the left of the contact  discontinuity

if (pm .le.  pl) then 
  ! left rarefaction 
  shl  =  ul - cl 

  if (s .le. shl ) then 
    ! sampled point is left data state 
    rho = rhol
    u  = ul
    p = pl
  else
   cml = cl*(pm/pl)**g1
   stl = um - cml 
   if (s .gt. stl) then 
    ! sampled point is star left state 
    rho = rhol*(pm/pl)**(1.d0/gamma)
    u = um 
    p = pm
  else 
    ! sampled point is inside left fan 
    u = g5*(cl + g7*ul + s)
    c = g5*(cl + g7*(ul - s))
    rho = rhol*(c/cl)**g4 
    p = pl*(c/cl)**g3 
   endif
  endif 
else
  ! left shock 
  pml  = pm/pl 
  sl = ul - cl*sqrt(g2*pml + g1)

  if (s .le. sl) then
   ! sampled point is left data state 
   rho = rhol
   u  = ul
   p = pl
  else 
    ! sampled point is star left state 
    rho = rhol*(pml + g6)/(pml*g6 + 1.d0)
    u = um
    p = pm 
  endif 
 endif  

else 
 ! sampling point lies to the right of the contact discontinuity 
 if (pm .gt. pr) then

 ! right shock 
 pmr = pm/pr 
 sr = ur +  cr*sqrt(g2*pmr + g1) 

 if (s .ge. sr)  then 
  ! sampled point is right data state 
  rho = rhor
   u  = ur
   p = pr
 else 
  ! sampled point is star right state 
  rho = rhor*(pmr + g6)/(pmr*g6 + 1.d0)
  u = um
  p = pm 
 endif 
else 
! right rarefaction 

shr = ur + cr 

            if (s .ge. shr) then 
              ! right data state 
               rho = rhor
               u  = ur
               p = pr 
             else 
              cmr = cr*(pm/pr)**g1 
              str = um + cmr 
                if (s .le. str) then
                  ! star right state 
                  rho = rhor*(pm/pr)**(1.d0/gamma)
                  u  = um 
                  p = pm  
                else 
                   ! sampled point is inside left fan 
                   u = g5*(-cr + g7*ur + s)
                   c = g5*(cr - g7*(ur - s))
                   rho = rhor*(c/cr)**g4
                   p = pr*(c/cr)**g3
                endif
            endif 

endif
endif 

RETURN
END SUBROUTINE sample
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END PROGRAM main
