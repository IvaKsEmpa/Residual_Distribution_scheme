MODULE scheme
  USE precisions
  USE GlobalParam
  USE ModuleInterface
  IMPLICIT NONE  
CONTAINS

  !----------------------------------------------------------------------------------------------------
  SUBROUTINE INITIALISATION(var,mesh, dx)
    IMPLICIT NONE
    TYPE(variable), INTENT(inout):: var
    TYPE(maillage), INTENT(inout):: mesh
    INTEGER :: ix, i, jseg, jt
    REAL (KIND = DP), INTENT(in):: Dx
    REAL (KIND = DP) :: press,  ro , x0

    INTEGER:: nu(2)
    ALLOCATE(mesh%x(mesh%ns))
    x0 =0.d0
    DO ix = 1, mesh%ns
       mesh%x(ix) = x0 + REAL(ix-1)*dx
    ENDDO
    ALLOCATE(mesh%elem(2,mesh%nt))
    DO ix =1, mesh%nt
       ! element [i,i+1]
       mesh%elem(1,ix)=ix
       mesh%elem(2,ix)=ix+1
    ENDDO
    ALLOCATE(mesh%nutv(2,Mesh%nseg)) ! from each "segment"[cell interface], what are the
    ! neighboring elements
    DO jt=1, mesh%nt
       jseg=mesh%elem(1,jt)
       mesh%nutv(2,jseg)=jt ! element jt is on right of edge jseg
       jseg=mesh%elem(2,jt)
       mesh%nutv(1,jseg)=jt! element jt is on left of edge jseg+1
    ENDDO
    mesh%nutv(1,1)=0
    mesh%nutv(2,mesh%nseg)=0
    ALLOCATE(mesh%edge(2,mesh%nt)) ! what are the edges of element jt
    DO jt=1, mesh%nt
       mesh%edge(1,jt)=mesh%elem(1,jt) ! in 1D, there is identification between edges and vertices
       mesh%edge(2,jt)=mesh%elem(2,jt)
    ENDDO

    DO ix = 1, mesh%nt
       !----------------------------------------------------------------------------------------------------
       nu=mesh%elem(:,ix)
       DO i=1,2

          IF (COND_LIM == 1) THEN      ! SOD 
             IF ( mesh%x(nu(i)) .LE. 0.5d0*Lx) THEN  
                var%prim(ix)%v(ro_pv) =1.d0  
                var%prim(ix)%v(press_pv)  = 1.d0 
             ELSE
                var%prim(ix)%v(ro_pv) =  0.125d0  
                var%prim(ix)%v(press_pv)  = 0.10d0
             ENDIF
             var%v(nu(i))  = 0.d0
             !----------------------------------------------------------------------------------------------------
             press = var%prim(ix)%v(press_pv); ro  = var%prim(ix)%v(ro_pv)
             var%prim(ix)%v(sound_pv) = SoundPolGas(ro, press)
             var%prim(ix)%v(ein_pv)   = InternalEn(ro, press)
          ELSEIF (COND_LIM == 2) THEN      ! 123-problem
             IF ( mesh%x(nu(i)) .LE. 0.5d0*(Lx+x0)) THEN  
  
                var%v(nu(i))  = -2.d0
             ELSE
                  var%v(nu(i))  = 2.d0
             ENDIF
               var%prim(ix)%v(ro_pv) =1.d0  
                var%prim(ix)%v(press_pv)  = 0.4d0 
             !----------------------------------------------------------------------------------------------------
             press = var%prim(ix)%v(press_pv); ro  = var%prim(ix)%v(ro_pv)
             var%prim(ix)%v(sound_pv) = SoundPolGas(ro, press)
             var%prim(ix)%v(ein_pv)   = InternalEn(ro, press)

          ELSEIF (COND_LIM == 3) THEN      ! strong shock
             IF ( mesh%x(nu(i)) .LE. 0.5d0*(Lx)) THEN  
  
                var%prim(ix)%v(press_pv)  = 1000.0d0
             ELSE
              var%prim(ix)%v(press_pv)  = 0.01d0
                 
             ENDIF
               var%prim(ix)%v(ro_pv) =1.d0  
               var%v(nu(i))  = 0.d0
             !----------------------------------------------------------------------------------------------------
             press = var%prim(ix)%v(press_pv); ro  = var%prim(ix)%v(ro_pv)
             var%prim(ix)%v(sound_pv) = SoundPolGas(ro, press)
             var%prim(ix)%v(ein_pv)   = InternalEn(ro, press)
          ELSEIF (COND_LIM == 4) THEN      ! Interaction of blast waves.
             IF ( mesh%x(nu(i)) .LE. 0.1d0) THEN  
  
                var%prim(ix)%v(press_pv)  = 1000.0d0
             ELSEIF ( mesh%x(nu(i)) .gt. 0.1d0 .and. mesh%x(nu(i)) .LE. 0.9d0) THEN
              var%prim(ix)%v(press_pv)  = 0.01d0
             ELSEIF ( mesh%x(nu(i)) .gt. 0.9d0 .and. mesh%x(nu(i)) .LE. 1.0d0) THEN
              var%prim(ix)%v(press_pv)  = 100.0d0
                 
             ENDIF
               var%prim(ix)%v(ro_pv) =1.d0  
               var%v(nu(i))  = 1.d0
             !----------------------------------------------------------------------------------------------------
             press = var%prim(ix)%v(press_pv); ro  = var%prim(ix)%v(ro_pv)
             var%prim(ix)%v(sound_pv) = SoundPolGas(ro, press)
             var%prim(ix)%v(ein_pv)   = InternalEn(ro, press)
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE INITIALISATION
  !----------------------------------------------------------------------------------------------------
  SUBROUTINE lim_cond(var,mesh)
    IMPLICIT NONE
    !INTEGER ::  ix , iv
    TYPE(variable), INTENT(inout) :: var
    TYPE(maillage), INTENT(in):: mesh

    IF (COND_LIM == 1) THEN  
      ! DO iv = 1, Nv_Prim
          var%prim(1)%v=var%prim(2)%v ! 1st and 2nd element same
          var%v(1)=var%v(2) ! 1st and second points same

          var%prim(mesh%nt)%v=var%prim(mesh%nt-1)%v
          var%v(mesh%ns)=var%v(mesh%ns-1)

       !END DO
    ENDIF
    RETURN
  END SUBROUTINE lim_cond
  !----------------------------------------------------------------------------------------------------
  SUBROUTINE shema(var, mesh, FLUX, lambda, rostar, pstar, einstar, alpha)
    IMPLICIT NONE
    INTEGER :: jt, ix
    TYPE(maillage), INTENT(in):: mesh
    TYPE(variable), INTENT(inout):: var
    REAL (KIND=DP), INTENT(in)::rostar(:),pstar(:), einstar(:),alpha(:), flux(:,:)
    !REAL(dp),DIMENSION(:), INTENT(in):: pente
    REAL(dp), INTENT(in)::  lambda

    REAL(Kind=dp)::rotilde(mesh%nseg) !, rotilde_0(mesh%nseg)
    REAL(dp):: vit(2), restot, ps(2), fl(3,mesh%nt), e(2), res(2),&
         & resE, ro, ro_n, rot(2) !,rot_n(2)
    REAL(dp):: res_u(2,mesh%nt), v_spe(2), usq, v(2), v_n(2), vbar, zl, zr
    INTEGER:: nu(2),jseg, l, ne(2)
    REAL(Kind=dp)::correctionU,correctionE, beta
    !--------------------------------------------------------------------------

    var%res_v=0.
    DO jt=1, mesh%nt
      ne=mesh%edge(:,jt)
      DO l=1,3
        fl(l,jt)=flux(l,ne(2))-flux(l,ne(1))
      ENDDO
    ENDDO
    fl(:,1)=0.; fl(:, mesh%nt)=0.
    !--------------------------------------------------------------------------
    ! update density
    var%prim_n%v(ro_pv)  = var%Prim(:)%v(ro_pv) - lambda*fl(1,:) 
  DO jt=1, mesh%nt
   if ( var%prim_n(jt)%v(ro_pv) .ne.  var%prim_n(jt)%v(ro_pv)) then 
    print*, 'var%prim_n%v(ro_pv)', var%prim_n%v(ro_pv)
    stop
   endif 
 ENDDO
    ! update velocity
    DO jseg=2,mesh%nseg-1
      ne=mesh%nutv(:,jseg)
      rotilde(jseg) = 0.5*SUM( var%prim_n(ne)%v(ro_pv) )

    ENDDO
    rotilde(mesh%nseg) = rotilde(mesh%nseg-1)
    rotilde(1)    = rotilde(2)

    DO jt=1, mesh%nt

       nu=mesh%elem(:,jt)
       ne=mesh%edge(:,jt)
       
       ro=var%prim(jt)%v(ro_pv)
       ro_n=var%prim_n(jt)%v(ro_pv)

       beta= ro_n*MAXVAL(alpha(ne))  ! to have the same dimension, we myltiply by \ro here 
       vit=var%v(nu)
       ps=pstar(ne)
       vbar=0.5*SUM(vit)

       

       Restot= ro*vbar*(vit(2)-vit(1))+ ps(2)-ps(1) ! fl(2,jt) !
       Res=restot*0.5+ beta*( vit-vbar)

       correctionU=( fl(2,jt)- SUM(res)-vbar*fl(1,jt))

       Res= Res + 0.5*correctionU 
       res_u(:,jt)=res ! stor residual on velocity for further use

       var%res_v(nu)=Var%res_v(nu)+res
       !-------------------------------------------------------------------------------   
    ENDDO
    ! here we update the velocity
    var%res_v(1:2)=0.; var%res_v(Mesh%ns-1:Mesh%ns)=0.

! DO jt= 1, mesh%nt
!        !----------------------------------------------------------------------------------------------------
!        nu=mesh%elem(:,jt)
!  DO ix=1,2
     var%v_n=var%v-lambda*var%res_v/rotilde
!    ! var%v_n=var%v-lambda*fl(2,:)/rotilde
!  ENDDO
! ENDDO 
     
  DO jt=1, 2

   if ( var%v_n(jt) .ne.  var%v_n(jt)) then 
    print*, 'var%v_n(jt)', var%v_n(jt)
   endif 
ENDDO

    ! update energy

    DO jt=2, mesh%nt-1
       nu=mesh%elem(:,jt)
       ne=mesh%edge(:,jt)

       beta=maxval(alpha(ne))

       res=res_u(:,jt)
       v=var%v(nu)
       v_n=var%v_n(nu)
       e=einstar(ne)
       ps=pstar(ne)

       vbar=SUM(v)*0.5
       ro=var%prim(jt)%v(ro_pv)
       ro_n=var%prim_n(jt)%v(ro_pv)
       rot=rotilde(nu)

       v_spe= 0.5*(ro+ro_n)/rot * 0.5* (v+v_n)
       usq=0.25*SUM(v*v+v_n*v_n) !0.5*SUM(v*v+v_n*v_n)
       !
      ! ! sol 1
         ! ResE = fl(3,jt)
       !! sol 2
        zl = var%prim(jt)%v(ro_pv)*SoundPolGas(var%prim(jt)%v(ro_pv), var%prim(jt)%v(press_pv))
        zr = var%prim(jt+1)%v(ro_pv)*SoundPolGas(var%prim(jt+1)%v(ro_pv), var%prim(jt+1)%v(press_pv))
        ! ResE= vbar* (e(2)-e(1))+ ( var%prim(jt)%v(ein_pv)+var%prim(jt)%v(press_pv) )* (v(2)-v(1))!+beta*( var%prim(jt)%v(ein_pv)-0.5*sum(e))
        ResE= vbar* (e(2)-e(1))+(0.5D0*sum(e+ps)-zl*zr/(zl + zr)*(v(2)-v(1)))*(v(2)-v(1))!+beta*( var%prim(jt)%v(ein_pv)-0.5*sum(e))
        !!!ResE= vbar* (e(2)-e(1))+ 0.5*sum(e+pstar(ne))* (v(2)-v(1))
        correctionE=fl(3,jt)-ResE - ( 0.5*usq*fl(1,jt)+ SUM(v_spe*res) )
        ResE= ResE+correctionE

      ! ! sol 1
            !  var%Prim_n(jt)%v(ein_pv) = Var%prim(jt)%v(ein_pv)  - lambda*ResE- 0.25* sum(ro_n*v_n*v_n-ro*v*v)
       !! sol 2
       var%Prim_n(jt)%v(ein_pv) = Var%prim(jt)%v(ein_pv)  - lambda*ResE
      if (var%Prim_n(jt)%v(ein_pv) .ne. var%Prim_n(jt)%v(ein_pv)) then 
       print*, 'var%Prim_n(jt)%v(ein_pv)', var%Prim_n(jt)%v(ein_pv)
       stop
      endif 
    
    ENDDO

    !var%Prim_n(:)%v(ein_pv) = Var%prim(jt)%v(ein_pv)  - lambda*fl(3,:) 
    var%prim_n(1)       = var%prim_n(2)
    var%prim_n(Mesh%nt) = var%prim_n(Mesh%nt-1)
    ! update
    DO l=1, Nv_prim
       var%prim(:)%v(l)=var%prim_n%v(l)
    ENDDO
    var%v=var%v_n
    CALL state(var%prim)
  END SUBROUTINE shema
  !----------------------------------------------------------------------------------------------------
  SUBROUTINE shema2(var, mesh, FLUX, lambda, rostar, pstar, einstar, alpha)
    IMPLICIT NONE
    INTEGER :: jt, ix
    TYPE(maillage), INTENT(in):: mesh
    TYPE(variable), INTENT(inout):: var
    REAL (KIND=DP), INTENT(in)::rostar(:),pstar(:), einstar(:),alpha(:), flux(:,:)
    !REAL(dp),DIMENSION(:), INTENT(in):: pente
    REAL(dp), INTENT(in)::  lambda

    REAL(Kind=dp)::rotilde(mesh%nseg) !, rotilde_0(mesh%nseg)
    REAL(dp):: vit(2), restot, ps(2), fl(3,mesh%nt), e(2), res(2),&
         & resE, ro, ro_n, rot(2) !,rot_n(2)
    REAL(dp):: res_u(2,mesh%nt), v_spe(2), usq, v(2), v_n(2), vbar, zl, zr
    INTEGER:: nu(2),jseg, l, ne(2)
    REAL(Kind=dp)::correctionU,correctionE, beta
    !--------------------------------------------------------------------------
! c'est pareille pour l'instant que shema 1 (ordre 1). il faut monter l'ordre ici
    var%res_v=0.
    DO jt=1, mesh%nt
      ne=mesh%edge(:,jt)
      DO l=1,3
        fl(l,jt)=flux(l,ne(2))-flux(l,ne(1))
      ENDDO
    ENDDO
    fl(:,1)=0.; fl(:, mesh%nt)=0.
    !--------------------------------------------------------------------------
    ! update density
    var%prim_n%v(ro_pv)  = var%Prim(:)%v(ro_pv) - lambda*fl(1,:) 
  DO jt=1, mesh%nt
   if ( var%prim_n(jt)%v(ro_pv) .ne.  var%prim_n(jt)%v(ro_pv)) then 
    print*, 'var%prim_n%v(ro_pv)', var%prim_n%v(ro_pv)
    stop
   endif 
 ENDDO
    ! update velocity
    DO jseg=2,mesh%nseg-1
      ne=mesh%nutv(:,jseg)
      rotilde(jseg) = 0.5*SUM( var%prim_n(ne)%v(ro_pv) )

    ENDDO
    rotilde(mesh%nseg) = rotilde(mesh%nseg-1)
    rotilde(1)    = rotilde(2)

    DO jt=1, mesh%nt

       nu=mesh%elem(:,jt)
       ne=mesh%edge(:,jt)
       
       ro=var%prim(jt)%v(ro_pv)
       ro_n=var%prim_n(jt)%v(ro_pv)

       beta= ro_n*MAXVAL(alpha(ne))  ! to have the same dimension, we myltiply by \ro here 
       vit=var%v(nu)
       ps=pstar(ne)
       vbar=0.5*SUM(vit)

       

       Restot= ro*vbar*(vit(2)-vit(1))+ ps(2)-ps(1) ! fl(2,jt) !
       Res=restot*0.5+ beta*( vit-vbar)

       correctionU=( fl(2,jt)- SUM(res)-vbar*fl(1,jt))

       Res= Res + 0.5*correctionU 
       res_u(:,jt)=res ! stor residual on velocity for further use

       var%res_v(nu)=Var%res_v(nu)+res
       !-------------------------------------------------------------------------------   
    ENDDO
    ! here we update the velocity
    var%res_v(1:2)=0.; var%res_v(Mesh%ns-1:Mesh%ns)=0.

! DO jt= 1, mesh%nt
!        !----------------------------------------------------------------------------------------------------
!        nu=mesh%elem(:,jt)
!  DO ix=1,2
     var%v_n=var%v-lambda*var%res_v/rotilde
!    ! var%v_n=var%v-lambda*fl(2,:)/rotilde
!  ENDDO
! ENDDO 
     
  DO jt=1, 2

   if ( var%v_n(jt) .ne.  var%v_n(jt)) then 
    print*, 'var%v_n(jt)', var%v_n(jt)
   endif 
ENDDO

    ! update energy

    DO jt=2, mesh%nt-1
       nu=mesh%elem(:,jt)
       ne=mesh%edge(:,jt)

       beta=maxval(alpha(ne))

       res=res_u(:,jt)
       v=var%v(nu)
       v_n=var%v_n(nu)
       e=einstar(ne)
       ps=pstar(ne)

       vbar=SUM(v)*0.5
       ro=var%prim(jt)%v(ro_pv)
       ro_n=var%prim_n(jt)%v(ro_pv)
       rot=rotilde(nu)

       v_spe= 0.5*(ro+ro_n)/rot * 0.5* (v+v_n)
       usq=0.25*SUM(v*v+v_n*v_n) !0.5*SUM(v*v+v_n*v_n)
       !
      ! ! sol 1
         ! ResE = fl(3,jt)
       !! sol 2
        zl = var%prim(jt)%v(ro_pv)*SoundPolGas(var%prim(jt)%v(ro_pv), var%prim(jt)%v(press_pv))
        zr = var%prim(jt+1)%v(ro_pv)*SoundPolGas(var%prim(jt+1)%v(ro_pv), var%prim(jt+1)%v(press_pv))
        ! ResE= vbar* (e(2)-e(1))+ ( var%prim(jt)%v(ein_pv)+var%prim(jt)%v(press_pv) )* (v(2)-v(1))!+beta*( var%prim(jt)%v(ein_pv)-0.5*sum(e))
        ResE= vbar* (e(2)-e(1))+(0.5D0*sum(e+ps)-zl*zr/(zl + zr)*(v(2)-v(1)))*(v(2)-v(1))!+beta*( var%prim(jt)%v(ein_pv)-0.5*sum(e))
        !!!ResE= vbar* (e(2)-e(1))+ 0.5*sum(e+pstar(ne))* (v(2)-v(1))
        correctionE=fl(3,jt)-ResE - ( 0.5*usq*fl(1,jt)+ SUM(v_spe*res) )
        ResE= ResE+correctionE

      ! ! sol 1
            !  var%Prim_n(jt)%v(ein_pv) = Var%prim(jt)%v(ein_pv)  - lambda*ResE- 0.25* sum(ro_n*v_n*v_n-ro*v*v)
       !! sol 2
       var%Prim_n(jt)%v(ein_pv) = Var%prim(jt)%v(ein_pv)  - lambda*ResE
      if (var%Prim_n(jt)%v(ein_pv) .ne. var%Prim_n(jt)%v(ein_pv)) then 
       print*, 'var%Prim_n(jt)%v(ein_pv)', var%Prim_n(jt)%v(ein_pv)
       stop
      endif 
    
    ENDDO

    !var%Prim_n(:)%v(ein_pv) = Var%prim(jt)%v(ein_pv)  - lambda*fl(3,:) 
    var%prim_n(1)       = var%prim_n(2)
    var%prim_n(Mesh%nt) = var%prim_n(Mesh%nt-1)
    ! update
    DO l=1, Nv_prim
       var%prim(:)%v(l)=var%prim_n%v(l)
    ENDDO
    var%v=var%v_n
    CALL state(var%prim)
  END SUBROUTINE shema2
  !-------------------------------------------------------------------------------
  SUBROUTINE HLLC(dx, var, mesh, FLUX,  pente, rostar, pstar, einstar, alpha)
    IMPLICIT NONE
  
    TYPE(maillage), INTENT(in):: mesh
    TYPE(variable), INTENT(in):: var
    REAL(dp), INTENT(in):: dx
    REAL(dp), INTENT(out):: FLUX(3,mesh%nseg)
    REAL (KIND=DP),INTENT(out):: rostar(Mesh%nseg), pstar(mesh%nseg), &
         & einstar(mesh%nseg), alpha(mesh%nseg)
    REAL(KIND=dp),INTENT(in):: pente(mesh%nseg)

    REAL(Kind=dp):: ul, ur, rol, ror, ml, mr, pressl, pressr
    REAL(Kind=dp):: einl, einr
    REAL(Kind=dp):: ustar,  sr, sl, cl, cr, EL, ER, Estar
    INTEGER:: jseg, jt1, jt2, is
    rostar = 0.d0; pstar = 0.d0;  einstar = 0.d0
    alpha = 0.d0

    DO jseg=2, mesh%nseg-1
       jt1=mesh%nutv(1,jseg); jt2=mesh%nutv(2,jseg); is=jseg

       ! solution of the Riemann problem
       ul=var%v(is)  - 0.5D0*dx*PENTE(is)
       ur=var%v(is)  + 0.5D0*dx*PENTE(is)

       rol=var%prim(jt1)%v(ro_pv);          ror=var%prim(jt2)%v(ro_pv)
       pressl=var%prim(jt1)%v(press_pv);    pressr=var%prim(jt2)%v(press_pv)
       einl = var%prim(jt1)%v(ein_pv) ;     einr = var%prim(jt2)%v(ein_pv)
       EL = TotalEn(rol, ul, einl); ER = TotalEn(ror, ur, einr) 
       !! sound speed
       cl=SoundPolGas(rol, pressl); cr=SoundPolGas(ror, pressr)
       ! davis
       sr=ur+cr; IF((ul+cl)> sr) sr = ul+cl
       Sl=ul-cl; IF ((ur-cr) < sl) sl = ur - cr
       alpha(jseg) =DMAX1( dabs(sr), dabs(sl))  ! max(abs(ur)+cr,abs(ul)+cl)!
       ml=rol*(ul-sl)
       mr=ror*(ur-sr)
       ustar = var%v(jseg)

       IF (ustar.GE.0.d0) THEN
          IF (sl.GE.0.d0) THEN
             rostar(jseg)=rol
             pstar(jseg) = pressl
             einstar(jseg) = einl
             FLUX(1,jseg) = rol*ul
             FLUX(2,jseg) = rol*ul**2 + pressl
             FLUX(3,jseg) = (El + pressl)*ul 
          ELSE 
             rostar(jseg)=ml/(ustar-sl)
             pstar(jseg)= ( (ur-ul)*mr*ml-mr*pressl+ml*pressr )/(ml-mr)
             Estar= ( El*(ul - sl)+pressl*ul-pstar(jseg)*ustar )/(ustar - sl)
             einstar(jseg)=  Estar - 0.5d0*rostar(jseg)*ustar**2
             FLUX(1,jseg) = rostar(jseg)*ustar
             FLUX(2,jseg) = rostar(jseg)*ustar**2 + pstar(jseg)
             FLUX(3,jseg) = (Estar + pstar(jseg))*ustar  
          ENDIF
       ELSE
          IF (sr.GE.0.d0) THEN
             rostar(jseg)=mr/(ustar-sr)
             pstar(jseg)= ((ur-ul)*mr*ml-mr*pressl+ml*pressr )/(ml-mr)
             Estar= ( Er*(ur - sr)+pressr*ur-pstar(jseg)*ustar )/(ustar - sr)
             einstar(jseg)=  Estar - 0.5d0*rostar(jseg)*ustar**2 

             FLUX(1,jseg) = rostar(jseg)*ustar
             FLUX(2,jseg) = rostar(jseg)*ustar**2 + pstar(jseg)
             FLUX(3,jseg) = (Estar + pstar(jseg))*ustar 
          ELSE
             rostar(jseg)=ror
             pstar(jseg) = pressr
             einstar(jseg) = einr
             FLUX(1,jseg) = ror*ur
             FLUX(2,jseg) = ror*ur**2 + pressr
             FLUX(3,jseg) = (Er + pressr)*ur
          ENDIF
       ENDIF

    ENDDO
    alpha(1) =  alpha(2)
    flux(:,1)=flux(:,2)
    alpha(mesh%nseg)=alpha(Mesh%nseg-1)
    flux(:,mesh%nseg)=flux(:,mesh%nseg-1)

  END SUBROUTINE HLLC

  !----------------------------------------------------------------------------------------------------

  SUBROUTINE HLL(dx, var, mesh, FLUX,  pente, rostar, pstar, einstar, alpha)
    IMPLICIT NONE
    INTEGER:: k, ix, i1, i2
    TYPE(maillage), INTENT(in):: mesh
    TYPE(variable), INTENT(in):: var
    REAL(dp), INTENT(in):: dx
    REAL(dp), INTENT(out):: FLUX(3,mesh%nseg)
    REAL (KIND=DP),INTENT(out):: rostar(Mesh%nseg), pstar(mesh%nseg), &
         & einstar(mesh%nseg), alpha(mesh%nseg)
    REAL(KIND=dp),INTENT(in):: pente(mesh%nseg)

    REAL(Kind=dp):: ul, ur, rol, ror, ml, mr, pressl, pressr
    REAL(Kind=dp):: lambda, einl, einr, FLL(3), FRR(3), ConsL(3), ConsR(3)
    REAL(Kind=dp):: ustar,  sr, sl, cl, cr, EL, ER, Estar
    INTEGER:: jseg, jt1, jt2, is
    rostar = 0.d0; pstar = 0.d0;  einstar = 0.d0
    alpha = 0.d0

    DO jseg=2, mesh%nseg-1
       jt1=mesh%nutv(1,jseg); jt2=mesh%nutv(2,jseg); is=jseg

       ! solution of the Riemann problem
       ul=var%v(is)  - 0.5D0*dx*PENTE(is)
       ur=var%v(is)  + 0.5D0*dx*PENTE(is)

       rol=var%prim(jt1)%v(ro_pv);          ror=var%prim(jt2)%v(ro_pv)
       pressl=var%prim(jt1)%v(press_pv);    pressr=var%prim(jt2)%v(press_pv)
       einl = var%prim(jt1)%v(ein_pv) ;     einr = var%prim(jt2)%v(ein_pv)
       EL = TotalEn(rol, ul, einl); ER = TotalEn(ror, ur, einr) 
       !! sound speed
       cl=SoundPolGas(rol, pressl); cr=SoundPolGas(ror, pressr)
       ! davis
       sr=ur+cr; IF((ul+cl)> sr) sr = ul+cl
       Sl=ul-cl; IF ((ur-cr) < sl) sl = ur - cr
       alpha(jseg) =DMAX1( dabs(sr), dabs(sl)) 
       ml=rol*(ul-sl)
       mr=ror*(ur-sr)
       ustar = var%v(jseg)

       !IF (ustar.GE.0.d0) THEN
          IF (sl.GE.0.d0) THEN
             rostar(jseg)=rol
             pstar(jseg) = pressl
             einstar(jseg) = einl
             FLUX(1,jseg) = rol*ul
             FLUX(2,jseg) = rol*ul**2 + pressl
             FLUX(3,jseg) = (El + pressl)*ul 
          ELSEIF (sl.le.0.d0 .and. sr.ge.0.d0) THEN 
             rostar(jseg)=(sr*ror - sl*rol + rol*ul - ror*ur)/(sr - sl) ! ml/(ustar-sl) !
             Estar=  (sr*Er - sl*El + (El+pressl)*ul - (Er+pressr)*ur)/(sr - sl) !( El*(ul - sl)+pressl*ul-pstar(jseg)*ustar )/(ustar - sl) !
             einstar(jseg)=  Estar - 0.5d0*rostar(jseg)*ustar**2
             pstar(jseg)=  ( (ur-ul)*mr*ml-mr*pressl+ml*pressr )/(ml-mr) ! Pressure(rostar(jseg),einstar(jseg)) !
            ! FLUX(1,jseg) = rostar(jseg)*ustar
            ! FLUX(2,jseg) = rostar(jseg)*ustar**2 + pstar(jseg)
            ! FLUX(3,jseg) = (Estar + pstar(jseg))*ustar  

             Frr(1) = ror*ur
             Frr(2) = ror*ur**2 + pressr
             Frr(3) = (Er + pressr)*ur

             FlL(1) = rol*ul
             FlL(2) = rol*ul**2 + pressl
             FlL(3) = (El + pressl)*ul 

             Consr(1) = ror
             Consr(2) = ror*ur
             Consr(3) = Er

             ConsL(1) = rol
             ConsL(2) = rol*ul
             ConsL(3) = El 

        FLUX(1,jseg) = (sr*Fll(1) - sl*Frr(1)  + sl*sr*(ConsR(1) - ConsL(1)) )/(sr - sl)
        FLUX(2,jseg) = (sr*Fll(2) - sl*Frr(2)  + sl*sr*(ConsR(2) - ConsL(2)) )/(sr - sl)
        FLUX(3,jseg) = (sr*Fll(3) - sl*Frr(3)  + sl*sr*(ConsR(3) - ConsL(3)) )/(sr - sl)
            
 !else 
       ELSEIF (sr.le.0.d0) THEN
             rostar(jseg)=ror
             pstar(jseg) = pressr
             einstar(jseg) = einr
             FLUX(1,jseg) = ror*ur
             FLUX(2,jseg) = ror*ur**2 + pressr
             FLUX(3,jseg) = (Er + pressr)*ur
          ENDIF
      ! ENDIF

    ENDDO
    alpha(1) =  alpha(2)
    flux(:,1)=flux(:,2)
    alpha(mesh%nseg)=alpha(Mesh%nseg-1)
    flux(:,mesh%nseg)=flux(:,mesh%nseg-1)

  END SUBROUTINE HLL
!--------------------------------------------------------------------------------------------------------------------------------------------------------------------
  SUBROUTINE PENTE_U( DX,Var, mesh, PENTE)
    IMPLICIT NONE
    INTEGER ::  jt1, jt2, iseg ! is,

    TYPE(maillage), INTENT(in):: mesh
    TYPE(variable), INTENT(in):: Var
    REAL (KIND = DP) :: DX, pent, pente11, pente12
    REAL (KIND = DP), INTENT(out) ::  PENTE(:)

    pente = 0.d0
    DO iseg=2, mesh%nseg-1 ! the first and last egde are special, boundaries
       jt1=mesh%nutv(1,iseg); jt2=mesh%nutv(2,iseg)
       PENTE11 = (Var%v(mesh%elem(2,jt1))-Var%v(mesh%elem(1,jt1)))/dx
       PENTE12=  (Var%v(mesh%elem(2,jt2))-Var%v(mesh%elem(1,jt2)))/dx
       pent = superbee(pente11,pente12)
       PENTE(iseg)=pent
    END DO
    PENTE(mesh%nseg) = PENTE( mesh%nseg-1) 
    PENTE(1) = PENTE(2)


  END SUBROUTINE PENTE_U
  !--------------------------------------------------------------------------------

  FUNCTION vleer(s1,s2) RESULT(slim)
    REAL(kind=dp) :: s1,s2, slim,prod   ! ss1,ss2,
    prod=s1*s2
    IF(prod>0.d0.AND.dabs(s1+s2)>1.d-6) THEN
       slim=2.d0*prod/(s1+s2)
    ELSE
       slim=0.d0
    ENDIF
  END FUNCTION vleer
  !--------------------------------------------------------------------------------
  FUNCTION minmod(s1,s2) RESULT(slim)
    REAL(kind=dp) :: s1,s2, slim,ss1,ss2

    IF(s1*s2 > 1.0E-08) THEN
       ss1=dabs(s1)
       ss2=dabs(s2)
       slim=dmin1(ss1,ss2)
       IF(s1 < 0.0)slim=-slim
    ELSE
       slim=0.0
    ENDIF
  END FUNCTION minmod
  !-------------------------------------------------------------------------------

  FUNCTION superbee(s1,s2) RESULT(slim)
    REAL(kind=dp) :: s1,s2, slim,ss1,ss2  
    IF(s1*s2.GT.0.d0) THEN
       ss1=dabs(s1)
       ss2=dabs(s2)
       slim=dmax1(dmin1(2.d0*ss1,ss2),dmin1(ss1,2.d0*ss2))
       IF(s1.LT.0.d0)slim=-slim
    ELSE
       slim=0.d0
    ENDIF
  END FUNCTION superbee

  !-------------------------------------------------------------------

  SUBROUTINE valbada(s1,s2,slim)
    IMPLICIT NONE
    REAL(kind=dp) :: s1,s2,slim  

    IF(s1.NE.0.d0.OR.s2.NE.0.d0) THEN
       slim=s1*s2*(s1+s2)/(s1*s1+s2*s2)
    ELSE
       slim=0.d0
    ENDIF
    RETURN
  END SUBROUTINE valbada
  !-------------------------------------------------------------------
  SUBROUTINE vmc(s1,s2,slim)
    IMPLICIT NONE
    REAL(kind=dp) :: s1,s2,slim,ss1,ss2   

    IF(s1*s2.GT.0.d0) THEN
       ss1=dabs(s1)
       ss2=dabs(s2)
       slim=dmin1(2.d0*ss1,2.d0*ss2,0.5d0*(ss1+ss2))
       IF(s1.LT.0.d0)slim=-slim
    ELSE
       slim=0.d0
    ENDIF
    RETURN
  END SUBROUTINE vmc
  !--------------------------------------------------------------------
  SUBROUTINE state(prim)
    IMPLICIT NONE
    TYPE(variable_ther), DIMENSION(:), INTENT(inout):: prim
    INTEGER:: nt, jt
    REAL(dp):: ro, p, ein,c
    nt=SIZE(prim)
    DO jt=1, nt
       ein=prim(jt)%v(ein_pv)
       ro=prim(jt)%v(ro_pv)
       p=pressure(ro, ein)
       c=soundPolGas(ro,p)
       prim(jt)%v(press_pv)=p
       prim(jt)%v(sound_pv)=c
    ENDDO
  END SUBROUTINE state
  ! ----------------------------------------------------------------------
END MODULE scheme
