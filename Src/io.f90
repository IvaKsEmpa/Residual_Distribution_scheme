MODULE InOutputs
  USE GlobalParam
  USE precisions
  IMPLICIT NONE
CONTAINS
  SUBROUTINE LECTURE_DONNEES()
    IMPLICIT NONE
    OPEN(UNIT=21, FILE = './Src/data.inp', STATUS = 'OLD')
    READ(21,*) cond_lim, solver !COND LIM :1 box/2 absorbtion/3 batteur/4 jump;
    READ(21,*) order
    READ(21,*) Nx    ! NUMBER OF CELLS
    READ(21,*) Lx 
    READ(21,*) TIMEOUT         ! OUTPUT TIME
    READ(21,*) iterfinal       ! Iteration final
    READ(21,*) CFL          
    READ(21,*) ImpreE, ImpreF, period_time ! POUR IMPRIMER LES FICHIER SUR ECRAN ET DANS LE FICHIER
    CLOSE(21)
    RETURN
  END SUBROUTINE LECTURE_DONNEES
  !----------------------------------------------------------------------------------------------------
  SUBROUTINE Ecriture_donnees(var,mesh)
    IMPLICIT NONE
    TYPE(variable):: var
    TYPE(maillage):: mesh
    INTEGER :: jt, nu(2)
    REAL (KIND = DP):: X,v

    OPEN(MyUnit,FILE = './resu/blast3200.out')
    REWIND(MyUnit)
    DO jt=1, Mesh%nt

       nu=mesh%elem(:,jt)
       x=0.5*sum(mesh%x(nu))
       v=0.5*sum(var%v(nu))
       WRITE(MyUnit,10) x, var%prim(jt)%v(ro_pv),   v , var%prim(jt)%v(press_pv), var%prim(jt)%v(ein_pv)/var%prim(jt)%v(ro_pv)

    END DO
    WRITE(MyUnit,*)
    CLOSE(MyUnit)
10  FORMAT((5(f20.13,1X)))
  END SUBROUTINE Ecriture_donnees
END MODULE InOutputs
! ----------------------------------------------------------------------
