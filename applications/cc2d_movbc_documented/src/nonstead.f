************************************************************************
* This file contains the main nonsteady loop to solve the NS equations.
* 
* In:
*  MFILE  - file handle to write output to
*  MSHOW  - level of output
*  ITFILM - the start number of NS iter. for output; 
*           is incremented in the loop
*  IFILEN - the start number for output with GMV,AVS,...; is 
*           incremented in the loop.
* 
* Return values can be found in the appropriate COMMON blocks.
* Before the iteration the caller should set all values in the 
* /NSCOUN/ Common block to 0. These values are incremented during
* the iteration to provide the caller with information it.
* 
* The routine initializes all necessary variables in the Common blocks
* by itself, using the calculated values in the initialization
* routine INIT1.
* That means: To start the same computation twice, only the start 
* solution vector on the finest level has to be the same 
* (e.g. zero + implementation of boundary values).
*
* Before this routine can be called, the following steps must be
* performed:
* - The memory management must be initialized (ZINIT)
* - The CC2D.DAT-parameters must be read into the COMMON-block
*   variables (RDOUT/RDDAT)
* - Any geometry information necessary for the computation must be
*   initialized
* - The solver structures must be initialized. This can be done
*   by INIT1 (for the first run) or ININSO (for all further runs)
* - The output channels must be initialized (INPWFS)
* - The time stepping must be initialized (INNSTT)
************************************************************************

      SUBROUTINE NONSTL (MFILE,MSHOW,ITFILM,IFILEN)

      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicelem.inc'
      
      INCLUDE 'cns.inc'
      INCLUDE 'cnsparfrac.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgtria.inc'
      INCLUDE 'cmgadr.inc'
      
      INCLUDE 'cinidat.inc'

C externals

C *** Parametrization of the domain
      DOUBLE PRECISION PARX,PARY,TMAX
      EXTERNAL PARX,PARY,TMAX
      
C *** Coefficient of exact solution
      DOUBLE PRECISION UE,PE,UEX,UEY
      EXTERNAL UE,PE,UEX,UEY
      
C *** definition of finite elements
      EXTERNAL E030,E031,EM30,EM31,E010
      
C parameters

      INTEGER MFILE,MSHOW,IFILEN,ITFILM
      
C local variables

      INTEGER ISETLV
      INTEGER LRHS,LTM0,LTM3,NNITR,ITNSR,NITNSL,LTMLO,IFRSTH,NMG
      INTEGER IUP,INDMAX,ITYP,ICUBER,IFPOST

      DOUBLE PRECISION TTTSUM,TTT0,TTT1,TIMEST,TTTOTH
      DOUBLE PRECISION TIMO31,TIMO32,TIMO11,TIMO12,TSTEPO,TIMNSH,THETAH
      DOUBLE PRECISION TSTEPN,RELU21,RELUM1,RELP21,RELPM1,DSXN,RELU20
      DOUBLE PRECISION RELUM0,RELP20,RELPM0,RELT1,EPSAD,CTPROJ,HTFACT
      DOUBLE PRECISION DRELTN,EX1,EX2,DFWVI,DAWVI,RELU2,RELP2,RELUM
      DOUBLE PRECISION RELPM,RELT

      LOGICAL BSTOP,BNL1,BNL2

      TIMEST=TIMENS
      CALL ZTIME(TTT0)
      TTTSUM = TTT0

C     Do we have a stationary or instationary simulation?

      IF (ISTAT.NE.0) THEN
      
C       Our simulation is instationary. In this case we have to allocate
C       some more vectors that contain information during the time steps.
      
        ISETLV=2
        ILEV=NLMAX
        CALL  SETLEV (ISETLV)

C       In case the right hand side is not homogenuous (IRHS <> 0 from the
C       DAT file), we have to calculate it - so we have to
C       allocate a RHS-vector - a local, not a COMMON-block variable.

        LRHS=0
        IF (IRHS.GE.1) THEN
          CALL ZNEW(NUP,1,LRHS  ,'DRHS  ')
          IF (IER.NE.0) GOTO 99998
          CALL LCP1(DWORK(KF1),DWORK(L(LRHS)),NUP)
        ENDIF

C       Allocate LTMP0 - a temporary vector that holds the solution
C       in the beginning of the current time step.

        CALL ZNEW(NUP,-1,LTM0  ,'DMT0  ')
        IF (IER.NE.0) GOTO 99998

C       In case that extrapolation replaces the nonlinear
C       iteration, we need additional auxiliary vectors.
C
C       Normally, the convective term in the Navier stokes equation
C                      u^(n+1) * grad(u^(n+1))
C       makes it necessary to perform nonlinear iteration with a
C       maximum of INLMAX iterations to calculate an approximation
C       to u^(n+1) from u^n for the next time step. This treatment
C       of the nonlinearity can be simplified by constant or linear
C       extrapolation in time. Using constant interpolation replaces
C       this term by
C                      u^n * grad(u^(n+1))
C       while linear extrapolation replaces it by
C                      [ 2 u^n + u^(n-1) ] * grad(u^(n+1))
C       where the latter term is 2nd order in time.
C       With both approaches, a nonlinear iteration is no more
C       necessary, thus speeding up the computation - with the
C       disadvantage of changing the equation, so obtaining a solution
C       to another equation than Navier-Stokes.
C
C       Extrapolation in time is active as soon as INLMIN=INLMAX.
C       If INLMIN=INLMAX=-1, the constant extrapolation is active,
C       while INLMIN=INLMAX=1 activates the linear extrapolation,
C       which is also indicated by the station flag ITEXL<>0.
C
C       For linear extrapolation we allocate a vector DTML,
C       which saves the old solution u^(n-1) of the previous time step.
C       LTML is a COMMON-block variable and updated in NSDEF
C       from u^(n-1) to u^n, while the solution vector in DU is
C       updated from u^n to u^(n+1).

        IF (ITEXL.NE.0) THEN
        
          CALL ZNEW(NUP,-1,LTML ,'DMTL  ')
          IF (IER.NE.0) GOTO 99998
          
C         If adaptive time-stepping technique 3 is used (prediction with
C         repetition and time-error handling) we need additional vectors. 
C         In that case, DTML takes a backup ???
          
          IF (ABS(IADTIM).EQ.3) THEN
          
            CALL ZNEW(NUP,-1,LTMLO,'DMTLO ')
            
            IF (IER.NE.0) GOTO 99998
            
C           Copy the current (starting) solution into DTML allocated 
C           above (not DTMLO as one would expect!).
C           That way the start vector can be restored.
C
C           In the very first step of the algorithm, if adaptive
C           time stepping technique 3 is used, the content of
C           DTML is copied to DTMLO. So by this we make sure there
C           is something to copy into DTMLO.

            CALL LCP1(DWORK(KU1),DWORK(L(LTML)),NUP)
            
          ENDIF ! (ABS(IADTIM).EQ.3)
          
        ENDIF ! (ITEXL.NE.0)

C       For adaptive time-stepping we need a vector DTM3 that
C       saves the result of the predictor-step.

        IF (IADTIM.NE.0) THEN
          CALL ZNEW(NUP,-1,LTM3  ,'DMT3  ')
          IF (IER.NE.0) GOTO 99998
        ENDIF
        
      ENDIF ! (ISTAT.NE.0)

C     To summarize the use of the vectors allocated above:
C
C     DU      - current solution vector u^n, updated to u^(n+1)
C               during NSDEF
C     DTM0    - holds a copy of u^n so that it can be restores;
C               must be restored for recalculation in case of too
C               big time-step and when switching from prediction
C               step to standard step-size.
C     DTML    - holds u^(n-1) for linear extrapolation;
C               COMMON-block variable, updated in NSDEF to u^n
C               when DU is updated from u^n to u^(n+1)
C     DTMLO   - Holds a copy of DTML in case that it has to be
C               restored because of recalculation. Only
C               if IADTIM=+-3.
C     DTM3    - Receives the calculated solution of the predictor
C               step in a time-dependent calculation

      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0

C=======================================================================
C *** Begin of nonsteady loop
C=======================================================================

C***********************************************************************
C *** NNITR = max number of repetitions for |IADAPT|=3

      NNITR=IREPIT

C***********************************************************************

C     Perform a loop over all macro-timesteps - regardless of whether
C     the simulation is steady or time-dependent

      DO ITNS=1,NITNS
      
      ITNSR=0

      IF (ISTAT.NE.0) THEN
      
C       If the simulation is time-dependent, copy the current solution
C       into DTM0 - a backup vector of the content of DU. In case
C       of an error, DU can be restored that way to repeat the
C       calculation.
      
        CALL ZTIME(TTT0)
        CALL LCP1(DWORK(KU1),DWORK(L(LTM0)),NUP)
        
C       In case of linear extrapolation with adaptive time-stepping 
C       technique 3 (prediction with
C       repetition and time-error handling), back up the current
C       start vector in DTMLO. So we can restore it if we have to
C       repeat the current time-step.
        
        IF ((ITEXL.NE.0).AND.(ABS(IADTIM).EQ.3)) THEN
        
C         Make a backup of the start vector of the current time-step
C         (currently in DTML) into the backup vector LTMLO. That way
C         we can restore the start vector if time-step control forces
C         us to do so.
        
          CALL LCP1(DWORK(L(LTML)),DWORK(L(LTMLO)),NUP)
          
C         Make a backup of the time-step weights that is used
C         for linear extrapolation in NSDEF.
C         The TIMLxx-variables are COMMON-block variables.
          
          TIMO31=TIML31
          TIMO32=TIML32
          TIMO11=TIML11
          TIMO12=TIML12
          
        ENDIF
        
        CALL ZTIME(TTT1)
        TTLC=TTLC+TTT1-TTT0
        
      ENDIF ! (ISTAT.NE.0)

C     Here the real calculation part begins. The following mark
C     is used to jump back and repeat the calculation in case
C     of an error.

110   CONTINUE

C     TSTEP (parameter from the DAT-file) defines the current
C     time step. 
C     Make a backup of TSTEP in TSTEPO - the overall time-step size - 
C     and multiply TSTEP with 3.

      TSTEPO=    TSTEP
      TSTEP =3D0*TSTEPO
      TIMNSH=TIMENS

      IF (MSHOW.GE.2) WRITE(MTERM,*)
      IF (MSHOW.GE.1) WRITE(MFILE,*)

C     IADTIM is = 0 if either a stationary simulation should be
C     performed (set this way in the initialisation routine),
C     or if the user switched it off in the DAT file.

      IF (IADTIM.NE.0) THEN

C=======================================================================
C *** first time step: 1 step (predictor step)
C=======================================================================

C       At this point we know, that any kind of adaptive time-stepping
C       should be performed. At first we will calculate one large step
C       with time-step size 3*TSTEP. This is the predictor-step
C       in the adaptive time-stepping algorithm.
C
C       This calculation branch is only called if a really 
C       time-dependent simulation should be performed. In case of a
C       stationary simulation, IADTIM=0 and so this code is not
C       executed. In that case the solution will be calculated later.
C
C       NITNSL counts the current number of substeps in the timestep.
C       Since we calculate only one large step, we set NITNSL to 1.

        NITNSL=1

        CALL ZTIME(TTT0)
        
C       The variable ITEXL is = 0 normally, except for when linear
C       extrapolation in time (INLMIN=INLMAX=1) is used.
C       In this case this variable indicates where we are
C       inside of the algorithm.
C
C       ITEXL describes a status machine that corresponds to the
C       content of LTML as well as TIML11,TIML12,TIML31,TIML32,
C       and defines the behaviour of the extrapolation handling.
C       The values are defined as following:
C        = -1 : Predictor-step with no substeps, i.e. time-step
C               with step-size 3xTSTEP. Use TIM1x for adaptive time-
C               stepping.
C               No information about previous time-step available,
C               i.e. LTML does not define u^(n-1)
C        =  1 : Anywhere during the simulation, predictor-step
C               with step-size 3xTSTEP. Use TIM1x for adaptive time-
C               stepping.
C               LTML identifies the solution u^(n-1)
C        = -3 : Standard macro time step consisting of 3 substeps, 
C               i.e. time-step with step-size TSTEP. 
C               Use TIM3x for adaptive time-stepping.
C               No information about previous time-step available,
C               i.e. LTML does not define u^(n-1).
C        =  3 : Anywhere during the simulation, small/standard time-step
C               with step-size TSTEP. 
C               Use TIM3x for adaptive time-stepping.
C               LTML identifies the solution u^(n-1).
C
C       When linear extrapolation should be used, ITEXL must be = 1
C       on entry on this routine and will be = 1 when this routine ends.
C       The variable is set before the call of MGSTP/NSDEF 
C       appropriately to define whether information about u^(n-1) is 
C       available or not and whether we are computing a macro step
C       or not.
        
        IF (ITEXL.NE.0) THEN
        
          IF (ITNS.EQ.1) THEN
      
C           We are in the very first step of the simulation.
C           This is marked by ITEXL=-1. From the second step on
C           we reset ITEXL back to 1.
        
            ITEXL =-1
            
C           As long as ITEXL < 0, NSDEF will not perform any
C           extrapolation.
C
C           In the very first step at all, copy the starting vector
C           to the array DTML. In later steps this is updated to 
C           always hold the iterate u^(n-1) of the previous time
C           step.
          
            CALL LCP1(DWORK(KU1),DWORK(L(LTML)),NUP)
          
          ELSE
        
C           Reset ITEXL to 1 to mark that we are no more in the
C           very first step of the time-dependent simulation.
C
C           This will switch NSDEF to use TIML1x to use for
C           linear extrapolation.
        
            ITEXL = 1
            
C           Copy the extrapolation weight of the last small-scale
C           substep to TIML11 to use it as coefficient in the
C           extrapolation.
            
            TIML11=TIML32
            
          ENDIF ! (ITNS.EQ.1)
          
        ENDIF ! (ITEXL.NE.0)

C       Hack the use of the Fractional Step Theta scheme:
C       If Fractional Step is active (IFRSTP=1 from the DAT-file),
C       switch it off to use a standard scheme to perform
C       a simple (i.e. "quick-and-dirty") calculation. Remember,
C       this is only a predictor step - it's not so important how
C       accurate it is, i.e. if we use an accurate scheme or not!
C
C       Save the old status and restore it later after the solver
C       is finished.

        IFRSTH=0
        
        IF (IFRSTP.EQ.1) THEN
        
C         Deactivate fractional step, activate Crank-Nicolson:
C         Set Theta-weight to THETA = 1/2. Save the old THETA-value
C         (coming from the DAT-file) and restore it later.
        
          THETAH=THETA
          THETA =0.5D0
          IFRSTH=1
          IFRSTP=0
          
        END IF

        IF (MSHOW.GE.2) WRITE(MTERM,2)
        IF (MSHOW.GE.2) WRITE(MTERM,2002) ITNS,TIMENS,TSTEP
        IF (MSHOW.GE.2) WRITE(MTERM,2)
        IF (MSHOW.GE.1) WRITE(MFILE,2)
        IF (MSHOW.GE.1) WRITE(MFILE,2002) ITNS,TIMENS,TSTEP
        IF (MSHOW.GE.1) WRITE(MFILE,2)
        CALL ZTIME(TTT1)
        TTLC=TTLC+TTT1-TTT0

C       Start the calculation routine to calculate one 
C       macro-time-step with stepsize 3xTSTEP.

        BSTOP=.FALSE.
        BNL1 =.FALSE.
        CALL MGSTP(MFILE,MSHOW,NITNSL,LRHS,BSTOP,BNL1,NMG)
        NMGU=NMGU+NMG

        NSUBST=NSUBST+NITNSL
        
C       Result of the calculation:
C         DU   -> solution vector, updated from u^n to u^(n+1)
C         DRHS -> updated RHS-vector
C       If linear extrapolation is active:
C         DTML -> solution of previous time step;
C                 updated from u^(n-1) to u^n
C
C       Reactivate the FS Theta-Scheme in case we switched it off before.

        IF (IFRSTH.EQ.1) THEN
          THETA =THETAH
          IFRSTP=1
        ENDIF

        CALL ZTIME(TTT0)
        ISETLV=2
        ILEV=NLEV
        CALL  SETLEV (ISETLV)
        CALL ZTIME(TTT1)
        TTLC=TTLC+TTT1-TTT0

C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================

C       The linear solver broke down for any reason (BSTOP=true). 
C       Perhaps the time-step was wrong?
C       The parameter IREPIT from the DAT-file (where NNITR holds a
C       local copy from) grants up to IREPIT repetitions of the
C       calculation. So try to repeat the calculation.

        IF ((BSTOP).AND.(ITNSR.LE.NNITR)) THEN
        
          CALL ZTIME(TTT0)
          ITNSR=ITNSR+1
          
C         Restore the start vector for a new run from DTM0 back to DU
          
          CALL LCP1(DWORK(L(LTM0)),DWORK(KU1),NUP)
          
          IF ((ITEXL.NE.0).AND.(ABS(IADTIM).EQ.3)) THEN

C           Restore u^(n-1) in DTML for linear extrapolation,
C           i.e. copy DTMLO back to DTML. 
C           Then restart the current time-step.

            CALL LCP1(DWORK(L(LTMLO)),DWORK(L(LTML)),NUP)
            
C           Restore the time-step weights used in adaptive
C           time stepping when linear extrapolation is active.
            
            TIML31=TIMO31
            TIML32=TIMO32
            TIML11=TIMO11
            TIML12=TIMO12
            
          ENDIF
          
C         Reduce the time-step size to 1/SQRT(DTFACT), depending 
C         on the parameter DTFACT in the DAT-file.
          
          TSTEPN=TSTEPO/SQRT(DTFACT)
          TSTEP=TSTEPN
          TIMENS=TIMNSH
          
          CALL ZTIME(TTT1)
          TTLC=TTLC+TTT1-TTT0
          
C         Jump back and repeat the calculation.
          
          GOTO 110
          
        ENDIF ! ((BSTOP).AND.(ITNSR.LE.NNITR))
        
        CALL ZTIME(TTT0)
        
C       Copy the current solution vector from DU to DTM3. 
C       DTM3 always holds the solution of the end of the time-stepping
C       procedure. Since we performed one large step, thus "jumping
C       over" the substeps, we can directly copy the solution vector
C       into the final solution vector.
C
C       If in the copy command there is a REAL() around the 2nd DWORK, the
C       behavior of the residuals in the solver is more like in the case of
C       real matrices. This is the only reason why this loop is not solved
C       by LCP1...
      
        DO IUP=1,NUP
          DWORK(L(LTM3)+IUP-1)=DWORK(KU1+IUP-1)
        END DO
 
        CALL ZTIME(TTT1)
        TTLC=TTLC+TTT1-TTT0

C       Reset the simulation time to the beginning of the time-step

        TIMENS=TIMNSH

      END IF ! (IADTIM.NE.0)

C=======================================================================
C *** second time step: 3 steps instationary or 1 step stationary
C=======================================================================

C     Now we return to the calculation of the "real" solution with the
C     standard time-step size. In case of a time-dependent simulation
C     we now perform 3 steps with standard time-step size TSTEP.
C     in case of a stationary simulation we only calculate once.

      IF (ISTAT.NE.0) THEN
      
        CALL ZTIME(TTT0)
        
C       Copy back the starting vector of the time step from DTM0
C       to DU, so we can begin again with the smaller time-step.
        
        CALL LCP1(DWORK(L(LTM0)),DWORK(KU1),NUP)
        
C       Now we perform 3 substeps to calculate the macrostep, 
C       so set NITNSL to 3:
        
        NITNSL=3

C       Normally we don't use linear extrapolation, so ITEXL=0. 

        IF (ITEXL.NE.0) THEN

          IF (ITNS.EQ.1) THEN

C           If we are using linear extrapolation, in the
C           very first step of the simulation set ITEXL to -3, which
C           deactivates any linear extrapolation (since there is
C           no solution u^(n-1) which could be used for that).
C
C           Furthermore in the very first time-step of the simulation,
C           save the starting vector into DTML, the backup of the 
C           start-vector in the current time-step. This will be
C           our u^0 in the linear extrapolation in the 2nd time
C           step, when we also have u^1 and want to extrapolate
C           to get an approximation u~ to u^2.

            ITEXL=-3
            CALL LCP1(DWORK(KU1),DWORK(L(LTML)),NUP)
            
          ELSE
          
C           In all later steps, set ITEXL to 3.
C           This will switch NSDEF to use TIML1x to use for
C           linear extrapolation.

            ITEXL= 3
          
          ENDIF
          
        ENDIF

C       Reduce the time-step size to 1/3. Remember that TSTEP was
C       multiplied by 3 at the beginning of the subroutine in case
C       of a possible large predictor-step - so this command reduces
C       TSTEP to the value in the DAT-file.

        TSTEP=TSTEP/DBLE(NITNSL)
        
        IF (MSHOW.GE.2) WRITE(MTERM,2)
        IF (MSHOW.GE.2) WRITE(MTERM,2001) ITNS,TIMENS,TSTEP
        IF (MSHOW.GE.2) WRITE(MTERM,2)
        IF (MSHOW.GE.1) WRITE(MFILE,2)
        IF (MSHOW.GE.1) WRITE(MFILE,2001) ITNS,TIMENS,TSTEP
        IF (MSHOW.GE.1) WRITE(MFILE,2)
        CALL ZTIME(TTT1)
        TTLC=TTLC+TTT1-TTT0
        
      ELSE
      
C       No, this is a stationary simulation - so only calculate
C       one solution by setting NITNSL to 1.
C       As TSTEP is not used we set it to 1 for sure.
      
        NITNSL=1
        TSTEP=1D0
        
      ENDIF ! (ISTAT.NE.0)

C     BUG !!! ??? Is the extrapolation information (TIMLxx) 
C                 not restored for the new calculation ???

C     Now calculate the macro time-step calculation routine.
C     Depending on NITNSL (=1 or 3) this will calculate the
C     solution by 1 substep (in case of a stationary simulation)
C     or by 3 substeps (in case of a time-dependent simulation).

      BSTOP=.FALSE.
      BNL2 =.FALSE.
      CALL MGSTP(MFILE,MSHOW,NITNSL,LRHS,BSTOP,BNL2,NMG)
      NMGU=NMGU+NMG
      IF (IER.NE.0) GOTO 99998

C       Result of the calculation:
C         DU   -> solution vector, updated from u^n to u^(n+1)
C         RHS  -> updated RHS-vector
C       If linear extrapolation is active:
C         DTML -> solution of previous time step;
C                 updated from u^(n-1) to u^n
      
C     Reset the time-step size to the value that it had in the
C     corrector step. This has only an effect in a time-dependent
C     simulation - in a stationary simulation or in a simulation
C     without adaptive time-stepping, the result is reset to the
C     value of the DAT-file directly afterwards.
      
      TSTEP=TSTEP*DBLE(NITNSL)
      NSUBST=NSUBST+NITNSL

C     In case of a stationary simulation or a simulation without
C     adaptive time-stepping, that's all here:
C     Reset TSTEP and jump out of the time-step handling.

      IF ((ISTAT.EQ.0).OR.(IADTIM.EQ.0)) THEN
        TSTEP=TSTEPO
        GOTO 900
      ENDIF

      CALL ZTIME(TTT0)
      ISETLV=2
      ILEV=NLEV
      CALL  SETLEV (ISETLV)

      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0

C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================

C     The solver broke down for any reason, or the nonlinear substeps
C     don't converge, or...
C     Perhaps the time-step was wrong?
C     The parameter IREPIT from the DAT-file (where NNITR holds a
C     local copy from) grants up to IREPIT repetitions of the
C     calculation. So try to repeat the calculation.

      IF (((BSTOP).OR.(.NOT.BNL2)).AND.(ITNSR.LE.NNITR)) THEN
      
        IF (((.NOT.BSTOP).AND.(.NOT.BNL2)).AND.(ABS(IADTIM).EQ.1)) 
     *        GOTO 130
     
        CALL ZTIME(TTT0)
        ITNSR=ITNSR+1

C       Restore the solution vector from the backup DTM0
         
        CALL LCP1(DWORK(L(LTM0)),DWORK(KU1),NUP)
        
        IF ((ITEXL.NE.0).AND.(ABS(IADTIM).EQ.3)) THEN

C         Restore u^(n-1) in DTML for linear extrapolation,
C         i.e. copy DTMLO back to DTML. 
C         Then restart the current time-step.

          CALL LCP1(DWORK(L(LTMLO)),DWORK(L(LTML)),NUP)
          
C         Restore the time-step weights used in adaptive
C         time stepping when linear extrapolation is active.
C         Otherwise the content of TIMLxx is ignored.
          
          TIML31=TIMO31
          TIML32=TIMO32
          TIML11=TIMO11
          TIML12=TIMO12
        ENDIF
        
C       If the nonlinear iteration broke down, reduce the
C       time-step to 1/SQRT(DTFACT) with DTFACT from the DAT-file.
C       If the linear solver broke down (BSTOP=1),
C       reduce the time-step by 1/DTFACT.
        
        IF (.NOT.BNL2) TSTEPN=TSTEPO/SQRT(DTFACT)
        IF (BSTOP)     TSTEPN=TSTEPO/DTFACT
        
        TSTEP=TSTEPN
        TIMENS=TIMNSH
        CALL ZTIME(TTT1)
        TTLC=TTLC+TTT1-TTT0
        
C       Repeat the calculation:
        
        GOTO 110
        
      ENDIF ! (((BSTOP).OR.(.NOT.BNL2)).AND.(ITNSR.LE.NNITR))
      
C     Local jump-mark to jump out of the previous IF-command...      
      
130   CONTINUE

C=======================================================================
C     Time step control
C=======================================================================

      CALL ZTIME(TTT0)
      ISETLV=2
      ILEV=NLEV
      CALL  SETLEV (ISETLV)

C     DU contains the result of the calculation with the standard
C     time-steps, DTM3 contains the result of the predictor-step
C     above. Subtract both solutions and save the result in DD

      DO IUP=1,NUP
        DWORK(L(LD1)+IUP-1)=DWORK(KU1+IUP-1)-DWORK(L(LTM3)+IUP-1)
      END DO

C     Calculate some residual information. This will help us later
C     to calculate the time-step.

      RELU21=0D0
      RELUM1=0D0
      RELP21=0D0
      RELPM1=0D0

      IF ((IEPSAD.GE.1).OR.(IEPSAD.EQ.-1)) THEN
        CALL LL21(DWORK(KU1)   ,2*NU,DSXN)
        RELU20=DSXN/(SQRT(DBLE(2*NU)))
        IF (RELU20.LE.1D0) RELU20=1D0
        CALL LL21(DWORK(L(LD1)),2*NU,DSXN)
        RELU21=DSXN/(SQRT(DBLE(2*NU))*RELU20)
      ENDIF

      IF ((IEPSAD.GE.1).OR.(IEPSAD.EQ.-2)) THEN
       CALL LLI1(DWORK(KU1   ),2*NU,DSXN,INDMAX)
       RELUM0=DSXN
       IF (RELUM0.LE.1D0) RELUM0=1D0
       CALL LLI1(DWORK(L(LD1)),2*NU,DSXN,INDMAX)
       RELUM1=DSXN/RELUM0
      ENDIF

      IF ((IEPSAD.GE.1).OR.(IEPSAD.EQ.-3)) THEN
       CALL LL21(DWORK(KU1   +2*NU),NP,DSXN)
       RELP20=DSXN/(SQRT(DBLE(NP)))
       IF (RELP20.LE.1D0) RELP20=1D0
       CALL LL21(DWORK(L(LD1)+2*NU),NP,DSXN)
       RELP21=DSXN/(SQRT(DBLE(NP))*RELP20)
      ENDIF

      IF ((IEPSAD.GE.1).OR.(IEPSAD.EQ.-4)) THEN
       CALL LLI1(DWORK(KU1   +2*NU),NP,DSXN,INDMAX)
       RELPM0=DSXN
       IF (RELPM0.LE.1D0) RELPM0=1D0
       CALL LLI1(DWORK(L(LD1)+2*NU),NP,DSXN,INDMAX)
       RELPM1=DSXN/RELPM0
      ENDIF

C     Depending on IEPSAD, set up RELT1=J(.), our error functional:

      IF (ABS(IEPSAD).EQ.1)  RELT1=RELU21
      IF (ABS(IEPSAD).EQ.2)  RELT1=RELUM1
      IF (ABS(IEPSAD).EQ.3)  RELT1=RELP21
      IF (ABS(IEPSAD).EQ.4)  RELT1=RELPM1
      IF (     IEPSAD.EQ.5)  RELT1=MAX(RELU21,RELP21)
      IF (     IEPSAD.EQ.6)  RELT1=MAX(RELUM1,RELPM1)
      IF (     IEPSAD.EQ.7)  RELT1=MAX(RELU21,RELUM1,RELP21,RELPM1)
      IF (     IEPSAD.EQ.8)  RELT1=MIN(RELU21,RELUM1,RELP21,RELPM1)

      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0

C=======================================================================

C     Now let's calculate the time-step - if we are allowed to...

      IF (ABS(IADTIM).GE.1) THEN

C       Adaptive time-stepping is active.
C       Make the necessary adjustments for the time-step
      
        CALL ZTIME(TTT0)
        
C       Using the values in TIMExx, EPSADx and IADIN, adjust (i.e.
C       calculate) a new stopping criterion EPSAD.
        
        CALL CRITAD(TIMEIN,TIMENS,TIMEST,EPSADI,EPSADL,EPSAD,IADIN)

        IF ((IFRSTP.NE.1).AND.(THETA.EQ.1D0)) THEN

C         Time approximation of 1st order; used Euler one-step
C         scheme. The error can be represented as:
C           J(v_k) - J(v) = k e(v) + O(k^2)

          TSTEPN=TSTEPO*2D0*EPSAD/RELT1
          CTPROJ=2D0

        ELSE

C         Time approximation of 1st order; used for Fractional
C         step Theta scheme and Crank Nicolson. The error can be 
C         represented as:
C           J(v_k) - J(v) = k^2 e(v) + O(k^4)

          TSTEPN=TSTEPO*SQRT(8D0*EPSAD/RELT1)
          CTPROJ=8D0

        ENDIF

C       When the nonlinear solver broke down and we are in a time
C       stepping mode that allowes repetition of the time step,
C       bound the time step: TSTEPN must result in the interval
C       [ TSTEP/SQRT(DTFACT) .. TSTEP ]:

        IF ((ABS(IADTIM).GT.1)) THEN
          IF ((TSTEPN.GT.TSTEPO).AND.(.NOT.BNL1)) TSTEPN=TSTEPO
          IF ((TSTEPN.LT.TSTEPO/SQRT(DTFACT)).AND.(.NOT.BNL1)) 
     *        TSTEPN=TSTEPO/SQRT(DTFACT)
        ENDIF

C       TSTEPN must not be smaller than TSTEPO/DTFACT

        IF  (TSTEPN.LT.TSTEPO/DTFACT) TSTEPN=TSTEPO/DTFACT

C       Use another upper bound if the solution was calculated
C       by repetition of a time step:

        IF (ITNSR.GT.0) THEN
          HTFACT=DTFACT**(1D0/DBLE(ITNSR+1))
          IF (TSTEPN.GT.TSTEPO*HTFACT) TSTEPN=TSTEPO*HTFACT
        ELSE       
          IF (TSTEPN.GT.TSTEPO*DTFACT) TSTEPN=TSTEPO*DTFACT       
        ENDIF       

C       Bound the time step to the interval DTMIN..DTMAX

        IF (TSTEPN.LT.DTMIN) TSTEPN=DTMIN       
        IF (TSTEPN.GT.DTMAX) TSTEPN=DTMAX       

        ITYP=0
        DRELTN=TSTEPN/TSTEPO
        TSTEP =TSTEPN
        
C       In adaptive time stepping control technique 3 (prediction with
C       repetition and time-error control), if the calculated
C       residual DRELTN is too large, set ITYP=1 to indicate that
C       the time-step has to be repeated.
        
        IF ((ABS(IADTIM).EQ.3).AND.(DRELTN.LT.EPSADU)) ITYP=1

        IF (MSHOW.GE.2) WRITE(MTERM,3)
        IF (MSHOW.GE.1) WRITE(MFILE,3)
        IF (MSHOW.GE.2) WRITE(MTERM,10002) TSTEPO,
     *     RELU21/CTPROJ,RELUM1/CTPROJ,RELP21/CTPROJ,RELPM1/CTPROJ
        IF (MSHOW.GE.1) WRITE(MFILE,10002) TSTEPO,
     *     RELU21/CTPROJ,RELUM1/CTPROJ,RELP21/CTPROJ,RELPM1/CTPROJ

        IF (MSHOW.GE.2) WRITE(MTERM,3001) ITYP,ITNS,TSTEPN,TSTEPO
        IF (MSHOW.GE.1) WRITE(MFILE,3001) ITYP,ITNS,TSTEPN,TSTEPO
        IF (MSHOW.GE.2) WRITE(MTERM,3)
        IF (MSHOW.GE.2) WRITE(MTERM,*)
        IF (MSHOW.GE.1) WRITE(MFILE,3)
        IF (MSHOW.GE.1) WRITE(MFILE,*)

        CALL ZTIME(TTT1)
        TTLC=TTLC+TTT1-TTT0

C       The parameter IREPIT from the DAT-file (where NNITR holds a
C       local copy from) grants up to IREPIT repetitions of the
C       calculation. If ITYP=1, we have to repeat the current 
C       time-step - as long as there are repetitions left.

        IF ((ITYP.EQ.1).AND.(ITNSR.LE.NNITR)) THEN
          
          CALL ZTIME(TTT0)
          ITNSR=ITNSR+1
          
C         Restore the solution vector from the backup DTM0
          
          CALL LCP1(DWORK(L(LTM0)),DWORK(KU1),NUP)
          
          IF ((ITEXL.NE.0).AND.(ABS(IADTIM).EQ.3)) THEN

C           Restore u^(n-1) in DTML for linear extrapolation,
C           i.e. copy DTMLO back to DTML. 
C           Then restart the current time-step.

            CALL LCP1(DWORK(L(LTMLO)),DWORK(L(LTML)),NUP)
            
C           Restore the time-step weights used in adaptive
C           time stepping when linear extrapolation is active.
C           Otherwise the content of TIMLxx is ignored.
            
            TIML31=TIMO31
            TIML32=TIMO32
            TIML11=TIMO11
            TIML12=TIMO12
          END IF
          
C         Restore the simulation time to the beginning of the macrostep
          
          TIMENS=TIMNSH
          CALL ZTIME(TTT1)
          
          TTLC=TTLC+TTT1-TTT0
          
C         Jump back to repeat the time-step
          
          GOTO 110
          
        END IF ! ((ITYP.EQ.1).AND.(ITNSR.LE.NNITR))

C=======================================================================
C     Extrapolation step
C=======================================================================
      
        CALL ZTIME(TTT0)
       
C       Set up extrapolation weights for the Theta scheme.
C       For 2nd order time accuracy (Fractional step, Crank Nicolson)
C       we use different weights than for 1st order Euler...
       
        IF ((IFRSTP.NE.1).AND.(THETA.EQ.1D0)) THEN
          EX1= 4D0/3D0
          EX2=-1D0/3D0
        ELSE
          EX1= 9D0/8D0
          EX2=-1D0/8D0
        END IF

        IF (((IADTIM.LT.-1).AND.(BNL1)).OR.(IADTIM.EQ.-1)) THEN
       
C         Calculate the new solution vector by weighted combination
C         of the standard-time-step solution DU and the predictor-step
C         in DTM3. 
       
          DO IUP=1,NUP
            DWORK(KU1+IUP-1)= EX1*DWORK(KU1+IUP-1)
     *                       +EX2*DWORK(L(LTM3)+IUP-1)
          END DO
        
        ENDIF

        CALL ZTIME(TTT1)
        TTLC=TTLC+TTT1-TTT0

      END IF ! ABS(IADTIM).GE.1)

C=======================================================================

C     The following jump mark is used for stationary simulations or
C     simulations without adaptive time-stepping. These types
C     of simulations skip the calculation of the time-steps.

900   CONTINUE

      CALL ZTIME(TTT0)
      ISETLV=2
      ILEV=NLEV
      CALL  SETLEV (ISETLV)

C=======================================================================
C     Error analysis for postprocessing
C=======================================================================

      IF (IERANA.GE.1)  THEN
       ICUBER=IERANA
       IF (IELT.EQ.0) 
     *  CALL ELPQU(DWORK(KU1),DWORK(KU2),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             DWORK(L(LCORVG)),E031,ICUBER,UE,UEX,UEY,MFILE)
       IF (IELT.EQ.1) 
     *  CALL ELPQU(DWORK(KU1),DWORK(KU2),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             DWORK(L(LCORVG)),E030,ICUBER,UE,UEX,UEY,MFILE)
       IF (IELT.EQ.2) 
     *  CALL ELPQN(DWORK(KU1),DWORK(KU2),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             DWORK(L(LCORVG)),EM31,ICUBER,UE,UEX,UEY,MFILE)
       IF (IELT.EQ.3) 
     *  CALL ELPQN(DWORK(KU1),DWORK(KU2),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             DWORK(L(LCORVG)),EM30,ICUBER,UE,UEX,UEY,MFILE)

       CALL ELPQP(DWORK(KP),DWORK(L(LD1)+2*NU),KWORK(L(LVERT)),
     *            KWORK(L(LMID)),DWORK(L(LCORVG)),E010,ICUBER,PE,MFILE)

       IF (IERANA.NE.1) THEN
        ICUBER=1
        IF (IELT.EQ.0) 
     *   CALL ELPQU(DWORK(KU1),DWORK(KU2),KWORK(L(LVERT)),
     *              KWORK(L(LMID)),DWORK(L(LCORVG)),E031,ICUBER,
     *              UE,UEX,UEY,MFILE)
        IF (IELT.EQ.1) 
     *   CALL ELPQU(DWORK(KU1),DWORK(KU2),KWORK(L(LVERT)),
     *              KWORK(L(LMID)),DWORK(L(LCORVG)),E030,ICUBER,
     *              UE,UEX,UEY,MFILE)
        IF (IELT.EQ.2) 
     *   CALL ELPQN(DWORK(KU1),DWORK(KU2),KWORK(L(LVERT)),
     *              KWORK(L(LMID)),DWORK(L(LCORVG)),EM31,ICUBER,
     *              UE,UEX,UEY,MFILE)
        IF (IELT.EQ.3) 
     *   CALL ELPQN(DWORK(KU1),DWORK(KU2),KWORK(L(LVERT)),
     *              KWORK(L(LMID)),DWORK(L(LCORVG)),EM30,ICUBER,
     *              UE,UEX,UEY,MFILE)

        CALL ELPQP(DWORK(KP),DWORK(L(LD1)+2*NU),KWORK(L(LVERT)),
     *             KWORK(L(LMID)),DWORK(L(LCORVG)),E010,ICUBER,PE,MFILE)
       ENDIF

      ENDIF

C     Call postprocessing routine; calculate drag, lift, write GMV,...

      IFPOST=1
      CALL FPOST(IFPOST,IFILEN,ITFILM,UE,MSHOW,DFWVI,DAWVI,0,'')

      CALL ZTIME(TTT1)
      TTPOST=TTPOST+TTT1-TTT0

C=======================================================================

C     If the current simulation is of stationary type, stop the
C     calculation here. Jump out of the time-stepping loop.

      IF (ISTAT.EQ.0) GOTO 999

C=======================================================================

C     Here we definitely have a time-dependent calculation.
C     Remember, in DTM0 there is still the start solution from where we
C     calculated the next time step. Substract that from the new
C     solution in DU, save the result in DD and calculate some norms
C     according to the IEPSAD-parameter in the DAT-file.

      CALL ZTIME(TTT0)
      CALL LCP1(DWORK(L(LTM0)),DWORK(L(LD1)),NUP)
      CALL LLC1(DWORK(KU1)    ,DWORK(L(LD1)),NUP,-1D0,1D0)

C     Calculate some residual information for output

      RELU2=0D0
      RELP2=0D0
      RELUM=0D0
      RELPM=0D0

      CALL LL21(DWORK(L(LD1)),2*NU,DSXN)
      RELU2=DSXN/(SQRT(DBLE(2*NU))*3D0*TSTEPO)
      CALL LL21(DWORK(L(LD1)+2*NU),NP,DSXN)
      RELP2=DSXN/(SQRT(DBLE(NP))*3D0*TSTEPO)

      IF ((ABS(IEPSAD).EQ.2).OR.(IEPSAD.GE.5)) THEN
        CALL LLI1(DWORK(L(LD1)),2*NU,DSXN,INDMAX)
        RELUM=DSXN/(3D0*TSTEPO)
      ENDIF

      IF ((ABS(IEPSAD).EQ.4).OR.(IEPSAD.GE.5)) THEN
        CALL LLI1(DWORK(L(LD1)+2*NU),NP,DSXN,INDMAX)
        RELPM=DSXN/(3D0*TSTEPO)
      ENDIF

      IF (ABS(IEPSAD).EQ.1)  RELT=RELU2
      IF (ABS(IEPSAD).EQ.2)  RELT=RELUM
      IF (ABS(IEPSAD).EQ.3)  RELT=RELP2
      IF (ABS(IEPSAD).EQ.4)  RELT=RELPM
      IF (    IEPSAD .EQ.5)  RELT=MAX(RELU2,RELP2)
      IF (    IEPSAD .EQ.6)  RELT=MAX(RELUM,RELPM)
      IF (    IEPSAD .EQ.7)  RELT=MAX(RELU2,RELUM,RELP2,RELPM)
      IF (    IEPSAD .EQ.8)  RELT=MIN(RELU2,RELUM,RELP2,RELPM)

      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0

      IF (MSHOW.GE.2) WRITE(MTERM,2)
      IF (MSHOW.GE.2) WRITE(MTERM,20001) ITNS,NSUBST,TIMENS,RELU2,RELP2,
     *                                   RELT
      IF (MSHOW.GE.2) WRITE(MTERM,2)
      IF (MSHOW.GE.2) WRITE(MTERM,*)

      IF (MSHOW.GE.1) WRITE(MFILE,2)
      IF (MSHOW.GE.0) WRITE(MFILE,20001) ITNS,NSUBST,TIMENS,RELU2,RELP2,
     *                                   RELT
      IF (MSHOW.GE.1) WRITE(MFILE,2)
      IF (MSHOW.GE.1) WRITE(MFILE,*)

C     Stop the calculation if the time-derivative is too low
C     (-> simulation got stationary), or if the maximum simulation
C     time is reached.

      IF ((RELT.LE.EPSNS).OR.(TIMENS.GE.TIMEMX)) GOTO 999

      CALL ZTIME(TTT1)
      TTTOTH=TTT1-TTTSUM
      IF (MSHOW.GE.3) WRITE(MTERM,*) 'total time : ', TTTOTH
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'total time : ', TTTOTH
      IF (MSHOW.GE.3) WRITE(MTERM,*) 'mavec time : ', TTADF
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'mavec time : ', TTADF
      IF (MSHOW.GE.3) WRITE(MTERM,*) 'konv. time : ', TTUPW
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'konv. time : ', TTUPW
      IF (MSHOW.GE.3) WRITE(MTERM,*) 'bdry  time : ', TTBDR
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'bdry  time : ', TTBDR
      IF (MSHOW.GE.3) WRITE(MTERM,*) 'LC    time : ', TTLC
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'LC    time : ', TTLC
      IF (MSHOW.GE.3) WRITE(MTERM,*) 'MG    time : ', TTMG
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'MG    time : ', TTMG

C     That's it with the time-dependent, adaptive time-stepping
C     simulation loop. Proceed with the next macro-timestep.

      END DO ! ITNS=1,NITNS

C     At this point the time-dependent simulation is completed.
C     Only the time-dependent simulation will reach this point,
C     as the stationary simulation jumps out of the caluculation
C     above!
C     Reset ITNS to NITNS as at the end of the simulation,
C     the DO-loop brings ITNS to NITNS+1...

      ITNS=NITNS
      
C Release the memory used:

      IF (ISTAT.NE.0) THEN
        
        ISETLV=2
        ILEV=NLMAX
        CALL SETLEV (ISETLV)
        
        IF (IRHS.GE.1) THEN
          CALL ZDISP(0,LRHS  ,'DRHS  ')
          IF (IER.NE.0) GOTO 99998
        ENDIF

        CALL ZDISP(0,LTM0  ,'DMT0  ')
        IF (IER.NE.0) GOTO 99998

        IF (ITEXL.NE.0) THEN
          CALL ZDISP(0,LTML ,'DMTL  ')
          IF (IER.NE.0) GOTO 99998
          IF (ABS(IADTIM).EQ.3) THEN
            CALL ZDISP(0,LTMLO,'DMTLO ')
            IF (IER.NE.0) GOTO 99998
          ENDIF
        ENDIF

        IF (IADTIM.NE.0) THEN
          CALL ZDISP(0,LTM3  ,'DMT3  ')
          IF (IER.NE.0) GOTO 99998
        ENDIF
      ENDIF
      
999   GOTO 99999

C=======================================================================
C     Error case
C=======================================================================
99998 WRITE(MTERM,*) 'IER', IER
      WRITE(MTERM,*) 'IN SUBROUTINE ',SUB

99999 RETURN

   1  FORMAT(80('-'))
   2  FORMAT(80('$'))
   3  FORMAT(80('@'))
1000  FORMAT (6E12.5)
1001  FORMAT(' IT DIV-L2',3X,'RHOMGP')
1003  FORMAT(I3,2(D9.2))
2001  FORMAT('MACRO STEP ',I4,' AT TIME = ',D10.3,
     *       ' WITH 3 STEPS: DT3 = ',D10.3)
2002  FORMAT('MACRO STEP ',I4,' AT TIME = ',D10.3,
     *       ' WITH 1 STEP : DT1 = ',D10.3)
3001  FORMAT('CHOICE ',I2,'(',I4,')',1X,'   ---- NEW DT = ',D10.3,
     *       ' -- OLD DT = ',D10.3)
10002 FORMAT ('OLD DT=',D9.2,
     *        1X,'U(L2)=',D9.2,1X,'U(MX)=',D9.2,
     *        1X,'P(L2)=',D9.2,1X,'P(MX)=',D9.2)
10003 FORMAT ('OLD DT=',D9.2,
     *        1X,'FWREL=',D9.2,1X,'AWREL=',D9.2,
     *        1X,'FW3  =',D9.2,1X,'FA3  =',D9.2)
20001 FORMAT ('#',I4,1X,'(',I4,')',1X,'TIME=',D10.3,1X,'RELU(L2)=',
     *        D9.2,1X,'RELP(L2)=',D9.2,1X,'REL=',D9.2)

      END
