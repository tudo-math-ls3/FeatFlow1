************************************************************************
      
      SUBROUTINE MGSTP(MFILE,MSHOW,NITNSL,LRHS,BSTOP,BNLEND,NMG)
      
************************************************************************
*
*   Purpose: - performs 1 macro (time) step of coupled method
*
* This routine calculates the solution u^(n+1) for the next time step
* based on the solution u^n in the current time step. 
*
* For stationary simulations, based on u^0, u^1 is calculated which
* is the results of the complete simulation.
* For true time-dependent simulation, u^(n+1) is calculated by u^n
* in NITNSL substeps with the current Theta-Scheme (which is configured
* by the NSPAR/NSFRAC COMMON-block).
*
* In:
*  MFILE  - file handle to write output to
*  MSHOW  - level of output
*  NITNSL - Number of substeps to compute the solution.
*           =1: use 1-step algorithm (e.g. stationary simulation or
*               predictor-step)
*           =3: use 3-step algorithm (Theta-scheme)
*  LRHS   - Handle to right-hand-side vector.
*           Can be 0 for homogenuous problems.
*           Must be <> 0 for nonsteady inhomogenuous problems.
*
* From COMMON-blocks:
*  DU     - u^n
*  DTML   - if linear extrapolation is used: receives u^(n-1);
*           otherwise ignored.
*  ITEXL  - <> 0: use linear extrapolation; in this case ITEXL
*           describes the current state of the algortithm
*           (predictor step or not, first step at all or not,...)
*  
* Out:
*  BSTOP  - =true: linear solver broke down
*  BNLEND - =false: nonlinear iteration broke down
*  NMG    - Number of multigrid iterations needed during calculation
*  DRHS   - the RHS-vector identified by LRHS is updated according
*           to the information in the next time-step 
*
* Updated variables in COMMON-blocks:
*  DU     - receives u^(n+1)
*  DTML   - if linear extrapolation is used: receives u^n
*
* Return:
*  NMG - number of MG steps that are used during the iteration.
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT NONE
      
C include all necessary COMMON blocks
      
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
      
C constants

      INTEGER NBLOCF
      PARAMETER (NBLOCF=2)
      LOGICAL BCONF(NBLOCF),BSNGLF(NBLOCF)

      DIMENSION ARRDF(NBLOCF)
      CHARACTER ARRDF*6

      INTEGER KFN(NBLOCF)
      
      DATA BCONF  /.FALSE.,.FALSE./,ARRDF/'DF1   ','DF2   '/
      DATA BSNGLF /.FALSE.,.FALSE./,KFN/1,1/
      
      SAVE BCONF, BSNGLF, ARRDF

C parameters

      INTEGER MFILE,MSHOW,NITNSL,LRHS,NMG
      LOGICAL BSTOP,BNLEND
      
C externals

C *** Parametrization of the domain
      DOUBLE PRECISION PARX,PARY,TMAX
      EXTERNAL PARX,PARY,TMAX
      
C *** Coefficient of stiffness matrix, right hand side, exact solution
      DOUBLE PRECISION COEFFN,RHS,UE,PE,UEX,UEY,NEUDAT
      EXTERNAL COEFFN,RHS,UE,PE,UEX,UEY,NEUDAT
      
C *** definition of finite elements
      EXTERNAL E030,E031,EM31,EM30

C local variables

      INTEGER KF(NNAB,NBLOCF),LF(NBLOCF)

      INTEGER ITNSL,ITMOD,ISETLV,NMBD,ICLRF,LF1,LF2,NMG1
      DOUBLE PRECISION TTT0,TTT1,TSTEPH,TSTEPB,TMSTEP,TRSTEP,TOSTEP
      DOUBLE PRECISION DSXN,RELU2,RELP2

C=======================================================================

C     Initialize descriptors for RHS-vector.
C     This is used in the assembling of additional RHS terms from the
C     current time step.

      KF(1,1)=1
      KF(1,2)=1

      NMG = 0

C=======================================================================
C *** Begin of nonsteady loop
C=======================================================================

C     Perform at most NITNSL iterations to calculate the iterate of the
C     next time step. For stationary simulations, the loop is canceled
C     below. In time-dependent simulations this performs NITNSL
C     substeps to calculate u^(n+1) from u^n.

      DO ITNSL=1,NITNSL
      
        CALL ZTIME(TTT0)

C       Set ITMOD to 1,2 or 3, depending on ITNSL. ITMOD identifies the
C       current substep in the Fractional step Theta scheme.

        ITMOD=MOD(ITNSL-1,3)+1

        ISETLV=2
        ILEV=NLEV
        CALL  SETLEV (ISETLV)

C       Should we use Fractional-Step?
C       IFRSTP is a COMMON-block variable from the DAT-file.
C       Adjust time indicator according to the current step in the
C       FS-Theta-Scheme.
C
C       TSTEPH denotes the actual length of the substep.
C       For 1-step schemes this is just TSTEP, and 3 steps with
C       length TSTEP=TSTEPH gives the macrostep with length 3xTSTEP.
C       For the FS Theta scheme the macrostep has still length
C       3xTSTEP, but the substeps are not equally distributed
C       by substep-length TSTEP. Instead the 1st and 3rd substep
C       have length 
C            Theta * steplength of macrostep   
C       while the second step has length
C            Theta' * steplength of macrostep
C       All three steps again sum up to the macroste length 3xTSTEP.

        IF (IFRSTP.EQ.1) THEN
        
          IF (ITMOD.EQ.2) THEN
            TSTEPH=3D0*TSTEP*THETAP
          ELSE
            TSTEPH=3D0*TSTEP*THETA
          ENDIF
          
        ELSE

          TSTEPH=TSTEP
          
        ENDIF ! (IFRSTP.EQ.1)

C       Adjust the current simulation time according to the Theta-scheme - 
C       except for in stationary simulations, there we can reset TIMENS to 0
C       as it's not used.

        TIMENS=TIMENS+TSTEPH
        IF (ISTAT.EQ.0) TIMENS=0D0     

        CALL ZTIME(TTT1)
        TTLC=TTLC+TTT1-TTT0

C=======================================================================

C       Now comes the big difference between stationary and instationary
C       simulation. In instationary simulations we have to prepare
C       a lot of stuff for the next iterate:

        IF (ISTAT.NE.0) THEN

          CALL ZTIME(TTT0)
          
C         Initialize Theta-Scheme weights, depending on whether we
C         use Fractional Step or a 1-step-scheme (see below).
          
          IF (IFRSTP.EQ.1) THEN
          
            IF (ITMOD.EQ.2) THEN
              TSTEPB= 3D0*TSTEP*THETAP
              THSTEP=-3D0*TSTEP*FALPHA*THETAP
              TMSTEP= 3D0*TSTEP*FALPHA*THETA
              TRSTEP= 3D0*TSTEP*FBETA *THETAP
            ELSE
              TSTEPB= 3D0*TSTEP*THETA
              THSTEP=-3D0*TSTEP*FBETA*THETA
              TMSTEP= 3D0*TSTEP*FALPHA*THETA
              TRSTEP= 3D0*TSTEP*FALPHA*THETA
            ENDIF !(ITMOD.EQ.2)
            
          ELSE
            TSTEPB= TSTEP
            THSTEP= TSTEP*(THETA-1D0)
            TMSTEP= TSTEP*THETA
            TRSTEP= TSTEP*THETA
          ENDIF ! (IFRSTP.EQ.1)    

C         If we use linear extrapolation, in all steps from
C         the second macrostep on (ITEXL > 0)  we have to
C         prepare extrapolation weights.
C         Remember:
C           ITEXL = -1: first macro-step, prediction step
C           ITEXL = -3: first macro-step, standard
C           ITEXL =  1: subsequent macro-step, prediction step,
C                       or: stationary simulation
C           ITEXL =  3: subsequent macro-step, standard

C         The linear extrapolation approximates the velocity u^(l+1) with
C         u~ computed by (cf. p. 205 (187f), Turek's book)
C
C         u~ = [ dT2 / (dT2-dT1) ] * u(t-dT1)  +  [ -dT1 / (dT2-dT1) ] * u(t-dT2)
C         
C         The TIMLxy-parameter are used as:
C
C           TIMLx1 = dT1       
C           TIMLx2 = (dT2-dT1) 
C
C         Set TIMLx2 = TIMLx1 to remember the previous time 
C         step =>   dT1 -> dT2-dT1
C         Set TIMLx1 to the new stepsize TSTEPH.
C         (So TIMLxy is used as a queue in y. The new step is put into
C         TIMLx1 while the old step goes to TIMLx2.)

          IF (ITEXL.EQ.1) THEN
            TIML12=TIML11
            TIML11=TSTEPH
          ENDIF

          IF (ITEXL.EQ.3) THEN
            TIML32=TIML31
            TIML31=TSTEPH
          ENDIF       

          CALL ZTIME(TTT1)
          TTLC=TTLC+TTT1-TTT0

C***********************************************************************
C         The Nav.St.-Equation (cf. p. 38 (43), Turek's book)
C
C             u_t + N(u)u + grad(p) = f,   div(u)=0
C
C         is discretised in time
C
C             [ I + Theta k N(u)]u + k grad(p) = g,   div(u)=0
C
C         and space
C
C             [ M + Theta k N(u_h)]u_h + k B p_h = g_h,   B^T u_h = 0
C
C         resulting in the abstract system
C
C             S(u_h)u_h  +  k B p_h  =  g_h,   B^T u_h = 0
C
C         with
C
C             S(u) = [alpha*M + Theta_1*nu*k*L + Theta_2*k*K(u)] u
C
C         Here, alpha=0 for stationationary and alpha=1 for
C         instationary problems,L=KST1=Laplace-Matrix.
C         The time-step constant k=TSTEPB in front of the pressure 
C         B p_h in the abstract system above can be neglected by 
C         scaling the pressure in every time-step.
C***********************************************************************

C         At first we scale the pressure by the step-length.
C         That's all for the handling of the pressure.

          CALL LLC1(DWORK(KP),DWORK(KP),NP,0D0,TSTEPH)

C***********************************************************************

          CALL ZTIME(TTT0)
          
          IF (IBDR.GE.2) THEN
          
C           We have time-dependent Dirichlet/Neumann boundary 
C           conditions! (IBDR=2 from the DAT-file). Implement them!
          
            DO ILEV=NLMIN,NLMAX
              NVT =KNVT(ILEV)
              NMT =KNMT(ILEV)
              NEL =KNEL(ILEV)
              NVBD=KNVBD(ILEV)

C             Restore nodal-property vector KNPR back to original.
C             The backup was stored in KLNPRO in the initialization 
C             phase.    

              CALL XLCP3(KLNPRO(ILEV),KLNPR(ILEV),NVT+NMT)
              IF (IER.NE.0) GOTO 99999

              NMBD =0
              INEUM=0
              KNMBD(ILEV)=0
              CALL ZDISP(0,KLMBD(ILEV),'KMBD  ')
              CALL ZDISP(0,KLDBD(ILEV),'DDBD  ')
              IF (IER.NE.0) GOTO 99999
              CALL ZNEW (NMT,3,KLMBD(ILEV),'KMBD  ')
              CALL ZNEW (NMT,1,KLDBD(ILEV),'DDBD  ')
              IF (IER.NE.0) GOTO 99999

C             Implement boundary information; calculate NMBD.
C             This modifies KNPR and KMBD according to the information
C             what is Dirichlet-boundary and what not.
C             Furthermore this routine returns in INEUM whether there
C             are Neumann components at all in our problem.
C             Remark: INEUM is a COMMON-block variable, so available
C             also to other routines that are called later.

              CALL BDRNEU(KWORK(L(KLMBD(ILEV))),KWORK(L(KLVBD(ILEV))),
     *              KWORK(L(KLEBD(ILEV))),KWORK(L(KLVERT(ILEV))),
     *              KWORK(L(KLMID(ILEV))),KWORK(L(KLNPR(ILEV))),
     *              DWORK(L(KLDBD(ILEV))),DWORK(L(KLMBDP(ILEV))),
     *              DWORK(L(KLCVG(ILEV))),KWORK(L(KLADJ(ILEV))),
     *              KNVBD(ILEV),KNVT(ILEV),KNEL(ILEV),
     *              NMBD,INEUM)
              KNMBD(ILEV)=NMBD

            END DO ! ILEV
       
          ENDIF ! (IBDR.GE.2)
      
          CALL ZTIME(TTT1)
          TTBDR=TTBDR+TTT1-TTT0
      
C***********************************************************************

          ISETLV=2
          ILEV=NLEV
          CALL SETLEV(ISETLV)

          CALL ZTIME(TTT0)
         
C         We are still in the handling of the instationary case here!
C         Now we set up the right hand side of the Theta-scheme.
C         This is based on the right-hand-side vector(s) of the equation
C         (identified by the handle LRHS) and the nonlinear term
C         depending on the current solution (see page 165/166 (151/152) 
C         in Turek's CFD-book):  
C
C         [I - theta_1 k N(u^n)]u^n + Theta_2 k f^(n+1) + Theta_3 k f^n
C
C         with k=Delta(t) being the step size of the substep. 
C         The RHS of the Theta-scheme is written into the arrays
C         identified by KFx.  
C
C         Classical 1-step scheme
C         -----------------------
C         In the classical 1-step scheme, the RHS has the form 
C         (cf. p.152):  
C
C         [I - (1-Theta) k N(u^n)]u^n + Theta k f^(n+1) + (1-Theta) k f^n
C            ^^^^^^^^^^^^^              ^^^^^^^           ^^^^^^^^^^^
C              THSTEP                   TRSTEP             -THSTEP
C
C         In case that we have steaty (in-)homogenuous RHS, we have
C         f^n = f^(n+1) = RHS and there fore the latter both terms
C         reduce:
C
C           Theta k f^(n+1) + (1-Theta) k f^n
C         = Theta k f^(n+1) + (1-Theta) k f^(n+1)
C         = k f^(n+1)
C         =    k      RHS
C           ^^^^^^^
C           TSTEPB
C
C         Fractional step Theta-scheme
C         ----------------------------
C         This is a little bit more complicated due to the different
C         substeps. As RHS we (again) sum together:
C
C         [I + THSTEP N(u^n)]u^n + TRSTEP f^(n+1) + -THSTEP f^n
C
C         But this time the step-variables have different values
C         depending on the step
C
C             Step:       1                2                  3
C         Variable:
C            THSTEP   -beta*Theta*K  alpha*Theta'*K      -beta*Theta*K
C            TRSTEP   alpha*Theta*K  -beta*Theta'*K      alpha*Theta*K
C
C         or in formulas, respectively:
C
C             Step:       1                2                  3
C         Variable:
C            THSTEP   -beta*Theta*K  alpha*(1-2Theta)*K   -beta*Theta*K
C            TRSTEP   alpha*Theta*K  -beta*(1-2Theta)*K   alpha*Theta*K
C
C         with:   alpha = FALPHA = (1-2Theta)/(1-Theta)
C                 beta  = FBETA  = Theta/(1-Theta)
C                 Theta = 1-2/sqrt(2)
C                 Theta'= (1-2Theta)
C
C         WARNING: The documentation on page 165 (151) in Turek's book 
C          on the Theta Scheme is different as implemented here!!!
C          This version of the Theta Scheme has never been published,
C          but was derived independently of Glowinski by Rannacher/
C          Turek! It makes more sense than the Glowinski-scheme
C          as the coefficients before f^n+1 and f^n add to 1 and
C          the RHS f^(n+1-Theta) is not used in both, the second
C          and the third step!
C
C         If we have homogenuous RHS, simply force it to 0.
         
          IF (IRHS.EQ.0) THEN
            CALL LCL1(DWORK(KF1),NUP)
          ENDIF

C         If we have steady inhomogenuous RHS, copy the calculated
C         RHS-vector from DRHS to DF1, weighted by the coefficients
C         of the Theta-scheme.    
C         The RHS is normally build in the INIT1 with the help
C         of the INDATxD.F-file.

          IF (IRHS.EQ.1) THEN

C           The RHS reduces as stated above. We only have to multiply by
C           
C                      TSTEPB  =  TRSTEP + (-THSTEP)  =  k

            CALL LLC1(DWORK(L(LRHS)),DWORK(KF1),NUP,TSTEPB,0D0)
            
          ENDIF
           
          CALL ZTIME(TTT1)
          TTLC=TTLC+TTT1-TTT0

C         If we have nonsteady inhomogenuous RHS, we have to calculate 
C         that. 

          IF (IRHS.EQ.2) THEN
          
C           We go "from back to front" when adding the terms
C           to the RHS in the Theta scheme (cf. p. 164 (150), Turek's 
C           book). Take the current RHS f^n, weight it according 
C           to the Theta-scheme with (1-Theta)*k and write it into DF1.

            CALL ZTIME(TTT0)
            CALL LLC1(DWORK(L(LRHS)),DWORK(KF1),NUP,-THSTEP,0D0)
            CALL ZTIME(TTT1)
            TTLC=TTLC+TTT1-TTT0

            CALL ZTIME(TTT0)
            
C           Implement pressure-drop boundary conditions into
C           RHS-vectors DF=(DF1,DF2,DP)
            
            IF (IBDR.GE.1) THEN
              TIMENS=TIMENS-TSTEPH
              CALL PDSET  (KWORK(L(KLVBD(NLEV))),KWORK(L(KLEBD(NLEV))),
     *                KWORK(L(KLVERT(NLEV))),KWORK(L(KLMID(NLEV))),
     *                KWORK(L(KLNPR(NLEV))),DWORK(L(LCORVG)),
     *                DWORK(L(KLMBDP(NLEV))),DWORK(KF1),DWORK(KF2),
     *                -THSTEP,KNVBD(NLEV))
              TIMENS=TIMENS+TSTEPH
            ENDIF ! (IBDR.GE.1)
            
            CALL ZTIME(TTT1)
            TTBDR=TTBDR+TTT1-TTT0

C           Assemble the RHS-vector f^(n+1) of the next time-step, 
C           obtain two new handles for it.    
            
            CALL ZTIME(TTT0)
            LF(1)=0
            LF(2)=0
            IF (IELT.EQ.0) 
     *        CALL  XVB0 (LF,NU,NBLOCF,ICLRF,E031,RHS,BCONF,KF,KFN,
     *                    ICUBF,ARRDF,BSNGLF)
            IF (IELT.EQ.1) 
     *        CALL  XVB0 (LF,NU,NBLOCF,ICLRF,E030,RHS,BCONF,KF,KFN,
     *                    ICUBF,ARRDF,BSNGLF)
            IF (IELT.EQ.2) 
     *        CALL  XVBM0(LF,NU,NBLOCF,ICLRF,EM31,RHS,BCONF,KF,KFN,
     *                    ICUBF,ARRDF,BSNGLF)
            IF (IELT.EQ.3) 
     *        CALL  XVBM0(LF,NU,NBLOCF,ICLRF,EM30,RHS,BCONF,KF,KFN,
     *                    ICUBF,ARRDF,BSNGLF)
            IF (IER.NE.0) GOTO 99999
            LF1=LF(1)
            LF2=LF(2)
            CALL ZTIME(TTT1)
            TTLC=TTLC+TTT1-TTT0

C           Implement pressure drop conditions into these newly assembled
C           RHS vectors.

            CALL ZTIME(TTT0)
            IF (IBDR.GE.1) THEN
              CALL PDSET (KWORK(L(KLVBD(NLEV))),KWORK(L(KLEBD(NLEV))),
     *               KWORK(L(KLVERT(NLEV))),KWORK(L(KLMID(NLEV))),
     *               KWORK(L(KLNPR(NLEV))),DWORK(L(LCORVG)),
     *               DWORK(L(KLMBDP(NLEV))),DWORK(L(LF1)),DWORK(L(LF2)),
     *               TRSTEP,KNVBD(NLEV))
            ENDIF ! (IBDR.GE.1)
    
            CALL ZTIME(TTT1)
            TTBDR=TTBDR+TTT1-TTT0

C           Add the newly assembled RHS-vectors to the RHS-vectors from
C           above, weighted according to the Theta-scheme with Theta*k.

            CALL ZTIME(TTT0)
            CALL  LLC1 (DWORK(L(LF1)),DWORK(KF1),NU,TRSTEP,1D0)
            CALL  LLC1 (DWORK(L(LF2)),DWORK(KF2),NU,TRSTEP,1D0)
            CALL ZTIME(TTT1)
            TTLC=TTLC+TTT1-TTT0

C           Implement Dirichlet boundary conditions into RHS-vector
            
            CALL ZTIME(TTT0)
            CALL BDRSET (DWORK(KF1),DWORK(KF2),KWORK(L(LNPR)),
     *               NVT,PARX,PARY,UE,
     *               DWORK(L(LCORVG)),KWORK(L(LADJ)),NEL,
     *               KWORK(L(LMID)),KWORK(L(LVERT)))
            CALL BDRSET (DWORK(KF1),DWORK(KF2),KWORK(L(LNPR)),
     *               NVT,PARX,PARY,UE,
     *               DWORK(L(LCORVG)),KWORK(L(LADJ)),NEL,
     *               KWORK(L(LMID)),KWORK(L(LVERT)))
            CALL ZTIME(TTT1)
            TTBDR=TTBDR+TTT1-TTT0

C           Copy the new RHS-vector to DRHS, overwriting the old one.

            CALL ZTIME(TTT0)
            CALL  LCP1 (DWORK(L(LF1)), DWORK(L(LRHS)),   NU)
            CALL  LCP1 (DWORK(L(LF2)), DWORK(L(LRHS)+NU),NU)

C           Now we can release the memory of the temporary RHS again.

            CALL  ZDISP (0,LF1,ARRDF(1))
            CALL  ZDISP (0,LF2,ARRDF(2))
            IF (IER.NE.0) GOTO 99999
            CALL ZTIME(TTT1)
            TTLC=TTLC+TTT1-TTT0
          ENDIF ! (IRHS.EQ.2)

C***********************************************************************

          CALL ZTIME(TTT0)
          ISETLV=2
          ILEV=NLEV
          CALL  SETLEV (ISETLV)

C         In case we have inhomogenuous boundary conditions (steady or
C         non-steady), (re-)implement pressure-drop conditions.

          IF (IBDR.GE.1) THEN
          
            IF ((IFRSTP.EQ.1).AND.(ITMOD.EQ.2)) TIMENS=TIMENS-TSTEPH
            
            CALL PDSET (KWORK(L(KLVBD(NLEV))),KWORK(L(KLEBD(NLEV))),
     *              KWORK(L(KLVERT(NLEV))),KWORK(L(KLMID(NLEV))),
     *              KWORK(L(KLNPR(NLEV))),DWORK(L(LCORVG)),
     *              DWORK(L(KLMBDP(NLEV))),DWORK(KF1),DWORK(KF2),
     *              TSTEPB,KNVBD(NLEV))
            
            IF ((IFRSTP.EQ.1).AND.(ITMOD.EQ.2)) TIMENS=TIMENS+TSTEPH
            
          ENDIF ! (IBDR.GE.1)
          CALL ZTIME(TTT1)
          TTBDR=TTBDR+TTT1-TTT0

C======================================================================

C         Set up linear part of RHS-vector of Theta-scheme as well as
C         linear part of nonlinear system matrix. Use the current
C         vector in KFx as source and overwrite it.
C         -->   KA1 = [  M*I  +  THSTEP * (-nu * Laplace(.))  ]
C         -->   KFx = [ M*u^n + THSTEP*(-nu*Laplace(u^n)) ]  +  KFx
C                                                              ^^^^^
C                    Remember:  KFx = TRSTEP*f^(n+1) + -THSTEP*RHS
C                                   = TSTEPB*RHS  (in steady case)
C                               THSTEP = -(1-Theta)*k

          CALL ZTIME(TTT0)
          CALL XMADF1(KM1,KST1,KA1,KCOLA,KLDA,KF1,KF2,KU1,KU2,
     *                NA,NU,THSTEP)
          CALL ZTIME(TTT1)
          TTADF=TTADF+TTT1-TTT0
        
C         Now treat the convective part.
C
C         Modify the current RHS-vector KF that way, that
C             THSTEP * u grad(u)
C         is added:
C
C             KFx = KFx + THSTEP * u grad(u)
C                       ^^^^^^^^^^^^^^^^^^^^
C                 = [ M*I + THSTEP*N(u) ]u + TRSTEP*f^(n+1) - THSTEP*f^n
C
C         But because upwinding/streamline diffusion always subtracts the
C         term THSTEP * u grad(u) (in order to build a defect vector), 
C         we use a trick and switch the sign of THSTEP - then the term is
C         added!
C
C         The matrix is not to be modified (we don't need it), so we
C         call the assembling routine with IDEF=2.
C
C         There's a special trick if streamline diffusion is used, as
C         these routines also handle the case that IPRECA=4!
C         If IPRECA=4 and IMASS=1, we also include the term with mass
C         matrix into the RHS. For a proper handling, we therefore set
C         the factor in front of the mass matrix DCMASS=-1, so
C         (-M*u) is subtracted from the RHS, so M*u is added to the RHS.

          IF (THSTEP.NE.0D0) THEN

            CALL ZTIME(TTT0)
            
            THSTEP=-THSTEP
      
            IF (IUPW.EQ.1) THEN
              CALL GUPWD (DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KU1),DWORK(KU2),
     *               DWORK(KF1),DWORK(KF2),DWORK(KA1),KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),2)
            ELSE
              IF (IELT.EQ.0) 
     *          CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(KF1),DWORK(KF2),DWORK(KA1),NA,KWORK(KCOLA),
     *                KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *                DWORK(L(LCORVG)),E031,COEFFN,2,-1D0)
              IF (IELT.EQ.1) 
     *          CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(KF1),DWORK(KF2),DWORK(KA1),NA,KWORK(KCOLA),
     *                KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *                DWORK(L(LCORVG)),E030,COEFFN,2,-1D0)
              IF (IELT.EQ.2) 
     *          CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(KF1),DWORK(KF2),DWORK(KA1),NA,KWORK(KCOLA),
     *                KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *                DWORK(L(LCORVG)),EM31,COEFFN,2,-1D0)
              IF (IELT.EQ.3) 
     *          CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(KF1),DWORK(KF2),DWORK(KA1),NA,KWORK(KCOLA),
     *                KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *                DWORK(L(LCORVG)),EM30,COEFFN,2,-1D0)
            ENDIF ! (IUPW.EQ.1)
          
            THSTEP=-THSTEP
          
            CALL ZTIME(TTT1)
            TTUPW=TTUPW+TTT1-TTT0
          ELSE

C           THSTEP=0D0 - This is the case for backward Euler only.
C           (THSTEP = c(1-THETA) = 0  <->  THETA=1)

            IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN

C             The handling is slightly different because we have
C             to set THSTEP=1 during the upwinding in this case,
C             and we must use DCMASS=0.

              CALL ZTIME(TTT0)

C             Use a step length of 1D0 in call to the creation of the
C             convective part. Save the old value in TOSTEP and restore
C             later.
              
              TOSTEP=THSTEP
              THSTEP=1D0
              
              IF (IELT.EQ.0) 
     *          CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(KF1),DWORK(KF2),DWORK(KA1),NA,KWORK(KCOLA),
     *                KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *                DWORK(L(LCORVG)),E031,COEFFN,2,0D0)
              IF (IELT.EQ.1) 
     *          CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(KF1),DWORK(KF2),DWORK(KA1),NA,KWORK(KCOLA),
     *                KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *                DWORK(L(LCORVG)),E030,COEFFN,2,0D0)
              IF (IELT.EQ.2) 
     *          CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(KF1),DWORK(KF2),DWORK(KA1),NA,KWORK(KCOLA),
     *                KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *                DWORK(L(LCORVG)),EM31,COEFFN,2,0D0)
              IF (IELT.EQ.3) 
     *          CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(KF1),DWORK(KF2),DWORK(KA1),NA,KWORK(KCOLA),
     *                KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *                DWORK(L(LCORVG)),EM30,COEFFN,2,0D0)
            ENDIF ! ((IPRECA.EQ.4).AND.(IMASS.EQ.1))
            
            THSTEP=TOSTEP
            
            CALL ZTIME(TTT1)
            TTUPW=TTUPW+TTT1-TTT0
          ENDIF ! (THSTEP.NE.0D0)

C         Incorporate Dirichlet boundary conditions into current
C         solution vector.

          CALL ZTIME(TTT0)
          CALL BDRSET (DWORK(KU1),DWORK(KU2),KWORK(L(LNPR)),
     *              NVT,PARX,PARY,UE,
     *              DWORK(L(LCORVG)),KWORK(L(LADJ)),NEL,
     *              KWORK(L(LMID)),KWORK(L(LVERT)))
          CALL BDRSET (DWORK(KF1),DWORK(KF2),KWORK(L(LNPR)),
     *              NVT,PARX,PARY,UE,
     *              DWORK(L(LCORVG)),KWORK(L(LADJ)),NEL,
     *              KWORK(L(LMID)),KWORK(L(LVERT)))
          CALL ZTIME(TTT1)
          TTBDR=TTBDR+TTT1-TTT0

C         Switch THSTEP to TMSTEP - the "matrix step length weight".
C
C         This is the weight in the matrix on the left hand side
C         of the equation! As example to identify it, we consider the
C         equation in the first step of the Fractional-Step Theta-
C         scheme (cf. p. 167 (153), Turek's book):
C
C            [ I + alpha k Theta A_(n+theta) ] u_(n+Theta)  =  RHS
C                  ^^^^^^^^^^^^^ ^^^^^^^^^^^                   ^^^
C                     TMSTEP     Created in NSDEF  what we created
C                                                  above
C         As NSDEF also uses THSTEP as coefficient before the matrix,
C         we switch THSTEP to TMSTEP, what is the correct coefficient
C         for the system matrix.

          THSTEP = TMSTEP

C         That's all for the instationary case. Now we can turn to
C         solve the nonlinear problem for both, stationary
C         and instationary...

        ENDIF ! (ISTAT.NE.0)

************************************************************************
C *** fixed point defect correction for stationary NS equations
C***********************************************************************

C       Call NSDEF and perform fixed point defect correction
C       to calculate the next iterate of this substep.
C       This will update the solution vector in DU1/DU2/DP.

        CALL  NSDEF (MFILE,MSHOW,BSTOP,BNLEND,NMG1)
        
C       sum up the number of multigrid steps

        NMG=NMG+NMG1
        IF (IER.NE.0) GOTO 99999

        IF ((.NOT.BNLEND).AND.(ABS(IADTIM).GT.1)) RETURN

C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================
        IF (BSTOP) RETURN

        CALL ZTIME(TTT0)
        ILEV=NLEV
        ISETLV=2
        CALL  SETLEV (ISETLV)

C***********************************************************************

C       Scale back the pressure by the time-step.
C       Remember: p was scaled down by the time-step above,
C       so we have now to scale by the inverse of the time-step
C       to scale it back.

        CALL LLC1(DWORK(KP),DWORK(KP),NP,0D0,1D0/TSTEPH)

C***********************************************************************

        CALL LL21(DWORK(KU1),2*NU,DSXN)
        RELU2=DSXN/SQRT(DBLE(NU))
        CALL LL21(DWORK(KU1+2*NU),NP,DSXN)
        RELP2=DSXN/SQRT(DBLE(NP))
        THSTEP=TSTEPH

        CALL ZTIME(TTT1)
        TTLC=TTLC+TTT1-TTT0

        IF (MSHOW.GE.2) WRITE(MTERM,3)
        IF (MSHOW.GE.2) WRITE(MTERM,10002) ITNSL,ITNS,TIMENS,RELU2,RELP2
        IF (MSHOW.GE.2) WRITE(MTERM,3)
        IF (MSHOW.GE.2) WRITE(MTERM,*)

        IF (MSHOW.GE.1) WRITE(MFILE,3)
        IF (MSHOW.GE.1) WRITE(MFILE,10002) ITNSL,ITNS,TIMENS,RELU2,RELP2
        IF (MSHOW.GE.1) WRITE(MFILE,3)
        IF (MSHOW.GE.1) WRITE(MFILE,*)

C======================================================================
C *** Return if stationary calculation !!!
C=======================================================================

        IF (ISTAT.EQ.0) RETURN

      END DO

   1  FORMAT(80('-'))
   3  FORMAT(80('+'))
1000  FORMAT (6E12.5)
1001  FORMAT(' IT DIV-L2',3X,'RHOMGP')
1003  FORMAT(I4,2(D9.2))
10002 FORMAT ('#',I4,'(',I4,')',1X,'TIME=',D10.3,2X,'NORM(U)=',
     *        D14.7,2X,'NORM(P)=',D14.7)

99999 END
