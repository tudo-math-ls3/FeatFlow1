************************************************************************
      SUBROUTINE  XOPTCN (KU1,KU2,KP,KU1OLD,KU2OLD,KPOLD,KF1,KF2,KFP,
     *                    KD1,KD2,KDP,KAUX1,KAUX2,KAUXP,KA1,KCOLA,KLDA,
     *                    KB1,KB2,KCOLB,KLDB,KST1,KM1,NA,NU,NP,
     *                    DELU,DELP,OMEGA,KMBD,NMBD,INEUM)
************************************************************************
* Computes the new OMEGA-parameter for the defect correction.
* Updates the KUx-arrays with the defect correction approach:
*
*   u^(l+1) = u^l + OMEGA * Y
*
* with Y=(KU1,KU2,KP) the solution from the Oseen equation and
* u^l=(KU1OLD,KU2OLD,KPOLD) the old solution.
* Calculates the norms of the changes into DELxx.
*
*    Input:
*      - vectors U,UOLD,F
*      - OMEGA for the new calculation of the nonlinear block A
*      - matrices and pointers A,...,ST,NA
*      - number of equations NU
*    Output:
*      - updated vector U
*      - optimal relaxation parameter OMEGA
*      - maximum relative changes  DELU 
*      - note that D is changed to  D=K*V with correction V=U-UOLD 
*
* In:
*   KU1OLD,
*   KU2OLD,
*   KPOLD   - pointers to array [1..*] of double
*             vector u^l of the previous time step
*   KU1,
*   KU2,
*   KP      - pointers to array [1..*] of double
*             update vector Y for the new time step, calculated
*             by solving the Oseen equation
*   KF1,
*   KF2,
*   KFP     - pointers to array [1..*] of double
*             right hand side of the Oseen equation
*   KD1,
*   KD2,
*   KDP     - pointers to array [1..*] of double
*             Defect vector = RHS of the Oseen equation
*
*   KAUX1  - pointer to array [1..NU] of double
*   KAUX2  - pointer to array [1..NU] of double
*   KAUXP  - pointer to array [1..NP] of double
*            temporary vector used for some computations
*   KA1,
*   KCOLA,
*   KLDA,
*   KB1,KB2,
*   KCOLB, 
*   KLDB   - current system matrix + structure   
*
*   KM     - mass matrix
*   NA,
*   NU,
*   NP     - size of matrix and vectors
*
* In (from COMMON-blocks):
*   OMGMIN - minimal OMEGA, from DAT-file
*   OMGMAX - maximal OMEGA, from DAT-file
*
* Out:
*   KU1,
*   KU2,
*   KP     - new solution vectors, updated with the defect correction
*            approach
*   OMEGA  - new relaxation parameter OMEGA for defect correction
*   DELU,
*   DELP,
*   DELT   - Norm of new solution vectors / relative change of
*            solution vectors - see below at the end of the file.
*
*   KAUX1,  
*   KAUX2  - is overwritten with some temporary data (cf. COMMON block)
*   KD1,
*   KD2,
*   KDP    - These vectors are overwritten with some temporary data
*
*   KA1    - is rebuild at the point u^l-omegaold*Y
************************************************************************
      IMPLICIT NONE
      
C common blocks

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgpar.inc'
      
      INCLUDE 'cns.inc'
      INCLUDE 'cnsparfrac.inc'

C *** user COMMON blocks

      INCLUDE 'cinidat.inc'

C parameters

      INTEGER KU1,KU2,KP,KU1OLD,KU2OLD,KPOLD,KF1,KF2,KFP,
     *        KD1,KD2,KDP,KAUX1,KAUX2,KAUXP,KA1,KCOLA,KLDA,
     *        KB1,KB2,KCOLB,KLDB,KST1,KM1,NA,NU,NP,
     *        KMBD(*),NMBD,INEUM
      DOUBLE PRECISION DELU,DELP,OMEGA

C local variables

      DOUBLE PRECISION TTT0, TTT1, AA1, AA2
      DOUBLE PRECISION SKV1, SKV2, SSS, DELU1, DELU2
      DOUBLE PRECISION DELT1, DELT2, DELT
      INTEGER INDU1, INDU2, INDP, INDT1, INDT2, INDT


C-----------------------------------------------------------------------
C     E X T E R N A L S
C-----------------------------------------------------------------------
C *** Coefficient of stiffness matrix
      DOUBLE PRECISION COEFFN
      EXTERNAL COEFFN
C *** definition of finite elements
      EXTERNAL E030,E031,EM31,EM30
      EXTERNAL MTMUL

C-----------------------------------------------------------------------

      IF (OMGMIN.EQ.OMGMAX) THEN
      
        IF (OMGMIN.LT.0D0) THEN

C         If OMGMIN=OMGMAX<0, we have nothing to do - no 
C         relative changes are calculated. Cancel here.

          OMEGA=ABS(OMGMIN)
          DELU=0D0
          DELP =0D0
          GOTO 99999
        ELSE

C         For OMGMIN=OMGMAX>0, there's no choice for OMEGA - 
C         prescribe it directly and skip its calculation.

          OMEGA=OMGMIN
          GOTO 999
          
        ENDIF
        
      ENDIF

C***********************************************************************

C     The first big part now is to calculate a new OMEGA parameter
C     with OMGMIN < OMEGA < OMGMAX.
C
C     The defect correction for a problem like T(u)u=f has the form
C
C           u^(l+1)  =  u^l  -  OMEGA * C * ( T(u^l)u^l - f )
C
C     with an appropriate preconditioner C (see below).
C     In our case, this iteration system can be written as:
C
C     (KU1)    (KU1)               ( [ KST1        KB1] (KU1)   (KF1) )
C     (KU2) := (KU2) - OMEGA * C * ( [       KST1  KB2] (KU2) - (KF2) )
C     (KP)     (KP)                ( [ KB1^T KB2^T  0 ] (KP)    (KFP) )
C
C                                  |----------------------------------|
C                                           = (KD1,KD2,KDP)^T
C                              |--------------------------------------|
C                                          = -Y = -(KU1,KU2,KP)
C
C     with KST1=KST1(KU1,KU2,KP) and Y being the solution from
C     the Oseen equation with 
C
C                      [ KST1        KB1]
C        C = T(u^l) =  [       KST1  KB2]
C                      [ KB1^T KB2^T  0 ]
C
C     The parameter OMEGA is calculated as the result of the 1D
C     minimization problem:
C
C       OMEGA = min_omega || T(u^l-omega*Y)*(u^l-omega*Y) - f ||_E
C
C               < T(u^l-omegaold*Y)Y , f - T(u^l-omegaold*Y)u^l >
C            ~= -------------------------------------------------
C                  < T(u^l-omegaold*Y)Y , T(u^l-omegaold*Y)Y >
C
C     when choosing omegaold=previous omega, which is a good choice
C     as one can see by linearization (see p. 187 (170), Turek's book).
C
C     in the Euclidian norm ||.||_E to the Euclidian scalar 
C     product <.,.>.

C=======================================================================
C *** Calculate on AUX1/2 the linearization point: AUX = UOLD+OMEGA*Y
C=======================================================================

C     We calculate UOLD+OMEGA*Y instead of UOLD-OMEGA*Y, because the
C     "-"-sign is already incorporated in Y -- Y is calculated
C     as residuum "b-Ax" rather than as a defect "Ax-b"...

      CALL ZTIME(TTT0)
      AA1=1.0D0
      AA2=OMEGA
      CALL LCP1 (DWORK(KU1   ),DWORK(KAUX1),NU)
      CALL LCP1 (DWORK(KU2   ),DWORK(KAUX2),NU)
      CALL LCP1 (DWORK(KP    ),DWORK(KAUXP),NP)

      CALL LLC1 (DWORK(KU1OLD),DWORK(KAUX1),NU,AA1,AA2)
      CALL LLC1 (DWORK(KU2OLD),DWORK(KAUX2),NU,AA1,AA2)
      CALL LLC1 (DWORK(KPOLD) ,DWORK(KAUXP),NP,AA1,AA2)
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0

C=======================================================================
C *** First term of scalar product in the nominator
C
C     Calculate the new nonlinear block A at the
C     point AUX = u^l-omegaold*Y
C=======================================================================

C     Construct the linear part of the nonlinear matrix on the lower
C     level. As there is no defect vector to be calculated (we
C     obtained it already by restriction), we use a routine here
C     that only builds the KA1-matrix.

      CALL ZTIME(TTT0)
      CALL XMADF3(KM1,KST1,KA1,KCOLA,KLDA,NA,NU,THSTEP,ISTAT)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0

C     Up to now our system matrix has consisted only of linear
C     terms:
C
C         KA1 = [  M*I  +  THSTEP * (-nu * Laplace(.))  ]
C
C     but this is not enough... the nonlinear part is still missing.
C     So we have to add the following term to KA1:
C
C         THSTEP * u grad(.)
C
C     what will finally result in the system matrix
C
C         KA1 = [  M*I  +  THSTEP * (-nu * Laplace(.))  ] + THSTEP * u grad(.)
C             = [  M*I  +  THSTEP * N(u) ]
C
C     where u=KAUX in this case. Don't perform any correction of defect
C     vectors (IDEF=0), as we only want to have the matrix modified.

      CALL ZTIME(TTT0)
      IF (IUPW.EQ.1) THEN
       CALL GUPWD (DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX1),DWORK(KAUX2),
     *             1D0,0D0,DWORK(KAUX1),DWORK(KAUX2),
     *             DWORK(KAUX1),DWORK(KAUX2),DWORK(KA1),KWORK(KCOLA),
     *             KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             DWORK(L(LCORVG)),0)
      ELSE
       IF (IELT.EQ.0) 
     *  CALL SUPWDG(DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX1),DWORK(KAUX2),
     *             1D0,0D0,DWORK(KAUX1),DWORK(KAUX2),
     *             DWORK(KAUX1),DWORK(KAUX2),DWORK(KA1),NA,KWORK(KCOLA),
     *             KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             DWORK(L(LCORVG)),E031,COEFFN,0,1D0)
       IF (IELT.EQ.1) 
     *  CALL SUPWDG(DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX1),DWORK(KAUX2),
     *             1D0,0D0,DWORK(KAUX1),DWORK(KAUX2),
     *             DWORK(KAUX1),DWORK(KAUX2),DWORK(KA1),NA,KWORK(KCOLA),
     *             KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             DWORK(L(LCORVG)),E030,COEFFN,0,1D0)
       IF (IELT.EQ.2) 
     *  CALL SUPWNP(DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX1),DWORK(KAUX2),
     *             1D0,0D0,DWORK(KAUX1),DWORK(KAUX2),
     *             DWORK(KAUX1),DWORK(KAUX2),DWORK(KA1),NA,KWORK(KCOLA),
     *             KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             DWORK(L(LCORVG)),EM31,COEFFN,0,1D0)
       IF (IELT.EQ.3) 
     *  CALL SUPWNP(DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX1),DWORK(KAUX2),
     *             1D0,0D0,DWORK(KAUX1),DWORK(KAUX2),
     *             DWORK(KAUX1),DWORK(KAUX2),DWORK(KA1),NA,KWORK(KCOLA),
     *             KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             DWORK(L(LCORVG)),EM30,COEFFN,0,1D0)
      ENDIF
      CALL ZTIME(TTT1)
      TTUPW=TTUPW+TTT1-TTT0

C     Incorporate Dirichlet-boundary into the matrix. Replace all
C     rows corresponding to Dirichlet nodes by unit vectors.

      CALL ZTIME(TTT0)
      CALL BDRYA (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),KMBD,NMBD)
      CALL ZTIME(TTT1)
      TTBDR=TTBDR+TTT1-TTT0

C=======================================================================
C *** Second term of the scalar product in the nominator
C     Calculate the defect  D = F-T*UOLD,
C     store it in KD1,KD2,KDP, overwriting the old defect
C=======================================================================

      CALL ZTIME(TTT0)
      CALL LCP1 (DWORK(KF1),DWORK(KD1),NU)
      CALL LCP1 (DWORK(KF2),DWORK(KD2),NU)
      CALL LCP1 (DWORK(KFP),DWORK(KDP),NP)
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0

      CALL ZTIME(TTT0)
      CALL MTMUL(DWORK(KD1),DWORK(KD2),DWORK(KDP),
     *             DWORK(KU1OLD),DWORK(KU2OLD),DWORK(KPOLD),-1D0,1D0,
     *             DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *             DWORK(KB1),DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *             NU,NP,KMBD,NMBD,INEUM)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0

C=======================================================================
C *** All terms in the fraction:
C     Calculate the value  AUX = T*Y
C=======================================================================

      CALL ZTIME(TTT0)
      CALL MTMUL(DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUXP),
     *             DWORK(KU1),DWORK(KU2),DWORK(KP),1D0,0D0,
     *             DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *             DWORK(KB1),DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *             NU,NP,KMBD,NMBD,INEUM)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0

C=======================================================================
C *** Calculation of the fraction terms.
C     Calculate nominator:    SKV1:= (T*Y,D)   = (AUX,D)
C     Calculate denominator:  SKV2:= (T*Y,T*Y) = (AUX,AUX)
C=======================================================================

      CALL ZTIME(TTT0)

      CALL LSP1 (DWORK(KAUX1),DWORK(KD1)  ,NU,SKV1)
      CALL LSP1 (DWORK(KAUX2),DWORK(KD2)  ,NU,SSS)
      SKV1=SKV1+SSS
      CALL LSP1 (DWORK(KAUXP),DWORK(KDP)  ,NP,SSS)
      SKV1=SKV1+SSS
      
      CALL LSP1 (DWORK(KAUX1),DWORK(KAUX1)  ,NU,SKV2)
      CALL LSP1 (DWORK(KAUX2),DWORK(KAUX2)  ,NU,SSS)
      SKV2=SKV2+SSS
      CALL LSP1 (DWORK(KAUXP),DWORK(KAUXP)  ,NP,SSS)
      SKV2=SKV2+SSS

      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0

      IF (SKV2.LT. 1.0D-40) THEN
        WRITE(MTERM,*) 'ERROR in OPTCOR: SKVKV is nearly zero'
        STOP
      ENDIF

C     Ok, we have the nominator and the denominator. Divide them
C     by each other to calculate the new OMEGA.

      OMEGA=SKV1/SKV2
      
C     And make sure it's in the allowed range:

      IF (OMEGA.LT.ABS(OMGMIN)) OMEGA=ABS(OMGMIN)
      IF (OMEGA.GT.ABS(OMGMAX)) OMEGA=ABS(OMGMAX)

C     That's it, we have our new Omega.

C=======================================================================
C *** Update of the solution vector.
C     Calculate  Unew  :=  u^l + Omega*Y  =  UOLD+OMEGA*Y 
C=======================================================================

999   CALL ZTIME(TTT0)

C=======================================================================
C *** Calculate maximum changes   DELU
C=======================================================================

C     Calculate the maximum norm of Y=(KU1, KU2, KP) and save the
C     it individually in DELU1, DELU2, DELP.
C     (INDxxx is the index of the component giving rise to the maximum)

      CALL  LLI1 (DWORK(KU1),NU,DELU1,INDU1)
      CALL  LLI1 (DWORK(KU2),NU,DELU2,INDU2)
      CALL  LLI1 (DWORK(KP) ,NP,DELP ,INDP )

C=======================================================================
C *** Defect correction.
C     Update the solution   U := UOLD + OMEGA*Y
C=======================================================================

      AA1= 1.D0
      AA2= OMEGA
      CALL  LLC1 (DWORK(KU1OLD),DWORK(KU1),NU,AA1,AA2)
      CALL  LLC1 (DWORK(KU2OLD),DWORK(KU2),NU,AA1,AA2)
      CALL  LLC1 (DWORK(KPOLD) ,DWORK(KP) ,NP,AA1,AA2)

C     This gives us our new solution (KU1, KU2, KP), overwriting our
C     old correction vector Y - this one is now history...

C=======================================================================
C *** relative maximum changes   
C=======================================================================

C     This simply calculates some postprocessing values of the relative
C     change in the solution.
C
C     Maximum norm of solution vector:
C
C                  || YP ||_max       || Pnew - Pold ||_max
C       DELP := ------------------- = ---------------------
C                  || P ||_max           || Pnew ||_max
C
C
C     Relative change of solution vector:
C
C                || (Y1,Y2) ||_max    || Unew - Uold ||_max
C       DELU := ------------------- = --------------------- 
C               || (KU1,KU2) ||_max       || Unew ||_max

      CALL  LLI1 (DWORK(KU1),NU,DELT1,INDT1)
      CALL  LLI1 (DWORK(KU2),NU,DELT2,INDT2)
      DELT=MAX(DELT1,DELT2)
      IF (ABS(DELT).LT.1D-8) DELT=1D0
      DELU=MAX(DELU1,DELU2)/DELT

      CALL  LLI1 (DWORK(KP) ,NP,DELT,INDT)
      IF (ABS(DELT).LT.1D-8) DELT=1D0
      DELP=DELP/DELT

      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0

99999 END
