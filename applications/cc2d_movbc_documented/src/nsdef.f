************************************************************************

      SUBROUTINE  NSDEF (MFILE,MSHOW,BSTOP,BNLEND,NMG)  

************************************************************************
*   Purpose: - solver for the stationary incompressible Navier Stokes
*              equations via fixed point defect correction plus
*              multigrid for linear Oseen problems
*            - nonlinear version:
*                   - fixed point defect correction as outer iteration
*                   - mg as solver for the linear auxiliary Oseen
*                     problems
*                   - nonlinear parameter optimization for the 
*                     correction from the linear solver step
*
* In:
*  MFILE  - file handle to write output to
*  MSHOW  - level of output
*
* Out:
*  BSTOP  - =true: linear solver broke down
*  BNLEND - =false: nonlinear iteration broke down
*  NMG - number of MG steps used 
*
* In (from COMMON blocks):
*   DU      - solution vector
*   KF1,
*   KF2,KP  - starting addtes of right-hand-side vector in DWORK
*   LD1     - handle to U/V/P-vector. Temporary vector that holds
*             the defect during the nonlinear iteration.
*   LTML    - handle to solution vector of the previous time step;
*             in the first time step: handle t a vector that
*             contains a copy of DU.
*   ITEXL  - status identifier for linear extrapolation
*   TIMxy  - coefficients for the linear extrapolation formula
*   THSTEP - double. Theta-scheme identifier.
*            =0: Forward Euler, =1: Backward Euler, =0.5: Crank Nicolson
*            For Fractional-Step theta scheme this is set to different
*            values, depending on the current step in the scheme.
*            For stationary simulations this parameter must be set
*            to 1 to include the full Laplacian matrix into the
*            system matrix:
*                    [alpha*M + THETA*KST1 + ... ] 
*   ISTAT  - whether the simulation is stationary
*   INEUM  - whether there are boundary components in our problem
*   OMGINI - Initial OMEGA-Value - from the DAT-file
*   LU1OLD - if <> 0, this must be a handle of
*               array[1..NUP] of double
*            and receives a backup of the solution vector on the finest
*            level before the iteration.
*
* Updated variables in COMMON-blocks:
*   DU    - new solution vector
*
* If linear extrapolation is active:
*   ITEXL = ABS(ITEXL) 
*         - to indicate that LTML identifies u^(n-1) when this routine is
*           finished
*   LTML  - The vector identified by LTML is updated to u^(n-1)
************************************************************************

C=======================================================================
C     Declarations
C=======================================================================

      IMPLICIT NONE
      
C main COMMON blocks
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgtria.inc'
      INCLUDE 'cmgadr.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      
      INCLUDE 'cns.inc'
      INCLUDE 'cnsparfrac.inc'
      
      INCLUDE 'cinidat.inc'

C parameters
      
      INTEGER MFILE, MSHOW, NMG
      LOGICAL BSTOP, BNLEND

C *** Arrays for multigrid modul M010 
      INTEGER KOFFX(NNLEV),KOFFD(NNLEV),KOFFB(NNLEV),KNEQ(NNLEV),
     *        KIT(NNLEV),KIT0(NNLEV)

C externals

C *** Coefficient of stiffness matrix
      EXTERNAL COEFFN
C *** definition of finite elements
      EXTERNAL E030,E031,EM30,EM31
C *** Multigrid components
      EXTERNAL  YAX,YPROL,YREST,YSM,YEX,YEXA,YDBC,YSTEP,I000

C local variables

      DOUBLE PRECISION TTT0,TTT1,DTIM1,DTIM2,A1L,A2L
      DOUBLE PRECISION RESOLD,RES0,RES,EPSRES,OMEGA
      CHARACTER CFILE*60
      
      INTEGER ISETLV,IDEFUP,INL,I1,KU1F,KU2F,KVERTF,KMIDF,KADJF
      INTEGER IRELMG,ITMG,IDEFMG,LDEF
      DOUBLE PRECISION RESU, RESDIV, RESOUT, RHOASM, DELP
      DOUBLE PRECISION RHOLMG,RHO,DELU
      LOGICAL BMG

C=======================================================================
C     Initialization
C=======================================================================

      NMG = 0

C     Initialization of the offset arrays KOFFX,KOFFB,KOFFD and KNEQ  
C     for the multigrid solver

      DO ILEV=NLMIN,NLMAX
        KOFFX(ILEV)=L(KLUP(ILEV))-1
        KOFFB(ILEV)=L(KLF12P(ILEV))-1
        KOFFD(ILEV)=L(KLAUX(ILEV))-1
        KNEQ(ILEV)=KNUP(ILEV)
      END DO
      
C     The offset of the RHS-vector on the finest level will be modified
C     to point to our defect vector, identified by the handle LD1.
C     That is because the KLF12P-handle on the finest level points
C     currently to the RHS of the nonlinear equation, but we use
C     another RHS. More precisely, we use a defect vector as RHS,
C     where the MG-solver (solving an Oseen equatiion) is a
C     preconditioner for the nonlinear iteration.

C=======================================================================
C     First generation of the nonlinear block A on level NLMAX
C     and generation of the defect vectors.
C=======================================================================

      CALL ZTIME(TTT0)
      BSTOP =.FALSE.
      BNLEND=.FALSE.

      ISETLV=2
      ILEV=NLMAX
      CALL SETLEV (ISETLV)

C     If INLMAX>1, we work with defect correction and really solve the
C     nonlinear problem with at most INLMAX nonlinear iterations.
C     If INLMAX=1, we use the linear extrapolation technique, thus
C     skipping all nonlinear iterations. In this case we only have to
C     solve a linear problem, which is later on used by the caller in
C     to extrapolate the solution for the next time step (in the
C     instationary case).

      IF (INLMAX.GT.1) THEN
        IDEFUP=1
      ELSE
        IDEFUP=0
      ENDIF

C     Copy the current RHS to the vector KD1, which is used to 
C     hold the defect during the nonlinear iteration.

      IF (INLMAX.GT.1) THEN
        CALL LCP1 (DWORK(KF1),DWORK(L(LD1)),NU)
        CALL LCP1 (DWORK(KF2),DWORK(L(LD1)+NU),NU)
        CALL LCP1 (DWORK(KFP),DWORK(L(LD1)+2*NU),NP)
      ENDIF
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0

C     Now calculate the linear part of the system matrix as well as
C     the defect vector for velocity/pressure arising from the linear
C     parts:
C
C      -->   KA1 = [  ALPHA*M*I  +  THSTEP * (-nu * Laplace(.))  ]
C      -->   (KD1) = (KD1)                                               ( KU1 )
C            (KD2) = (KD2) - [ ALPHA*M*u^n + THSTEP*(-nu*Laplace(u^n)) ] ( KU2 )
C            (KDP) = (KDP)                                               ( KP  )
C
C     with ALPHA=ISTAT. The matrix KA1 will later also serve as a
C     preconditioner in the Oseen equation. The nonlinear part
C     is subtracted later...

      CALL ZTIME(TTT0)
      CALL XMADF2(KM1,KST1,KA1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *            L(LD1),L(LD1)+NU,L(LD1)+2*NU,KU1,KU2,KP,NA,NU,NP,
     *            KWORK(L(KLMBD(NLEV))),KNMBD(NLEV),INEUM,THSTEP,ISTAT)
     
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0

      CALL ZTIME(TTT0)
      
C     Check if we have Stokes flow or Navier Stokes flow.
C
C     When calculating Navier Stokes,
C
C         du/dt + nu*Laplace(u) + u*grad(u) + grad(p) = f
C
C     we have to treat the nonlinearity u*grad(u). In contrast when
C     treating Stokes flow ("slow flow"),
C
C         du/dt + nu*Laplace(u) + grad(p) = f
C
C     the nonlinear term is neglected and so we have nothing to do
C     there to treat it.
C
C     Only exception: Stokes Flow with IPRECA=4=build always everything
C     without relying on precalculated data in the KST1-matrices.
      
      IF ((ISTOK.NE.1).OR.((ISTOK.EQ.1).AND.(IPRECA.EQ.4))) THEN

C       Treat the nonlinearity with linear extrapolation or
C       by standard handling:

        IF (ITEXL.NE.0) THEN
        
C         Linear extrapolation handling (cf. p. 204f (187f), 
C         Turek's book).
C
C         When using linear extrapolation, the velocity u~ for the 
C         nonlinearity is approximated by the weighted mean
C
C         u~ = [ dT2 / (dT2-dT1) ] * u(t-dT1)  +  [ -dT1 / (dT2-dT1) ] * u(t-dT2)
C
C         We make the following setting to approximate the
C         velocity u^(l+1): 
C
C         u(t-dT1) = u^n
C         u(t-dT2) = u^(n-1)
C         
C         DTIM1 = dT1       / TSTEP = TIMLx1
C         DTIM2 = (dT2-dT1) / TSTEP = TIMLx2
C
C         (both DTIM1 and DTIM2 are normalized by TSTEP, therefore
C         there is no step size appearent in TIMLxy!)
C
C         With this setting, the above formula results in
C
C         u~ = [ dT2 / (dT2-dT1) ] * u(t-dT1)  +  [ -dT1 / (dT2-dT1) ] * u(t-dT2)
C
C            = [ (DTIM1+DTIM2) / DTIM2 ] * u^n  +  [ DTIM1 / DTIM2 ] * u^(n-1)
C
C            = A2L * DTML  +  A1L * DU
        
          IF (ITEXL.EQ.1) THEN
            DTIM2=TIML12
            DTIM1=TIML11
          ENDIF       
          IF (ITEXL.EQ.3) THEN
            DTIM2=TIML32
            DTIM1=TIML31
          ENDIF   
              
C         In the very first step (ITEXL < 0) deactivate extrapolation - 
C         since we don't have a solution of the previous time step:

          IF (ITEXL.LT.0) THEN
            
            A1L  =0D0
            A2L  =1D0
            
C           When this routine ends, we have a solution u^(n-1) saved
C           in LTML that can be used for the next extrapolation step.
C           So remove the "-"-sign from ITEXL to indicate that we have
C           u^(n-1) in the next time step.
            
            ITEXL=ABS(ITEXL)        
            
          ELSE
            
            A1L  =-DTIM1/DTIM2
            A2L  =(DTIM1+DTIM2)/DTIM2
            
          ENDIF

          IF (IUPW.EQ.1) THEN
            CALL GUPWD(DWORK(L(LTML)),DWORK(L(LTML)+NU),
     *               DWORK(KU1),DWORK(KU2),A1L,A2L,
     *               DWORK(KU1),DWORK(KU2),
     *               DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *               DWORK(KA1),KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),IDEFUP)
          ELSE
            IF (IELT.EQ.0) 
     *        CALL SUPWDG(DWORK(L(LTML)),DWORK(L(LTML)+NU),
     *                 DWORK(KU1),DWORK(KU2),A1L,A2L,
     *                 DWORK(KU1),DWORK(KU2),
     *                 DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                 DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                 KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                 E031,COEFFN,IDEFUP,1D0)
            IF (IELT.EQ.1) 
     *        CALL SUPWDG(DWORK(L(LTML)),DWORK(L(LTML)+NU),
     *                 DWORK(KU1),DWORK(KU2),A1L,A2L,
     *                 DWORK(KU1),DWORK(KU2),
     *                 DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                 DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                 KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                 E030,COEFFN,IDEFUP,1D0)
            IF (IELT.EQ.2) 
     *        CALL SUPWNP(DWORK(L(LTML)),DWORK(L(LTML)+NU),
     *                 DWORK(KU1),DWORK(KU2),A1L,A2L,
     *                 DWORK(KU1),DWORK(KU2),
     *                 DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                 DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                 KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                 EM31,COEFFN,IDEFUP,1D0)
            IF (IELT.EQ.3) 
     *        CALL SUPWNP(DWORK(L(LTML)),DWORK(L(LTML)+NU),
     *                 DWORK(KU1),DWORK(KU2),A1L,A2L,
     *                 DWORK(KU1),DWORK(KU2),
     *                 DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                 DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                 KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                 EM30,COEFFN,IDEFUP,1D0)
          ENDIF

          IF (ITEXL.EQ.3) CALL LCP1(DWORK(KU1),DWORK(L(LTML)),NUP)

C         Make a linear combination of u^(n-1) and u^n 
C         to create DU := u^(n+1).
C         If ITEXL < 0, There will be A1L=0, A2L=1, so
C         there will be no real combination.
C         The (ITEXL.GT.0)-IF is nonsense because ITEXL is
C         set to > 0 above!

          IF (ITEXL.GT.0) CALL LLC1(DWORK(L(LTML)),DWORK(KU1),NUP,
     *                              A1L,A2L)
        ELSE

C         Standard handling routine without linear extrapolation
C
C         Add the convective part u*grad(u) to the system matrix KA1 -
C         the linear Oseen matrix that will later serve as a
C         preconditioner.
C
C         Up to now our system matrix has consisted only of linear
C         terms:
C
C             KA1 = [  M*I  +  THSTEP * (-nu * Laplace(.))  ]
C
C         but this is not enough... the nonlinear part is still missing.
C         So we have to add the following term to KA1:
C
C             THSTEP * u grad(.)
C
C         what will finally result in the system matrix
C
C             KA1 = [  M*I  +  THSTEP * (-nu * Laplace(.))  ] + THSTEP * u grad(.)
C                 = [  M*I  +  THSTEP * N(u) ]
C
C         When we perform nonlinear iteration (i.e. no linear
C         extrapolation, IDEFUP=1), add the corresponding
C         nonlinearity to the defect vectors, what is still missing.
C
C         This overall procedure then results in a matrix for the
C         "linearizes Navier Stokes" equation, i.e. the Oseen equation.
C         This matrix is later used as a preconditioner for the defect
C         vector (cf. p. 184 (168), Turek's book).
C
C         Additionally to setting up the nonlinear part, also update
C         the defect vector; subtract (Nonliner Part)*(Solution) to
C         incorporate the nonlinear defect according to
C
C        -->   (KD1) = (KD1)                 ( KU1 )
C              (KD2) = (KD2) - [ u grad(.) ] ( KU2 )
C              (KDP) = (KDP)                 ( KP  )
        
          IF (IUPW.EQ.1) THEN
            CALL GUPWD (DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                 1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                 DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                 DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *                 KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                 IDEFUP)
          ELSE
            IF (IELT.EQ.0) 
     *        CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                  1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                  DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                  E031,COEFFN,IDEFUP,1D0)
            IF (IELT.EQ.1) 
     *        CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                  1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                  DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                  E030,COEFFN,IDEFUP,1D0)
            IF (IELT.EQ.2) 
     *        CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                  1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                  DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                  EM31,COEFFN,IDEFUP,1D0)
            IF (IELT.EQ.3) 
     *        CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                  1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                  DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                  DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                  EM30,COEFFN,IDEFUP,1D0)

          ENDIF ! (IUPW.EQ.1)

        ENDIF ! (ITEXL.NE.0)

C       That's it for the nonlinearity. Now we come to the real
C       iteration.

      ENDIF ! ((ISTOK.NE.1).OR.((ISTOK.EQ.1).AND.(IPRECA.EQ.4)))
        
      CALL ZTIME(TTT1)
      TTUPW=TTUPW+TTT1-TTT0
   
C=======================================================================
C     Treat Dirichlet-nodes:
C=======================================================================

C     Replace all rows in the system matrix corresponding to Dirichlet
C     nodes by unit vectors:
   
      CALL ZTIME(TTT0)
      CALL BDRYA (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *            KWORK(L(KLMBD(ILEV))),KNMBD(ILEV))

C     And if we want to perform real nonlinear iteration,
C     set the entries in the defect vector corresponding to Dirichlet
C     nodes to 0.

      IF (INLMAX.GT.1) THEN
       CALL BDRY0 (DWORK(L(LD1)),DWORK(L(LD1)+NU),KWORK(L(KLMBD(ILEV))),
     *             KNMBD(ILEV))
      ENDIF
      CALL ZTIME(TTT1)
      TTBDR=TTBDR+TTT1-TTT0

C=======================================================================
C     Calculation of initial defects
C=======================================================================

      CALL ZTIME(TTT0)
      
C     When we want to perform nonlinear iteration, calculate the norm of
C     the initial residuum. The residuum itself was calculated
C     in the assembling-routines above (XMADF2, ...)
      
      IF (INLMAX.GT.1) THEN
      
C       Calculate the norm of the residual. Results of this routine
C       are saved into RESU,RESDIV.
      
        CALL  RESDFK(DWORK(KU1),DWORK(KU2),DWORK(KP),
     *               DWORK(L(LD1)),DWORK(L(LD1)+NU),DWORK(L(LD1)+2*NU),
     *               DWORK(KF1),DWORK(KF2),DWORK(KFP),NU,NP,RESU,RESDIV)
     
        RESOLD=SQRT(RESU*RESU+RESDIV*RESDIV)
        RES0  =RESOLD
        RES   =RESOLD
        EPSRES=DMPD*RESOLD
        
      ELSE
      
        RES   =0D0
        RESOLD=0D0
        EPSRES=0D0
        
      ENDIF

C=======================================================================
C     Solution
C=======================================================================

      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)

      IF (MSHOW.GE.2) WRITE(MTERM,1001)
      IF (MSHOW.GE.1) WRITE(MFILE,1001)
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)
      INL=0
      IF (MSHOW.GE.2) WRITE(MTERM,1002)  INL,RESU,RESDIV,RESOLD
      IF (MSHOW.GE.1) WRITE(MFILE,1002)  INL,RESU,RESDIV,RESOLD
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)

C     Initialize the OMEGA-parameter

      OMEGA=OMGINI

      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0

C=======================================================================
C *** Loop of nonlinear iteration
C=======================================================================

C     Now we start with the actual nonlinear iteration to solve the
C     problem. Out problem is to find a solution u that fulfills
C     (cf p. 178 (163) in Turek's book):
C
C         T(u)u = f
C
C     for an appropriate operator T(.) and a solution u. In our context,
C     T(.) is the discrete Navier-Stokes operator and u=(u,v,p) is the
C     velocity/pressure solution vector.
C
C     The general approach for the nonlinear here is the defect
C     correction:
C
C           u^(l+1)  =  u^l  -  OMEGA * C * ( T(u^l)u^l - f )
C
C     with l=1..INLMAX for a start vector u^1 and C being an 
C     appropriate preconditioner. (In fact, as preconditioner
C     we choose C=T(u^l)^{-1} here!) The parameter OMEGA
C     can be prescribed by the user by setting OMGMIN=OMGMAX=anything,
C     or can be calculated automatically when OMGMIN < OMGMAX with
C     the "adaptive fixed point defect correction" method.
C
C     Perform a maximum of INLMAX nonlinear iterations:

      DO INL=1,INLMAX
      
        NNONL=NNONL+1

C       Make a backup of the solution vector to DDC so that we can 
C       use it later in our parameter estimation for the correction.
C       The IF-clause is nonsense, LU1OLD must always be <> 0 -
C       otherwise there would also be an error when calling XOPTCN!

        IF (LU1OLD.NE.0) THEN
          CALL ZTIME(TTT0)
          CALL  LCP1 (DWORK(KU1),DWORK(L(LU1OLD)),NU)
          CALL  LCP1 (DWORK(KU2),DWORK(L(LU1OLD)+NU),NU)
          CALL  LCP1 (DWORK(KP ),DWORK(L(LU1OLD)+2*NU),NP)
          CALL ZTIME(TTT1)
          TTLC=TTLC+TTT1-TTT0
        ENDIF

C=======================================================================
C *** Generate the A blocks for all coarse levels
C=======================================================================

        IF (NLMAX.GT.NLMIN)  THEN
        
C         Loop through all levels, looking from the lower level
C         to the upper level always.
        
          DO ILEV=NLMAX-1,NLMIN,-1
            CALL ZTIME(TTT0)
            ISETLV=2
            CALL  SETLEV (ISETLV)

C           Restrict the solution from the finer level to the coarser
C           one. The solution vector of the finest level always exists
C           as that is the iteration vector of the algorithm.
C           Now transported this to the lower levels, so we have a
C           representation of the solution on each level.

            I1=ILEV+1
            KU1F=L(KLUP(I1))
            KU2F=KU1F+KNU(I1)
            KVERTF=L(KLVERT(I1))
            KMIDF =L(KLMID (I1))
            KADJF =L(KLADJ (I1))
            CALL  RESTRU (DWORK(KU1),DWORK(KU2),DWORK(KU1F),DWORK(KU2F),
     *                    KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *                    NU,NP,NVT,KWORK(KVERTF),KWORK(KMIDF),
     *                    KWORK(KADJF),KNU(I1),KNP(I1),KNVT(I1) )

            CALL ZTIME(TTT1)
            TTLC=TTLC+TTT1-TTT0

C           Construct the linear part of the nonlinear matrix on the lower
C           level. Use the restricted solution from above to construct
C           the nonlinearity.
C           There is no defect vector to be calculated, we use a 
C           routine here that only builds the nonlinearity of the 
C           KA1-matrix.

            CALL ZTIME(TTT0)
            CALL XMADF3(KM1,KST1,KA1,KCOLA,KLDA,NA,NU,THSTEP,ISTAT)
            CALL ZTIME(TTT1)
            TTADF=TTADF+TTT1-TTT0

C           In the case that we have to calculate Navier Stokes (ISTOK=0)
C           or that we have a Stokes calculation with IPRECA=4,
C           add the convective part u*grad(u) to the system matrix KA1.
C
C           Up to now our system matrix has consisted only of linear
C           terms:
C
C               KA1 = [  M*I  +  THSTEP * (-nu * Laplace(.))  ]
C
C           but this is not enough... the nonlinear part is still missing.
C           So we have to add the following term to KA1:
C
C               THSTEP * u grad(.)
C
C           what will finally result in the system matrix
C
C               KA1 = [  M*I  +  THSTEP * (-nu * Laplace(.))  ] + THSTEP * u grad(.)
C                   = [  M*I  +  THSTEP * N(u) ]

            CALL ZTIME(TTT0)
            IF ((ISTOK.NE.1).OR.((ISTOK.EQ.1).AND.(IPRECA.EQ.4))) THEN
              IF (IUPW.EQ.1) THEN
               CALL GUPWD (DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                     1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                     DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                     DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *                     KWORK(L(LVERT)),KWORK(L(LMID)),
     *                     DWORK(L(LCORVG)),0)
              ELSE
               IF (IELT.EQ.0) 
     *          CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                      1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                      DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                      DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                      KWORK(L(LVERT)),KWORK(L(LMID)),
     *                      DWORK(L(LCORVG)),E031,COEFFN,0,1D0)
               IF (IELT.EQ.1) 
     *          CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                      1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                      DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                      DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                      KWORK(L(LVERT)),KWORK(L(LMID)),
     *                      DWORK(L(LCORVG)),E030,COEFFN,0,1D0)
               IF (IELT.EQ.2) 
     *          CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                      1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                      DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                      DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                      KWORK(L(LVERT)),KWORK(L(LMID)),
     *                      DWORK(L(LCORVG)),EM31,COEFFN,0,1D0)
               IF (IELT.EQ.3) 
     *          CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                      1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                      DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                      DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                      KWORK(L(LVERT)),KWORK(L(LMID)),
     *                      DWORK(L(LCORVG)),EM30,COEFFN,0,1D0)
              ENDIF
            ENDIF
            CALL ZTIME(TTT1)
            TTUPW=TTUPW+TTT1-TTT0

C           *MK1*: Call matrix restriction routine to rebuild the matrix
C           entrys belonging to anisotropic cells. This procedure is 
C           applied to all matrices except for the finest level

            IF (ILEV.NE.NLMAX) THEN
              CALL ZTIME(TTT0)
              
              CALL IMPMRS (ILEV)

              CALL ZTIME(TTT1)
              TTUPW=TTUPW+TTT1-TTT0
            END IF

C           Because the coarse grid matrix is constructed with the 
C           help of the fine grid matrix, the implementation of the 
C           boundary conditions into the matrix can't be done here! 
C           This is done later in a second step...  

          END DO

C         Implement Dirichlet boundary conditions the matrices on all
C         levels. Replace the rows corresponding to Dirichlet nodes
C         by unit vectors.

          DO ILEV=NLMAX-1,NLMIN,-1
            CALL ZTIME(TTT0)
            ISETLV=2
            CALL SETLEV (ISETLV)
            CALL ZTIME(TTT0)
            CALL BDRYA (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *                  KWORK(L(KLMBD(ILEV))),KNMBD(ILEV))

            CALL ZTIME(TTT1)
            TTBDR=TTBDR+TTT1-TTT0
          END DO

        ENDIF

C=======================================================================
C     Ok, small summary - what have we calculated up to now and where?
C     (cf. p. 179 (164) in Turek's book)
C
C     KD1,KD2,KDP  ->  currently contains d^l = T(u^l)-f,
C                      was calculated with XMADF2 at the beginning.
C     KF1,KF2,KFP  ->  Right hand side of the equation, was
C                      provided by the caller
C     KST1,KB1,KB2 ->  Matrix for the Oseen equation = linearized
C                      Navier Stokes
C                      (denoted as T~(u^l) in the book).
C                      This is now used as a preconditioner for
C                      the defect(=update) vector.
C
C     Remember, the basic nonlinear iteration, given by
C
C       u^(l+1)  =  u^l  -  OMEGA * T~(u^l)^(-1) * ( T(u^l)u^l - f )
C
C     - with the preconditioner C = T~(u^l)^{-1} - can be decomposed 
C     into three steps:
C
C      1.) Calculate the defect (done above):
C             d^l = T(u^l) u^l - f
C
C      2.) Solve an auxiliary problem (Oseen equation):
C             T~(u^l) y^l = d^l
C
C      3.) Update the solution:
C             u^(l+1) = u^l - Omega * y^l
C
C     So we are now at step 2 - we must solve the Oseen equation
C     to apply the preconditioner T~(u^l)^{-1} to the defect vector.
C
C     Remark: The precontidioner T~(u^l) was set up the same way
C      as the system matrix T(u^l) that was used for the defect
C      correction. Theoretically (cf. p. 168) the preconditioner
C      can be
C
C       [ ALPHA*M + nu*L + K~(u^l)                          B1 ] ^ {-1}
C   C = [                          ALPHA*M + nu*L + K~(u^l) B2 ]
C       [        B1^T                      B2^T             0  ]
C
C      with K~ not necessarily build up the same way as K in the system
C      matrix
C
C       [ ALPHA*M + nu*L + K(u^l)                          B1 ] ^ {-1}
C   T = [                         ALPHA*M + nu*L + K~(u^l) B2 ]
C       [        B1^T                     B2^T             0  ]
C
C       [ KA1       B1 ]
C     = [      KA1  B1 ]
C       [ B1^T B2^T 0  ]
C
C      For T, K must be as exact as possible, so theat the approximate
C      solution converges to the correct solution - but not for C.
C      Here we theoretically could use another scheme for setting up
C      the convective part.
C      Nevertheless, K~ := K is the natural choice - and as we already
C      built it for the defect correction, we can also use it here,
C      also to save some time.
C=======================================================================

C=======================================================================
C *** Perform MG for Oseen equations
C=======================================================================

C       Set finest level NLMAX to update NU,NP,...

        ISETLV=2
        ILEV=NLMAX
        CALL  SETLEV (ISETLV)

C       At this point, the linearized nonlinear system matrix
C       is set up, so we have an stationary Oseen-equation here:
C
C       [ ALPHA*M + THETA*(-nu Laplace(.)) + u^n*grad(.) ] y = (KD1,KD2,KDP)
C                                            ^^^^^^^^^^^   ^
C                                          Linearization   to calc.
C
C       This is just a linear system of the form Au=b, which we can
C       solve with our multigrid solver - so call that.
C
C       Note that solving the linear system will destroy our current
C       solution. Fortunately we normally make a backup of that 
C       (in the if-clause (LU1OLD.NE.0) above), so we can restore it
C       later for the defect correction.
C       To make the starting vector of the MG iteration not represent
C       a totally wrong vector, we fill it with 0:

        CALL LCL1 (DWORK(L(KLUP(NLMAX))),KNEQ(NLMAX))
        
C       Make a backup of the RHS of the nonlinear system and replace
C       it by the defect vector. Restore the RHS of the nonlinear
C       system later.

        CALL ZNEW (2*NU+NP,-1,LDEF,'DDEF  ')
        IF (IER.NE.0) GOTO 99999
  
        CALL LCP1(DWORK(KF1)   ,DWORK(L(LDEF)),NU)
        CALL LCP1(DWORK(KF2)   ,DWORK(L(LDEF)+NU),NU)
        CALL LCP1(DWORK(KFP)   ,DWORK(L(LDEF)+2*NU),NP)
        CALL LCP1(DWORK(L(LD1)),DWORK(KF1)    ,NU)
        CALL LCP1(DWORK(L(LD1)+NU),DWORK(KF1+NU) ,NU)
        CALL LCP1(DWORK(L(LD1)+2*NU),DWORK(KF1+2*NU),NP)

C       When ILMAX=ILMIN, suppress check of residuals for convergence
C       and suppress output of MG iteration; always perform ILMIN=ILMAX
C       MG-iterations, regardless of the defect vectror.

        IRELMG=1
        ITMG=ILMIN
        IF (ILMAX.GT.ILMIN) THEN
          IDEFMG=1
        ELSE
          IDEFMG=0
        ENDIF

C       Solve the Oseen equation to obtain ( Y1, Y2, YP )^T.
C       This is saved in the solution vectors KU1,KU2,KP on the
C       finest level.

C      CALL  M011 (DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
C     *            KNEQ,ILMAX,ITMG,DMPMG,EPSMG,
C     *            YAX,YPROL,YREST,YSM,YSM,YEX,YEXA,YDBC,YSTEP,
C     *            KIT0,KIT,IRELMG,IDEFMG,RHOLMG,BMG)
        MT = MT-1
        CALL  M011 (DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
     *            KNEQ,ILMAX,ITMG,DMPMG,EPSMG,RESOUT,
     *            YAX,YPROL,YREST,YSM,YSM,YEX,YEXA,YDBC,YSTEP,
     *            .FALSE.,I000,
     *            KIT0,KIT,IRELMG,IDEFMG,RHOLMG,RHOASM,BMG,.TRUE.)
        MT = MT+1

C       This overwrote KU1/KU2/KUP with the update vector,
C       arising from solving Oseen with the defect as a RHS.
C       The vector Y=(KU1,KU2,KP) has to be added to our solution 
C       vector. The current solution vector is still saved in the
C       array of LU1OLD/LU2OLD/LPOLD.

C *MK*!!!: Change .TRUE.->.FALSE in the above line of call to M011 !!!

        NMG=NMG+ITMG

C       That's it, we have our new iterate u^new inside of the
C       nonlinear iteration. If INLMAX=1, this is even our
C       solution vector u^(n+1).

C=======================================================================
C    End of multigrid iteration
C=======================================================================

C       Set finest level NLMAX to update NU,NP,...

        ISETLV=2
        ILEV=NLMAX
        CALL  SETLEV (ISETLV)
        
C       Restore the RHS-vector on the finest level and release
C       the memory used

        CALL LCP1(DWORK(L(LDEF)),DWORK(KF1),NU)
        CALL LCP1(DWORK(L(LDEF)+NU),DWORK(KF2),NU)
        CALL LCP1(DWORK(L(LDEF)+2*NU),DWORK(KFP),NP)

        CALL ZDISP(0,LDEF,'DDEF  ')
        IF (IER.NE.0) GOTO 99999

C=======================================================================
C *** Calculate the optimal correction
C=======================================================================

C       Starting with the OMEGA-parameter, calculate a new OMEGA
C       for the defect correction and add the update-vector
C       Y=(KU1,KU2,KP) to the current solution (KU1OLD,KU2OLD,KPOLD)
C       to obtain a new iterate (KU1,KU2,KP) of the nonlinear iteration.
C       So the following call updates KU1,KU2,KP and OMEGA.
C       As a side effect (!) the system matrix KA1 is rebuild.

        CALL  ZTIME(TTT0)
        CALL  XOPTCN (KU1,KU2,KP,L(LU1OLD),L(LU1OLD)+NU,L(LU1OLD)+2*NU,
     *              KF1,KF2,KFP,
     *              L(LD1),L(LD1)+NU,L(LD1)+2*NU,KAUX1,KAUX2,KAUXP,KA1,
     *              KCOLA,
     *              KLDA,KB1,KB2,KCOLB,KLDB,KST1,KM1,NA,NU,NP,DELU,
     *              DELP,OMEGA,KWORK(L(KLMBD(ILEV))),KNMBD(ILEV),INEUM)

C       The LD1,LD2,LDP-vectors are now written with undefined data.
C       Now we have our new iterate (KU1,KU2,KP)=u^(l+1) - hopefully...

C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================

        IF (BMG.AND.(ABS(OMEGA).GE.1D-2).AND.(ISTAT.EQ.1)) THEN
          WRITE (MTERM,'(A)') 'Calculation canceled, MG-solver broke '//
     *                        'down'
          WRITE (MFILE,'(A)') 'Calculation canceled, MG-solver broke '//
     *                        'down'
          BSTOP=.TRUE.
        END IF

C=======================================================================
C *** Calculation of defects
C=======================================================================

C       We only calculate defects if we perform real nonlinear
C       iteration...

        IF (INLMAX.GT.1) THEN
        
C         Copy the RHS of the system to the defect vector, identified
C         by the handle LD1:
        
          CALL ZTIME(TTT0)
          CALL LCP1 (DWORK(KF1),DWORK(L(LD1)),NU)
          CALL LCP1 (DWORK(KF2),DWORK(L(LD1)+NU),NU)
          CALL LCP1 (DWORK(KFP),DWORK(L(LD1)+2*NU),NP)
          CALL ZTIME(TTT1)
          TTLC=TTLC+TTT1-TTT0

C         Build the linear part of the system matrix into KST1.
C         Calculate the linear part of the defect vector into LD1:

          CALL ZTIME(TTT0)
          CALL XMADF2(KM1,KST1,KA1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                L(LD1),L(LD1)+NU,L(LD1)+2*NU,KU1,KU2,KP,NA,NU,NP,
     *                KWORK(L(KLMBD(NLEV))),KNMBD(NLEV),INEUM,THSTEP,
     *                ISTAT)
          CALL ZTIME(TTT1)
          TTADF=TTADF+TTT1-TTT0

C         Incorporate the nonlinearity into the system matrix:
C
C             KA1 = [  ALPHA*M*I  +  THSTEP * (-nu * Laplace(.))  ] + THSTEP * u grad(.)
C                 = [  ALPHA*M*I  +  THSTEP * N(u) ]
C
C         (with ALPHA=ISTAT)
C
C         Incorporate the nonlinearity into the defect vector.

          CALL  ZTIME(TTT0)
          IF ((ISTOK.NE.1).OR.((ISTOK.EQ.1).AND.(IPRECA.EQ.4))) THEN
            IF (IUPW.EQ.1) THEN
              CALL GUPWD (DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                    1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                    DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                    DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *                    KWORK(L(LVERT)),KWORK(L(LMID)),
     *                    DWORK(L(LCORVG)),IDEFUP)
            ELSE
              IF (IELT.EQ.0) 
     *          CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                      1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                      DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                      DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                      KWORK(L(LVERT)),KWORK(L(LMID)),
     *                      DWORK(L(LCORVG)),E031,COEFFN,IDEFUP,1D0)
              IF (IELT.EQ.1) 
     *          CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                      1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                      DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                      DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                      KWORK(L(LVERT)),KWORK(L(LMID)),
     *                      DWORK(L(LCORVG)),E030,COEFFN,IDEFUP,1D0)
              IF (IELT.EQ.2) 
     *          CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                      1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                      DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                      DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                      KWORK(L(LVERT)),KWORK(L(LMID)),
     *                      DWORK(L(LCORVG)),EM31,COEFFN,IDEFUP,1D0)
              IF (IELT.EQ.3) 
     *          CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                      1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                      DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                      DWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                      KWORK(L(LVERT)),KWORK(L(LMID)),
     *                      DWORK(L(LCORVG)),EM30,COEFFN,IDEFUP,1D0)
            ENDIF
          ENDIF
          
          CALL ZTIME(TTT1)
          TTUPW=TTUPW+TTT1-TTT0

C         Incorporate Dirichlet conditions into the matrix and the
C         defect vector.

          CALL ZTIME(TTT0)
          CALL BDRYA  (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *                 KWORK(L(KLMBD(ILEV))),KNMBD(ILEV))

          CALL BDRY0 (DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                KWORK(L(KLMBD(ILEV))),KNMBD(ILEV))
          CALL ZTIME(TTT1)
          TTBDR=TTBDR+TTT1-TTT0

C=======================================================================
C *** Calculation of defects and norms
C=======================================================================

C         Use RESDFK to compute the norms RESU and RESDIV from
C         the solution, the RHS and the residuum:
C
C                       || (D1,D2) ||_E 
C         RESU = -----------------------------
C                max ( ||F1||_E , ||F2||_E )
C
C                     || P ||_E
C         RESDIV = ----------------
C                  || (U1,U2) ||_E
C
C         RES = || (RESU,RESDIV) ||_E


          CALL  ZTIME(TTT0)
          CALL  RESDFK(DWORK(KU1),DWORK(KU2),DWORK(KP),
     *                 DWORK(L(LD1)),DWORK(L(LD1)+NU),
     *                 DWORK(L(LD1)+2*NU),DWORK(KF1),DWORK(KF2),
     *                 DWORK(KFP),NU,NP,RESU,RESDIV)

          RES = SQRT(RESU*RESU+RESDIV*RESDIV)

C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================

          IF ( RES/RESOLD.GT.1D2) BSTOP=.TRUE.

          RESOLD=RES
          RHO   =(RES/RES0)**(1D0/DBLE(INL))
          CALL ZTIME(TTT1)
          TTLC=TTLC+TTT1-TTT0
          
        ELSE
        
C         In case of linear extrapolation, the norms are not
C         checked for divergence.
        
          RES   =0D0
          RESU  =0D0
          RESDIV=0D0
          RESOLD=0D0
          RHO   =0D0
          
        END IF ! (INLMAX.GT.1)

C=======================================================================
C *** Control of terminating the nonlinear iteration
C=======================================================================

C       Check if the nonlinear iteration can prematurely terminate.
C       In that case sei BNLENT=true - a later IF-case will leave
C       the iteration in this case.

        IF ((DELU.LE.EPSUR).AND.(DELP.LE.EPSPR)   .AND.
     *      (RESU.LE.EPSD) .AND.(RESDIV.LE.EPSDIV).AND.
     *      (RES.LE.EPSRES).AND.(INL.GE.INLMIN)) BNLEND=.TRUE.
        IF ((DELU.LE.EPSUR).AND.(DELP.LE.EPSPR)   .AND.
     *      (INLMIN.EQ.INLMAX).AND.(INLMIN.EQ.1)) BNLEND=.TRUE.

        IF (MSHOW.GE.2) 
     *    WRITE(MTERM,1003) INL,DELU,DELP,RESU,RESDIV,RES,RHO,
     *                      OMEGA,RHOLMG
        IF (MSHOW.GE.1) 
     *    WRITE(MFILE,1003) INL,DELU,DELP,RESU,RESDIV,RES,RHO,
     *                      OMEGA,RHOLMG
C
C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================

C       Cancel the iteration here if the MG-solver broke down above.

        IF (BSTOP) GOTO 221

C=======================================================================
C *** Autosave
C=======================================================================

        IF (IAUSAV.NE.0) THEN
          IF (MOD(INL,IAUSAV).EQ.0) THEN      
            CALL ZTIME(TTT0)
            CFILE='#data/#AUTOSAV   '
            CALL  OF0 (39,CFILE,0)
            CALL  OWA1 (DWORK(KU1),'DU12P ',NUP, 39,0)
            REWIND(39)
            CLOSE(39)
            CALL ZTIME(TTT1)
            TTPOST=TTPOST+TTT1-TTT0
          ENDIF
        ENDIF

C=======================================================================
C *** Return if BNLEND=true
C=======================================================================

C       Leave the nonlinear iteration if BNLEND=true

        IF ((BNLEND).OR.((INLMIN.EQ.INLMAX).AND.(INLMIN.EQ.1))
     *              .OR.(ABS(OMEGA).LT.1D-1)) GOTO 221

C       That's it with the nonlinear iteration

      END DO ! INL

C=======================================================================
C *** End of the nonlinear loop
C=======================================================================

221   IF (MSHOW.GE.4) THEN
      
        WRITE(MTERM,*) NNONL,NMG
        WRITE(MTERM,*) ' MULTIGRID COMPONENTS [in percent]:',TTMG
        WRITE(MTERM,*) ' smoothing     :', 1.D2*TTS/TTMG
        WRITE(MTERM,*) ' solver        :', 1.D2*TTE/TTMG
        WRITE(MTERM,*) ' defect calc.  :', 1.D2*TTD/TTMG
        WRITE(MTERM,*) ' prolongation  :', 1.D2*TTP/TTMG
        WRITE(MTERM,*) ' restriction   :', 1.D2*TTR/TTMG
        WRITE(MTERM,1)
      ENDIF

      IF (MSHOW.GE.3) THEN            
        WRITE(MFILE,*) NNONL,NMG
        WRITE(MFILE,*) ' MULTIGRID COMPONENTS [in percent]:',TTMG
        WRITE(MFILE,*) ' smoothing     :', 1.D2*TTS/TTMG
        WRITE(MFILE,*) ' solver        :', 1.D2*TTE/TTMG
        WRITE(MFILE,*) ' defect calc.  :', 1.D2*TTD/TTMG
        WRITE(MFILE,*) ' prolongation  :', 1.D2*TTP/TTMG
        WRITE(MFILE,*) ' restriction   :', 1.D2*TTR/TTMG
        WRITE(MFILE,1)
      ENDIF

C     In case that the maximum number of nonlinear steps were done,
C     INL = INLMAX+1 because of the DO-loop - reset INL=INLMAX
C     for output purposes.

      IF (INL.GT.INLMAX) INL=INLMAX

   1  FORMAT(80('-'))
1001  FORMAT(' IT RELU',5X,'RELP',5X,'DEF-U',4X,'DEF-DIV',2X,
     *       'DEF-TOT',2X,'RHONL ',3X,'OMEGNL',3X,'RHOMG')
1002  FORMAT(I3,18X,3(D9.2))
1003  FORMAT(I3,8(D9.2))

99999 END
