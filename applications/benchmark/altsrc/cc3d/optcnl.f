************************************************************************
      SUBROUTINE  XOPTCN (KU1,KU2,KU3,KP,KU1OLD,KU2OLD,KU3OLD,
     *                    KPOLD,KF1,KF2,KF3,KFP,KD1,KD2,KD3,KDP,
     *                    KAUX1,KAUX2,KAUX3,KAUXP,KA1,KCOLA,KLDA,
     *                    KB1,KB2,KB3,KCOLB,KLDB,KST1,KM1,NA,NU,NP,
     *                    DELU,DELP,OMEGA,KABD,NABD,INEUM)
************************************************************************
*    Purpose: - Computes for given vectors U and UOLD the optimal 
*               weighted correction OMEGA*V with V:=U-UOLD
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
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)      
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION KABD(*)
      DIMENSION VWORK(1),KWORK(1)
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGIEL/  KLINT(NNLEV)
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      COMMON /NSTIME/ TTGRID,TTPOST,TTADF,TTUPW,TTBDR,TTLC
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
      SAVE 
C
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C     E X T E R N A L S
C-----------------------------------------------------------------------
C *** Coefficient of stiffness matrix
      EXTERNAL COEFFN
C *** definition of finite elements
      EXTERNAL E030,E031,EM31,EM30
C
      EXTERNAL MATML
C
C-----------------------------------------------------------------------
C
      IF (OMGMIN.EQ.OMGMAX) THEN
       IF (OMGMIN.LT.0D0) THEN
        OMEGA=ABS(OMGMIN)
        DELU=0D0
        DELP=0D0
        GOTO 99999
       ELSE
        OMEGA=OMGMIN
        GOTO 999
       ENDIF
      ENDIF
C
C=======================================================================
C *** Calculate on AUX1/2 the linearization point: UOLD+OMEGA*(U-UOLD)
C=======================================================================
C
      CALL ZTIME(TTT0)
      CALL  LCP1 (DWORK(KU1),DWORK(KAUX1),NU)
      CALL  LCP1 (DWORK(KU2),DWORK(KAUX2),NU)
      CALL  LCP1 (DWORK(KU3),DWORK(KAUX3),NU)
      AA1=1.0D0-OMEGA
      AA2=OMEGA
      CALL  LLC1 (DWORK(KU1OLD),DWORK(KAUX1),NU,AA1,AA2)
      CALL  LLC1 (DWORK(KU2OLD),DWORK(KAUX2),NU,AA1,AA2)
      CALL  LLC1 (DWORK(KU3OLD),DWORK(KAUX3),NU,AA1,AA2)
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C=======================================================================
C *** Calculate the new nonlinear block A at the point AUX
C=======================================================================
C
      CALL ZTIME(TTT0)
      CALL XMADF3(KM1,KST1,KA1,KCOLA,KLDA,NA,NU,THSTEP,ISTAT)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      IF (IUPW.EQ.1) THEN
       CALL GUPWD  (DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),
     *              DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),1D0,0D0,
     *              DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),
     *              DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),VWORK(KA1),
     *              KWORK(KCOLA),KWORK(KLDA),KWORK(L(LVERT)),
     *              KWORK(L(LAREA)),DWORK(L(LCORVG)),0)
      ELSE
       IF (IELT.EQ.0) 
     *  CALL SUPWDG(DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),
     *              DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),1D0,0D0,
     *              DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),
     *              DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),
     *              VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *              KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *              KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),
     *              E031,COEFFN,0,1D0)
       IF (IELT.EQ.1) 
     *  CALL SUPWDG(DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),
     *              DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),1D0,0D0,
     *              DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),
     *              DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),
     *              VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *              KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *              KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),
     *              E030,COEFFN,0,1D0)
       IF (IELT.EQ.2) 
     *  CALL SUPWNP(DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),
     *              DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),1D0,0D0,
     *              DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),
     *              DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),
     *              VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *              KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *              KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),
     *              EM31,COEFFN,0,1D0)
       IF (IELT.EQ.3) 
     *  CALL SUPWNP(DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),
     *              DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),1D0,0D0,
     *              DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),
     *              DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),
     *              VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *              KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *              KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),
     *              EM30,COEFFN,0,1D0)
      ENDIF
      CALL ZTIME(TTT1)
      TTUPW=TTUPW+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      CALL BDRYA (VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),KABD,NABD)
      CALL ZTIME(TTT1)
      TTBDR=TTBDR+TTT1-TTT0
C
C *** Calculate the new defect  D=F-K*U
      CALL ZTIME(TTT0)
      CALL LCP1 (DWORK(KF1),DWORK(KD1),NU)
      CALL LCP1 (DWORK(KF2),DWORK(KD2),NU)
      CALL LCP1 (DWORK(KF3),DWORK(KD3),NU)
      CALL LCP1 (DWORK(KFP),DWORK(KDP),NP)
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      CALL  MATML(DWORK(KD1),DWORK(KD2),DWORK(KD3),DWORK(KDP),
     *             DWORK(KU1),DWORK(KU2),DWORK(KU3),DWORK(KP),-1D0,1D0,
     *             VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *             VWORK(KB1),VWORK(KB2),VWORK(KB3),KWORK(KCOLB),
     *             KWORK(KLDB),NU,NP,KABD,NABD,INEUM)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
C=======================================================================
C *** Calculate on AUX the old defect  AUX=F-K*UOLD
C=======================================================================
C
C *** Calculate the new defect  D=F-K*U
      CALL ZTIME(TTT0)
      CALL LCP1 (DWORK(KF1),DWORK(KAUX1),NU)
      CALL LCP1 (DWORK(KF2),DWORK(KAUX2),NU)
      CALL LCP1 (DWORK(KF3),DWORK(KAUX3),NU)
      CALL LCP1 (DWORK(KFP),DWORK(KAUXP),NP)
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      CALL  MATML(DWORK(KAUX1),DWORK(KAUX2),DWORK(KAUX3),DWORK(KAUXP),
     *             DWORK(KU1OLD),DWORK(KU2OLD),DWORK(KU3OLD),
     *             DWORK(KPOLD),-1D0,1D0,VWORK(KA1),KWORK(KCOLA),
     *             KWORK(KLDA),VWORK(KB1),VWORK(KB2),VWORK(KB3),
     *             KWORK(KCOLB),KWORK(KLDB),NU,NP,KABD,NABD,INEUM)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
C=======================================================================
C *** Calculate   D:=AUX-D  (this is K*V with V=U-UOLD)
C=======================================================================
C
      CALL ZTIME(TTT0)
      AA1= 1.D0
      AA2=-1.D0
      CALL  LLC1 (DWORK(KAUX1),DWORK(KD1),NU,AA1,AA2)
      CALL  LLC1 (DWORK(KAUX2),DWORK(KD2),NU,AA1,AA2)
      CALL  LLC1 (DWORK(KAUX3),DWORK(KD3),NU,AA1,AA2)
      CALL  LLC1 (DWORK(KAUXP),DWORK(KDP),NP,AA1,AA2)
C
C=======================================================================
C *** Calculate   SKVKV:= (K*V,K*V) = (D,D)
C=======================================================================
C
      CALL  LSP1 (DWORK(KD1),DWORK(KD1),NU, SKVKV)
      CALL  LSP1 (DWORK(KD2),DWORK(KD2),NU, SSS)
      CALL  LSP1 (DWORK(KD3),DWORK(KD3),NU, SSK)
      SKVKV=SKVKV+SSS+SSK
      CALL  LSP1 (DWORK(KDP),DWORK(KDP),NP, SSS)
      SKVKV=SKVKV+SSS
C
C=======================================================================
C *** Calculate   SKVD:= (K*V,D) = (D,AUX)
C=======================================================================
C
      CALL  LSP1 (DWORK(KAUX1),DWORK(KD1),NU, SKVD)
      CALL  LSP1 (DWORK(KAUX2),DWORK(KD2),NU, SSS)
      CALL  LSP1 (DWORK(KAUX3),DWORK(KD3),NU, SSK)
      SKVD=SKVD+SSS+SSK
      CALL  LSP1 (DWORK(KAUXP),DWORK(KDP),NP, SSS)
      SKVD=SKVD+SSS
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C=======================================================================
C *** Calculate the optimal relaxation parameter OMEGA
C=======================================================================
C
      IF (SKVKV.LT. 1.0D-40) THEN
          WRITE(MTERM,*) 'ERROR in OPTCOR: SKVKV is nearly zero'
          STOP
      ENDIF
C
      OMEGA=SKVD/SKVKV
      IF (OMEGA.LT.ABS(OMGMIN)) OMEGA=ABS(OMGMIN)
      IF (OMEGA.GT.ABS(OMGMAX)) OMEGA=ABS(OMGMAX)
C
C=======================================================================
C *** Calculate the optimal correction  U:=OMEGA*(U-UOLD)
C=======================================================================
C
999   CALL ZTIME(TTT0)
      AA1= -OMEGA
      AA2=  OMEGA
      CALL  LLC1 (DWORK(KU1OLD),DWORK(KU1),NU,AA1,AA2)
      CALL  LLC1 (DWORK(KU2OLD),DWORK(KU2),NU,AA1,AA2)
      CALL  LLC1 (DWORK(KU3OLD),DWORK(KU3),NU,AA1,AA2)
      CALL  LLC1 (DWORK(KPOLD) ,DWORK(KP) ,NP,AA1,AA2)
C
C=======================================================================
C *** Calculate maximum changes   DELU
C=======================================================================
C
      CALL  LLI1 (DWORK(KU1),NU,DELU1,INDU1)
      CALL  LLI1 (DWORK(KU2),NU,DELU2,INDU2)
      CALL  LLI1 (DWORK(KU3),NU,DELU3,INDU3)
      CALL  LLI1 (DWORK(KP) ,NP,DELP ,INDP )
C
C=======================================================================
C *** Update the solution   U:=UOLD+U
C=======================================================================
C
      AA1= 1.D0
      AA2= 1.D0
      CALL  LLC1 (DWORK(KU1OLD),DWORK(KU1),NU,AA1,AA2)
      CALL  LLC1 (DWORK(KU2OLD),DWORK(KU2),NU,AA1,AA2)
      CALL  LLC1 (DWORK(KU3OLD),DWORK(KU3),NU,AA1,AA2)
      CALL  LLC1 (DWORK(KPOLD) ,DWORK(KP) ,NP,AA1,AA2)
C
C=======================================================================
C *** relative maximum changes   DELU
C=======================================================================
C
      CALL  LLI1 (DWORK(KU1),NU,DELT1,INDT1)
      CALL  LLI1 (DWORK(KU2),NU,DELT2,INDT2)
      CALL  LLI1 (DWORK(KU3),NU,DELT3,INDT3)
      DELT=MAX(DELT1,DELT2,DELT3)
      IF (ABS(DELT).LT.1D-8) DELT=1D0
      DELU=MAX(DELU1,DELU2,DELU3)/DELT
C
      CALL  LLI1 (DWORK(KP) ,NP,DELT,INDT)
      IF (ABS(DELT).LT.1D-8) DELT=1D0
      DELP=DELP/DELT
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C
99999 END
