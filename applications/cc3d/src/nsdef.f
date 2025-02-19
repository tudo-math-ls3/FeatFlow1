************************************************************************
      SUBROUTINE  NSDEF (MFILE,MSHOW,BSTOP,BNLEND)  
************************************************************************
*   Purpose: - solver for the stationary incompressible Navier Stokes
*              equations via fixed point defect correction plus
*              imultigrid for linear Oseen problems
*            - nonlinear version:
*                   - fixed point defect correction as outer iteration
*                   - mg as solver for the linear auxiliary Oseen
*                     problems
*                   - nonlinear parameter optimization for the 
*                     correction from the linear solver step
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299,NNLEV=9,NNAB=21,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120,CFILE*60
C
C *** Arrays for multigrid modul M011 
      DIMENSION  KOFFX(NNLEV),KOFFD(NNLEV),KOFFB(NNLEV),KNEQ(NNLEV),
     *           KIT(NNLEV),KIT0(NNLEV)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
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
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
      COMMON /NSTIME/ TTGRID,TTPOST,TTADF,TTUPW,TTBDR,TTLC
      COMMON /NSCOUN/ NNONL,NMG
      COMMON /NSEXL/  ITEXL,LTML,TIML11,TIML12,TIML31,TIML32
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTIME/ TTMG,TTS,TTE,TTD,TTP,TTR,IMTIME
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
      COMMON /MGIEL/  KLINT(NNLEV)
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
      SAVE
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C     E X T E R N A L S
C-----------------------------------------------------------------------
C *** Coefficient of stiffness matrix
      EXTERNAL COEFFN
C *** definition for linear matrix generation
      EXTERNAL E030,E031,EM30,EM31
C *** Multigrid components
      EXTERNAL  YAX,YPROL,YREST,YSM,YEX,YEXA,YDBC,YSTEP
C
C=======================================================================
C     Initialization
C=======================================================================
C
C *** Initialization of the offset arrays KOFFX,KOFFB,KOFFD and KNEQ
      DO 10 ILEV=NLMIN,NLMAX
      KOFFX(ILEV)=L(KLUP(ILEV))-1
      KOFFB(ILEV)=L(KLF12P(ILEV))-1
      KOFFD(ILEV)=L(KLAUX(ILEV))-1
      KNEQ(ILEV)=KNUP(ILEV)
10    CONTINUE
C
C=======================================================================
C     First generation of the nonlinear block A on level NLMAX
C=======================================================================
C
      CALL ZTIME(TTT0)
      BSTOP =.FALSE.
      BNLEND=.FALSE.
C
      ISETLV=2
      ILEV=NLMAX
      CALL  SETLEV (ISETLV)
C
      IF (INLMAX.GT.1) THEN
       IDEFUP=1
      ELSE
       IDEFUP=0
      ENDIF
C
      IF (INLMAX.GT.1) THEN
       CALL LCP1 (DWORK(KF1),DWORK(L(LD1)),NU)
       CALL LCP1 (DWORK(KF2),DWORK(L(LD2)),NU)
       CALL LCP1 (DWORK(KF3),DWORK(L(LD3)),NU)
       CALL LCP1 (DWORK(KFP),DWORK(L(LDP)),NP)
      ENDIF
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      CALL XMADF2(KM1,KST1,KA1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *            L(LD1),L(LD2),L(LD3),L(LDP),KU1,KU2,KU3,KP,NA,NU,NP,
     *            KWORK(L(KLABD(NLEV))),KNABD(NLEV),INEUM,THSTEP,ISTAT)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
C
      CALL ZTIME(TTT0)
      IF ((ISTOK.NE.1).OR.((ISTOK.EQ.1).AND.(IPRECA.EQ.4))) THEN
C
       IF (ITEXL.NE.0) THEN
        IF (ITEXL.EQ.1) THEN
         DTIM2=TIML12
         DTIM1=TIML11
        ENDIF       
        IF (ITEXL.EQ.3) THEN
         DTIM2=TIML32
         DTIM1=TIML31
        ENDIF       
        IF (ITEXL.LT.0) THEN
         A1L  =0D0
         A2L  =1D0
         ITEXL=ABS(ITEXL)        
        ELSE
         A1L  =-DTIM1/DTIM2
         A2L  =(DTIM1+DTIM2)/DTIM2
        ENDIF
C
        IF (IUPW.EQ.1) THEN
         CALL GUPWD(DWORK(L(LTML)),DWORK(L(LTML)+NU),
     *              DWORK(L(LTML)+2*NU),
     *              DWORK(KU1),DWORK(KU2),DWORK(KU3),A1L,A2L,
     *              DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *              DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *              VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),KWORK(L(LVERT)),
     *              KWORK(L(LAREA)),DWORK(L(LCORVG)),IDEFUP)
        ELSE
         IF (IELT.EQ.0) 
     *    CALL SUPWDG(DWORK(L(LTML)),DWORK(L(LTML)+NU),
     *                DWORK(L(LTML)+2*NU),
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),A1L,A2L,
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *                KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),E031,
     *                COEFFN,IDEFUP,1D0)
         IF (IELT.EQ.1) 
     *    CALL SUPWDG(DWORK(L(LTML)),DWORK(L(LTML)+NU),
     *                DWORK(L(LTML)+2*NU),
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),A1L,A2L,
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *                KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),E030,
     *                COEFFN,IDEFUP,1D0)
         IF (IELT.EQ.2) 
     *    CALL SUPWNP(DWORK(L(LTML)),DWORK(L(LTML)+NU),
     *                DWORK(L(LTML)+2*NU),
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),A1L,A2L,
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *                KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),EM31,
     *                COEFFN,IDEFUP,1D0)
         IF (IELT.EQ.3) 
     *    CALL SUPWNP(DWORK(L(LTML)),DWORK(L(LTML)+NU),
     *                DWORK(L(LTML)+2*NU),
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),A1L,A2L,
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *                KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),EM30,
     *                COEFFN,IDEFUP,1D0)
        ENDIF
C
        IF (ITEXL.EQ.3) CALL LCP1(DWORK(KU1),DWORK(L(LTML)),NUP)
        IF (ITEXL.GT.0) CALL LLC1(DWORK(L(LTML)),DWORK(KU1),NUP,
     *                            A1L,A2L)
       ELSE
        IF (IUPW.EQ.1) THEN
         CALL GUPWD(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *              DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *              DWORK(KU1),DWORK(KU2),DWORK(KU3), 
     *              DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *              VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),KWORK(L(LVERT)),
     *              KWORK(L(LAREA)),DWORK(L(LCORVG)),IDEFUP)
        ELSE
         IF (IELT.EQ.0) 
     *    CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *                KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),E031,
     *                COEFFN,IDEFUP,1D0)
         IF (IELT.EQ.1) 
     *    CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *                KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),E030,
     *                COEFFN,IDEFUP,1D0)
         IF (IELT.EQ.2) 
     *    CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *                KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),EM31,
     *                COEFFN,IDEFUP,1D0)
         IF (IELT.EQ.3) 
     *    CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *                KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),EM30,
     *                COEFFN,IDEFUP,1D0)
        ENDIF
       ENDIF
C
      ENDIF
      CALL ZTIME(TTT1)
      TTUPW=TTUPW+TTT1-TTT0
C
C
      CALL ZTIME(TTT0)
      CALL BDRYA  (VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *             KWORK(L(KLABD(ILEV))),KNABD(ILEV))
C
      IF (INLMAX.GT.1) THEN
      CALL BDRY0 (DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *            KWORK(L(KLABD(ILEV))),KNABD(ILEV))
      ENDIF
      CALL ZTIME(TTT1)
      TTBDR=TTBDR+TTT1-TTT0
C
C=======================================================================
C     Calculation of initial defects
C=======================================================================
C
      CALL ZTIME(TTT0)
      IF (INLMAX.GT.1) THEN
      CALL  RESDFK(DWORK(KU1),DWORK(KU2),DWORK(KU3),DWORK(KP),
     *             DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *             DWORK(L(LDP)),DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *             DWORK(KFP),NU,NP,RESU,RESDIV)
      RESOLD=DSQRT(RESU*RESU+RESDIV*RESDIV)
      RES0=RESOLD
      RES =RESOLD
       EPSRES=DMPD*RESOLD
      ELSE
       RES   =0D0
       RESOLD=0D0
       EPSRES=0D0
      ENDIF
C
C=======================================================================
C     Solution
C=======================================================================
C
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)
C
      IF (MSHOW.GE.2) WRITE(MTERM,1001)
      IF (MSHOW.GE.1) WRITE(MFILE,1001)
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)
      INL=0
      IF (MSHOW.GE.2) WRITE(MTERM,1002)  INL,RESU,RESDIV,RESOLD
      IF (MSHOW.GE.1) WRITE(MFILE,1002)  INL,RESU,RESDIV,RESOLD
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)
C
      OMEGA=OMGINI
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C=======================================================================
C *** Loop of nonlinear iteration
C=======================================================================
C
      DO 222  INL=1,INLMAX
      NNONL=NNONL+1
C
      IF (LU1OLD.NE.0) THEN
       CALL ZTIME (TTT0)
       CALL  LCP1 (DWORK(KU1),DWORK(L(LU1OLD)),NU)
       CALL  LCP1 (DWORK(KU2),DWORK(L(LU2OLD)),NU)
       CALL  LCP1 (DWORK(KU3),DWORK(L(LU3OLD)),NU)
       CALL  LCP1 (DWORK(KP ),DWORK(L(LPOLD )),NP)
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
      ENDIF
C
C=======================================================================
C *** Generate the A blocks for all coarse levels
C=======================================================================
C
      IF (NLMAX.GT.NLMIN)  THEN
        DO 22  ILEV=NLMAX-1,NLMIN,-1
        CALL ZTIME(TTT0)
        ISETLV=2
        CALL  SETLEV (ISETLV)
C
        I1=ILEV+1
        KU1F=L(KLUP(I1))
        KU2F=KU1F+KNU(I1)
        KU3F=KU2F+KNU(I1)
        KVERTF=L(KLVERT(I1))
        KAREAF =L(KLAREA (I1))
        KADJF =L(KLADJ (I1))
        CALL  RESTRU(DWORK(KU1),DWORK(KU2),DWORK(KU3),  
     *               DWORK(KU1F),DWORK(KU2F),DWORK(KU3F),
     *               KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LADJ)),
     *               NU,NP,NVT,KWORK(KVERTF),KWORK(KAREAF),KWORK(KADJF),
     *               KNU(I1),KNP(I1),KNVT(I1),VWORK(L(KLVOL(ILEV))))
        CALL ZTIME(TTT1)
        TTLC=TTLC+TTT1-TTT0
C
        CALL ZTIME(TTT0)
        CALL XMADF3(KM1,KST1,KA1,KCOLA,KLDA,NA,NU,THSTEP,ISTAT)
        CALL ZTIME(TTT1)
        TTADF=TTADF+TTT1-TTT0
C
        CALL ZTIME(TTT0)
        IF ((ISTOK.NE.1).OR.((ISTOK.EQ.1).AND.(IPRECA.EQ.4))) THEN
         IF (IUPW.EQ.1) THEN
          CALL GUPWD(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *               DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *               DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *               DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *               VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *               KWORK(L(LVERT)),KWORK(L(LAREA)),DWORK(L(LCORVG)),0)
         ELSE
          IF (IELT.EQ.0) 
     *     CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                 DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *                 DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                 DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *                 VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                 KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *                 KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),E031,
     *                 COEFFN,0,1D0)
          IF (IELT.EQ.1) 
     *     CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                 DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *                 DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                 DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *                 VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                 KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *                 KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),E030,
     *                 COEFFN,0,1D0)
          IF (IELT.EQ.2) 
     *     CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                 DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *                 DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                 DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *                 VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                 KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *                 KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),EM31,
     *                 COEFFN,0,1D0)
          IF (IELT.EQ.3) 
     *     CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                 DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *                 DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                 DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *                 VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                 KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *                 KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),EM30,
     *                 COEFFN,0,1D0)
         ENDIF
        ENDIF
        CALL ZTIME(TTT1)
        TTUPW=TTUPW+TTT1-TTT0
C
        CALL ZTIME(TTT0)
        CALL  BDRYA (VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *               KWORK(L(KLABD(ILEV))),KNABD(ILEV))
C
        CALL ZTIME(TTT1)
        TTBDR=TTBDR+TTT1-TTT0
C
  22    CONTINUE
      ENDIF
C
C=======================================================================
C *** Perform MG for Oseen equations
C=======================================================================
C
      CALL ZTIME(TTT0)
C
      ISETLV=2
      ILEV=NLMAX
      CALL  SETLEV (ISETLV)
C
      CALL ZNEW (NUP,-1,LDEF,'DDEF  ')
      IF (IER.NE.0) GOTO 99999
C
      CALL LCP1(DWORK(KF1)   ,DWORK(L(LDEF)),NUP)
      CALL LCP1(DWORK(L(LD1)),DWORK(KF1)    ,NUP)
C
      CALL LCL1 (DWORK(KU1),NUP)
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C
      IRELMG=1
      ITMG=ILMIN
      IF (ILMAX.GT.ILMIN) THEN
       IDEFMG=1
      ELSE
       IDEFMG=0
      ENDIF
C
      CALL  M011 (DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
     *            KNEQ,ILMAX,ITMG,DMPMG,EPSMG,
     *            YAX,YPROL,YREST,YSM,YSM,YEX,YEXA,YDBC,YSTEP,
     *            KIT0,KIT,IRELMG,IDEFMG,RHOLMG,BMG)
      NMG=NMG+ITMG
C
C=======================================================================
C    End of multigrid iteration
C=======================================================================
C
C *** Set finest level NLMAX
      ISETLV=2
      ILEV=NLMAX
      CALL  SETLEV (ISETLV)
C
      CALL  ZTIME(TTT0)
C
      CALL LCP1(DWORK(L(LDEF)),DWORK(KF1),NUP)
      CALL ZDISP(0,LDEF,'DDEF  ')
      IF (IER.NE.0) GOTO 99999
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C
C=======================================================================
C *** Calculate the optimal correction
C=======================================================================
C
      CALL  ZTIME(TTT0)
      CALL  XOPTCN(KU1,KU2,KU3,KP,L(LU1OLD),L(LU2OLD),L(LU3OLD),
     *             L(LPOLD),KF1,KF2,KF3,KFP,L(LD1),L(LD2),L(LD3),L(LDP),
     *             KAUX1,KAUX2,KAUX3,KAUXP,KA1,KCOLA,KLDA,KB1,KB2,KB3,
     *             KCOLB,KLDB,KST1,KM1,NA,NU,NP,DELU,DELP,OMEGA,
     *             KWORK(L(KLABD(ILEV))),KNABD(ILEV),INEUM)
C
C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================
       IF (BMG.AND.(ABS(OMEGA).GE.1D-2).AND.(ISTAT.EQ.1)) BSTOP=.TRUE.
C
C=======================================================================
C *** Calculation of defects
C=======================================================================
C
      IF (INLMAX.GT.1) THEN
       CALL  ZTIME(TTT0)
       CALL LCP1 (DWORK(KF1),DWORK(L(LD1)),NU)
       CALL LCP1 (DWORK(KF2),DWORK(L(LD2)),NU)
       CALL LCP1 (DWORK(KF3),DWORK(L(LD3)),NU)
       CALL LCP1 (DWORK(KFP),DWORK(L(LDP)),NP)
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
C
       CALL  ZTIME(TTT0)
       CALL XMADF2(KM1,KST1,KA1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *             L(LD1),L(LD2),L(LD3),L(LDP),KU1,KU2,KU3,KP,NA,NU,NP,
     *             KWORK(L(KLABD(NLEV))),KNABD(NLEV),INEUM,THSTEP,ISTAT)
       CALL ZTIME(TTT1)
       TTADF=TTADF+TTT1-TTT0
C
       CALL  ZTIME(TTT0)
       IF ((ISTOK.NE.1).OR.((ISTOK.EQ.1).AND.(IPRECA.EQ.4))) THEN
        IF (IUPW.EQ.1) THEN
         CALL GUPWD(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *              DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *              DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *              DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *              VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),KWORK(L(LVERT)),
     *              KWORK(L(LAREA)),DWORK(L(LCORVG)),IDEFUP)
        ELSE
         IF (IELT.EQ.0) 
     *    CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *                KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),E031,
     *                COEFFN,IDEFUP,1D0)
         IF (IELT.EQ.1) 
     *    CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *                KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),E030,
     *                COEFFN,IDEFUP,1D0)
         IF (IELT.EQ.2) 
     *    CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *                KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),EM31,
     *                COEFFN,IDEFUP,1D0)
         IF (IELT.EQ.3) 
     *    CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *                DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *                DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *                KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),EM30,
     *                COEFFN,IDEFUP,1D0)
        ENDIF
       ENDIF
       CALL ZTIME(TTT1)
       TTUPW=TTUPW+TTT1-TTT0
C
       CALL ZTIME(TTT0)
       CALL BDRYA  (VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *              KWORK(L(KLABD(ILEV))),KNABD(ILEV))
C
       CALL BDRY0 (DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *             KWORK(L(KLABD(ILEV))),KNABD(ILEV))
       CALL ZTIME(TTT1)
       TTBDR=TTBDR+TTT1-TTT0
C
C=======================================================================
C *** Calculation of defects and norms
C=======================================================================
C
       CALL  ZTIME(TTT0)
       CALL  RESDFK(DWORK(KU1),DWORK(KU2),DWORK(KU3),DWORK(KP),
     *              DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LD3)),
     *              DWORK(L(LDP)),DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *              DWORK(KFP),NU,NP,RESU,RESDIV)
C
       RES   =SQRT(RESU*RESU+RESDIV*RESDIV)
C
C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================
       IF ( RES/RESOLD.GT.1D2) BSTOP=.TRUE.
C
       RESOLD=RES
       RHO   =(RES/RES0)**(1D0/DBLE(INL))
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
      ELSE
       RES   =0D0
       RESU  =0D0
       RESDIV=0D0
       RESOLD=0D0
       RHO   =0D0
      ENDIF
C
C=======================================================================
C *** Control of terminating the nonlinear iteration
C=======================================================================
C
      IF ((DELU.LE.EPSUR).AND.(DELP.LE.EPSPR)   .AND.
     *    (RESU.LE.EPSD) .AND.(RESDIV.LE.EPSDIV).AND.
     *    (RES.LE.EPSRES).AND.(INL.GE.INLMIN)) BNLEND=.TRUE.
      IF ((DELU.LE.EPSUR).AND.(DELP.LE.EPSPR)   .AND.
     *    (INLMIN.EQ.INLMAX).AND.(INLMIN.EQ.1)) BNLEND=.TRUE.
C
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,1003) INL,DELU,DELP,RESU,RESDIV,RES,RHO,OMEGA,RHOLMG
      IF (MSHOW.GE.1) 
     * WRITE(MFILE,1003) INL,DELU,DELP,RESU,RESDIV,RES,RHO,OMEGA,RHOLMG
C
C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================
       IF (BSTOP) GOTO 221
C
C=======================================================================
C *** Autosave
C=======================================================================
C
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
C
C=======================================================================
C *** Return if BNLEND=true
C=======================================================================
      IF ((BNLEND).OR.((INLMIN.EQ.INLMAX).AND.(INLMIN.EQ.1))
     *            .OR.(ABS(OMEGA).LT.1D-1)) GOTO 221
C
C
C
222   CONTINUE
C
C=======================================================================
C *** End of the nonlinear loop
C=======================================================================
C
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
C
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
C
      IF (INL.GT.INLMAX)  INL=INLMAX
C
   1  FORMAT(80('-'))
1001  FORMAT(' IT RELU',5X,'RELP',5X,'DEF-U',4X,'DEF-DIV',2X,
     *       'DEF-TOT',2X,'RHONL ',3X,'OMEGNL',3X,'RHOMG')
1002  FORMAT(I3,18X,3(D9.2))
1003  FORMAT(I3,8(D9.2))
C
C
C
99999 END
