************************************************************************
      SUBROUTINE  NSDEF (MFILE,MSHOW,BSTOP,BNLEND)  
************************************************************************
*   Purpose: - solver for the stationary Burgers-equation
*              equations via fixed point defect correction plus
*              multigrid for linear problems
*            - nonlinear version:
*                   - fixed point defect correction as outer iteration
*                   - mg as solver for the linear auxiliary
*                     problems
*                   - nonlinear parameter optimization for the 
*                     correction from the linear solver step
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299, NNLEV=9, NNAB=21, NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120,CFILE*60
C
C *** Arrays for multigrid modul M010 
      DIMENSION  KOFFX(NNLEV),KOFFD(NNLEV),KOFFB(NNLEV),KNEQ(NNLEV),
     *           KIT(NNLEV),KIT0(NNLEV)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
      COMMON /NSTIME/ TTGRID,TTPOST,TTADF,TTUPW,TTBDR,TTLC,TTILU,
     *                TTMGU,TTSU,TTEU,TTDU,TTPU,TTRU,
     *                TTMGP,TTSP,TTEP,TTDP,TTPP,TTRP
      COMMON /NSCOUN/ NNONL,NMGU,NMGP
      COMMON /NSEXL/  ITEXL,LTML,TIML11,TIML12,TIML31,TIML32
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTIME/ TTMG,TTS,TTE,TTD,TTP,TTR,IMTIME
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
      COMMON /MGILUU/ ISORTU,KLAILU(NNLEV),KLTRA1(NNLEV),KLTRA2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
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
C *** definition of finite elements
      EXTERNAL E030,E031,EM30,EM31
C *** Multigrid components
      EXTERNAL  YAXU,YPROLU,YRESTU,YSMU,YEXU,YEXAU,YDBCU,YSTEPU
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
      ENDIF
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      CALL XMADF2(KM1,KMASS1,KST1,KA1,KCOLA,KLDA,L(LD1),L(LD2),KU1,KU2,
     *            NA,NU,KWORK(L(KLMBD(NLEV))),KNMBD(NLEV),THSTEP)
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
     *              DWORK(KU1),DWORK(KU2),A1L,A2L,
     *              DWORK(KU1),DWORK(KU2),
     *              DWORK(L(LD1)),DWORK(L(LD2)),VWORK(KA1),KWORK(KCOLA),
     *              KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *              DWORK(L(LCORVG)),IDEFUP,2)
        ELSE
         IF (IELT.EQ.0) 
     *    CALL SUPWDG(DWORK(L(LTML)),DWORK(L(LTML)+NU),
     *                DWORK(KU1),DWORK(KU2),A1L,A2L,
     *                DWORK(KU1),DWORK(KU2),
     *                DWORK(L(LD1)),DWORK(L(LD2)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                E031,COEFFN,IDEFUP,2,1D0)
         IF (IELT.EQ.1) 
     *    CALL SUPWDG(DWORK(L(LTML)),DWORK(L(LTML)+NU),
     *                DWORK(KU1),DWORK(KU2),A1L,A2L,
     *                DWORK(KU1),DWORK(KU2),
     *                DWORK(L(LD1)),DWORK(L(LD2)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                E030,COEFFN,IDEFUP,2,1D0)
         IF (IELT.EQ.2) 
     *    CALL SUPWNP(DWORK(L(LTML)),DWORK(L(LTML)+NU),
     *                DWORK(KU1),DWORK(KU2),A1L,A2L,
     *                DWORK(KU1),DWORK(KU2),
     *                DWORK(L(LD1)),DWORK(L(LD2)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                EM31,COEFFN,IDEFUP,2,1D0)
         IF (IELT.EQ.3) 
     *    CALL SUPWNP(DWORK(L(LTML)),DWORK(L(LTML)+NU),
     *                DWORK(KU1),DWORK(KU2),A1L,A2L,
     *                DWORK(KU1),DWORK(KU2),
     *                DWORK(L(LD1)),DWORK(L(LD2)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                EM30,COEFFN,IDEFUP,2,1D0)
        ENDIF
C
        IF (ITEXL.EQ.3) CALL LCP1(DWORK(KU1),DWORK(L(LTML)),2*NU)
        IF (ITEXL.GT.0) CALL LLC1(DWORK(L(LTML)),DWORK(KU1),2*NU,
     *                            A1L,A2L)
       ELSE
        IF (IUPW.EQ.1) THEN
         CALL GUPWD (DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KU1),DWORK(KU2),
     *               DWORK(L(LD1)),DWORK(L(LD2)),
     *               VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *               KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *               IDEFUP,2)
        ELSE
         IF (IELT.EQ.0) 
     *    CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(L(LD1)),DWORK(L(LD2)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                E031,COEFFN,IDEFUP,2,1D0)
         IF (IELT.EQ.1) 
     *    CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(L(LD1)),DWORK(L(LD2)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                E030,COEFFN,IDEFUP,2,1D0)
         IF (IELT.EQ.2) 
     *    CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(L(LD1)),DWORK(L(LD2)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                EM31,COEFFN,IDEFUP,2,1D0)
         IF (IELT.EQ.3) 
     *    CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(L(LD1)),DWORK(L(LD2)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                EM30,COEFFN,IDEFUP,2,1D0)
        ENDIF
       ENDIF
C
      ENDIF
      CALL ZTIME(TTT1)
      TTUPW=TTUPW+TTT1-TTT0
C
C
      CALL ZTIME(TTT0)
      CALL BDRYA (VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *            KWORK(L(KLMBD(ILEV))),KNMBD(ILEV))
C
      IF (INLMAX.GT.1) THEN
       CALL BDRY0 (DWORK(L(LD1)),DWORK(L(LD2)),KWORK(L(KLMBD(ILEV))),
     *             KNMBD(ILEV))
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
       CALL  RESDFK(DWORK(L(LD1)),DWORK(L(LD2)),DWORK(KF1),DWORK(KF2),
     *              NU,RESU1,RESU2)
       RESOLD=MAX(RESU1,RESU2)
       RES0=RESOLD
       RES =RESOLD
       EPSRES=DMPUD*RESOLD
      ELSE
       RES   =0D0
       RESU1 =0D0
       RESU2 =0D0
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
      IF (MSHOW.GE.2) WRITE(MTERM,1002)  INL,RESU1,RESU2
      IF (MSHOW.GE.1) WRITE(MFILE,1002)  INL,RESU1,RESU2
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)
C
      IF ((INLMAX.GT.1).AND.(RES.LE.1D-12)) THEN
       CALL LCL1(DWORK(KU1),2*NU)
       RETURN
      ENDIF
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
       CALL ZTIME(TTT0)
       CALL  LCP1 (DWORK(KU1),DWORK(L(LU1OLD)),NU)
       CALL  LCP1 (DWORK(KU2),DWORK(L(LU2OLD)),NU)
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
        KVERTF=L(KLVERT(I1))
        KMIDF =L(KLMID (I1))
        KADJF =L(KLADJ (I1))
        CALL  RESTRU (DWORK(KU1),DWORK(KU2), DWORK(KU1F),DWORK(KU2F),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *                NU,NP,NVT,KWORK(KVERTF),KWORK(KMIDF),KWORK(KADJF),
     *                KNU(I1),KNP(I1),KNVT(I1),2)
C
        CALL ZTIME(TTT1)
        TTLC=TTLC+TTT1-TTT0
C
        CALL ZTIME(TTT0)
        CALL XMADF3(KM1,KMASS1,KST1,KA1,KCOLA,KLDA,NA,NU,THSTEP)
        CALL ZTIME(TTT1)
        TTADF=TTADF+TTT1-TTT0
C
        CALL ZTIME(TTT0)
        IF ((ISTOK.NE.1).OR.((ISTOK.EQ.1).AND.(IPRECA.EQ.4))) THEN
         IF (IUPW.EQ.1) THEN
          CALL GUPWD (DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(L(LD1)),DWORK(L(LD2)),
     *                VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                0,2)
         ELSE
          IF (IELT.EQ.0) 
     *     CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                 1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                 DWORK(L(LD1)),DWORK(L(LD2)),
     *                 VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                 KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                 E031,COEFFN,0,2,1D0)
          IF (IELT.EQ.1) 
     *      CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                 1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                 DWORK(L(LD1)),DWORK(L(LD2)),
     *                 VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                 KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                 E030,COEFFN,0,2,1D0)
          IF (IELT.EQ.2) 
     *     CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                 1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                 DWORK(L(LD1)),DWORK(L(LD2)),
     *                 VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                 KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                 EM31,COEFFN,0,2,1D0)
          IF (IELT.EQ.3) 
     *     CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                 1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                 DWORK(L(LD1)),DWORK(L(LD2)),
     *                 VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                 KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                 EM30,COEFFN,0,2,1D0)
         ENDIF
        ENDIF
        CALL ZTIME(TTT1)
        TTUPW=TTUPW+TTT1-TTT0
C
        CALL ZTIME(TTT0)
        CALL  BDRYA (VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *               KWORK(L(KLMBD(ILEV))),KNMBD(ILEV))
C
        CALL ZTIME(TTT1)
        TTBDR=TTBDR+TTT1-TTT0
C
22      CONTINUE
      ENDIF
C
C=======================================================================
C     matrix sorting ISORTU > 0
C=======================================================================
C
      CALL ZTIME(TTT0)
      IF ((ISORTU.GT.0).OR.(ISMU.EQ.4).OR.(ISLU.EQ.4)) THEN
C
       DO 11  ILEV=NLMIN,NLMAX
C
       IF (ISORTU.GT.0) THEN
        CALL ZNEW(KNA(ILEV)  ,-2,LAH  ,'VAH   ')
        CALL ZNEW(KNA(ILEV)  ,-3,LCOLH,'KCOLH ')
        CALL ZNEW(KNU(ILEV)+1,-3,LLDH ,'KLDH  ')
        IF (IER.NE.0) GOTO 99999
C
        CALL ZCPY(KLA(ILEV)   ,'VA    ',LAH  ,'VAH   ')
        CALL ZCPY(KLCOLA(ILEV),'KCOLA ',LCOLH,'KCOLH ')
        CALL ZCPY(KLLDA(ILEV) ,'KLDA  ',LLDH, 'KLDH  ')
        IF (IER.NE.0) GOTO 99999
C
        CALL MTSRTV(VWORK(L(KLA(ILEV)))   ,VWORK(L(LAH)),
     *              KWORK(L(KLCOLA(ILEV))),KWORK(L(LCOLH)),
     *              KWORK(L(KLLDA(ILEV))) ,KWORK(L(LLDH)),
     *              KWORK(L(KLTRA1(ILEV))),KWORK(L(KLTRA2(ILEV))),
     *              KNU(ILEV))
C
        CALL ZDISP(0,LLDH ,'KLDH  ')
        CALL ZDISP(0,LCOLH,'KCOLH ')
        CALL ZDISP(0,LAH  ,'DAH   ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF ((ISMU.EQ.4).OR.(ISLU.EQ.3).OR.(ISLU.EQ.4)) THEN
        CALL ZNEW(KNA(ILEV),-2,KLAILU(ILEV),'VVAILU')
        IF (IER.NE.0) GOTO 99999
        CALL ZCPY(KLA(ILEV),'VA    ',KLAILU(ILEV),'VVAILU')
        IF (IER.NE.0) GOTO 99999
C
        TOLILU=1D-12
        ALPILU=0.0D0
        INDILU=1
        CALL IFD27(VWORK(L(KLAILU(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *             KWORK(L(KLLDA(ILEV))),KNU(ILEV),INDILU,ALPILU,TOLILU)
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
11     CONTINUE
C
      ENDIF
C
      CALL ZTIME(TTT1)
      TTILU=TTILU+TTT1-TTT0
C
C=======================================================================
C *** Initialization of the offset arrays KOFFX,KOFFB,KOFFD and KNEQ
C=======================================================================
C
      DO 12  ILEV=NLMIN,NLMAX
      KOFFX(ILEV)=L(KLUP  (ILEV))-1
      KOFFB(ILEV)=L(KLF12P(ILEV))-1
      KOFFD(ILEV)=L(KLAUX (ILEV))-1
      KNEQ (ILEV)=KNU(ILEV)
      KPRSM(ILEV)=NSMU*NSMUFA**(NLMAX-ILEV)
      KPOSM(ILEV)=NSMU*NSMUFA**(NLMAX-ILEV)
12    CONTINUE
C
      ICYCLE=ICYCU
      IRELMG=1
      ITMG =ILMINU
      EPSUMG=1D99
      IF (ILMAXU.GT.ILMINU) THEN
       IDEFMG=1
      ELSE
       IDEFMG=0
      ENDIF
c	do 19 i=1,30
c19    write (*,*) i,vwork(ka1-1+i),'u'
c
      CALL  M011 (DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
     *            KNEQ,ILMAXU,ITMG,DMPUMG,EPSUMG,DEFUMG,
     *            YAXU,YPROLU,YRESTU,YSMU,YSMU,YEXU,YEXAU,YDBCU,YSTEPU,
     *            KIT0,KIT,IRELMG,IDEFMG,RHOMG1,BMGU1)
      NMGU=NMGU+ITMG
C
      TTMGU=TTMGU+TTMG
      TTSU=TTSU+TTS
      TTEU=TTEU+TTE
      TTDU=TTDU+TTD
      TTPU=TTPU+TTP
      TTRU=TTRU+TTR
C
C
      DO 13  ILEV=NLMIN,NLMAX
      KOFFX(ILEV)=L(KLUP  (ILEV))-1+KNU(ILEV)
      KOFFB(ILEV)=L(KLF12P(ILEV))-1+KNU(ILEV)
      KOFFD(ILEV)=L(KLAUX (ILEV))-1+KNU(ILEV)
      KNEQ (ILEV)=KNU(ILEV)
      KPRSM(ILEV)=NSMU*NSMUFA**(NLMAX-ILEV)
      KPOSM(ILEV)=NSMU*NSMUFA**(NLMAX-ILEV)
13    CONTINUE
C
      ICYCLE=ICYCU
      IRELMG=1
      ITMG =ILMINU
      EPSUMG=1D99
      IF (ILMAXU.GT.ILMINU) THEN
       IDEFMG=1
      ELSE
       IDEFMG=0
      ENDIF
      CALL  M011 (DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
     *            KNEQ,ILMAXU,ITMG,DMPUMG,EPSUMG,DEFUMG,
     *            YAXU,YPROLU,YRESTU,YSMU,YSMU,YEXU,YEXAU,YDBCU,YSTEPU,
     *            KIT0,KIT,IRELMG,IDEFMG,RHOMG2,BMGU2)
      NMGU=NMGU+ITMG
C
      TTMGU=TTMGU+TTMG
      TTSU=TTSU+TTS
      TTEU=TTEU+TTE
      TTDU=TTDU+TTD
      TTPU=TTPU+TTP
      TTRU=TTRU+TTR
C
      RHOLMG=MAX(RHOMG1,RHOMG2)
C
C=======================================================================
C     matrix resorting ISORTU > 0
C=======================================================================
C
      CALL ZTIME(TTT0)
      IF ((ISORTU.GT.0).OR.(ISMU.EQ.4).OR.(ISLU.EQ.3).OR.
     *    (ISLU.EQ.4)) THEN
C
       DO 14  ILEV=NLMIN,NLMAX
C
       IF (ISORTU.GT.0) THEN
        CALL ZNEW(KNA(ILEV)  ,-3,LCOLH,'KCOLH ')
        CALL ZNEW(KNU(ILEV)+1,-3,LLDH ,'KLDH  ')
        IF (IER.NE.0) GOTO 99999
C
        CALL ZCPY(KLCOLA(ILEV),'KCOLA ',LCOLH,'KCOLH ')
        CALL ZCPY(KLLDA(ILEV) ,'KLDA  ',LLDH, 'KLDH  ')
        IF (IER.NE.0) GOTO 99999
C
        CALL MTSRTR(KWORK(L(KLCOLA(ILEV))),KWORK(L(LCOLH)),
     *              KWORK(L(KLLDA(ILEV))) ,KWORK(L(LLDH)),
     *              KWORK(L(KLTRA1(ILEV))),KWORK(L(KLTRA2(ILEV))),
     *              KNU(ILEV))
C
         CALL ZDISP(0,LLDH ,'KLDH  ')
         CALL ZDISP(0,LCOLH,'KCOLH ')
         IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF ((ISMU.EQ.4).OR.(ISLU.EQ.3).OR.(ISLU.EQ.4)) THEN
        CALL ZDISP(0,KLAILU(ILEV) ,'VVAILU  ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
14     CONTINUE
C
      ENDIF
C
      CALL ZTIME(TTT1)
      TTILU=TTILU+TTT1-TTT0
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
C=======================================================================
C *** Calculate the optimal correction
C=======================================================================
C
      CALL  ZTIME(TTT0)
      CALL  XOPTCN (KU1,KU2,L(LU1OLD),L(LU2OLD),KF1,KF2,L(LD1),L(LD2),
     *              KAUX1,KAUX2,KA1,KCOLA,KLDA,KST1,KM1,KMASS1,NA,NU,
     *              DELU1,DELU2,OMEGA,KWORK(L(KLMBD(ILEV))),
     *              KNMBD(ILEV),INEUM)
C
C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================
       IF ((BMGU1.OR.BMGU2).AND.(ABS(OMEGA).GE.1D-2)) BSTOP=.TRUE.
C
C=======================================================================
C *** Calculation of defects
C=======================================================================
C
      IF (INLMAX.GT.1) THEN
       CALL  ZTIME(TTT0)
       CALL LCP1 (DWORK(KF1),DWORK(L(LD1)),NU)
       CALL LCP1 (DWORK(KF2),DWORK(L(LD2)),NU)
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
C
       CALL  ZTIME(TTT0)
       CALL XMADF2(KM1,KMASS1,KST1,KA1,KCOLA,KLDA,L(LD1),L(LD2),KU1,KU2,
     *             NA,NU,KWORK(L(KLMBD(NLEV))),KNMBD(NLEV),THSTEP)
       CALL ZTIME(TTT1)
       TTADF=TTADF+TTT1-TTT0
C
       CALL  ZTIME(TTT0)
       IF ((ISTOK.NE.1).OR.((ISTOK.EQ.1).AND.(IPRECA.EQ.4))) THEN
        IF (IUPW.EQ.1) THEN
         CALL GUPWD (DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KU1),DWORK(KU2),
     *               DWORK(L(LD1)),DWORK(L(LD2)),
     *               VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *               KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *               IDEFUP,2)
        ELSE
         IF (IELT.EQ.0) 
     *    CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(L(LD1)),DWORK(L(LD2)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                E031,COEFFN,IDEFUP,2,1D0)
         IF (IELT.EQ.1) 
     *    CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(L(LD1)),DWORK(L(LD2)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                E030,COEFFN,IDEFUP,2,1D0)
         IF (IELT.EQ.2) 
     *    CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(L(LD1)),DWORK(L(LD2)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                EM31,COEFFN,IDEFUP,2,1D0)
         IF (IELT.EQ.3) 
     *    CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(L(LD1)),DWORK(L(LD2)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                EM30,COEFFN,IDEFUP,2,1D0)
        ENDIF
       ENDIF
       CALL ZTIME(TTT1)
       TTUPW=TTUPW+TTT1-TTT0
C
       CALL ZTIME(TTT0)
       CALL BDRYA  (VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *              KWORK(L(KLMBD(ILEV))),KNMBD(ILEV))
C
       CALL BDRY0 (DWORK(L(LD1)),DWORK(L(LD2)),
     *             KWORK(L(KLMBD(ILEV))),KNMBD(ILEV))
       CALL ZTIME(TTT1)
       TTBDR=TTBDR+TTT1-TTT0
C
C=======================================================================
C *** Calculation of defects and norms
C=======================================================================
C
       CALL  ZTIME(TTT0)
       CALL  RESDFK(DWORK(L(LD1)),DWORK(L(LD2)),DWORK(KF1),DWORK(KF2),
     *              NU,RESU1,RESU2)
C
       RES=MAX(RESU1,RESU2)
C
C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================
       IF ( RES/RESOLD.GT.1D2) BSTOP=.TRUE.
C
       RESOLD=RES
       RHO   =(RES/RES0)**(1D0/DBLE(INL))
       DELU  =MAX(DELU1,DELU2)
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
      ELSE
       RES   =0D0
       RESU1 =0D0
       RESU2 =0D0
       RESOLD=0D0
       RHO   =0D0
       DELU  =MAX(DELU1,DELU2)
      ENDIF
C
C=======================================================================
C *** Control of terminating the nonlinear iteration
C=======================================================================
C
      IF ((DELU.LE.EPSUR).AND.(RES.LE.EPSUD).AND.(RES.LE.EPSRES).AND.
     *    (INL.GE.INLMIN))  BNLEND=.TRUE.
      IF ((DELU.LE.EPSUR).AND.(INLMIN.EQ.INLMAX).AND.(INLMIN.EQ.1))
     *                      BNLEND=.TRUE.
C
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,1003)  INL,DELU1,DELU2,RESU1,RESU2,RHO,OMEGA,RHOMG1,
     *                   RHOMG2
      IF (MSHOW.GE.1) 
     * WRITE(MFILE,1003)  INL,DELU1,DELU2,RESU1,RESU2,RHO,OMEGA,RHOMG1,
     *                   RHOMG2
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
       CALL  OWA1 (DWORK(KU1),'DU12P ',NUP,54,0)
       REWIND(54)
       CLOSE(54)
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
      WRITE(MTERM,*)
      WRITE(MTERM,*) ' U-MULTIGRID COMPONENTS [in percent]:',TTMGU
      WRITE(MTERM,*) ' smoothing     :', 1.D2*TTSU/TTMGU
      WRITE(MTERM,*) ' solver        :', 1.D2*TTEU/TTMGU
      WRITE(MTERM,*) ' defect calc.  :', 1.D2*TTDU/TTMGU
      WRITE(MTERM,*) ' prolongation  :', 1.D2*TTPU/TTMGU
      WRITE(MTERM,*) ' restriction   :', 1.D2*TTRU/TTMGU
      WRITE(MTERM,1)
      ENDIF
C
      IF (MSHOW.GE.3) THEN
      WRITE(MFILE,*)
      WRITE(MFILE,*) ' U-MULTIGRID COMPONENTS [in percent]:',TTMGU
      WRITE(MFILE,*) ' smoothing     :', 1.D2*TTSU/TTMGU
      WRITE(MFILE,*) ' solver        :', 1.D2*TTEU/TTMGU
      WRITE(MFILE,*) ' defect calc.  :', 1.D2*TTDU/TTMGU
      WRITE(MFILE,*) ' prolongation  :', 1.D2*TTPU/TTMGU
      WRITE(MFILE,*) ' restriction   :', 1.D2*TTRU/TTMGU
      WRITE(MFILE,1)
      ENDIF
C
      IF (INL.GT.INLMAX)  INL=INLMAX
C
   1  FORMAT(80('-'))
1001  FORMAT(' IT REL-U1',3X,'REL-U2',3X,'DEF-U1',3X,'DEF-U2',3X,
     *       'RHONL ',3X,'OMEGNL',3X,'RHOMG1',3X,'RHOMG2')
1002  FORMAT(I3,18X,2(D9.2))
1003  FORMAT(I3,9(D9.2))
C
C
C
99999 END
