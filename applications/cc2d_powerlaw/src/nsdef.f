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
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120,CFILE*60
C
C *** Arrays for multigrid modul M010 
      DIMENSION  KOFFX(NNLEV),KOFFD(NNLEV),KOFFB(NNLEV),  KNEQ(NNLEV),
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
c
      CHARACTER CPARM1*60,CMESH1*60,CFILE1*60,CSTART*60,CSOL*60
      COMMON /FILES/ IMESH1,MMESH1,CPARM1,CMESH1,MFILE1,CFILE1,
     *               ISTART,MSTART,CSTART,ISOL,MSOL,CSOL
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
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
      COMMON /NSTIME/ TTGRID,TTPOST,TTADF,TTUPW,TTBDR,TTLC
      COMMON /NSCOUN/ NNONL,NMG
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
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
c
      INCLUDE 'jump.inc'
      INCLUDE 'bouss.inc'
      INCLUDE 'block.inc'
      SAVE
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C     E X T E R N A L S
C-----------------------------------------------------------------------
C *** Parametrization of the domain
      EXTERNAL PARX,PARY,TMAX
C *** Coefficient of stiffness matrix, right hand side, exact solution
      EXTERNAL COEFFN,RHS,UE,PE,UEX,UEY
C *** definition of finite elements
      EXTERNAL E030,E031,EM30,EM31
C *** Multigrid components
c$$$      EXTERNAL  YAX,YPROL,YREST,YSM,YEX,YEXA,YDBC,YSTEP
      EXTERNAL  YAX,YPROL2,YREST2,YPROL,YREST,YSM,YEX,YEXA,YDBC,YSTEP
      EXTERNAL  YSM1
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
C
C=======================================================================
C     First generation of the nonlinear block A on level NLMAX
C=======================================================================
C
      RESOLD=1.0D0
c
      CALL ZTIME(TTT0)
      BSTOP =.FALSE.
      BNLEND=.FALSE.
C
      ISETLV=2
      ILEV=NLMAX
      CALL  SETLEV (ISETLV)
c
      IF (INLMAX.GT.1) THEN
       IDEFUP=1
      ELSE
       IDEFUP=0
      ENDIF
C
c$$$C=======================================================================
c$$$          CALL GARLI3 (ILEV,NEL,NVT,NMT,
c$$$     *             DWORK(L(LCORVG)), KWORK(L(LVERT)),KWORK(L(LMID)),
c$$$     *             DWORK(KU1),DWORK(KU2),95)
c$$$C=======================================================================
c$$$C=======================================================================
c$$$          CALL GARLI3 (ILEV,NEL,NVT,NMT,
c$$$     *             DWORK(L(LCORVG)), KWORK(L(LVERT)),KWORK(L(LMID)),
c$$$     *             DWORK(KF1),DWORK(KF2),96)
c$$$C=======================================================================
c
        CALL BDRSET (DWORK(Ku1),DWORK(Ku2),DWORK(Kf1),
     *               DWORK(Kf2),KWORK(L(KLMBD(ILEV))),
     *               DWORK(L(KLDBD(ILEV))),KWORK(L(LNPR)),
     *               KNMBD(ILEV),NVT,PARX,PARY,UE)

      IF (INLMAX.GT.1) THEN
       CALL LCP1 (DWORK(KF1),DWORK(L(LD1)),NU)
       CALL LCP1 (DWORK(KF2),DWORK(L(LD2)),NU)
       CALL LCP1 (DWORK(KFP),DWORK(L(LDP)),NP)
      ENDIF
c
      IF(RESOLD .LT. DTOLERNEW .AND. ILEV .EQ. NLEV)IUSENEWT=2
      IF(INEWTON.EQ.1.AND.INL.GE.IFIXMIN .AND. IUSENEWT .EQ. 2 
     *          .OR. INL.GE.IFIXMAX)IUSENEWT=1
c

c          write (*,*) (i,dwork(ku1-1+i),DWORK(ku2-1+i),'\n',i=1,nmt)
 
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      CALL XMADF2(KM1,KST1,KA1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *            L(LD1),L(LD2),L(LDP),KU1,KU2,KP,NA,NU,NP,
     *            KWORK(L(KLMBD(NLEV))),KNMBD(NLEV),INEUM,THSTEP,ISTAT)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
      CALL ZTIME(TTT0)
c
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
      if (igrad.eq.2) igrad=0
      if (ipreco.eq.2)ipreco=0
C
         CALL STABIL (DWORK(L(LTML)),DWORK(L(LTML)+NU),
     *                DWORK(KU1),DWORK(KU2),A1L,A2L,
     *                DWORK(KU1),DWORK(KU2),
     *                DWORK(L(LD1)),DWORK(L(LD2)),
     *                VWORK(KA1),NA,KWORK(KCOLA),
     *                KWORK(KLDA),
     *          VWORK(KB1),VWORK(KB2),NB,
     *          VWORK(KB1+NB),VWORK(KB2+NB),KWORK(KCOLB),KWORK(KLDB),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),
     *                KWORK(L(LADJ)),KWORK(L(LMEL)), 
     *                DWORK(L(LCORVG)),E031,COEFFN,IDEFUP,1D0)
C
        IF (ITEXL.EQ.3) CALL LCP1(DWORK(KU1),DWORK(L(LTML)),NUP)
        IF (ITEXL.GT.0) CALL LLC1(DWORK(L(LTML)),DWORK(KU1),NUP,
     *                            A1L,A2L)
       ELSE
C
        CALL STABIL (DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(L(LD1)),DWORK(L(LD2)),
     *                VWORK(KA1),NA,KWORK(KCOLA),
     *                KWORK(KLDA),
     *          VWORK(KB1),VWORK(KB2),NB,
     *          VWORK(KB1+NB),VWORK(KB2+NB),KWORK(KCOLB),KWORK(KLDB),
     *          KWORK(L(LVERT)),KWORK(L(LMID)),
     *                KWORK(L(LADJ)),KWORK(L(LMEL)),
     *                DWORK(L(LCORVG)),E031,COEFFN,IDEFUP,1D0) 
C
       ENDIF
      ENDIF
      CALL ZTIME(TTT1)
      TTUPW=TTUPW+TTT1-TTT0
C
c         PRINT*,'VWORK(KB1)',(VWORK(KB1+I-1),I=1,5)
c         PRINT*,'VWORK(KB1+NB)',(VWORK(KB1+NB+I-1),I=1,5)
C
      CALL ZTIME(TTT0)
      CALL BDRYA  (VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *             KWORK(L(KLMBD(ILEV))),KNMBD(ILEV))
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
      CALL  RESDFC(DWORK(KU1),DWORK(KU2),DWORK(KP),
     *             DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LDP)),
     *             RESU,RESDIV,VWORK(KM1),VWORK(L(KLAREA(ILEV))))
      IF ((ABS(ICHDEF).EQ.1).OR.(ABS(ICHDEF).EQ.3)) THEN
       RESOLD=SQRT(RESU*RESU+RESDIV*RESDIV)
      ELSE
       RESOLD=MAX(RESU,RESDIV)
      ENDIF
      RES0  =RESOLD
      RES   =RESOLD
      EPSRES=DMPD*RESOLD
c
c
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
      IF (IGRAD.EQ.1) THEN
         IF (MSHOW.GE.2) WRITE(MTERM,1002)  INL,RESU,RESDIV,RESOLD
         IF (MSHOW.GE.1) WRITE(MFILE,1002)  INL,RESU,RESDIV,RESOLD
      ELSE
         IF (MSHOW.GE.2) WRITE(MTERM,10021)  INL,RESU,RESDIV,RESOLD,
     *    IGRAD, IPRECO
         IF (MSHOW.GE.1) WRITE(MFILE,10021)  INL,RESU,RESDIV,RESOLD,
     *    IGRAD, IPRECO
      ENDIF
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
      IZBV2=INL
C
      IF (LU1OLD.NE.0) THEN
       CALL ZTIME(TTT0)
       CALL  LCP1 (DWORK(KU1),DWORK(L(LU1OLD)),NU)
       CALL  LCP1 (DWORK(KU2),DWORK(L(LU2OLD)),NU)
       CALL  LCP1 (DWORK(KP ),DWORK(L(LPOLD )),NP)
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
      ENDIF

C=======================================================================
C *** Generate the A blocks for all coarse levels
C=======================================================================
C
      IF (NLMAX.GT.NLMIN)  THEN
       DO 22  ILEV=NLMAX-1,NLMIN,-1
        CALL ZTIME(TTT0)
        ISETLV=2
        CALL  SETLEV (ISETLV)
c
        I1=ILEV+1
        KU1F=L(KLUP(I1))
        KU2F=KU1F+KNU(I1)

        IF  (MODINT.ge.0) THEN
c          IF (INL.EQ.1) write (*,*) 'klassische Rest/Prol'
        KVERTF=L(KLVERT(I1))
        KMIDF =L(KLMID (I1))
        KADJF =L(KLADJ (I1))
        CALL  RESTRU (DWORK(KU1),DWORK(KU2),DWORK(KU1F),DWORK(KU2F),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *                NU,NP,NVT,KWORK(KVERTF),KWORK(KMIDF),KWORK(KADJF),
     *                KNU(I1),KNP(I1),KNVT(I1) )
        ELSE !   (MODINT.ge.0)
c          IF (INL.EQ.1) write (*,*) 'modifizierte Rest/Prol'
         LLV2=L(KLVERT(ILEV+1))
         LLV1=L(KLVERT(ILEV))
         LLA2=L(KLADJ(ILEV+1))
         LLA1=L(KLADJ(ILEV))
         LLM2=L(KLMID(ILEV+1))
         LLM1=L(KLMID(ILEV))
         NVT2=KNVT(ILEV+1)
         NVT1=KNVT(ILEV)
         NEL2=KNEL(ILEV+1)
         NEL1=KNEL(ILEV)
         NMT2=KNMT(ILEV+1)
         NMT1=KNMT(ILEV)
C
         LLAR2=L(KLAREA(ILEV+1))
         LLAR1=L(KLAREA(ILEV))

         LASP1=L(KLASPR(ILEV))
C
         CALL RESTU2(DWORK(KU1F),DWORK(KU1),KWORK(LLV2),KWORK(LLV1),
     *               KWORK(LLM2),KWORK(LLM1),KWORK(LLA2),KWORK(LLA1),
     *               DWORK(L(LCORVG)),VWORK(LLAR1),DWORK(LASP1),
     *               NVT2,NVT1,NEL2,NEL1,NMT1,NMT2)
C
         CALL RESTU2(DWORK(KU2F),DWORK(KU2),KWORK(LLV2),KWORK(LLV1),
     *               KWORK(LLM2),KWORK(LLM1),KWORK(LLA2),KWORK(LLA1),
     *               DWORK(L(LCORVG)),VWORK(LLAR1),DWORK(LASP1),
     *               NVT2,NVT1,NEL2,NEL1,NMT1,NMT2)
c
        ENDIF!   (MODINT.ge.0)
c

        CALL ZTIME(TTT1)
        TTLC=TTLC+TTT1-TTT0
C
        CALL ZTIME(TTT0)
        CALL XMADF3(KM1,KST1,KA1,KCXSOLA,KLDA,NA,NU,THSTEP,ISTAT)
        CALL ZTIME(TTT1)
        TTADF=TTADF+TTT1-TTT0
C
        CALL ZTIME(TTT0)
        IF ((ISTOK.NE.1).OR.((ISTOK.EQ.1).AND.(IPRECA.EQ.4))) THEN
C
        CALL STABIL (DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(L(LD1)),DWORK(L(LD2)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *          VWORK(KB1),VWORK(KB2),NB,
     *          VWORK(KB1+NB),VWORK(KB2+NB),KWORK(KCOLB),KWORK(KLDB),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),
     *                KWORK(L(LADJ)),KWORK(L(LMEL)),
     *                DWORK(L(LCORVG)),E031,COEFFN,0,1D0)
C
        ENDIF
        CALL ZTIME(TTT1)
        TTUPW=TTUPW+TTT1-TTT0
C
        CALL ZTIME(TTT0)
c
        CALL  BDRYA (VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *               KWORK(L(KLMBD(ILEV))),KNMBD(ILEV))
C
4444    IF  (MODINT.le.-2) THEN
           IF (ILEV.NE.NLMAX) THEN
C
         IF (INL.EQ.1)  WRITE (*,*) 'Modifizierte Matrix-HOSSA!'
         write(*,*) 'ist noch nicht umgeschrieben!!!'
         goto 1277
c$$$           KA2   =L(KLA(ILEV+1))
c$$$           KCOLA2=L(KLCOLA(ILEV+1))
c$$$           KLDA2 =L(KLLDA(ILEV+1))
c$$$C
c$$$           CALL  MAREST(KWORK(LLV1),KWORK(LLV2),KWORK(LLM1),KWORK(LLM2),
c$$$     *               KWORK(LLA1),KWORK(LLA2),VWORK(KA1),VWORK(KA2),
c$$$     *               KWORK(KLDA),KWORK(KLDA2),KWORK(KCOLA),
c$$$     *               KWORK(KCOLA2),DWORK(L(LCORVG)),VWORK(LLAR1),
c$$$     *               NEL1,NEL2,NVT1,NVT2)
c$$$C

           ENDIF
        ENDIF

 1277    continue

        CALL ZTIME(TTT1)
        TTBDR=TTBDR+TTT1-TTT0
C
22      CONTINUE
      ENDIF
C

C=======================================================================
C
      DO 23  ILEV=NLMIN,NLMAX
C
      CALL ZTIME(TTT0)
      ISETLV=2
      CALL  SETLEV (ISETLV)
C
      CALL ZTIME(TTT0)
      CALL  BDRYA (VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *             KWORK(L(KLMBD(ILEV))),KNMBD(ILEV))
C
      CALL ZTIME(TTT1)
      TTBDR=TTBDR+TTT1-TTT0
C
23    CONTINUE
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

      CALL  RESDFC(DWORK(KU1),DWORK(KU2),DWORK(KP),
     *             DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LDP)),
     *             RESU,RESDIV,VWORK(KM1),VWORK(L(KLAREA(ILEV))))
c---------------------------------------------------
c     Auf dem groben Level wird die Matrix aufgehoben
c---------------------------------------------------
      BMATCL=.TRUE. !Solange 'true' wird die GG-Matrix neu berechnet
      IF ((ISL.eq.3).OR.(BLOCSL)) THEN
         CALL ZNEW(KNU(NLMIN)*KNU(NLMIN),1,LBUMU,'BUMatU')!Fuer u
         CALL ZNEW(KNU(NLMIN),3,LBUPiU,'pivot1')
         CALL ZNEW(KNP(NLMIN)*KNP(NLMIN),1,LBUMP,'BUMatP')!Fuer p
         CALL ZNEW(KNP(NLMIN),3,LBUPiP,'pivot2')
         ENDIF
c---------------------------------------------------
c---------------------------------------------------
c     Steuert, ob die LSK-Matrix gespeichert wird, oder nicht
c---------------------------------------------------
         DO 1210 IL =1,NNLEV
            BSKBLD(IL)=.TRUE.
 1210       CONTINUE
c
c$$$C=======================================================================
c$$$          CALL GARLI3 (ILEV,NEL,NVT,NMT,
c$$$     *             DWORK(L(LCORVG)), KWORK(L(LVERT)),KWORK(L(LMID)),
c$$$     *             DWORK(L(LD1)),DWORK(L(LD2)),40+inl)
c$$$C=======================================================================
c$$$       CALL BDRY0 (DWORK(L(LD1)),DWORK(L(LD2)),KWORK(L(KLMBD(ILEV))),
c$$$     *             KNMBD(ILEV))
c$$$C=======================================================================
c$$$          CALL GARLI3 (ILEV,NEL,NVT,NMT,
c$$$     *             DWORK(L(LCORVG)), KWORK(L(LVERT)),KWORK(L(LMID)),
c$$$     *             DWORK(ku1),DWORK(ku2),40+inl)
c$$$C=======================================================================
C=======================================================================
C     MG-Vorkonditionierer. Wahl: Vanca-Glaetter oder gepatcht
C=======================================================================
      IF (IBLOCK.LE.0) THEN  
C=======================================================================

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
      IRELMG=1
      ITMG=ILMIN
      IF (ILMAX.GT.ILMIN) THEN
       IDEFMG=1
      ELSE
       IDEFMG=0
      ENDIF
C
      KAREA=L(KLAREA(NLMAX))
      NP=KNP(NLMAX)
      CALL prtavg(DWORK(KP),VWORK(KAREA),NP,INEUM)
        IF (MODINT.ge.0) THEN
      CALL  M011 (DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
     *            KNEQ,ILMAX,ITMG,DMPMG,EPSMG,
     *            YAX,YPROL,YREST,YSM1,YSM1,YEX,YEXA,YDBC,YSTEP,
     *            KIT0,KIT,IRELMG,IDEFMG,RHOLMG,BMG)
c
        ELSE !   (MODINT.ge.0)
      CALL  M011 (DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
     *            KNEQ,ILMAX,ITMG,DMPMG,EPSMG,
     *            YAX,YPROL2,YREST2,YSM1,YSM1,YEX,YEXA,YDBC,YSTEP,
     *            KIT0,KIT,IRELMG,IDEFMG,RHOLMG,BMG)
         ENDIF!  (MODINT.ge.0)
      NMG=NMG+ITMG
C
      CALL prtavg(DWORK(KP),VWORK(KAREA),NP,INEUM)
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
      ELSE !gehoert zu IF (IBLOCK.LE.0)
                CALL ZNEW (NUP,-1,LDEF,'DDEF  ')
         IF (IER.NE.0) GOTO 99999
C
         CALL LCP1(DWORK(KF1)   ,DWORK(L(LDEF)),NUP)
         CALL LCP1(DWORK(L(LD1)),DWORK(KF1)    ,NUP)
C
         CALL LCL1 (DWORK(KU1),NUP)
c

         CALL XOSEE (MFILE,MSHOW,BMG,RHOLMG,INL)
C *** Set finest level NLMAX
      ISETLV=2
      ILEV=NLMAX
      CALL  SETLEV (ISETLV)
C
      CALL  ZTIME(TTT0)
         CALL LCP1(DWORK(L(LDEF)),DWORK(KF1),NUP)
         CALL ZDISP(0,LDEF,'DDEF  ')
         IF (IER.NE.0) GOTO 99999
c
      ENDIF !gehoert zu IF (IBLOCK.LE.0)
C=======================================================================
C    End of multigrid iteration
C=======================================================================
c---------------------------------------------------
      IF ((ISL.eq.3).OR.(ISL.EQ.10)) THEN
         CALL ZDISP(0,LBUMU,'BUMatU')!Fuer u
         CALL ZDISP(0,LBUPiU,'pivot1')
         CALL ZDISP(0,LBUMP,'BUMatP')!Fuer p
         CALL ZDISP(0,LBUPiP,'pivot2')
         ENDIF
c---------------------------------------------------
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C
C=======================================================================
C *** Calculate the optimal correction
C=======================================================================
C
      CALL  ZTIME(TTT0)
      CALL  XOPTCN (KU1,KU2,KP,L(LU1OLD),L(LU2OLD),L(LPOLD),KF1,KF2,KFP,
     *              L(LD1),L(LD2),L(LDP),KAUX1,KAUX2,KAUXP,KA1,KCOLA,
     *              KLDA,KB1,KB2,KCOLB,KLDB,KST1,KM1,NA,NU,NP,DELU,
     *              DELP,OMEGA,KWORK(L(KLMBD(ILEV))),KNMBD(ILEV),INEUM)
c
      CALL prtavg(DWORK(KP),VWORK(KAREA),NP,INEUM)
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
       CALL LCP1 (DWORK(KP ),DWORK(L(LPOLD )),NP)
       CALL LCP1 (DWORK(KF1),DWORK(L(LD1)),NU)
       CALL LCP1 (DWORK(KF2),DWORK(L(LD2)),NU)
       CALL LCP1 (DWORK(KFP),DWORK(L(LDP)),NP)
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
C
       CALL  ZTIME(TTT0)
       CALL XMADF2(KM1,KST1,KA1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *             L(LD1),L(LD2),L(LDP),KU1,KU2,KP,NA,NU,NP,
     *             KWORK(L(KLMBD(NLEV))),KNMBD(NLEV),INEUM,THSTEP,ISTAT)
       CALL ZTIME(TTT1)
       TTADF=TTADF+TTT1-TTT0
C
      IF(RESOLD .LT. DTOLERNEW .AND. ILEV .EQ. NLEV)IUSENEWT=2
      IF(INEWTON.EQ.1.AND.INL.GE.IFIXMIN .AND. IUSENEWT .EQ. 2 
     *          .OR. INL.GE.IFIXMAX)IUSENEWT=1
c
 

       CALL  ZTIME(TTT0)
       IF ((ISTOK.NE.1).OR.((ISTOK.EQ.1).AND.(IPRECA.EQ.4))) THEN
        CALL STABIL (DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(L(LD1)),DWORK(L(LD2)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *          VWORK(KB1),VWORK(KB2),NB,
     *          VWORK(KB1+NB),VWORK(KB2+NB),KWORK(KCOLB),KWORK(KLDB),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),
     *                KWORK(L(LADJ)),KWORK(L(LMEL)),
     *                DWORK(L(LCORVG)),E031,COEFFN,IDEFUP,1D0)
C
C
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

c$$$C=======================================================================
c$$$          CALL GARLI3 (ILEV,NEL,NVT,NMT,
c$$$     *             DWORK(L(LCORVG)), KWORK(L(LVERT)),KWORK(L(LMID)),
c$$$     *             DWORK(ku1),DWORK(ku2),50+inl)
c$$$          CALL GARLI3 (ILEV,NEL,NVT,NMT,
c$$$     *             DWORK(L(LCORVG)), KWORK(L(LVERT)),KWORK(L(LMID)),
c$$$     *             DWORK(l(ld1)),DWORK(l(ld2)),50+inl)
c$$$
c$$$         CALL BDRSET (DWORK(Ku1),DWORK(Ku2),DWORK(Ku1),
c$$$     *               DWORK(Ku2),KWORK(L(KLMBD(ILEV))),
c$$$     *               DWORK(L(KLDBD(ILEV))),KWORK(L(LNPR)),
c$$$     *               KNMBD(ILEV),NVT,PARX,PARY,UE)
c$$$        CALL BDRY0 (DWORK(L(LD1)),DWORK(L(LD2)),KWORK(L(KLMBD(ILEV))),
c$$$     *             KNMBD(ILEV))
c$$$C=======================================================================
c
       CALL  ZTIME(TTT0)
       CALL  RESDFC(DWORK(KU1),DWORK(KU2),DWORK(KP),
     *              DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LDP)),
     *              RESU,RESDIV,VWORK(KM1),VWORK(L(KLAREA(ILEV))))
       IF ((ABS(ICHDEF).EQ.1).OR.(ABS(ICHDEF).EQ.3)) THEN
        RES=SQRT(RESU*RESU+RESDIV*RESDIV)
       ELSE
        RES=MAX(RESU,RESDIV)
       ENDIF

c        write (*,*) 'residuum',res
C
c$$$
c$$$       CALL  ZTIME(TTT0)
c$$$       CALL  RESDFC(DWORK(KU1),DWORK(KU2),DWORK(KP),
c$$$     *              DWORK(L(LD1)),DWORK(L(LD2)),DWORK(L(LDP)),
c$$$     *              RESU,RESDIV,VWORK(KM1),VWORK(L(KLAREA(ILEV))))
c$$$       IF ((ABS(ICHDEF).EQ.1).OR.(ABS(ICHDEF).EQ.3)) THEN
c$$$        RESi=SQRT(RESU*RESU+RESDIV*RESDIV)
c$$$       ELSE
c$$$        RESi=MAX(RESU,RESDIV)
c$$$       ENDIF
c$$$c
c$$$       write (*,*) 'residuumNEU',resi,resu,resdiv
C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================
       RRHO= RES/RESOLD
       IF ( RES/RESOLD.GT.1D19) BSTOP=.TRUE.
C
      IF(IUSENEWT .EQ. 1 .and. RRHO .gt. 1d-3)THEN
         DMPD=(RRHO**2)
         DMPMG=(RRHO**2)
         EPSMG=1d-2
         DMPSL=1d-2
         EPSSL=(RRHO**2)
         IF(DMPD .LE. 1D-4)DMPD=1D-4
         IF(DMPMG .LE. 1D-4)DMPMG=1D-4
         IF(EPSSL .LE. 1D-4)EPSSL=1D-4
      ENDIF  
C
       RHOOLD=RHO
       RESOLD=RES
       RHO   =(RES/RES0)**(1D0/DBLE(INL))
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
      ELSE
       CALL  ZTIME(TTT0)
       CALL LCP1 (DWORK(KFP),DWORK(L(LDP)),NP)
       CALL XDFKD (-1D0,1D0)
       CALL  RESDFC(DWORK(KU1),DWORK(KU2),DWORK(KP),
     *              DWORK(KF1),DWORK(KF2),DWORK(KFP),
     *              RESU,RESDIV,VWORK(KM1),VWORK(L(KLAREA(ILEV))))
       CALL LCP1 (DWORK(L(LDP)),DWORK(KFP),NP)
c
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
     * WRITE(MTERM,1003) INL,DELU,DELP,RESU,RESDIV,RES,RHO,RRHO,RHOLMG
c$     * WRITE(MTERM,1003) INL,DELU,DELP,RESU,RESDIV,RES,RHO,OMEGA,RHOLMG
      IF (MSHOW.GE.1) 
     * WRITE(MFILE,1003) INL,DELU,DELP,RESU,RESDIV,RES,RHO,RRHO,RHOLMG
C$     * WRITE(MFILE,1003) INL,DELU,DELP,RESU,RESDIV,RES,RHO,OMEGA,RHOLMG
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
     *       'DEF-TOT',2X,'RHONL ',3X,'RHO',3X,'RHOMG')
1002  FORMAT(I3,18X,3(D9.2),5X,'Gradientenform')
10021 FORMAT(I3,18X,3(D9.2),5X,'Deformationstensor',2(I3))
1003  FORMAT(I3,8(D9.2))
C
C
C
99999 END
