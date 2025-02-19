************************************************************************
      SUBROUTINE PROSTP (MFILE,MSHOW,IPROJ,NITNSL,LRHS,BSTOP,BNLEND)
************************************************************************
*
*   Purpose: - performs 1 macro step of discrete projection method
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299,NNLEV=9,NNAB=21,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NBLOCF=2)
      DIMENSION BCONF(NBLOCF),ARRDF(NBLOCF),BSNGLF(NBLOCF)
      DIMENSION KFN(NBLOCF),KF(NNAB,NBLOCF),LF(NBLOCF)
      CHARACTER ARRDF*6
      DATA BCONF /.FALSE.,.FALSE./,ARRDF/'DF1   ','DF2   '/
      DATA BSNGLF/.FALSE.,.FALSE./,KFN/1,1/
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
      COMMON /NSSAV/  INSAV,INSAVN
      COMMON /NSSAVF/ DTFILM,DTFILO,DTAVS,DTAVSO,DTGMV,DTGMVO,
     *                IFUSAV,IFPSAV,IFXSAV,IAVS,IGMV,IFINIT
      COMMON /NSADAT/ TIMEMX,DTMIN,DTMAX,DTFACT,TIMEIN,EPSADI,EPSADL,
     *                EPSADU,IEPSAD,IADIN,IREPIT,IADTIM,PRDIF1,PRDIF2
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
      COMMON /MGPROJ/ KLC(NNLEV),KLCOLC(NNLEV),KLLDC(NNLEV),KNC(NNLEV)
      COMMON /MGILUP/ ISORTP,KLCILU(NNLEV),KLTRC1(NNLEV),KLTRC2(NNLEV)
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
C *** Parametrization of the domain
      EXTERNAL PARX,PARY,TMAX
C *** Coefficient of stiffness matrix, right hand side, exact solution
      EXTERNAL COEFFN,RHS,UE,PE,UEX,UEY
C *** definition of finite elements
      EXTERNAL E030,E031,EM31,EM30
C *** Multigrid components
      EXTERNAL  YAXP,YPROLP,YRESTP,YSMP,YEXP,YEXAP,YDBCP,YSTEPP
C
C=======================================================================
C
      KF(1,1)=1
      KF(1,2)=1
C
C=======================================================================
C *** Begin of nonsteady loop
C=======================================================================
C
      DO 100 ITNSL=1,NITNSL
      CALL ZTIME(TTT0)
C
      ITMOD=MOD(ITNSL-1,3)+1
C
      ISETLV=2
      ILEV=NLEV
      CALL  SETLEV (ISETLV)
C
      IF (IFRSTP.EQ.1) THEN
       IF (ITMOD.EQ.2) THEN
        TSTEPH=3D0*TSTEP*THETAP
       ELSE
        TSTEPH=3D0*TSTEP*THETA
       ENDIF
      ELSE
       TSTEPH=TSTEP
      ENDIF
      TIMENS=TIMENS+TSTEPH
C
      IF (IFRSTP.EQ.1) THEN
       IF (ITMOD.EQ.2) THEN
        TSTEPB= 3D0*TSTEP*THETAP
        THSTEP=-3D0*TSTEP*FALPHA*THETAP
        TMSTEP= 3D0*TSTEP*FALPHA*THETA
        TRSTEP= 3D0*TSTEP*FBETA *THETAP
       ELSE
        TSTEPB= 3D0*TSTEP*THETA
        THSTEP=-3D0*TSTEP*FBETA *THETA
        TMSTEP= 3D0*TSTEP*FALPHA*THETA
        TRSTEP= 3D0*TSTEP*FALPHA*THETA
       ENDIF
      ELSE
       TSTEPB= TSTEP
       THSTEP= TSTEP*(THETA-1D0)
       TMSTEP= TSTEP*THETA
       TRSTEP= TSTEP*THETA
      ENDIF
C
      IF (ITEXL.EQ.1) THEN
       TIML12=TIML11
       TIML11=TSTEPH
      ENDIF
C
      IF (ITEXL.EQ.3) THEN
       TIML32=TIML31
       TIML31=TSTEPH
      ENDIF       
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C=======================================================================
C     Variable Neumann parts
C=======================================================================
C
      IF (IBDR.GE.2) THEN
       DO 10 ILEV=NLMIN,NLMAX
C
       ISETLV=2
       CALL SETLEV(ISETLV)
C
       CALL ZTIME(TTT0)
       NVT =KNVT(ILEV)
       NMT =KNMT(ILEV)
       NEL =KNEL(ILEV)
       NVBD=KNVBD(ILEV)
C
       CALL XLCP3(KLNPRO(ILEV),KLNPR(ILEV),NVT+NMT)
       IF (IER.NE.0) GOTO 99999
C
       NMBD =0
       INEUM=0
       KNMBD(ILEV)=0
       CALL ZDISP(0,KLMBD(ILEV),'KMBD  ')
       CALL ZDISP(0,KLDBD(ILEV),'DDBD  ')
       IF (IER.NE.0) GOTO 99999
       CALL ZNEW (NVBD,3,KLMBD(ILEV),'KMBD  ')
       CALL ZNEW (NVBD,1,KLDBD(ILEV),'DDBD  ')
       IF (IER.NE.0) GOTO 99999
C
       CALL BDRNEU(KWORK(L(KLMBD(ILEV))),KWORK(L(KLVBD(ILEV))),
     *             KWORK(L(KLEBD(ILEV))),KWORK(L(KLVERT(ILEV))),
     *             KWORK(L(KLMID(ILEV))),KWORK(L(KLNPR(ILEV))),
     *             DWORK(L(KLDBD(ILEV))),DWORK(L(KLMBDP(ILEV))),
     *             NMBD,INEUM)
       KNMBD(ILEV)=NMBD
       CALL ZTIME(TTT1)
       TTBDR=TTBDR+TTT1-TTT0
C
C
C
       ISETLV=2
       CALL SETLEV(ISETLV)
       CALL ZTIME(TTT0)
C
C=======================================================================
C     matrix resorting for ISORTP > 0
C=======================================================================
C
       IF (ISORTP.GT.0) THEN
        CALL ZNEW(KNC(ILEV)  ,-3,LCOLH,'KCOLH ')
        CALL ZNEW(KNP(ILEV)+1,-3,LLDH ,'KLDH  ')
        IF (IER.NE.0) GOTO 99999
C
        CALL ZCPY(KLCOLC(ILEV),'KCOLC ',LCOLH,'KCOLH ')
        CALL ZCPY(KLLDC(ILEV) ,'KLDC  ',LLDH, 'KLDH  ')
        IF (IER.NE.0) GOTO 99999
C
        CALL MTSRTR(KWORK(L(KLCOLC(ILEV))),KWORK(L(LCOLH)),
     *              KWORK(L(KLLDC(ILEV))) ,KWORK(L(LLDH)),
     *              KWORK(L(KLTRC1(ILEV))),KWORK(L(KLTRC2(ILEV))),
     *              KNP(ILEV))
C
         CALL ZDISP(0,LLDH ,'KLDH  ')
         CALL ZDISP(0,LCOLH,'KCOLH ')
         IF (IER.NE.0) GOTO 99999
       ENDIF
C
C=======================================================================
C     P matrix generation 
C=======================================================================
C
       CALL XLCL1(KLC(ILEV),KNC(ILEV))
       IF (IER.NE.0) GOTO 99999
C
       CALL PROJMA(DWORK(L(KLC(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *             KWORK(L(KLLDC(ILEV))),KWORK(L(LNPR)),
     *             KWORK(L(LMID)),KWORK(L(LADJ)),
     *             DWORK(KM1),DWORK(KB1),DWORK(KB2),
     *             KWORK(KCOLB),KWORK(KLDB),NEL,NVT,NMT)
C
C=======================================================================
C     ILU(C) matrix generation + grid sorting for P
C=======================================================================
C
       IF ((ISMP.EQ.4).OR.(ISLP.EQ.3).OR.(ISLP.EQ.4).OR.
     *     (ISORTP.GT.0)) THEN
C
C=======================================================================
C     matrix sorting for ISORTP > 0
C=======================================================================
C
        IF (ISORTP.GT.0) THEN
C
         CALL ZNEW(KNC(ILEV)   ,-1,LCH  ,'DCH   ')
         CALL ZNEW(KNC(ILEV)   ,-3,LCOLH,'KCOLH ')
         CALL ZNEW(KNEL(ILEV)+1,-3,LLDH ,'KLDH  ')
         IF (IER.NE.0) GOTO 99999
C
         CALL ZCPY(KLC(ILEV)   ,'CC    ',LCH,  'DCH   ')
         CALL ZCPY(KLCOLC(ILEV),'KCOLC ',LCOLH,'KCOLH ')
         CALL ZCPY(KLLDC(ILEV) ,'KLDC  ',LLDH, 'KLDH  ')
         IF (IER.NE.0) GOTO 99999
C
         CALL MTSRTD(DWORK(L(KLC(ILEV)))   ,DWORK(L(LCH)),
     *               KWORK(L(KLCOLC(ILEV))),KWORK(L(LCOLH)),
     *               KWORK(L(KLLDC(ILEV))) ,KWORK(L(LLDH)),
     *               KWORK(L(KLTRC1(ILEV))),KWORK(L(KLTRC2(ILEV))),
     *               KNEL(ILEV))
C
         CALL ZDISP(0,LLDH ,'KLDH  ')
         CALL ZDISP(0,LCOLH,'KCOLH ')
         CALL ZDISP(0,LCH  ,'DCH   ')
         IF (IER.NE.0) GOTO 99999
        ENDIF
C
C=======================================================================
C     ILU(C) matrix
C=======================================================================
C
        IF ((ISMP.EQ.4).OR.(ISLP.EQ.3).OR.(ISLP.EQ.4)) THEN
         CALL ZCTYPE(1,KLCILU(ILEV),'DDCILU')
         IF (IER.NE.0) GOTO 99999
         CALL ZCPY(KLC(ILEV),'CC    ',KLCILU(ILEV),'DDCILU')
         IF (IER.NE.0) GOTO 99999
C
         TOLILU=1D-12
         ALPILU=0D0
         INDILU=1
         CALL IFD17(DWORK(L(KLCILU(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *              KWORK(L(KLLDC(ILEV))),NEL,INDILU,ALPILU,TOLILU)
         IF (IER.NE.0) GOTO 99999
C
         CALL ZCTYPE(2,KLCILU(ILEV),'DDCILU')
         IF (IER.NE.0) GOTO 99999
        ENDIF
C
       ENDIF
C
       CALL ZTIME(TTT1)
       TTADF=TTADF+TTT1-TTT0
C
10     CONTINUE
      ENDIF
C      
C***********************************************************************
C
      ISETLV=2
      ILEV=NLEV
      CALL SETLEV(ISETLV)
C
      CALL ZTIME(TTT0)
      IF (IRHS.EQ.0) THEN
       CALL LCL1(DWORK(KF1),NUP)
      ENDIF
C
      IF (IRHS.EQ.1) THEN
       CALL LLC1(DWORK(L(LRHS)),DWORK(KF1),NUP,TSTEPB,0D0)
      ENDIF
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      IF (IRHS.EQ.2) THEN
       CALL ZTIME(TTT0)
       CALL LLC1(DWORK(L(LRHS)),DWORK(KF1),NUP,-THSTEP,0D0)
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
C
       CALL ZTIME(TTT0)
       IF (IBDR.GE.1) THEN
        TIMENS=TIMENS-TSTEPH
        CALL PDSET  (KWORK(L(KLVBD(NLEV))),KWORK(L(KLEBD(NLEV))),
     *               KWORK(L(KLVERT(NLEV))),KWORK(L(KLMID(NLEV))),
     *               KWORK(L(KLNPR(NLEV))),DWORK(L(LCORVG)),
     *               DWORK(L(KLMBDP(NLEV))),DWORK(KF1),DWORK(KF2),
     *               -THSTEP)
        TIMENS=TIMENS+TSTEPH
       ENDIF
       CALL ZTIME(TTT1)
       TTBDR=TTBDR+TTT1-TTT0
C
       CALL ZTIME(TTT0)
       LF(1)=0
       LF(2)=0
       IF (IELT.EQ.0) 
     *  CALL  XVB0 (LF,NU,NBLOCF,ICLRF,E031,RHS,BCONF,KF,KFN,ICUBF,
     *              ARRDF,BSNGLF)
       IF (IELT.EQ.1) 
     *  CALL  XVB0 (LF,NU,NBLOCF,ICLRF,E030,RHS,BCONF,KF,KFN,ICUBF,
     *              ARRDF,BSNGLF)
       IF (IELT.EQ.2) 
     *  CALL  XVBM0(LF,NU,NBLOCF,ICLRF,EM31,RHS,BCONF,KF,KFN,ICUBF,
     *              ARRDF,BSNGLF)
       IF (IELT.EQ.3) 
     *  CALL  XVBM0(LF,NU,NBLOCF,ICLRF,EM30,RHS,BCONF,KF,KFN,ICUBF,
     *              ARRDF,BSNGLF)
       IF (IER.NE.0) GOTO 99999
       LF1=LF(1)
       LF2=LF(2)
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
C
       CALL ZTIME(TTT0)
       IF (IBDR.GE.1) THEN
        CALL PDSET  (KWORK(L(KLVBD(NLEV))),KWORK(L(KLEBD(NLEV))),
     *               KWORK(L(KLVERT(NLEV))),KWORK(L(KLMID(NLEV))),
     *               KWORK(L(KLNPR(NLEV))),DWORK(L(LCORVG)),
     *               DWORK(L(KLMBDP(NLEV))),DWORK(L(LF1)),DWORK(L(LF2)),
     *               TRSTEP)
       ENDIF
       CALL ZTIME(TTT1)
       TTBDR=TTBDR+TTT1-TTT0
C
       CALL ZTIME(TTT0)
       CALL  LLC1 (DWORK(L(LF1)),DWORK(KF1),NU,TRSTEP,1D0)
       CALL  LLC1 (DWORK(L(LF2)),DWORK(KF2),NU,TRSTEP,1D0)
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
C
       CALL ZTIME(TTT0)
       CALL BDRSET (DWORK(KF1),DWORK(KF2),DWORK(KF1),
     *              DWORK(KF2),KWORK(L(KLMBD(ILEV))),
     *              DWORK(L(KLDBD(ILEV))),KWORK(L(LNPR)),
     *              KNMBD(ILEV),NVT,PARX,PARY,UE,2)
       CALL ZTIME(TTT1)
       TTBDR=TTBDR+TTT1-TTT0
C
       CALL ZTIME(TTT0)
       CALL  LCP1 (DWORK(L(LF1)), DWORK(L(LRHS)),   NU)
       CALL  LCP1 (DWORK(L(LF2)), DWORK(L(LRHS)+NU),NU)
C
       CALL  ZDISP (0,LF1,ARRDF(1))
       CALL  ZDISP (0,LF2,ARRDF(2))
       IF (IER.NE.0) GOTO 99999
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
      ENDIF
C
C***********************************************************************
C
      CALL ZTIME(TTT0)
      ISETLV=2
      ILEV=NLEV
      CALL  SETLEV (ISETLV)
C
      IF (IBDR.GE.1) THEN
       IF ((IFRSTP.EQ.1).AND.(ITMOD.EQ.2)) TIMENS=TIMENS-TSTEPH
       CALL PDSET  (KWORK(L(KLVBD(NLEV))),KWORK(L(KLEBD(NLEV))),
     *              KWORK(L(KLVERT(NLEV))),KWORK(L(KLMID(NLEV))),
     *              KWORK(L(KLNPR(NLEV))),DWORK(L(LCORVG)),
     *              DWORK(L(KLMBDP(NLEV))),DWORK(KF1),DWORK(KF2),TSTEPB)
       IF ((IFRSTP.EQ.1).AND.(ITMOD.EQ.2)) TIMENS=TIMENS+TSTEPH
      ENDIF
      CALL ZTIME(TTT1)
      TTBDR=TTBDR+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      CALL XMADF1(KM1,KMASS1,KST1,KA1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *            KF1,KF2,KU1,KU2,KP,NA,NU,THSTEP,TSTEPH,IPROJ)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
C
      IF (THSTEP.NE.0D0) THEN
C
       CALL ZTIME(TTT0)
       THSTEP=-THSTEP
       IF (IUPW.EQ.1) THEN
        CALL GUPWD (DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *              1D0,0D0,DWORK(KU1),DWORK(KU2),
     *              DWORK(KF1),DWORK(KF2),VWORK(KA1),KWORK(KCOLA),
     *              KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *              DWORK(L(LCORVG)),2,2)
       ELSE
        IF (IELT.EQ.0) 
     *   CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KU1),DWORK(KU2),
     *               DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),E031,COEFFN,2,2,-1D0)
        IF (IELT.EQ.1) 
     *   CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KU1),DWORK(KU2),
     *               DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),E030,COEFFN,2,2,-1D0)
        IF (IELT.EQ.2) 
     *   CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KU1),DWORK(KU2),
     *               DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),EM31,COEFFN,2,2,-1D0)
        IF (IELT.EQ.3) 
     *   CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KU1),DWORK(KU2),
     *               DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),EM30,COEFFN,2,2,-1D0)
       ENDIF
       THSTEP=-THSTEP
       CALL ZTIME(TTT1)
       TTUPW=TTUPW+TTT1-TTT0
      ELSE
       IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
        CALL ZTIME(TTT0)
        TOSTEP=THSTEP
        THSTEP=1D0
        IF (IELT.EQ.0) 
     *   CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KU1),DWORK(KU2),
     *               DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),E031,COEFFN,2,2,0D0)
        IF (IELT.EQ.1) 
     *   CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KU1),DWORK(KU2),
     *               DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),E030,COEFFN,2,2,0D0)
        IF (IELT.EQ.2) 
     *   CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KU1),DWORK(KU2),
     *               DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),EM31,COEFFN,2,2,0D0)
        IF (IELT.EQ.3) 
     *   CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *               1D0,0D0,DWORK(KU1),DWORK(KU2),
     *               DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),EM30,COEFFN,2,2,0D0)
       ENDIF
       THSTEP=TOSTEP
       CALL ZTIME(TTT1)
       TTUPW=TTUPW+TTT1-TTT0
      ENDIF
C
      CALL ZTIME(TTT0)
      CALL BDRSET (DWORK(KU1),DWORK(KU2),DWORK(KF1),
     *             DWORK(KF2),KWORK(L(KLMBD(ILEV))),
     *             DWORK(L(KLDBD(ILEV))),KWORK(L(LNPR)),
     *             KNMBD(ILEV),NVT,PARX,PARY,UE,2)
      CALL ZTIME(TTT1)
      TTBDR=TTBDR+TTT1-TTT0
C
      THSTEP=TMSTEP
C
C***********************************************************************
C *** fixed point defect correction for stationary Burgers-equ.
C***********************************************************************
C
      CALL  NSDEF (MFILE,MSHOW,BSTOP,BNLEND)
      IF (IER.NE.0) GOTO 99999
C
      IF ((.NOT.BNLEND).AND.(ABS(IADTIM).GT.1)) RETURN
C
C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================
      IF (BSTOP) RETURN
C
      ILEV=NLEV
      ISETLV=2
      CALL  SETLEV (ISETLV)
C
C***********************************************************************
C *** Calculation of QR = 1/K B^T U~
C***********************************************************************
C
      CALL ZTIME(TTT0)
      CALL XDFKD(KB1,KB2,KCOLB,KLDB,KFP,KU1,KU2,NU,NP,TSTEPH)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      CALL LL21 (DWORK(KU1),NU,DNRMU1)
      CALL LL21 (DWORK(KU2),NU,DNRMU2)
      DNORMU=SQRT(DNRMU1*DNRMU1+DNRMU2*DNRMU2)
      IF (ABS(DNORMU).LT.1D-8) DNORMU=1D0
C
C***********************************************************************
C *** Solution of B^T M^-1 B P = FP
C***********************************************************************
C
      CALL LCP1(DWORK(KP)      ,DWORK(L(LDP))  ,NP)
      CALL LCP1(DWORK(L(LPOLD)),DWORK(KP)      ,NP)
      CALL LCP1(DWORK(L(LDP))  ,DWORK(L(LPOLD)),NP)
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      DO 30 ILEV=NLMIN,NLMAX
      KOFFX(ILEV)=L(KLUP  (ILEV))-1+2*KNU(ILEV)
      KOFFB(ILEV)=L(KLF12P(ILEV))-1+2*KNU(ILEV)
      KOFFD(ILEV)=L(KLAUX (ILEV))-1+2*KNU(ILEV)
      KNEQ (ILEV)=KNP(ILEV)
      KPRSM(ILEV)=NSMP*NSMPFA**(NLMAX-ILEV)
      KPOSM(ILEV)=NSMP*NSMPFA**(NLMAX-ILEV)
30    CONTINUE
C
      ICYCLE=ICYCP
      IRELMG=1
      ITMGP =ILMINP
      IF (ILMAXP.GT.ILMINP) THEN
       IDEFMG=1
      ELSE
       IDEFMG=0
      ENDIF
      CALL  M011 (DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
     *            KNEQ,ILMAXP,ITMGP,DMPPMG,DNORMU*EPSP/TSTEPH,DEFPMG,
     *            YAXP,YPROLP,YRESTP,YSMP,YSMP,YEXP,YEXAP,YDBCP,YSTEPP,
     *            KIT0,KIT,IRELMG,IDEFMG,RHOMGP,BMGP)
      NMGP=NMGP+ITMGP
C
      TTMGP=TTMGP+TTMG
      TTSP =TTSP+TTS
      TTEP =TTEP+TTE
      TTDP =TTDP+TTD
      TTPP =TTPP+TTP
      TTRP =TTRP+TTR
C
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)
      IF (MSHOW.GE.2) WRITE(MTERM,1001)
      IF (MSHOW.GE.1) WRITE(MFILE,1001)
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)
      IF (MSHOW.GE.2) WRITE(MTERM,1003) ITMGP,TSTEPH*DEFPMG/DNORMU,
     *                                  RHOMGP
      IF (MSHOW.GE.1) WRITE(MFILE,1003) ITMGP,TSTEPH*DEFPMG/DNORMU,
     *                                  RHOMGP
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)
C
      IF (MSHOW.GE.4) THEN
      WRITE(MTERM,*)
      WRITE(MTERM,*) ' P-MULTIGRID COMPONENTS [in percent]:',TTMGP
      WRITE(MTERM,*) ' smoothing     :', 1.D2*TTSP/TTMGP
      WRITE(MTERM,*) ' solver        :', 1.D2*TTEP/TTMGP
      WRITE(MTERM,*) ' defect calc.  :', 1.D2*TTDP/TTMGP
      WRITE(MTERM,*) ' prolongation  :', 1.D2*TTPP/TTMGP
      WRITE(MTERM,*) ' restriction   :', 1.D2*TTRP/TTMGP
      WRITE(MTERM,1)
      ENDIF
C
      IF (MSHOW.GE.3) THEN
      WRITE(MFILE,*)
      WRITE(MFILE,*) ' P-MULTIGRID COMPONENTS [in percent]:',TTMGP
      WRITE(MFILE,*) ' smoothing     :', 1.D2*TTSP/TTMGP
      WRITE(MFILE,*) ' solver        :', 1.D2*TTEP/TTMGP
      WRITE(MFILE,*) ' defect calc.  :', 1.D2*TTDP/TTMGP
      WRITE(MFILE,*) ' prolongation  :', 1.D2*TTPP/TTMGP
      WRITE(MFILE,*) ' restriction   :', 1.D2*TTRP/TTMGP
      WRITE(MFILE,1)
      ENDIF
C
      IF (IER.GT.0) IER=0
      IF (IER.NE.0) GOTO 99999
C
C***********************************************************************
C
      ISETLV=2
      ILEV=NLMAX
      CALL  SETLEV (ISETLV)
C
C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================
      IF (BMGP) THEN
       BSTOP=.TRUE.
       RETURN
      ENDIF
C
C***********************************************************************
C *** Update of U = U~ - k B P
C***********************************************************************
C
      CALL ZTIME(TTT0)
      CALL XDFKG(KM1,KB1,KB2,KCOLB,KLDB,L(LD1),L(LD2),KU1,KU2,KP,NU,NP,
     *           KWORK(L(KLMBD(NLEV))),KNMBD(NLEV),TSTEPH)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C ************************************************************************
C *** Update of P(n+1)=P(n) + alphap P for second order
C***********************************************************************
C
      CALL ZTIME(TTT0)
C
      IF (IPROJ.EQ.1) THEN
       ALPHAP=PRDIF1
       ALPHAT=1D0
      ELSE
       ALPHAP=PRDIF1
       ALPHAT=0D0
       IF (IPROJ.EQ.-1) IPROJ=1
       IF (IPROJ.LT.-1) IPROJ=IPROJ+1
      ENDIF
C
      CALL LCP1(DWORK(KFP)     ,DWORK(KAUXP)   ,NP)
C
      CALL LCP1(DWORK(L(LPOLD)),DWORK(KFP)     ,NP)
      CALL LCP1(DWORK(KP)      ,DWORK(L(LDP))  ,NP)
      CALL LLC1(DWORK(L(LPOLD)),DWORK(KP)      ,NP,ALPHAT,ALPHAP)
      CALL LCP1(DWORK(L(LDP))  ,DWORK(L(LPOLD)),NP)
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      APDIV=0D0
      APDIV=PRDIF2*TMSTEP
      IF (APDIV.NE.0D0) 
     *    CALL XDFKDV(KP,KAUXP,VWORK(L(KLAREA(NLEV))),NP,APDIV)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C ************************************************************************
C
      CALL ZTIME(TTT0)
      CALL LCP1(DWORK(KFP),DWORK(L(LDP)),NP)
      CALL LLC1(DWORK(KP) ,DWORK(L(LDP)),NP,-1D0,1D0)
C
      CALL LL21(DWORK(L(LD1)+2*NU),NP,DSXN)
      RELP2=DSXN/(SQRT(DBLE(NP))*TSTEPH)
      CALL LLI1(DWORK(L(LDP)),NP,DSXN,INDMAX)
      RELPM=DSXN/TSTEPH
      THSTEP=TSTEPH
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      IF (MSHOW.GE.2) WRITE(MTERM,3)
      IF (MSHOW.GE.2) WRITE(MTERM,10002) ITNSL,ITNS,
     *                                   TIMENS,RELP2,RELPM
      IF (MSHOW.GE.2) WRITE(MTERM,3)
      IF (MSHOW.GE.2) WRITE(MTERM,*)
C
      IF (MSHOW.GE.1) WRITE(MFILE,3)
      IF (MSHOW.GE.1) WRITE(MFILE,10002) ITNSL,ITNS,
     *                                   TIMENS,RELP2,RELPM
      IF (MSHOW.GE.1) WRITE(MFILE,3)
      IF (MSHOW.GE.1) WRITE(MFILE,*)
C
100   CONTINUE
C
      GOTO 99999
C
C
C
   1  FORMAT(80('-'))
   3  FORMAT(80('+'))
1000  FORMAT (6E12.5)
1001  FORMAT('  IT DIV-L2',3X,'RHOMGP')
1003  FORMAT(I4,2(D9.2))
10002 FORMAT ('#',I4,'(',I4,')',1X,'TIME=',D10.3,1X,'REL2(P)=',
     *        D10.3,1X,'RELM(P)=',D10.3)
C
C
C
99999 END
