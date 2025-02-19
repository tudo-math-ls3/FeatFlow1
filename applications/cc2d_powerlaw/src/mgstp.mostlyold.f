************************************************************************
      SUBROUTINE MGSTP(MFILE,MSHOW,NITNSL,LRHS,BSTOP,BNLEND)
************************************************************************
*
*   Purpose: - performs 1 macro step of coupled method
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299, NNLEV=9,NNAB=21,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NBLOCF=2)
      DIMENSION BCONF(NBLOCF),ARRDF(NBLOCF),BSNGLF(NBLOCF)
      DIMENSION KFN(NBLOCF),KF(NNAB,NBLOCF),LF(NBLOCF)
      CHARACTER ARRDF*6
      DATA BCONF  /.FALSE.,.FALSE./,ARRDF/'DF1   ','DF2   '/
      DATA BSNGLF /.FALSE.,.FALSE./,KFN/1,1/
C
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
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
      COMMON /NSSAV/  INSAV,INSAVN
      COMMON /NSSAVF/ DTFILM,DTFILO,DTAVS,DTAVSO,DTGMV,DTGMVO,
     *                IFUSAV,IFPSAV,IFXSAV,IAVS,IGMV,IFINIT
      COMMON /NSADAT/ TIMEMX,DTMIN,DTMAX,DTFACT,TIMEIN,EPSADI,EPSADL,
     *                EPSADU,IEPSAD,IADIN,IREPIT,IADTIM
      COMMON /NSTIME/ TTGRID,TTPOST,TTADF,TTUPW,TTBDR,TTLC
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
      IF (ISTAT.EQ.0) TIMENS=0D0     
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C=======================================================================
C
      IF (ISTAT.NE.0) THEN
C
       CALL ZTIME(TTT0)
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
C***********************************************************************
C
       CALL LLC1(DWORK(KP),DWORK(KP),NP,0D0,TSTEPH)
C
C***********************************************************************
C
       CALL ZTIME(TTT0)
       IF (IBDR.GE.2) THEN
        DO 10 ILEV=NLMIN,NLMAX
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
     *              KWORK(L(KLEBD(ILEV))),KWORK(L(KLVERT(ILEV))),
     *              KWORK(L(KLMID(ILEV))),KWORK(L(KLNPR(ILEV))),
     *              DWORK(L(KLDBD(ILEV))),DWORK(L(KLMBDP(ILEV))),
     *              NMBD,INEUM)
        KNMBD(ILEV)=NMBD
C
10      CONTINUE
       ENDIF
C      
       CALL ZTIME(TTT1)
       TTBDR=TTBDR+TTT1-TTT0
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
     *                KWORK(L(KLVERT(NLEV))),KWORK(L(KLMID(NLEV))),
     *                KWORK(L(KLNPR(NLEV))),DWORK(L(LCORVG)),
     *                DWORK(L(KLMBDP(NLEV))),DWORK(KF1),DWORK(KF2),
     *                -THSTEP)
         TIMENS=TIMENS+TSTEPH
        ENDIF
        CALL ZTIME(TTT1)
        TTBDR=TTBDR+TTT1-TTT0
C
        CALL ZTIME(TTT0)
        LF(1)=0
        LF(2)=0
        IF (IELT.EQ.0) 
     *   CALL  XVB0 (LF,NU,NBLOCF,ICLRF,E031,RHS,BCONF,KF,KFN,ICUBF,
     *               ARRDF,BSNGLF)
        IF (IELT.EQ.1) 
     *   CALL  XVB0 (LF,NU,NBLOCF,ICLRF,E030,RHS,BCONF,KF,KFN,ICUBF,
     *               ARRDF,BSNGLF)
        IF (IELT.EQ.2) 
     *   CALL  XVBM0(LF,NU,NBLOCF,ICLRF,EM31,RHS,BCONF,KF,KFN,ICUBF,
     *               ARRDF,BSNGLF)
        IF (IELT.EQ.3) 
     *   CALL  XVBM0(LF,NU,NBLOCF,ICLRF,EM30,RHS,BCONF,KF,KFN,ICUBF,
     *               ARRDF,BSNGLF)
        IF (IER.NE.0) GOTO 99999
        LF1=LF(1)
        LF2=LF(2)
        CALL ZTIME(TTT1)
        TTLC=TTLC+TTT1-TTT0
C
        CALL ZTIME(TTT0)
        IF (IBDR.GE.1) THEN
         CALL PDSET (KWORK(L(KLVBD(NLEV))),KWORK(L(KLEBD(NLEV))),
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
     *               DWORK(KF2),KWORK(L(KLMBD(ILEV))),
     *               DWORK(L(KLDBD(ILEV))),KWORK(L(LNPR)),
     *               KNMBD(ILEV),NVT,PARX,PARY,UE)
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
        CALL PDSET (KWORK(L(KLVBD(NLEV))),KWORK(L(KLEBD(NLEV))),
     *              KWORK(L(KLVERT(NLEV))),KWORK(L(KLMID(NLEV))),
     *              KWORK(L(KLNPR(NLEV))),DWORK(L(LCORVG)),
     *              DWORK(L(KLMBDP(NLEV))),DWORK(KF1),DWORK(KF2),TSTEPB)
        IF ((IFRSTP.EQ.1).AND.(ITMOD.EQ.2)) TIMENS=TIMENS+TSTEPH
       ENDIF
       CALL ZTIME(TTT1)
       TTBDR=TTBDR+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      CALL XMADF1(KM1,KST1,KA1,KCOLA,KLDA,KF1,KF2,KU1,KU2,NA,NU,THSTEP)
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
ccc       DO 20 INU=1,NU
ccc       WRITE(6,*) INU,VWORK(KA1+KWORK(KLDA+INU-1)-1),
ccc     *            VWORK(KST1+KWORK(KLDA+INU-1)-1),VWORK(KM1+INU-1),
ccc     *            DWORK(KU1+INU-1),DWORK(KU2+INU-1),
ccc     *            DWORK(KF1+INU-1),DWORK(KF2+INU-1)
ccc 20    CONTINUE
C
C
       IF (THSTEP.NE.0D0) THEN
C
        CALL ZTIME(TTT0)
        THSTEP=-THSTEP
c
        write (*,*) 'simmada'
        CALL STABIL (DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
     *                KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *                DWORK(L(LCORVG)),E031,COEFFN,2,-1D0)
c
c$$$        IF (IUPW.EQ.1) THEN
c$$$         CALL GUPWD (DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
c$$$     *               1D0,0D0,DWORK(KU1),DWORK(KU2),
c$$$     *               DWORK(KF1),DWORK(KF2),VWORK(KA1),KWORK(KCOLA),
c$$$     *               KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
c$$$     *               DWORK(L(LCORVG)),2)
c$$$        ELSE
c$$$         IF (IELT.EQ.0) 
c$$$     *    CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
c$$$     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
c$$$     *                DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
c$$$     *                KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
c$$$     *                DWORK(L(LCORVG)),E031,COEFFN,2,-1D0)
c$$$         IF (IELT.EQ.1) 
c$$$     *    CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
c$$$     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
c$$$     *                DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
c$$$     *                KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
c$$$     *                DWORK(L(LCORVG)),E030,COEFFN,2,-1D0)
c$$$         IF (IELT.EQ.2) 
c$$$     *    CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
c$$$     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
c$$$     *                DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
c$$$     *                KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
c$$$     *                DWORK(L(LCORVG)),EM31,COEFFN,2,-1D0)
c$$$         IF (IELT.EQ.3) 
c$$$     *    CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
c$$$     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
c$$$     *                DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
c$$$     *                KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
c$$$     *                DWORK(L(LCORVG)),EM30,COEFFN,2,-1D0)
c$$$        ENDIF
        THSTEP=-THSTEP
        CALL ZTIME(TTT1)
        TTUPW=TTUPW+TTT1-TTT0
       ELSE
        IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
         CALL ZTIME(TTT0)
         TOSTEP=THSTEP
         THSTEP=1D0
         IF (IELT.EQ.0) 
     *    CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
     *                KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),E031,COEFFN,2,0D0)
         IF (IELT.EQ.1) 
     *    CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
     *                KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *                DWORK(L(LCORVG)),E030,COEFFN,2,0D0)
         IF (IELT.EQ.2) 
     *    CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
     *                KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *                DWORK(L(LCORVG)),EM31,COEFFN,2,0D0)
         IF (IELT.EQ.3) 
     *    CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(KF1),DWORK(KF2),VWORK(KA1),NA,KWORK(KCOLA),
     *                KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *                DWORK(L(LCORVG)),EM30,COEFFN,2,0D0)
        ENDIF
        THSTEP=TOSTEP
        CALL ZTIME(TTT1)
        TTUPW=TTUPW+TTT1-TTT0
       ENDIF
C
       CALL ZTIME(TTT0)
       CALL BDRSET (DWORK(KU1),DWORK(KU2),DWORK(KF1),
     *              DWORK(KF2),KWORK(L(KLMBD(ILEV))),
     *              DWORK(L(KLDBD(ILEV))),KWORK(L(LNPR)),
     *              KNMBD(ILEV),NVT,PARX,PARY,UE)
       CALL ZTIME(TTT1)
       TTBDR=TTBDR+TTT1-TTT0
C
       THSTEP=TMSTEP
C
      ENDIF
C
      stop
************************************************************************
C *** fixed point defect correction for stationary NS equations
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
      CALL ZTIME(TTT0)
      ILEV=NLEV
      ISETLV=2
      CALL  SETLEV (ISETLV)
C
C***********************************************************************
C
      CALL LLC1(DWORK(KP),DWORK(KP),NP,0D0,1D0/TSTEPH)
C
C***********************************************************************
C
      CALL LL21(DWORK(KU1),2*NU,DSXN)
      RELU2=DSXN/SQRT(DBLE(NU))
      CALL LL21(DWORK(KU1+2*NU),NP,DSXN)
      RELP2=DSXN/SQRT(DBLE(NP))
      THSTEP=TSTEPH
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      IF (MSHOW.GE.2) WRITE(MTERM,3)
      IF (MSHOW.GE.2) WRITE(MTERM,10002) ITNSL,ITNS,TIMENS,RELU2,RELP2
      IF (MSHOW.GE.2) WRITE(MTERM,3)
      IF (MSHOW.GE.2) WRITE(MTERM,*)
C
      IF (MSHOW.GE.1) WRITE(MFILE,3)
      IF (MSHOW.GE.1) WRITE(MFILE,10002) ITNSL,ITNS,TIMENS,RELU2,RELP2
      IF (MSHOW.GE.1) WRITE(MFILE,3)
      IF (MSHOW.GE.1) WRITE(MFILE,*)
C
C=======================================================================
C *** Return if stationary calculation !!!
C=======================================================================
      IF (ISTAT.EQ.0) RETURN
C
100   CONTINUE
C
C
C
   1  FORMAT(80('-'))
   3  FORMAT(80('+'))
1000  FORMAT (6E12.5)
1001  FORMAT(' IT DIV-L2',3X,'RHOMGP')
1003  FORMAT(I4,2(D9.2))
10002 FORMAT ('#',I4,'(',I4,')',1X,'TIME=',D10.3,2X,'NORM(U)=',
     *        D14.7,2X,'NORM(P)=',D14.7)
C
C
C
99999 END
