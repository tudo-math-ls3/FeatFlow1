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
      PARAMETER (NNARR=299, NNLEV=9, NNAB=21,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NBLOCF=3)
      DIMENSION BCONF(NBLOCF),ARRDF(NBLOCF),BSNGLF(NBLOCF)
      DIMENSION KFN(NBLOCF),KF(NNAB,NBLOCF),LF(NBLOCF)
      CHARACTER ARRDF*6
      DATA BCONF  /.FALSE.,.FALSE.,.FALSE./,ARRDF/'DF1   ','DF2   ',
     *            'DF3   '/
      DATA BSNGLF /.FALSE.,.FALSE.,.FALSE./,KFN/1,1,1/
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
C
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
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
      COMMON /MGIEL/  KLINT(NNLEV)
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C     E X T E R N A L S
C-----------------------------------------------------------------------
C *** Parametrization of the domain
      EXTERNAL PARX,PARY,PARZ
C *** Coefficient of stiffness matrix, right hand side, exact solution
      EXTERNAL COEFFN,RHS,UE,PE,UEX,UEY,UEZ
C *** definition of finite elements
      EXTERNAL E030,E031,EM31,EM30
C
C=======================================================================
C
      ILINT=0
C
      KF(1,1)=1
      KF(1,2)=1
      KF(1,3)=1
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
       NVT=KNVT(ILEV)
       NAT=KNAT(ILEV)
       NEL=KNEL(ILEV)
       NET=0
C
       CALL XLCP3(KLNPRO(ILEV),KLNPR(ILEV),NVT+NAT)
       IF (IER.NE.0) GOTO 99999
C
       NABD =0
       INEUM=0
       KNABD(ILEV)=0
       CALL ZDISP(0,KLABD(ILEV),'KABD  ')
       CALL ZDISP(0,KELBD(ILEV),'KELBD ')
       IF (IER.NE.0) GOTO 99999
       CALL ZNEW(NAT,3,KLABD(ILEV),'KABD  ')
       CALL ZNEW(NAT,3,KELBD(ILEV),'KELBD ')
       IF (IER.NE.0) GOTO 99999
       CALL BDRNEU(KWORK(L(KLABD(ILEV))),KWORK(L(KLNPR(ILEV))),
     *             DWORK(L(KLCVG(ILEV))),INEUM,KWORK(L(KLAREA(ILEV))),
     *             KWORK(L(KLADJ(ILEV))),KWORK(L(KLVERT(ILEV))),
     *             KWORK(L(KELBD(ILEV))),ilev,nlmax)
       CALL ZDISP(NABD,KLABD(ILEV),'KABD  ')
       CALL ZDISP(NABD,KELBD(ILEV),'KELBD ')
       IF (IER.NE.0) GOTO 99999
       KNABD(ILEV)=NABD
C
10     CONTINUE
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
         CALL PDSET(DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *              KWORK(L(KLABD(NLEV))),KNABD(NLEV),DWORK(L(LCORVG)),
     *              KWORK(L(KELBD(NLEV))),KWORK(L(LVERT)),
     *              KWORK(L(LAREA)),-TSTEPH)
         TIMENS=TIMENS+TSTEPH
        ENDIF
        CALL ZTIME(TTT1)
        TTBDR=TTBDR+TTT1-TTT0
C
        CALL ZTIME(TTT0)
        LF(1)=0
        LF(2)=0
        LF(3)=0
       IF (IELT.EQ.0)
     *  CALL  XVB0 (LF,NU,NBLOCF,ICLRF,E031,
     *              RHS,BCONF,KF,KFN,ICUBF,ILINT,BSNGLF,ARRDF)
       IF (IELT.EQ.1) 
     *  CALL  XVB0 (LF,NU,NBLOCF,ICLRF,E030,
     *              RHS,BCONF,KF,KFN,ICUBF,ILINT,BSNGLF,ARRDF)
       IF (IELT.EQ.2) 
     *  CALL  XVBM0(LF,NU,NBLOCF,ICLRF,EM31,
     *              RHS,BCONF,KF,KFN,ICUBF,ILINT,BSNGLF,ARRDF)
       IF (IELT.EQ.3) 
     *  CALL  XVBM0(LF,NU,NBLOCF,ICLRF,EM30,
     *              RHS,BCONF,KF,KFN,ICUBF,ILINT,BSNGLF,ARRDF)
        IF (IER.NE.0) GOTO 99999
        LF1=LF(1)
        LF2=LF(2)
        LF3=LF(3)
        CALL ZTIME(TTT1)
        TTLC=TTLC+TTT1-TTT0
C
        CALL ZTIME(TTT0)
        IF (IBDR.GE.1) THEN
         CALL PDSET(DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *              KWORK(L(KLABD(NLEV))),KNABD(NLEV),DWORK(L(LCORVG)),
     *              KWORK(L(KELBD(NLEV))),KWORK(L(LVERT)),
     *              KWORK(L(LAREA)),TRSTEP)
        ENDIF
        CALL ZTIME(TTT1)
        TTBDR=TTBDR+TTT1-TTT0
C
        CALL ZTIME(TTT0)
        CALL  LLC1 (DWORK(L(LF1)),DWORK(KF1),NU,TRSTEP,1D0)
        CALL  LLC1 (DWORK(L(LF2)),DWORK(KF2),NU,TRSTEP,1D0)
        CALL  LLC1 (DWORK(L(LF3)),DWORK(KF3),NU,TRSTEP,1D0)
        CALL ZTIME(TTT1)
        TTLC=TTLC+TTT1-TTT0
C
        CALL ZTIME(TTT0)
        CALL BDRSET(DWORK(KU1),DWORK(KU2),DWORK(KU3),DWORK(KF1),
     *              DWORK(KF2),DWORK(KF3),KWORK(L(KLABD(ILEV))),
     *              KNABD(ILEV),DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *              KWORK(L(LAREA)),KWORK(L(KELBD(ILEV))),UE)
        CALL ZTIME(TTT1)
        TTBDR=TTBDR+TTT1-TTT0
C
        CALL ZTIME(TTT0)
        CALL  LCP1 (DWORK(L(LF1)), DWORK(L(LRHS)),   NU)
        CALL  LCP1 (DWORK(L(LF2)), DWORK(L(LRHS)+NU),NU)
        CALL  LCP1 (DWORK(L(LF3)), DWORK(L(LRHS)+NU+NU),NU)
C
        CALL  ZDISP (0,LF1,ARRDF(1))
        CALL  ZDISP (0,LF2,ARRDF(2))
        CALL  ZDISP (0,LF3,ARRDF(3))
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
        CALL PDSET(DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *             KWORK(L(KLABD(NLEV))),KNABD(NLEV),DWORK(L(LCORVG)),
     *             KWORK(L(KELBD(NLEV))),KWORK(L(LVERT)),
     *             KWORK(L(LAREA)),TSTEPB)
        IF ((IFRSTP.EQ.1).AND.(ITMOD.EQ.2)) TIMENS=TIMENS+TSTEPH
       ENDIF
       CALL ZTIME(TTT1)
       TTBDR=TTBDR+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      CALL XMADF1(KM1,KST1,KA1,KCOLA,KLDA,KF1,KF2,KF3,
     *            KU1,KU2,KU3,NA,NU,THSTEP)
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
        IF (IUPW.EQ.1) THEN
        CALL  GUPWD(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *              DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *              DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *              DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *              VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *              KWORK(L(LVERT)),KWORK(L(LAREA)),DWORK(L(LCORVG)),2)
        ELSE
         IF (IELT.EQ.0) 
     *   CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *               DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *               DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *               DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *               VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *               KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *               KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),E031,COEFFN,
     *               2,-1D0)
         IF (IELT.EQ.1) 
     *   CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *               DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *               DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *               DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *               VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *               KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *               KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),E030,COEFFN,
     *               2,-1D0)
         IF (IELT.EQ.2) 
     *   CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *               DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *               DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *               DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *               VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *               KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *               KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),EM31,COEFFN,
     *               2,-1D0)
         IF (IELT.EQ.3) 
     *   CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *               DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *               DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *               DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *               VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *               KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *               KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),EM30,COEFFN,
     *               2,-1D0)
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
     *   CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *               DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *               DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *               DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *               VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *               KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *               KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),E031,COEFFN,
     *               2,0D0)
        IF (IELT.EQ.1) 
     *   CALL SUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *               DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *               DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *               DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *               VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *               KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *               KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),E030,COEFFN,
     *               2,0D0)
        IF (IELT.EQ.2) 
     *   CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *               DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *               DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *               DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *               VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *               KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *               KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),EM31,COEFFN,
     *               2,0D0)
        IF (IELT.EQ.3) 
     *   CALL SUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *               DWORK(KU1),DWORK(KU2),DWORK(KU3),1D0,0D0,
     *               DWORK(KU1),DWORK(KU2),DWORK(KU3),
     *               DWORK(KF1),DWORK(KF2),DWORK(KF3),
     *               VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *               KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *               KWORK(L(KLINT(NLEV))),DWORK(L(LCORVG)),EM30,COEFFN,
     *               2,0D0)
        ENDIF
        THSTEP=TOSTEP
        CALL ZTIME(TTT1)
        TTUPW=TTUPW+TTT1-TTT0
       ENDIF
C
       CALL ZTIME(TTT0)
       CALL BDRSET (DWORK(KU1),DWORK(KU2),DWORK(KU3),DWORK(KF1),
     *              DWORK(KF2),DWORK(KF3),KWORK(L(KLABD(ILEV))),
     *              KNABD(ILEV),DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *              KWORK(L(LAREA)),KWORK(L(KELBD(ILEV))),UE)
       CALL ZTIME(TTT1)
       TTBDR=TTBDR+TTT1-TTT0
C
       THSTEP=TMSTEP
C
      ENDIF
C
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
      CALL LL21(DWORK(KU1),3*NU,DSXN)
      RELU2=DSXN/SQRT(DBLE(NU))
      CALL LL21(DWORK(KU1+3*NU),NP,DSXN)
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
