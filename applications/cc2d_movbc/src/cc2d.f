************************************************************************
      PROGRAM  CC2D_MOVBC  
************************************************************************
*
*   Explanation: works similar as standard CC2D 
*
*   Differences: bndry.f contains position (and movement) of "fictitious
*                boundary components" 
*
*                indat2d.f contains position (and movement) of 
*                "fictitious boundary components" and additionally the 
*                Dirichlet values on the surface
*
*                INPR=-2 is set for such "fictitious b.c.'s"
*                INPR=11 is set for the corresponding nodes for 
*                        visualization reasons only
*
*   Perform the given example (stat. flow around circle) 
*   and modify w.r.t. your application !!!
*
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299,NNLEV=9,NNAB=21)
      INCLUDE 'cc2d.inc'
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CDATA*60
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
      COMMON /NSADAT/ TIMEMX,DTMIN,DTMAX,DTFACT,TIMEIN,EPSADI,EPSADL,
     *                EPSADU,IEPSAD,IADIN,IREPIT,IADTIM
      COMMON /NSSAV/  INSAV,INSAVN
      COMMON /NSSAVF/ DTFILM,DTFILO,DTAVS,DTAVSO,DTGMV,DTGMVO,
     *                IFUSAV,IFPSAV,IFXSAV,IAVS,IGMV,IFINIT
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
C *** Coefficient of exact solution
      EXTERNAL UE,PE,UEX,UEY
C *** definition of finite elements
      EXTERNAL E030,E031,EM30,EM31,E010
      EXTERNAL ZVALUE1
C=======================================================================
C     Initialization
C=======================================================================
      CALL ZTIME(TTTSUM)
C
      CALL ZTIME(TTT0)
      CALL ZINIT(NNWORK,'feat.msg','#data/CC2D.ERR','#data/CC2D.PRT',
     *                             '#data/CC2D.SYS','#data/CC2D.TRC') 
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      MDATA=79
      CDATA='#data/cc2d.dat'
      CALL  INIT1 (MDATA,CDATA,MFILE,MSHOW,IWORKG,IWMAXG)
      IWORKI=IWORK
      IWMAXI=IWMAX
C
      IFILEN=0
      ITFILM=IFINIT-1
      NSUBST=0
      NNONL =0
      NMG   =0
C
      TIMEST=TIMENS
      CALL ZTIME(TTT0)
C
      IF (ISTAT.NE.0) THEN
       ISETLV=2
       ILEV=NLMAX
       CALL  SETLEV (ISETLV)
C
       LRHS=0
       IF (IRHS.GE.1) THEN
        CALL ZNEW(NUP,1,LRHS  ,'DRHS  ')
        IF (IER.NE.0) GOTO 99998
        CALL LCP1(DWORK(KF1),DWORK(L(LRHS)),NUP)
       ENDIF
C
       CALL ZNEW(NUP,-1,LTM0  ,'DMT0  ')
       IF (IER.NE.0) GOTO 99998
C
       IF (ITEXL.NE.0) THEN
        CALL ZNEW(NUP,-1,LTML ,'DMTL  ')
        IF (IER.NE.0) GOTO 99998
        IF (ABS(IADTIM).EQ.3) THEN
         CALL ZNEW(NUP,-1,LTMLO,'DMTLO ')
         IF (IER.NE.0) GOTO 99998
         CALL LCP1(DWORK(KU1),DWORK(L(LTML)),NUP)
        ENDIF
       ENDIF
C
       IF (IADTIM.NE.0) THEN
        CALL ZNEW(NUP,-2,LTM3  ,'DMT3  ')
        IF (IER.NE.0) GOTO 99998
       ENDIF
      ENDIF
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C=======================================================================
C *** Begin of nonsteady loop
C=======================================================================
C
C***********************************************************************
C *** NNITR = max number of repetitions for |IADAPT|=3
      NNITR=IREPIT
C***********************************************************************
C
      DO 100 ITNS=1,NITNS
      ITNSR=0
C
      IF (ISTAT.NE.0) THEN
       CALL ZTIME(TTT0)
       CALL LCP1(DWORK(KU1),DWORK(L(LTM0)),NUP)
       IF ((ITEXL.NE.0).AND.(ABS(IADTIM).EQ.3)) THEN
        CALL LCP1(DWORK(L(LTML)),DWORK(L(LTMLO)),NUP)
        TIMO31=TIML31
        TIMO32=TIML32
        TIMO11=TIML11
        TIMO12=TIML12
       ENDIF
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
      ENDIF
C
110   TSTEPO=    TSTEP
      TSTEP =3D0*TSTEPO
      TIMNSH=TIMENS
C
      IF (MSHOW.GE.2) WRITE(MTERM,*)
      IF (MSHOW.GE.1) WRITE(MFILE,*)
C
      IF (IADTIM.EQ.0) GOTO 120
C
C=======================================================================
C *** first time step: 1 step
C=======================================================================
C
      NITNSL=1
C
      CALL ZTIME(TTT0)
      IF (ITEXL.NE.0) THEN
       IF (ITNS.EQ.1) THEN
        ITEXL =-1
        CALL LCP1(DWORK(KU1),DWORK(L(LTML)),NUP)
       ELSE
        ITEXL = 1
        TIML11=TIML32
       ENDIF
      ENDIF
C
      IFRSTH=0
      IF ((NITNSL.EQ.1).AND.(IFRSTP.EQ.1)) THEN
       THETAH=THETA
       THETA =0.5D0
       IFRSTH=1
       IFRSTP=0
      ENDIF
C
      TSTEP=TSTEP/DBLE(NITNSL)
C
      IF (MSHOW.GE.2) WRITE(MTERM,2)
      IF (MSHOW.GE.2) WRITE(MTERM,2002) ITNS,TIMENS,TSTEP
      IF (MSHOW.GE.2) WRITE(MTERM,2)
      IF (MSHOW.GE.1) WRITE(MFILE,2)
      IF (MSHOW.GE.1) WRITE(MFILE,2002) ITNS,TIMENS,TSTEP
      IF (MSHOW.GE.1) WRITE(MFILE,2)
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      BSTOP=.FALSE.
      BNL1 =.FALSE.
      CALL MGSTP(MFILE,MSHOW,NITNSL,LRHS,BSTOP,BNL1)
      TSTEP=TSTEP*DBLE(NITNSL)
      NSUBST=NSUBST+NITNSL
C
      IF ((NITNSL.EQ.1).AND.(IFRSTH.EQ.1)) THEN
       THETA =THETAH
       IFRSTP=1
      ENDIF
C
      CALL ZTIME(TTT0)
      ISETLV=2
      ILEV=NLEV
      CALL  SETLEV (ISETLV)
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================
      IF ((BSTOP).AND.(ITNSR.LE.NNITR)) THEN
       CALL ZTIME(TTT0)
       ITNSR=ITNSR+1
       CALL LCP1(DWORK(L(LTM0)),DWORK(KU1),NUP)
       IF ((ITEXL.NE.0).AND.(ABS(IADTIM).EQ.3)) THEN
        CALL LCP1(DWORK(L(LTMLO)),DWORK(L(LTML)),NUP)
        TIML31=TIMO31
        TIML32=TIMO32
        TIML11=TIMO11
        TIML12=TIMO12
       ENDIF
       TSTEPN=TSTEPO/SQRT(DTFACT)
       TSTEP=TSTEPN
       TIMENS=TIMNSH
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
       GOTO 110
      ENDIF
C
C
      CALL ZTIME(TTT0)
      DO 112 IUP=1,NUP
112   VWORK(L(LTM3)+IUP-1)=REAL(DWORK(KU1+IUP-1))
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      TIMENS=TIMNSH
C
C=======================================================================
C *** second time step: 3 steps
C=======================================================================
C
120   IF (ISTAT.NE.0) THEN
       CALL ZTIME(TTT0)
       CALL LCP1(DWORK(L(LTM0)),DWORK(KU1),NUP)
       NITNSL=3
C
       IF (ITEXL.NE.0) THEN
        IF (ITNS.EQ.1) THEN
         ITEXL=-3
         CALL LCP1(DWORK(KU1),DWORK(L(LTML)),NUP)
        ELSE
         ITEXL= 3
        ENDIF
       ENDIF
C
       TSTEP=TSTEP/DBLE(NITNSL)
       IF (MSHOW.GE.2) WRITE(MTERM,2)
       IF (MSHOW.GE.2) WRITE(MTERM,2001) ITNS,TIMENS,TSTEP
       IF (MSHOW.GE.2) WRITE(MTERM,2)
       IF (MSHOW.GE.1) WRITE(MFILE,2)
       IF (MSHOW.GE.1) WRITE(MFILE,2001) ITNS,TIMENS,TSTEP
       IF (MSHOW.GE.1) WRITE(MFILE,2)
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
      ELSE
       NITNSL=1
       TSTEP=1D0
      ENDIF
C
      BSTOP=.FALSE.
      BNL2 =.FALSE.
      CALL MGSTP(MFILE,MSHOW,NITNSL,LRHS,BSTOP,BNL2)
      IF (IER.NE.0) GOTO 99998
      TSTEP=TSTEP*DBLE(NITNSL)
      NSUBST=NSUBST+NITNSL
C
      IF ((ISTAT.EQ.0).OR.(IADTIM.EQ.0)) THEN
       TSTEP=TSTEPO
       GOTO 900
      ENDIF
C
      CALL ZTIME(TTT0)
      ISETLV=2
      ILEV=NLEV
      CALL  SETLEV (ISETLV)
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C=======================================================================
C *** Unexpected STOP !!!
C=======================================================================
      IF (((BSTOP).OR.(.NOT.BNL2)).AND.(ITNSR.LE.NNITR)) THEN
       IF (((.NOT.BSTOP).AND.(.NOT.BNL2)).AND.(ABS(IADTIM).EQ.1)) 
     *        GOTO 130
       CALL ZTIME(TTT0)
       ITNSR=ITNSR+1
       CALL LCP1(DWORK(L(LTM0)),DWORK(KU1),NUP)
       IF ((ITEXL.NE.0).AND.(ABS(IADTIM).EQ.3)) THEN
        CALL LCP1(DWORK(L(LTMLO)),DWORK(L(LTML)),NUP)
        TIML31=TIMO31
        TIML32=TIMO32
        TIML11=TIMO11
        TIML12=TIMO12
       ENDIF
       IF (.NOT.BNL2) TSTEPN=TSTEPO/SQRT(DTFACT)
       IF (BSTOP)     TSTEPN=TSTEPO/DTFACT
       TSTEP=TSTEPN
       TIMENS=TIMNSH
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
       GOTO 110
      ENDIF
C
C=======================================================================
C     Time step control
C=======================================================================
C
130   CALL ZTIME(TTT0)
      ISETLV=2
      ILEV=NLEV
      CALL  SETLEV (ISETLV)
C
      DO 132 IUP=1,NUP
132   DWORK(L(LD1)+IUP-1)=DWORK(KU1+IUP-1)-DBLE(VWORK(L(LTM3)+IUP-1))
C
C
      RELU21=0D0
      RELUM1=0D0
      RELP21=0D0
      RELPM1=0D0
C
      IF ((IEPSAD.GE.1).OR.(IEPSAD.EQ.-1)) THEN
       CALL LL21(DWORK(KU1)   ,2*NU,DSXN)
       RELU20=DSXN/(SQRT(DBLE(2*NU)))
       IF (RELU20.LE.1D0) RELU20=1D0
       CALL LL21(DWORK(L(LD1)),2*NU,DSXN)
       RELU21=DSXN/(SQRT(DBLE(2*NU))*RELU20)
      ENDIF
C
      IF ((IEPSAD.GE.1).OR.(IEPSAD.EQ.-2)) THEN
       CALL LLI1(DWORK(KU1   ),2*NU,DSXN,INDMAX)
       RELUM0=DSXN
       IF (RELUM0.LE.1D0) RELUM0=1D0
       CALL LLI1(DWORK(L(LD1)),2*NU,DSXN,INDMAX)
       RELUM1=DSXN/RELUM0
      ENDIF
C
      IF ((IEPSAD.GE.1).OR.(IEPSAD.EQ.-3)) THEN
       CALL LL21(DWORK(KU1   +2*NU),NP,DSXN)
       RELP20=DSXN/(SQRT(DBLE(NP)))
       IF (RELP20.LE.1D0) RELP20=1D0
       CALL LL21(DWORK(L(LD1)+2*NU),NP,DSXN)
       RELP21=DSXN/(SQRT(DBLE(NP))*RELP20)
      ENDIF
C
      IF ((IEPSAD.GE.1).OR.(IEPSAD.EQ.-4)) THEN
       CALL LLI1(DWORK(KU1   +2*NU),NP,DSXN,INDMAX)
       RELPM0=DSXN
       IF (RELPM0.LE.1D0) RELPM0=1D0
       CALL LLI1(DWORK(L(LD1)+2*NU),NP,DSXN,INDMAX)
       RELPM1=DSXN/RELPM0
      ENDIF
C
C
      IF (ABS(IEPSAD).EQ.1)  RELT1=RELU21
      IF (ABS(IEPSAD).EQ.2)  RELT1=RELUM1
      IF (ABS(IEPSAD).EQ.3)  RELT1=RELP21
      IF (ABS(IEPSAD).EQ.4)  RELT1=RELPM1
      IF (     IEPSAD.EQ.5)  RELT1=MAX(RELU21,RELP21)
      IF (     IEPSAD.EQ.6)  RELT1=MAX(RELUM1,RELPM1)
      IF (     IEPSAD.EQ.7)  RELT1=MAX(RELU21,RELUM1,RELP21,RELPM1)
      IF (     IEPSAD.EQ.8)  RELT1=MIN(RELU21,RELUM1,RELP21,RELPM1)
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C=======================================================================
C
      IF ((ABS(IADTIM).GE.1)) THEN
C
       CALL ZTIME(TTT0)
       CALL CRITAD(TIMEIN,TIMENS,TIMEST,EPSADI,EPSADL,EPSAD,IADIN)
C
      IF ((IFRSTP.NE.1).AND.(THETA.EQ.1D0)) THEN
        TSTEPN=TSTEPO*2D0*EPSAD/RELT1
        CTPROJ=2D0
       ELSE
        TSTEPN=TSTEPO*SQRT(8D0*EPSAD/RELT1)
        CTPROJ=8D0
       ENDIF
C
       IF ((ABS(IADTIM).GT.1)) THEN
        IF ((TSTEPN.GT.TSTEPO).AND.(.NOT.BNL1)) TSTEPN=TSTEPO
        IF ((TSTEPN.LT.TSTEPO/SQRT(DTFACT)).AND.(.NOT.BNL1)) 
     *       TSTEPN=TSTEPO/SQRT(DTFACT)
       ENDIF
C
       IF  (TSTEPN.LT.TSTEPO/DTFACT) TSTEPN=TSTEPO/DTFACT
C
       IF (ITNSR.GT.0) THEN
        HTFACT=DTFACT**(1D0/DBLE(ITNSR+1))
        IF (TSTEPN.GT.TSTEPO*HTFACT) TSTEPN=TSTEPO*HTFACT
       ELSE       
        IF (TSTEPN.GT.TSTEPO*DTFACT) TSTEPN=TSTEPO*DTFACT       
       ENDIF       
       IF (TSTEPN.LT.DTMIN) TSTEPN=DTMIN       
       IF (TSTEPN.GT.DTMAX) TSTEPN=DTMAX       
C
       ITYP=0
       DRELTN=TSTEPN/TSTEPO
       TSTEP =TSTEPN
       IF ((ABS(IADTIM).EQ.3).AND.(DRELTN.LT.EPSADU)) ITYP=1
C
      IF (MSHOW.GE.2) WRITE(MTERM,3)
      IF (MSHOW.GE.1) WRITE(MFILE,3)
      IF (MSHOW.GE.2) WRITE(MTERM,10002) TSTEPO,
     *    RELU21/CTPROJ,RELUM1/CTPROJ,RELP21/CTPROJ,RELPM1/CTPROJ
      IF (MSHOW.GE.1) WRITE(MFILE,10002) TSTEPO,
     *    RELU21/CTPROJ,RELUM1/CTPROJ,RELP21/CTPROJ,RELPM1/CTPROJ
C
      IF (MSHOW.GE.2) WRITE(MTERM,3001) ITYP,ITNS,TSTEPN,TSTEPO
      IF (MSHOW.GE.1) WRITE(MFILE,3001) ITYP,ITNS,TSTEPN,TSTEPO
      IF (MSHOW.GE.2) WRITE(MTERM,3)
      IF (MSHOW.GE.2) WRITE(MTERM,*)
      IF (MSHOW.GE.1) WRITE(MFILE,3)
      IF (MSHOW.GE.1) WRITE(MFILE,*)
C
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
C
C
       IF ((ITYP.EQ.1).AND.(ITNSR.LE.NNITR)) THEN
        CALL ZTIME(TTT0)
        ITNSR=ITNSR+1
        CALL LCP1(DWORK(L(LTM0)),DWORK(KU1),NUP)
        IF ((ITEXL.NE.0).AND.(ABS(IADTIM).EQ.3)) THEN
         CALL LCP1(DWORK(L(LTMLO)),DWORK(L(LTML)),NUP)
         TIML31=TIMO31
         TIML32=TIMO32
         TIML11=TIMO11
         TIML12=TIMO12
        ENDIF
        TIMENS=TIMNSH
        CALL ZTIME(TTT1)
        TTLC=TTLC+TTT1-TTT0
        GOTO 110
       ENDIF
C
C=======================================================================
C     Extrapolation step
C=======================================================================
C      
       CALL ZTIME(TTT0)
       IF ((IFRSTP.NE.1).AND.(THETA.EQ.1D0)) THEN
        EX1= 4D0/3D0
        EX2=-1D0/3D0
       ELSE
        EX1= 9D0/8D0
        EX2=-1D0/8D0
       ENDIF
C
       IF (((IADTIM.LT.-1).AND.(BNL1)).OR.(IADTIM.EQ.-1)) THEN
        DO 142 IUP=1,NUP
142     DWORK(KU1+IUP-1)= EX1*DWORK(KU1+IUP-1)
     *                   +EX2*DBLE(VWORK(L(LTM3)+IUP-1))
       ENDIF
C
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
C
      ENDIF
C
C=======================================================================
C
900   CALL ZTIME(TTT0)
      ISETLV=2
      ILEV=NLEV
      CALL  SETLEV (ISETLV)
C
C=======================================================================
C     Error handling
C=======================================================================
C
      IF (IERANA.GE.1)  THEN
       ICUBER=IERANA
       IF (IELT.EQ.0) 
     *  CALL ELPQU(DWORK(KU1),DWORK(KU2),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             DWORK(L(LCORVG)),E031,ICUBER,UE,UEX,UEY,MFILE)
       IF (IELT.EQ.1) 
     *  CALL ELPQU(DWORK(KU1),DWORK(KU2),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             DWORK(L(LCORVG)),E030,ICUBER,UE,UEX,UEY,MFILE)
       IF (IELT.EQ.2) 
     *  CALL ELPQN(DWORK(KU1),DWORK(KU2),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             DWORK(L(LCORVG)),EM31,ICUBER,UE,UEX,UEY,MFILE)
       IF (IELT.EQ.3) 
     *  CALL ELPQN(DWORK(KU1),DWORK(KU2),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             DWORK(L(LCORVG)),EM30,ICUBER,UE,UEX,UEY,MFILE)
C
       CALL ELPQP(DWORK(KP),DWORK(L(LDP)),KWORK(L(LVERT)),
     *            KWORK(L(LMID)),DWORK(L(LCORVG)),E010,ICUBER,PE,MFILE)
C
       IF (IERANA.NE.1) THEN
        ICUBER=1
        IF (IELT.EQ.0) 
     *   CALL ELPQU(DWORK(KU1),DWORK(KU2),KWORK(L(LVERT)),
     *              KWORK(L(LMID)),DWORK(L(LCORVG)),E031,ICUBER,
     *              UE,UEX,UEY,MFILE)
        IF (IELT.EQ.1) 
     *   CALL ELPQU(DWORK(KU1),DWORK(KU2),KWORK(L(LVERT)),
     *              KWORK(L(LMID)),DWORK(L(LCORVG)),E030,ICUBER,
     *              UE,UEX,UEY,MFILE)
        IF (IELT.EQ.2) 
     *   CALL ELPQN(DWORK(KU1),DWORK(KU2),KWORK(L(LVERT)),
     *              KWORK(L(LMID)),DWORK(L(LCORVG)),EM31,ICUBER,
     *              UE,UEX,UEY,MFILE)
        IF (IELT.EQ.3) 
     *   CALL ELPQN(DWORK(KU1),DWORK(KU2),KWORK(L(LVERT)),
     *              KWORK(L(LMID)),DWORK(L(LCORVG)),EM30,ICUBER,
     *              UE,UEX,UEY,MFILE)
C
        CALL ELPQP(DWORK(KP),DWORK(L(LDP)),KWORK(L(LVERT)),
     *             KWORK(L(LMID)),DWORK(L(LCORVG)),E010,ICUBER,PE,MFILE)
       ENDIF
C
      ENDIF
C
      IFPOST=1
      CALL FPOST(IFPOST,IFILEN,ITFILM,UE,MFILE,MSHOW)
C
      CALL ZTIME(TTT1)
      TTPOST=TTPOST+TTT1-TTT0
C
C=======================================================================
C
      IF (ISTAT.EQ.0) GOTO 999
C
C=======================================================================
C
      CALL ZTIME(TTT0)
      CALL LCP1(DWORK(L(LTM0)),DWORK(L(LD1)),NUP)
      CALL LLC1(DWORK(KU1)    ,DWORK(L(LD1)),NUP,-1D0,1D0)
C
      RELU2=0D0
      RELP2=0D0
      RELUM=0D0
      RELPM=0D0
C
      CALL LL21(DWORK(L(LD1)),2*NU,DSXN)
      RELU2=DSXN/(SQRT(DBLE(2*NU))*3D0*TSTEPO)
      CALL LL21(DWORK(L(LD1)+2*NU),NP,DSXN)
      RELP2=DSXN/(SQRT(DBLE(NP))*3D0*TSTEPO)
C
      IF ((ABS(IEPSAD).EQ.2).OR.(IEPSAD.GE.5)) THEN
       CALL LLI1(DWORK(L(LD1)),2*NU,DSXN,INDMAX)
       RELUM=DSXN/(3D0*TSTEPO)
      ENDIF
C
      IF ((ABS(IEPSAD).EQ.4).OR.(IEPSAD.GE.5)) THEN
       CALL LLI1(DWORK(L(LD1)+2*NU),NP,DSXN,INDMAX)
       RELPM=DSXN/(3D0*TSTEPO)
      ENDIF
C
      IF (ABS(IEPSAD).EQ.1)  RELT=RELU2
      IF (ABS(IEPSAD).EQ.2)  RELT=RELUM
      IF (ABS(IEPSAD).EQ.3)  RELT=RELP2
      IF (ABS(IEPSAD).EQ.4)  RELT=RELPM
      IF (    IEPSAD .EQ.5)  RELT=MAX(RELU2,RELP2)
      IF (    IEPSAD .EQ.6)  RELT=MAX(RELUM,RELPM)
      IF (    IEPSAD .EQ.7)  RELT=MAX(RELU2,RELUM,RELP2,RELPM)
      IF (    IEPSAD .EQ.8)  RELT=MIN(RELU2,RELUM,RELP2,RELPM)
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      IF (MSHOW.GE.2) WRITE(MTERM,2)
      IF (MSHOW.GE.2) WRITE(MTERM,20001) ITNS,NSUBST,TIMENS,RELU2,RELP2,
     *                                   RELT
      IF (MSHOW.GE.2) WRITE(MTERM,2)
      IF (MSHOW.GE.2) WRITE(MTERM,*)
C
      IF (MSHOW.GE.1) WRITE(MFILE,2)
      IF (MSHOW.GE.0) WRITE(MFILE,20001) ITNS,NSUBST,TIMENS,RELU2,RELP2,
     *                                   RELT
      IF (MSHOW.GE.1) WRITE(MFILE,2)
      IF (MSHOW.GE.1) WRITE(MFILE,*)
C
      IF ((RELT.LE.EPSNS).OR.(TIMENS.GE.TIMEMX)) GOTO 999
C
      CALL ZTIME(TTT1)
      TTTOTH=TTT1-TTTSUM
      IF (MSHOW.GE.3) WRITE(MTERM,*) 'total time : ', TTTOTH
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'total time : ', TTTOTH
      IF (MSHOW.GE.3) WRITE(MTERM,*) 'mavec time : ', TTADF
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'mavec time : ', TTADF
      IF (MSHOW.GE.3) WRITE(MTERM,*) 'konv. time : ', TTUPW
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'konv. time : ', TTUPW
      IF (MSHOW.GE.3) WRITE(MTERM,*) 'bdry  time : ', TTBDR
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'bdry  time : ', TTBDR
      IF (MSHOW.GE.3) WRITE(MTERM,*) 'LC    time : ', TTLC
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'LC    time : ', TTLC
      IF (MSHOW.GE.3) WRITE(MTERM,*) 'MG    time : ', TTMG
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'MG    time : ', TTMG
C
100   CONTINUE
C
      ITNS=NITNS
C
C=======================================================================
C *** Postprocessing
C=======================================================================
C
999   CALL ZTIME(TTT0)
C
      IFPOST=0
      CALL FPOST(IFPOST,IFILEN,ITFILM,UE,MFILE,MSHOW)
C
      CALL ZTIME(TTT1)
      TTPOST=TTPOST+TTT1-TTT0
C
C=======================================================================
C     Statistics
C=======================================================================
C
      IF (MSHOW.GE.2) WRITE (MTERM,*)
      IF (MSHOW.GE.0) WRITE (MFILE,*)
      IF (MSHOW.GE.2) WRITE (MTERM,*) 'STATISTICS :'
      IF (MSHOW.GE.0) WRITE (MFILE,*) 'STATISTICS :'
C
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'NWORK :      ',NWORK
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'NWORK :      ',NWORK
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'IWORKG:      ',IWORKG
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'IWORKG:      ',IWORKG
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'IWMAXG:      ',IWMAXG
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'IWMAXG:      ',IWMAXG
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'IWORKI:      ',IWORKI
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'IWORKI:      ',IWORKI
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'IWMAXI:      ',IWMAXI
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'IWMAXI:      ',IWMAXI
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'IWORK :      ',IWORK
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'IWORK :      ',IWORK
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'IWMAX :      ',IWMAX
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'IWMAX :      ',IWMAX
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.0) WRITE(MFILE,1)
C
      CALL ZTIME(TTT1)
      TTTSUM=TTT1-TTTSUM
      TTTH  =TTGRID+TTPOST+TTADF+TTUPW+TTBDR+TTLC+TTMG
      TTLIN =TTADF+TTUPW+TTBDR+TTLC
C
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'total time : ', TTTSUM
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'total time : ', TTTSUM
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'appr. time : ', TTTH
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'appr. time : ', TTTH
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'grid  time : ', TTGRID
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'grid  time : ', TTGRID
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'post  time : ', TTPOST
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'post  time : ', TTPOST
C
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'lin.  time : ', TTLIN
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'lin.  time : ', TTLIN
      IF (MSHOW.GE.2) WRITE(MTERM,*) '-> mavec time : ', TTADF
      IF (MSHOW.GE.0) WRITE(MFILE,*) '-> mavec time : ', TTADF
      IF (MSHOW.GE.2) WRITE(MTERM,*) '-> konv. time : ', TTUPW
      IF (MSHOW.GE.0) WRITE(MFILE,*) '-> konv. time : ', TTUPW
      IF (MSHOW.GE.2) WRITE(MTERM,*) '-> bdry  time : ', TTBDR
      IF (MSHOW.GE.0) WRITE(MFILE,*) '-> bdry  time : ', TTBDR
      IF (MSHOW.GE.2) WRITE(MTERM,*) '-> LC    time : ', TTLC
      IF (MSHOW.GE.0) WRITE(MFILE,*) '-> LC    time : ', TTLC
C
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'mg    time : ', TTMG
      IF (MSHOW.GE.1) WRITE(MFILE,*) 'mg    time : ', TTMG
C
      IF (MSHOW.GE.2) WRITE(MTERM,*) '#substeps  : ', NSUBST
      IF (MSHOW.GE.0) WRITE(MFILE,*) '#substeps  : ', NSUBST
      IF (MSHOW.GE.2) WRITE(MTERM,*) '#nonlinear : ', NNONL
      IF (MSHOW.GE.1) WRITE(MFILE,*) '#nonlinear : ', NNONL
      IF (MSHOW.GE.2) WRITE(MTERM,*) '#mg        : ', NMG
      IF (MSHOW.GE.1) WRITE(MFILE,*) '#mg        : ', NMG
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)
C
      IF (MSHOW.GE.2) THEN
      WRITE(MTERM,*)
      WRITE(MTERM,*) ' MULTIGRID COMPONENTS [in percent]:'
      WRITE(MTERM,*) ' smoothing     :', 1.D2*TTS/TTMG
      WRITE(MTERM,*) ' solver        :', 1.D2*TTE/TTMG
      WRITE(MTERM,*) ' defect calc.  :', 1.D2*TTD/TTMG
      WRITE(MTERM,*) ' prolongation  :', 1.D2*TTP/TTMG
      WRITE(MTERM,*) ' restriction   :', 1.D2*TTR/TTMG
      WRITE(MTERM,1)
      ENDIF
C
      IF (MSHOW.GE.0) THEN
      WRITE(MFILE,*)
      WRITE(MFILE,*) ' MULTIGRID COMPONENTS [in percent]:'
      WRITE(MFILE,*) ' smoothing     :', 1.D2*TTS/TTMG
      WRITE(MFILE,*) ' solver        :', 1.D2*TTE/TTMG
      WRITE(MFILE,*) ' defect calc.  :', 1.D2*TTD/TTMG
      WRITE(MFILE,*) ' prolongation  :', 1.D2*TTP/TTMG
      WRITE(MFILE,*) ' restriction   :', 1.D2*TTR/TTMG
      WRITE(MFILE,1)
      ENDIF
C
      GOTO 99999
C=======================================================================
C     Error case
C=======================================================================
99998 WRITE(MTERM,*) 'IER', IER
      WRITE(MTERM,*) 'IN SUBROUTINE ',SUB
C
99999 CLOSE(MFILE)
C
C
C
   1  FORMAT(80('-'))
   2  FORMAT(80('$'))
   3  FORMAT(80('@'))
1000  FORMAT (6E12.5)
1001  FORMAT(' IT DIV-L2',3X,'RHOMGP')
1003  FORMAT(I3,2(D9.2))
2001  FORMAT('MACRO STEP ',I4,' AT TIME = ',D10.3,
     *       ' WITH 3 STEPS: DT3 = ',D10.3)
2002  FORMAT('MACRO STEP ',I4,' AT TIME = ',D10.3,
     *       ' WITH 1 STEP : DT1 = ',D10.3)
3001  FORMAT('CHOICE ',I2,'(',I4,')',1X,'   ---- NEW DT = ',D10.3,
     *       ' -- OLD DT = ',D10.3)
10002 FORMAT ('OLD DT=',D9.2,
     *        1X,'U(L2)=',D9.2,1X,'U(MX)=',D9.2,
     *        1X,'P(L2)=',D9.2,1X,'P(MX)=',D9.2)
10003 FORMAT ('OLD DT=',D9.2,
     *        1X,'FWREL=',D9.2,1X,'AWREL=',D9.2,
     *        1X,'FW3  =',D9.2,1X,'FA3  =',D9.2)
20001 FORMAT ('#',I4,1X,'(',I4,')',1X,'TIME=',D10.3,1X,'RELU(L2)=',
     *        D9.2,1X,'RELP(L2)=',D9.2,1X,'REL=',D9.2)
C
C
C
      END
