************************************************************************
      SUBROUTINE FPOST(ITYP,IFILEN,ITFILM,UE,MFILE,MSHOW)
************************************************************************
*   Purpose: - performs the nonsteady postprocess:
*                   - output for film generation
*                   - output for AVS
*                   - output for GMV
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CXX(10)*15,CFILE*15
C
      DATA CXX/'#ns/#DX1       ','#ns/#DX2       ','#ns/#DX3       ',
     *         '#ns/#DX4       ','#ns/#DX5       ','#ns/#DX6       ',
     *         '#ns/#DX7       ','#ns/#DX8       ','#ns/#DX9       ',
     *         '#ns/#DX10      '/
C
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
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
      COMMON /NSPTS/  KPU(2),KPP(4),KPX(4),KPI(2),DPI(2,2),DPF(2)
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
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
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
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
      CHARACTER CPARM1*60,CMESH1*60,CFILE1*60,CSTART*60,CSOL*60
      COMMON /FILES/ IMESH1,MMESH1,CPARM1,CMESH1,MFILE1,CFILE1,
     *               ISTART,MSTART,CSTART,ISOL,MSOL,CSOL
C
      INCLUDE 'jump.inc'
      INCLUDE 'bouss.inc'
      EXTERNAL DVISCO,DVISCOPRIM
      SAVE 
C-----------------------------------------------------------------------
C
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
C *** Coefficient of exact solution and fem basis functions
      EXTERNAL UE,EM30,EM31
C=======================================================================
      SUB='FPOST '
C
C
C
      NGRAPH=MAX(IAVS,IGMV)
C
      ILEV=NLEV
      ISETLV=2
      CALL  SETLEV (ISETLV)
C
      CALL PTSDAT(TIMENS,DBLE(NY))
C
C=======================================================================
C *** write the solution vector if ABS(ISOL)=1
C=======================================================================
C
      IF ((ITYP.EQ.0).AND.(ABS(ISOL).EQ.1)) THEN
       IF (ISOL.EQ.1) THEN
        IFMTS=0
       ELSE
        IFMTS=1
       ENDIF
       CALL  OF0 (MSOL,CSOL,IFMTS)
       CALL  OWA1 (DWORK(KU1),'DU12P ',NUP,MSOL,IFMTS)
       CLOSE(MSOL)
      ENDIF
C
C=======================================================================
C *** write unformatted time dep. solution vector
C=======================================================================
C
      IF ((ITYP.EQ.1).AND.(INSAV.GT.0)) THEN
      IF (MOD(ITNS,INSAV).EQ.0) THEN
       IFILEN=IFILEN+1
       ITWX=MOD(IFILEN+INSAVN-1,INSAVN)+1
       CALL  OF0 (39,CXX(ITWX),0)
       CALL  OWA1 (DWORK(KU1),'DU12P ',NUP,39,0)
       REWIND(39)
       CLOSE (39)
       IF (IER.NE.0) GOTO 99999
      ENDIF
      ENDIF
C
C=======================================================================
C
      IF ((ITYP.EQ.1).AND.((TIMENS-DTFILO).GE.DTFILM)) THEN
       ITFILM=ITFILM+1
       WRITE(40,*) REAL(TIMENS)
      ENDIF
C
C=======================================================================
C *** write velocities 
C=======================================================================
C
       KMOV1=L(LD1)
       KMOV2=KMOV1+NVT
       KAUXM=L(KLAUX(NLEV))
C
       CALL  INTUVD(DWORK(KU1),DWORK(KU2),DWORK(KMOV1),DWORK(KMOV2),
     *              DWORK(KAUXM),NVT,NEL,NVBD,
     *              KWORK(L(LMID)),KWORK(L(LVERT)),KWORK(L(LVBD)),
     *              KWORK(L(KLMBD(NLEV))),DWORK(L(LCORVG)),UE)
C
       IF ((ITYP.EQ.1).AND.(IFUSAV.GT.0)
     *                .AND.((TIMENS-DTFILO).GE.DTFILM)) THEN
        CFILE='#film/#DU      '
        IF ((ITFILM.GE.0).AND.(ITFILM.LT.10)) 
     *       WRITE(CFILE(10:10),'(I1.1)') ITFILM
        IF ((ITFILM.GE.10).AND.(ITFILM.LT.100)) 
     *       WRITE(CFILE(10:11),'(I2.2)') ITFILM
        IF ((ITFILM.GE.100).AND.(ITFILM.LT.1000)) 
     *       WRITE(CFILE(10:12),'(I3.3)') ITFILM
        IF ((ITFILM.GE.1000).AND.(ITFILM.LT.10000)) 
     *       WRITE(CFILE(10:13),'(I4.4)') ITFILM
        IF (ITFILM.GE.10000) STOP
C
        CALL OF0 (39,CFILE,0)
        IF (IER.NE.0) GOTO 99998
        CALL  OWA1 (DWORK(KMOV1),'DUF   ',KNVT(IFUSAV),39,0)
        CALL  OWA1 (DWORK(KMOV2),'DUF   ',KNVT(IFUSAV),39,0)
        REWIND(39)
        CLOSE (39)
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (((ITYP.EQ.0).AND.(IAVS.GT.0)).OR.
     *     ((ITYP.EQ.0).AND.(IGMV.GT.0)).OR.
     *     ((ITYP.EQ.1).AND.(IAVS.GT.0)
     *                 .AND.((TIMENS-DTAVSO).GE.DTAVS)) .OR.
     *     ((ITYP.EQ.1).AND.(IGMV.GT.0)
     *                 .AND.((TIMENS-DTGMVO).GE.DTGMV))) THEN
        CALL ZNEW (2*KNVT(NGRAPH),-2,LAVSU,'AVSU  ')
        IF (IER.NE.0) GOTO 99998
        KAVS1=L(LAVSU)
        KAVS2=KAVS1+KNVT(NGRAPH)
        DO 100 IVTA=1,KNVT(NGRAPH)
        VWORK(KAVS1+IVTA-1)=REAL(DWORK(KMOV1+IVTA-1))
100     VWORK(KAVS2+IVTA-1)=REAL(DWORK(KMOV2+IVTA-1))
       ENDIF
C
      IF (ITYP.EQ.1) THEN
       IF (MSHOW.GE.2) 
     *  WRITE(MTERM,1001) DWORK(KMOV1+KPU(1)-1),DWORK(KMOV2+KPU(1)-1),
     *                    DWORK(KMOV1+KPU(2)-1),DWORK(KMOV2+KPU(2)-1)
       IF (MSHOW.GE.1) 
     *  WRITE(MFILE,1001) DWORK(KMOV1+KPU(1)-1),DWORK(KMOV2+KPU(1)-1),
     *                    DWORK(KMOV1+KPU(2)-1),DWORK(KMOV2+KPU(2)-1)
       WRITE(41,*) REAL(TIMENS),REAL(DWORK(KMOV1+KPU(1)-1))
       WRITE(42,*) REAL(TIMENS),REAL(DWORK(KMOV2+KPU(1)-1))
       WRITE(43,*) REAL(TIMENS),REAL(DWORK(KMOV1+KPU(2)-1))
       WRITE(44,*) REAL(TIMENS),REAL(DWORK(KMOV2+KPU(2)-1))
      ENDIF
C
C=======================================================================
C *** write PRESSURE 
C=======================================================================
C

       KPL  =L(LD1)
       KAUXM=L(LD2)
       LAREA=KLAREA(NLEV)
C
       CALL  INTPV (DWORK(KP),DWORK(KPL),DWORK(KAUXM),VWORK(L(LAREA)),
     *              KWORK(L(LVERT)))
       IF (IER.NE.0) GOTO 99999
C
       IF ((ITYP.EQ.1).AND.(IFPSAV.GT.0)
     *                .AND.((TIMENS-DTFILO).GE.DTFILM)) THEN
        CFILE='#film/#DP      '
        IF ((ITFILM.GE.0).AND.(ITFILM.LT.10)) 
     *       WRITE(CFILE(10:10),'(I1.1)') ITFILM
        IF ((ITFILM.GE.10).AND.(ITFILM.LT.100)) 
     *       WRITE(CFILE(10:11),'(I2.2)') ITFILM
        IF ((ITFILM.GE.100).AND.(ITFILM.LT.1000)) 
     *       WRITE(CFILE(10:12),'(I3.3)') ITFILM
        IF ((ITFILM.GE.1000).AND.(ITFILM.LT.10000)) 
     *       WRITE(CFILE(10:13),'(I4.4)') ITFILM
        IF (ITFILM.GE.10000) STOP
C
        CALL OF0 (39,CFILE,0)
        IF (IER.NE.0) GOTO 99998
        CALL  OWA1 (DWORK(KPL),'DPL   ',KNVT(IFPSAV),39,0)
        REWIND(39)
        CLOSE (39)
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (((ITYP.EQ.0).AND.(IAVS.GT.0)).OR.
     *     ((ITYP.EQ.0).AND.(IGMV.GT.0)).OR.
     *     ((ITYP.EQ.1).AND.(IAVS.GT.0)
     *                 .AND.((TIMENS-DTAVSO).GE.DTAVS)) .OR.
     *     ((ITYP.EQ.1).AND.(IGMV.GT.0)
     *                 .AND.((TIMENS-DTGMVO).GE.DTGMV))) THEN
        CALL ZNEW (KNVT(NGRAPH),-2,LAVSP,'AVSP  ')
        IF (IER.NE.0) GOTO 99998
        KAVSP=L(LAVSP)
        DO 200 IVTA=1,KNVT(NGRAPH)
200     VWORK(KAVSP+IVTA-1)=REAL(DWORK(KPL+IVTA-1))
       ENDIF
C
      IF (ITYP.EQ.1) THEN
       IF (MSHOW.GE.2) 
     *  WRITE(MTERM,1002) DWORK(KPL+KPP(1)-1),DWORK(KPL+KPP(2)-1),
     *                    DWORK(KPL+KPP(3)-1),DWORK(KPL+KPP(4)-1)
       IF (MSHOW.GE.1) 
     *  WRITE(MFILE,1002) DWORK(KPL+KPP(1)-1),DWORK(KPL+KPP(2)-1),
     *                    DWORK(KPL+KPP(3)-1),DWORK(KPL+KPP(4)-1)
       WRITE(45,*) REAL(TIMENS),REAL(DWORK(KPL+KPP(1)-1))
       WRITE(46,*) REAL(TIMENS),REAL(DWORK(KPL+KPP(2)-1))
       WRITE(47,*) REAL(TIMENS),REAL(DWORK(KPL+KPP(3)-1))
       WRITE(48,*) REAL(TIMENS),REAL(DWORK(KPL+KPP(4)-1))
      ENDIF
C
       CALL  BDPRES(DWORK(KPL),KWORK(L(LVERT)),KWORK(L(LNPR)),
     *              KWORK(L(LVBD)),KWORK(L(LMM)),
     *              DWORK(L(LCORVG)),DWORK(L(LVBDP)),
     *              KNVBD(NLEV),P5,P6)
C
      IF (ITYP.EQ.1) THEN
       IF (MSHOW.GE.2) 
     *  WRITE(MTERM,1003) P5,P6
       IF (MSHOW.GE.1) 
     *  WRITE(MFILE,1003) P5,P6
       WRITE(49,*) REAL(TIMENS),REAL(P5)
       WRITE(50,*) REAL(TIMENS),REAL(P6)
      ENDIF
C
C
C
C=======================================================================
C *** write lift and drag as volume integral 
C=======================================================================
C

c        KNY=L(KLNY(NLEV))
      DO I=1,NA
         VWORK(0*NA+KA1+I-1)=0.0D0
         VWORK(1*NA+KA1+I-1)=0.0D0
         VWORK(2*NA+KA1+I-1)=0.0D0
         VWORK(3*NA+KA1+I-1)=0.0D0
      ENDDO   
      ISTOKOLD=ISTOK
      IJUMPOLD=IJUMP
      IUSENEWTOLD=IUSENEWT
      DNPWOLD=DNPW
      DNNSOLD=DNNS
      DNPWSTARTold=DNPWSTART
      DNPOLD=DNP
      THSTEPOLD=THSTEP
      IMASMATRIXOLD=IMASMATRIX
      THSTEP=1.0D0
      IMASMATRIX=0
      DNNS=0.0d0
      DNPW=0.0d0
      DNPWSTART=0.0d0
      DNP=0.0D0
      IUSENEWT=0
      IJUMP=0
      ISTOK=1
      
        CALL STABIL   (DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(L(LD1)),DWORK(L(LD2)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *          VWORK(KB1),VWORK(KB2),NB,
     *          VWORK(KB1+NB),VWORK(KB2+NB),KWORK(KCOLB),KWORK(KLDB),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),
     *                KWORK(L(LADJ)),KWORK(L(LMEL)),
     *                DWORK(L(LCORVG)),EM30,COEFFN,0,1D0)
C
       IUSENEWT=IUSENEWTOLD
       ISTOK=ISTOKOLD
       IJUMP=IJUMPOLD
       DNPW=DNPWOLD 
       DNNS=DNNSOLD
       DNP=DNPOLD
       DNPWSTART=DNPWSTARTold
       THSTEP=THSTEPOLD
       IMASMATRIX=1
c
c$$$       print*,'A1',(VWORK(0*NA+KA1+I-1),I=1,5)
c$$$       print*,'A1',(VWORK(1*NA+KA1+I-1),I=1,5)
c$$$       print*,'A1',(VWORK(2*NA+KA1+I-1),I=1,5)
c$$$       print*,'A4',(VWORK(3*NA+KA1+I-1),I=1,5)
         CALL ZNEW (NU,-1,lvd1,'vd1  ')
         CALL ZNEW (NU,-1,lvd2,'vd2  ')
         CALL ZNEW (NU,-1,lvl1,'vl1  ')
         CALL ZNEW (NU,-1,lvl2,'vl2  ')
         CALL ZNEW (NU,-1,ldau1,'dau1 ')
         CALL ZNEW (NU,-1,ldau2,'dau2 ')
         CALL ZNEW (NP,-1,lDVDBU,'lDVDBU')
         CALL ZNEW (NP,-1,lDVLBU,'lDVLBU')
C
         IF ((IELT.EQ.0).OR.(IELT.EQ.2)) 
     *   CALL VIBDFORC(DWORK(KU1),DWORK(KU2),DWORK(KP),
     *      KWORK(L(LVERT)),VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *      VWORK(KB1),VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *      KWORK(L(LNPR)),KWORK(L(LMID)),KWORK(L(LVBD)),KWORK(L(LEBD)),
     *      KWORK(L(LMM)),DWORK(L(LCORVG)),EM31,DFWVI,DAWVI,
     *      dviappp,dvifppp,dviauuu,dvifuuu,
     *      dwork(l(lVD1)),dwork(l(lVD2)),
     *      dwork(l(lVl1)),dwork(l(lVl2)),
     *      dwork(l(ldau1)),dwork(l(ldau2)),
     *      dwork(l(lDVDBU)),dwork(l(lDVLBU))) 
       IF ((IELT.EQ.1).OR.(IELT.EQ.3)) 
     *   CALL VIBDFORC(DWORK(KU1),DWORK(KU2),DWORK(KP),
     *      KWORK(L(LVERT)),VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *      VWORK(KB1),VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *      KWORK(L(LNPR)),KWORK(L(LMID)),KWORK(L(LVBD)),KWORK(L(LEBD)),
     *      KWORK(L(LMM)),DWORK(L(LCORVG)),EM30,DFWVI,DAWVI,
     *      dviappp,dvifppp,dviauuu,dvifuuu,
     *      dwork(l(lVD1)),dwork(l(lVD2)),
     *      dwork(l(lVl1)),dwork(l(lVl2)),
     *      dwork(l(ldau1)),dwork(l(ldau2)),
     *      dwork(l(lDVDBU)),dwork(l(lDVLBU))) 
C
         CALL ZDISP (0,lvd1,'vd1  ')
         CALL ZDISP (0,lvd2,'vd2  ')
         CALL ZDISP (0,lvl1,'vl1  ')
         CALL ZDISP (0,lvl2,'vl2  ')
         CALL ZDISP (0,ldau1,'dau1 ')
         CALL ZDISP (0,ldau2,'dau2 ')
         CALL ZDISP (0,lDVDBU,'lDVDBU')
         CALL ZDISP (0,lDVLBU,'lDVLBU')
         IF (ITYP.EQ.1) THEN
            IF (MSHOW.GE.2) THEN
               WRITE(MTERM,1004) DAWVI,DFWVI
            ENDIF     
            
            IF (MSHOW.GE.1) THEN
               WRITE(MFILE,1004) DAWVI,DFWVI
            ENDIF 
            
c            WRITE(51,*) REAL(TIMENS),REAL(DFW),REAL(DAWVI)
c            WRITE(52,*) REAL(TIMENS),REAL(DAW),REAL(DFWVI)
         ENDIF
C
C=======================================================================
C *** write lift and drag 
C=======================================================================
C
       IF ((IELT.EQ.0).OR.(IELT.EQ.2)) 
     *  CALL  BDFORC(DWORK(KU1),DWORK(KU2),DWORK(KP),KWORK(L(LVERT)),
     *               KWORK(L(LMID)),KWORK(L(LVBD)),KWORK(L(LEBD)),
     *               KWORK(L(LMM)),DWORK(L(LCORVG)),EM31,DFW,DAW)
C
       IF ((IELT.EQ.1).OR.(IELT.EQ.3)) 
     *  CALL  BDFORC(DWORK(KU1),DWORK(KU2),DWORK(KP),KWORK(L(LVERT)),
     *               KWORK(L(LMID)),KWORK(L(LVBD)),KWORK(L(LEBD)),
     *               KWORK(L(LMM)),DWORK(L(LCORVG)),EM30,DFW,DAW)
C
      IF (ITYP.EQ.1) THEN
       IF (MSHOW.GE.2) 
     *  WRITE(MTERM,1004) DFW,DAW
       IF (MSHOW.GE.1) 
     *  WRITE(MFILE,1004) DFW,DAW
c       WRITE(51,*) REAL(TIMENS),REAL(DFW)
c       WRITE(52,*) REAL(TIMENS),REAL(DAW)
      ENDIF
C
C=======================================================================
C *** write lift and drag 
C=======================================================================
C
       IF ((IELT.EQ.0).OR.(IELT.EQ.2)) 
     *  CALL  VIBF_STG(DWORK(KU1),DWORK(KU2),DWORK(KP),KWORK(L(LVERT)),
     *   KWORK(L(LNPR)),KWORK(L(LMID)),DWORK(L(LCORVG)),EM31,DFWI,DAWI)
C
       IF ((IELT.EQ.1).OR.(IELT.EQ.3)) 
     *  CALL  VIBF_STG(DWORK(KU1),DWORK(KU2),DWORK(KP),KWORK(L(LVERT)),
     *   KWORK(L(LNPR)),KWORK(L(LMID)),DWORK(L(LCORVG)),EM30,DFWI,DAWI)
C
      IF (ITYP.EQ.1) THEN
       IF (MSHOW.GE.2) 
     *  WRITE(MTERM,1004) DFWI,DAWI
       IF (MSHOW.GE.1) 
     *  WRITE(MFILE,1004) DFWI,DAWI
       WRITE(51,*) REAL(TIMENS),REAL(DFW),REAL(DFWI)
       WRITE(52,*) REAL(TIMENS),REAL(DAW),REAL(DAWI)
      ENDIF
C
C=======================================================================
C *** write streamfunction 
C=======================================================================
C
       KISO =L(LD1)
       KVIND=L(LD2)
C
       CALL  U2ISO (DWORK(L(LCORVG)),KWORK(L(LVERT)),KWORK(L(LMID)),
     *              KWORK(L(LADJ)),DWORK(KVIND),DWORK(KISO),
     *              DWORK(KU1),DWORK(KU2))
       IF (IER.NE.0) GOTO 99999
C
       IF ((ITYP.EQ.1).AND.(IFXSAV.GT.0)
     *                .AND.((TIMENS-DTFILO).GE.DTFILM)) THEN
        CFILE='#film/#DX      '
        IF ((ITFILM.GE.0).AND.(ITFILM.LT.10)) 
     *       WRITE(CFILE(10:10),'(I1.1)') ITFILM
        IF ((ITFILM.GE.10).AND.(ITFILM.LT.100)) 
     *       WRITE(CFILE(10:11),'(I2.2)') ITFILM
        IF ((ITFILM.GE.100).AND.(ITFILM.LT.1000)) 
     *       WRITE(CFILE(10:12),'(I3.3)') ITFILM
        IF ((ITFILM.GE.1000).AND.(ITFILM.LT.10000)) 
     *       WRITE(CFILE(10:13),'(I4.4)') ITFILM
        IF (ITFILM.GE.10000) STOP
C
        CALL OF0 (39,CFILE,0)
        IF (IER.NE.0) GOTO 99998
        CALL  OWA1 (DWORK(KISO),'DISO  ',KNVT(IFXSAV),39,0)
        REWIND(39)
        CLOSE (39)
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (((ITYP.EQ.0).AND.(IAVS.GT.0)).OR.
     *     ((ITYP.EQ.0).AND.(IGMV.GT.0)).OR.
     *     ((ITYP.EQ.1).AND.(IAVS.GT.0)
     *                 .AND.((TIMENS-DTAVSO).GE.DTAVS)) .OR.
     *     ((ITYP.EQ.1).AND.(IGMV.GT.0)
     *                 .AND.((TIMENS-DTGMVO).GE.DTGMV))) THEN
        CALL ZNEW (KNVT(NGRAPH),-2,LAVSI,'AVSI  ')
        IF (IER.NE.0) GOTO 99998
        KAVSI=L(LAVSI)
        DO 300 IVTA=1,KNVT(NGRAPH)
300     VWORK(KAVSI+IVTA-1)=REAL(DWORK(KISO+IVTA-1))
       ENDIF
C
      IF (ITYP.EQ.1) THEN
       IF (MSHOW.GE.2) 
     *  WRITE(MTERM,1005) DWORK(KISO+KPX(1)-1)-DWORK(KISO+KPX(2)-1)
       IF (MSHOW.GE.1) 
     *  WRITE(MFILE,1005) DWORK(KISO+KPX(1)-1)-DWORK(KISO+KPX(2)-1)
       WRITE(53,*) REAL(TIMENS),REAL( DWORK(KISO+KPX(1)-1)
     *                               -DWORK(KISO+KPX(2)-1))
      ENDIF
C
C=======================================================================
C
      IF ((ITYP.EQ.1).AND.((TIMENS-DTFILO).GE.DTFILM)) THEN
       DTFILO=TIMENS
      ENDIF
C
C=======================================================================
C *** write viscosity 
C=======================================================================
C
       KNY=L(KLNY(NLEV))
       KNYL  =L(LD1)
       KAUXM =L(LD2)
       LAREA=KLAREA(NLEV)
C
       CALL  INTPV (DWORK(KNY),DWORK(KNYL),DWORK(KAUXM),
     *              VWORK(L(LAREA)),KWORK(L(LVERT)))
       IF (IER.NE.0) GOTO 99999
C
       IF (((ITYP.EQ.0).AND.(IAVS.GT.0)).OR.
     *     ((ITYP.EQ.0).AND.(IGMV.GT.0)).OR.
     *     ((ITYP.EQ.1).AND.(IAVS.GT.0)
     *                 .AND.((TIMENS-DTAVSO).GE.DTAVS)) .OR.
     *     ((ITYP.EQ.1).AND.(IGMV.GT.0)
     *                 .AND.((TIMENS-DTGMVO).GE.DTGMV))) THEN
        CALL ZNEW (KNVT(NGRAPH),-2,LAVSNY,'AVSNY ')
        IF (IER.NE.0) GOTO 99998
        KAVSNY=L(LAVSNY)
        DO 400 IVTA=1,KNVT(NGRAPH)
           DWORK(KNYL+IVTA-1)=DWORK(KNYL+IVTA-1)*NY
           VWORK(KAVSNY+IVTA-1)=REAL(DWORK(KNYL+IVTA-1))
 400    CONTINUE

       ENDIF
C
CCC
C=======================================================================
C *** write norm of grad function 
C=======================================================================
c
C
       CALL ZNEW(NEL,1,LG,'DNG')
       KG=L(LG)
C
       IF (IELT.EQ.0 .OR. IELT.EQ.2) 
     * CALL NGRAD(DWORK(KU1),DWORK(KU2),DWORK(KG),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),
     *                 DWORK(L(LCORVG)),EM31)
       IF (IELT.EQ.1 .OR. IELT.EQ.3) 
     * CALL NGRAD(DWORK(KU1),DWORK(KU2),DWORK(KG),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),
     *                 DWORK(L(LCORVG)),EM30)
c
c
       KGG= L(LD1)
       KAUXM=L(LD2)
       LAREA=KLAREA(NLEV)
C
       CALL  INTPV (DWORK(KG),DWORK(KGG),DWORK(KAUXM),
     *             VWORK(L(LAREA)),KWORK(L(LVERT)))
       IF (IER.NE.0) GOTO 99999
C
       CALL ZDISP(0,LG,'DNG')
C
c
       IF (((ITYP.EQ.0).AND.(IAVS.GT.0)).OR.
     *     ((ITYP.EQ.0).AND.(IGMV.GT.0)).OR.
     *     ((ITYP.EQ.1).AND.(IAVS.GT.0)
     *                 .AND.((TIMENS-DTAVSO).GE.DTAVS)) .OR.
     *     ((ITYP.EQ.1).AND.(IGMV.GT.0)
     *                 .AND.((TIMENS-DTGMVO).GE.DTGMV))) THEN
        CALL ZNEW (KNVT(NGRAPH),-2,LAVSG,'AVSG  ')
        IF (IER.NE.0) GOTO 99998
        KAVSG=L(LAVSG)
        DO  IVTA=1,KNVT(NGRAPH)
              VWORK(KAVSG+IVTA-1)=REAL(DWORK(KGG+IVTA-1))
        ENDDO
       ENDIF
C
c
C=======================================================================
C *** write TENSOR 
C=======================================================================
c
C
       CALL ZNEW(NEL,1,LT,'DT')
       KT=L(LT)
       KNY=L(KLNY(NLEV))
C
       IF (IELT.EQ.0 .OR. IELT.EQ.2) 
     * CALL TENSOR(DWORK(KU1),DWORK(KU2),DWORK(KP),DWORK(KNY),
     *              DWORK(KT),KWORK(L(LVERT)),KWORK(L(LMID)),
     *                 DWORK(L(LCORVG)),EM31)
       IF (IELT.EQ.1 .OR. IELT.EQ.3) 
     * CALL TENSOR(DWORK(KU1),DWORK(KU2),DWORK(KP),DWORK(KNY),
     *              DWORK(KT),KWORK(L(LVERT)),KWORK(L(LMID)),
     *                 DWORK(L(LCORVG)),EM30)

c
C
       KTT= L(LD1)
       KAUXM=L(LD2)
       LAREA=KLAREA(NLEV)
C
       CALL  INTPV (DWORK(KT),DWORK(KTT),DWORK(KAUXM),
     *             VWORK(L(LAREA)),KWORK(L(LVERT)))
       IF (IER.NE.0) GOTO 99999
C
       CALL ZDISP(0,LT,'DT')
c
       IF (((ITYP.EQ.0).AND.(IAVS.GT.0)).OR.
     *     ((ITYP.EQ.0).AND.(IGMV.GT.0)).OR.
     *     ((ITYP.EQ.1).AND.(IAVS.GT.0)
     *                 .AND.((TIMENS-DTAVSO).GE.DTAVS)) .OR.
     *     ((ITYP.EQ.1).AND.(IGMV.GT.0)
     *                 .AND.((TIMENS-DTGMVO).GE.DTGMV))) THEN
        CALL ZNEW (KNVT(NGRAPH),-2,LAVST,'AVST  ')
        IF (IER.NE.0) GOTO 99998
        KAVST=L(LAVST)
        DO  IVTA=1,KNVT(NGRAPH)
              VWORK(KAVST+IVTA-1)=REAL(DWORK(KTT+IVTA-1))
        ENDDO
       ENDIF
C
C=======================================================================
C *** write TENSOR 
C=======================================================================
c
C
       CALL ZNEW(NEL,1,LT1,'DT1')
       KT1=L(LT1)
       KNY=L(KLNY(NLEV))
C
       IF (IELT.EQ.0 .OR. IELT.EQ.2) 
     * CALL TENSOR1(DWORK(KU1),DWORK(KU2),DWORK(KP),DWORK(KNY),
     *              DWORK(KT1),KWORK(L(LVERT)),KWORK(L(LMID)),
     *                 DWORK(L(LCORVG)),EM31)
       IF (IELT.EQ.1 .OR. IELT.EQ.3) 
     * CALL TENSOR1(DWORK(KU1),DWORK(KU2),DWORK(KP),DWORK(KNY),
     *              DWORK(KT1),KWORK(L(LVERT)),KWORK(L(LMID)),
     *                 DWORK(L(LCORVG)),EM30)

c
C
       KT1T1= L(LD1)
       KAUXM=L(LD2)
       LAREA=KLAREA(NLEV)
C
       CALL  INTPV (DWORK(KT1),DWORK(KT1T1),DWORK(KAUXM),
     *             VWORK(L(LAREA)),KWORK(L(LVERT)))
       IF (IER.NE.0) GOTO 99999
C
       CALL ZDISP(0,LT1,'DT1')
c
       IF (((ITYP.EQ.0).AND.(IAVS.GT.0)).OR.
     *     ((ITYP.EQ.0).AND.(IGMV.GT.0)).OR.
     *     ((ITYP.EQ.1).AND.(IAVS.GT.0)
     *                 .AND.((TIMENS-DTAVSO).GE.DTAVS)) .OR.
     *     ((ITYP.EQ.1).AND.(IGMV.GT.0)
     *                 .AND.((TIMENS-DTGMVO).GE.DTGMV))) THEN
        CALL ZNEW (KNVT(NGRAPH),-2,LAVST1,'AVST1  ')
        IF (IER.NE.0) GOTO 99998
        KAVST1=L(LAVST1)
        DO  IVTA=1,KNVT(NGRAPH)
              VWORK(KAVST1+IVTA-1)=REAL(DWORK(KT1T1+IVTA-1))
        ENDDO
       ENDIF
C
C=======================================================================
C *** write TENSOR 
C=======================================================================
c
C
       CALL ZNEW(NEL,1,LT2,'DT2')
       KT2=L(LT2)
       KNY=L(KLNY(NLEV))
C
       IF (IELT.EQ.0 .OR. IELT.EQ.2) 
     * CALL TENSOR2(DWORK(KU1),DWORK(KU2),DWORK(KP),DWORK(KNY),
     *              DWORK(KT2),KWORK(L(LVERT)),KWORK(L(LMID)),
     *                 DWORK(L(LCORVG)),EM31)
       IF (IELT.EQ.1 .OR. IELT.EQ.3) 
     * CALL TENSOR2(DWORK(KU1),DWORK(KU2),DWORK(KP),DWORK(KNY),
     *              DWORK(KT2),KWORK(L(LVERT)),KWORK(L(LMID)),
     *                 DWORK(L(LCORVG)),EM30)

c
C
       KT2T2= L(LD1)
       KAUXM=L(LD2)
       LAREA=KLAREA(NLEV)
C
       CALL  INTPV (DWORK(KT2),DWORK(KT2T2),DWORK(KAUXM),
     *             VWORK(L(LAREA)),KWORK(L(LVERT)))
       IF (IER.NE.0) GOTO 99999
C
       CALL ZDISP(0,LT2,'DT2')
c
       IF (((ITYP.EQ.0).AND.(IAVS.GT.0)).OR.
     *     ((ITYP.EQ.0).AND.(IGMV.GT.0)).OR.
     *     ((ITYP.EQ.1).AND.(IAVS.GT.0)
     *                 .AND.((TIMENS-DTAVSO).GE.DTAVS)) .OR.
     *     ((ITYP.EQ.1).AND.(IGMV.GT.0)
     *                 .AND.((TIMENS-DTGMVO).GE.DTGMV))) THEN
        CALL ZNEW (KNVT(NGRAPH),-2,LAVST2,'AVST2  ')
        IF (IER.NE.0) GOTO 99998
        KAVST2=L(LAVST2)
        DO  IVTA=1,KNVT(NGRAPH)
              VWORK(KAVST2+IVTA-1)=REAL(DWORK(KT2T2+IVTA-1))
        ENDDO
       ENDIF
C

c
C=======================================================================
C *** write AVS and/or GMV
C=======================================================================
C
      IF (((ITYP.EQ.0).AND.(IAVS.GT.0)).OR.
     *    ((ITYP.EQ.0).AND.(IGMV.GT.0)).OR.
     *    ((ITYP.EQ.1).AND.(IAVS.GT.0)
     *                .AND.((TIMENS-DTAVSO).GE.DTAVS)) .OR.
     *    ((ITYP.EQ.1).AND.(IGMV.GT.0)
     *                .AND.((TIMENS-DTGMVO).GE.DTGMV))) THEN
C
       IF (((ITYP.EQ.0).AND.(IAVS.GT.0)).OR.
     *     ((ITYP.EQ.1).AND.(IAVS.GT.0)
     *                 .AND.((TIMENS-DTAVSO).GE.DTAVS))) THEN
       CFILE='#avs/u.        '
       IF ((ITNS+IFINIT.GE.0).AND.(ITNS+IFINIT.LT.10)) 
     *      WRITE(CFILE(8:8),'(I1.1)') ITNS+IFINIT
       IF ((ITNS+IFINIT.GE.10).AND.(ITNS+IFINIT.LT.100)) 
     *      WRITE(CFILE(8:9),'(I2.2)') ITNS+IFINIT
       IF ((ITNS+IFINIT.GE.100).AND.(ITNS+IFINIT.LT.1000)) 
     *      WRITE(CFILE(8:10),'(I3.3)') ITNS+IFINIT
       IF ((ITNS+IFINIT.GE.1000).AND.(ITNS+IFINIT.LT.10000)) 
     *      WRITE(CFILE(8:11),'(I4.4)') ITNS+IFINIT
       IF (ITNS+IFINIT.GE.10000) STOP
C
       CALL XAVS2D(39,CFILE,KNEL(IAVS),KNVT(IAVS),
     *             KWORK(L(KLVERT(IAVS))),DWORK(L(LCORVG)),
     *             VWORK(KAVS1),VWORK(KAVS2),VWORK(KAVSP),
     *             VWORK(KAVSI))
       DTAVSO=TIMENS
       ENDIF
C
       IF (((ITYP.EQ.0).AND.(IGMV.GT.0)).OR.
     *     ((ITYP.EQ.1).AND.(IGMV.GT.0)
     *                 .AND.((TIMENS-DTGMVO).GE.DTGMV))) THEN
       CFILE='#gmv/u.        '
       IF ((ITNS+IFINIT.GE.0).AND.(ITNS+IFINIT.LT.10)) 
     *      WRITE(CFILE(8:8),'(I1.1)') ITNS+IFINIT
       IF ((ITNS+IFINIT.GE.10).AND.(ITNS+IFINIT.LT.100)) 
     *      WRITE(CFILE(8:9),'(I2.2)') ITNS+IFINIT
       IF ((ITNS+IFINIT.GE.100).AND.(ITNS+IFINIT.LT.1000)) 
     *      WRITE(CFILE(8:10),'(I3.3)') ITNS+IFINIT
       IF ((ITNS+IFINIT.GE.1000).AND.(ITNS+IFINIT.LT.10000)) 
     *      WRITE(CFILE(8:11),'(I4.4)') ITNS+IFINIT
       IF (ITNS+IFINIT.GE.10000) STOP
C
       CALL XGMV2D(39,CFILE,KNEL(IGMV),KNVT(IGMV),
     *             KWORK(L(KLVERT(IGMV))),DWORK(L(LCORVG)),
     *             VWORK(KAVS1),VWORK(KAVS2),VWORK(KAVSP),
     *             VWORK(KAVSI),VWORK(KAVSNY),VWORK(KAVSG),
     *             VWORK(KAVST),VWORK(KAVST1),VWORK(KAVST2),TIMENS)
       DTGMVO=TIMENS
       ENDIF
C
       CALL ZDISP(0,LAVSI,'AVSI  ')
       CALL ZDISP(0,LAVSP,'AVSP  ')
       CALL ZDISP(0,LAVSU,'AVSU  ')
       CALL ZDISP(0,LAVSNY,'AVSNY ')
       CALL ZDISP(0,LAVSG,'AVSG  ')
       CALL ZDISP(0,LAVST,'AVST  ')
       CALL ZDISP(0,LAVST1,'AVST1  ')
       CALL ZDISP(0,LAVST2,'AVST2  ')
       IF (IER.NE.0) GOTO 99998
      ENDIF
C
C
      GOTO 99999
C
C=======================================================================
C     Error case
C=======================================================================
99998 WRITE(MTERM,*) 'IER', IER
      WRITE(MTERM,*) 'IN SUBROUTINE ',SUB
C
1000  FORMAT (6E12.5)
1001  FORMAT ('P(VELO) ',4(D12.5))
1002  FORMAT ('P(PRES) ',4(D12.5))
1003  FORMAT ('I(PRES) ',2(D12.5))
1004  FORMAT ('I(FORCE)',2(D12.5))
1005  FORMAT ('P(FLUX) ',1(D12.5))
c
c
c
99999 END

