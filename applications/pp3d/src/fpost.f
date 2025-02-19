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
C
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CXX(10)*15,CFILE*15
      DATA CXX/'#ns/DX1        ','#ns/DX2        ','#ns/DX3        ',
     *         '#ns/DX4        ','#ns/DX5        ','#ns/DX6        ',
     *         '#ns/DX7        ','#ns/DX8        ','#ns/DX9        ',
     *         '#ns/DX10       '/
C
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
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
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
      COMMON /NSSAV/  INSAV,INSAVN
      COMMON /NSSAVF/ DTFILM,DTFILO,DTAVS,DTAVSO,DTGMV,DTGMVO,
     *                IFUSAV,IFPSAV,IFXSAV,IAVS,IGMV,IFINIT
      COMMON /NSPTS/  KPU(2),KPP(4)
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
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
      COMMON /MGBDRY/ INEUM,LELBD,KELBD(NNLEV),KLNPRO(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C
C *** user COMMON block
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
      CHARACTER CPARM1*60,CMESH1*60,CFILE1*60,CSTART*60,CSOL*60
      COMMON /FILES/ IMESH1,MMESH1,CPARM1,CMESH1,MFILE1,CFILE1,
     *               ISTART,MSTART,CSTART,ISOL,MSOL,CSOL
C
      SAVE 
C-----------------------------------------------------------------------
C
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
C *** Coefficient of exact solution and element basis functions
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
      IF ((TIMENS-DTFILO).GE.DTFILM) THEN
       ITFILM=ITFILM+1
       WRITE(40,*) REAL(TIMENS),REAL(TIMENS)
      ENDIF
C
C=======================================================================
C *** write velocities 
C=======================================================================
C
      KMOV1=L(LD1)
      KMOV2=KMOV1+NVT
      KMOV3=KMOV2+NVT
      KAUXM=L(KLAUX(NLEV))
C
      CALL  INTUVD(DWORK(KU1)  ,DWORK(KU2)  ,DWORK(KU3)  ,        
     *             DWORK(KMOV1),DWORK(KMOV2),DWORK(KMOV3),
     *             DWORK(KAUXM), NVT, NEL, NVBD,NABD,KWORK(L(LAREA)),
     *             KWORK(L(LVERT)),KWORK(L(LNPR)),KWORK(L(LVBD)),
     *             KWORK(L(KLABD(NLEV))),DWORK(L(LCORVG)),
     *             UE,INEUM,KWORK(L(KELBD(NLEV))))
C
      IF ((ITYP.EQ.1).AND.(IFUSAV.GT.0)
     *               .AND.((TIMENS-DTFILO).GE.DTFILM)) THEN
       CFILE='#film/#DU      '
       IF ((ITFILM.GE.0).AND.(ITFILM.LT.10)) 
     *      WRITE(CFILE(10:10),'(I1.1)') ITFILM
       IF ((ITFILM.GE.10).AND.(ITFILM.LT.100)) 
     *      WRITE(CFILE(10:11),'(I2.2)') ITFILM
       IF ((ITFILM.GE.100).AND.(ITFILM.LT.1000)) 
     *      WRITE(CFILE(10:12),'(I3.3)') ITFILM
       IF ((ITFILM.GE.1000).AND.(ITFILM.LT.10000)) 
     *      WRITE(CFILE(10:13),'(I4.4)') ITFILM
       IF (ITFILM.GE.10000) STOP
C
       CALL OF0 (39,CFILE,0)
       IF (IER.NE.0) GOTO 99998
       CALL  OWA1 (DWORK(KMOV1),'DUF   ',KNVT(IFUSAV),39,0)
       CALL  OWA1 (DWORK(KMOV2),'DUF   ',KNVT(IFUSAV),39,0)
       CALL  OWA1 (DWORK(KMOV3),'DUF   ',KNVT(IFUSAV),39,0)
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
         CALL ZNEW (3*KNVT(NGRAPH),-2,LAVSU,'AVSU  ')
         IF (IER.NE.0) GOTO 99998
         KAVS1=L(LAVSU)
         KAVS2=KAVS1+KNVT(NGRAPH)
         KAVS3=KAVS2+KNVT(NGRAPH)
         DO 100 IVTA=1,KNVT(NGRAPH)
         VWORK(KAVS1+IVTA-1)=REAL(DWORK(KMOV1+IVTA-1))
         VWORK(KAVS2+IVTA-1)=REAL(DWORK(KMOV2+IVTA-1))
100      VWORK(KAVS3+IVTA-1)=REAL(DWORK(KMOV3+IVTA-1))
        ENDIF
C
      IF (ITYP.EQ.1) THEN
       IF (MSHOW.GE.2) WRITE(MTERM,1001) 
     *  DWORK(KMOV1+KPU(1)-1),DWORK(KMOV2+KPU(1)-1),
     *  DWORK(KMOV3+KPU(1)-1),DWORK(KMOV1+KPU(2)-1),
     *  DWORK(KMOV2+KPU(2)-1),DWORK(KMOV3+KPU(2)-1)
       IF (MSHOW.GE.1) WRITE(MFILE,1001) 
     *  DWORK(KMOV1+KPU(1)-1),DWORK(KMOV2+KPU(1)-1),
     *  DWORK(KMOV3+KPU(1)-1),DWORK(KMOV1+KPU(2)-1),
     *  DWORK(KMOV2+KPU(2)-1),DWORK(KMOV3+KPU(2)-1)
       WRITE(41,2000) REAL(TIMENS),REAL(DWORK(KMOV1+KPU(1)-1))
       WRITE(42,2000) REAL(TIMENS),REAL(DWORK(KMOV2+KPU(1)-1))
       WRITE(43,2000) REAL(TIMENS),REAL(DWORK(KMOV3+KPU(1)-1))
       WRITE(44,2000) REAL(TIMENS),REAL(DWORK(KMOV1+KPU(2)-1))
       WRITE(45,2000) REAL(TIMENS),REAL(DWORK(KMOV2+KPU(2)-1))
       WRITE(46,2000) REAL(TIMENS),REAL(DWORK(KMOV3+KPU(2)-1))
      ENDIF
C
C=======================================================================
C *** write PRESSURE 
C=======================================================================
C
      KPL  =L(LD1)
      KAUXM=L(LD2)
      LVOL=KLVOL(NLEV)
C
      CALL  INTPVD(DWORK(KP),DWORK(KPL),DWORK(KAUXM),
     *             VWORK(L(KLVOL(NLEV))),KWORK(L(LVERT)))
      IF (IER.NE.0) GOTO 99999
C
      IF ((ITYP.EQ.1).AND.(IFPSAV.GT.0)
     *               .AND.((TIMENS-DTFILO).GE.DTFILM)) THEN
       CFILE='#film/#DP      '
       IF ((ITFILM.GE.0).AND.(ITFILM.LT.10)) 
     *      WRITE(CFILE(10:10),'(I1.1)') ITFILM
       IF ((ITFILM.GE.10).AND.(ITFILM.LT.100)) 
     *      WRITE(CFILE(10:11),'(I2.2)') ITFILM
       IF ((ITFILM.GE.100).AND.(ITFILM.LT.1000)) 
     *      WRITE(CFILE(10:12),'(I3.3)') ITFILM
       IF ((ITFILM.GE.1000).AND.(ITFILM.LT.10000)) 
     *      WRITE(CFILE(10:13),'(I4.4)') ITFILM
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
       WRITE(47,2000) REAL(TIMENS),REAL(DWORK(KPL+KPP(1)-1))
       WRITE(48,2000) REAL(TIMENS),REAL(DWORK(KPL+KPP(2)-1))
       WRITE(49,2000) REAL(TIMENS),REAL(DWORK(KPL+KPP(3)-1))
       WRITE(50,2000) REAL(TIMENS),REAL(DWORK(KPL+KPP(4)-1))
      ENDIF
C
      CALL  BDPRES(DWORK(KPL),KWORK(L(LVERT)),KWORK(L(LAREA)),
     *             KWORK(L(KLABD(NLEV))),KWORK(L(KELBD(NLEV))),
     *             DWORK(L(LCORVG)),KNABD(NLEV),P5,P6)
C
      IF (ITYP.EQ.1) THEN
       IF (MSHOW.GE.2) 
     *  WRITE(MTERM,1003) P5,P6
       IF (MSHOW.GE.1) 
     *  WRITE(MFILE,1003) P5,P6
       WRITE(51,2000) REAL(TIMENS),REAL(P5)
       WRITE(52,2000) REAL(TIMENS),REAL(P6)
      ENDIF
C
C=======================================================================
C *** write lift and drag 
C=======================================================================
C
       IF ((IELT.EQ.0).OR.(IELT.EQ.2)) 
     *  CALL  BDFORC(DWORK(KU1),DWORK(KU2),DWORK(KU3),DWORK(KPL),
     *               KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *               KWORK(L(KLABD(NLEV))),KWORK(L(KELBD(NLEV))),
     *               DWORK(L(LCORVG)),EM31,DFW,DAW)
C
       IF ((IELT.EQ.1).OR.(IELT.EQ.3)) 
     *  CALL  BDFORC(DWORK(KU1),DWORK(KU2),DWORK(KU3),DWORK(KPL),
     *               KWORK(L(LVERT)),KWORK(L(LAREA)),KWORK(L(LEDGE)),
     *               KWORK(L(KLABD(NLEV))),KWORK(L(KELBD(NLEV))),
     *               DWORK(L(LCORVG)),EM30,DFW,DAW)
C
      IF (ITYP.EQ.1) THEN
       IF (MSHOW.GE.2) 
     *  WRITE(MTERM,1004) DFW,DAW
       IF (MSHOW.GE.1) 
     *  WRITE(MFILE,1004) DFW,DAW
       WRITE(53,2000) REAL(TIMENS),REAL(DFW)
       WRITE(54,2000) REAL(TIMENS),REAL(DAW)
      ENDIF
C
C=======================================================================
C
      IF ((ITYP.EQ.1).AND.((TIMENS-DTFILO).GE.DTFILM)) THEN
       DTFILO=TIMENS
      ENDIF
C
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
     *    ((ITYP.EQ.1).AND.(IAVS.GT.0)
     *                .AND.((TIMENS-DTAVSO).GE.DTAVS))) THEN
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
       CALL XAVS3D(39,CFILE,KNEL(IAVS),KNVT(IAVS),
     *             KWORK(L(KLVERT(IAVS))),DWORK(L(LCORVG)),
     *             VWORK(KAVS1),VWORK(KAVS2),VWORK(KAVS3),VWORK(KAVSP))
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
       CALL XGMV3D(39,CFILE,KNEL(IGMV),KNVT(IGMV),
     *             KWORK(L(KLVERT(IGMV))),DWORK(L(LCORVG)),
     *             VWORK(KAVS1),VWORK(KAVS2),VWORK(KAVS3),VWORK(KAVSP),
     *             TIMENS)
       DTGMVO=TIMENS
       ENDIF
C
       CALL ZDISP(0,LAVSP,'AVSP  ')
       CALL ZDISP(0,LAVSU,'AVSU  ')
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
1001  FORMAT ('P(VELO) ',6(D12.5))
1002  FORMAT ('P(PRES) ',4(D12.5))
1003  FORMAT ('I(PRES) ',2(D12.5))
1004  FORMAT ('I(FORCE)',2(D12.5))
2000  FORMAT(2E12.5,120(' '))
C
C
C
99999 END
