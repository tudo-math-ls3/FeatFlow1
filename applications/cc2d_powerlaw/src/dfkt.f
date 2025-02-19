************************************************************************
      SUBROUTINE XMADF1(KM1,KST1,KA1,KCOLA,KLDA,KF1,KF2,KU1,KU2,
     *                  NA,NU,THSTEP,ISTAT)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)      
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CFILST*12,CFILM*12,CFILE*12
      DIMENSION CFILST(NNLEV),CFILM(NNLEV)
      DATA CFILST/'#ns/ST1     ','#ns/ST2     ','#ns/ST3     ',
     *            '#ns/ST4     ','#ns/ST5     ','#ns/ST6     ',
     *            '#ns/ST7     ','#ns/ST8     ','#ns/ST9     '/
      DATA CFILM /'#ns/MA1     ','#ns/MA2     ','#ns/MA3     ',
     *            '#ns/MA4     ','#ns/MA5     ','#ns/MA6     ',
     *            '#ns/MA7     ','#ns/MA8     ','#ns/MA9     '/
C
      DIMENSION VWORK(1),KWORK(1)
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
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
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE
C-----------------------------------------------------------------------
C
      IF (IPRECA.EQ.0) THEN
C
       IF (IMASS.EQ.1) THEN
        IF (THSTEP.NE.0D0) THEN
         CALL LLC2(VWORK(KST1)  ,VWORK(KA1),NA,THSTEP,0D0)
         CALL LLC2(VWORK(KM1),VWORK(KA1),NA,1D0,1D0)
        ELSE
         CALL LCP2(VWORK(KM1),VWORK(KA1),NA)
        ENDIF
       ELSE
        IF (THSTEP.NE.0D0) THEN
         CALL LLC2(VWORK(KST1),VWORK(KA1),NA,THSTEP,0D0)
         DO 100 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
100      VWORK(KA1+INUA)=VWORK(KA1+INUA)+VWORK(KM1+INU-1)
        ELSE
         CALL LCL2(VWORK(KA1),NA)
         DO 110 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
110      VWORK(KA1+INUA)=VWORK(KM1+INU-1)
        ENDIF
       ENDIF
C
       CALL  LAX37 (VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *              DWORK(KU1),DWORK(KF1),1D0,1D0)
       CALL  LAX37 (VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *              DWORK(KU2),DWORK(KF2),1D0,1D0)
C
      ENDIF
C
C=======================================================================
C
      IF (IPRECA.EQ.1) THEN
C
       IF (IMASS.EQ.0) THEN
        IF (THSTEP.NE.0D0) THEN
         CALL  LAX17 (DWORK(KST1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                DWORK(KU1),DWORK(KF1),THSTEP,1D0)
         CALL  LAX17 (DWORK(KST1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                DWORK(KU2),DWORK(KF2),THSTEP,1D0)
C
         DO 200 INU=1,NU
         DWORK(KF1+INU-1)= DWORK(KF1+INU-1)
     *                    +DWORK(KM1+INU-1)*DWORK(KU1+INU-1)
         DWORK(KF2+INU-1)= DWORK(KF2+INU-1)
     *                    +DWORK(KM1+INU-1)*DWORK(KU2+INU-1)
200      CONTINUE
C
         DO 210 INA=1,NA
210      VWORK(KA1+INA-1)=REAL(THSTEP*DWORK(KST1+INA-1))
C
         DO 220 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
220      VWORK(KA1+INUA)=VWORK(KA1+INUA)+REAL(DWORK(KM1+INU-1))
        ELSE
         DO 230 INU=1,NU
         DWORK(KF1+INU-1)= DWORK(KF1+INU-1)
     *                    +DWORK(KM1+INU-1)*DWORK(KU1+INU-1)
         DWORK(KF2+INU-1)= DWORK(KF2+INU-1)
     *                    +DWORK(KM1+INU-1)*DWORK(KU2+INU-1)
230      CONTINUE
C
         CALL LCL2(VWORK(KA1),NA)
C
         DO 240 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
240      VWORK(KA1+INUA)=REAL(DWORK(KM1+INU-1))
        ENDIF
       ELSE
        IF (THSTEP.NE.0D0) THEN
         CALL  LAX17 (DWORK(KST1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                DWORK(KU1),DWORK(KF1),THSTEP,1D0)
         CALL  LAX17 (DWORK(KST1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                DWORK(KU2),DWORK(KF2),THSTEP,1D0)
         CALL  LAX17 (DWORK(KM1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                DWORK(KU1),DWORK(KF1),1D0,1D0)
         CALL  LAX17 (DWORK(KM1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                DWORK(KU2),DWORK(KF2),1D0,1D0)
C
         DO 250 INA=1,NA
250      VWORK(KA1+INA-1)=REAL(THSTEP*DWORK(KST1+INA-1)
     *                               +DWORK(KM1+INA-1))
        ELSE
         CALL  LAX17 (DWORK(KM1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                DWORK(KU1),DWORK(KF1),1D0,1D0)
         CALL  LAX17 (DWORK(KM1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                DWORK(KU2),DWORK(KF2),1D0,1D0)
C
         CALL LCL2(VWORK(KA1),NA)
C
         DO 260 INA=1,NA
260      VWORK(KA1+INA-1)=REAL(DWORK(KM1+INA-1))
        ENDIF
       ENDIF
C
      ENDIF
C
C=======================================================================
C
      IF (IPRECA.EQ.2) THEN
C
       IF (THSTEP.NE.0D0) THEN
        CALL  OF0 (59,CFILST(ILEV),0)
        CFILE='STMAT '
        CALL  ORA2 (VWORK(KA1),CFILE,59,0)
        REWIND(59)
        CLOSE (59)
        IF (IER.NE.0) RETURN
        CALL LLC2(VWORK(KA1),VWORK(KA1),NA,THSTEP,0D0)
       ELSE
        CALL LCL2(VWORK(KA1),NA)
       ENDIF
C
       IF (IMASS.EQ.1) THEN
        CALL  OF0 (59,CFILM(ILEV),0)
        CFILE='MASMAT'
        CALL  ORALC2 (VWORK(KA1),1D0,CFILE,59,0)
        REWIND(59)
        CLOSE (59)
        IF (IER.NE.0) RETURN
       ELSE
        IF (THSTEP.NE.0D0) THEN
         DO 300 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
300      VWORK(KA1+INUA)=VWORK(KA1+INUA)+VWORK(KM1+INU-1)
        ELSE
         DO 310 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
310      VWORK(KA1+INUA)=VWORK(KM1+INU-1)
        ENDIF
       ENDIF
C
       CALL  LAX37 (VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *              DWORK(KU1),DWORK(KF1),1D0,1D0)
       CALL  LAX37 (VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *              DWORK(KU2),DWORK(KF2),1D0,1D0)
C
      ENDIF
C
C=======================================================================
C
      IF (IPRECA.EQ.3) THEN
C
       CALL ZCTYPE(1,KLA(ILEV),'VA    ')
       IF (IER.NE.0) GOTO 99998
       KA1=L(KLA(ILEV))
C
       IF (THSTEP.NE.0D0) THEN
        CALL  OF0 (59,CFILST(ILEV),0)
        CFILE='STMAT '
        CALL  ORA1 (DWORK(KA1),CFILE,59,0)
        REWIND(59)
        CLOSE (59)
        IF (IER.NE.0) RETURN
        CALL LLC1(DWORK(KA1),DWORK(KA1),NA,THSTEP,0D0)
       ELSE
        CALL LCL1(DWORK(KA1),NA)
       ENDIF
C
       IF (IMASS.EQ.1) THEN
        CALL  OF0 (59,CFILM(ILEV),0)
        CFILE='MASMAT'
        CALL  ORALC1 (DWORK(KA1),1D0,CFILE,59,0)
        REWIND(59)
        CLOSE (59)
        IF (IER.NE.0) RETURN
       ELSE
        IF (THSTEP.NE.0D0) THEN
         DO 400 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
400      DWORK(KA1+INUA)=DWORK(KA1+INUA)+DWORK(KM1+INU-1)
        ELSE
         DO 410 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
410      DWORK(KA1+INUA)=DWORK(KM1+INU-1)
        ENDIF
       ENDIF
C
       CALL  LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *              DWORK(KU1),DWORK(KF1),1D0,1D0)
       CALL  LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *              DWORK(KU2),DWORK(KF2),1D0,1D0)
C
       CALL ZCTYPE(2,KLA(ILEV),'VA    ')
       IF (IER.NE.0) GOTO 99998
       KA1=L(KLA(ILEV))
C
      ENDIF
C
C=======================================================================
C
      IF (IPRECA.EQ.4) THEN
C
       IF (IMASS.EQ.0) THEN

        IF (ISTAT.NE.0) THEN
        DO 510 INU=1,NU
        DWORK(KF1+INU-1)= DWORK(KF1+INU-1)
     *                   +DWORK(KM1+INU-1)*DWORK(KU1+INU-1)
        DWORK(KF2+INU-1)= DWORK(KF2+INU-1)
     *                   +DWORK(KM1+INU-1)*DWORK(KU2+INU-1)
510     CONTINUE
        ENDIF
C
        CALL LCL2(VWORK(KA1),5*NA) !CHANGEDEV !
c        CALL LCL2(VWORK(KA1),4*NA) !CHANGEDEV !
C
        IF (ISTAT.NE.0) THEN
        DO 520 INU=1,NU
        INUA=KWORK(KLDA+INU-1)-1
        VWORK(KA1+INUA)=REAL(DWORK(KM1+INU-1))
520     VWORK(KA1+3*NA+INUA)=REAL(DWORK(KM1+INU-1)) !CHANGEDEV !
      ENDIF
      ELSE
        CALL LCL2(VWORK(KA1),5*NA) !CHANGEDEV !
c        CALL LCL2(VWORK(KA1),4*NA) !CHANGEDEV !
       ENDIF
C
      ENDIF
C
C
C
99998 END
C
C
C
************************************************************************
      SUBROUTINE XMADF2(KM1,KST1,KA1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                  KD1,KD2,KDP,KU1,KU2,KP,NA,NU,NP,KMBD,NMBD,INEUM,
     *                  THSTEP,ISTAT)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)      
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CFILST*12,CFILM*12,CFILE*12
      DIMENSION CFILST(NNLEV),CFILM(NNLEV)
      DATA CFILST/'#ns/ST1     ','#ns/ST2     ','#ns/ST3     ',
     *            '#ns/ST4     ','#ns/ST5     ','#ns/ST6     ',
     *            '#ns/ST7     ','#ns/ST8     ','#ns/ST9     '/
      DATA CFILM /'#ns/MA1     ','#ns/MA2     ','#ns/MA3     ',
     *            '#ns/MA4     ','#ns/MA5     ','#ns/MA6     ',
     *            '#ns/MA7     ','#ns/MA8     ','#ns/MA9     '/
C
      DIMENSION KMBD(*)
      DIMENSION VWORK(1),KWORK(1)
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
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
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE
C-----------------------------------------------------------------------
C
      IF (IPRECA.EQ.0) THEN
C

       IF (ISTAT.EQ.1) THEN
        IF (IMASS.EQ.1) THEN
         IF (THSTEP.NE.0D0) THEN
          CALL LLC2(VWORK(KST1),VWORK(KA1),NA,THSTEP,0D0)
          CALL LLC2(VWORK(KM1 ),VWORK(KA1),NA,1D0,1D0)
         ELSE
          CALL LCP2(VWORK(KM1),VWORK(KA1),NA)
         ENDIF
        ELSE
         IF (THSTEP.NE.0D0) THEN
          CALL LLC2(VWORK(KST1),VWORK(KA1),NA,THSTEP,0D0)
          DO 100 INU=1,NU
          INUA=KWORK(KLDA+INU-1)-1
100       VWORK(KA1+INUA)=VWORK(KA1+INUA)+VWORK(KM1+INU-1)
         ELSE
          CALL LCL2(VWORK(KA1),NA)
          DO 110 INU=1,NU
          INUA=KWORK(KLDA+INU-1)-1
110       VWORK(KA1+INUA)=VWORK(KM1+INU-1)
         ENDIF
        ENDIF
       ELSE
        CALL LLC2(VWORK(KST1),VWORK(KA1),NA,THSTEP,0D0)
       ENDIF

C
       IF (INLMAX.GT.1) THEN
        CALL  LAX37 (VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *               DWORK(KU1),DWORK(KD1),-1D0,1D0)
        CALL  LAX37 (VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *               DWORK(KU2),DWORK(KD2),-1D0,1D0)
       ENDIF
C

      ENDIF
C
C=======================================================================
C
      IF (IPRECA.EQ.1) THEN
C
      IF (ISTAT.EQ.1) THEN
       IF (IMASS.EQ.0) THEN
        IF (THSTEP.NE.0D0) THEN
         IF (INLMAX.GT.1) THEN
          CALL  LAX17 (DWORK(KST1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DWORK(KU1),DWORK(KD1),-THSTEP,1D0)
          CALL  LAX17 (DWORK(KST1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DWORK(KU2),DWORK(KD2),-THSTEP,1D0)
C
          DO 200 INU=1,NU
          DWORK(KD1+INU-1)= DWORK(KD1+INU-1)
     *                     -DWORK(KM1+INU-1)*DWORK(KU1+INU-1)
          DWORK(KD2+INU-1)= DWORK(KD2+INU-1)
     *                     -DWORK(KM1+INU-1)*DWORK(KU2+INU-1)
200       CONTINUE
         ENDIF
C
         DO 210 INA=1,NA
210      VWORK(KA1+INA-1)=REAL(THSTEP*DWORK(KST1+INA-1))
C
         DO 220 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
220      VWORK(KA1+INUA)=VWORK(KA1+INUA)+REAL(DWORK(KM1+INU-1))
        ELSE
         IF (INLMAX.GT.1) THEN
          DO 230 INU=1,NU
          DWORK(KD1+INU-1)= DWORK(KD1+INU-1)
     *                     -DWORK(KM1+INU-1)*DWORK(KU1+INU-1)
          DWORK(KD2+INU-1)= DWORK(KD2+INU-1)
     *                     -DWORK(KM1+INU-1)*DWORK(KU2+INU-1)
230       CONTINUE
         ENDIF
C
         CALL LCL2(VWORK(KA1),NA)
C
         DO 240 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
240      VWORK(KA1+INUA)=REAL(DWORK(KM1+INU-1))
        ENDIF
       ELSE
        IF (THSTEP.NE.0D0) THEN
         IF (INLMAX.GT.1) THEN
          CALL  LAX17 (DWORK(KST1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DWORK(KU1),DWORK(KD1),-THSTEP,1D0)
          CALL  LAX17 (DWORK(KST1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DWORK(KU2),DWORK(KD2),-THSTEP,1D0)
          CALL  LAX17 (DWORK(KM1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DWORK(KU1),DWORK(KD1),-1D0,1D0)
          CALL  LAX17 (DWORK(KM1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DWORK(KU2),DWORK(KD2),-1D0,1D0)
         ENDIF
C
         DO 250 INA=1,NA
250      VWORK(KA1+INA-1)=REAL(THSTEP*DWORK(KST1+INA-1)
     *                               +DWORK(KM1+INA-1))
        ELSE
         IF (INLMAX.GT.1) THEN
          CALL  LAX17 (DWORK(KM1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DWORK(KU1),DWORK(KD1),-1D0,1D0)
          CALL  LAX17 (DWORK(KM1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DWORK(KU2),DWORK(KD2),-1D0,1D0)
         ENDIF
C
         CALL LCL2(VWORK(KA1),NA)
C
         DO 260 INA=1,NA
260      VWORK(KA1+INA-1)=REAL(DWORK(KM1+INA-1))
        ENDIF
       ENDIF
      ELSE
       IF (INLMAX.GT.1) THEN
        CALL  LAX17 (DWORK(KST1),KWORK(KCOLA),KWORK(KLDA),NU,
     *               DWORK(KU1),DWORK(KD1),-THSTEP,1D0)
        CALL  LAX17 (DWORK(KST1),KWORK(KCOLA),KWORK(KLDA),NU,
     *               DWORK(KU2),DWORK(KD2),-THSTEP,1D0)
       ENDIF
       DO 270 INA=1,NA
270    VWORK(KA1+INA-1)=REAL(THSTEP*DWORK(KST1+INA-1))
      ENDIF
C
      ENDIF
C
C=======================================================================
C
      IF (IPRECA.EQ.2) THEN
C
      IF (ISTAT.EQ.1) THEN
       IF (THSTEP.NE.0D0) THEN
        CALL  OF0 (59,CFILST(ILEV),0)
        CFILE='STMAT '
        CALL  ORA2 (VWORK(KA1),CFILE,59,0)
        REWIND(59)
        CLOSE (59)
        IF (IER.NE.0) RETURN
        CALL LLC2(VWORK(KA1),VWORK(KA1),NA,THSTEP,0D0)
       ELSE
        CALL LCL2(VWORK(KA1),NA)
       ENDIF
C
       IF (IMASS.EQ.1) THEN
        CALL  OF0 (59,CFILM(ILEV),0)
        CFILE='MASMAT'
        CALL  ORALC2 (VWORK(KA1),1D0,CFILE,59,0)
        REWIND(59)
        CLOSE (59)
        IF (IER.NE.0) RETURN
       ELSE
        IF (THSTEP.NE.0D0) THEN
         DO 300 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
300      VWORK(KA1+INUA)=VWORK(KA1+INUA)+VWORK(KM1+INU-1)
        ELSE
         DO 310 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
310      VWORK(KA1+INUA)=VWORK(KM1+INU-1)
        ENDIF
       ENDIF
      ELSE
       CALL  OF0 (59,CFILST(ILEV),0)
       CFILE='STMAT '
       CALL  ORA2 (VWORK(KA1),CFILE,59,0)
       REWIND(59)
       CLOSE (59)
       IF (IER.NE.0) RETURN
       CALL LLC2(VWORK(KA1),VWORK(KA1),NA,THSTEP,0D0)
      ENDIF
C
      IF (INLMAX.GT.1) THEN
       CALL  LAX37 (VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *              DWORK(KU1),DWORK(KD1),-1D0,1D0)
       CALL  LAX37 (VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *              DWORK(KU2),DWORK(KD2),-1D0,1D0)
      ENDIF
C
      ENDIF
C
C=======================================================================
C
      IF (IPRECA.EQ.3) THEN
C
       CALL ZCTYPE(1,KLA(ILEV),'VA    ')
       IF (IER.NE.0) GOTO 99998
       KA1=L(KLA(ILEV))
C
       IF (THSTEP.NE.0D0) THEN
        CALL  OF0 (59,CFILST(ILEV),0)
        CFILE='STMAT '
        CALL  ORA1 (DWORK(KA1),CFILE,59,0)
        REWIND(59)
        CLOSE (59)
        IF (IER.NE.0) RETURN
        CALL LLC1(DWORK(KA1),DWORK(KA1),NA,THSTEP,0D0)
       ELSE
        CALL LCL1(DWORK(KA1),NA)
       ENDIF
C
       IF (ISTAT.EQ.1) THEN
        IF (IMASS.EQ.1) THEN
         CALL  OF0 (59,CFILM(ILEV),0)
         CFILE='MASMAT'
         CALL  ORALC1 (DWORK(KA1),1D0,CFILE,59,0)
         REWIND(59)
         CLOSE (59)
         IF (IER.NE.0) RETURN
        ELSE
         IF (THSTEP.NE.0D0) THEN
          DO 400 INU=1,NU
          INUA=KWORK(KLDA+INU-1)-1
400       DWORK(KA1+INUA)=DWORK(KA1+INUA)+DWORK(KM1+INU-1)
         ELSE
          DO 410 INU=1,NU
          INUA=KWORK(KLDA+INU-1)-1
410       DWORK(KA1+INUA)=DWORK(KM1+INU-1)
         ENDIF
        ENDIF
       ENDIF
C
       IF (INLMAX.GT.1) THEN
        CALL  LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *               DWORK(KU1),DWORK(KD1),-1D0,1D0)
        CALL  LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *               DWORK(KU2),DWORK(KD2),-1D0,1D0)
       ENDIF
C
       CALL ZCTYPE(2,KLA(ILEV),'VA    ')
       IF (IER.NE.0) GOTO 99998
       KA1=L(KLA(ILEV))
C
      ENDIF
C
C=======================================================================
C
      IF (IPRECA.EQ.4) THEN
C
      IF (ISTAT.EQ.1) THEN
       IF (IMASS.EQ.0) THEN
        IF (INLMAX.GT.1) THEN
         DO 510 INU=1,NU
         DWORK(KD1+INU-1)= DWORK(KD1+INU-1)
     *                    -DWORK(KM1+INU-1)*DWORK(KU1+INU-1)
         DWORK(KD2+INU-1)= DWORK(KD2+INU-1)
     *                    -DWORK(KM1+INU-1)*DWORK(KU2+INU-1)
510      CONTINUE
        ENDIF
C
        CALL LCL2(VWORK(KA1),5*NA) !CHANGEDEV !
c        CALL LCL2(VWORK(KA1),4*NA) !CHANGEDEV !
C
        DO 520 INU=1,NU
        INUA=KWORK(KLDA+INU-1)-1
        VWORK(KA1+INUA)=REAL(DWORK(KM1+INU-1))
        VWORK(KA1+3*NA+INUA)=REAL(DWORK(KM1+INU-1))
 520    continue
       ELSE         !imass
        CALL LCL2(VWORK(KA1),5*NA)
c        CALL LCL2(VWORK(KA1),4*NA)
       ENDIF        !imass
      ELSE      !istat
       CALL LCL2(VWORK(KA1),5*NA)
c       CALL LCL2(VWORK(KA1),4*NA)
      ENDIF     !istat
C
      ENDIF !ipreca.eq.4
C
C=======================================================================
C
      IF (IPRECB.EQ.0) THEN
       IF (INEUM.EQ.1) CALL  BDRDEF(DWORK(KD1),KMBD,NMBD,0.5D0)
       CALL  LAX39 (VWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,
     *              DWORK(KP),DWORK(KD1),-1D0,1D0)
       IF (INEUM.EQ.1) CALL  BDRDEF(DWORK(KD1),KMBD,NMBD,2.0D0)
       IF (INEUM.EQ.1) CALL  BDRDEF(DWORK(KD2),KMBD,NMBD,0.5D0)
       CALL  LAX39 (VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,
     *              DWORK(KP),DWORK(KD2),-1D0,1D0)
       IF (INEUM.EQ.1) CALL  BDRDEF(DWORK(KD2),KMBD,NMBD,2.0D0)
      ENDIF
C
      IF (IPRECB.EQ.1) THEN
       IF (INEUM.EQ.1) CALL  BDRDEF(DWORK(KD1),KMBD,NMBD,0.5D0)
       CALL  LAX19 (DWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,
     *              DWORK(KP),DWORK(KD1),-1D0,1D0)
       IF (INEUM.EQ.1) CALL  BDRDEF(DWORK(KD1),KMBD,NMBD,2.0D0)
       IF (INEUM.EQ.1) CALL  BDRDEF(DWORK(KD2),KMBD,NMBD,0.5D0)
       CALL  LAX19 (DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,
     *              DWORK(KP),DWORK(KD2),-1D0,1D0)
       IF (INEUM.EQ.1) CALL  BDRDEF(DWORK(KD2),KMBD,NMBD,2.0D0)
      ENDIF
C
      IF (IPRECB.EQ.2) THEN
       CALL BMUL1 (KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *             DWORK(L(LCORVG)),DWORK(KP),DWORK(KD1),DWORK(KD2),
     *             NEL,NVT,NMT,-1D0,1D0)
      ENDIF
C
C=======================================================================
C
      IF (INLMAX.GT.1) CALL  BDRY0 (DWORK(KD1),DWORK(KD2),KMBD,NMBD)
C
C=======================================================================
C
      IF (IPRECB.EQ.0) THEN
       CALL  LTX39 (VWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU1),DWORK(KDP),-1D0,1D0)
       CALL  LTX39 (VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU2),DWORK(KDP),-1D0,1D0)
      ENDIF
C
      IF (IPRECB.EQ.1) THEN
       CALL  LTX19 (DWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU1),DWORK(KDP),-1D0,1D0)
       CALL  LTX19 (DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU2),DWORK(KDP),-1D0,1D0)
      ENDIF
C
      IF (IPRECB.EQ.2) THEN
       WRITE(6,*) 'BTMUL ???'
       CALL BTMUL1(KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *             DWORK(L(LCORVG)),DWORK(KU1),DWORK(KU2),DWORK(KDP),
     *             NEL,NVT,NMT,-1D0)
      ENDIF
C
C
C
99998 END
C
C
C
************************************************************************
      SUBROUTINE XMADF3(KM1,KST1,KA1,KCOLA,KLDA,NA,NU,THSTEP,ISTAT)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)      
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CFILST*12,CFILM*12,CFILE*12
      DIMENSION CFILST(NNLEV),CFILM(NNLEV)
      DATA CFILST/'#ns/ST1     ','#ns/ST2     ','#ns/ST3     ',
     *            '#ns/ST4     ','#ns/ST5     ','#ns/ST6     ',
     *            '#ns/ST7     ','#ns/ST8     ','#ns/ST9     '/
      DATA CFILM /'#ns/MA1     ','#ns/MA2     ','#ns/MA3     ',
     *            '#ns/MA4     ','#ns/MA5     ','#ns/MA6     ',
     *            '#ns/MA7     ','#ns/MA8     ','#ns/MA9     '/
C
      DIMENSION VWORK(1),KWORK(1)
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
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
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE
C-----------------------------------------------------------------------
C
      IF (IPRECA.EQ.0) THEN
C
      IF (ISTAT.EQ.1) THEN
       IF (IMASS.EQ.1) THEN
        IF (THSTEP.NE.0D0) THEN
         CALL LLC2(VWORK(KST1)  ,VWORK(KA1),NA,THSTEP,0D0)
         CALL LLC2(VWORK(KM1),VWORK(KA1),NA,1D0,1D0)
        ELSE
         CALL LCP2(VWORK(KM1),VWORK(KA1),NA)
        ENDIF
       ELSE
        IF (THSTEP.NE.0D0) THEN
         CALL LLC2(VWORK(KST1),VWORK(KA1),NA,THSTEP,0D0)
         DO 100 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
100      VWORK(KA1+INUA)=VWORK(KA1+INUA)+VWORK(KM1+INU-1)
        ELSE
         CALL LCL2(VWORK(KA1),NA)
         DO 110 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
110      VWORK(KA1+INUA)=VWORK(KM1+INU-1)
        ENDIF
       ENDIF
      ELSE
       CALL LLC2(VWORK(KST1),VWORK(KA1),NA,THSTEP,0D0)
      ENDIF
C
      ENDIF
C
C=======================================================================
C
      IF (IPRECA.EQ.1) THEN
C
       IF (IMASS.EQ.0) THEN
        IF (THSTEP.NE.0D0) THEN
         DO 210 INA=1,NA
210      VWORK(KA1+INA-1)=REAL(THSTEP*DWORK(KST1+INA-1))
C
         DO 220 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
220      VWORK(KA1+INUA)=VWORK(KA1+INUA)+REAL(DWORK(KM1+INU-1))
        ELSE
         CALL LCL2(VWORK(KA1),NA)
         DO 240 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
240      VWORK(KA1+INUA)=REAL(DWORK(KM1+INU-1))
        ENDIF
       ELSE
        IF (THSTEP.NE.0D0) THEN
         DO 250 INA=1,NA
250      VWORK(KA1+INA-1)=REAL(THSTEP*DWORK(KST1+INA-1)
     *                               +DWORK(KM1+INA-1))
        ELSE
         CALL LCL2(VWORK(KA1),NA)
         DO 260 INA=1,NA
260      VWORK(KA1+INA-1)=REAL(DWORK(KM1+INA-1))
        ENDIF
       ENDIF
C
      ENDIF
C
C=======================================================================
C
      IF (IPRECA.EQ.2) THEN
C
       IF (THSTEP.NE.0D0) THEN
        CALL  OF0 (59,CFILST(ILEV),0)
        CFILE='STMAT '
        CALL  ORA2 (VWORK(KA1),CFILE,59,0)
        REWIND(59)
        CLOSE (59)
        IF (IER.NE.0) RETURN
        CALL LLC2(VWORK(KA1),VWORK(KA1),NA,THSTEP,0D0)
       ELSE
        CALL LCL2(VWORK(KA1),NA)
       ENDIF
C
       IF (IMASS.EQ.1) THEN
        CALL  OF0 (59,CFILM(ILEV),0)
        CFILE='MASMAT'
        CALL  ORALC2 (VWORK(KA1),1D0,CFILE,59,0)
        REWIND(59)
        CLOSE (59)
        IF (IER.NE.0) RETURN
       ELSE
        IF (THSTEP.NE.0D0) THEN
         DO 300 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
300      VWORK(KA1+INUA)=VWORK(KA1+INUA)+REAL(DWORK(KM1+INU-1))
        ELSE
         DO 310 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
310      VWORK(KA1+INUA)=REAL(DWORK(KM1+INU-1))
        ENDIF
       ENDIF
C
      ENDIF
C
C=======================================================================
C
      IF (IPRECA.EQ.3) THEN
C
       CALL ZCTYPE(1,KLA(ILEV),'VA    ')
       IF (IER.NE.0) GOTO 99998
       KA1=L(KLA(ILEV))
C
       IF (THSTEP.NE.0D0) THEN
        CALL  OF0 (59,CFILST(ILEV),0)
        CFILE='STMAT '
        CALL  ORA1 (DWORK(KA1),CFILE,59,0)
        REWIND(59)
        CLOSE (59)
        IF (IER.NE.0) RETURN
        CALL LLC1(DWORK(KA1),DWORK(KA1),NA,THSTEP,0D0)
       ELSE
        CALL LCL1(DWORK(KA1),NA)
       ENDIF
C
       IF (IMASS.EQ.1) THEN
        CALL  OF0 (59,CFILM(ILEV),0)
        CFILE='MASMAT'
        CALL  ORALC1 (DWORK(KA1),1D0,CFILE,59,0)
        REWIND(59)
        CLOSE (59)
        IF (IER.NE.0) RETURN
       ELSE
        IF (THSTEP.NE.0D0) THEN
         DO 400 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
400      DWORK(KA1+INUA)=DWORK(KA1+INUA)+DWORK(KM1+INU-1)
        ELSE
         DO 410 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
410      DWORK(KA1+INUA)=DWORK(KM1+INU-1)
        ENDIF
       ENDIF
C
       CALL ZCTYPE(2,KLA(ILEV),'VA    ')
       IF (IER.NE.0) GOTO 99998
       KA1=L(KLA(ILEV))
C
      ENDIF
C
C=======================================================================
C
      IF (IPRECA.EQ.4) THEN
C
       IF (IMASS.EQ.0) THEN
        CALL LCL2(VWORK(KA1),5*NA) !CHANGEDEV !
c        CALL LCL2(VWORK(KA1),4*NA) !CHANGEDEV !
        IF (ISTAT.NE.0) THEN
          DO 510 INU=1,NU
          INUA=KWORK(KLDA+INU-1)-1
          VWORK(KA1+INUA)=REAL(DWORK(KM1+INU-1))
          VWORK(KA1+3*NA+INUA)=REAL(DWORK(KM1+INU-1))
 510      continue
        ENDIF
       ELSE
        CALL LCL2(VWORK(KA1),5*NA)
c        CALL LCL2(VWORK(KA1),4*NA)
       ENDIF
C
      ENDIF
C
C=======================================================================
C
99998 END
C
C
C
************************************************************************
      SUBROUTINE RESDFK(U1,U2,P,D1,D2,DP,F1,F2,FP,NU,NP,RESU,RESDIV)
************************************************************************
*    Purpose:  Computes the norms RESU, RESDIV
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z)
      PARAMETER (NNLEV=9)
      DIMENSION U1(*),U2(*),P(*),D1(*),D2(*),DP(*),F1(*),F2(*),FP(*)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      INCLUDE 'block.inc'
      SAVE 
C
C-----------------------------------------------------------------------
C     Compute the relative l2-norms  RESU,RESDIV
C-----------------------------------------------------------------------
C
      CALL LL21 (F1,NU,RESF1)
      CALL LL21 (F2,NU,RESF2)
      RESF=MAX(RESF1,RESF2)
      IF (ABS(RESF).LT.1D-8) RESF=1D0
C
      CALL LL21 (D1,NU,RESU1)
      CALL LL21 (D2,NU,RESU2)
c$$$      RESU=SQRT(RESU1*RESU1+RESU2*RESU2)/RESF
       RESU=MAX(RESu1,RESu2)/SQRT(DBLE(NU))
C
C
      CALL LL21 (U1,NU,DNRMU1)
      CALL LL21 (U2,NU,DNRMU2)
      DNORMU=SQRT(DNRMU1*DNRMU1+DNRMU2*DNRMU2)
      IF (ABS(DNORMU).LT.1D-8) DNORMU=1D0
C
      CALL LL21 (DP,NP,RESDIV)
       RESDIV=RESDIV/SQRT(DBLE(NP))
c      RESDIV=RESDIV/DNORMU
C

      END
C
C
C
************************************************************************
      SUBROUTINE RESDFC(U1,U2,P,D1,D2,DP,RESU,RESDIV,VMASSU,VMASSP)
************************************************************************
C
C-----------------------------------------------------------------------
*    Purpose:  Computes the residuals
C
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'
C-----------------------------------------------------------------------
C
*-----------------------------------------------------------------------
      DIMENSION U1(*),U2(*),P(*),D1(*),D2(*),DP(*)
      DIMENSION VMASSU(*),VMASSP(*)
C
C
C
      IF (ABS(ICHDEF).EQ.1) THEN
       CALL LL21 (D1,NU,RES1)
       CALL LL21 (D2,NU,RES2)
       RESU=MAX(RES1,RES2)/SQRT(DBLE(NU))
C
       CALL LL21 (DP,NP,RESDIV)
       RESDIV=RESDIV/SQRT(DBLE(NP))
      ENDIF
C
C
      IF (ABS(ICHDEF).EQ.2) THEN
       CALL LLI1 (D1,NU,RES1,IND1)
       CALL LLI1 (D2,NU,RES2,IND2)
       RESU=MAX(RES1,RES2)
C
       CALL LLI1 (DP,NP,RESDIV,INDP)
      ENDIF
C
C
      IF (ABS(ICHDEF).EQ.3) THEN
c         write (*,*) ilev, nu, '\n',(DBLE(VMASSU(i)),
c     *    dble(vwork(l(klm(ilev))-1+i)),
c     *    D1(I),D2(I),'\n', i=1,25)
c         write (*,*) ilev, nu, '\n'
       RES1=0D0
       RES2=0D0
       DO 300 IU=1,NU
       RES1=RES1+D1(IU)**2/DBLE(VMASSU(IU))
       RES2=RES2+D2(IU)**2/DBLE(VMASSU(IU))
300    CONTINUE
       RES1=SQRT(RES1)
       RES2=SQRT(RES2)
       RESU=MAX(RES1,RES2)
C
       RESDIV=0D0
       DO 310 IP=1,NP
       RESDIV=RESDIV+DP(IP)**2/DBLE(VMASSP(IP))
310    CONTINUE
       RESDIV=SQRT(RESDIV)
      ENDIF
C
C
      IF (ABS(ICHDEF).EQ.4) THEN
       RES1=0D0
       RES2=0D0
       DO 400 IU=1,NU
       RES1=MAX(RES1,ABS(D1(IU)/DBLE(VMASSU(IU))))
       RES2=MAX(RES2,ABS(D2(IU)/DBLE(VMASSU(IU))))
400    CONTINUE
       RESU=MAX(RES1,RES2)
C
       RESDIV=0D0
       DO 410 IP=1,NP
       RESDIV=MAX(RESDIV,ABS(DP(IP)/DBLE(VMASSP(IP))))
410    CONTINUE
      ENDIF
C
C
C
      END
C
************************************************************************
      SUBROUTINE RESDU(DX,NX,ILEVU,RES)
************************************************************************
C
C-----------------------------------------------------------------------
C Purpose:   
C
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'
C-----------------------------------------------------------------------
C
      DIMENSION DX(*)
C
C
C
      IF (ABS(ICHDEF).EQ.1) THEN
       CALL LL21 (DX,NX,RES)
       RES=RES/SQRT(DBLE(NX))
      ENDIF
C
C
      IF (ABS(ICHDEF).EQ.2) THEN
       CALL LLI1 (DX,NX,RES,INDU)
      ENDIF
C
C
      IF (ABS(ICHDEF).EQ.3) THEN
       RES=0D0
       DO 300 IU=1,NX
       RES=RES+DX(IU)**2/DBLE(VWORK(KM1+IU-1))
300    CONTINUE
       RES=SQRT(RES)
      ENDIF
C
C
      IF (ABS(ICHDEF).EQ.4) THEN
       RES=0D0
       DO 400 IU=1,NX
       RES=MAX(RES,ABS(DX(IU)/DBLE(VWORK(KM1+IU-1))))
400    CONTINUE
      ENDIF
C
C
C
      END
c
************************************************************************
      SUBROUTINE RESDP(DX,NX,ILEVP,RES)
************************************************************************
C
C-----------------------------------------------------------------------
C Purpose:   
C
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'
C-----------------------------------------------------------------------
C
      DIMENSION DX(*)
C
C
C
      IF (ABS(ICHDEF).EQ.1) THEN
       CALL LL21 (DX,NX,RES)
       RES=RES/SQRT(DBLE(NX))
      ENDIF
C
C
      IF (ABS(ICHDEF).EQ.2) THEN
       CALL LLI1 (DX,NX,RES,INDP)
      ENDIF
C
C
      IF (ABS(ICHDEF).EQ.3) THEN
       RES=0D0
       DO 300 IP=1,NX
       RES=RES+DX(IP)**2/DBLE(VWORK(L(KLAREA(ILEVP))+IP-1))
300    CONTINUE
       RES=SQRT(RES)
      ENDIF
C
C
      IF (ABS(ICHDEF).EQ.4) THEN
       RES=0D0
       DO 400 IP=1,NX
       RES=MAX(RES,ABS(DX(IP)/DBLE(VWORK(L(KLAREA(ILEVP))+IP-1))))
400    CONTINUE
      ENDIF
C
C
C
      END
************************************************************************
      SUBROUTINE XDFKD(A1,A2)
************************************************************************
C
C-----------------------------------------------------------------------
C Purpose:   
C
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'
C-----------------------------------------------------------------------
C
C
C
ccc      IF (A2.EQ.0D0) WRITE(6,*) 'A2=0 !!!!!    XDFKD'
C
      IF (IPRECB.EQ.0) THEN
       CALL  LTX39 (VWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU1),DWORK(KFP),1D0/A1,A2)
       CALL  LTX39 (VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU2),DWORK(KFP),1D0/A1,1D0)
      ENDIF
C
C=======================================================================
C
      IF (IPRECB.EQ.1) THEN
       CALL  LTX19 (DWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU1),DWORK(KFP),1D0/A1,A2)
       CALL  LTX19 (DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU2),DWORK(KFP),1D0/A1,1D0)
      ENDIF
C
C=======================================================================
C
      IF (IPRECB.EQ.2) THEN
      IF (A2.EQ.0D0) WRITE(6,*) 'A2=0 !!!!!    ?????'
       CALL BTMUL1(KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *             DWORK(L(LCORVG)),DWORK(KU1),DWORK(KU2),DWORK(KFP),
     *             NEL,NVT,NMT,1D0/A1)
      ENDIF
C
C
C
      END
C
C
C





