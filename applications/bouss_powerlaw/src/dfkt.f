************************************************************************
      SUBROUTINE XMADF1(KM1,KMASS1,KST1,KA1,KCOLA,KLDA,KB1,KB2,
     *                  KCOLB,KLDB,KF1,KF2,KU1,KU2,KP,NA,NU,
     *                  THSTEP,TSTEPH,IPROJ)
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
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
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
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
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
         CALL LLC2(VWORK(KMASS1),VWORK(KA1),NA,1D0,1D0)
        ELSE
         CALL LCP2(VWORK(KMASS1),VWORK(KA1),NA)
        ENDIF
       ELSE
        IF (THSTEP.NE.0D0) THEN
         CALL LLC2(VWORK(KST1),VWORK(KA1),NA,THSTEP,0D0)
         DO 100 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
100      VWORK(KA1+INUA)=VWORK(KA1+INUA)+REAL(DWORK(KM1+INU-1))
        ELSE
         CALL LCL2(VWORK(KA1),NA)
         DO 110 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
110      VWORK(KA1+INUA)=REAL(DWORK(KM1+INU-1))
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
         CALL  LAX17 (DWORK(KMASS1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                DWORK(KU1),DWORK(KF1),1D0,1D0)
         CALL  LAX17 (DWORK(KMASS1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                DWORK(KU2),DWORK(KF2),1D0,1D0)
C
         DO 250 INA=1,NA
250      VWORK(KA1+INA-1)=REAL(THSTEP*DWORK(KST1+INA-1)
     *                               +DWORK(KMASS1+INA-1))
        ELSE
         CALL  LAX17 (DWORK(KMASS1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                DWORK(KU1),DWORK(KF1),1D0,1D0)
         CALL  LAX17 (DWORK(KMASS1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                DWORK(KU2),DWORK(KF2),1D0,1D0)
C
         CALL LCL2(VWORK(KA1),NA)
C
         DO 260 INA=1,NA
260      VWORK(KA1+INA-1)=REAL(DWORK(KMASS1+INA-1))
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
        DO 510 INU=1,NU
        DWORK(KF1+INU-1)= DWORK(KF1+INU-1)
     *                   +DWORK(KM1+INU-1)*DWORK(KU1+INU-1)
        DWORK(KF2+INU-1)= DWORK(KF2+INU-1)
     *                   +DWORK(KM1+INU-1)*DWORK(KU2+INU-1)
510     CONTINUE
C
        CALL LCL2(VWORK(KA1),4*NA)
C
        DO 520 INU=1,NU
           INUA=KWORK(KLDA+INU-1)-1
           VWORK(KA1+INUA)=REAL(DWORK(KM1+INU-1))
           VWORK(3*NA+KA1+INUA)=REAL(DWORK(KM1+INU-1))
520     continue
       ELSE
        CALL LCL2(VWORK(KA1),4*NA)
       ENDIF
C
      ENDIF
C
C=======================================================================
C
       IF (IPROJ.EQ.1) THEN
        IF (IPRECB.EQ.0) THEN
         CALL  LAX39 (VWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,
     *                DWORK(KP),DWORK(KF1),-TSTEPH,1D0)
         CALL  LAX39 (VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,
     *                DWORK(KP),DWORK(KF2),-TSTEPH,1D0)
        ENDIF
        IF (IPRECB.EQ.1) THEN
         CALL  LAX19 (DWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,
     *                DWORK(KP),DWORK(KF1),-TSTEPH,1D0)
         CALL  LAX19 (DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,
     *                DWORK(KP),DWORK(KF2),-TSTEPH,1D0)
        ENDIF
        IF (IPRECB.EQ.2) THEN
         CALL BMUL1 (KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *               DWORK(L(LCORVG)),DWORK(KP),DWORK(KF1),DWORK(KF2),
     *               NEL,NVT,NMT,-TSTEPH,1D0)
        ENDIF
       ENDIF
C
C
99998 END
C
C
C
************************************************************************
      SUBROUTINE XMADF2(KM1,KMASS1,KST1,KA1,KCOLA,KLDA,KD1,KD2,KU1,KU2,
     *                  NA,NU,KMBD,NMBD,THSTEP)
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
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
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
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
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
         CALL LLC2(VWORK(KMASS1),VWORK(KA1),NA,1D0,1D0)
        ELSE
         CALL LCP2(VWORK(KMASS1),VWORK(KA1),NA)
        ENDIF
       ELSE
        IF (THSTEP.NE.0D0) THEN
         CALL LLC2(VWORK(KST1),VWORK(KA1),NA,THSTEP,0D0)
         DO 100 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
100      VWORK(KA1+INUA)=VWORK(KA1+INUA)+REAL(DWORK(KM1+INU-1))
        ELSE
         CALL LCL2(VWORK(KA1),NA)
         DO 110 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
110      VWORK(KA1+INUA)=REAL(DWORK(KM1+INU-1))
        ENDIF
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
          CALL  LAX17 (DWORK(KMASS1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DWORK(KU1),DWORK(KD1),-1D0,1D0)
          CALL  LAX17 (DWORK(KMASS1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DWORK(KU2),DWORK(KD2),-1D0,1D0)
         ENDIF
C
         DO 250 INA=1,NA
250      VWORK(KA1+INA-1)=REAL(THSTEP*DWORK(KST1+INA-1)
     *                               +DWORK(KMASS1+INA-1))
        ELSE
         IF (INLMAX.GT.1) THEN
          CALL  LAX17 (DWORK(KMASS1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DWORK(KU1),DWORK(KD1),-1D0,1D0)
          CALL  LAX17 (DWORK(KMASS1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DWORK(KU2),DWORK(KD2),-1D0,1D0)
         ENDIF
C
         CALL LCL2(VWORK(KA1),NA)
C
         DO 260 INA=1,NA
260      VWORK(KA1+INA-1)=REAL(DWORK(KMASS1+INA-1))
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
        CALL LCL2(VWORK(KA1),4*NA)
C
        DO 520 INU=1,NU
        INUA=KWORK(KLDA+INU-1)-1
        VWORK(KA1+INUA)=REAL(DWORK(KM1+INU-1))
        VWORK(3*NA+KA1+INUA)=REAL(DWORK(KM1+INU-1))
520     continue
       ELSE
        CALL LCL2(VWORK(KA1),4*NA)
       ENDIF
C
      ENDIF
C
C=======================================================================
C
      IF (INLMAX.GT.1) CALL  BDRY0 (DWORK(KD1),DWORK(KD2),KMBD,NMBD)
C
C=======================================================================
C
C
99998 END
C
C
C
************************************************************************
      SUBROUTINE XMADF3(KM1,KMASS1,KST1,KA1,KCOLA,KLDA,NA,NU,THSTEP)
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
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
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
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
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
         CALL LLC2(VWORK(KMASS1),VWORK(KA1),NA,1D0,1D0)
        ELSE
         CALL LCP2(VWORK(KMASS1),VWORK(KA1),NA)
        ENDIF
       ELSE
        IF (THSTEP.NE.0D0) THEN
         CALL LLC2(VWORK(KST1),VWORK(KA1),NA,THSTEP,0D0)
         DO 100 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
100      VWORK(KA1+INUA)=VWORK(KA1+INUA)+REAL(DWORK(KM1+INU-1))
        ELSE
         CALL LCL2(VWORK(KA1),NA)
         DO 110 INU=1,NU
         INUA=KWORK(KLDA+INU-1)-1
110      VWORK(KA1+INUA)=REAL(DWORK(KM1+INU-1))
        ENDIF
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
     *                               +DWORK(KMASS1+INA-1))
        ELSE
         CALL LCL2(VWORK(KA1),NA)
         DO 260 INA=1,NA
260      VWORK(KA1+INA-1)=REAL(DWORK(KMASS1+INA-1))
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
        CALL LCL2(VWORK(KA1),4*NA)
        DO 510 INU=1,NU
        INUA=KWORK(KLDA+INU-1)-1
        VWORK(KA1+INUA)=REAL(DWORK(KM1+INU-1))
        VWORK(3*NA+KA1+INUA)=REAL(DWORK(KM1+INU-1))
510     continue
       ELSE
        CALL LCL2(VWORK(KA1),4*NA)
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
      SUBROUTINE XDFKD(KB1,KB2,KCOLB,KLDB,KFP,KU1,KU2,NU,NP,TSTEPH)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)      
      PARAMETER (NNARR=299,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
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
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE
C-----------------------------------------------------------------------
C
      IF (IPRECB.EQ.0) THEN
       CALL  LTX39 (VWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU1),DWORK(KFP),1D0/TSTEPH,0D0)
       CALL  LTX39 (VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU2),DWORK(KFP),1D0/TSTEPH,1D0)
      ENDIF
C
C=======================================================================
C
      IF (IPRECB.EQ.1) THEN
       CALL  LTX19 (DWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU1),DWORK(KFP),1D0/TSTEPH,0D0)
       CALL  LTX19 (DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU2),DWORK(KFP),1D0/TSTEPH,1D0)
      ENDIF
C
C=======================================================================
C
      IF (IPRECB.EQ.2) THEN
       CALL BTMUL1(KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *             DWORK(L(LCORVG)),DWORK(KU1),DWORK(KU2),DWORK(KFP),
     *             NEL,NVT,NMT,1D0/TSTEPH)
      ENDIF
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE XDFKG(KM1,KB1,KB2,KCOLB,KLDB,KD1,KD2,KU1,KU2,KP,NU,NP,
     *                 KMBD,NMBD,TSTEPH)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)      
      PARAMETER (NNARR=299,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION KMBD(*)
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
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE
C-----------------------------------------------------------------------
C
      IF (IPRECB.EQ.0) THEN
       CALL  LAX39 (VWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,
     *              DWORK(KP),DWORK(KD1),TSTEPH,0D0)
       CALL  LAX39 (VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,
     *              DWORK(KP),DWORK(KD2),TSTEPH,0D0)
      ENDIF
C
C=======================================================================
C
      IF (IPRECB.EQ.1) THEN
       CALL  LAX19 (DWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,
     *              DWORK(KP),DWORK(KD1),TSTEPH,0D0)
       CALL  LAX19 (DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,
     *              DWORK(KP),DWORK(KD2),TSTEPH,0D0)
      ENDIF
C
C=======================================================================
C
      IF (IPRECB.EQ.2) THEN
       CALL BMUL1 (KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *             DWORK(L(LCORVG)),DWORK(KP),DWORK(KD1),DWORK(KD2),
     *             NEL,NVT,NMT,TSTEPH,0D0)
      ENDIF
C
      CALL  BDRY0 (DWORK(KD1),DWORK(KD2),KMBD,NMBD)
C
      DO 10 IU=1,NU
      DWORK(KU1+IU-1)=DWORK(KU1+IU-1)-DWORK(KD1+IU-1)/DWORK(KM1+IU-1)
      DWORK(KU2+IU-1)=DWORK(KU2+IU-1)-DWORK(KD2+IU-1)/DWORK(KM1+IU-1)
10    CONTINUE
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE RESDFK(D1,D2,F1,F2,NU,RESU1,RESU2)
************************************************************************
*    Purpose:  Computes the norms RESU
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z)
      DIMENSION D1(*),D2(*),F1(*),F2(*)
      SAVE 
C
      RESU1 =0D0
      RESU2 =0D0
C-----------------------------------------------------------------------
C     Compute the relative l2-norms  RESU1,RESU2
C-----------------------------------------------------------------------
      CALL LL21 (F1,NU,RESF1)
      CALL LL21 (F2,NU,RESF2)
      RESF=MAX(RESF1,RESF2)
      IF (ABS(RESF).LT.1D-8) RESF=1D0
C
      CALL LL21 (D1,NU,RESU1)
      CALL LL21 (D2,NU,RESU2)
      RESU1=RESU1/RESF
      RESU2=RESU2/RESF
C
      END
C
C
C
************************************************************************
      SUBROUTINE XDFKDV(KP,KFP,VMP,NP,A1,locny)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)      
      PARAMETER (NNARR=299,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VMP(*)
      DIMENSION VWORK(1),KWORK(1)
      DOUBLE PRECISION locny
      dimension locny(*)
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
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE
C-----------------------------------------------------------------------
C
      DO 10 IP=1,NP
      DWORK(KP+IP-1)=DWORK(KP+IP-1)+
     *       A1*NY*locny(IP)*DWORK(KFP+IP-1)/DBLE(VMP(IP))
10    CONTINUE
C
      END






