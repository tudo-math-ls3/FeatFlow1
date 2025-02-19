************************************************************************
      SUBROUTINE XMDF1T(KM1,KMASS1,KST1,KA1,KCOLA,KLDA,
     *                  KFT,KT,NA,NU,
     *                  THSTEP,TSTEPH)
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
     *              DWORK(KT),DWORK(KFT),1D0,1D0)
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
     *                DWORK(KT),DWORK(KFT),THSTEP,1D0)
C
         DO 200 INU=1,NU
         DWORK(KFT+INU-1)= DWORK(KFT+INU-1)
     *                    +DWORK(KM1+INU-1)*DWORK(KT+INU-1)
 
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
         DWORK(KFT+INU-1)= DWORK(KFT+INU-1)
     *                    +DWORK(KM1+INU-1)*DWORK(KT+INU-1)
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
     *                DWORK(KT),DWORK(KFT),THSTEP,1D0)
         CALL  LAX17 (DWORK(KMASS1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                DWORK(KT),DWORK(KFT),1D0,1D0)
C
         DO 250 INA=1,NA
250      VWORK(KA1+INA-1)=REAL(THSTEP*DWORK(KST1+INA-1)
     *                               +DWORK(KMASS1+INA-1))
        ELSE
         CALL  LAX17 (DWORK(KMASS1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                DWORK(KT),DWORK(KFT),1D0,1D0)
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
     *              DWORK(KT),DWORK(KFT),1D0,1D0)
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
     *              DWORK(KT),DWORK(KFT),1D0,1D0)
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
        DWORK(KFT+INU-1)= DWORK(KFT+INU-1)
     *                   +DWORK(KM1+INU-1)*DWORK(KT+INU-1)
510     CONTINUE
C
        CALL LCL2(VWORK(KA1),NA)
C
        DO 520 INU=1,NU
        INUA=KWORK(KLDA+INU-1)-1
520     VWORK(KA1+INUA)=REAL(DWORK(KM1+INU-1))
       ELSE
        CALL LCL2(VWORK(KA1),NA)
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
      SUBROUTINE XMDF2T(KM1,KMASS1,KST1,KA1,KCOLA,KLDA,KDT,KT,
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
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
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
     *               DWORK(KT),DWORK(KDT),-1D0,1D0)
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
     *                 DWORK(KT),DWORK(KDT),-THSTEP,1D0)
C
          DO 200 INU=1,NU
          DWORK(KDT+INU-1)= DWORK(KDT+INU-1)
     *                     -DWORK(KM1+INU-1)*DWORK(KT+INU-1)
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
          DWORK(KDT+INU-1)= DWORK(KDT+INU-1)
     *                     -DWORK(KM1+INU-1)*DWORK(KT+INU-1)
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
     *                 DWORK(KT),DWORK(KDT),-THSTEP,1D0)
          CALL  LAX17 (DWORK(KMASS1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DWORK(KT),DWORK(KDT),-1D0,1D0)
         ENDIF
C
         DO 250 INA=1,NA
250      VWORK(KA1+INA-1)=REAL(THSTEP*DWORK(KST1+INA-1)
     *                               +DWORK(KMASS1+INA-1))
        ELSE
         IF (INLMAX.GT.1) THEN
          CALL  LAX17 (DWORK(KMASS1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DWORK(KT),DWORK(KDT),-1D0,1D0)
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
     *               DWORK(KT),DWORK(KDT),-1D0,1D0)
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
     *               DWORK(KT),DWORK(KDT),-1D0,1D0)
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
         DWORK(KDT+INU-1)= DWORK(KDT+INU-1)
     *                    -DWORK(KM1+INU-1)*DWORK(KT+INU-1)
510      CONTINUE
        ENDIF
C
        CALL LCL2(VWORK(KA1),NA)
C
        DO 520 INU=1,NU
        INUA=KWORK(KLDA+INU-1)-1
520     VWORK(KA1+INUA)=REAL(DWORK(KM1+INU-1))
       ELSE
        CALL LCL2(VWORK(KA1),NA)
       ENDIF
C
      ENDIF
C
C=======================================================================
C
c      IF (INLMAX.GT.1) CALL  BDRY0T (DWORK(KDT),KMBD,NMBD)
       CALL  BDRY0T (DWORK(KDT),KMBD,NMBD,
     *            KWORK(L(KLNPR(ILEV))),DWORK(L(KLDBD(ILEV))))
C
C=======================================================================
C
C
99998 END

