      SUBROUTINE CHEC3A(DCORVG,DVBDP,KVBD,KVERT,KADJ,KNPR,KMM,KNELM,
     *                  ICH3A)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=4,NNELM=1000)
      DIMENSION DCORVG(2,*),DVBDP(*)
      DIMENSION KVERT(NNVE,*),KADJ(NNVE,*),KNPR(*),KMM(2,*),KVBD(*),
     *          KNELM(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
C
      COMMON /ELEMOD/ ITYPEL,NELMOD,DELMA,DELMB,KELMOD(NNELM)
C
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='CHEC3A'
      IF (ICHECK.GE.997) CALL OTRC('CHEC3A','05/15/95')
C
C
C
      ICH3A=0
      DO 100 IEL=1,NEL
C
      ITYP=KNELM(IEL)
      IF (ITYP.NE.-1) GOTO 100
C
      IADJ1=KADJ(1,IEL)
      IADJ2=KADJ(2,IEL)
      IADJ3=KADJ(3,IEL)
      IADJ4=KADJ(4,IEL)
C
      NADJ0=0
      IF (IADJ1.EQ.0) NADJ0=NADJ0+1
      IF (IADJ2.EQ.0) NADJ0=NADJ0+1
      IF (IADJ3.EQ.0) NADJ0=NADJ0+1
      IF (IADJ4.EQ.0) NADJ0=NADJ0+1
C
      IF (NADJ0.NE.1) THEN
       WRITE(6,*) 'ELEMENT ',IEL,' CANCELLED FROM LIST IN CHEC3A'
       KNELM(IEL)=0
       ICH3A=ICH3A+1
      ENDIF
C
100   CONTINUE
C
      END
C
C *****************************************************************
C
      SUBROUTINE CHEC3B(DCORVG,DVBDP,KVBD,KVERT,KADJ,KNPR,KMM,KNELM,
     *                  ICH3B)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=4,NNELM=1000)
      DIMENSION DCORVG(2,*),DVBDP(*)
      DIMENSION KVERT(NNVE,*),KADJ(NNVE,*),KNPR(*),KMM(2,*),KVBD(*),
     *          KNELM(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
C
      COMMON /ELEMOD/ ITYPEL,NELMOD,DELMA,DELMB,KELMOD(NNELM)
C
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='CHEC3B'
      IF (ICHECK.GE.997) CALL OTRC('CHEC3B','05/15/95')
C
C
C
      ICH3B=0
      DO 100 IEL=1,NEL
C
      ITYP=KNELM(IEL)
      IF (ITYP.EQ.-1) THEN
       ICH3B=ICH3B+1
       GOTO 100
      ENDIF
C
      IADJ1=KADJ(1,IEL)
      IADJ2=KADJ(2,IEL)
      IADJ3=KADJ(3,IEL)
      IADJ4=KADJ(4,IEL)
C
      ITYP1=0
      IF (IADJ1.NE.0) ITYP1=KNELM(IADJ1)
      IF (ITYP1.EQ.-1) THEN
       DO 110 JADJ=1,4
       IADJ=KADJ(JADJ,IADJ1)
       IF (IADJ.EQ.0) GOTO 112
110    CONTINUE
       WRITE(6,*) 'ERROR 1 IN CHEC3B'
C
112    JADJ1=MOD(JADJ,4)+1
       IADJH1=KADJ(JADJ1,IADJ1)
       JADJ3=MOD(JADJ1+1,4)+1
       IADJH3=KADJ(JADJ3,IADJ1)
       IF ((IADJH1.NE.IEL).AND.(IADJH3.NE.IEL)) ITYP1=0
      ELSE
       ITYP1=0
      ENDIF
C
      ITYP2=0
      IF (IADJ2.NE.0) ITYP2=KNELM(IADJ2)
      IF (ITYP2.EQ.-1) THEN
       DO 120 JADJ=1,4
       IADJ=KADJ(JADJ,IADJ2)
       IF (IADJ.EQ.0) GOTO 122
120    CONTINUE
       WRITE(6,*) 'ERROR 2 IN CHEC3B'
C
122    JADJ1=MOD(JADJ,4)+1
       IADJH1=KADJ(JADJ1,IADJ2)
       JADJ3=MOD(JADJ1+1,4)+1
       IADJH3=KADJ(JADJ3,IADJ2)
       IF ((IADJH1.NE.IEL).AND.(IADJH3.NE.IEL)) ITYP2=0
      ELSE
       ITYP2=0
      ENDIF
C
      ITYP3=0
      IF (IADJ3.NE.0) ITYP3=KNELM(IADJ3)
      IF (ITYP3.EQ.-1) THEN
       DO 130 JADJ=1,4
       IADJ=KADJ(JADJ,IADJ3)
       IF (IADJ.EQ.0) GOTO 132
130    CONTINUE
       WRITE(6,*) 'ERROR 3 IN CHEC3B'
C
132    JADJ1=MOD(JADJ,4)+1
       IADJH1=KADJ(JADJ1,IADJ3)
       JADJ3=MOD(JADJ1+1,4)+1
       IADJH3=KADJ(JADJ3,IADJ3)
       IF ((IADJH1.NE.IEL).AND.(IADJH3.NE.IEL)) ITYP3=0
      ELSE
       ITYP3=0
      ENDIF
C
      ITYP4=0
      IF (IADJ4.NE.0) ITYP4=KNELM(IADJ4)
      IF (ITYP4.EQ.-1) THEN
       DO 140 JADJ=1,4
       IADJ=KADJ(JADJ,IADJ4)
       IF (IADJ.EQ.0) GOTO 142
140    CONTINUE
       WRITE(6,*) 'ERROR 4 IN CHEC3B'
C
142    JADJ1=MOD(JADJ,4)+1
       IADJH1=KADJ(JADJ1,IADJ4)
       JADJ3=MOD(JADJ1+1,4)+1
       IADJH3=KADJ(JADJ3,IADJ4)
       IF ((IADJH1.NE.IEL).AND.(IADJH3.NE.IEL)) ITYP4=0
      ELSE
       ITYP4=0
      ENDIF
C
      ITYP=ITYP1+ITYP2+ITYP3+ITYP4
C
      NADJ0=0
      IF (IADJ1.EQ.0) NADJ0=NADJ0+1
      IF (IADJ2.EQ.0) NADJ0=NADJ0+1
      IF (IADJ3.EQ.0) NADJ0=NADJ0+1
      IF (IADJ4.EQ.0) NADJ0=NADJ0+1
C
      IF (ITYP.EQ.-1) THEN
       IF (NADJ0.EQ.1) THEN
        KNELM(IEL)=-2
        ICH3B=ICH3B+2
       ELSE
        KNELM(IEL)=-3
C        ICH3B=ICH3B+2
        WRITE(6,*) 'IEL= ',IEL, 
     *             ' : REFINEMENT PROCEDURE NOT YET IMPLEMENTED !!!'
       ENDIF
      ENDIF
C
      IF (ITYP.EQ.-2) THEN
       IF (((ITYP1+ITYP3).EQ.-2).OR.((ITYP2+ITYP4).EQ.-2)) THEN
        KNELM(IEL)=-1
        ICH3B=ICH3B+2
       ELSE
        IF (NADJ0.EQ.2) THEN 
         KNELM(IEL)=-5
         ICH3B=ICH3B+3
        ELSE
         KNELM(IEL)=-4
         ICH3B=ICH3B+2
        ENDIF
       ENDIF
      ENDIF
C
      IF ((ITYP.EQ.-3).OR.(ITYP.EQ.-4)) THEN
       WRITE(6,*) 'ERROR 5 IN CHEC3B'
      ENDIF
C
100   CONTINUE
C
      END
C
C *****************************************************************
C
      SUBROUTINE CREF3A(DCORVG,DVBDP,KVBD,KVERT,KADJ,KNPR,KMM,KNELM,
     *                  KINDEL,NELOLD,NVTOLD)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=4,NNELM=1000)
      DIMENSION DCORVG(2,*),DVBDP(*)
      DIMENSION KVERT(NNVE,*),KADJ(NNVE,*),KNPR(*),KMM(2,*),KVBD(*),
     *          KNELM(*),KINDEL(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
C
      COMMON /ELEMOD/ ITYPEL,NELMOD,DELMA,DELMB,KELMOD(NNELM)
C
      EXTERNAL PARX,PARY,TMAX
C
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='CREF3A'
      IF (ICHECK.GE.997) CALL OTRC('CREF3A','05/15/95')
C
C
C
      INDIEL=0
      DO 100 IEL=1,NELOLD
C
      ITYP=KNELM(IEL)
      IF (ITYP.NE.-1) GOTO 100
C
      INDIEL=INDIEL+1
      KINDEL(IEL)=INDIEL
C
      DO 105 JADJ=1,4
      IADJ=KADJ(JADJ,IEL)
      IF (IADJ.EQ.0) GOTO 107
105   CONTINUE
      WRITE(6,*) 'ERROR 1 IN CREF3A'
C
107   JADJ1=MOD(JADJ,4)+1
      IADJ1=KADJ(JADJ1,IEL)
      JADJ3=MOD(JADJ1+1,4)+1
      IADJ3=KADJ(JADJ3,IEL)
C
      ITYP1=0
      ITYP3=0
      IF (IADJ1.NE.0) ITYP1=KNELM(IADJ1)
      IF (IADJ3.NE.0) ITYP3=KNELM(IADJ3)
C
C
      IVT1=KVERT(1,IEL)
      IVT2=KVERT(2,IEL)
      IVT3=KVERT(3,IEL)
      IVT4=KVERT(4,IEL)
C
      IF (JADJ1.EQ.1) THEN
       IVTH1=IVT1
       IVTH2=IVT2
       IVTH3=IVT3
       IVTH4=IVT4
      ENDIF
C
      IF (JADJ1.EQ.2) THEN
       IVTH1=IVT2
       IVTH2=IVT3
       IVTH3=IVT4
       IVTH4=IVT1
      ENDIF
C
      IF (JADJ1.EQ.3) THEN
       IVTH1=IVT3
       IVTH2=IVT4
       IVTH3=IVT1
       IVTH4=IVT2
      ENDIF
C
      IF (JADJ1.EQ.4) THEN
       IVTH1=IVT4
       IVTH2=IVT1
       IVTH3=IVT2
       IVTH4=IVT3
      ENDIF
C
      INPR1=KNPR(IVTH1)
      INPR2=KNPR(IVTH2)
      INPR3=KNPR(IVTH3)
      INPR4=KNPR(IVTH4)
C
      PX1=DCORVG(1,IVTH1)
      PX2=DCORVG(1,IVTH2)
      PX3=DCORVG(1,IVTH3)
      PX4=DCORVG(1,IVTH4)
      PY1=DCORVG(2,IVTH1)
      PY2=DCORVG(2,IVTH2)
      PY3=DCORVG(2,IVTH3)
      PY4=DCORVG(2,IVTH4)
C
      IF ((ITYP1.EQ.-1).AND.(KINDEL(IADJ1).LT.INDIEL)
     *                 .AND.(KINDEL(IADJ1).GT.0)) THEN
       DO 110 JADJH=1,4
       IADJH=KADJ(JADJH,IADJ1)
       IF (IADJH.EQ.0) GOTO 112
110    CONTINUE
       WRITE(6,*) 'ERROR 2 IN CREF3A'
C
112    JADJH1=MOD(JADJH,4)+1
       IADJH1=KADJ(JADJH1,IADJ1)
       JADJH3=MOD(JADJH1+1,4)+1
       IADJH3=KADJ(JADJH3,IADJ1)
C
       IF (IADJH1.EQ.IEL) THEN
        IVTS1=KVERT(2,IADJ1)
       ELSE
        IVTS1=KVERT(3,IADJ1)
       ENDIF
      ELSE
       IVTS1=NVT+1
       CALL XCORS(PX1,PX2,PY1,PY2,PXS1,PXS2,PYS1,PYS2,
     *            INPR1,INPR2,INPRS1,DELMA,DELMB)
       DCORVG(1,IVTS1)=PXS1
       DCORVG(2,IVTS1)=PYS1
       KNPR(IVTS1)=INPRS1
       NVT=NVT+1
      ENDIF
C
      IF ((ITYP3.EQ.-1).AND.(KINDEL(IADJ3).LT.INDIEL)
     *                 .AND.(KINDEL(IADJ3).GT.0)) THEN
       DO 120 JADJH=1,4
       IADJH=KADJ(JADJH,IADJ3)
       IF (IADJH.EQ.0) GOTO 122
120    CONTINUE
       WRITE(6,*) 'ERROR 3 IN CREF3A'
C
122    JADJH1=MOD(JADJH,4)+1
       IADJH1=KADJ(JADJH1,IADJ3)
       JADJH3=MOD(JADJH1+1,4)+1
       IADJH3=KADJ(JADJH3,IADJ3)
C
       IF (IADJH1.EQ.IEL) THEN
        IVTS3=KVERT(2,IADJ3)
       ELSE
        IVTS3=KVERT(3,IADJ3)
       ENDIF
      ELSE
       IVTS3=NVT+1
       CALL XCORS(PX3,PX4,PY3,PY4,PXS5,PXS6,PYS5,PYS6,
     *            INPR3,INPR4,INPRS3,1D0-DELMB,1D0-DELMA)
       DCORVG(1,IVTS3)=PXS6
       DCORVG(2,IVTS3)=PYS6
       KNPR(IVTS3)=INPRS3
       NVT=NVT+1
      ENDIF
C
C
      KVERT(1,IEL)=IVTH1
      KVERT(2,IEL)=IVTS1
      KVERT(3,IEL)=IVTS3
      KVERT(4,IEL)=IVTH4
C
      KVERT(1,NEL+1)=IVTH3
      KVERT(2,NEL+1)=IVTS3
      KVERT(3,NEL+1)=IVTS1
      KVERT(4,NEL+1)=IVTH2
C
      NEL=NEL+1
C
100   CONTINUE
C
99999 END
C
C *****************************************************************
C
      SUBROUTINE CREF3B(DCORVG,DVBDP,KVBD,KVERT,KADJ,KNPR,KMM,KNELM,
     *                  KINDEL,NELOLD,NVTOLD)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=4,NNELM=1000)
      DIMENSION DCORVG(2,*),DVBDP(*)
      DIMENSION KVERT(NNVE,*),KADJ(NNVE,*),KNPR(*),KMM(2,*),KVBD(*),
     *          KNELM(*),KINDEL(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
C
      COMMON /ELEMOD/ ITYPEL,NELMOD,DELMA,DELMB,KELMOD(NNELM)
C
      EXTERNAL PARX,PARY,TMAX
C
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='CREF3B'
      IF (ICHECK.GE.997) CALL OTRC('CREF3B','05/15/95')
C
C
C
      DO 100 IEL=1,NELOLD
C
      ITYP=KNELM(IEL)
      IF (ITYP.NE.-2) GOTO 100
C
102   DO 110 JADJ1=1,4
      IADJ1=KADJ(JADJ1,IEL)
      IF (KNELM(IADJ1).EQ.-1) GOTO 112
110   CONTINUE
      WRITE(6,*) 'ERROR 1 IN CREF3E'
C
112   DO 120 JADJH=1,4
      IADJH=KADJ(JADJH,IADJ1)
      IF (IADJH.EQ.0) GOTO 122
120   CONTINUE
      WRITE(6,*) 'ERROR 2 IN CREF3E'
C
122   JADJH1=MOD(JADJH,4)+1
      IADJH1=KADJ(JADJH1,IADJ1)
      JADJH3=MOD(JADJH1+1,4)+1
      IADJH3=KADJ(JADJH3,IADJ1)
C
      IF (IADJH1.EQ.IEL) THEN
       ICORM=2
       IVTS1=KVERT(2,IADJ1)
      ELSE
       ICORM=1
       IVTS1=KVERT(3,IADJ1)
      ENDIF
C
C
      IVT1=KVERT(1,IEL)
      IVT2=KVERT(2,IEL)
      IVT3=KVERT(3,IEL)
      IVT4=KVERT(4,IEL)
C
      IF (JADJ1.EQ.1) THEN
       IVTH1=IVT1
       IVTH2=IVT2
       IVTH3=IVT3
       IVTH4=IVT4
      ENDIF
C
      IF (JADJ1.EQ.2) THEN
       IVTH1=IVT2
       IVTH2=IVT3
       IVTH3=IVT4
       IVTH4=IVT1
      ENDIF
C
      IF (JADJ1.EQ.3) THEN
       IVTH1=IVT3
       IVTH2=IVT4
       IVTH3=IVT1
       IVTH4=IVT2
      ENDIF
C
      IF (JADJ1.EQ.4) THEN
       IVTH1=IVT4
       IVTH2=IVT1
       IVTH3=IVT2
       IVTH4=IVT3
      ENDIF
C
      INPR1=KNPR(IVTH1)
      INPR2=KNPR(IVTH2)
      INPR3=KNPR(IVTH3)
      INPR4=KNPR(IVTH4)
C
      PX1=DCORVG(1,IVTH1)
      PX2=DCORVG(1,IVTH2)
      PX3=DCORVG(1,IVTH3)
      PX4=DCORVG(1,IVTH4)
      PY1=DCORVG(2,IVTH1)
      PY2=DCORVG(2,IVTH2)
      PY3=DCORVG(2,IVTH3)
      PY4=DCORVG(2,IVTH4)
C
      IF (ICORM.EQ.1) THEN
       CALL XCORI(PX1,PX2,PX3,PX4,PY1,PY2,PY3,PY4,
     *            PXN1,PXN2,PXN3,PXN4,PYN1,PYN2,PYN3,PYN4,
     *            DELMA,DELMB)
       IVTN=NVT+1
       DCORVG(1,IVTN)=PXN1
       DCORVG(2,IVTN)=PYN1
      ELSE
       CALL XCORI(PX2,PX3,PX4,PX1,PY2,PY3,PY4,PY1,
     *            PXN1,PXN2,PXN3,PXN4,PYN1,PYN2,PYN3,PYN4,
     *            DELMA,DELMB)
       IVTN=NVT+1
       DCORVG(1,IVTN)=PXN1
       DCORVG(2,IVTN)=PYN1
      ENDIF
C
      KNPR(IVTN)=0
      NVT=NVT+1
C
C
      IF (ICORM.EQ.1) THEN
       IVTS =NVT+1
       IVTSH=NVT+2
       CALL XCORS(PX1,PX4,PY1,PY4,PXS1,PXS2,PYS1,PYS2,
     *            INPR1,INPR4,INPRS1,DELMA,1D0-DELMA)
       DCORVG(1,IVTS)=PXS1
       DCORVG(2,IVTS)=PYS1
       KNPR(IVTS)=INPRS1
       NVT=NVT+1
       CALL XCORSP(DCORVG,DVBDP,KVBD,KMM,NBCT,NVBD,IEL,
     *             IVTH1,IVTH4,IVTS,IVTSH,INPR1,INPR4,
     *             DELMA,1D0-DELMA,PARX,PARY,TMAX)
       NVBD=NVBD-1
      ELSE
       IVTS= NVT+1
       IVTSH=NVT+2
       CALL XCORS(PX2,PX3,PY2,PY3,PXS1,PXS2,PYS1,PYS2,
     *            INPR2,INPR3,INPRS2,DELMA,1D0-DELMA)
       DCORVG(1,IVTS)=PXS1
       DCORVG(2,IVTS)=PYS1
       KNPR(IVTS)=INPRS2
       NVT=NVT+1
       CALL XCORSP(DCORVG,DVBDP,KVBD,KMM,NBCT,NVBD,IEL,
     *             IVTH2,IVTH3,IVTS,IVTSH,INPR2,INPR3,
     *             DELMA,1D0-DELMA,PARX,PARY,TMAX)
       NVBD=NVBD-1
      ENDIF
C
C
      IF (ICORM.EQ.1) THEN
       KVERT(1,IEL)=IVTH1
       KVERT(2,IEL)=IVTS1
       KVERT(3,IEL)=IVTN
       KVERT(4,IEL)=IVTS
C
       KVERT(1,NEL+1)=IVTH2
       KVERT(2,NEL+1)=IVTH3
       KVERT(3,NEL+1)=IVTN
       KVERT(4,NEL+1)=IVTS1
C
       KVERT(1,NEL+2)=IVTH3
       KVERT(2,NEL+2)=IVTH4
       KVERT(3,NEL+2)=IVTS
       KVERT(4,NEL+2)=IVTN
      ELSE
       KVERT(1,IEL)=IVTH2
       KVERT(2,IEL)=IVTS
       KVERT(3,IEL)=IVTN
       KVERT(4,IEL)=IVTS1
C
       KVERT(1,NEL+1)=IVTH3
       KVERT(2,NEL+1)=IVTH4
       KVERT(3,NEL+1)=IVTN
       KVERT(4,NEL+1)=IVTS
C
       KVERT(1,NEL+2)=IVTH4
       KVERT(2,NEL+2)=IVTH1
       KVERT(3,NEL+2)=IVTS1
       KVERT(4,NEL+2)=IVTN
      ENDIF
C
C
C
      NEL=NEL+2
C
100   CONTINUE
C
C
C
99999 END
C
C *****************************************************************
C
      SUBROUTINE CREF3D(DCORVG,DVBDP,KVBD,KVERT,KADJ,KNPR,KMM,KNELM,
     *                  KINDEL,NELOLD,NVTOLD)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=4,NNELM=1000)
      DIMENSION DCORVG(2,*),DVBDP(*)
      DIMENSION KVERT(NNVE,*),KADJ(NNVE,*),KNPR(*),KMM(2,*),KVBD(*),
     *          KNELM(*),KINDEL(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
C
      COMMON /ELEMOD/ ITYPEL,NELMOD,DELMA,DELMB,KELMOD(NNELM)
C
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='CREF3D'
      IF (ICHECK.GE.997) CALL OTRC('CREF3D','05/15/95')
C
C
C
      DO 100 IEL=1,NELOLD
C
      ITYP=KNELM(IEL)
      IF (ITYP.NE.-4) GOTO 100
C
      JADS1=1
102   DO 110 JADJ1=JADS1,4
      IADJ1=KADJ(JADJ1,IEL)
      IF (KNELM(IADJ1).EQ.-1) GOTO 112
110   CONTINUE
      WRITE(6,*) 'ERROR 1 IN CREF3D'
C
112   DO 120 JADJH=1,4
      IADJH=KADJ(JADJH,IADJ1)
      IF (IADJH.EQ.0) GOTO 122
120   CONTINUE
      WRITE(6,*) 'ERROR 2 IN CREF3D'
C
122   JADJH1=MOD(JADJH,4)+1
      IADJH1=KADJ(JADJH1,IADJ1)
      JADJH3=MOD(JADJH1+1,4)+1
      IADJH3=KADJ(JADJH3,IADJ1)
C
      IF (IADJH1.EQ.IEL) THEN
       ICORM=1
       IVTS1=KVERT(2,IADJ1)
      ELSE
       ICORM=2
       IVTS1=KVERT(3,IADJ1)
      ENDIF
C
C
      JADJ2=MOD(JADJ1,4)+1
      IADJ2=KADJ(JADJ2,IEL)
C
      IF (KNELM(IADJ2).NE.-1) THEN
       JADS1=JADJ1+1
       GOTO 102
      ENDIF
C
      DO 130 JADJH=1,4
      IADJH=KADJ(JADJH,IADJ2)
      IF (IADJH.EQ.0) GOTO 132
130   CONTINUE
      WRITE(6,*) 'ERROR 3 IN CREF3D'
C
132   JADJH1=MOD(JADJH,4)+1
      IADJH1=KADJ(JADJH1,IADJ2)
      JADJH3=MOD(JADJH1+1,4)+1
      IADJH3=KADJ(JADJH3,IADJ2)
C
      IF (IADJH1.EQ.IEL) THEN
       IVTS4=KVERT(2,IADJ2)
      ELSE
       IVTS4=KVERT(3,IADJ2)
      ENDIF
C
C
      IVT1=KVERT(1,IEL)
      IVT2=KVERT(2,IEL)
      IVT3=KVERT(3,IEL)
      IVT4=KVERT(4,IEL)
C
      IF (JADJ1.EQ.1) THEN
       IVTH1=IVT2
       IVTH2=IVT3
       IVTH3=IVT4
       IVTH4=IVT1
      ENDIF
C
      IF (JADJ1.EQ.2) THEN
       IVTH1=IVT3
       IVTH2=IVT4
       IVTH3=IVT1
       IVTH4=IVT2
      ENDIF
C
      IF (JADJ1.EQ.3) THEN
       IVTH1=IVT4
       IVTH2=IVT1
       IVTH3=IVT2
       IVTH4=IVT3
      ENDIF
C
      IF (JADJ1.EQ.4) THEN
       IVTH1=IVT1
       IVTH2=IVT2
       IVTH3=IVT3
       IVTH4=IVT4
      ENDIF
C
      PX1=DCORVG(1,IVTH1)
      PX2=DCORVG(1,IVTH2)
      PX3=DCORVG(1,IVTH3)
      PX4=DCORVG(1,IVTH4)
      PY1=DCORVG(2,IVTH1)
      PY2=DCORVG(2,IVTH2)
      PY3=DCORVG(2,IVTH3)
      PY4=DCORVG(2,IVTH4)
C
      IF (ICORM.EQ.1) THEN
       CALL XCORI(PX1,PX2,PX3,PX4,PY1,PY2,PY3,PY4,
     *            PXN1,PXN2,PXN3,PXN4,PYN1,PYN2,PYN3,PYN4,
     *            DELMA,DELMB)
      ELSE
       CALL XCORI(PX1,PX2,PX3,PX4,PY1,PY2,PY3,PY4,
     *            PXN1,PXN2,PXN3,PXN4,PYN1,PYN2,PYN3,PYN4,
     *            1D0-DELMB,1D0-DELMA)
      ENDIF
C
      IVTN1=NVT+1
      NVT=NVT+1
C
      DCORVG(1,IVTN1)=PXN1
      DCORVG(2,IVTN1)=PYN1
      KNPR(IVTN1)=0
C
      KVERT(1,IEL)=IVTH1
      KVERT(2,IEL)=IVTS4
      KVERT(3,IEL)=IVTN1
      KVERT(4,IEL)=IVTS1
C
      KVERT(1,NEL+1)=IVTS4
      KVERT(2,NEL+1)=IVTH2
      KVERT(3,NEL+1)=IVTH3
      KVERT(4,NEL+1)=IVTN1
C
      KVERT(1,NEL+2)=IVTH3
      KVERT(2,NEL+2)=IVTH4
      KVERT(3,NEL+2)=IVTS1
      KVERT(4,NEL+2)=IVTN1
C
      NEL=NEL+2
C
100   CONTINUE
C
99999 END
C
C *****************************************************************
C
      SUBROUTINE CREF3E(DCORVG,DVBDP,KVBD,KVERT,KADJ,KNPR,KMM,KNELM,
     *                  KINDEL,NELOLD,NVTOLD)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=4,NNELM=1000)
      DIMENSION DCORVG(2,*),DVBDP(*)
      DIMENSION KVERT(NNVE,*),KADJ(NNVE,*),KNPR(*),KMM(2,*),KVBD(*),
     *          KNELM(*),KINDEL(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
C
      COMMON /ELEMOD/ ITYPEL,NELMOD,DELMA,DELMB,KELMOD(NNELM)
C
      EXTERNAL PARX,PARY,TMAX
C
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='CREF3E'
      IF (ICHECK.GE.997) CALL OTRC('CREF3E','05/15/95')
C
C
C
      DO 100 IEL=1,NELOLD
C
      ITYP=KNELM(IEL)
      IF (ITYP.NE.-5) GOTO 100
C
      JADS1=1
102   DO 110 JADJ1=JADS1,4
      IADJ1=KADJ(JADJ1,IEL)
      IF (KNELM(IADJ1).EQ.-1) GOTO 112
110   CONTINUE
      WRITE(6,*) 'ERROR 1 IN CREF3E'
C
112   DO 120 JADJH=1,4
      IADJH=KADJ(JADJH,IADJ1)
      IF (IADJH.EQ.0) GOTO 122
120   CONTINUE
      WRITE(6,*) 'ERROR 2 IN CREF3E'
C
122   JADJH1=MOD(JADJH,4)+1
      IADJH1=KADJ(JADJH1,IADJ1)
      JADJH3=MOD(JADJH1+1,4)+1
      IADJH3=KADJ(JADJH3,IADJ1)
C
      IF (IADJH1.EQ.IEL) THEN
       ICORM=1
       IVTS2=KVERT(2,IADJ1)
      ELSE
       ICORM=2
       IVTS2=KVERT(3,IADJ1)
      ENDIF
C
C
      JADJ2=MOD(JADJ1,4)+1
      IADJ2=KADJ(JADJ2,IEL)
C
      IF (KNELM(IADJ2).NE.-1) THEN
       JADS1=JADJ1+1
       GOTO 102
      ENDIF
C
      DO 130 JADJH=1,4
      IADJH=KADJ(JADJH,IADJ2)
      IF (IADJH.EQ.0) GOTO 132
130   CONTINUE
      WRITE(6,*) 'ERROR 3 IN CREF3E'
C
132   JADJH1=MOD(JADJH,4)+1
      IADJH1=KADJ(JADJH1,IADJ2)
      JADJH3=MOD(JADJH1+1,4)+1
      IADJH3=KADJ(JADJH3,IADJ2)
C
      IF (IADJH1.EQ.IEL) THEN
       IVTS3=KVERT(2,IADJ2)
      ELSE
       IVTS3=KVERT(3,IADJ2)
      ENDIF
C
C
      IVT1=KVERT(1,IEL)
      IVT2=KVERT(2,IEL)
      IVT3=KVERT(3,IEL)
      IVT4=KVERT(4,IEL)
C
      IF (JADJ1.EQ.1) THEN
       IVTH1=IVT2
       IVTH2=IVT3
       IVTH3=IVT4
       IVTH4=IVT1
      ENDIF
C
      IF (JADJ1.EQ.2) THEN
       IVTH1=IVT3
       IVTH2=IVT4
       IVTH3=IVT1
       IVTH4=IVT2
      ENDIF
C
      IF (JADJ1.EQ.3) THEN
       IVTH1=IVT4
       IVTH2=IVT1
       IVTH3=IVT2
       IVTH4=IVT3
      ENDIF
C
      IF (JADJ1.EQ.4) THEN
       IVTH1=IVT1
       IVTH2=IVT2
       IVTH3=IVT3
       IVTH4=IVT4
      ENDIF
C
      INPR1=KNPR(IVTH1)
      INPR2=KNPR(IVTH2)
      INPR3=KNPR(IVTH3)
      INPR4=KNPR(IVTH4)
C
      PX1=DCORVG(1,IVTH1)
      PX2=DCORVG(1,IVTH2)
      PX3=DCORVG(1,IVTH3)
      PX4=DCORVG(1,IVTH4)
      PY1=DCORVG(2,IVTH1)
      PY2=DCORVG(2,IVTH2)
      PY3=DCORVG(2,IVTH3)
      PY4=DCORVG(2,IVTH4)
C
      IF (ICORM.EQ.1) THEN
       CALL XCORI(PX1,PX2,PX3,PX4,PY1,PY2,PY3,PY4,
     *            PXN1,PXN2,PXN3,PXN4,PYN1,PYN2,PYN3,PYN4,
     *            DELMA,DELMB)
      ELSE
       CALL XCORI(PX1,PX2,PX3,PX4,PY1,PY2,PY3,PY4,
     *            PXN1,PXN2,PXN3,PXN4,PYN1,PYN2,PYN3,PYN4,
     *            1D0-DELMB,1D0-DELMA)
      ENDIF
C
C
C
      IVTN3=NVT+1
      DCORVG(1,IVTN3)=PXN3
      DCORVG(2,IVTN3)=PYN3
C
      KNPR(IVTN3)=0
C
      KVERT(1,IEL)=IVTH1
      KVERT(2,IEL)=IVTS3
      KVERT(3,IEL)=IVTN3
      KVERT(4,IEL)=IVTS2
C
      NVT=NVT+1
C
C
      IVTS5=NVT+1
      IVTS6=NVT+2
      CALL XCORS(PX3,PX2,PY3,PY2,PXS5,PXS6,PYS5,PYS6,
     *           INPR3,INPR2,INPRS3,DELMA,DELMB)
      DCORVG(1,IVTS5)=PXS5
      DCORVG(2,IVTS5)=PYS5
      KNPR(IVTS5)=INPRS3
      NVT=NVT+1
      CALL XCORSP(DCORVG,DVBDP,KVBD,KMM,NBCT,NVBD,IEL,
     *            IVTH3,IVTH2,IVTS5,IVTS6,INPR3,INPR2,
     *            DELMA,DELMB,PARX,PARY,TMAX)
      NVBD=NVBD-1
C
      IVTS7=NVT+1
      IVTS8=NVT+2
      CALL XCORS(PX3,PX4,PY3,PY4,PXS7,PXS8,PYS7,PYS8,
     *           INPR3,INPR4,INPRS4,DELMA,DELMB)
      DCORVG(1,IVTS7)=PXS7
      DCORVG(2,IVTS7)=PYS7
      KNPR(IVTS7)=INPRS4
      NVT=NVT+1
      CALL XCORSP(DCORVG,DVBDP,KVBD,KMM,NBCT,NVBD,IEL,
     *            IVTH3,IVTH4,IVTS7,IVTS8,INPR3,INPR4,
     *            DELMA,DELMB,PARX,PARY,TMAX)
      NVBD=NVBD-1
C
C
      KVERT(1,NEL+1)=IVTH2
      KVERT(2,NEL+1)=IVTS5
      KVERT(3,NEL+1)=IVTN3
      KVERT(4,NEL+1)=IVTS3
C
      KVERT(1,NEL+2)=IVTH3
      KVERT(2,NEL+2)=IVTS7
      KVERT(3,NEL+2)=IVTN3
      KVERT(4,NEL+2)=IVTS5
C
      KVERT(1,NEL+3)=IVTH4
      KVERT(2,NEL+3)=IVTS2
      KVERT(3,NEL+3)=IVTN3
      KVERT(4,NEL+3)=IVTS7
C
C
      NEL=NEL+3
C
100   CONTINUE
C
C
C
99999 END
C
C
C
