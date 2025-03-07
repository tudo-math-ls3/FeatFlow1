      SUBROUTINE XCORI(PX1,PX2,PX3,PX4,PY1,PY2,PY3,PY4,
     *                 PXN1,PXN2,PXN3,PXN4,PYN1,PYN2,PYN3,PYN4,
     *                 DNELM1,DNELM2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
      A11=ABS(DNELM1)
      A21=1D0-A11
C
      A12=ABS(DNELM2)
      A22=1D0-A12
C
      PXN1= (1D0-A11)*(1D0-A11)*PX1+(1D0-A11)*(1D0-A21)*PX2
     *     +(1D0-A21)*(1D0-A21)*PX3+(1D0-A21)*(1D0-A11)*PX4
      PXN2= (1D0-A11)*(1D0-A12)*PX1+(1D0-A11)*(1D0-A22)*PX2
     *     +(1D0-A21)*(1D0-A22)*PX3+(1D0-A21)*(1D0-A12)*PX4
      PXN3= (1D0-A12)*(1D0-A12)*PX1+(1D0-A12)*(1D0-A22)*PX2
     *     +(1D0-A22)*(1D0-A22)*PX3+(1D0-A22)*(1D0-A12)*PX4
      PXN4= (1D0-A12)*(1D0-A11)*PX1+(1D0-A12)*(1D0-A21)*PX2
     *     +(1D0-A22)*(1D0-A21)*PX3+(1D0-A22)*(1D0-A11)*PX4
      PYN1= (1D0-A11)*(1D0-A11)*PY1+(1D0-A11)*(1D0-A21)*PY2
     *     +(1D0-A21)*(1D0-A21)*PY3+(1D0-A21)*(1D0-A11)*PY4
      PYN2= (1D0-A11)*(1D0-A12)*PY1+(1D0-A11)*(1D0-A22)*PY2
     *     +(1D0-A21)*(1D0-A22)*PY3+(1D0-A21)*(1D0-A12)*PY4
      PYN3= (1D0-A12)*(1D0-A12)*PY1+(1D0-A12)*(1D0-A22)*PY2
     *     +(1D0-A22)*(1D0-A22)*PY3+(1D0-A22)*(1D0-A12)*PY4
      PYN4= (1D0-A12)*(1D0-A11)*PY1+(1D0-A12)*(1D0-A21)*PY2
     *     +(1D0-A22)*(1D0-A21)*PY3+(1D0-A22)*(1D0-A11)*PY4
C
      END
C
C *****************************************************************
C
      SUBROUTINE XCORS(PX1,PX2,PY1,PY2,PXN1,PXN2,PYN1,PYN2,
     *                 INPR1,INPR2,INPR,DNELM1,DNELM2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
C
C
      A11=ABS(DNELM1)
c      IF (A11.LE.1D-5  ) A11=1D-5
c      IF (A11.GE.0.49D0) A11=0.49D0
      A21=1D0-A11
C
      A12=ABS(DNELM2)
c      IF (A12.LE.1D-4  ) A12=1D-4
c      IF (A12.GE.0.99D0) A12=0.99
      A22=1D0-A12
C
C
C
      PXN1=(1D0-A11)*PX1+(1D0-A21)*PX2
      PXN2=(1D0-A12)*PX1+(1D0-A22)*PX2
      PYN1=(1D0-A11)*PY1+(1D0-A21)*PY2
      PYN2=(1D0-A12)*PY1+(1D0-A22)*PY2
C
      IF (INPR1.EQ.INPR2) THEN
       INPR=INPR1
      ELSE
       INPR=0
      ENDIF
C
      END
C
C *****************************************************************
C
      SUBROUTINE XCORSP(DCORVG,DVBDP,KVBD,KMM,NBCT,NVBD,IEL,
     *                  IVT1,IVT2,IVTS1,IVTS2,INPR1,INPR2,
     *                  DNELM1,DNELM2,PARX,PARY,TMAX)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KMM(2,*),KVBD(*),DCORVG(2,*),DVBDP(*)
C
      A11=DNELM1
      A21=1D0-A11
C
      A12=DNELM2
      A22=1D0-A12
C
      IMM1=KMM(1,INPR1)
      IMM2=KMM(2,INPR1)
      IVBMM1=0
      IVBMM2=0
      INMMB=0
      DO 10 IVBD=1,NVBD
      IF (KVBD(IVBD).EQ.IMM1) THEN
       IVBMM1=IVBD
       INMMB=INMMB+1
      ENDIF
      IF (KVBD(IVBD).EQ.IMM2) THEN
       IVBMM2=IVBD
       INMMB=INMMB+1
      ENDIF
      IF (INMMB.EQ.2) GOTO 99
10    CONTINUE
C
99    ICVBD=0
C
      DO 100 IVBD=1,NVBD
      IVT=KVBD(IVBD)
C
      IF (IVT.EQ.IVT1) THEN
       ICVBD=ICVBD+1
       IVBD1=IVBD
      ENDIF
C
      IF (IVT.EQ.IVT2) THEN
       ICVBD=ICVBD+1
       IVBD2=IVBD
      ENDIF
C
      IF (ICVBD.GE.2) GOTO 110
100   CONTINUE
C
110   DPAR1=DVBDP(IVBD1)
      DPAR2=DVBDP(IVBD2)
C
      IF ((IVT1.EQ.KMM(1,INPR1)).AND.
     *    (DPAR2.GT.0.5D0*TMAX(INPR2))) DPAR1=TMAX(INPR1)
      IF ((IVT2.EQ.KMM(1,INPR2)).AND.
     *    (DPAR1.GT.0.5D0*TMAX(INPR1))) DPAR2=TMAX(INPR2)
C
      KVBD(NVBD+1)=IVTS1
      KVBD(NVBD+2)=IVTS2      
      DVBDP(NVBD+1)=(1D0-A11)*DPAR1+(1D0-A21)*DPAR2
      DVBDP(NVBD+2)=(1D0-A12)*DPAR1+(1D0-A22)*DPAR2
C
      DCORVG(1,IVTS1)=PARX(DVBDP(NVBD+1),INPR1)
      DCORVG(2,IVTS1)=PARY(DVBDP(NVBD+1),INPR1)
      DCORVG(1,IVTS2)=PARX(DVBDP(NVBD+2),INPR2)
      DCORVG(2,IVTS2)=PARY(DVBDP(NVBD+2),INPR2)
C
      IF (DVBDP(NVBD+1).GT.DVBDP(IVBMM2)) KMM(2,INPR1)=IVTS1
      IF (DVBDP(NVBD+2).GT.DVBDP(IVBMM2)) KMM(2,INPR1)=IVTS2
      IF (DVBDP(NVBD+1).LT.DVBDP(IVBMM1)) KMM(1,INPR1)=IVTS1
      IF (DVBDP(NVBD+2).LT.DVBDP(IVBMM1)) KMM(1,INPR1)=IVTS2
C
      NVBD=NVBD+2
C
      END
C
C *****************************************************************
C
      SUBROUTINE XCORM(PX1,PX2,PX3,PX4,PY1,PY2,PY3,PY4,
     *                 PXN1,PXN2,PYN1,PYN2,DNELM1,DNELM2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
      A11=ABS(DNELM1)
      A21=1D0-A11
C
      A12=ABS(DNELM2)
      A22=1D0-A12
C
      PXN1=0.5D0*( (1D0-A11)*PX1+(1D0-A21)*PX2
     *            +(1D0-A21)*PX3+(1D0-A11)*PX4)
      PXN2=0.5D0*( (1D0-A12)*PX1+(1D0-A22)*PX2
     *            +(1D0-A22)*PX3+(1D0-A12)*PX4)
      PYN1=0.5D0*( (1D0-A11)*PY1+(1D0-A21)*PY2
     *            +(1D0-A21)*PY3+(1D0-A11)*PY4)
      PYN2=0.5D0*( (1D0-A12)*PY1+(1D0-A22)*PY2
     *            +(1D0-A22)*PY3+(1D0-A12)*PY4)
C
      END
C
C *****************************************************************
C
      SUBROUTINE CHCOOR(DCORVG,KVERT,KMID,KADJ,KNPR,NEL,NVT)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DCORVG(2,*),KVERT(4,*),KMID(4,*),KADJ(4,*),KNPR(*)
C
C
C
      DO 100 IEL=1,NEL/4
C
      IEL1=IEL
      IEL2=KADJ(2,IEL1)
      IEL3=KADJ(2,IEL2)
      IEL4=KADJ(2,IEL3)
C
      IADJ1=KADJ(1,IEL1)
      IADJ2=KADJ(1,IEL2)
      IADJ3=KADJ(1,IEL3)
      IADJ4=KADJ(1,IEL4)
C
      IVT1=KVERT(2,IEL1)
      IVT2=KVERT(2,IEL2)
      IVT3=KVERT(2,IEL3)
      IVT4=KVERT(2,IEL4)
C
      NADJ0=0
      IF (IADJ1.EQ.0) NADJ0=NADJ0+1
      IF (IADJ2.EQ.0) NADJ0=NADJ0+1
      IF (IADJ3.EQ.0) NADJ0=NADJ0+1
      IF (IADJ4.EQ.0) NADJ0=NADJ0+1
C
C
      IF (NADJ0.EQ.0) GOTO 100
C
C
      PX1=DCORVG(1,IVT1)
      PX2=DCORVG(1,IVT2)
      PX3=DCORVG(1,IVT3)
      PX4=DCORVG(1,IVT4)
C
      PY1=DCORVG(2,IVT1)
      PY2=DCORVG(2,IVT2)
      PY3=DCORVG(2,IVT3)
      PY4=DCORVG(2,IVT4)
C
      IVTM=KVERT(3,IEL)
      PXM=DCORVG(1,IVTM)
      PYM=DCORVG(2,IVTM)
C
C
      IF (NADJ0.EQ.1) THEN
       IF ((IADJ1.EQ.0).OR.(IADJ3.EQ.0)) THEN
        PX=0.5D0*(PX1+PX3)
        PY=0.5D0*(PY1+PY3)
       ELSE
        PX=0.5D0*(PX2+PX4)
        PY=0.5D0*(PY2+PY4)
       ENDIF
       DCORVG(1,IVTM)=PX
       DCORVG(2,IVTM)=PY
      ENDIF
C
C      
      IF (NADJ0.GT.1) THEN
       PX=0.25D0*(PX1+PX2+PX3+PX4)
       PY=0.25D0*(PY1+PY2+PY3+PY4)
       DCORVG(1,IVTM)=PX
       DCORVG(2,IVTM)=PY
      ENDIF
C     
100   CONTINUE
C
C
C
      END
