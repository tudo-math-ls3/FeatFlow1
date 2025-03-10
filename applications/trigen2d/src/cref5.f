      SUBROUTINE CREF5A(DCORVG,DVBDP,KVBD,KVERT,KADJ,KNPR,KMM,KNELM,
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
      SUB='CREF5A'
      IF (ICHECK.GE.997) CALL OTRC('CREF5A','05/15/95')
C
C
C
      DO 100 IEL=1,NELOLD
C
      ITYP=KNELM(IEL)
      IF (ITYP.NE.10) GOTO 100
C
      IVT1=KVERT(1,IEL)
      IVT2=KVERT(2,IEL)
      IVT3=KVERT(3,IEL)
      IVT4=KVERT(4,IEL)
C
      PX1=DCORVG(1,IVT1)
      PX2=DCORVG(1,IVT2)
      PX3=DCORVG(1,IVT3)
      PX4=DCORVG(1,IVT4)
      PY1=DCORVG(2,IVT1)
      PY2=DCORVG(2,IVT2)
      PY3=DCORVG(2,IVT3)
      PY4=DCORVG(2,IVT4)
C
      CALL XCORI(PX1,PX2,PX3,PX4,PY1,PY2,PY3,PY4,
     *           PXN1,PXN2,PXN3,PXN4,PYN1,PYN2,PYN3,PYN4,
     *           DELMA,DELMB)
C
      DCORVG(1,NVT+1)=PXN1
      DCORVG(2,NVT+1)=PYN1
      DCORVG(1,NVT+2)=PXN2
      DCORVG(2,NVT+2)=PYN2
      DCORVG(1,NVT+3)=PXN3
      DCORVG(2,NVT+3)=PYN3
      DCORVG(1,NVT+4)=PXN4
      DCORVG(2,NVT+4)=PYN4
C
      KNPR(NVT+1)=0
      KNPR(NVT+2)=0
      KNPR(NVT+3)=0
      KNPR(NVT+4)=0
C
      KVERT(1,IEL)=NVT+1
      KVERT(2,IEL)=NVT+2
      KVERT(3,IEL)=NVT+3
      KVERT(4,IEL)=NVT+4
C
      KVERT(1,NEL+1)=IVT1
      KVERT(2,NEL+1)=IVT2
      KVERT(3,NEL+1)=NVT+2
      KVERT(4,NEL+1)=NVT+1
C
      KVERT(1,NEL+2)=IVT2
      KVERT(2,NEL+2)=IVT3
      KVERT(3,NEL+2)=NVT+3
      KVERT(4,NEL+2)=NVT+2
C
      KVERT(1,NEL+3)=IVT3
      KVERT(2,NEL+3)=IVT4
      KVERT(3,NEL+3)=NVT+4
      KVERT(4,NEL+3)=NVT+3
C
      KVERT(1,NEL+4)=IVT4
      KVERT(2,NEL+4)=IVT1
      KVERT(3,NEL+4)=NVT+1
      KVERT(4,NEL+4)=NVT+4
C
      NEL=NEL+4
      NVT=NVT+4
C
100   CONTINUE
C
C
C
99999 END
