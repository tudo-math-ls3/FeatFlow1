      DOUBLE PRECISION FUNCTION PARX(T1,T2,T3,IBCT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARX=T1
99999 END
C
C
C
      DOUBLE PRECISION FUNCTION PARY(T1,T2,T3,IBCT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARY=T2
99999 END
C
C
C
      DOUBLE PRECISION FUNCTION PARZ(T1,T2,T3,IBCT)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARZ=T3
      GOTO 99999
99999 END
C
C
C
************************************************************************
      SUBROUTINE   TRPARV (DCORVG,KNPR,KVEL,NVT,NVEL)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DCORVG(3,*),KNPR(*),KVEL(NVEL,*)
C
C-----------------------------------------------------------------------
C
C      RETURN
C
      DO 10 IVT=1,NVT
      INPR=KNPR(IVT)
      IF (INPR.EQ.0) GOTO 10
C
      IEL=KVEL(4,IVT)
      PX=DCORVG(1,IVT)
      PY=DCORVG(2,IVT)
      PZ=DCORVG(3,IVT)
C
      IF ((ABS(PZ-0.0D0).GT.1D-8).AND.(ABS(PZ-0.41D0).GT.1D-8).AND.
     *    (ABS(PX-0.0D0).GT.1D-8).AND.(ABS(PX-2.50D0).GT.1D-8).AND.
     *    (ABS(PY-0.0D0).GT.1D-8).AND.(ABS(PY-0.41D0).GT.1D-8)) THEN
       PXM=0.50D0
       PYM=0.20D0
       RAD=0.05D0
       DL=SQRT((PX-PXM)**2+(PY-PYM)**2)
       DCORVG(1,IVT)=PXM+RAD/DL*(PX-PXM)
       DCORVG(2,IVT)=PYM+RAD/DL*(PY-PYM)
       GOTO 10
      ENDIF
C
C
      IF ((ABS(PZ-0.0D0).LE.1D-8).OR.(ABS(PZ-0.41D0).LE.1D-8)) THEN
       PXM=0.50D0
       PYM=0.20D0
       RAD=0.05D0
       RADH=0.05D0+1D-8
       IF ((ABS(PX-PXM).LE.RADH).AND.(ABS(PY-PYM).LE.RADH).AND.
     *     (IEL.EQ.0)) THEN
        DL=SQRT((PX-PXM)**2+(PY-PYM)**2)
        DCORVG(1,IVT)=PXM+RAD/DL*(PX-PXM)
        DCORVG(2,IVT)=PYM+RAD/DL*(PY-PYM)
        GOTO 10
       ENDIF
      ENDIF
C
10    CONTINUE     
C
C
C
      END
