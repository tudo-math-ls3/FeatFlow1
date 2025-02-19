      DOUBLE PRECISION FUNCTION PARX(T1,T2,T3,IBCT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARX=T1
99999 END
C
C
C
      DOUBLE PRECISION FUNCTION PARY(T1,T2,T3,IBCT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARY=T2
99999 END
C
C
C
      DOUBLE PRECISION FUNCTION PARZ(T1,T2,T3,IBCT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
      PARAMETER (PI=3.1415926535897931D0)
C
C-----------------------------------------------------------------------
C
c      RETURN
C
      DO 10 IVT=1,NVT
      INPR=KNPR(IVT)
      IF (INPR.EQ.0) GOTO 10
C
      IF (NVEL.GE.3) THEN
       IEL=KVEL(3,IVT)
      ELSE
       IEL=0
      ENDIF
C
      PX=DCORVG(1,IVT)
      PY=DCORVG(2,IVT)
      PZ=DCORVG(3,IVT)
C
      IF ((ABS(PZ-3D0).GT.1D-8).AND.(ABS(PZ+3D0).GT.1D-8)) THEN
       PXM=0D0
       PYM=0D0
       RAD=SQRT(2D0)
       DL=SQRT((PX-PXM)**2+(PY-PYM)**2)
          DCORVG(1,IVT)=PXM+RAD/DL*(PX-PXM)
          DCORVG(2,IVT)=PYM+RAD/DL*(PY-PYM)
       GOTO 11
      ENDIF
C
      IF ((ABS(PZ-3D0).LE.1D-8).OR.(ABS(PZ+3D0).LE.1D-8)) THEN
       PXM=0D0
       PYM=0D0
       RAD=SQRT(2D0)
       IF (IEL.EQ.0) THEN
        DL=SQRT((PX-PXM)**2+(PY-PYM)**2)
           DCORVG(1,IVT)=PXM+RAD/DL*(PX-PXM)
           DCORVG(2,IVT)=PYM+RAD/DL*(PY-PYM)
        GOTO 11
       ENDIF
      ENDIF
C
C
C
11    IF ((ABS(PZ-3D0).GT.1D-8).AND.(ABS(PZ+3D0).GT.1D-8)) THEN
       PXM=0D0
       PYM=0D0
       if(sqrt(PX**2+PY**2).gt.0d0) then
         PXH=PX/SQRT(PX**2+PY**2)
         IF (PXH.GT. 1D0) PXH= 1D0
         IF (PXH.LT.-1D0) PXH=-1D0
         THETA=ACOS(PXH)
         IF (PY.LT.0D0) THETA=2D0*PI-THETA
         THETA=THETA+0.125D0*(PZ+3D0)*PI
         CFAC=0.5D0
         RAD=SQRT(2D0*(1D0+CFAC*ABS(SIN(2D0*THETA)))**2)
         DL=SQRT((PX-PXM)**2+(PY-PYM)**2)
         DCORVG(1,IVT)=PXM+RAD/DL*(PX-PXM)
         DCORVG(2,IVT)=PYM+RAD/DL*(PY-PYM)
       endif
       GOTO 10
      ENDIF
C
      IF ((ABS(PZ-3D0).LE.1D-8).OR.(ABS(PZ+3D0).LE.1D-8)) THEN
       PXM=0D0
       PYM=0D0
       if(sqrt(PX**2+PY**2).gt.0d0) then
         PXH=PX/SQRT(PX**2+PY**2)
         IF (PXH.GT. 1D0) PXH= 1D0
         IF (PXH.LT.-1D0) PXH=-1D0
         THETA=ACOS(PXH)
         IF (PY.LT.0D0) THETA=2D0*PI-THETA
         THETA=THETA+0.125D0*(PZ+3D0)*PI
         CFAC=0.5D0
         RAD=SQRT(2D0*(1D0+CFAC*ABS(SIN(2D0*THETA)))**2)
         IF (IEL.EQ.0) THEN
           DL=SQRT((PX-PXM)**2+(PY-PYM)**2)
           DCORVG(1,IVT)=PXM+RAD/DL*(PX-PXM)
           DCORVG(2,IVT)=PYM+RAD/DL*(PY-PYM)
           GOTO 10
         ENDIF
       endif
      ENDIF
C
10    CONTINUE     
C
C
C
      END
