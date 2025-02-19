************************************************************************
      DOUBLE PRECISION FUNCTION FDATIN(ITYP,IBLOC,X,Y,Z,TIMENS,RE)
*
*     Prescribed data for files coeff.f and bndry.f
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (PI=3.1415926535897931D0)
C     
C
C
      FDATIN=0D0
C
C
C
C=======================================================================
C *** Case 1: Velocity boundary values and/or exact solution
C=======================================================================
C
      IF (ITYP.EQ.1) THEN
C
       IF (IBLOC.EQ.1) THEN
        IF (X.EQ.0D0) 
     *      FDATIN=16D0*2.25D0/(0.41D0**4)*Y*(0.41D0-Y)*Z*(0.41D0-Z)
       ENDIF
C
       IF (IBLOC.EQ.2) THEN
        IF (X.EQ.0.0D0) FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.3) THEN
        IF (X.EQ.0.0D0) FDATIN=0D0
       ENDIF
C
      ENDIF
C
C
C=======================================================================
C *** Case 2: Velocity x-derivative of exact solution
C=======================================================================
C
      IF (ITYP.EQ.2) THEN
C
       IF (IBLOC.EQ.1) THEN
        IF (X.EQ.0.0D0) FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.2) THEN
        IF (X.EQ.0.0D0) FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.3) THEN
        IF (X.EQ.0.0D0) FDATIN=0D0
       ENDIF
      ENDIF
C
C
C=======================================================================
C *** Case 3: Velocity y-derivative of exact solution
C=======================================================================
C
      IF (ITYP.EQ.3) THEN
C
       IF (IBLOC.EQ.1) THEN
        FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.2) THEN
        FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.3) THEN
        FDATIN=0D0
       ENDIF
C
      ENDIF
C
C
C=======================================================================
C *** Case 4: Velocity z-derivative of exact solution
C=======================================================================
C
      IF (ITYP.EQ.4) THEN
C
       IF (IBLOC.EQ.1) THEN
        FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.2) THEN
        FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.3) THEN
        FDATIN=0D0
       ENDIF
C
      ENDIF
C
C
C=======================================================================
C *** Case 5: Exact pressure solution
C=======================================================================
C
      IF (ITYP.EQ.5) THEN
C
       FDATIN=0D0
C
      ENDIF
C
C
C=======================================================================
C *** Case 6: Right hand side for momentum equation
C=======================================================================
C
      IF (ITYP.EQ.6) THEN
C
       IF (IBLOC.EQ.1) THEN
        FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.2) THEN
        FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.3) THEN
        FDATIN=0D0
       ENDIF
C
      ENDIF
C
C
C=======================================================================
C *** Case 7: Right hand side for continuity equation
C=======================================================================
C
      IF (ITYP.EQ.7) THEN
C
       FDATIN=0D0
C
      ENDIF
C
C
C=======================================================================
C *** Case 8: Mean pressure values
C=======================================================================
C
      IF (ITYP.EQ.8) THEN
C
       FDATIN=0D0
C
      ENDIF
C
C
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE NEUDAT(IEL,INPR,PX,PY,PZ,TIMENS,IFLAG)
*
*     Neumann-boundary part
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C     
C=======================================================================
C *** Set Neumann-boundary parts
C=======================================================================
C
      IF (PX.EQ.2.5D0) THEN
       IFLAG=1      
      ENDIF
C
C
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE BDPDAT(IEL,INPR,PX,PY,PZ,TIMENS,IFLAG1,IFLAG2)
*
*     Pressure integral boundary part
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C     
C=======================================================================
C *** Set pressure integral boundary parts
C=======================================================================
C
      IF (((PZ.GT.0.00D0).AND.(PZ.LT.0.41D0)).AND.
     *    ((PX.GE.0.45D0).AND.(PX.LE.0.55D0)).AND.
     *    ((PY.GE.0.15D0).AND.(PY.LE.0.25D0))) THEN
       IFLAG1=1      
      ENDIF
C
      IF (((PZ.GT.0.00D0).AND.(PZ.LT.0.41D0)).AND.
     *    ((PX.GE.0.45D0).AND.(PX.LE.0.50D0)).AND.
     *    ((PY.GE.0.15D0).AND.(PY.LE.0.25D0))) THEN
       IFLAG2=1      
      ENDIF
C
C
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE BDFDAT(IEL,INPR,PX,PY,PZ,TIMENS,DNU,IFLAG,DPF1,DPF2)
*
*     lift and drag data
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C     
C=======================================================================
C *** Parameters for lift (DFW) and drag (DAW) in bdforc.f (INPR=2)
C ***
C *** dfw=2 int_s [dpf(1) dut/dn n_y - p n_x] ds / dpf(2)
C *** daw=2 int_s [dpf(1) dut/dn n_x + p n_y] ds / dpf(2)
C ***
C=======================================================================
C
      RHO  =1.0D0
      DIST =0.041D0
      UMEAN=1.0D0
C
      DPF1=RHO*DNU
      DPF2=RHO*DIST*UMEAN**2
C
C=======================================================================
C
      IF (((PZ.GT.0.00D0).AND.(PZ.LT.0.41D0)).AND.
     *    ((PX.GE.0.45D0).AND.(PX.LE.0.55D0)).AND.
     *    ((PY.GE.0.15D0).AND.(PY.LE.0.25D0))) THEN
       IFLAG=1      
      ENDIF
C
C
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE PTSDAT(TIMENS,DNU)
*
*     Data for Point-output (for fpost)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      COMMON /NSPTS/  KPU(2),KPP(4)
      SAVE
C
C     
C
C=======================================================================
C *** Points for velocity and pressure
C=======================================================================
C
      KPU(1)=3421
      KPU(2)=677
C
      KPP(1)=687
      KPP(2)=690
      KPP(3)=3498
      KPP(4)=677
C
C
C
99999 END
