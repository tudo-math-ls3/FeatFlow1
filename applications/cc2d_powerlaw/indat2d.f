************************************************************************
      DOUBLE PRECISION FUNCTION FDATIN(ITYP,IBLOC,X,Y,TIMENS,RE)
*
*     Prescribed data for files coeff.f and bndry.f
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (PI=3.1415926535897931D0)
c      EXTERNAL  f_u,f_v
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
         IF(Y.GE.0.999999D0) FDATIN=1.0d0
CCCC  4.0D0*X*(1-X)
       ENDIF

       IF (IBLOC.EQ.2) THEN
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
      ENDIF
C
C
C=======================================================================
C *** Case 4: Exact pressure solution
C=======================================================================
C
      IF (ITYP.EQ.4) THEN
C
       FDATIN=0D0
C
      ENDIF
C
C
C=======================================================================
C *** Case 5: Right hand side for momentum equation
C=======================================================================
C
      IF (ITYP.EQ.5) THEN
C
       IF (IBLOC.EQ.1) THEN
        FDATIN=0D0
       ENDIF
C
       IF (IBLOC.EQ.2) THEN
        FDATIN=0D0
       ENDIF
C
      ENDIF
C
C
C=======================================================================
C *** Case 6: Right hand side for continuity equation
C=======================================================================
C
      IF (ITYP.EQ.6) THEN
C
       FDATIN=0D0
C
      ENDIF
C
C
C=======================================================================
C *** Case 7: Mean pressure values
C=======================================================================
C
      IF (ITYP.EQ.7) THEN
       DPAR=X
       INPR=IBLOC
C
       IF ((DPAR.GT.1D0).AND.(DPAR.LT.2D0).AND.(INPR.EQ.1)) THEN
        FDATIN=0D0
       ENDIF
C
       IF ((DPAR.GT.3D0).AND.(DPAR.LT.4D0).AND.(INPR.EQ.1)) THEN
        FDATIN=0D0
       ENDIF
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
      SUBROUTINE NEUDAT(INPART,INPRN,DPARN1,DPARN2,TIMENS)
*
*     Neumann-boundary part
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C     
C
C
C=======================================================================
C *** Case 0: Set number of Neumann-boundary parts
C=======================================================================
C
C      IF (INPART.EQ.0) THEN
C       INPART=1      
C      ENDIF
C
C
C=======================================================================
C *** Case <>0: Specify Neumann-boundary parts
C=======================================================================
C
C      IF (INPART.GT.0) THEN
C
C       IF (INPART.EQ.1) THEN
C        DPARN1=1D0
C        DPARN2=2D0
C        INPRN =1
C       ENDIF
C       
C      ENDIF
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
*     Data for Point-output (for fpost and bdpres and bdforc)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      COMMON /NSPTS/  KPU(2),KPP(4),KPX(4),KPI(2),DPI(2,2),DPF(2)
      SAVE
C
      EXTERNAL UE
C
C
C
C=======================================================================
C *** Points for velocity, pressure and flux-tracing
C=======================================================================
C
      KPU(1)=2
      KPU(2)=12
C
      KPP(1)=2
      KPP(2)=12
      KPP(3)=24
      KPP(4)=12
C
      KPX(1)=2
      KPX(2)=12
C
C
C=======================================================================
C *** Parameters for 2 integral pressures in bdpres.f
C=======================================================================
C
      KPI(1)  =0
      DPI(1,1)=0D0
      DPI(2,1)=1D0
C
      KPI(2)  =0
      DPI(1,2)=0D0
      DPI(2,2)=1D0
C
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
      DIST =1.0D0
      UMEAN=1.0D0
C
C      DPF(1)=RHO*DNU
C      DPF(2)=RHO*DIST*UMEAN**2
C
      DPF(1)=1d0
      DPF(2)=1d0
C
C

C
99999 END
