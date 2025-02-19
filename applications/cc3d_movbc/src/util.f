************************************************************************
      SUBROUTINE   INTUV (DU1,DU2,DU3,V1,V2,V3,KAUX,NVT,NEL,NVBD,NABD,
     *                    KAREA,KVERT,KNPR,KVBD,KABD,DCORVG,
     *                    UE,INEUM,KELBD)
************************************************************************
*    Purpose:  - Interpolates the solution vector (DU1,DU2,DU3) to
*                the (REAL) vector (V1,V2,V3) of dimension NVT with
*                values in the vertices
*              - the values of vertices at the boundary are computed
*                via the exact function UE
*
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=8,NNAE=6,NNLEV=9)
C
      REAL  V1,V2,V3
      PARAMETER (A1=2D0/3D0,A2=1D0/3D0)
      DIMENSION DU1(*),DU2(*),DU3(*),V1(*),V2(*),V3(*),KAUX(*),
     *          KVERT(NNVE,*),KNPR(*),KAREA(6,*),DCORVG(3,*),
     *          KVBD(*),KABD(*),KELBD(*)
C
      SAVE
C-----------------------------------------------------------------------
C *** Zero initialization of (V1,V2,V3,AUX)
      CALL  LCL2 (V1,NVT)
      CALL  LCL2 (V2,NVT)
      CALL  LCL2 (V3,NVT)
      CALL  LCL3 (KAUX,NVT)
C
C
      DO 10 IEL=1,NEL
C
      IA1=KAREA(1,IEL)
      IA2=KAREA(2,IEL)
      IA3=KAREA(3,IEL)
      IA4=KAREA(4,IEL)
      IA5=KAREA(5,IEL)
      IA6=KAREA(6,IEL)
C
      IV1=KVERT(1,IEL)
      IV2=KVERT(2,IEL)
      IV3=KVERT(3,IEL)
      IV4=KVERT(4,IEL)
      IV5=KVERT(5,IEL)
      IV6=KVERT(6,IEL)
      IV7=KVERT(7,IEL)
      IV8=KVERT(8,IEL)
C
      DUH1=DU1(IA1)
      DUH2=DU1(IA2)
      DUH3=DU1(IA3)
      DUH4=DU1(IA4)
      DUH5=DU1(IA5)
      DUH6=DU1(IA6)
C
      DVH1=DU2(IA1)
      DVH2=DU2(IA2)
      DVH3=DU2(IA3)
      DVH4=DU2(IA4)
      DVH5=DU2(IA5)
      DVH6=DU2(IA6)
C
      DWH1=DU3(IA1)
      DWH2=DU3(IA2)
      DWH3=DU3(IA3)
      DWH4=DU3(IA4)
      DWH5=DU3(IA5)
      DWH6=DU3(IA6)
C
      KAUX(IV1)=KAUX(IV1)+1
      KAUX(IV2)=KAUX(IV2)+1
      KAUX(IV3)=KAUX(IV3)+1
      KAUX(IV4)=KAUX(IV4)+1
      KAUX(IV5)=KAUX(IV5)+1
      KAUX(IV6)=KAUX(IV6)+1
      KAUX(IV7)=KAUX(IV7)+1
      KAUX(IV8)=KAUX(IV8)+1
C
C
C
C
      R1=DUH1+DUH2+DUH3
      R2=DUH1+DUH2+DUH5
      R3=DUH1+DUH3+DUH4
      R4=DUH1+DUH4+DUH5
      R5=DUH2+DUH3+DUH6
      R6=DUH2+DUH5+DUH6
      R7=DUH3+DUH4+DUH6
      R8=DUH4+DUH5+DUH6
C
      V1(IV1)=V1(IV1) + A1*R2 - A2*R7
      V1(IV2)=V1(IV2) + A1*R1 - A2*R8
      V1(IV3)=V1(IV3) + A1*R3 - A2*R6
      V1(IV4)=V1(IV4) + A1*R4 - A2*R5
      V1(IV5)=V1(IV5) + A1*R6 - A2*R3
      V1(IV6)=V1(IV6) + A1*R5 - A2*R4
      V1(IV7)=V1(IV7) + A1*R7 - A2*R2
      V1(IV8)=V1(IV8) + A1*R8 - A2*R1
C
      S1=DVH1+DVH2+DVH3
      S2=DVH1+DVH2+DVH5
      S3=DVH1+DVH3+DVH4
      S4=DVH1+DVH4+DVH5
      S5=DVH2+DVH3+DVH6
      S6=DVH2+DVH5+DVH6
      S7=DVH3+DVH4+DVH6
      S8=DVH4+DVH5+DVH6
C
      V2(IV1)=V2(IV1) + A1*S2 - A2*S7
      V2(IV2)=V2(IV2) + A1*S1 - A2*S8
      V2(IV3)=V2(IV3) + A1*S3 - A2*S6
      V2(IV4)=V2(IV4) + A1*S4 - A2*S5
      V2(IV5)=V2(IV5) + A1*S6 - A2*S3
      V2(IV6)=V2(IV6) + A1*S5 - A2*S4
      V2(IV7)=V2(IV7) + A1*S7 - A2*S2
      V2(IV8)=V2(IV8) + A1*S8 - A2*S1
C
      T1=DWH1+DWH2+DWH3
      T2=DWH1+DWH2+DWH5
      T3=DWH1+DWH3+DWH4
      T4=DWH1+DWH4+DWH5
      T5=DWH2+DWH3+DWH6
      T6=DWH2+DWH5+DWH6
      T7=DWH3+DWH4+DWH6
      T8=DWH4+DWH5+DWH6
C
      V3(IV1)=V3(IV1) + A1*T2 - A2*T7
      V3(IV2)=V3(IV2) + A1*T1 - A2*T8
      V3(IV3)=V3(IV3) + A1*T3 - A2*T6
      V3(IV4)=V3(IV4) + A1*T4 - A2*T5
      V3(IV5)=V3(IV5) + A1*T6 - A2*T3
      V3(IV6)=V3(IV6) + A1*T5 - A2*T4
      V3(IV7)=V3(IV7) + A1*T7 - A2*T2
      V3(IV8)=V3(IV8) + A1*T8 - A2*T1
C
10    CONTINUE
C
      DO 20  IV=1,NVT
      V1(IV)=V1(IV)/REAL(KAUX(IV))
      V2(IV)=V2(IV)/REAL(KAUX(IV))
      V3(IV)=V3(IV)/REAL(KAUX(IV))
 20   CONTINUE
C
19    CONTINUE
      IF (INEUM.EQ.0) THEN
      DO 30  IV=1,NVT
      INPR=KNPR(IV)
      IF (INPR.EQ.0) GOTO 30
      X=DCORVG(1,IV)
      Y=DCORVG(2,IV)
      Z=DCORVG(3,IV)
      V1(IV)=UE(X,Y,Z,1)
      V2(IV)=UE(X,Y,Z,2)
      V3(IV)=UE(X,Y,Z,3)
  30  CONTINUE
      ENDIF
C
C
      IF (INEUM.EQ.1) THEN
      DO 31 IAT=1,NABD
      IABD=KABD(IAT)
      IELBD=KELBD(IAT)
        IF ((IELBD.GT.0).AND.(IABD.GT.0)) THEN
         DO 32 IA=1,6
         IAH=KAREA(IA,IELBD)
         IF (IAH.EQ.IABD) THEN
          IF (IA.EQ.1) THEN
          IVT1=KVERT(1,IELBD)
          IVT2=KVERT(2,IELBD)
          IVT3=KVERT(3,IELBD)
          IVT4=KVERT(4,IELBD)
          ENDIF
          IF (IA.EQ.2) THEN
          IVT1=KVERT(1,IELBD)
          IVT2=KVERT(2,IELBD)
          IVT3=KVERT(6,IELBD)
          IVT4=KVERT(5,IELBD)
          ENDIF
          IF (IA.EQ.3) THEN
          IVT1=KVERT(2,IELBD)
          IVT2=KVERT(3,IELBD)
          IVT3=KVERT(6,IELBD)
          IVT4=KVERT(7,IELBD)
          ENDIF
          IF (IA.EQ.4) THEN
          IVT1=KVERT(3,IELBD)
          IVT2=KVERT(4,IELBD)
          IVT3=KVERT(8,IELBD)
          IVT4=KVERT(7,IELBD)
          ENDIF
          IF (IA.EQ.5) THEN
          IVT1=KVERT(4,IELBD)
          IVT2=KVERT(1,IELBD)
          IVT3=KVERT(5,IELBD)
          IVT4=KVERT(8,IELBD)
          ENDIF
          IF (IA.EQ.6) THEN
          IVT1=KVERT(5,IELBD)
          IVT2=KVERT(6,IELBD)
          IVT3=KVERT(7,IELBD)
          IVT4=KVERT(8,IELBD)
          ENDIF
           X1=DCORVG(1,IVT1)
           Y1=DCORVG(2,IVT1)
           Z1=DCORVG(3,IVT1)
           V1(IVT1)=UE(X1,Y1,Z1,1)
           V2(IVT1)=UE(X1,Y1,Z1,2)
           V3(IVT1)=UE(X1,Y1,Z1,3)
           X2=DCORVG(1,IVT2)
           Y2=DCORVG(2,IVT2)
           Z2=DCORVG(3,IVT2)
           V1(IVT2)=UE(X2,Y2,Z2,1)
           V2(IVT2)=UE(X2,Y2,Z2,2)
           V3(IVT2)=UE(X2,Y2,Z2,3)
           X3=DCORVG(1,IVT3)
           Y3=DCORVG(2,IVT3)
           Z3=DCORVG(3,IVT3)
           V1(IVT3)=UE(X3,Y3,Z3,1)
           V2(IVT3)=UE(X3,Y3,Z3,2)
           V3(IVT3)=UE(X3,Y3,Z3,3)
           X4=DCORVG(1,IVT4)
           Y4=DCORVG(2,IVT4)
           Z4=DCORVG(3,IVT4)
           V1(IVT4)=UE(X4,Y4,Z4,1)
           V2(IVT4)=UE(X4,Y4,Z4,2)
           V3(IVT4)=UE(X4,Y4,Z4,3)
          ENDIF
32       CONTINUE
         ENDIF
31     CONTINUE
      ENDIF
C
      do 333 II=1,NVT
      INPR=KNPR(II)
      IF (INPR.EQ.11) THEN
      V1(II)=0D0
      V2(II)=0D0
      V3(II)=0D0
      ENDIF
333   CONTINUE
C
      END
c
c
c
c
************************************************************************
      SUBROUTINE INTUVD (DU1,DU2,DU3,DL1,DL2,DL3,DAUX,NVT,NEL,NVBD,NABD,
     *                   KAREA,KVERT,KNPR,KVBD,KABD,DCORVG,
     *                   UE,INEUM,KELBD)
************************************************************************
*    Purpose:  - Interpolates the solution vector (DU1,DU2,DU3) to
*                the vector (DL1,DL2,DL3) of dimension NVT with
*                values in the vertices
*              - the values of vertices at the boundary are computed
*                via the exact function UE
*
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=8,NNAE=6,NNLEV=9)
C
      PARAMETER (A1=2D0/3D0,A2=1D0/3D0)
      DIMENSION DU1(*),DU2(*),DU3(*),DL1(*),DL2(*),DL3(*),DAUX(*),
     *          KVERT(NNVE,*),KNPR(*),KAREA(6,*),DCORVG(3,*),
     *          KVBD(*),KABD(*),KELBD(*)
C
      SAVE
C-----------------------------------------------------------------------
      DO 1 IVT=1,NVT
      DL1 (IVT)=0D0
      DL2 (IVT)=0D0
      DL3 (IVT)=0D0
1     DAUX(IVT)=0D0
C
C
      DO 10 IEL=1,NEL
C
      IA1=KAREA(1,IEL)
      IA2=KAREA(2,IEL)
      IA3=KAREA(3,IEL)
      IA4=KAREA(4,IEL)
      IA5=KAREA(5,IEL)
      IA6=KAREA(6,IEL)
C
      IV1=KVERT(1,IEL)
      IV2=KVERT(2,IEL)
      IV3=KVERT(3,IEL)
      IV4=KVERT(4,IEL)
      IV5=KVERT(5,IEL)
      IV6=KVERT(6,IEL)
      IV7=KVERT(7,IEL)
      IV8=KVERT(8,IEL)
C
      DUH1=DU1(IA1)
      DUH2=DU1(IA2)
      DUH3=DU1(IA3)
      DUH4=DU1(IA4)
      DUH5=DU1(IA5)
      DUH6=DU1(IA6)
C
      DVH1=DU2(IA1)
      DVH2=DU2(IA2)
      DVH3=DU2(IA3)
      DVH4=DU2(IA4)
      DVH5=DU2(IA5)
      DVH6=DU2(IA6)
C
      DWH1=DU3(IA1)
      DWH2=DU3(IA2)
      DWH3=DU3(IA3)
      DWH4=DU3(IA4)
      DWH5=DU3(IA5)
      DWH6=DU3(IA6)
C
      DAUX(IV1)=DAUX(IV1)+1D0
      DAUX(IV2)=DAUX(IV2)+1D0
      DAUX(IV3)=DAUX(IV3)+1D0
      DAUX(IV4)=DAUX(IV4)+1D0
      DAUX(IV5)=DAUX(IV5)+1D0
      DAUX(IV6)=DAUX(IV6)+1D0
      DAUX(IV7)=DAUX(IV7)+1D0
      DAUX(IV8)=DAUX(IV8)+1D0
C
C
C
C
      R1=DUH1+DUH2+DUH3
      R2=DUH1+DUH2+DUH5
      R3=DUH1+DUH3+DUH4
      R4=DUH1+DUH4+DUH5
      R5=DUH2+DUH3+DUH6
      R6=DUH2+DUH5+DUH6
      R7=DUH3+DUH4+DUH6
      R8=DUH4+DUH5+DUH6
C
      DL1(IV1)=DL1(IV1) + A1*R2 - A2*R7
      DL1(IV2)=DL1(IV2) + A1*R1 - A2*R8
      DL1(IV3)=DL1(IV3) + A1*R3 - A2*R6
      DL1(IV4)=DL1(IV4) + A1*R4 - A2*R5
      DL1(IV5)=DL1(IV5) + A1*R6 - A2*R3
      DL1(IV6)=DL1(IV6) + A1*R5 - A2*R4
      DL1(IV7)=DL1(IV7) + A1*R7 - A2*R2
      DL1(IV8)=DL1(IV8) + A1*R8 - A2*R1
C
      S1=DVH1+DVH2+DVH3
      S2=DVH1+DVH2+DVH5
      S3=DVH1+DVH3+DVH4
      S4=DVH1+DVH4+DVH5
      S5=DVH2+DVH3+DVH6
      S6=DVH2+DVH5+DVH6
      S7=DVH3+DVH4+DVH6
      S8=DVH4+DVH5+DVH6
C
      DL2(IV1)=DL2(IV1) + A1*S2 - A2*S7
      DL2(IV2)=DL2(IV2) + A1*S1 - A2*S8
      DL2(IV3)=DL2(IV3) + A1*S3 - A2*S6
      DL2(IV4)=DL2(IV4) + A1*S4 - A2*S5
      DL2(IV5)=DL2(IV5) + A1*S6 - A2*S3
      DL2(IV6)=DL2(IV6) + A1*S5 - A2*S4
      DL2(IV7)=DL2(IV7) + A1*S7 - A2*S2
      DL2(IV8)=DL2(IV8) + A1*S8 - A2*S1
C
      T1=DWH1+DWH2+DWH3
      T2=DWH1+DWH2+DWH5
      T3=DWH1+DWH3+DWH4
      T4=DWH1+DWH4+DWH5
      T5=DWH2+DWH3+DWH6
      T6=DWH2+DWH5+DWH6
      T7=DWH3+DWH4+DWH6
      T8=DWH4+DWH5+DWH6
C
      DL3(IV1)=DL3(IV1) + A1*T2 - A2*T7
      DL3(IV2)=DL3(IV2) + A1*T1 - A2*T8
      DL3(IV3)=DL3(IV3) + A1*T3 - A2*T6
      DL3(IV4)=DL3(IV4) + A1*T4 - A2*T5
      DL3(IV5)=DL3(IV5) + A1*T6 - A2*T3
      DL3(IV6)=DL3(IV6) + A1*T5 - A2*T4
      DL3(IV7)=DL3(IV7) + A1*T7 - A2*T2
      DL3(IV8)=DL3(IV8) + A1*T8 - A2*T1
C
10    CONTINUE
C
      DO 20  IV=1,NVT
      DL1(IV)=DL1(IV)/DAUX(IV)
      DL2(IV)=DL2(IV)/DAUX(IV)
      DL3(IV)=DL3(IV)/DAUX(IV)
 20   CONTINUE
C
19    CONTINUE
      IF (INEUM.EQ.0) THEN
      DO 30  IV=1,NVT
      INPR=KNPR(IV)
      IF (INPR.EQ.0) GOTO 30
      X=DCORVG(1,IV)
      Y=DCORVG(2,IV)
      Z=DCORVG(3,IV)
      DL1(IV)=UE(X,Y,Z,1)
      DL2(IV)=UE(X,Y,Z,2)
      DL3(IV)=UE(X,Y,Z,3)
  30  CONTINUE
      ENDIF
C
C
      IF (INEUM.EQ.1) THEN
      DO 31 IAT=1,NABD
      IABD=KABD(IAT)
      IELBD=KELBD(IAT)
        IF ((IELBD.GT.0).AND.(IABD.GT.0)) THEN
         DO 32 IA=1,6
         IAH=KAREA(IA,IELBD)
         IF (IAH.EQ.IABD) THEN
          IF (IA.EQ.1) THEN
          IVT1=KVERT(1,IELBD)
          IVT2=KVERT(2,IELBD)
          IVT3=KVERT(3,IELBD)
          IVT4=KVERT(4,IELBD)
          ENDIF
          IF (IA.EQ.2) THEN
          IVT1=KVERT(1,IELBD)
          IVT2=KVERT(2,IELBD)
          IVT3=KVERT(6,IELBD)
          IVT4=KVERT(5,IELBD)
          ENDIF
          IF (IA.EQ.3) THEN
          IVT1=KVERT(2,IELBD)
          IVT2=KVERT(3,IELBD)
          IVT3=KVERT(6,IELBD)
          IVT4=KVERT(7,IELBD)
          ENDIF
          IF (IA.EQ.4) THEN
          IVT1=KVERT(3,IELBD)
          IVT2=KVERT(4,IELBD)
          IVT3=KVERT(8,IELBD)
          IVT4=KVERT(7,IELBD)
          ENDIF
          IF (IA.EQ.5) THEN
          IVT1=KVERT(4,IELBD)
          IVT2=KVERT(1,IELBD)
          IVT3=KVERT(5,IELBD)
          IVT4=KVERT(8,IELBD)
          ENDIF
          IF (IA.EQ.6) THEN
          IVT1=KVERT(5,IELBD)
          IVT2=KVERT(6,IELBD)
          IVT3=KVERT(7,IELBD)
          IVT4=KVERT(8,IELBD)
          ENDIF
           X1=DCORVG(1,IVT1)
           Y1=DCORVG(2,IVT1)
           Z1=DCORVG(3,IVT1)
           DL1(IVT1)=UE(X1,Y1,Z1,1)
           DL2(IVT1)=UE(X1,Y1,Z1,2)
           DL3(IVT1)=UE(X1,Y1,Z1,3)
           X2=DCORVG(1,IVT2)
           Y2=DCORVG(2,IVT2)
           Z2=DCORVG(3,IVT2)
           DL1(IVT2)=UE(X2,Y2,Z2,1)
           DL2(IVT2)=UE(X2,Y2,Z2,2)
           DL3(IVT2)=UE(X2,Y2,Z2,3)
           X3=DCORVG(1,IVT3)
           Y3=DCORVG(2,IVT3)
           Z3=DCORVG(3,IVT3)
           DL1(IVT3)=UE(X3,Y3,Z3,1)
           DL2(IVT3)=UE(X3,Y3,Z3,2)
           DL3(IVT3)=UE(X3,Y3,Z3,3)
           X4=DCORVG(1,IVT4)
           Y4=DCORVG(2,IVT4)
           Z4=DCORVG(3,IVT4)
           DL1(IVT4)=UE(X4,Y4,Z4,1)
           DL2(IVT4)=UE(X4,Y4,Z4,2)
           DL3(IVT4)=UE(X4,Y4,Z4,3)
          ENDIF
32       CONTINUE
         ENDIF
31     CONTINUE
      ENDIF
C
      do 334 II=1,NVT
      INPR=KNPR(II)
      IF (INPR.EQ.11) THEN
      DL1(II)=0D0
      DL2(II)=0D0
      DL3(II)=0D0
      ENDIF
334   CONTINUE
C
      END
c
c
c
c
************************************************************************
      SUBROUTINE   INTPV (DP,VPL,VAUX,AVOL,KVERT,KNPR)
************************************************************************
*    Purpose:  - Interpolates the solution pressure DP to
*                the (REAL) vector VPL of dimension NVT with
*                values in the vertices
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      REAL AVOL
C
      PARAMETER (NNVE=8)
C
      REAL  VPL,VAUX
      DIMENSION DP(*),VPL(*),VAUX(*),KVERT(NNVE,*),AVOL(*),KNPR(*)
C
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
C
      SAVE
C-----------------------------------------------------------------------
      DO 10 IEL=1,NEL
C
      VPIEL=REAL(DP(IEL))
      VAVOL=AVOL(IEL)
C
      IV1=KVERT(1,IEL)
      IV2=KVERT(2,IEL)
      IV3=KVERT(3,IEL)
      IV4=KVERT(4,IEL)
      IV5=KVERT(5,IEL)
      IV6=KVERT(6,IEL)
      IV7=KVERT(7,IEL)
      IV8=KVERT(8,IEL)
C
      VPL(IV1)=VPL(IV1)+0.125E0*VAVOL*VPIEL
      VPL(IV2)=VPL(IV2)+0.125E0*VAVOL*VPIEL
      VPL(IV3)=VPL(IV3)+0.125E0*VAVOL*VPIEL
      VPL(IV4)=VPL(IV4)+0.125E0*VAVOL*VPIEL
      VPL(IV5)=VPL(IV5)+0.125E0*VAVOL*VPIEL
      VPL(IV6)=VPL(IV6)+0.125E0*VAVOL*VPIEL
      VPL(IV7)=VPL(IV7)+0.125E0*VAVOL*VPIEL
      VPL(IV8)=VPL(IV8)+0.125E0*VAVOL*VPIEL
C
      VAUX(IV1)=VAUX(IV1)+0.125E0*VAVOL
      VAUX(IV2)=VAUX(IV2)+0.125E0*VAVOL
      VAUX(IV3)=VAUX(IV3)+0.125E0*VAVOL
      VAUX(IV4)=VAUX(IV4)+0.125E0*VAVOL
      VAUX(IV5)=VAUX(IV5)+0.125E0*VAVOL
      VAUX(IV6)=VAUX(IV6)+0.125E0*VAVOL
      VAUX(IV7)=VAUX(IV7)+0.125E0*VAVOL
      VAUX(IV8)=VAUX(IV8)+0.125E0*VAVOL
C
10    CONTINUE
C
C
      DO 20 IVT=1,NVT
20    VPL(IVT)=VPL(IVT)/VAUX(IVT)
C
      do 335 II=1,NVT
      INPR=KNPR(II)
      IF (INPR.EQ.11) THEN
      VPL(II)=0D0
      ENDIF
335   CONTINUE
C
C      VPH=0E0
C      DO 30 IVT=1,NVT
C30    VPH=VPH+VPL(IVT)
C      VMWP=VPH/REAL(NVT)
C
C
C      DO 40 IVT=1,NVT
C40    VPL(IVT)=VPL(IVT)-VMWP
C
C
      END
c
c
************************************************************************
      SUBROUTINE  INTPVD (DP,DPL,DAUX,AVOL,KVERT,KNPR)
************************************************************************
*    Purpose:  - Interpolates the solution pressure DP to
*                the (REAL) vector VPL of dimension NVT with
*                values in the vertices
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL AVOL
      PARAMETER (NNVE=8)
C
      DIMENSION DP(*),DPL(*),DAUX(*),KVERT(NNVE,*),AVOL(*),KNPR(*)
C
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
C
      SAVE
C-----------------------------------------------------------------------
      DO 1 IVT=1,NVT
      DPL (IVT)=0D0
1     DAUX(IVT)=0D0
C
C
      DO 10 IEL=1,NEL
C
      DPIEL=DP(IEL)
      DAVOL=DBLE(AVOL(IEL))
C
      IV1=KVERT(1,IEL)
      IV2=KVERT(2,IEL)
      IV3=KVERT(3,IEL)
      IV4=KVERT(4,IEL)
      IV5=KVERT(5,IEL)
      IV6=KVERT(6,IEL)
      IV7=KVERT(7,IEL)
      IV8=KVERT(8,IEL)
C
      DPL(IV1)=DPL(IV1)+0.125D0*DAVOL*DPIEL
      DPL(IV2)=DPL(IV2)+0.125D0*DAVOL*DPIEL
      DPL(IV3)=DPL(IV3)+0.125D0*DAVOL*DPIEL
      DPL(IV4)=DPL(IV4)+0.125D0*DAVOL*DPIEL
      DPL(IV5)=DPL(IV5)+0.125D0*DAVOL*DPIEL
      DPL(IV6)=DPL(IV6)+0.125D0*DAVOL*DPIEL
      DPL(IV7)=DPL(IV7)+0.125D0*DAVOL*DPIEL
      DPL(IV8)=DPL(IV8)+0.125D0*DAVOL*DPIEL
C
      DAUX(IV1)=DAUX(IV1)+0.125D0*DAVOL
      DAUX(IV2)=DAUX(IV2)+0.125D0*DAVOL
      DAUX(IV3)=DAUX(IV3)+0.125D0*DAVOL
      DAUX(IV4)=DAUX(IV4)+0.125D0*DAVOL
      DAUX(IV5)=DAUX(IV5)+0.125D0*DAVOL
      DAUX(IV6)=DAUX(IV6)+0.125D0*DAVOL
      DAUX(IV7)=DAUX(IV7)+0.125D0*DAVOL
      DAUX(IV8)=DAUX(IV8)+0.125D0*DAVOL
C
10    CONTINUE
C
C
      DO 20 IVT=1,NVT
20    DPL(IVT)=DPL(IVT)/DAUX(IVT)
C
c      do 335 II=1,NVT
c      INPR=KNPR(II)
c      IF (INPR.EQ.11) THEN
c      DPL(II)=0D0
c      ENDIF
c335   CONTINUE
C
C
C      VPH=0E0
C      DO 30 IVT=1,NVT
C30    VPH=VPH+VPL(IVT)
C      VMWP=VPH/REAL(NVT)
C
C
C      DO 40 IVT=1,NVT
C40    VPL(IVT)=VPL(IVT)-VMWP
C
C
      END
c
cc
c
************************************************************************
      SUBROUTINE   SETARE  (AVOL,NEL,KVERT,DCORVG)
************************************************************************
*
*   Purpose: - writes on  AVOL(IEL)  the VOLUME of the element IEL,
*              IEL=1,...,NEL
*            - writes on  AVOL(NEL+1) the sum of all  AVOL(IEL)
*            - KVERT,DCORVG are the usual FEAT arrays
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      REAL AVOL
C
      PARAMETER (NNVE=8)
      PARAMETER (A1=1D0/6D0)
C
      DIMENSION  AVOL(*),KVERT(NNVE,*),DCORVG(3,*)
C=======================================================================
      SUM=0.D0
      DO  11  IEL=1,NEL
C
      I1=KVERT(1,IEL)
      I2=KVERT(2,IEL)
      I3=KVERT(3,IEL)
      I4=KVERT(4,IEL)
      I5=KVERT(5,IEL)
      I6=KVERT(6,IEL)
      I7=KVERT(7,IEL)
      I8=KVERT(8,IEL)
C
      X1=DCORVG(1,I1)
      X2=DCORVG(1,I2)
      X3=DCORVG(1,I3)
      X4=DCORVG(1,I4)
      X5=DCORVG(1,I5)
      X6=DCORVG(1,I6)
      X7=DCORVG(1,I7)
      X8=DCORVG(1,I8)
C
      Y1=DCORVG(2,I1)
      Y2=DCORVG(2,I2)
      Y3=DCORVG(2,I3)
      Y4=DCORVG(2,I4)
      Y5=DCORVG(2,I5)
      Y6=DCORVG(2,I6)
      Y7=DCORVG(2,I7)
      Y8=DCORVG(2,I8)
C
      Z1=DCORVG(3,I1)
      Z2=DCORVG(3,I2)
      Z3=DCORVG(3,I3)
      Z4=DCORVG(3,I4)
      Z5=DCORVG(3,I5)
      Z6=DCORVG(3,I6)
      Z7=DCORVG(3,I7)
      Z8=DCORVG(3,I8)
C
      AAA=A1*((DABS((X4-X1)*(Y4-Y3)*(Z4-Z8)+(Y4-Y1)*
     *       (Z4-Z3)*(X4-X8)+(Z4-Z1)*(X4-X3)*(Y4-Y8)-
     *       (X4-X8)*(Y4-Y3)*(Z4-Z1)-(Y4-Y8)*(Z4-Z3)*
     *       (X4-X1)-(Z4-Z8)*(X4-X3)*(Y4-Y1)))+
     *       (DABS((X2-X3)*(Y2-Y1)*(Z2-Z6)+(Y2-Y3)*
     *       (Z2-Z1)*(X2-X6)+(Z2-Z3)*(X2-X1)*(Y2-Y6)-
     *       (X2-X6)*(Y2-Y1)*(Z2-Z3)-(Y2-Y6)*(Z2-Z1)*
     *       (X2-X3)-(Z2-Z6)*(X2-X1)*(Y2-Y3)))+
     *       (DABS((X5-X8)*(Y5-Y6)*(Z5-Z1)+(Y5-Y8)*
     *       (Z5-Z6)*(X5-X1)+(Z5-Z8)*(X5-X6)*(Y5-Y1)-
     *       (X5-X1)*(Y5-Y6)*(Z5-Z8)-(Y5-Y1)*(Z5-Z6)*
     *       (X5-X8)-(Z5-Z1)*(X5-X6)*(Y5-Y8)))+
     *       (DABS((X7-X6)*(Y7-Y8)*(Z7-Z3)+(Y7-Y6)*
     *       (Z7-Z8)*(X7-X3)+(Z7-Z6)*(X7-X8)*(Y7-Y3)-
     *       (X7-X3)*(Y7-Y8)*(Z7-Z6)-(Y7-Y3)*(Z7-Z8)*
     *       (X7-X6)-(Z7-Z3)*(X7-X8)*(Y7-Y6)))+
     *       (DABS((X1-X3)*(Y1-Y8)*(Z1-Z6)+(Y1-Y3)*
     *       (Z1-Z8)*(X1-X6)+(Z1-Z3)*(X1-X8)*(Y1-Y6)-
     *       (X1-X6)*(Y1-Y8)*(Z1-Z3)-(Y1-Y6)*(Z1-Z8)*
     *       (X1-X3)-(Z1-Z6)*(X1-X8)*(Y1-Y3))))
      AVOL(IEL)=REAL(AAA)
      SUM=SUM+AAA
  11  CONTINUE
C
      AVOL(NEL+1)=REAL(SUM)
C
      END
c
c
c
************************************************************************
      SUBROUTINE    TOL20A  (P,AVOL, NEL,INEUM)
************************************************************************
*
*    Purpose: - Transforms the vector P into the space L2_0
*             - uses the vector AVOL with the VOLUME of all elements
*               and on AVOL(NEL+1) the sum of all AVOL(IEL)
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      REAL AVOL
C
      DIMENSION P(*),AVOL(*)
C
      IF (INEUM.EQ.1) RETURN
C
      PINT=0.D0
      DO 2  IEL=1,NEL
  2   PINT=PINT+P(IEL)*DBLE(AVOL(IEL))
C
      C=PINT/DBLE(AVOL(NEL+1))
C
      DO 3  IEL=1,NEL
  3   P(IEL)=P(IEL)-C
C
      END
c
c
c
c
************************************************************************
      SUBROUTINE   SETLEV (ISETLV)
************************************************************************
*
*   Purpose:  sets all data for current level ILEV (from /MGPAR/)
*            
*   Input:
*   -------
*     ILEV        - current level number (from /MGPAR/)
*     ISETLV >=1  - update of /TRIAA/,TRIAD/ 
*            >=2  - update of /LEVDIM/,/ADRFLD/
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299,NNLEV=9,  NNWORK=1)
C
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C *** Standard COMMON blocks
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
C
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
C
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** COMMON blocks for multigrid data management
C
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
C
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNET(NNLEV),
     *                KNAT(NNLEV),KNVE(NNLEV),KNEE(NNLEV),
     *                KNAE(NNLEV),KNVEL(NNLEV),KNEEL(NNLEV),
     *                KNVED(NNLEV),KNVAR(NNLEV),KNEAR(NNLEV),
     *                KNBCT(NNLEV),KNVBD(NNLEV),KNEBD(NNLEV),
     *                KNABD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLCAG(NNLEV),
     *                KLVERT(NNLEV),KLEDGE(NNLEV),KLAREA(NNLEV),
     *                KLADJ(NNLEV),KLVEL(NNLEV),KLEEL(NNLEV),
     *                KLAEL(NNLEV),KLVED(NNLEV),KLAED(NNLEV),
     *                KLVAR(NNLEV),KLEAR(NNLEV),KLEVE(NNLEV),
     *                KLAVE(NNLEV),KLNPR(NNLEV),KLBCT(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLABD(NNLEV)
C
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
C
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,  NU,NP,  NUP
C
      COMMON /ADRFLD/ KA1,KST1,
     *                KM1,KCOLA,KLDA,KB1,KB2,KB3,KCOLB,KLDB,
     *                KU1,KU2,KU3,KP,KF1,KF2,KF3,KFP,KAUX1,KAUX2,
     *                KAUX3,KAUXP
C
C
      SAVE 
C-----------------------------------------------------------------------
C *** elementary check
      IF (ILEV.LT.NLMIN .OR. ILEV.GT.NLMAX .OR. ISETLV.LT.1  .OR.
     *    ISETLV.GT.2 )  THEN
        WRITE(MTERM,*) 'ERROR in SETLEV: ILEV or ISETLV is wrong'
        STOP
      ENDIF
C
C *** update of /TRIAD/,/TRIAA/
C
      NEL =KNEL(ILEV)    
      NVT =KNVT(ILEV)
      NET =KNET(ILEV)
      NAT =KNAT(ILEV)
      NVE =KNVE(ILEV)
      NEE =KNEE(ILEV)
      NAE =KNAE(ILEV)
      NVEL =KNVEL(ILEV)
      NEEL =KNEEL(ILEV)
      NVED =KNVED(ILEV)
      NVAR =KNVAR(ILEV)
      NEAR =KNEAR(ILEV)
      NBCT =KNBCT(ILEV)
      NVBD =KNVBD(ILEV)
      NEBD =KNEBD(ILEV)
      NABD =KNABD(ILEV)
C
      LCORVG =KLCVG(ILEV) 
      LCORMG =KLCMG(ILEV) 
      LCORAG =KLCAG(ILEV) 
      LVERT  =KLVERT(ILEV)
      LEDGE  =KLEDGE(ILEV) 
      LAREA  =KLAREA(ILEV) 
      LADJ   =KLADJ(ILEV) 
      LVEL   =KLVEL(ILEV) 
      LEEL   =KLEEL(ILEV)
      LAEL   =KLAEL(ILEV)   
      LVED   =KLVED(ILEV) 
      LAED   =KLAED(ILEV) 
      LVAR   =KLVAR(ILEV) 
      LEAR   =KLEAR(ILEV) 
      LEVE   =KLEVE(ILEV) 
      LAVE   =KLAVE(ILEV) 
      LNPR   =KLNPR(ILEV) 
      LBCT   =KLBCT(ILEV) 
      LVBD   =KLVBD(ILEV) 
      LEBD   =KLEBD(ILEV) 
      LABD   =KLABD(ILEV)
C
C *** update of /LEVDIM/,/ADRFLD/ if  ISETLV=2
C
      IF (ISETLV.EQ.2)  THEN
C
         NA =KNA  (ILEV)
         NB =KNB  (ILEV)
         NU =KNU  (ILEV)
         NP =KNP  (ILEV)
         NUP=KNUP (ILEV)
C
         KA1  =L(KLA    (ILEV))
         KST1 =L(KLST   (ILEV))
         KM1  =L(KLM    (ILEV))
         KCOLA=L(KLCOLA (ILEV))
         KLDA =L(KLLDA  (ILEV))
         KB1  =L(KLB1   (ILEV))
         KB2  =L(KLB2   (ILEV))
         KB3  =L(KLB3   (ILEV))
         KCOLB=L(KLCOLB (ILEV))
         KLDB =L(KLLDB  (ILEV))
         KU1  =L(KLUP   (ILEV))
         KU2  =KU1+NU
         KU3  =KU2+NU
         KP   =KU3+NU
         KF1  =L(KLF12P (ILEV))
         KF2  =KF1+NU
         KF3  =KF2+NU
         KFP  =KF3+NU
         KAUX1=L(KLAUX  (ILEV))
         KAUX2=KAUX1+NU
         KAUX3=KAUX2+NU
         KAUXP=KAUX3+NU
C
      ENDIF
C
      END
c
c
************************************************************************
      SUBROUTINE CRITAD(TIMEIN,TIMENS,TIMEST,EPSADI,EPSADL,EPSAD,IADIN)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (PI=3.1415926535897931D0)
C
      EPSAD=EPSADL
C
      IF (TIMEIN.GT.0) THEN
       TDIFF=TIMENS-TIMEST
C
       IF (IADIN.EQ.0) THEN
        IF (TDIFF.LE.TIMEIN) THEN
         EPSAD=EPSADI
        ELSE
         EPSAD=EPSADL
        ENDIF
       ENDIF
C
       IF (IADIN.EQ.1) THEN
        IF (TDIFF.LE.TIMEIN) THEN
         EPSAD=EPSADI+TDIFF/TIMEIN*(EPSADL-EPSADI)
        ELSE
         EPSAD=EPSADL
        ENDIF
       ENDIF
C
       IF (IADIN.EQ.2) THEN
        IF (TDIFF.LE.TIMEIN) THEN
         EPSAD=EPSADI**(1D0-TDIFF/TIMEIN)*EPSADL**(TDIFF/TIMEIN)
        ELSE
         EPSAD=EPSADL
        ENDIF
       ENDIF
C
      ENDIF
C
C
C     
      END
C
C
C
************************************************************************
      SUBROUTINE C2N2DM (DPC,DPL,KAREA,KADJ,NEL,NAT,NVT,IPAR)
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION DPC(*),DPL(*),KAREA(NNAE,*),KADJ(NNAE,*)
C
C-----------------------------------------------------------------------
C
C *** constant to nonconforming = 0
      IF (IPAR.EQ.0) THEN
C
       DO 10 IEL=1,NEL
       DPH  =DPC(IEL)
C
       DO 20 IVE=1,6
       IADJ=KADJ(IVE,IEL)
       IAREA=KAREA(IVE,IEL)
C
       IF (IADJ.EQ.0)   DPL(IAREA)=DPH
       IF (IADJ.GT.IEL) DPL(IAREA)=0.5D0*(DPH+DPC(IADJ))
C
20     CONTINUE
10     CONTINUE
C
      ELSE
C
       DO 110 IEL=1,NEL
       DPC(IEL)=0.16666666D0*( DPL(KAREA(1,IEL))+
     *                  DPL(KAREA(2,IEL))
     *                  +DPL(KAREA(3,IEL))+DPL(KAREA(4,IEL))
     *                  +DPL(KAREA(5,IEL))+DPL(KAREA(6,IEL)))
110    CONTINUE
C
      ENDIF      
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE EM30(XI1,XI2,XI3,IPAR)
************************************************************************
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNAE=6)
      PARAMETER (NNDIM=3,NNCOF=10)
      PARAMETER (Q2=1D0/3D0)
      PARAMETER (Q8=-0.5625D0,Q9=0.5208333D0)
      DIMENSION DXM(6),DYM(6),DZM(6),A(6,6),F(6),CKH(6),CK(6,6)
      DIMENSION P1X(6),P1Y(6),P1Z(6)
      DIMENSION P2X(6),P2Y(6),P2Z(6)
      DIMENSION P3X(6),P3Y(6),P3Z(6)
      DIMENSION P4X(6),P4Y(6),P4Z(6)
      DIMENSION P5X(6),P5Y(6),P5Z(6)
      DIMENSION P6X(6),P6Y(6),P6Z(6)
      DIMENSION P7X(6),P7Y(6),P7Z(6)
      DIMENSION P8X(6),P8Y(6),P8Z(6)
      DIMENSION PP1(6),PP2(6),PP3(6),PP4(6),PP5(6)
      DIMENSION PS1(6),PS2(6),PQ1(6),PQ2(6),PQ3(6)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,
     *                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),
     *                IEL,NDIM
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COFBAS/ COB(NNBAS,NNCOF)
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/,/CHAR/,/CUB/
C
      F1(X,Y,Z,CA1,CB1,CC1,CA2,CB2,CC2,CA3,CB3,CC3,
     *   CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)=1D0
      F2(X,Y,Z,CA1,CB1,CC1,CA2,CB2,CC2,CA3,CB3,CC3,
     *   CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)=
     *                   CA1*X  +CB1*Y +CC1*Z
      F3(X,Y,Z,CA1,CB1,CC1,CA2,CB2,CC2,CA3,CB3,CC3,
     *   CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)=
     *                   CA2*X  +CB2*Y +CC2*Z
      F4(X,Y,Z,CA1,CB1,CC1,CA2,CB2,CC2,CA3,CB3,CC3,
     *   CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)=
     *                   CA3*X  +CB3*Y +CC3*Z
      F5(X,Y,Z,CA1,CB1,CC1,CA2,CB2,CC2,CA3,CB3,CC3,
     *   CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)=
     *  CD1*X*X+CD2*Y*Y+CD3*Z*Z+CD4*X*Y+CD5*X*Z+CD6*Y*Z
      F6(X,Y,Z,CA1,CB1,CC1,CA2,CB2,CC2,CA3,CB3,CC3,
     *   CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)=
     *  CE1*X*X+CE2*Y*Y+CE3*Z*Z+CE4*X*Y+CE5*X*Z+CE6*Y*Z
C
C
      SUB='EM30'
      IF (ICHECK.GE.998) CALL OTRC('EM30  ','27/02/95')
C
C
C *** Dummy call
      IF (IPAR.EQ.-1) THEN
       IER=0
       IPAR=30
       GOTO 99999
      ENDIF
C
C
      IF (IPAR.EQ.-2) THEN
       DXM(1)=0.25D0*(DX(1)+DX(2)+DX(3)+DX(4))
       DYM(1)=0.25D0*(DY(1)+DY(2)+DY(3)+DY(4))
       DZM(1)=0.25D0*(DZ(1)+DZ(2)+DZ(3)+DZ(4))
C
       DXM(2)=0.25D0*(DX(1)+DX(2)+DX(6)+DX(5))
       DYM(2)=0.25D0*(DY(1)+DY(2)+DY(6)+DY(5))
       DZM(2)=0.25D0*(DZ(1)+DZ(2)+DZ(6)+DZ(5))
C
       DXM(3)=0.25D0*(DX(2)+DX(3)+DX(7)+DX(6))
       DYM(3)=0.25D0*(DY(2)+DY(3)+DY(7)+DY(6))
       DZM(3)=0.25D0*(DZ(2)+DZ(3)+DZ(7)+DZ(6))
C
       DXM(4)=0.25D0*(DX(3)+DX(7)+DX(8)+DX(4))
       DYM(4)=0.25D0*(DY(3)+DY(7)+DY(8)+DY(4))
       DZM(4)=0.25D0*(DZ(3)+DZ(7)+DZ(8)+DZ(4))
C
       DXM(5)=0.25D0*(DX(1)+DX(4)+DX(8)+DX(5))
       DYM(5)=0.25D0*(DY(1)+DY(4)+DY(8)+DY(5))
       DZM(5)=0.25D0*(DZ(1)+DZ(4)+DZ(8)+DZ(5))
C
       DXM(6)=0.25D0*(DX(5)+DX(6)+DX(7)+DX(8))
       DYM(6)=0.25D0*(DY(5)+DY(6)+DY(7)+DY(8))
       DZM(6)=0.25D0*(DZ(5)+DZ(6)+DZ(7)+DZ(8))
C
C
       CA1=(DXM(1)-DXM(6))/SQRT((DXM(1)-DXM(6))**2+(DYM(1)-DYM(6))**2+
     *     (DZM(1)-DZM(6))**2)
       CB1=(DYM(1)-DYM(6))/SQRT((DXM(1)-DXM(6))**2+(DYM(1)-DYM(6))**2+
     *     (DZM(1)-DZM(6))**2)
       CC1=(DZM(1)-DZM(6))/SQRT((DXM(1)-DXM(6))**2+(DYM(1)-DYM(6))**2+
     *     (DZM(1)-DZM(6))**2)
       CA2=(DXM(2)-DXM(4))/SQRT((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2+
     *     (DZM(2)-DZM(4))**2)
       CB2=(DYM(2)-DYM(4))/SQRT((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2+
     *     (DZM(2)-DZM(4))**2)
       CC2=(DZM(2)-DZM(4))/SQRT((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2+
     *     (DZM(2)-DZM(4))**2)
       CA3=(DXM(3)-DXM(5))/SQRT((DXM(3)-DXM(5))**2+(DYM(3)-DYM(5))**2+
     *     (DZM(3)-DZM(5))**2)
       CB3=(DYM(3)-DYM(5))/SQRT((DXM(3)-DXM(5))**2+(DYM(3)-DYM(5))**2+
     *     (DZM(3)-DZM(5))**2)
       CC3=(DZM(3)-DZM(5))/SQRT((DXM(3)-DXM(5))**2+(DYM(3)-DYM(5))**2+
     *     (DZM(3)-DZM(5))**2)
       CD1=CA1**2-CA2**2
       CD2=CB1**2-CB2**2
       CD3=CC1**2-CC2**2
       CD4=2D0*CA1*CB1-2D0*CA2*CB2
       CD5=2D0*CA1*CC1-2D0*CA2*CC2
       CD6=2D0*CB1*CC1-2D0*CB2*CC2
       CE1=CA2**2-CA3**2
       CE2=CB2**2-CB3**2
       CE3=CC2**2-CC3**2
       CE4=2D0*CA2*CB2-2D0*CA3*CB3
       CE5=2D0*CA2*CC2-2D0*CA3*CC3
       CE6=2D0*CB2*CC2-2D0*CB3*CC3
C
       P1X(1)=Q2*(DX(1)+DX(2)+DX(3))
       P1Y(1)=Q2*(DY(1)+DY(2)+DY(3))
       P1Z(1)=Q2*(DZ(1)+DZ(2)+DZ(3))
       P2X(1)=0.6D0*DX(1)+0.2D0*DX(2)+0.2D0*DX(3)
       P2Y(1)=0.6D0*DY(1)+0.2D0*DY(2)+0.2D0*DY(3)
       P2Z(1)=0.6D0*DZ(1)+0.2D0*DZ(2)+0.2D0*DZ(3)
       P3X(1)=0.2D0*DX(1)+0.6D0*DX(2)+0.2D0*DX(3)
       P3Y(1)=0.2D0*DY(1)+0.6D0*DY(2)+0.2D0*DY(3)
       P3Z(1)=0.2D0*DZ(1)+0.6D0*DZ(2)+0.2D0*DZ(3)
       P4X(1)=0.2D0*DX(1)+0.2D0*DX(2)+0.6D0*DX(3)
       P4Y(1)=0.2D0*DY(1)+0.2D0*DY(2)+0.6D0*DY(3)
       P4Z(1)=0.2D0*DZ(1)+0.2D0*DZ(2)+0.6D0*DZ(3)
       P5X(1)=Q2*(DX(1)+DX(3)+DX(4))
       P5Y(1)=Q2*(DY(1)+DY(3)+DY(4))
       P5Z(1)=Q2*(DZ(1)+DZ(3)+DZ(4))
       P6X(1)=0.6D0*DX(1)+0.2D0*DX(3)+0.2D0*DX(4)
       P6Y(1)=0.6D0*DY(1)+0.2D0*DY(3)+0.2D0*DY(4)
       P6Z(1)=0.6D0*DZ(1)+0.2D0*DZ(3)+0.2D0*DZ(4)
       P7X(1)=0.2D0*DX(1)+0.6D0*DX(3)+0.2D0*DX(4)
       P7Y(1)=0.2D0*DY(1)+0.6D0*DY(3)+0.2D0*DY(4)
       P7Z(1)=0.2D0*DZ(1)+0.6D0*DZ(3)+0.2D0*DZ(4)
       P8X(1)=0.2D0*DX(1)+0.2D0*DX(3)+0.6D0*DX(4)
       P8Y(1)=0.2D0*DY(1)+0.2D0*DY(3)+0.6D0*DY(4)
       P8Z(1)=0.2D0*DZ(1)+0.2D0*DZ(3)+0.6D0*DZ(4)
C ***  das sind die 8 Kubaturpunkte auf der ersten flaeche
       PP1(1)=DSQRT((DX(1)-DX(2))**2+(DY(1)-DY(2))**2+
     *       (DZ(1)-DZ(2))**2)
       PP2(1)=DSQRT((DX(2)-DX(3))**2+(DY(2)-DY(3))**2+
     *       (DZ(2)-DZ(3))**2)
       PP3(1)=DSQRT((DX(1)-DX(3))**2+(DY(1)-DY(3))**2+
     *       (DZ(1)-DZ(3))**2)
       PP4(1)=DSQRT((DX(3)-DX(4))**2+(DY(3)-DY(4))**2+
     *       (DZ(3)-DZ(4))**2)
       PP5(1)=DSQRT((DX(1)-DX(4))**2+(DY(1)-DY(4))**2+
     *       (DZ(1)-DZ(4))**2)
       PS1(1)=(PP1(1)+PP2(1)+PP3(1))/2D0
       PS2(1)=(PP3(1)+PP4(1)+PP5(1))/2D0
       PQ1(1)=DSQRT(PS1(1)*(PS1(1)-PP1(1))*(PS1(1)-
     *       PP2(1))*(PS1(1)-PP3(1)))
       PQ2(1)=DSQRT(PS2(1)*(PS2(1)-PP3(1))*(PS2(1)-
     *       PP4(1))*(PS2(1)-PP5(1)))
       PQ3(1)=PQ1(1)+PQ2(1)
C
       P1X(2)=Q2*(DX(1)+DX(2)+DX(6))
       P1Y(2)=Q2*(DY(1)+DY(2)+DY(6))
       P1Z(2)=Q2*(DZ(1)+DZ(2)+DZ(6))
       P2X(2)=0.6D0*DX(1)+0.2D0*DX(2)+0.2D0*DX(6)
       P2Y(2)=0.6D0*DY(1)+0.2D0*DY(2)+0.2D0*DY(6)
       P2Z(2)=0.6D0*DZ(1)+0.2D0*DZ(2)+0.2D0*DZ(6)
       P3X(2)=0.2D0*DX(1)+0.6D0*DX(2)+0.2D0*DX(6)
       P3Y(2)=0.2D0*DY(1)+0.6D0*DY(2)+0.2D0*DY(6)
       P3Z(2)=0.2D0*DZ(1)+0.6D0*DZ(2)+0.2D0*DZ(6)
       P4X(2)=0.2D0*DX(1)+0.2D0*DX(2)+0.6D0*DX(6)
       P4Y(2)=0.2D0*DY(1)+0.2D0*DY(2)+0.6D0*DY(6)
       P4Z(2)=0.2D0*DZ(1)+0.2D0*DZ(2)+0.6D0*DZ(6)
       P5X(2)=Q2*(DX(1)+DX(6)+DX(5))
       P5Y(2)=Q2*(DY(1)+DY(6)+DY(5))
       P5Z(2)=Q2*(DZ(1)+DZ(6)+DZ(5))
       P6X(2)=0.6D0*DX(1)+0.2D0*DX(6)+0.2D0*DX(5)
       P6Y(2)=0.6D0*DY(1)+0.2D0*DY(6)+0.2D0*DY(5)
       P6Z(2)=0.6D0*DZ(1)+0.2D0*DZ(6)+0.2D0*DZ(5)
       P7X(2)=0.2D0*DX(1)+0.6D0*DX(6)+0.2D0*DX(5)
       P7Y(2)=0.2D0*DY(1)+0.6D0*DY(6)+0.2D0*DY(5)
       P7Z(2)=0.2D0*DZ(1)+0.6D0*DZ(6)+0.2D0*DZ(5)
       P8X(2)=0.2D0*DX(1)+0.2D0*DX(6)+0.6D0*DX(5)
       P8Y(2)=0.2D0*DY(1)+0.2D0*DY(6)+0.6D0*DY(5)
       P8Z(2)=0.2D0*DZ(1)+0.2D0*DZ(6)+0.6D0*DZ(5)
C ***  das sind die 8 Kubaturpunkte auf der zweiten flaeche
       PP1(2)=DSQRT((DX(1)-DX(2))**2+(DY(1)-DY(2))**2+
     *       (DZ(1)-DZ(2))**2)
       PP2(2)=DSQRT((DX(2)-DX(6))**2+(DY(2)-DY(6))**2+
     *       (DZ(2)-DZ(6))**2)
       PP3(2)=DSQRT((DX(1)-DX(6))**2+(DY(1)-DY(6))**2+
     *       (DZ(1)-DZ(6))**2)
       PP4(2)=DSQRT((DX(5)-DX(6))**2+(DY(5)-DY(6))**2+
     *       (DZ(5)-DZ(6))**2)
       PP5(2)=DSQRT((DX(1)-DX(5))**2+(DY(1)-DY(5))**2+
     *       (DZ(1)-DZ(5))**2)
       PS1(2)=(PP1(2)+PP2(2)+PP3(2))/2D0
       PS2(2)=(PP3(2)+PP4(2)+PP5(2))/2D0
       PQ1(2)=DSQRT(PS1(2)*(PS1(2)-PP1(2))*(PS1(2)-
     *       PP2(2))*(PS1(2)-PP3(2)))
       PQ2(2)=DSQRT(PS2(2)*(PS2(2)-PP3(2))*(PS2(2)-
     *       PP4(2))*(PS2(2)-PP5(2)))
       PQ3(2)=PQ1(2)+PQ2(2)
C
C
       P1X(3)=Q2*(DX(2)+DX(3)+DX(7))
       P1Y(3)=Q2*(DY(2)+DY(3)+DY(7))
       P1Z(3)=Q2*(DZ(2)+DZ(3)+DZ(7))
       P2X(3)=0.6D0*DX(2)+0.2D0*DX(3)+0.2D0*DX(7)
       P2Y(3)=0.6D0*DY(2)+0.2D0*DY(3)+0.2D0*DY(7)
       P2Z(3)=0.6D0*DZ(2)+0.2D0*DZ(3)+0.2D0*DZ(7)
       P3X(3)=0.2D0*DX(2)+0.6D0*DX(3)+0.2D0*DX(7)
       P3Y(3)=0.2D0*DY(2)+0.6D0*DY(3)+0.2D0*DY(7)
       P3Z(3)=0.2D0*DZ(2)+0.6D0*DZ(3)+0.2D0*DZ(7)
       P4X(3)=0.2D0*DX(2)+0.2D0*DX(3)+0.6D0*DX(7)
       P4Y(3)=0.2D0*DY(2)+0.2D0*DY(3)+0.6D0*DY(7)
       P4Z(3)=0.2D0*DZ(2)+0.2D0*DZ(3)+0.6D0*DZ(7)
       P5X(3)=Q2*(DX(2)+DX(7)+DX(6))
       P5Y(3)=Q2*(DY(2)+DY(7)+DY(6))
       P5Z(3)=Q2*(DZ(2)+DZ(7)+DZ(6))
       P6X(3)=0.6D0*DX(2)+0.2D0*DX(7)+0.2D0*DX(6)
       P6Y(3)=0.6D0*DY(2)+0.2D0*DY(7)+0.2D0*DY(6)
       P6Z(3)=0.6D0*DZ(2)+0.2D0*DZ(7)+0.2D0*DZ(6)
       P7X(3)=0.2D0*DX(2)+0.6D0*DX(7)+0.2D0*DX(6)
       P7Y(3)=0.2D0*DY(2)+0.6D0*DY(7)+0.2D0*DY(6)
       P7Z(3)=0.2D0*DZ(2)+0.6D0*DZ(7)+0.2D0*DZ(6)
       P8X(3)=0.2D0*DX(2)+0.2D0*DX(7)+0.6D0*DX(6)
       P8Y(3)=0.2D0*DY(2)+0.2D0*DY(7)+0.6D0*DY(6)
       P8Z(3)=0.2D0*DZ(2)+0.2D0*DZ(7)+0.6D0*DZ(6)
C ***  das sind die 8 Kubaturpunkte auf der dritten flaeche
       PP1(3)=DSQRT((DX(2)-DX(3))**2+(DY(2)-DY(3))**2+
     *       (DZ(2)-DZ(3))**2)
       PP2(3)=DSQRT((DX(3)-DX(7))**2+(DY(3)-DY(7))**2+
     *       (DZ(3)-DZ(7))**2)
       PP3(3)=DSQRT((DX(2)-DX(7))**2+(DY(2)-DY(7))**2+
     *       (DZ(2)-DZ(7))**2)
       PP4(3)=DSQRT((DX(7)-DX(6))**2+(DY(7)-DY(6))**2+
     *       (DZ(7)-DZ(6))**2)
       PP5(3)=DSQRT((DX(2)-DX(6))**2+(DY(2)-DY(6))**2+
     *       (DZ(2)-DZ(6))**2)
       PS1(3)=(PP1(3)+PP2(3)+PP3(3))/2D0
       PS2(3)=(PP3(3)+PP4(3)+PP5(3))/2D0
       PQ1(3)=DSQRT(PS1(3)*(PS1(3)-PP1(3))*(PS1(3)-
     *       PP2(3))*(PS1(3)-PP3(3)))
       PQ2(3)=DSQRT(PS2(3)*(PS2(3)-PP3(3))*(PS2(3)-
     *       PP4(3))*(PS2(3)-PP5(3)))
       PQ3(3)=PQ1(3)+PQ2(3)
C
       P1X(4)=Q2*(DX(3)+DX(4)+DX(8))
       P1Y(4)=Q2*(DY(3)+DY(4)+DY(8))
       P1Z(4)=Q2*(DZ(3)+DZ(4)+DZ(8))
       P2X(4)=0.6D0*DX(3)+0.2D0*DX(4)+0.2D0*DX(8)
       P2Y(4)=0.6D0*DY(3)+0.2D0*DY(4)+0.2D0*DY(8)
       P2Z(4)=0.6D0*DZ(3)+0.2D0*DZ(4)+0.2D0*DZ(8)
       P3X(4)=0.2D0*DX(3)+0.6D0*DX(4)+0.2D0*DX(8)
       P3Y(4)=0.2D0*DY(3)+0.6D0*DY(4)+0.2D0*DY(8)
       P3Z(4)=0.2D0*DZ(3)+0.6D0*DZ(4)+0.2D0*DZ(8)
       P4X(4)=0.2D0*DX(3)+0.2D0*DX(4)+0.6D0*DX(8)
       P4Y(4)=0.2D0*DY(3)+0.2D0*DY(4)+0.6D0*DY(8)
       P4Z(4)=0.2D0*DZ(3)+0.2D0*DZ(4)+0.6D0*DZ(8)
       P5X(4)=Q2*(DX(3)+DX(8)+DX(7))
       P5Y(4)=Q2*(DY(3)+DY(8)+DY(7))
       P5Z(4)=Q2*(DZ(3)+DZ(8)+DZ(7))
       P6X(4)=0.6D0*DX(3)+0.2D0*DX(8)+0.2D0*DX(7)
       P6Y(4)=0.6D0*DY(3)+0.2D0*DY(8)+0.2D0*DY(7)
       P6Z(4)=0.6D0*DZ(3)+0.2D0*DZ(8)+0.2D0*DZ(7)
       P7X(4)=0.2D0*DX(3)+0.6D0*DX(8)+0.2D0*DX(7)
       P7Y(4)=0.2D0*DY(3)+0.6D0*DY(8)+0.2D0*DY(7)
       P7Z(4)=0.2D0*DZ(3)+0.6D0*DZ(8)+0.2D0*DZ(7)
       P8X(4)=0.2D0*DX(3)+0.2D0*DX(8)+0.6D0*DX(7)
       P8Y(4)=0.2D0*DY(3)+0.2D0*DY(8)+0.6D0*DY(7)
       P8Z(4)=0.2D0*DZ(3)+0.2D0*DZ(8)+0.6D0*DZ(7)
C ***  das sind die 8 Kubaturpunkte auf der vierten flaeche
       PP1(4)=DSQRT((DX(3)-DX(4))**2+(DY(3)-DY(4))**2+
     *       (DZ(3)-DZ(4))**2)
       PP2(4)=DSQRT((DX(4)-DX(8))**2+(DY(4)-DY(8))**2+
     *       (DZ(4)-DZ(8))**2)
       PP3(4)=DSQRT((DX(3)-DX(8))**2+(DY(3)-DY(8))**2+
     *       (DZ(3)-DZ(8))**2)
       PP4(4)=DSQRT((DX(7)-DX(8))**2+(DY(7)-DY(8))**2+
     *       (DZ(7)-DZ(8))**2)
       PP5(4)=DSQRT((DX(3)-DX(7))**2+(DY(3)-DY(7))**2+
     *       (DZ(3)-DZ(7))**2)
       PS1(4)=(PP1(4)+PP2(4)+PP3(4))/2D0
       PS2(4)=(PP3(4)+PP4(4)+PP5(4))/2D0
       PQ1(4)=DSQRT(PS1(4)*(PS1(4)-PP1(4))*(PS1(4)-
     *       PP2(4))*(PS1(4)-PP3(4)))
       PQ2(4)=DSQRT(PS2(4)*(PS2(4)-PP3(4))*(PS2(4)-
     *       PP4(4))*(PS2(4)-PP5(4)))
       PQ3(4)=PQ1(4)+PQ2(4)
C
       P1X(5)=Q2*(DX(1)+DX(4)+DX(8))
       P1Y(5)=Q2*(DY(1)+DY(4)+DY(8))
       P1Z(5)=Q2*(DZ(1)+DZ(4)+DZ(8))
       P2X(5)=0.6D0*DX(1)+0.2D0*DX(4)+0.2D0*DX(8)
       P2Y(5)=0.6D0*DY(1)+0.2D0*DY(4)+0.2D0*DY(8)
       P2Z(5)=0.6D0*DZ(1)+0.2D0*DZ(4)+0.2D0*DZ(8)
       P3X(5)=0.2D0*DX(1)+0.6D0*DX(4)+0.2D0*DX(8)
       P3Y(5)=0.2D0*DY(1)+0.6D0*DY(4)+0.2D0*DY(8)
       P3Z(5)=0.2D0*DZ(1)+0.6D0*DZ(4)+0.2D0*DZ(8)
       P4X(5)=0.2D0*DX(1)+0.2D0*DX(4)+0.6D0*DX(8)
       P4Y(5)=0.2D0*DY(1)+0.2D0*DY(4)+0.6D0*DY(8)
       P4Z(5)=0.2D0*DZ(1)+0.2D0*DZ(4)+0.6D0*DZ(8)
       P5X(5)=Q2*(DX(1)+DX(5)+DX(8))
       P5Y(5)=Q2*(DY(1)+DY(5)+DY(8))
       P5Z(5)=Q2*(DZ(1)+DZ(5)+DZ(8))
       P6X(5)=0.6D0*DX(1)+0.2D0*DX(5)+0.2D0*DX(8)
       P6Y(5)=0.6D0*DY(1)+0.2D0*DY(5)+0.2D0*DY(8)
       P6Z(5)=0.6D0*DZ(1)+0.2D0*DZ(5)+0.2D0*DZ(8)
       P7X(5)=0.2D0*DX(1)+0.6D0*DX(5)+0.2D0*DX(8)
       P7Y(5)=0.2D0*DY(1)+0.6D0*DY(5)+0.2D0*DY(8)
       P7Z(5)=0.2D0*DZ(1)+0.6D0*DZ(5)+0.2D0*DZ(8)
       P8X(5)=0.2D0*DX(1)+0.2D0*DX(5)+0.6D0*DX(8)
       P8Y(5)=0.2D0*DY(1)+0.2D0*DY(5)+0.6D0*DY(8)
       P8Z(5)=0.2D0*DZ(1)+0.2D0*DZ(5)+0.6D0*DZ(8)
C ***  das sind die 8 Kubaturpunkte auf der fuenften flaeche
       PP1(5)=DSQRT((DX(1)-DX(4))**2+(DY(1)-DY(4))**2+
     *       (DZ(1)-DZ(4))**2)
       PP2(5)=DSQRT((DX(4)-DX(8))**2+(DY(4)-DY(8))**2+
     *       (DZ(4)-DZ(8))**2)
       PP3(5)=DSQRT((DX(1)-DX(8))**2+(DY(1)-DY(8))**2+
     *       (DZ(1)-DZ(8))**2)
       PP4(5)=DSQRT((DX(1)-DX(5))**2+(DY(1)-DY(5))**2+
     *       (DZ(1)-DZ(5))**2)
       PP5(5)=DSQRT((DX(5)-DX(8))**2+(DY(5)-DY(8))**2+
     *       (DZ(5)-DZ(8))**2)
       PS1(5)=(PP1(5)+PP2(5)+PP3(5))/2D0
       PS2(5)=(PP3(5)+PP4(5)+PP5(5))/2D0
       PQ1(5)=DSQRT(PS1(5)*(PS1(5)-PP1(5))*(PS1(5)-
     *       PP2(5))*(PS1(5)-PP3(5)))
       PQ2(5)=DSQRT(PS2(5)*(PS2(5)-PP3(5))*(PS2(5)-
     *       PP4(5))*(PS2(5)-PP5(5)))
       PQ3(5)=PQ1(5)+PQ2(5)
C
       P1X(6)=Q2*(DX(5)+DX(6)+DX(7))
       P1Y(6)=Q2*(DY(5)+DY(6)+DY(7))
       P1Z(6)=Q2*(DZ(5)+DZ(6)+DZ(7))
       P2X(6)=0.6D0*DX(5)+0.2D0*DX(6)+0.2D0*DX(7)
       P2Y(6)=0.6D0*DY(5)+0.2D0*DY(6)+0.2D0*DY(7)
       P2Z(6)=0.6D0*DZ(5)+0.2D0*DZ(6)+0.2D0*DZ(7)
       P3X(6)=0.2D0*DX(5)+0.6D0*DX(6)+0.2D0*DX(7)
       P3Y(6)=0.2D0*DY(5)+0.6D0*DY(6)+0.2D0*DY(7)
       P3Z(6)=0.2D0*DZ(5)+0.6D0*DZ(6)+0.2D0*DZ(7)
       P4X(6)=0.2D0*DX(5)+0.2D0*DX(6)+0.6D0*DX(7)
       P4Y(6)=0.2D0*DY(5)+0.2D0*DY(6)+0.6D0*DY(7)
       P4Z(6)=0.2D0*DZ(5)+0.2D0*DZ(6)+0.6D0*DZ(7)
       P5X(6)=Q2*(DX(5)+DX(7)+DX(8))
       P5Y(6)=Q2*(DY(5)+DY(7)+DY(8))
       P5Z(6)=Q2*(DZ(5)+DZ(7)+DZ(8))
       P6X(6)=0.6D0*DX(5)+0.2D0*DX(7)+0.2D0*DX(8)
       P6Y(6)=0.6D0*DY(5)+0.2D0*DY(7)+0.2D0*DY(8)
       P6Z(6)=0.6D0*DZ(5)+0.2D0*DZ(7)+0.2D0*DZ(8)
       P7X(6)=0.2D0*DX(5)+0.6D0*DX(7)+0.2D0*DX(8)
       P7Y(6)=0.2D0*DY(5)+0.6D0*DY(7)+0.2D0*DY(8)
       P7Z(6)=0.2D0*DZ(5)+0.6D0*DZ(7)+0.2D0*DZ(8)
       P8X(6)=0.2D0*DX(5)+0.2D0*DX(7)+0.6D0*DX(8)
       P8Y(6)=0.2D0*DY(5)+0.2D0*DY(7)+0.6D0*DY(8)
       P8Z(6)=0.2D0*DZ(5)+0.2D0*DZ(7)+0.6D0*DZ(8)
C ***  das sind die 8 Kubaturpunkte auf der sechsten flaeche
       PP1(6)=DSQRT((DX(5)-DX(6))**2+(DY(5)-DY(6))**2+
     *       (DZ(5)-DZ(6))**2)
       PP2(6)=DSQRT((DX(6)-DX(7))**2+(DY(6)-DY(7))**2+
     *       (DZ(6)-DZ(7))**2)
       PP3(6)=DSQRT((DX(5)-DX(7))**2+(DY(5)-DY(7))**2+
     *       (DZ(5)-DZ(7))**2)
       PP4(6)=DSQRT((DX(7)-DX(8))**2+(DY(7)-DY(8))**2+
     *       (DZ(7)-DZ(8))**2)
       PP5(6)=DSQRT((DX(5)-DX(8))**2+(DY(5)-DY(8))**2+
     *       (DZ(5)-DZ(8))**2)
       PS1(6)=(PP1(6)+PP2(6)+PP3(6))/2D0
       PS2(6)=(PP3(6)+PP4(6)+PP5(6))/2D0
       PQ1(6)=DSQRT(PS1(6)*(PS1(6)-PP1(6))*(PS1(6)-
     *       PP2(6))*(PS1(6)-PP3(6)))
       PQ2(6)=DSQRT(PS2(6)*(PS2(6)-PP3(6))*(PS2(6)-
     *       PP4(6))*(PS2(6)-PP5(6)))
       PQ3(6)=PQ1(6)+PQ2(6)
C
************************************************************************
C
       DO 22 IA=1,6
       A(IA,1)=(PQ1(IA)/PQ3(IA))*(Q8*F1(P1X(IA),P1Y(IA),P1Z(IA),
     * CA1,CB1,CC1,CA2,CB2,CC2,
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F1(P2X(IA),P2Y(IA),P2Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F1(P3X(IA),P3Y(IA),P3Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F1(P4X(IA),P4Y(IA),P4Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6))+
     * (PQ2(IA)/PQ3(IA))*(Q8*F1(P5X(IA),P5Y(IA),P5Z(IA),
     * CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F1(P6X(IA),P6Y(IA),P6Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F1(P7X(IA),P7Y(IA),P7Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F1(P8X(IA),P8Y(IA),P8Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6))
C
       A(IA,2)=(PQ1(IA)/PQ3(IA))*(Q8*F2(P1X(IA),P1Y(IA),P1Z(IA),
     * CA1,CB1,CC1,CA2,CB2,CC2,
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F2(P2X(IA),P2Y(IA),P2Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F2(P3X(IA),P3Y(IA),P3Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F2(P4X(IA),P4Y(IA),P4Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6))+
     * (PQ2(IA)/PQ3(IA))*(Q8*F2(P5X(IA),P5Y(IA),P5Z(IA),
     * CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F2(P6X(IA),P6Y(IA),P6Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F2(P7X(IA),P7Y(IA),P7Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F2(P8X(IA),P8Y(IA),P8Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6))
C
       A(IA,3)=(PQ1(IA)/PQ3(IA))*(Q8*F3(P1X(IA),P1Y(IA),P1Z(IA),
     * CA1,CB1,CC1,CA2,CB2,CC2,
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F3(P2X(IA),P2Y(IA),P2Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F3(P3X(IA),P3Y(IA),P3Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F3(P4X(IA),P4Y(IA),P4Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6))+
     * (PQ2(IA)/PQ3(IA))*(Q8*F3(P5X(IA),P5Y(IA),P5Z(IA),
     * CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F3(P6X(IA),P6Y(IA),P6Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F3(P7X(IA),P7Y(IA),P7Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F3(P8X(IA),P8Y(IA),P8Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6))
C
       A(IA,4)=(PQ1(IA)/PQ3(IA))*(Q8*F4(P1X(IA),P1Y(IA),P1Z(IA),
     * CA1,CB1,CC1,CA2,CB2,CC2,
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F4(P2X(IA),P2Y(IA),P2Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F4(P3X(IA),P3Y(IA),P3Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F4(P4X(IA),P4Y(IA),P4Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6))+
     * (PQ2(IA)/PQ3(IA))*(Q8*F4(P5X(IA),P5Y(IA),P5Z(IA),
     * CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F4(P6X(IA),P6Y(IA),P6Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F4(P7X(IA),P7Y(IA),P7Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F4(P8X(IA),P8Y(IA),P8Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6))
C
       A(IA,5)=(PQ1(IA)/PQ3(IA))*(Q8*F5(P1X(IA),P1Y(IA),P1Z(IA),
     * CA1,CB1,CC1,CA2,CB2,CC2,
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F5(P2X(IA),P2Y(IA),P2Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F5(P3X(IA),P3Y(IA),P3Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F5(P4X(IA),P4Y(IA),P4Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6))+
     * (PQ2(IA)/PQ3(IA))*(Q8*F5(P5X(IA),P5Y(IA),P5Z(IA),
     * CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F5(P6X(IA),P6Y(IA),P6Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F5(P7X(IA),P7Y(IA),P7Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F5(P8X(IA),P8Y(IA),P8Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6))
C
       A(IA,6)=(PQ1(IA)/PQ3(IA))*(Q8*F6(P1X(IA),P1Y(IA),P1Z(IA),
     * CA1,CB1,CC1,CA2,CB2,CC2,
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F6(P2X(IA),P2Y(IA),P2Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F6(P3X(IA),P3Y(IA),P3Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F6(P4X(IA),P4Y(IA),P4Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6))+
     * (PQ2(IA)/PQ3(IA))*(Q8*F6(P5X(IA),P5Y(IA),P5Z(IA),
     * CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F6(P6X(IA),P6Y(IA),P6Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F6(P7X(IA),P7Y(IA),P7Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)+
     * Q9*F6(P8X(IA),P8Y(IA),P8Z(IA),CA1,CB1,CC1,CA2,CB2,CC2,       
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6))
C
22     CONTINUE
C
       CALL INVERT(A,F,CKH,0)       
C
       DO 24 IK1=1,6
       DO 24 IK2=1,6
24     CK(IK1,IK2)=A(IK2,IK1)
C
       DO 26 IK=1,6
       COB(IK,1)=CK(IK,6)*CE1+CK(IK,5)*CD1
       COB(IK,2)=CK(IK,6)*CE2+CK(IK,5)*CD2
       COB(IK,3)=CK(IK,6)*CE3+CK(IK,5)*CD3
       COB(IK,4)=CK(IK,6)*CE4+CK(IK,5)*CD4
       COB(IK,5)=CK(IK,6)*CE5+CK(IK,5)*CD5
       COB(IK,6)=CK(IK,6)*CE6+CK(IK,5)*CD6
       COB(IK,7)=CK(IK,4)*CA3+CK(IK,3)*CA2+CK(IK,2)*CA1
       COB(IK,8)=CK(IK,4)*CB3+CK(IK,3)*CB2+CK(IK,2)*CB1
       COB(IK,9)=CK(IK,4)*CC3+CK(IK,3)*CC2+CK(IK,2)*CC1
       COB(IK,10)=CK(IK,1)
26     CONTINUE
      ENDIF
C
C

      IER=-1
C *** No second order derivatives available
C *** Used for second order problems only
      IF (BDER(5).OR.BDER(6).OR.BDER(7)) GOTO 99999
C
      IER=0
************************************************************************
C
      IF (.NOT.BDER(1)) GOTO 101
C
C *** Function values
      DBAS(1,1,1)= COB(1,1)*XI1**2+COB(1,2)*XI2**2+COB(1,3)*XI3**2
     *          +COB(1,4)*XI1*XI2+COB(1,5)*XI1*XI3+COB(1,6)*XI2*XI3
     *          +COB(1,7)*XI1+COB(1,8)*XI2+COB(1,9)*XI3+COB(1,10)
      DBAS(1,2,1)= COB(2,1)*XI1**2+COB(2,2)*XI2**2+COB(2,3)*XI3**2
     *          +COB(2,4)*XI1*XI2+COB(2,5)*XI1*XI3+COB(2,6)*XI2*XI3
     *          +COB(2,7)*XI1+COB(2,8)*XI2+COB(2,9)*XI3+COB(2,10)
      DBAS(1,3,1)= COB(3,1)*XI1**2+COB(3,2)*XI2**2+COB(3,3)*XI3**2
     *          +COB(3,4)*XI1*XI2+COB(3,5)*XI1*XI3+COB(3,6)*XI2*XI3
     *          +COB(3,7)*XI1+COB(3,8)*XI2+COB(3,9)*XI3+COB(3,10)
      DBAS(1,4,1)= COB(4,1)*XI1**2+COB(4,2)*XI2**2+COB(4,3)*XI3**2
     *          +COB(4,4)*XI1*XI2+COB(4,5)*XI1*XI3+COB(4,6)*XI2*XI3
     *          +COB(4,7)*XI1+COB(4,8)*XI2+COB(4,9)*XI3+COB(4,10)
      DBAS(1,5,1)= COB(5,1)*XI1**2+COB(5,2)*XI2**2+COB(5,3)*XI3**2
     *          +COB(5,4)*XI1*XI2+COB(5,5)*XI1*XI3+COB(5,6)*XI2*XI3
     *          +COB(5,7)*XI1+COB(5,8)*XI2+COB(5,9)*XI3+COB(5,10)
      DBAS(1,6,1)= COB(6,1)*XI1**2+COB(6,2)*XI2**2+COB(6,3)*XI3**2
     *          +COB(6,4)*XI1*XI2+COB(6,5)*XI1*XI3+COB(6,6)*XI2*XI3
     *          +COB(6,7)*XI1+COB(6,8)*XI2+COB(6,9)*XI3+COB(6,10)
C
101   IF (.NOT.(BDER(2).OR.BDER(3).OR.(BDER(4)))) GOTO 99999
C
C *** First order derivatives
      IF (.NOT.BDER(2)) GOTO 102
      DBAS(1,1,2)= 2D0*COB(1,1)*XI1+COB(1,4)*XI2+COB(1,5)*XI3+COB(1,7)
      DBAS(1,2,2)= 2D0*COB(2,1)*XI1+COB(2,4)*XI2+COB(2,5)*XI3+COB(2,7)
      DBAS(1,3,2)= 2D0*COB(3,1)*XI1+COB(3,4)*XI2+COB(3,5)*XI3+COB(3,7)
      DBAS(1,4,2)= 2D0*COB(4,1)*XI1+COB(4,4)*XI2+COB(4,5)*XI3+COB(4,7)
      DBAS(1,5,2)= 2D0*COB(5,1)*XI1+COB(5,4)*XI2+COB(5,5)*XI3+COB(5,7)
      DBAS(1,6,2)= 2D0*COB(6,1)*XI1+COB(6,4)*XI2+COB(6,5)*XI3+COB(6,7)
C
102   IF (.NOT.BDER(3)) GOTO 103
      DBAS(1,1,3)= 2D0*COB(1,2)*XI2+COB(1,4)*XI1+COB(1,6)*XI3+COB(1,8)
      DBAS(1,2,3)= 2D0*COB(2,2)*XI2+COB(2,4)*XI1+COB(2,6)*XI3+COB(2,8)
      DBAS(1,3,3)= 2D0*COB(3,2)*XI2+COB(3,4)*XI1+COB(3,6)*XI3+COB(3,8)
      DBAS(1,4,3)= 2D0*COB(4,2)*XI2+COB(4,4)*XI1+COB(4,6)*XI3+COB(4,8)
      DBAS(1,5,3)= 2D0*COB(5,2)*XI2+COB(5,4)*XI1+COB(5,6)*XI3+COB(5,8)
      DBAS(1,6,3)= 2D0*COB(6,2)*XI2+COB(6,4)*XI1+COB(6,6)*XI3+COB(6,8)
C
103   IF (.NOT.BDER(4)) GOTO 99999
      DBAS(1,1,4)= 2D0*COB(1,3)*XI3+COB(1,5)*XI1+COB(1,6)*XI2+COB(1,9)
      DBAS(1,2,4)= 2D0*COB(2,3)*XI3+COB(2,5)*XI1+COB(2,6)*XI2+COB(2,9)
      DBAS(1,3,4)= 2D0*COB(3,3)*XI3+COB(3,5)*XI1+COB(3,6)*XI2+COB(3,9)
      DBAS(1,4,4)= 2D0*COB(4,3)*XI3+COB(4,5)*XI1+COB(4,6)*XI2+COB(4,9)
      DBAS(1,5,4)= 2D0*COB(5,3)*XI3+COB(5,5)*XI1+COB(5,6)*XI2+COB(5,9)
      DBAS(1,6,4)= 2D0*COB(6,3)*XI3+COB(6,5)*XI1+COB(6,6)*XI2+COB(6,9)
c      read(5,*) nnn
c      if (nnn.eq.1) then
c      if ((xi1.eq.0.5).and.(xi2.eq.0.25).and.(xi3.eq.0.25)) then
c      write(*,*)'xi1,xi2,xi3',xi1,xi2,xi3,ipar
c      endif
c      write(*,*)'dbas(1,1,1)=',dbas(1,1,1)
c      write(*,*)'dbas(1,2,1)=',dbas(1,2,1)
c      write(*,*)'dbas(1,3,1)=',dbas(1,3,1)
c      write(*,*)'dbas(1,4,1)=',dbas(1,4,1)
c      write(*,*)'dbas(1,5,1)=',dbas(1,5,1)
c      write(*,*)'dbas(1,6,1)=',dbas(1,6,1)
c      endif
C
99999 END
C
C
C
      SUBROUTINE EM31(XI1,XI2,XI3,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNAE=6)
      PARAMETER (NNDIM=3,NNCOF=10)
      DIMENSION DXM(6),DYM(6),DZM(6),A(6,6),F(6),CKH(6),CK(6,6)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,
     *                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),
     *                IEL,NDIM
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COFBAS/ COB(NNBAS,NNCOF)
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/,/CHAR/,/CUB/
C
      F1(X,Y,Z,CA1,CB1,CC1,CA2,CB2,CC2,CA3,CB3,CC3,
     *   CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)=1D0
      F2(X,Y,Z,CA1,CB1,CC1,CA2,CB2,CC2,CA3,CB3,CC3,
     *   CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)=
     *                   CA1*X  +CB1*Y +CC1*Z
      F3(X,Y,Z,CA1,CB1,CC1,CA2,CB2,CC2,CA3,CB3,CC3,
     *   CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)=
     *                   CA2*X  +CB2*Y +CC2*Z
      F4(X,Y,Z,CA1,CB1,CC1,CA2,CB2,CC2,CA3,CB3,CC3,
     *   CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)=
     *                   CA3*X  +CB3*Y +CC3*Z
      F5(X,Y,Z,CA1,CB1,CC1,CA2,CB2,CC2,CA3,CB3,CC3,
     *   CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)=
     *  CD1*X*X+CD2*Y*Y+CD3*Z*Z+CD4*X*Y+CD5*X*Z+CD6*Y*Z
      F6(X,Y,Z,CA1,CB1,CC1,CA2,CB2,CC2,CA3,CB3,CC3,
     *   CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)=
     *  CE1*X*X+CE2*Y*Y+CE3*Z*Z+CE4*X*Y+CE5*X*Z+CE6*Y*Z
C
C
      SUB='EM31'
      IF (ICHECK.GE.998) CALL OTRC('EM31  ','01/08/94')
C
C
C *** Dummy call
      IF (IPAR.EQ.-1) THEN
       IER=0
       IPAR=31
       GOTO 99999
      ENDIF
C
C
      IF (IPAR.EQ.-2) THEN
       DXM(1)=0.25D0*(DX(1)+DX(2)+DX(3)+DX(4))
       DYM(1)=0.25D0*(DY(1)+DY(2)+DY(3)+DY(4))
       DZM(1)=0.25D0*(DZ(1)+DZ(2)+DZ(3)+DZ(4))
C
       DXM(2)=0.25D0*(DX(1)+DX(2)+DX(6)+DX(5))
       DYM(2)=0.25D0*(DY(1)+DY(2)+DY(6)+DY(5))
       DZM(2)=0.25D0*(DZ(1)+DZ(2)+DZ(6)+DZ(5))
C
       DXM(3)=0.25D0*(DX(2)+DX(3)+DX(7)+DX(6))
       DYM(3)=0.25D0*(DY(2)+DY(3)+DY(7)+DY(6))
       DZM(3)=0.25D0*(DZ(2)+DZ(3)+DZ(7)+DZ(6))
C
       DXM(4)=0.25D0*(DX(3)+DX(7)+DX(8)+DX(4))
       DYM(4)=0.25D0*(DY(3)+DY(7)+DY(8)+DY(4))
       DZM(4)=0.25D0*(DZ(3)+DZ(7)+DZ(8)+DZ(4))
C
       DXM(5)=0.25D0*(DX(1)+DX(4)+DX(8)+DX(5))
       DYM(5)=0.25D0*(DY(1)+DY(4)+DY(8)+DY(5))
       DZM(5)=0.25D0*(DZ(1)+DZ(4)+DZ(8)+DZ(5))
C
       DXM(6)=0.25D0*(DX(5)+DX(6)+DX(7)+DX(8))
       DYM(6)=0.25D0*(DY(5)+DY(6)+DY(7)+DY(8))
       DZM(6)=0.25D0*(DZ(5)+DZ(6)+DZ(7)+DZ(8))
C
       CA1=(DXM(1)-DXM(6))/SQRT((DXM(1)-DXM(6))**2+(DYM(1)-DYM(6))**2+
     *     (DZM(1)-DZM(6))**2)
       CB1=(DYM(1)-DYM(6))/SQRT((DXM(1)-DXM(6))**2+(DYM(1)-DYM(6))**2+
     *     (DZM(1)-DZM(6))**2)
       CC1=(DZM(1)-DZM(6))/SQRT((DXM(1)-DXM(6))**2+(DYM(1)-DYM(6))**2+
     *     (DZM(1)-DZM(6))**2)
       CA2=(DXM(2)-DXM(4))/SQRT((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2+
     *     (DZM(2)-DZM(4))**2)
       CB2=(DYM(2)-DYM(4))/SQRT((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2+
     *     (DZM(2)-DZM(4))**2)
       CC2=(DZM(2)-DZM(4))/SQRT((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2+
     *     (DZM(2)-DZM(4))**2)
       CA3=(DXM(3)-DXM(5))/SQRT((DXM(3)-DXM(5))**2+(DYM(3)-DYM(5))**2+
     *     (DZM(3)-DZM(5))**2)
       CB3=(DYM(3)-DYM(5))/SQRT((DXM(3)-DXM(5))**2+(DYM(3)-DYM(5))**2+
     *     (DZM(3)-DZM(5))**2)
       CC3=(DZM(3)-DZM(5))/SQRT((DXM(3)-DXM(5))**2+(DYM(3)-DYM(5))**2+
     *     (DZM(3)-DZM(5))**2)
       CD1=CA1**2-CA2**2
       CD2=CB1**2-CB2**2
       CD3=CC1**2-CC2**2
       CD4=2D0*CA1*CB1-2D0*CA2*CB2
       CD5=2D0*CA1*CC1-2D0*CA2*CC2
       CD6=2D0*CB1*CC1-2D0*CB2*CC2
       CE1=CA2**2-CA3**2
       CE2=CB2**2-CB3**2
       CE3=CC2**2-CC3**2
       CE4=2D0*CA2*CB2-2D0*CA3*CB3
       CE5=2D0*CA2*CC2-2D0*CA3*CC3
       CE6=2D0*CB2*CC2-2D0*CB3*CC3
************************************************************************
C
       DO 22 IA=1,6
       A(IA,1)=F1(DXM(IA),DYM(IA),DZM(IA),CA1,CB1,CC1,CA2,CB2,CC2,
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)       
       A(IA,2)=F2(DXM(IA),DYM(IA),DZM(IA),CA1,CB1,CC1,CA2,CB2,CC2,
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)       
       A(IA,3)=F3(DXM(IA),DYM(IA),DZM(IA),CA1,CB1,CC1,CA2,CB2,CC2,
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)       
       A(IA,4)=F4(DXM(IA),DYM(IA),DZM(IA),CA1,CB1,CC1,CA2,CB2,CC2,
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)       
       A(IA,5)=F5(DXM(IA),DYM(IA),DZM(IA),CA1,CB1,CC1,CA2,CB2,CC2,
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)       
       A(IA,6)=F6(DXM(IA),DYM(IA),DZM(IA),CA1,CB1,CC1,CA2,CB2,CC2,
     * CA3,CB3,CC3,CD1,CD2,CD3,CD4,CD5,CD6,CE1,CE2,CE3,CE4,CE5,CE6)       
22     CONTINUE
C
       CALL INVERT(A,F,CKH,0)       
C
       DO 24 IK1=1,6
       DO 24 IK2=1,6
24     CK(IK1,IK2)=A(IK2,IK1)
C
       DO 26 IK=1,6
       COB(IK,1)=CK(IK,6)*CE1+CK(IK,5)*CD1
       COB(IK,2)=CK(IK,6)*CE2+CK(IK,5)*CD2
       COB(IK,3)=CK(IK,6)*CE3+CK(IK,5)*CD3
       COB(IK,4)=CK(IK,6)*CE4+CK(IK,5)*CD4
       COB(IK,5)=CK(IK,6)*CE5+CK(IK,5)*CD5
       COB(IK,6)=CK(IK,6)*CE6+CK(IK,5)*CD6
       COB(IK,7)=CK(IK,4)*CA3+CK(IK,3)*CA2+CK(IK,2)*CA1
       COB(IK,8)=CK(IK,4)*CB3+CK(IK,3)*CB2+CK(IK,2)*CB1
       COB(IK,9)=CK(IK,4)*CC3+CK(IK,3)*CC2+CK(IK,2)*CC1
       COB(IK,10)=CK(IK,1)
26     CONTINUE
      ENDIF
C
C

      IER=-1
C *** No second order derivatives available
C *** Used for second order problems only
      IF (BDER(5).OR.BDER(6).OR.BDER(7)) GOTO 99999
C
      IER=0
************************************************************************
C
      IF (.NOT.BDER(1)) GOTO 101
C
C *** Function values
      DBAS(1,1,1)= COB(1,1)*XI1**2+COB(1,2)*XI2**2+COB(1,3)*XI3**2
     *          +COB(1,4)*XI1*XI2+COB(1,5)*XI1*XI3+COB(1,6)*XI2*XI3
     *          +COB(1,7)*XI1+COB(1,8)*XI2+COB(1,9)*XI3+COB(1,10)
      DBAS(1,2,1)= COB(2,1)*XI1**2+COB(2,2)*XI2**2+COB(2,3)*XI3**2
     *          +COB(2,4)*XI1*XI2+COB(2,5)*XI1*XI3+COB(2,6)*XI2*XI3
     *          +COB(2,7)*XI1+COB(2,8)*XI2+COB(2,9)*XI3+COB(2,10)
      DBAS(1,3,1)= COB(3,1)*XI1**2+COB(3,2)*XI2**2+COB(3,3)*XI3**2
     *          +COB(3,4)*XI1*XI2+COB(3,5)*XI1*XI3+COB(3,6)*XI2*XI3
     *          +COB(3,7)*XI1+COB(3,8)*XI2+COB(3,9)*XI3+COB(3,10)
      DBAS(1,4,1)= COB(4,1)*XI1**2+COB(4,2)*XI2**2+COB(4,3)*XI3**2
     *          +COB(4,4)*XI1*XI2+COB(4,5)*XI1*XI3+COB(4,6)*XI2*XI3
     *          +COB(4,7)*XI1+COB(4,8)*XI2+COB(4,9)*XI3+COB(4,10)
      DBAS(1,5,1)= COB(5,1)*XI1**2+COB(5,2)*XI2**2+COB(5,3)*XI3**2
     *          +COB(5,4)*XI1*XI2+COB(5,5)*XI1*XI3+COB(5,6)*XI2*XI3
     *          +COB(5,7)*XI1+COB(5,8)*XI2+COB(5,9)*XI3+COB(5,10)
      DBAS(1,6,1)= COB(6,1)*XI1**2+COB(6,2)*XI2**2+COB(6,3)*XI3**2
     *          +COB(6,4)*XI1*XI2+COB(6,5)*XI1*XI3+COB(6,6)*XI2*XI3
     *          +COB(6,7)*XI1+COB(6,8)*XI2+COB(6,9)*XI3+COB(6,10)
C
101   IF (.NOT.(BDER(2).OR.BDER(3).OR.(BDER(4)))) GOTO 99999
C
C *** First order derivatives
      IF (.NOT.BDER(2)) GOTO 102
      DBAS(1,1,2)= 2D0*COB(1,1)*XI1+COB(1,4)*XI2+COB(1,5)*XI3+COB(1,7)
      DBAS(1,2,2)= 2D0*COB(2,1)*XI1+COB(2,4)*XI2+COB(2,5)*XI3+COB(2,7)
      DBAS(1,3,2)= 2D0*COB(3,1)*XI1+COB(3,4)*XI2+COB(3,5)*XI3+COB(3,7)
      DBAS(1,4,2)= 2D0*COB(4,1)*XI1+COB(4,4)*XI2+COB(4,5)*XI3+COB(4,7)
      DBAS(1,5,2)= 2D0*COB(5,1)*XI1+COB(5,4)*XI2+COB(5,5)*XI3+COB(5,7)
      DBAS(1,6,2)= 2D0*COB(6,1)*XI1+COB(6,4)*XI2+COB(6,5)*XI3+COB(6,7)
C
102   IF (.NOT.BDER(3)) GOTO 103
      DBAS(1,1,3)= 2D0*COB(1,2)*XI2+COB(1,4)*XI1+COB(1,6)*XI3+COB(1,8)
      DBAS(1,2,3)= 2D0*COB(2,2)*XI2+COB(2,4)*XI1+COB(2,6)*XI3+COB(2,8)
      DBAS(1,3,3)= 2D0*COB(3,2)*XI2+COB(3,4)*XI1+COB(3,6)*XI3+COB(3,8)
      DBAS(1,4,3)= 2D0*COB(4,2)*XI2+COB(4,4)*XI1+COB(4,6)*XI3+COB(4,8)
      DBAS(1,5,3)= 2D0*COB(5,2)*XI2+COB(5,4)*XI1+COB(5,6)*XI3+COB(5,8)
      DBAS(1,6,3)= 2D0*COB(6,2)*XI2+COB(6,4)*XI1+COB(6,6)*XI3+COB(6,8)
C
103   IF (.NOT.BDER(4)) GOTO 99999
      DBAS(1,1,4)= 2D0*COB(1,3)*XI3+COB(1,5)*XI1+COB(1,6)*XI2+COB(1,9)
      DBAS(1,2,4)= 2D0*COB(2,3)*XI3+COB(2,5)*XI1+COB(2,6)*XI2+COB(2,9)
      DBAS(1,3,4)= 2D0*COB(3,3)*XI3+COB(3,5)*XI1+COB(3,6)*XI2+COB(3,9)
      DBAS(1,4,4)= 2D0*COB(4,3)*XI3+COB(4,5)*XI1+COB(4,6)*XI2+COB(4,9)
      DBAS(1,5,4)= 2D0*COB(5,3)*XI3+COB(5,5)*XI1+COB(5,6)*XI2+COB(5,9)
      DBAS(1,6,4)= 2D0*COB(6,3)*XI3+COB(6,5)*XI1+COB(6,6)*XI2+COB(6,9)
c      read(5,*) nnn
c      if (nnn.eq.1) then
c      if ((xi1.eq.0.5).and.(xi2.eq.0.25).and.(xi3.eq.0.25)) then
c      write(*,*)'xi1,xi2,xi3',xi1,xi2,xi3,ipar
c      endif
c      write(*,*)'dbas(1,1,1)=',dbas(1,1,1)
c      write(*,*)'dbas(1,2,1)=',dbas(1,2,1)
c      write(*,*)'dbas(1,3,1)=',dbas(1,3,1)
c      write(*,*)'dbas(1,4,1)=',dbas(1,4,1)
c      write(*,*)'dbas(1,5,1)=',dbas(1,5,1)
c      write(*,*)'dbas(1,6,1)=',dbas(1,6,1)
c      endif
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE INVERT(A,F,X,IPAR)
      DOUBLE PRECISION A,B,F,X
      PARAMETER (NDIM=6)
      DIMENSION A(NDIM,NDIM),B(NDIM,NDIM),F(NDIM),X(NDIM),
     *          MERKX(NDIM),MERKY(NDIM)
C
C
      IF (IPAR.EQ.0) THEN
       CALL AUSTAU(NDIM,NDIM,A,B,MERKX,MERKY,IFEHL)
       DO 10 IA=1,NDIM
       DO 10 IB=1,NDIM
10     A(IA,IB)=B(IA,IB)
      ENDIF
C
      IF (IPAR.EQ.1) THEN
       DO 20 IA=1,NDIM
       X(IA)=0D0
       DO 22 IB=1,NDIM
       X(IA)=X(IA)+A(IA,IB)*F(IB)
22     CONTINUE
20     CONTINUE
      ENDIF
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE AUSTAU(NDIM,N,A,B,MERKX,MERKY,IFEHL)
      DOUBLE PRECISION A,B,HILF,PIVOT
      DIMENSION A(NDIM,N),B(NDIM,N),MERKX(N),MERKY(N)
C
      IFEHL=1
      DO 100 I=1,N
      MERKX(I)=0
      MERKY(I)=0
      DO 100 L=1,N
100   B(I,L)=A(I,L)
C
      DO 400 I=1,N
      PIVOT=0D0
      DO 200 IX=1,N
      IF (MERKX(IX).NE.0) GOTO 200
      DO 180 IY=1,N
      IF (MERKY(IY).NE.0) GOTO 180
      IF (ABS(B(IX,IY)).LE.ABS(PIVOT)) GOTO 180
      PIVOT=B(IX,IY)
      INDX=IX
      INDY=IY
180   CONTINUE
200   CONTINUE
C
      IF (ABS(PIVOT).LE.0.0) GOTO 770
      MERKX(INDX)=INDY
      MERKY(INDY)=INDX
      B(INDX,INDY)=1D0/PIVOT
      DO 300 L=1,N
      IF (L.EQ.INDX) GOTO 300
      DO 280 M=1,N
      IF (M.EQ.INDY) GOTO 280
      B(L,M)=B(L,M)-B(L,INDY)*B(INDX,M)/PIVOT
280   CONTINUE
300   CONTINUE
C
      DO 390 IX=1,N
      IF (IX.NE.INDX) B(IX,INDY)=B(IX,INDY)/PIVOT
390   CONTINUE
      DO 400 IY=1,N
      IF (IY.NE.INDY) B(INDX,IY)=-B(INDX,IY)/PIVOT
400   CONTINUE
C
C
      DO 500 I=2,N
      IX=I-1
      IF (MERKX(IX).EQ.IX) GOTO 500
      DO 450 J=1,N
      IY=J
      IF (MERKX(IY).EQ.IX) GOTO 460
450   CONTINUE
460   DO 490 K=1,N
      HILF=B(IX,K)
      B(IX,K)=B(IY,K)
490   B(IY,K)=HILF
      MERKX(IY)=MERKX(IX)
500   MERKX(IX)=IX
      DO 600 I=2,N
      IX=I-1
      IF (MERKY(IX).EQ.IX) GOTO 600
      DO 550 J=1,N
      IY=J
      IF (MERKY(IY).EQ.IX) GOTO 560
550   CONTINUE
560   DO 590 K=1,N
      HILF=B(K,IX)
      B(K,IX)=B(K,IY)
590   B(K,IY)=HILF
      MERKY(IY)=MERKY(IX)
600   MERKY(IX)=IX
C
      IFEHL=0
770   RETURN
      END
C
C
************************************************************************
      SUBROUTINE SETIEL(DCORVG,KVERT,KAREA,KINT,NEL,NINT0,NINT1,NINT2)
************************************************************************
*    Purpose:  set element type for ref.element transformation
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION DCORVG(3,*),KVERT(NNVE,*),KAREA(NNAE,*),KINT(*)
      DIMENSION DNXA(6),DNYA(6),DNZA(6)
C
C
C
      CRIT1=1D-2
      CRIT2=1D-2
      CRIT3=1D-2
C
      NINT0=0
      NINT1=0
      NINT2=0
C
C
C
      DO 10 IEL=1,NEL
      KINT(IEL)=0
C
      IV1=KVERT(1,IEL)
      IV2=KVERT(2,IEL)
      IV3=KVERT(3,IEL)
      IV4=KVERT(4,IEL)
      IV5=KVERT(5,IEL)
      IV6=KVERT(6,IEL)
      IV7=KVERT(7,IEL)
      IV8=KVERT(8,IEL)
C
      XV1=DCORVG(1,IV1)
      YV1=DCORVG(2,IV1)
      ZV1=DCORVG(3,IV1)
      XV2=DCORVG(1,IV2)
      YV2=DCORVG(2,IV2)
      ZV2=DCORVG(3,IV2)
      XV3=DCORVG(1,IV3)
      YV3=DCORVG(2,IV3)
      ZV3=DCORVG(3,IV3)
      XV4=DCORVG(1,IV4)
      YV4=DCORVG(2,IV4)
      ZV4=DCORVG(3,IV4)
      XV5=DCORVG(1,IV5)
      YV5=DCORVG(2,IV5)
      ZV5=DCORVG(3,IV5)
      XV6=DCORVG(1,IV6)
      YV6=DCORVG(2,IV6)
      ZV6=DCORVG(3,IV6)
      XV7=DCORVG(1,IV7)
      YV7=DCORVG(2,IV7)
      ZV7=DCORVG(3,IV7)
      XV8=DCORVG(1,IV8)
      YV8=DCORVG(2,IV8)
      ZV8=DCORVG(3,IV8)
C
C *** PXC,PYC,PZC are coordinates of the center of the element
      PXC=0.125D0*(XV1+XV2+XV3+XV4+XV5+XV6+XV7+XV8)
      PYC=0.125D0*(YV1+YV2+YV3+YV4+YV5+YV6+YV7+YV8)
      PZC=0.125D0*(ZV1+ZV2+ZV3+ZV4+ZV5+ZV6+ZV7+ZV8)
C
C
C      WRITE(6,*) '&&&&',IEL
C      WRITE(6,*) '&&&&',XV1,YV1,ZV1
C      WRITE(6,*) '&&&&',XV2,YV2,ZV2
C      WRITE(6,*) '&&&&',XV3,YV3,ZV3
C      WRITE(6,*) '&&&&',XV4,YV4,ZV4
C      WRITE(6,*) '&&&&',XV5,YV5,ZV5
C      WRITE(6,*) '&&&&',XV6,YV6,ZV6
C      WRITE(6,*) '&&&&',XV7,YV7,ZV7
C      WRITE(6,*) '&&&&',XV8,YV8,ZV8
C
      DO 20 IA=1,6
C
      IF (IA.EQ.1) THEN
       IVT1=KVERT(1,IEL)
       IVT2=KVERT(2,IEL)
       IVT3=KVERT(3,IEL)
       IVT4=KVERT(4,IEL)
       GOTO 30
      ENDIF
C
      IF (IA.EQ.2) THEN
       IVT1=KVERT(1,IEL)
       IVT2=KVERT(2,IEL)
       IVT3=KVERT(6,IEL)
       IVT4=KVERT(5,IEL)
       GOTO 30
      ENDIF
C
      IF (IA.EQ.3) THEN
       IVT1=KVERT(2,IEL)
       IVT2=KVERT(3,IEL)
       IVT3=KVERT(6,IEL)
       IVT4=KVERT(7,IEL)
       GOTO 30
      ENDIF
C
      IF (IA.EQ.4) THEN
       IVT1=KVERT(3,IEL)
       IVT2=KVERT(4,IEL)
       IVT3=KVERT(8,IEL)
       IVT4=KVERT(7,IEL)
       GOTO 30
      ENDIF
C
      IF (IA.EQ.5) THEN
       IVT1=KVERT(4,IEL)
       IVT2=KVERT(1,IEL)
       IVT3=KVERT(5,IEL)
       IVT4=KVERT(8,IEL)
       GOTO 30
      ENDIF
C
      IF (IA.EQ.6) THEN
       IVT1=KVERT(5,IEL)
       IVT2=KVERT(6,IEL)
       IVT3=KVERT(7,IEL)
       IVT4=KVERT(8,IEL)
       GOTO 30
      ENDIF
C
C
30    P1X=DCORVG(1,IVT1)
      P1Y=DCORVG(2,IVT1)
      P1Z=DCORVG(3,IVT1)
      P2X=DCORVG(1,IVT2)
      P2Y=DCORVG(2,IVT2)
      P2Z=DCORVG(3,IVT2)
      P3X=DCORVG(1,IVT3)
      P3Y=DCORVG(2,IVT3)
      P3Z=DCORVG(3,IVT3)
      P4X=DCORVG(1,IVT4)
      P4Y=DCORVG(2,IVT4)
      P4Z=DCORVG(3,IVT4)
C
      AX2=P2X-P1X
      AY2=P2Y-P1Y
      AZ2=P2Z-P1Z
      AX3=P3X-P1X
      AY3=P3Y-P1Y
      AZ3=P3Z-P1Z
C
      PXA=(P1X+P2X+P3X+P4X)*0.25D0
      PYA=(P1Y+P2Y+P3Y+P4Y)*0.25D0
      PZA=(P1Z+P2Z+P3Z+P4Z)*0.25D0
C
      AX =PXC-PXA
      AY =PYC-PYA
      AZ =PZC-PZA
C
      DNX=(AY3*AZ2)-(AZ3*AY2)
      DNY=(AZ3*AX2)-(AX3*AZ2)
      DNZ=(AX3*AY2)-(AY3*AX2)
      DNLEN=SQRT(DNX*DNX+DNY*DNY+DNZ*DNZ)
C
      DNX =DNX/DNLEN
      DNY =DNY/DNLEN
      DNZ =DNZ/DNLEN
C
      DHN=DNX*AX+DNY*AY+DNZ*AZ
      IF (DHN.LT.0D0) THEN
       DNX=-DNX
       DNY=-DNY
       DNZ=-DNZ
      ENDIF
C
      DNXA(IA)=DNX     
      DNYA(IA)=DNY    
      DNZA(IA)=DNZ
C
      IF (IA.GT.3) THEN
       DNXA(IA)=-DNXA(IA)     
       DNYA(IA)=-DNYA(IA)    
       DNZA(IA)=-DNZA(IA)
      ENDIF     
C
C      WRITE(6,*) IA,DNXA(IA),DNYA(IA),DNZA(IA)
C
20    CONTINUE
C
C
      IF (ABS(DNXA(6)).GT.CRIT1) THEN
       HX1=DNXA(1)/DNXA(6)
      ELSE
       HX1=1D0
      ENDIF
C
      IF (ABS(DNYA(6)).GT.CRIT1) THEN
       HY1=DNYA(1)/DNYA(6)
      ELSE
       HY1=1D0
      ENDIF
C
      IF (ABS(DNZA(6)).GT.CRIT1) THEN
       HZ1=DNZA(1)/DNZA(6)
      ELSE
       HZ1=1D0
      ENDIF
C
C
      IF (ABS(DNXA(4)).GT.CRIT1) THEN
       HX2=DNXA(2)/DNXA(4)
      ELSE
       HX2=1D0
      ENDIF
C
      IF (ABS(DNYA(4)).GT.CRIT1) THEN
       HY2=DNYA(2)/DNYA(4)
      ELSE
       HY2=1D0
      ENDIF
C
      IF (ABS(DNZA(4)).GT.CRIT1) THEN
       HZ2=DNZA(2)/DNZA(4)
      ELSE
       HZ2=1D0
      ENDIF
C
C
      IF (ABS(DNXA(5)).GT.CRIT1) THEN
       HX3=DNXA(3)/DNXA(5)
      ELSE
       HX3=1D0
      ENDIF
C
      IF (ABS(DNYA(5)).GT.CRIT1) THEN
       HY3=DNYA(3)/DNYA(5)
      ELSE
       HY3=1D0
      ENDIF
C
      IF (ABS(DNZA(5)).GT.CRIT1) THEN
       HZ3=DNZA(3)/DNZA(5)
      ELSE
       HZ3=1D0
      ENDIF
C
C      WRITE(6,*) HX1,HY1,HZ1
C      WRITE(6,*) HX2,HY2,HZ2
C      WRITE(6,*) HX3,HY3,HZ3
C
      IH1=0
      IH2=0
      IH3=0
      IF ((ABS(HX1-1D0).LE.CRIT2).AND.(ABS(HY1-1D0).LE.CRIT2).AND.
     *    (ABS(HZ1-1D0).LE.CRIT2)) IH1=1
      IF ((ABS(HX2-1D0).LE.CRIT2).AND.(ABS(HY2-1D0).LE.CRIT2).AND.
     *    (ABS(HZ2-1D0).LE.CRIT2)) IH2=1
      IF ((ABS(HX3-1D0).LE.CRIT2).AND.(ABS(HY3-1D0).LE.CRIT2).AND.
     *    (ABS(HZ3-1D0).LE.CRIT2)) IH3=1
      IH=IH1*IH2*IH3
C
C      WRITE(6,*) '%%%',IEL,IH1,IH2,IH3,IH
C
      IF (IH.EQ.1) THEN
       H1=DNXA(1)*DNXA(2)+DNYA(1)*DNYA(2)+DNZA(1)*DNZA(2)
       H2=DNXA(1)*DNXA(3)+DNYA(1)*DNYA(3)+DNZA(1)*DNZA(3)
       H3=DNXA(2)*DNXA(3)+DNYA(2)*DNYA(3)+DNZA(2)*DNZA(3)
       IF ((ABS(H1).LE.CRIT3).AND.(ABS(H2).LE.CRIT3).AND.
     *     (ABS(H3).LE.CRIT3)) THEN
        ILINT=2
        NINT2=NINT2+1
       ELSE
        ILINT=1
        NINT1=NINT1+1
       ENDIF
      ELSE
       ILINT=0
       NINT0=NINT0+1
      ENDIF
C
C      WRITE(6,*) '%%%',IEL,ILINT,H1,H2,H3
C
      KINT(IEL)=ILINT
C
10    CONTINUE
C
      END
C
C
C
************************************************************************
      SUBROUTINE GRDIST(DCORVG,KNPR,NVT,DIEPS)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DCORVG(2,*),KNPR(*)
C
      H=1D0/(SQRT(DBLE(NVT))-1)
      HDIST=DIEPS*H
C
      DO 1 IVT=1,NVT
      IF (KNPR(IVT).EQ.0) THEN
       DCORVG(1,IVT)=DCORVG(1,IVT)+DBLE((-1)**MOD(IVT,17))*HDIST
       DCORVG(2,IVT)=DCORVG(2,IVT)+DBLE((-1)**IVT)*HDIST
      ENDIF
1     CONTINUE
C
      END
C
C
C
************************************************************************
      SUBROUTINE CHCOOR(DCORVG,KVERT,KAREA,KADJ,KNPR,NEL,NVT)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DCORVG(3,*),KVERT(8,*),KAREA(6,*),KADJ(6,*),KNPR(*)
C
C
C
      DO 100 IEL=1,NEL/8
C
      IEL1=IEL
      IEL2=KADJ(3,IEL1)
      IEL3=KADJ(3,IEL2)
      IEL4=KADJ(3,IEL3)
      IEL5=KADJ(6,IEL1)
      IEL6=KADJ(3,IEL5)
      IEL7=KADJ(3,IEL6)
      IEL8=KADJ(3,IEL7)
C
      IADJ1=KADJ(1,IEL1)
      IADJ2=KADJ(2,IEL1)
      IADJ3=KADJ(2,IEL2)
      IADJ4=KADJ(2,IEL3)
      IADJ5=KADJ(2,IEL4)
      IADJ6=KADJ(1,IEL5)
C
      IVT1=KVERT(3,IEL1)
      IVT2=KVERT(6,IEL1)
      IVT3=KVERT(6,IEL2)
      IVT4=KVERT(6,IEL3)
      IVT5=KVERT(6,IEL4)
      IVT6=KVERT(3,IEL5)
C
      JVT1 =KVERT(2,IEL1)
      JVT2 =KVERT(2,IEL2)
      JVT3 =KVERT(2,IEL3)
      JVT4 =KVERT(2,IEL4)
      JVT5 =KVERT(5,IEL1)
      JVT6 =KVERT(5,IEL2)
      JVT7 =KVERT(5,IEL3)
      JVT8 =KVERT(5,IEL4)
      JVT9 =KVERT(2,IEL5)
      JVT10=KVERT(2,IEL6)
      JVT11=KVERT(2,IEL7)
      JVT12=KVERT(2,IEL8)
C
      NADJ0=0
      IF (IADJ1.EQ.0) NADJ0=NADJ0+1
      IF (IADJ2.EQ.0) NADJ0=NADJ0+1
      IF (IADJ3.EQ.0) NADJ0=NADJ0+1
      IF (IADJ4.EQ.0) NADJ0=NADJ0+1
      IF (IADJ5.EQ.0) NADJ0=NADJ0+1
      IF (IADJ6.EQ.0) NADJ0=NADJ0+1
C
C
      IF (NADJ0.EQ.0) GOTO 100
C
C
      QX1 =DCORVG(1,JVT1)
      QX2 =DCORVG(1,JVT2)
      QX3 =DCORVG(1,JVT3)
      QX4 =DCORVG(1,JVT4)
      QX5 =DCORVG(1,JVT5)
      QX6 =DCORVG(1,JVT6)
      QX7 =DCORVG(1,JVT7)
      QX8 =DCORVG(1,JVT8)
      QX9 =DCORVG(1,JVT9)
      QX10=DCORVG(1,JVT10)
      QX11=DCORVG(1,JVT11)
      QX12=DCORVG(1,JVT12)
C
      QY1 =DCORVG(2,JVT1)
      QY2 =DCORVG(2,JVT2)
      QY3 =DCORVG(2,JVT3)
      QY4 =DCORVG(2,JVT4)
      QY5 =DCORVG(2,JVT5)
      QY6 =DCORVG(2,JVT6)
      QY7 =DCORVG(2,JVT7)
      QY8 =DCORVG(2,JVT8)
      QY9 =DCORVG(2,JVT9)
      QY10=DCORVG(2,JVT10)
      QY11=DCORVG(2,JVT11)
      QY12=DCORVG(2,JVT12)
C
      QZ1 =DCORVG(3,JVT1)
      QZ2 =DCORVG(3,JVT2)
      QZ3 =DCORVG(3,JVT3)
      QZ4 =DCORVG(3,JVT4)
      QZ5 =DCORVG(3,JVT5)
      QZ6 =DCORVG(3,JVT6)
      QZ7 =DCORVG(3,JVT7)
      QZ8 =DCORVG(3,JVT8)
      QZ9 =DCORVG(3,JVT9)
      QZ10=DCORVG(3,JVT10)
      QZ11=DCORVG(3,JVT11)
      QZ12=DCORVG(3,JVT12)
C
C
      IF (NADJ0.EQ.1) THEN
       IF ((IADJ1.EQ.0).OR.(IADJ6.EQ.0)) THEN
        DCORVG(1,IVT2)=0.5D0*(QX1 +QX9)
        DCORVG(2,IVT2)=0.5D0*(QY1 +QY9)
        DCORVG(3,IVT2)=0.5D0*(QZ1 +QZ9)
        DCORVG(1,IVT3)=0.5D0*(QX2 +QX10)
        DCORVG(2,IVT3)=0.5D0*(QY2 +QY10)
        DCORVG(3,IVT3)=0.5D0*(QZ2 +QZ10)
        DCORVG(1,IVT4)=0.5D0*(QX3 +QX11)
        DCORVG(2,IVT4)=0.5D0*(QY3 +QY11)
        DCORVG(3,IVT4)=0.5D0*(QZ3 +QZ11)
        DCORVG(1,IVT5)=0.5D0*(QX4 +QX12)
        DCORVG(2,IVT5)=0.5D0*(QY4 +QY12)
        DCORVG(3,IVT5)=0.5D0*(QZ4 +QZ12)
       ENDIF
       IF ((IADJ2.EQ.0).OR.(IADJ4.EQ.0)) THEN
        DCORVG(1,IVT1)=0.5D0*(QX1 +QX3)
        DCORVG(2,IVT1)=0.5D0*(QY1 +QY3)
        DCORVG(3,IVT1)=0.5D0*(QZ1 +QZ3)
        DCORVG(1,IVT3)=0.5D0*(QX6 +QX7)
        DCORVG(2,IVT3)=0.5D0*(QY6 +QY7)
        DCORVG(3,IVT3)=0.5D0*(QZ6 +QZ7)
        DCORVG(1,IVT5)=0.5D0*(QX5 +QX8)
        DCORVG(2,IVT5)=0.5D0*(QY5 +QY8)
        DCORVG(3,IVT5)=0.5D0*(QZ5 +QZ8)
        DCORVG(1,IVT6)=0.5D0*(QX9 +QX11)
        DCORVG(2,IVT6)=0.5D0*(QY9 +QY11)
        DCORVG(3,IVT6)=0.5D0*(QZ9 +QZ11)
       ENDIF
       IF ((IADJ3.EQ.0).OR.(IADJ5.EQ.0)) THEN
        DCORVG(1,IVT1)=0.5D0*(QX2 +QX4)
        DCORVG(2,IVT1)=0.5D0*(QY2 +QY4)
        DCORVG(3,IVT1)=0.5D0*(QZ2 +QZ4)
        DCORVG(1,IVT2)=0.5D0*(QX5 +QX6)
        DCORVG(2,IVT2)=0.5D0*(QY5 +QY6)
        DCORVG(3,IVT2)=0.5D0*(QZ5 +QZ6)
        DCORVG(1,IVT4)=0.5D0*(QX7 +QX8)
        DCORVG(2,IVT4)=0.5D0*(QY7 +QY8)
        DCORVG(3,IVT4)=0.5D0*(QZ7 +QZ8)
        DCORVG(1,IVT6)=0.5D0*(QX10+QX12)
        DCORVG(2,IVT6)=0.5D0*(QY10+QY12)
        DCORVG(3,IVT6)=0.5D0*(QZ10+QZ12)
       ENDIF
      ENDIF
C
      IF (NADJ0.GT.1) THEN
       DCORVG(1,IVT1)=0.25D0*(QX1 +QX2 +QX3 +QX4)
       DCORVG(1,IVT2)=0.25D0*(QX1 +QX5 +QX6 +QX9)
       DCORVG(1,IVT3)=0.25D0*(QX2 +QX6 +QX7 +QX10)
       DCORVG(1,IVT4)=0.25D0*(QX3 +QX7 +QX8 +QX11)
       DCORVG(1,IVT5)=0.25D0*(QX4 +QX5 +QX8 +QX12)
       DCORVG(1,IVT6)=0.25D0*(QX9 +QX10+QX11+QX12)
       DCORVG(2,IVT1)=0.25D0*(QY1 +QY2 +QY3 +QY4)
       DCORVG(2,IVT2)=0.25D0*(QY1 +QY5 +QY6 +QY9)
       DCORVG(2,IVT3)=0.25D0*(QY2 +QY6 +QY7 +QY10)
       DCORVG(2,IVT4)=0.25D0*(QY3 +QY7 +QY8 +QY11)
       DCORVG(2,IVT5)=0.25D0*(QY4 +QY5 +QY8 +QY12)
       DCORVG(2,IVT6)=0.25D0*(QY9 +QY10+QY11+QY12)
       DCORVG(3,IVT1)=0.25D0*(QZ1 +QZ2 +QZ3 +QZ4)
       DCORVG(3,IVT2)=0.25D0*(QZ1 +QZ5 +QZ6 +QZ9)
       DCORVG(3,IVT3)=0.25D0*(QZ2 +QZ6 +QZ7 +QZ10)
       DCORVG(3,IVT4)=0.25D0*(QZ3 +QZ7 +QZ8 +QZ11)
       DCORVG(3,IVT5)=0.25D0*(QZ4 +QZ5 +QZ8 +QZ12)
       DCORVG(3,IVT6)=0.25D0*(QZ9 +QZ10+QZ11+QZ12)
      ENDIF
C
C
      IVTM=KVERT(7,IEL)
      PXM=DCORVG(1,IVTM)
      PYM=DCORVG(2,IVTM)
      PZM=DCORVG(3,IVTM)
C
      PX1=DCORVG(1,IVT1)
      PX2=DCORVG(1,IVT2)
      PX3=DCORVG(1,IVT3)
      PX4=DCORVG(1,IVT4)
      PX5=DCORVG(1,IVT5)
      PX6=DCORVG(1,IVT6)
C
      PY1=DCORVG(2,IVT1)
      PY2=DCORVG(2,IVT2)
      PY3=DCORVG(2,IVT3)
      PY4=DCORVG(2,IVT4)
      PY5=DCORVG(2,IVT5)
      PY6=DCORVG(2,IVT6)
C
      PZ1=DCORVG(3,IVT1)
      PZ2=DCORVG(3,IVT2)
      PZ3=DCORVG(3,IVT3)
      PZ4=DCORVG(3,IVT4)
      PZ5=DCORVG(3,IVT5)
      PZ6=DCORVG(3,IVT6)
C
      IF (NADJ0.EQ.1) THEN
       IF ((IADJ1.EQ.0).OR.(IADJ6.EQ.0)) THEN
        PX=0.5D0*(PX1+PX6)
        PY=0.5D0*(PY1+PY6)
        PZ=0.5D0*(PZ1+PZ6)
       ENDIF
       IF ((IADJ2.EQ.0).OR.(IADJ4.EQ.0)) THEN
        PX=0.5D0*(PX2+PX4)
        PY=0.5D0*(PY2+PY4)
        PZ=0.5D0*(PZ2+PZ4)
       ENDIF
       IF ((IADJ3.EQ.0).OR.(IADJ5.EQ.0)) THEN
        PX=0.5D0*(PX3+PX5)
        PY=0.5D0*(PY3+PY5)
        PZ=0.5D0*(PZ3+PZ5)
       ENDIF
       DCORVG(1,IVTM)=PX
       DCORVG(2,IVTM)=PY
       DCORVG(3,IVTM)=PZ
      ENDIF
C
C      
      IF (NADJ0.GT.1) THEN
       PX=1D0/6D0*(PX1+PX2+PX3+PX4+PX5+PX6)
       PY=1D0/6D0*(PY1+PY2+PY3+PY4+PY5+PY6)
       PZ=1D0/6D0*(PZ1+PZ2+PZ3+PZ4+PZ5+PZ6)
       DCORVG(1,IVTM)=PX
       DCORVG(2,IVTM)=PY
       DCORVG(3,IVTM)=PZ
      ENDIF
C
C
100   CONTINUE
C
C
C
      END
