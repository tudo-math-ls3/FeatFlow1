************************************************************************
      SUBROUTINE BDRNEU(KABD,KNPR,DCORVG,INEUM,KAREA,KADJ,KVERT,KELBD,
     *                  II,NLMAX)
************************************************************************
*    Purpose:  sets the DIRICHLET- and NEUMANN components
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION KABD(*),KELBD(*),KAREA(6,*),KADJ(6,*)
      DIMENSION KNPR(*),DCORVG(3,*),KVERT(8,*)
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      SAVE
C
C
      NABD =0
      INEUM=0
C
C
      DO 1 IEL=1,NEL
C
      DO 2 IAT=1,6
      IAR =KAREA(IAT,IEL)
      INPR=KNPR(NVT+IAR)
      IF (INPR.EQ.0) GOTO 2
C
      NABD=NABD+1
      KABD (NABD)=IAR
      KELBD(NABD)=IEL
C
C
      IVT1=KVERT(1,IEL)
      IVT2=KVERT(2,IEL)
      IVT3=KVERT(3,IEL)
      IVT4=KVERT(4,IEL)
      IVT5=KVERT(5,IEL)
      IVT6=KVERT(6,IEL)
      IVT7=KVERT(7,IEL)
      IVT8=KVERT(8,IEL)
C
      XV1=DCORVG(1,IVT1)
      YV1=DCORVG(2,IVT1)
      ZV1=DCORVG(3,IVT1)
      XV2=DCORVG(1,IVT2)
      YV2=DCORVG(2,IVT2)
      ZV2=DCORVG(3,IVT2)
      XV3=DCORVG(1,IVT3)
      YV3=DCORVG(2,IVT3)
      ZV3=DCORVG(3,IVT3)
      XV4=DCORVG(1,IVT4)
      YV4=DCORVG(2,IVT4)
      ZV4=DCORVG(3,IVT4)
      XV5=DCORVG(1,IVT5)
      YV5=DCORVG(2,IVT5)
      ZV5=DCORVG(3,IVT5)
      XV6=DCORVG(1,IVT6)
      YV6=DCORVG(2,IVT6)
      ZV6=DCORVG(3,IVT6)
      XV7=DCORVG(1,IVT7)
      YV7=DCORVG(2,IVT7)
      ZV7=DCORVG(3,IVT7)
      XV8=DCORVG(1,IVT8)
      YV8=DCORVG(2,IVT8)
      ZV8=DCORVG(3,IVT8)
C
      IF (IAT.EQ.1) THEN
       PX=(XV1+XV2+XV3+XV4)*0.25D0
       PY=(YV1+YV2+YV3+YV4)*0.25D0
       PZ=(ZV1+ZV2+ZV3+ZV4)*0.25D0
      ENDIF
C
      IF (IAT.EQ.2) THEN
       PX=(XV1+XV2+XV6+XV5)*0.25D0
       PY=(YV1+YV2+YV6+YV5)*0.25D0
       PZ=(ZV1+ZV2+ZV6+ZV5)*0.25D0
      ENDIF
C
      IF (IAT.EQ.3) THEN
       PX=(XV2+XV3+XV7+XV6)*0.25D0
       PY=(YV2+YV3+YV7+YV6)*0.25D0
       PZ=(ZV2+ZV3+ZV7+ZV6)*0.25D0
      ENDIF
C
      IF (IAT.EQ.4) THEN
       PX=(XV3+XV4+XV8+XV7)*0.25D0
       PY=(YV3+YV4+YV8+YV7)*0.25D0
       PZ=(ZV3+ZV4+ZV8+ZV7)*0.25D0
      ENDIF
C
      IF (IAT.EQ.5) THEN
       PX=(XV4+XV1+XV5+XV8)*0.25D0
       PY=(YV4+YV1+YV5+YV8)*0.25D0
       PZ=(ZV4+ZV1+ZV5+ZV8)*0.25D0
      ENDIF
C
      IF (IAT.EQ.6) THEN
       PX=(XV5+XV6+XV7+XV8)*0.25D0
       PY=(YV5+YV6+YV7+YV8)*0.25D0
       PZ=(ZV5+ZV6+ZV7+ZV8)*0.25D0
      ENDIF
C
      IFLAG=0
      CALL NEUDAT(IEL,INPR,PX,PY,PZ,TIMENS,IFLAG)
C
      IF (IFLAG.EQ.1) THEN
       INEUM=1
       KNPR(NVT+IAR)=0
       KABD(NABD) =-KABD(NABD)
       KELBD(NABD)=-KELBD(NABD)
      ENDIF
C
2     CONTINUE
C
1     CONTINUE
C
C
C
      DO 3 IEL=1,NEL
C
      DO 33 IVE=1,8
      IVT=KVERT(IVE,IEL)
      PX=DCORVG(1,IVT)
      PY=DCORVG(2,IVT)
      PZ=DCORVG(3,IVT)
      PXM=0D0
      PYM=0D0
      PZM=0D0
      DIST=SQRT((PX-PXM)**2+(PY-PYM)**2+(PZ-PZM)**2)
      IF (DIST.LE.0.5D0) THEN
       KNPR(IVT)=11
      ENDIF
33    CONTINUE
C
      DO 4 IAT=1,6
      IAR =KAREA(IAT,IEL)
      IF (KADJ(IAT,IEL).LT.IEL) THEN
C
      IVT1=KVERT(1,IEL)
      IVT2=KVERT(2,IEL)
      IVT3=KVERT(3,IEL)
      IVT4=KVERT(4,IEL)
      IVT5=KVERT(5,IEL)
      IVT6=KVERT(6,IEL)
      IVT7=KVERT(7,IEL)
      IVT8=KVERT(8,IEL)
C
      XV1=DCORVG(1,IVT1)
      YV1=DCORVG(2,IVT1)
      ZV1=DCORVG(3,IVT1)
      XV2=DCORVG(1,IVT2)
      YV2=DCORVG(2,IVT2)
      ZV2=DCORVG(3,IVT2)
      XV3=DCORVG(1,IVT3)
      YV3=DCORVG(2,IVT3)
      ZV3=DCORVG(3,IVT3)
      XV4=DCORVG(1,IVT4)
      YV4=DCORVG(2,IVT4)
      ZV4=DCORVG(3,IVT4)
      XV5=DCORVG(1,IVT5)
      YV5=DCORVG(2,IVT5)
      ZV5=DCORVG(3,IVT5)
      XV6=DCORVG(1,IVT6)
      YV6=DCORVG(2,IVT6)
      ZV6=DCORVG(3,IVT6)
      XV7=DCORVG(1,IVT7)
      YV7=DCORVG(2,IVT7)
      ZV7=DCORVG(3,IVT7)
      XV8=DCORVG(1,IVT8)
      YV8=DCORVG(2,IVT8)
      ZV8=DCORVG(3,IVT8)
C
      IF (IAT.EQ.1) THEN
       PX=(XV1+XV2+XV3+XV4)*0.25D0
       PY=(YV1+YV2+YV3+YV4)*0.25D0
       PZ=(ZV1+ZV2+ZV3+ZV4)*0.25D0
      ENDIF
C
      IF (IAT.EQ.2) THEN
       PX=(XV1+XV2+XV6+XV5)*0.25D0
       PY=(YV1+YV2+YV6+YV5)*0.25D0
       PZ=(ZV1+ZV2+ZV6+ZV5)*0.25D0
      ENDIF
C
      IF (IAT.EQ.3) THEN
       PX=(XV2+XV3+XV6+XV7)*0.25D0
       PY=(YV2+YV3+YV6+YV7)*0.25D0
       PZ=(ZV2+ZV3+ZV6+ZV7)*0.25D0
      ENDIF
C
      IF (IAT.EQ.4) THEN
       PX=(XV3+XV4+XV8+XV7)*0.25D0
       PY=(YV3+YV4+YV8+YV7)*0.25D0
       PZ=(ZV3+ZV4+ZV8+ZV7)*0.25D0
      ENDIF
C
      IF (IAT.EQ.5) THEN
       PX=(XV4+XV1+XV5+XV8)*0.25D0
       PY=(YV4+YV1+YV5+YV8)*0.25D0
       PZ=(ZV4+ZV1+ZV5+ZV8)*0.25D0
      ENDIF
C
      IF (IAT.EQ.6) THEN
       PX=(XV5+XV6+XV7+XV8)*0.25D0
       PY=(YV5+YV6+YV7+YV8)*0.25D0
       PZ=(ZV5+ZV6+ZV7+ZV8)*0.25D0
      ENDIF
C
C
      PXM=0D0
      PYM=0D0
      PZM=0D0
      DIST=SQRT((PX-PXM)**2+(PY-PYM)**2+(PZ-PZM)**2)
      IF (DIST.LE.0.5D0) THEN
       KNPR(IAR+NVT)=-1
       NABD=NABD+1
       KABD(NABD)=IAR
      ENDIF
C
      ENDIF
C
4     CONTINUE
C
3     CONTINUE
      END
C
************************************************************************
      SUBROUTINE BDRSET(DU1,DU2,DU3,DF1,DF2,DF3,KABD,NABD,DCORVG,
     *                  KVERT,KAREA,KELBD,UE)
************************************************************************
*    Purpose:  updates the solution vector (DU1,DU2,DU3) and the right 
*              hand side (DF1,DF2,DF3) for all DIRICHLET boundary nodes
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION DU1(*),DU2(*),DU3(*),DF1(*),DF2(*),DF3(*)
      DIMENSION KABD(*),KELBD(*),DCORVG(3,*)
      DIMENSION KVERT(8,*),KAREA(6,*)
      EXTERNAL UE
C
C
      DO 1 IAT=1,NABD
      IABD =KABD (IAT)
      IELBD=KELBD(IAT)
      IF (IABD.LT.0) GOTO 1
C
      IVT1=KVERT(1,IELBD)
      IVT2=KVERT(2,IELBD)
      IVT3=KVERT(3,IELBD)
      IVT4=KVERT(4,IELBD)
      IVT5=KVERT(5,IELBD)
      IVT6=KVERT(6,IELBD)
      IVT7=KVERT(7,IELBD)
      IVT8=KVERT(8,IELBD)
C
      XV1=DCORVG(1,IVT1)
      YV1=DCORVG(2,IVT1)
      ZV1=DCORVG(3,IVT1)
      XV2=DCORVG(1,IVT2)
      YV2=DCORVG(2,IVT2)
      ZV2=DCORVG(3,IVT2)
      XV3=DCORVG(1,IVT3)
      YV3=DCORVG(2,IVT3)
      ZV3=DCORVG(3,IVT3)
      XV4=DCORVG(1,IVT4)
      YV4=DCORVG(2,IVT4)
      ZV4=DCORVG(3,IVT4)
      XV5=DCORVG(1,IVT5)
      YV5=DCORVG(2,IVT5)
      ZV5=DCORVG(3,IVT5)
      XV6=DCORVG(1,IVT6)
      YV6=DCORVG(2,IVT6)
      ZV6=DCORVG(3,IVT6)
      XV7=DCORVG(1,IVT7)
      YV7=DCORVG(2,IVT7)
      ZV7=DCORVG(3,IVT7)
      XV8=DCORVG(1,IVT8)
      YV8=DCORVG(2,IVT8)
      ZV8=DCORVG(3,IVT8)
C
      IF (IABD.EQ.KAREA(1,IELBD)) THEN
       PX=(XV1+XV2+XV3+XV4)*0.25D0
       PY=(YV1+YV2+YV3+YV4)*0.25D0
       PZ=(ZV1+ZV2+ZV3+ZV4)*0.25D0
      ENDIF
C
      IF (IABD.EQ.KAREA(2,IELBD)) THEN
       PX=(XV1+XV2+XV6+XV5)*0.25D0
       PY=(YV1+YV2+YV6+YV5)*0.25D0
       PZ=(ZV1+ZV2+ZV6+ZV5)*0.25D0
      ENDIF
C
      IF (IABD.EQ.KAREA(3,IELBD)) THEN
       PX=(XV2+XV3+XV7+XV6)*0.25D0
       PY=(YV2+YV3+YV7+YV6)*0.25D0
       PZ=(ZV2+ZV3+ZV7+ZV6)*0.25D0
      ENDIF
C
      IF (IABD.EQ.KAREA(4,IELBD)) THEN
       PX=(XV3+XV4+XV8+XV7)*0.25D0
       PY=(YV3+YV4+YV8+YV7)*0.25D0
       PZ=(ZV3+ZV4+ZV8+ZV7)*0.25D0
      ENDIF
C
      IF (IABD.EQ.KAREA(5,IELBD)) THEN
       PX=(XV4+XV1+XV5+XV8)*0.25D0
       PY=(YV4+YV1+YV5+YV8)*0.25D0
       PZ=(ZV4+ZV1+ZV5+ZV8)*0.25D0
      ENDIF
C
      IF (IABD.EQ.KAREA(6,IELBD)) THEN
       PX=(XV5+XV6+XV7+XV8)*0.25D0
       PY=(YV5+YV6+YV7+YV8)*0.25D0
       PZ=(ZV5+ZV6+ZV7+ZV8)*0.25D0
      ENDIF
C
      U1=UE(PX,PY,PZ,1)
      U2=UE(PX,PY,PZ,2)
      U3=UE(PX,PY,PZ,3)
C
      DU1(IABD)= U1
      DU2(IABD)= U2
      DU3(IABD)= U3
      DF1(IABD)= U1 
      DF2(IABD)= U2 
      DF3(IABD)= U3
C
1     CONTINUE
C
      END
C
************************************************************************
      SUBROUTINE PDSET(DF1,DF2,DF3,KABD,NABD,DCORVG,KELBD,
     *                 KVERT,KAREA,TSTEPB)
************************************************************************
*    Purpose:  updates the right hand 
*              side (DF1,DF2,DF3) for all NEUMANN boundary nodes
*              only necEssary, if NEUMANN boundary implemented
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION DF1(*),DF2(*),DF3(*),KELBD(*),KVERT(8,*),KAREA(6,*)
      DIMENSION KABD(*),DCORVG(3,*)
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
      SAVE
C
C
      DO 10 IAT=1,NABD
      IABD =KABD(IAT)
      IELBD=KELBD(IAT)
C
      IF ((IABD.LT.0).AND.(IELBD.LT.0)) THEN
C
       IV1=KVERT(1,-IELBD)
       IV2=KVERT(2,-IELBD)
       IV3=KVERT(3,-IELBD)
       IV4=KVERT(4,-IELBD)
       IV5=KVERT(5,-IELBD)
       IV6=KVERT(6,-IELBD)
       IV7=KVERT(7,-IELBD)
       IV8=KVERT(8,-IELBD)
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
C ***  PXC,PYC,PZC are coordinates of the center of the element
       PXC=0.125D0*(XV1+XV2+XV3+XV4+XV5+XV6+XV7+XV8)
       PYC=0.125D0*(YV1+YV2+YV3+YV4+YV5+YV6+YV7+YV8)
       PZC=0.125D0*(ZV1+ZV2+ZV3+ZV4+ZV5+ZV6+ZV7+ZV8)
C
C
       IF (-IABD.EQ.KAREA(1,-IELBD)) THEN
        PX=(XV1+XV2+XV3+XV4)*0.25D0
        PY=(YV1+YV2+YV3+YV4)*0.25D0
        PZ=(ZV1+ZV2+ZV3+ZV4)*0.25D0
       ENDIF
C
       IF (-IABD.EQ.KAREA(2,-IELBD)) THEN
        PX=(XV1+XV2+XV6+XV5)*0.25D0
        PY=(YV1+YV2+YV6+YV5)*0.25D0
        PZ=(ZV1+ZV2+ZV6+ZV5)*0.25D0
       ENDIF
C
       IF (-IABD.EQ.KAREA(3,-IELBD)) THEN
        PX=(XV2+XV3+XV7+XV6)*0.25D0
        PY=(YV2+YV3+YV7+YV6)*0.25D0
        PZ=(ZV2+ZV3+ZV7+ZV6)*0.25D0
       ENDIF
C
       IF (-IABD.EQ.KAREA(4,-IELBD)) THEN
        PX=(XV3+XV4+XV8+XV7)*0.25D0
        PY=(YV3+YV4+YV8+YV7)*0.25D0
        PZ=(ZV3+ZV4+ZV8+ZV7)*0.25D0
       ENDIF
C
       IF (-IABD.EQ.KAREA(5,-IELBD)) THEN
        PX=(XV4+XV1+XV5+XV8)*0.25D0
        PY=(YV4+YV1+YV5+YV8)*0.25D0
        PZ=(ZV4+ZV1+ZV5+ZV8)*0.25D0
       ENDIF
C
       IF (-IABD.EQ.KAREA(6,-IELBD)) THEN
        PX=(XV5+XV6+XV7+XV8)*0.25D0
        PY=(YV5+YV6+YV7+YV8)*0.25D0
        PZ=(ZV5+ZV6+ZV7+ZV8)*0.25D0
       ENDIF
C
C
       IFLAG=0
       CALL NEUDAT(-IELBD,INPR,PX,PY,PZ,TIMENS,IFLAG)
C
       IF (IFLAG.EQ.1) THEN
C
        DO 20 IA=1,6
        IAH=KAREA(IA,-IELBD)
C
        IF (IAH.EQ.-IABD) THEN
C
         IF (IA.EQ.1) THEN
          IVT1=KVERT(1,-IELBD)
          IVT2=KVERT(2,-IELBD)
          IVT3=KVERT(3,-IELBD)
          IVT4=KVERT(4,-IELBD)
          GOTO 30
         ENDIF
C
         IF (IA.EQ.2) THEN
          IVT1=KVERT(1,-IELBD)
          IVT2=KVERT(2,-IELBD)
          IVT3=KVERT(6,-IELBD)
          IVT4=KVERT(5,-IELBD)
          GOTO 30
         ENDIF
C
         IF (IA.EQ.3) THEN
          IVT1=KVERT(2,-IELBD)
          IVT2=KVERT(3,-IELBD)
          IVT3=KVERT(7,-IELBD)
          IVT4=KVERT(6,-IELBD)
          GOTO 30
         ENDIF
C
         IF (IA.EQ.4) THEN
          IVT1=KVERT(3,-IELBD)
          IVT2=KVERT(4,-IELBD)
          IVT3=KVERT(8,-IELBD)
          IVT4=KVERT(7,-IELBD)
          GOTO 30
         ENDIF
C
         IF (IA.EQ.5) THEN
          IVT1=KVERT(4,-IELBD)
          IVT2=KVERT(1,-IELBD)
          IVT3=KVERT(5,-IELBD)
          IVT4=KVERT(8,-IELBD)
          GOTO 30
         ENDIF
C
         IF (IA.EQ.6) THEN
          IVT1=KVERT(5,-IELBD)
          IVT2=KVERT(6,-IELBD)
          IVT3=KVERT(7,-IELBD)
          IVT4=KVERT(8,-IELBD)
          GOTO 30
         ENDIF
C
30       P1X=DCORVG(1,IVT1)
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
         AX4=P4X-P1X
         AY4=P4Y-P1Y
         AZ4=P4Z-P1Z
C
         AX =PXC-PX
         AY =PYC-PY
         AZ =PZC-PZ
C
         DNX=(AY3*AZ2)-(AZ3*AY2)
         DNY=(AZ3*AX2)-(AX3*AZ2)
         DNZ=(AX3*AY2)-(AY3*AX2)
         DNAR1=SQRT(DNX*DNX+DNY*DNY+DNZ*DNZ)
C
         DNX=(AY4*AZ3)-(AZ4*AY3)
         DNY=(AZ4*AX3)-(AX4*AZ3)
         DNZ=(AX4*AY3)-(AY4*AX3)
         DNAR2=SQRT(DNX*DNX+DNY*DNY+DNZ*DNZ)
C
         DFAC=0.5D0*(DNAR1+DNAR2)/DNAR2
         DNX =DFAC*DNX
         DNY =DFAC*DNY
         DNZ =DFAC*DNZ
C
         DHN=DNX*AX+DNY*AY+DNZ*AZ
         IF (DHN.LT.0D0) THEN
          DNX=-DNX
          DNY=-DNY
          DNZ=-DNZ
         ENDIF      

         PDROP=FDATIN(8,0,PX,PY,PZ,TIMENS,RE)
C
         DF1(-IABD)=DF1(-IABD)+PDROP*DNX*TSTEPB
         DF2(-IABD)=DF2(-IABD)+PDROP*DNY*TSTEPB
         DF3(-IABD)=DF3(-IABD)+PDROP*DNZ*TSTEPB
        ENDIF
C
20      CONTINUE
C
       ENDIF
C
      ENDIF
C
10    CONTINUE
C
      END
C
************************************************************************
      SUBROUTINE    BDRDEF  (DX,KABD,NABD,A1)
************************************************************************
*    Purpose:  sets the NEUMANN-components of the vector DX to A1
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      DIMENSION DX(*)
      DIMENSION KABD(*)
C
C *** loop over all boundary vertices
C
      DO 1 IAT=1,NABD
      IABD=KABD(IAT)
      IF (IABD.GT.0) GOTO 1
      DX(-IABD)=A1*DX(-IABD)
1     CONTINUE
      END
C
C
************************************************************************
      SUBROUTINE BDRYA(VA,KCOL,KLD,KABD,NABD)
************************************************************************
*    Purpose:  updates the matrix entries for all DIRICHLET boundary
*              nodes
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL VA
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION VA(*),KCOL(*),KLD(*),KABD(*)
C
      DO 1 IAT=1,NABD
      IABD=KABD(IAT)
      IF (IABD.LT.0) GOTO 1
C
      VA(KLD(IABD))=1E0
      DO 2 ICOL=KLD(IABD)+1,KLD(IABD+1)-1
2     VA(ICOL)=0E0
C
1     CONTINUE
C
      END
C
************************************************************************
      SUBROUTINE BDRY0(D1,D2,D3,KABD,NABD)
************************************************************************
*    Purpose:  sets the DIRICHLET-components of the vector (D1,D2,D3) to
*              zero
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=8)
      DIMENSION D1(*),D2(*),D3(*)
      DIMENSION KABD(*)
C
      DO 1 IAT=1,NABD
      IABD=KABD(IAT)
      IF (IABD.LT.0) GOTO 1
C
      D1(IABD)=0D0
      D2(IABD)=0D0
      D3(IABD)=0D0
C
1     CONTINUE
C
      END
C
************************************************************************
      SUBROUTINE BDPRES(DP,KVERT,KAREA,KABD,KELBD,DCORVG,NABD,
     *                  PINT1,PINT2)
************************************************************************
*    Purpose:  Calculates integral boundary pressure
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      DIMENSION DP(*),KELBD(*),KABD(*),DCORVG(3,*),KVERT(8,*),KAREA(6,*)
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      SAVE
C
C
C
      DLEN1=0D0
      DLEN2=0D0
C
      PINT1=0D0
      PINT2=0D0
C
      DO 10 IAT=1,NABD
      IABD =KABD(IAT)
      IELBD=KELBD(IAT)
C
      IF ((IABD.LT.0).OR.(IELBD.LT.0)) GOTO 10
C
       IVT1=KVERT(1,IELBD)
       IVT2=KVERT(2,IELBD)
       IVT3=KVERT(3,IELBD)
       IVT4=KVERT(4,IELBD)
       IVT5=KVERT(5,IELBD)
       IVT6=KVERT(6,IELBD)
       IVT7=KVERT(7,IELBD)
       IVT8=KVERT(8,IELBD)
C
       XV1=DCORVG(1,IVT1)
       YV1=DCORVG(2,IVT1)
       ZV1=DCORVG(3,IVT1)
       XV2=DCORVG(1,IVT2)
       YV2=DCORVG(2,IVT2)
       ZV2=DCORVG(3,IVT2)
       XV3=DCORVG(1,IVT3)
       YV3=DCORVG(2,IVT3)
       ZV3=DCORVG(3,IVT3)
       XV4=DCORVG(1,IVT4)
       YV4=DCORVG(2,IVT4)
       ZV4=DCORVG(3,IVT4)
       XV5=DCORVG(1,IVT5)
       YV5=DCORVG(2,IVT5)
       ZV5=DCORVG(3,IVT5)
       XV6=DCORVG(1,IVT6)
       YV6=DCORVG(2,IVT6)
       ZV6=DCORVG(3,IVT6)
       XV7=DCORVG(1,IVT7)
       YV7=DCORVG(2,IVT7)
       ZV7=DCORVG(3,IVT7)
       XV8=DCORVG(1,IVT8)
       YV8=DCORVG(2,IVT8)
       ZV8=DCORVG(3,IVT8)
C
       IF (IABD.EQ.KAREA(1,IELBD)) THEN
        PX=(XV1+XV2+XV3+XV4)*0.25D0
        PY=(YV1+YV2+YV3+YV4)*0.25D0
        PZ=(ZV1+ZV2+ZV3+ZV4)*0.25D0
       ENDIF
C
       IF (IABD.EQ.KAREA(2,IELBD)) THEN
        PX=(XV1+XV2+XV6+XV5)*0.25D0
        PY=(YV1+YV2+YV6+YV5)*0.25D0
        PZ=(ZV1+ZV2+ZV6+ZV5)*0.25D0
       ENDIF
C
       IF (IABD.EQ.KAREA(3,IELBD)) THEN
        PX=(XV2+XV3+XV7+XV6)*0.25D0
        PY=(YV2+YV3+YV7+YV6)*0.25D0
        PZ=(ZV2+ZV3+ZV7+ZV6)*0.25D0
       ENDIF
C
       IF (IABD.EQ.KAREA(4,IELBD)) THEN
        PX=(XV3+XV4+XV8+XV7)*0.25D0
        PY=(YV3+YV4+YV8+YV7)*0.25D0
        PZ=(ZV3+ZV4+ZV8+ZV7)*0.25D0
       ENDIF
C
       IF (IABD.EQ.KAREA(5,IELBD)) THEN
        PX=(XV4+XV1+XV5+XV8)*0.25D0
        PY=(YV4+YV1+YV5+YV8)*0.25D0
        PZ=(ZV4+ZV1+ZV5+ZV8)*0.25D0
       ENDIF
C
       IF (IABD.EQ.KAREA(6,IELBD)) THEN
        PX=(XV5+XV6+XV7+XV8)*0.25D0
        PY=(YV5+YV6+YV7+YV8)*0.25D0
        PZ=(ZV5+ZV6+ZV7+ZV8)*0.25D0
       ENDIF
C
C
       IFLAG1=0
       IFLAG2=0
       CALL BDPDAT(IELBD,INPR,PX,PY,PZ,TIMENS,IFLAG1,IFLAG2)
C
       IF ((IFLAG1.EQ.1).OR.(IFLAG2.EQ.1)) THEN
C
        DO 20 IA=1,6
        IAH=KAREA(IA,IELBD)
C
        IF (IAH.EQ.IABD) THEN
C
         IF (IA.EQ.1) THEN
          IVT1=KVERT(1,IELBD)
          IVT2=KVERT(2,IELBD)
          IVT3=KVERT(3,IELBD)
          IVT4=KVERT(4,IELBD)
          GOTO 30
         ENDIF
C
         IF (IA.EQ.2) THEN
          IVT1=KVERT(1,IELBD)
          IVT2=KVERT(2,IELBD)
          IVT3=KVERT(6,IELBD)
          IVT4=KVERT(5,IELBD)
          GOTO 30
         ENDIF
C
         IF (IA.EQ.3) THEN
          IVT1=KVERT(2,IELBD)
          IVT2=KVERT(3,IELBD)
          IVT3=KVERT(7,IELBD)
          IVT4=KVERT(6,IELBD)
          GOTO 30
         ENDIF
C
         IF (IA.EQ.4) THEN
          IVT1=KVERT(3,IELBD)
          IVT2=KVERT(4,IELBD)
          IVT3=KVERT(8,IELBD)
          IVT4=KVERT(7,IELBD)
          GOTO 30
         ENDIF
C
         IF (IA.EQ.5) THEN
          IVT1=KVERT(4,IELBD)
          IVT2=KVERT(1,IELBD)
          IVT3=KVERT(5,IELBD)
          IVT4=KVERT(8,IELBD)
          GOTO 30
         ENDIF
C
         IF (IA.EQ.6) THEN
          IVT1=KVERT(5,IELBD)
          IVT2=KVERT(6,IELBD)
          IVT3=KVERT(7,IELBD)
          IVT4=KVERT(8,IELBD)
          GOTO 30
         ENDIF
C
30       P1X=DCORVG(1,IVT1)
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
         AX4=P4X-P1X
         AY4=P4Y-P1Y
         AZ4=P4Z-P1Z
C
         DNX=(AY3*AZ2)-(AZ3*AY2)
         DNY=(AZ3*AX2)-(AX3*AZ2)
         DNZ=(AX3*AY2)-(AY3*AX2)
         DNAR1=SQRT(DNX*DNX+DNY*DNY+DNZ*DNZ)
C
         DNX=(AY4*AZ3)-(AZ4*AY3)
         DNY=(AZ4*AX3)-(AX4*AZ3)
         DNZ=(AX4*AY3)-(AY4*AX3)
         DNAR2=SQRT(DNX*DNX+DNY*DNY+DNZ*DNZ)
C
         DAREAL=0.5D0*(DNAR1+DNAR2)
         PH=0.25D0*DAREAL*(DP(IVT1)+DP(IVT2)+DP(IVT3)+DP(IVT4))
C
C
         IF (IFLAG1.EQ.1) THEN
          DLEN1=DLEN1+DAREAL
          PINT1=PINT1+PH
         ENDIF
C         
         IF (IFLAG2.EQ.1) THEN
          DLEN2=DLEN2+DAREAL
          PINT2=PINT2+PH
         ENDIF
C
        ENDIF
C
20      CONTINUE
C
       ENDIF
C
10    CONTINUE
C
      IF (DLEN1.NE.0D0) PINT1=PINT1/DLEN1
      IF (DLEN2.NE.0D0) PINT2=PINT2/DLEN2
C
C
C
      END
C
************************************************************************
      SUBROUTINE BDFORC(DU1,DU2,DU3,DP,KVERT,KAREA,KEDGE,KABD,KELBD,
     *                  DCORVG,ELE,DFW,DAW)
************************************************************************
*    Purpose:  Calculates lift (DFW) and drag (DAW)
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNAB=21,NNDIM=3,NNCOF=10)
      DIMENSION DU1(*),DU2(*),DU3(*),DP(*),DCORVG(3,*)
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*),
     *          KABD(*),KELBD(*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DZ(NNVE),DJAC(3,3),DETJ,
     *                DBAS(NNDIM,NNBAS,NNDER),BDER(NNDER),KVE(NNVE),
     *                IEL,NDIM
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
C
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
      SAVE
C
C
C
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
      BDER(4)=.TRUE.
C
      IELTYP=-1
      NDIM=1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
C
      DFW  =0D0
      DAW  =0D0
C
      DO 10 IAT=1,NABD
      IABD =KABD(IAT)
      IELBD=KELBD(IAT)
C
      IF ((IABD.LT.0).OR.(IELBD.LT.0)) GOTO 10
C
       IVT1=KVERT(1,IELBD)
       IVT2=KVERT(2,IELBD)
       IVT3=KVERT(3,IELBD)
       IVT4=KVERT(4,IELBD)
       IVT5=KVERT(5,IELBD)
       IVT6=KVERT(6,IELBD)
       IVT7=KVERT(7,IELBD)
       IVT8=KVERT(8,IELBD)
C
       XV1=DCORVG(1,IVT1)
       YV1=DCORVG(2,IVT1)
       ZV1=DCORVG(3,IVT1)
       XV2=DCORVG(1,IVT2)
       YV2=DCORVG(2,IVT2)
       ZV2=DCORVG(3,IVT2)
       XV3=DCORVG(1,IVT3)
       YV3=DCORVG(2,IVT3)
       ZV3=DCORVG(3,IVT3)
       XV4=DCORVG(1,IVT4)
       YV4=DCORVG(2,IVT4)
       ZV4=DCORVG(3,IVT4)
       XV5=DCORVG(1,IVT5)
       YV5=DCORVG(2,IVT5)
       ZV5=DCORVG(3,IVT5)
       XV6=DCORVG(1,IVT6)
       YV6=DCORVG(2,IVT6)
       ZV6=DCORVG(3,IVT6)
       XV7=DCORVG(1,IVT7)
       YV7=DCORVG(2,IVT7)
       ZV7=DCORVG(3,IVT7)
       XV8=DCORVG(1,IVT8)
       YV8=DCORVG(2,IVT8)
       ZV8=DCORVG(3,IVT8)
C
C ***  PXC,PYC,PZC are coordinates of the center of the element
       PXC=0.125D0*(XV1+XV2+XV3+XV4+XV5+XV6+XV7+XV8)
       PYC=0.125D0*(YV1+YV2+YV3+YV4+YV5+YV6+YV7+YV8)
       PZC=0.125D0*(ZV1+ZV2+ZV3+ZV4+ZV5+ZV6+ZV7+ZV8)
C
C
       IF (IABD.EQ.KAREA(1,IELBD)) THEN
        PX=(XV1+XV2+XV3+XV4)*0.25D0
        PY=(YV1+YV2+YV3+YV4)*0.25D0
        PZ=(ZV1+ZV2+ZV3+ZV4)*0.25D0
       ENDIF
C
       IF (IABD.EQ.KAREA(2,IELBD)) THEN
        PX=(XV1+XV2+XV6+XV5)*0.25D0
        PY=(YV1+YV2+YV6+YV5)*0.25D0
        PZ=(ZV1+ZV2+ZV6+ZV5)*0.25D0
       ENDIF
C
       IF (IABD.EQ.KAREA(3,IELBD)) THEN
        PX=(XV2+XV3+XV7+XV6)*0.25D0
        PY=(YV2+YV3+YV7+YV6)*0.25D0
        PZ=(ZV2+ZV3+ZV7+ZV6)*0.25D0
       ENDIF
C
       IF (IABD.EQ.KAREA(4,IELBD)) THEN
        PX=(XV3+XV4+XV8+XV7)*0.25D0
        PY=(YV3+YV4+YV8+YV7)*0.25D0
        PZ=(ZV3+ZV4+ZV8+ZV7)*0.25D0
       ENDIF
C
       IF (IABD.EQ.KAREA(5,IELBD)) THEN
        PX=(XV4+XV1+XV5+XV8)*0.25D0
        PY=(YV4+YV1+YV5+YV8)*0.25D0
        PZ=(ZV4+ZV1+ZV5+ZV8)*0.25D0
       ENDIF
C
       IF (IABD.EQ.KAREA(6,IELBD)) THEN
        PX=(XV5+XV6+XV7+XV8)*0.25D0
        PY=(YV5+YV6+YV7+YV8)*0.25D0
        PZ=(ZV5+ZV6+ZV7+ZV8)*0.25D0
       ENDIF
C
C
       IFLAG=0
       CALL BDFDAT(IELBD,INPR,PX,PY,PZ,TIMENS,NY,IFLAG,DPF1,DPF2)
C
       IF (IFLAG.EQ.1) THEN
C
        DO 20 IA=1,6
        IAH=KAREA(IA,IELBD)
C
        IF (IAH.EQ.IABD) THEN
C
         IF (IA.EQ.1) THEN
          IVT1=KVERT(1,IELBD)
          IVT2=KVERT(2,IELBD)
          IVT3=KVERT(3,IELBD)
          IVT4=KVERT(4,IELBD)
          GOTO 30
         ENDIF
C
         IF (IA.EQ.2) THEN
          IVT1=KVERT(1,IELBD)
          IVT2=KVERT(2,IELBD)
          IVT3=KVERT(6,IELBD)
          IVT4=KVERT(5,IELBD)
          GOTO 30
         ENDIF
C
         IF (IA.EQ.3) THEN
          IVT1=KVERT(2,IELBD)
          IVT2=KVERT(3,IELBD)
          IVT3=KVERT(7,IELBD)
          IVT4=KVERT(6,IELBD)
          GOTO 30
         ENDIF
C
         IF (IA.EQ.4) THEN
          IVT1=KVERT(3,IELBD)
          IVT2=KVERT(4,IELBD)
          IVT3=KVERT(8,IELBD)
          IVT4=KVERT(7,IELBD)
          GOTO 30
         ENDIF
C
         IF (IA.EQ.5) THEN
          IVT1=KVERT(4,IELBD)
          IVT2=KVERT(1,IELBD)
          IVT3=KVERT(5,IELBD)
          IVT4=KVERT(8,IELBD)
          GOTO 30
         ENDIF
C
         IF (IA.EQ.6) THEN
          IVT1=KVERT(5,IELBD)
          IVT2=KVERT(6,IELBD)
          IVT3=KVERT(7,IELBD)
          IVT4=KVERT(8,IELBD)
          GOTO 30
         ENDIF
C
30       P1X=DCORVG(1,IVT1)
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
         AX4=P4X-P1X
         AY4=P4Y-P1Y
         AZ4=P4Z-P1Z
C
         AX =PXC-PX
         AY =PYC-PY
         AZ =PZC-PZ
C
         DNX=(AY3*AZ2)-(AZ3*AY2)
         DNY=(AZ3*AX2)-(AX3*AZ2)
         DNZ=(AX3*AY2)-(AY3*AX2)
         DNAR1=SQRT(DNX*DNX+DNY*DNY+DNZ*DNZ)
C
         DNX=(AY4*AZ3)-(AZ4*AY3)
         DNY=(AZ4*AX3)-(AX4*AZ3)
         DNZ=(AX4*AY3)-(AY4*AX3)
         DNAR2=SQRT(DNX*DNX+DNY*DNY+DNZ*DNZ)
C
         DNX =DNX/DNAR2
         DNY =DNY/DNAR2
         DNZ =DNZ/DNAR2
C
         DHN=DNX*AX+DNY*AY+DNZ*AZ
         IF (DHN.LT.0D0) THEN
          DNX=-DNX
          DNY=-DNY
          DNZ=-DNZ
         ENDIF
C
c         DNZ=0D0
C
         DTX= DNY
         DTY=-DNX     
         DTZ= DNZ     
C
         DAREAL=0.5D0*(DNAR1+DNAR2)
         DPCONT=0.25D0*(DP(IVT1)+DP(IVT2)+DP(IVT3)+DP(IVT4))
C
         IEL=IELBD
         CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
C
         DO 120 IVE=1,NVE
         JP=KVERT(IVE,IEL)
         KVE(IVE)=JP
         DX(IVE)=DCORVG(1,JP)
         DY(IVE)=DCORVG(2,JP)
         DZ(IVE)=DCORVG(3,JP)
120      CONTINUE
C
         XX=PX
         YY=PY
         ZZ=PZ
C
         CALL ELE(0D0,0D0,0D0,-2)
         CALL ELE(XX,YY,ZZ,-3)
C
         DUT=0
         DO 130 JDFL=1,IDFL
         DUT=DUT+DU1(KDFG(JDFL))*DBAS(1,KDFL(JDFL),2)*DTX*DNX
     *          +DU2(KDFG(JDFL))*DBAS(1,KDFL(JDFL),2)*DTY*DNX
     *          +DU1(KDFG(JDFL))*DBAS(1,KDFL(JDFL),3)*DTX*DNY
     *          +DU2(KDFG(JDFL))*DBAS(1,KDFL(JDFL),3)*DTY*DNY
130      CONTINUE
C
         DFW=DFW+DAREAL*(DPF1*DUT*DNY-DPCONT*DNX)
         DAW=DAW-DAREAL*(DPF1*DUT*DNX+DPCONT*DNY)
C
C
        ENDIF
C
20      CONTINUE
C
       ENDIF
C
10    CONTINUE
C
      IF (DPF2.NE.0D0) THEN
       DFW=2D0*DFW/DPF2      
       DAW=2D0*DAW/DPF2 
      ENDIF
C
C
C
99999 END
C

c#######################################################################
c#   new lift and drag evaluations .... 
c#######################################################################

************************************************************************
      subroutine lbdfnew(du1,du2,du3,dp,kvert,karea,kedge,kabd,kelbd,
     *    dcorvg,ele,dfwx,dfwy,dfwz)
************************************************************************
*       purpose:  calculates lift (dfw) and drag (daw)
*       general boundary integration
*       linear pressure
*-----------------------------------------------------------------------
      implicit double precision(a,c-h,o-u,w-z),logical(b)
c       
      parameter (nnbas=27,nnder=10,nncubp=36,nnve=8,nnee=12,nnae=6,
     *    nnab=21,nndim=3,nncof=10)
      dimension du1(*),du2(*),du3(*),dp(*),dcorvg(3,*)
      dimension kvert(nnve,*),karea(nnae,*),kedge(nnee,*)
      dimension kabd(*),kelbd(*)
      dimension kdfg(nnbas),kdfl(nnbas)
      character sub*6,fmt*15,cparam*120
c       
      common /output/ m,mt,mkeyb,mterm,merr,mprot,msys,mtrc,irecl8
      common /errctl/ ier,icheck
      common /char/   sub,fmt(3),cparam
      common /elem/   dx(nnve),dy(nnve),dz(nnve),djac(3,3),detj,
     *    dbas(nndim,nnbas,nnder),bder(nnder),kve(nnve),
     *    iel,ndim
      common /triad/  nel,nvt,net,nat,nve,nee,nae,nvel,neel,nved,
     *    nvar,near,nbct,nvbd,nebd,nabd
      common /cub/    dxi(nncubp,3),domega(nncubp),ncubp,icubp
      common /coaux1/ kdfg,kdfl,idfl
c       
      common /nspar/  tstep,theta,thstep,timens,epsns,nitns,itns,istat
c       
      double precision vrparm,ny
      dimension vrparm(100)
      equivalence (ny,vrparm)                          
      common /rparm/  ny,re,upsam,omgmin,omgmax,omgini,epsd,epsdiv,
     *    epsur,epspr,dmpd,dmpmg,epsmg,dmpsl,epssl,
     *    rlxsm,rlxsl,aminmg,amaxmg
      common /iparm/ iausav,ielt,istok,irhs,ibdr,ierana,imass,imassl,
     *    iupw,ipreca,iprecb,icubm,icuba,icubn,icubb,icubf,
     *    inlmin,inlmax,icyc,ilmin,ilmax,iint,ism,isl,
     *    nsm,nsl,nsmfac 
c       
      save
c       
c       
c       
      bder(1)=.true.
      bder(2)=.true.
      bder(3)=.true.
      bder(4)=.true.
c       
      ieltyp=-1
      ndim=1
      call ele(0d0,0d0,0d0,ieltyp)
      idfl=ndfl(ieltyp)
c       
      dfwx  =0d0
      dfwy  =0d0
      dfwz  =0d0
c       
      do iat=1,nabd
        iabd =kabd(iat)
        ielbd=kelbd(iat)
c       
        if ((iabd.le.0).or.(ielbd.le.0)) goto 10
c       
        ivt1=kvert(1,ielbd)
        ivt2=kvert(2,ielbd)
        ivt3=kvert(3,ielbd)
        ivt4=kvert(4,ielbd)
        ivt5=kvert(5,ielbd)
        ivt6=kvert(6,ielbd)
        ivt7=kvert(7,ielbd)
        ivt8=kvert(8,ielbd)
c       
        xv1=dcorvg(1,ivt1)
        yv1=dcorvg(2,ivt1)
        zv1=dcorvg(3,ivt1)
        xv2=dcorvg(1,ivt2)
        yv2=dcorvg(2,ivt2)
        zv2=dcorvg(3,ivt2)
        xv3=dcorvg(1,ivt3)
        yv3=dcorvg(2,ivt3)
        zv3=dcorvg(3,ivt3)
        xv4=dcorvg(1,ivt4)
        yv4=dcorvg(2,ivt4)
        zv4=dcorvg(3,ivt4)
        xv5=dcorvg(1,ivt5)
        yv5=dcorvg(2,ivt5)
        zv5=dcorvg(3,ivt5)
        xv6=dcorvg(1,ivt6)
        yv6=dcorvg(2,ivt6)
        zv6=dcorvg(3,ivt6)
        xv7=dcorvg(1,ivt7)
        yv7=dcorvg(2,ivt7)
        zv7=dcorvg(3,ivt7)
        xv8=dcorvg(1,ivt8)
        yv8=dcorvg(2,ivt8)
        zv8=dcorvg(3,ivt8)
c       
c       ***  pxc,pyc,pzc are coordinates of the center of the element
        pxc=0.125d0*(xv1+xv2+xv3+xv4+xv5+xv6+xv7+xv8)
        pyc=0.125d0*(yv1+yv2+yv3+yv4+yv5+yv6+yv7+yv8)
        pzc=0.125d0*(zv1+zv2+zv3+zv4+zv5+zv6+zv7+zv8)
c       
c       
        if (iabd.eq.karea(1,ielbd)) then
          ia=1
          px=(xv1+xv2+xv3+xv4)*0.25d0
          py=(yv1+yv2+yv3+yv4)*0.25d0
          pz=(zv1+zv2+zv3+zv4)*0.25d0
          ivt1=kvert(1,ielbd)
          ivt2=kvert(2,ielbd)
          ivt3=kvert(3,ielbd)
          ivt4=kvert(4,ielbd)
        endif
c       
        if (iabd.eq.karea(2,ielbd)) then
          ia=2
          px=(xv1+xv2+xv6+xv5)*0.25d0
          py=(yv1+yv2+yv6+yv5)*0.25d0
          pz=(zv1+zv2+zv6+zv5)*0.25d0
          ivt1=kvert(1,ielbd)
          ivt2=kvert(2,ielbd)
          ivt3=kvert(6,ielbd)
          ivt4=kvert(5,ielbd)
        endif
c       
        if (iabd.eq.karea(3,ielbd)) then
          ia=3
          px=(xv2+xv3+xv7+xv6)*0.25d0
          py=(yv2+yv3+yv7+yv6)*0.25d0
          pz=(zv2+zv3+zv7+zv6)*0.25d0
          ivt1=kvert(2,ielbd)
          ivt2=kvert(3,ielbd)
          ivt3=kvert(7,ielbd)
          ivt4=kvert(6,ielbd)
        endif
c       
        if (iabd.eq.karea(4,ielbd)) then
          ia=4
          px=(xv3+xv4+xv8+xv7)*0.25d0
          py=(yv3+yv4+yv8+yv7)*0.25d0
          pz=(zv3+zv4+zv8+zv7)*0.25d0
          ivt1=kvert(3,ielbd)
          ivt2=kvert(4,ielbd)
          ivt3=kvert(8,ielbd)
          ivt4=kvert(7,ielbd)
        endif
c       
        if (iabd.eq.karea(5,ielbd)) then
          ia=5
          px=(xv4+xv1+xv5+xv8)*0.25d0
          py=(yv4+yv1+yv5+yv8)*0.25d0
          pz=(zv4+zv1+zv5+zv8)*0.25d0
          ivt1=kvert(4,ielbd)
          ivt2=kvert(1,ielbd)
          ivt3=kvert(5,ielbd)
          ivt4=kvert(8,ielbd)
        endif
c       
        if (iabd.eq.karea(6,ielbd)) then
          ia=6
          px=(xv5+xv6+xv7+xv8)*0.25d0
          py=(yv5+yv6+yv7+yv8)*0.25d0
          pz=(zv5+zv6+zv7+zv8)*0.25d0
          ivt1=kvert(5,ielbd)
          ivt2=kvert(6,ielbd)
          ivt3=kvert(7,ielbd)
          ivt4=kvert(8,ielbd)
        endif
c       
c       
        iflag=0
        call bdfdat(ielbd,inpr,px,py,pz,timens,ny,iflag,dpf1,dpf2)
c       
        if (iflag.eq.1) then
c       
          iah=iabd
          
          p1x=dcorvg(1,ivt1)
          p1y=dcorvg(2,ivt1)
          p1z=dcorvg(3,ivt1)
          p2x=dcorvg(1,ivt2)
          p2y=dcorvg(2,ivt2)
          p2z=dcorvg(3,ivt2)
          p3x=dcorvg(1,ivt3)
          p3y=dcorvg(2,ivt3)
          p3z=dcorvg(3,ivt3)
          p4x=dcorvg(1,ivt4)
          p4y=dcorvg(2,ivt4)
          p4z=dcorvg(3,ivt4)
c       
          ax1=p2x-p1x
          ay1=p2y-p1y
          az1=p2z-p1z
          ax2=p3x-p2x
          ay2=p3y-p2y
          az2=p3z-p2z
          ax3=p4x-p3x
          ay3=p4y-p3y
          az3=p4z-p3z
          ax4=p1x-p4x
          ay4=p1y-p4y
          az4=p1z-p4z
c       
          ax =pxc-px
          ay =pyc-py
          az =pzc-pz
c       
          dn1x=(ay1*az4)-(az1*ay4)
          dn1y=(az1*ax4)-(ax1*az4)
          dn1z=(ax1*ay4)-(ay1*ax4)
          dn1ar=sqrt(dn1x*dn1x+dn1y*dn1y+dn1z*dn1z)
c
          dn2x=(ay2*az1)-(az2*ay1)
          dn2y=(az2*ax1)-(ax2*az1)
          dn2z=(ax2*ay1)-(ay2*ax1)
          dn2ar=sqrt(dn2x*dn2x+dn2y*dn2y+dn2z*dn2z)
c
          dn3x=(ay3*az2)-(az3*ay2)
          dn3y=(az3*ax2)-(ax3*az2)
          dn3z=(ax3*ay2)-(ay3*ax2)
          dn3ar=sqrt(dn3x*dn3x+dn3y*dn3y+dn3z*dn3z)
c
          dn4x=(ay4*az3)-(az4*ay3)
          dn4y=(az4*ax3)-(ax4*az3)
          dn4z=(ax4*ay3)-(ay4*ax3)
          dn4ar=sqrt(dn4x*dn4x+dn4y*dn4y+dn4z*dn4z)
c
          dnx=dn1x+dn2x+dn3x+dn4x
          dny=dn1y+dn2y+dn3y+dn4y
          dnz=dn1z+dn2z+dn3z+dn4z
          dnar=sqrt(dnx*dnx+dny*dny+dnz*dnz)
          dnx=dnx/dnar
          dny=dny/dnar
          dnz=dnz/dnar
c
          dhn=dnx*ax+dny*ay+dnz*az

          if (dhn.lt.0d0) then
            dnx=-dnx
            dny=-dny
            dnz=-dnz
          endif
c       
          dar=0.25d0*(dn1ar+dn2ar+dn3ar+dn4ar)
          dpc=0.25d0*(dp(ivt1)+dp(ivt2)+dp(ivt3)+dp(ivt4))
c       
          iel=ielbd
          call ndfgl(iel,1,ieltyp,kvert,kedge,karea,kdfg,kdfl)
c       
          do ive=1,nve
            jp=kvert(ive,iel)
            kve(ive)=jp
            dx(ive)=dcorvg(1,jp)
            dy(ive)=dcorvg(2,jp)
            dz(ive)=dcorvg(3,jp)
          enddo
c       
          xx=px
          yy=py
          zz=pz
c       
          call ele(0d0,0d0,0d0,-2)
          call ele(xx,yy,zz,-3)
c       
          du1x=0d0
          du1y=0d0
          du1z=0d0
          du2x=0d0
          du2y=0d0
          du2z=0d0
          du3x=0d0
          du3y=0d0
          du3z=0d0
          do jdfl=1,idfl
            dbx=dbas(1,kdfl(jdfl),2)
            dby=dbas(1,kdfl(jdfl),3)
            dbz=dbas(1,kdfl(jdfl),4)
            ig=kdfg(jdfl)
            du1x=du1x+du1(ig)*dbx
            du1y=du1y+du1(ig)*dby
            du1z=du1z+du1(ig)*dbz
            du2x=du2x+du2(ig)*dbx
            du2y=du2y+du2(ig)*dby
            du2z=du2z+du2(ig)*dbz
            du3x=du3x+du3(ig)*dbx
            du3y=du3y+du3(ig)*dby
            du3z=du3z+du3(ig)*dbz
          enddo
c       
          dfwx=dfwx+dar*(-dpc*dnx + dpf1*(du1x*dnx+du1y*dny+du1z*dnz))
          dfwy=dfwy+dar*(-dpc*dny + dpf1*(du2x*dnx+du2y*dny+du2z*dnz))
          dfwz=dfwz+dar*(-dpc*dnz + dpf1*(du3x*dnx+du3y*dny+du3z*dnz))
c       
        endif
c       
  10  continue
      enddo
c       
      if (dpf2.ne.0d0) then
        dfwx=2d0*dfwx/dpf2      
        dfwy=2d0*dfwy/dpf2      
        dfwz=2d0*dfwz/dpf2      
      endif
c       
99999 end
c
************************************************************************
      subroutine cbdfnew(du1,du2,du3,dp,kvert,karea,kedge,kabd,kelbd,
     *    dcorvg,ele,dfwx,dfwy,dfwz)
************************************************************************
*       purpose:  calculates lift (dfw) and drag (daw)
*       general boundary integration
*       constant pressure
*-----------------------------------------------------------------------
      implicit double precision(a,c-h,o-u,w-z),logical(b)
c       
      parameter (nnbas=27,nnder=10,nncubp=36,nnve=8,nnee=12,nnae=6,
     *    nnab=21,nndim=3,nncof=10)
      dimension du1(*),du2(*),du3(*),dp(*),dcorvg(3,*)
      dimension kvert(nnve,*),karea(nnae,*),kedge(nnee,*)
      dimension kabd(*),kelbd(*)
      dimension kdfg(nnbas),kdfl(nnbas)
      character sub*6,fmt*15,cparam*120
c       
      common /output/ m,mt,mkeyb,mterm,merr,mprot,msys,mtrc,irecl8
      common /errctl/ ier,icheck
      common /char/   sub,fmt(3),cparam
      common /elem/   dx(nnve),dy(nnve),dz(nnve),djac(3,3),detj,
     *    dbas(nndim,nnbas,nnder),bder(nnder),kve(nnve),
     *    iel,ndim
      common /triad/  nel,nvt,net,nat,nve,nee,nae,nvel,neel,nved,
     *    nvar,near,nbct,nvbd,nebd,nabd
      common /cub/    dxi(nncubp,3),domega(nncubp),ncubp,icubp
      common /coaux1/ kdfg,kdfl,idfl
c       
      common /nspar/  tstep,theta,thstep,timens,epsns,nitns,itns,istat
c       
      double precision vrparm,ny
      dimension vrparm(100)
      equivalence (ny,vrparm)                          
      common /rparm/  ny,re,upsam,omgmin,omgmax,omgini,epsd,epsdiv,
     *    epsur,epspr,dmpd,dmpmg,epsmg,dmpsl,epssl,
     *    rlxsm,rlxsl,aminmg,amaxmg
      common /iparm/ iausav,ielt,istok,irhs,ibdr,ierana,imass,imassl,
     *    iupw,ipreca,iprecb,icubm,icuba,icubn,icubb,icubf,
     *    inlmin,inlmax,icyc,ilmin,ilmax,iint,ism,isl,
     *    nsm,nsl,nsmfac 
c       
      save
c       
c       
c       
      bder(1)=.true.
      bder(2)=.true.
      bder(3)=.true.
      bder(4)=.true.
c       
      ieltyp=-1
      ndim=1
      call ele(0d0,0d0,0d0,ieltyp)
      idfl=ndfl(ieltyp)
c       
      dfwx  =0d0
      dfwy  =0d0
      dfwz  =0d0
c       
      do iat=1,nabd
        iabd =kabd(iat)
        ielbd=kelbd(iat)
c       
        if ((iabd.le.0).or.(ielbd.le.0)) goto 10
c       
        ivt1=kvert(1,ielbd)
        ivt2=kvert(2,ielbd)
        ivt3=kvert(3,ielbd)
        ivt4=kvert(4,ielbd)
        ivt5=kvert(5,ielbd)
        ivt6=kvert(6,ielbd)
        ivt7=kvert(7,ielbd)
        ivt8=kvert(8,ielbd)
c       
        xv1=dcorvg(1,ivt1)
        yv1=dcorvg(2,ivt1)
        zv1=dcorvg(3,ivt1)
        xv2=dcorvg(1,ivt2)
        yv2=dcorvg(2,ivt2)
        zv2=dcorvg(3,ivt2)
        xv3=dcorvg(1,ivt3)
        yv3=dcorvg(2,ivt3)
        zv3=dcorvg(3,ivt3)
        xv4=dcorvg(1,ivt4)
        yv4=dcorvg(2,ivt4)
        zv4=dcorvg(3,ivt4)
        xv5=dcorvg(1,ivt5)
        yv5=dcorvg(2,ivt5)
        zv5=dcorvg(3,ivt5)
        xv6=dcorvg(1,ivt6)
        yv6=dcorvg(2,ivt6)
        zv6=dcorvg(3,ivt6)
        xv7=dcorvg(1,ivt7)
        yv7=dcorvg(2,ivt7)
        zv7=dcorvg(3,ivt7)
        xv8=dcorvg(1,ivt8)
        yv8=dcorvg(2,ivt8)
        zv8=dcorvg(3,ivt8)
c       
c       ***  pxc,pyc,pzc are coordinates of the center of the element
        pxc=0.125d0*(xv1+xv2+xv3+xv4+xv5+xv6+xv7+xv8)
        pyc=0.125d0*(yv1+yv2+yv3+yv4+yv5+yv6+yv7+yv8)
        pzc=0.125d0*(zv1+zv2+zv3+zv4+zv5+zv6+zv7+zv8)
c       
c       
        if (iabd.eq.karea(1,ielbd)) then
          ia=1
          px=(xv1+xv2+xv3+xv4)*0.25d0
          py=(yv1+yv2+yv3+yv4)*0.25d0
          pz=(zv1+zv2+zv3+zv4)*0.25d0
          ivt1=kvert(1,ielbd)
          ivt2=kvert(2,ielbd)
          ivt3=kvert(3,ielbd)
          ivt4=kvert(4,ielbd)
        endif
c       
        if (iabd.eq.karea(2,ielbd)) then
          ia=2
          px=(xv1+xv2+xv6+xv5)*0.25d0
          py=(yv1+yv2+yv6+yv5)*0.25d0
          pz=(zv1+zv2+zv6+zv5)*0.25d0
          ivt1=kvert(1,ielbd)
          ivt2=kvert(2,ielbd)
          ivt3=kvert(6,ielbd)
          ivt4=kvert(5,ielbd)
        endif
c       
        if (iabd.eq.karea(3,ielbd)) then
          ia=3
          px=(xv2+xv3+xv7+xv6)*0.25d0
          py=(yv2+yv3+yv7+yv6)*0.25d0
          pz=(zv2+zv3+zv7+zv6)*0.25d0
          ivt1=kvert(2,ielbd)
          ivt2=kvert(3,ielbd)
          ivt3=kvert(7,ielbd)
          ivt4=kvert(6,ielbd)
        endif
c       
        if (iabd.eq.karea(4,ielbd)) then
          ia=4
          px=(xv3+xv4+xv8+xv7)*0.25d0
          py=(yv3+yv4+yv8+yv7)*0.25d0
          pz=(zv3+zv4+zv8+zv7)*0.25d0
          ivt1=kvert(3,ielbd)
          ivt2=kvert(4,ielbd)
          ivt3=kvert(8,ielbd)
          ivt4=kvert(7,ielbd)
        endif
c       
        if (iabd.eq.karea(5,ielbd)) then
          ia=5
          px=(xv4+xv1+xv5+xv8)*0.25d0
          py=(yv4+yv1+yv5+yv8)*0.25d0
          pz=(zv4+zv1+zv5+zv8)*0.25d0
          ivt1=kvert(4,ielbd)
          ivt2=kvert(1,ielbd)
          ivt3=kvert(5,ielbd)
          ivt4=kvert(8,ielbd)
        endif
c       
        if (iabd.eq.karea(6,ielbd)) then
          ia=6
          px=(xv5+xv6+xv7+xv8)*0.25d0
          py=(yv5+yv6+yv7+yv8)*0.25d0
          pz=(zv5+zv6+zv7+zv8)*0.25d0
          ivt1=kvert(5,ielbd)
          ivt2=kvert(6,ielbd)
          ivt3=kvert(7,ielbd)
          ivt4=kvert(8,ielbd)
        endif
c       
c       
        iflag=0
        call bdfdat(ielbd,inpr,px,py,pz,timens,ny,iflag,dpf1,dpf2)
c       
        if (iflag.eq.1) then
c       
          iah=iabd
          
          p1x=dcorvg(1,ivt1)
          p1y=dcorvg(2,ivt1)
          p1z=dcorvg(3,ivt1)
          p2x=dcorvg(1,ivt2)
          p2y=dcorvg(2,ivt2)
          p2z=dcorvg(3,ivt2)
          p3x=dcorvg(1,ivt3)
          p3y=dcorvg(2,ivt3)
          p3z=dcorvg(3,ivt3)
          p4x=dcorvg(1,ivt4)
          p4y=dcorvg(2,ivt4)
          p4z=dcorvg(3,ivt4)
c       
          ax1=p2x-p1x
          ay1=p2y-p1y
          az1=p2z-p1z
          ax2=p3x-p2x
          ay2=p3y-p2y
          az2=p3z-p2z
          ax3=p4x-p3x
          ay3=p4y-p3y
          az3=p4z-p3z
          ax4=p1x-p4x
          ay4=p1y-p4y
          az4=p1z-p4z
c       
          ax =pxc-px
          ay =pyc-py
          az =pzc-pz
c       
          dn1x=(ay1*az4)-(az1*ay4)
          dn1y=(az1*ax4)-(ax1*az4)
          dn1z=(ax1*ay4)-(ay1*ax4)
          dn1ar=sqrt(dn1x*dn1x+dn1y*dn1y+dn1z*dn1z)
c
          dn2x=(ay2*az1)-(az2*ay1)
          dn2y=(az2*ax1)-(ax2*az1)
          dn2z=(ax2*ay1)-(ay2*ax1)
          dn2ar=sqrt(dn2x*dn2x+dn2y*dn2y+dn2z*dn2z)
c
          dn3x=(ay3*az2)-(az3*ay2)
          dn3y=(az3*ax2)-(ax3*az2)
          dn3z=(ax3*ay2)-(ay3*ax2)
          dn3ar=sqrt(dn3x*dn3x+dn3y*dn3y+dn3z*dn3z)
c
          dn4x=(ay4*az3)-(az4*ay3)
          dn4y=(az4*ax3)-(ax4*az3)
          dn4z=(ax4*ay3)-(ay4*ax3)
          dn4ar=sqrt(dn4x*dn4x+dn4y*dn4y+dn4z*dn4z)
c
          dnx=dn1x+dn2x+dn3x+dn4x
          dny=dn1y+dn2y+dn3y+dn4y
          dnz=dn1z+dn2z+dn3z+dn4z
          dnar=sqrt(dnx*dnx+dny*dny+dnz*dnz)
          dnx=dnx/dnar
          dny=dny/dnar
          dnz=dnz/dnar
c
          dhn=dnx*ax+dny*ay+dnz*az
          if (dhn.lt.0d0) then
            dnx=-dnx
            dny=-dny
            dnz=-dnz
          endif
c       
          iel=ielbd

          dar=0.25d0*(dn1ar+dn2ar+dn3ar+dn4ar)
          dpc=dp(iel)
c       
          call ndfgl(iel,1,ieltyp,kvert,kedge,karea,kdfg,kdfl)
c       
          do ive=1,nve
            jp=kvert(ive,iel)
            kve(ive)=jp
            dx(ive)=dcorvg(1,jp)
            dy(ive)=dcorvg(2,jp)
            dz(ive)=dcorvg(3,jp)
          enddo
c       
          xx=px
          yy=py
          zz=pz
c       
          call ele(0d0,0d0,0d0,-2)
          call ele(xx,yy,zz,-3)
c       
          du1x=0d0
          du1y=0d0
          du1z=0d0
          du2x=0d0
          du2y=0d0
          du2z=0d0
          du3x=0d0
          du3y=0d0
          du3z=0d0
          do jdfl=1,idfl
            dbx=dbas(1,kdfl(jdfl),2)
            dby=dbas(1,kdfl(jdfl),3)
            dbz=dbas(1,kdfl(jdfl),4)
            ig=kdfg(jdfl)
            du1x=du1x+du1(ig)*dbx
            du1y=du1y+du1(ig)*dby
            du1z=du1z+du1(ig)*dbz
            du2x=du2x+du2(ig)*dbx
            du2y=du2y+du2(ig)*dby
            du2z=du2z+du2(ig)*dbz
            du3x=du3x+du3(ig)*dbx
            du3y=du3y+du3(ig)*dby
            du3z=du3z+du3(ig)*dbz
          enddo
c       
          dfwx=dfwx+dar*(-dpc*dnx + dpf1*(du1x*dnx+du1y*dny+du1z*dnz))
          dfwy=dfwy+dar*(-dpc*dny + dpf1*(du2x*dnx+du2y*dny+du2z*dnz))
          dfwz=dfwz+dar*(-dpc*dnz + dpf1*(du3x*dnx+du3y*dny+du3z*dnz))
c       
        endif
c       
  10  continue
      enddo
c       
      if (dpf2.ne.0d0) then
        dfwx=2d0*dfwx/dpf2      
        dfwy=2d0*dfwy/dpf2      
        dfwz=2d0*dfwz/dpf2      
      endif
c       
99999 end
c
************************************************************************
      subroutine bdfvol(du1,du2,du3,dp,kvert,karea,kedge,kabd,kelbd,
     *    dcorvg,ele,dfwx,dfwy,dfwz,da,ip)
************************************************************************
*       purpose: calculates lift (dfw) and drag (daw) 
*       volume integration
*       ip=0,1 ... for constant, linear pressure
*-----------------------------------------------------------------------
      implicit double precision (a,c-h,o-u,w-z),logical(b)
      character sub*6,fmt*15,cparam*120
c       
      parameter (nnbas=27,nnder=10,nncubp=36)
      parameter (nnve=8,nnee=12,nnae=6,nndim=3)
      dimension dcorvg(3,*),kvert(nnve,*),karea(nnae,*),kedge(nnee,*)
      dimension kabd(*),kelbd(*)
      dimension kdfg(nnbas),kdfl(nnbas)
      dimension kentry(nnbas,nnbas),dentry(nnbas,nnbas)
c
      dimension du1(*),du2(*),du3(*),dp(*),da(*)
c       
      common /output/ m,mt,mkeyb,mterm,merr,mprot,msys,mtrc,irecl8
      common /errctl/ ier,icheck
      common /char/   sub,fmt(3),cparam
      common /elem/   dx(nnve),dy(nnve),dz(nnve),djac(3,3),detj,
     *    dbas(nndim,nnbas,nnder),bder(nnder),kve(nnve),
     *    iel,ndim
      common /triad/  nel,nvt,net,nat,nve,nee,nae,nvel,neel,nved,
     *    nvar,near,nbct,nvbd,nebd,nabd
      common /cub/    dxi(nncubp,3),domega(nncubp),ncubp,icubp
      common /coaux1/ kdfg,kdfl,idfl

      double precision vrparm,ny
      dimension vrparm(100)
      equivalence (ny,vrparm)                          
      common /rparm/  ny,re,upsam,omgmin,omgmax,omgini,epsd,epsdiv,
     *    epsur,epspr,dmpd,dmpmg,epsmg,dmpsl,epssl,
     *    rlxsm,rlxsl,aminmg,amaxmg
      common /iparm/ iausav,ielt,istok,irhs,ibdr,ierana,imass,imassl,
     *    iupw,ipreca,iprecb,icubm,icuba,icubn,icubb,icubf,
     *    inlmin,inlmax,icyc,ilmin,ilmax,iint,ism,isl,
     *    nsm,nsl,nsmfac 

      save 
      external ele
c       
c       prepare vector dal such that it is 1 at the boundary where the
c       force should be computed and 0 otherwise 
      
      call  lcl1( da, nat) 
      
      do iat=1,nabd
        iabd =kabd(iat)
        ielbd=kelbd(iat)
c       
        if ((iabd.gt.0).and.(ielbd.gt.0)) then
c       
          ivt1=kvert(1,ielbd)
          ivt2=kvert(2,ielbd)
          ivt3=kvert(3,ielbd)
          ivt4=kvert(4,ielbd)
          ivt5=kvert(5,ielbd)
          ivt6=kvert(6,ielbd)
          ivt7=kvert(7,ielbd)
          ivt8=kvert(8,ielbd)
c       
          xv1=dcorvg(1,ivt1)
          yv1=dcorvg(2,ivt1)
          zv1=dcorvg(3,ivt1)
          xv2=dcorvg(1,ivt2)
          yv2=dcorvg(2,ivt2)
          zv2=dcorvg(3,ivt2)
          xv3=dcorvg(1,ivt3)
          yv3=dcorvg(2,ivt3)
          zv3=dcorvg(3,ivt3)
          xv4=dcorvg(1,ivt4)
          yv4=dcorvg(2,ivt4)
          zv4=dcorvg(3,ivt4)
          xv5=dcorvg(1,ivt5)
          yv5=dcorvg(2,ivt5)
          zv5=dcorvg(3,ivt5)
          xv6=dcorvg(1,ivt6)
          yv6=dcorvg(2,ivt6)
          zv6=dcorvg(3,ivt6)
          xv7=dcorvg(1,ivt7)
          yv7=dcorvg(2,ivt7)
          zv7=dcorvg(3,ivt7)
          xv8=dcorvg(1,ivt8)
          yv8=dcorvg(2,ivt8)
          zv8=dcorvg(3,ivt8)
c       
c       ***  pxc,pyc,pzc are coordinates of the center of the element
          pxc=0.125d0*(xv1+xv2+xv3+xv4+xv5+xv6+xv7+xv8)
          pyc=0.125d0*(yv1+yv2+yv3+yv4+yv5+yv6+yv7+yv8)
          pzc=0.125d0*(zv1+zv2+zv3+zv4+zv5+zv6+zv7+zv8)
c       
c       
          if (iabd.eq.karea(1,ielbd)) then
            ia=1
            px=(xv1+xv2+xv3+xv4)*0.25d0
            py=(yv1+yv2+yv3+yv4)*0.25d0
            pz=(zv1+zv2+zv3+zv4)*0.25d0
            ivt1=kvert(1,ielbd)
            ivt2=kvert(2,ielbd)
            ivt3=kvert(3,ielbd)
            ivt4=kvert(4,ielbd)
          endif
c       
          if (iabd.eq.karea(2,ielbd)) then
            ia=2
            px=(xv1+xv2+xv6+xv5)*0.25d0
            py=(yv1+yv2+yv6+yv5)*0.25d0
            pz=(zv1+zv2+zv6+zv5)*0.25d0
            ivt1=kvert(1,ielbd)
            ivt2=kvert(2,ielbd)
            ivt3=kvert(6,ielbd)
            ivt4=kvert(5,ielbd)
          endif
c       
          if (iabd.eq.karea(3,ielbd)) then
            ia=3
            px=(xv2+xv3+xv7+xv6)*0.25d0
            py=(yv2+yv3+yv7+yv6)*0.25d0
            pz=(zv2+zv3+zv7+zv6)*0.25d0
            ivt1=kvert(2,ielbd)
            ivt2=kvert(3,ielbd)
            ivt3=kvert(7,ielbd)
            ivt4=kvert(6,ielbd)
          endif
c       
          if (iabd.eq.karea(4,ielbd)) then
            ia=4
            px=(xv3+xv4+xv8+xv7)*0.25d0
            py=(yv3+yv4+yv8+yv7)*0.25d0
            pz=(zv3+zv4+zv8+zv7)*0.25d0
            ivt1=kvert(3,ielbd)
            ivt2=kvert(4,ielbd)
            ivt3=kvert(8,ielbd)
            ivt4=kvert(7,ielbd)
          endif
c       
          if (iabd.eq.karea(5,ielbd)) then
            ia=5
            px=(xv4+xv1+xv5+xv8)*0.25d0
            py=(yv4+yv1+yv5+yv8)*0.25d0
            pz=(zv4+zv1+zv5+zv8)*0.25d0
            ivt1=kvert(4,ielbd)
            ivt2=kvert(1,ielbd)
            ivt3=kvert(5,ielbd)
            ivt4=kvert(8,ielbd)
          endif
c       
          if (iabd.eq.karea(6,ielbd)) then
            ia=6
            px=(xv5+xv6+xv7+xv8)*0.25d0
            py=(yv5+yv6+yv7+yv8)*0.25d0
            pz=(zv5+zv6+zv7+zv8)*0.25d0
            ivt1=kvert(5,ielbd)
            ivt2=kvert(6,ielbd)
            ivt3=kvert(7,ielbd)
            ivt4=kvert(8,ielbd)
          endif
c       
          iflag=0
          call bdfdat(ielbd,inpr,px,py,pz,timens,ny,iflag,dpf1,dpf2)
c       
          if (iflag.eq.1) then
            da(iabd)=1.0d0
          endif
          
        endif
c       
      enddo
      
c       use the vector dal as a test function in the weak formulation
c       (only the difusion part is taken, should take even the
c       convection but in featflow it seems to be ok) ... 
c       this needs to be checked in timedependent case
      
      dfwx=0d0
      dfwy=0d0
      dfwz=0d0

c       ip=0 ... constant pressure , ip=1 ... linear pressure
      if(ip.eq.0) then
        call cstgforce(du1,du2,du3,dp,kvert,karea,kedge,kabd,kelbd,
     *      dcorvg,ele,dfwx,dfwy,dfwz,da,ielt,9,dpf1,dpf2)
      else
        call lstgforce(du1,du2,du3,dp,kvert,karea,kedge,kabd,kelbd,
     *      dcorvg,ele,dfwx,dfwy,dfwz,da,ielt,9,dpf1,dpf2)
      endif        
      if (dpf2.ne.0d0) then
        dfwx=2d0*dfwx/dpf2
        dfwy=2d0*dfwy/dpf2
        dfwz=2d0*dfwz/dpf2
      endif      
      
99999 end
      




c=============================================================================
c  the volume integration with constant pressure dp(1..nel)

      subroutine cstgforce(du1,du2,du3,dp,kvert,karea,kedge,kabd,kelbd,
     *    dcorvg,ele,dfwx,dfwy,dfwz,da,ielt,icub,dpf1,dpf2)
c       
      implicit double precision (a,c-h,o-u,w-z),logical(b)
      character sub*6,fmt*15,cparam*120
c       
      parameter (nnbas=27,nnder=10,nncubp=36)
      parameter (nnve=8,nnee=12,nnae=6,nndim=3)
      dimension dcorvg(3,*),kvert(nnve,*),karea(nnae,*),kedge(nnee,*)
      dimension kabd(*),kelbd(*)
      dimension kdfg(nnbas),kdfl(nnbas)
      dimension kentry(nnbas,nnbas),dentry(nnbas,nnbas)
c
      dimension du1(*),du2(*),du3(*),dp(*),da(*)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
c       
      common /output/ m,mt,mkeyb,mterm,merr,mprot,msys,mtrc,irecl8
      common /errctl/ ier,icheck
      common /char/   sub,fmt(3),cparam
      common /elem/   dx(nnve),dy(nnve),dz(nnve),djac(3,3),detj,
     *    dbas(nndim,nnbas,nnder),bder(nnder),kve(nnve),
     *    iel,ndim
      common /triad/  nel,nvt,net,nat,nve,nee,nae,nvel,neel,nved,
     *    nvar,near,nbct,nvbd,nebd,nabd
      common /cub/    dxi(nncubp,3),domega(nncubp),ncubp,icubp
      common /coaux1/ kdfg,kdfl,idfl

      save 
      external ele
c       

c       *** preparation - evaluation of parameters
      ier=0

c       *** which derivatives of basis functions are needed?
      do  i = 1,nnder
        bder(i)=.false.
      enddo
      
      bder(1)=.true.
      bder(2)=.true.
      bder(3)=.true.
      bder(4)=.true.
      
c       *** dummy call of ele sets number of element
      ieltyp=-1
      call ele(0d0,0d0,0d0,ieltyp)

      idfl=ndfl(ieltyp)
      call cb3h(icub)
      if (ier.ne.0) goto 99999
c       
c       *** dummy call - ele may save arithmetic operations
      icubp=icub
      call ele(0d0,0d0,0d0,-2)
      
c       *** loop over all elements
      do iel=1,nel
        call ndfgl(iel,1,ieltyp,kvert,kedge,karea,kdfg,kdfl)
        if (ier.lt.0) goto 99999
c       
c       *** evaluation of coordinates of the vertices
        do ive = 1, nve
          jp=kvert(ive,iel)
          kve(ive)=jp
          dx(ive)=dcorvg(1,jp)
          dy(ive)=dcorvg(2,jp)
          dz(ive)=dcorvg(3,jp)
        enddo
c       
       dj11=( dx(1)+dx(2)+dx(3)+dx(4)+dx(5)+dx(6)+dx(7)+dx(8))*q8
       dj12=( dy(1)+dy(2)+dy(3)+dy(4)+dy(5)+dy(6)+dy(7)+dy(8))*q8
       dj13=( dz(1)+dz(2)+dz(3)+dz(4)+dz(5)+dz(6)+dz(7)+dz(8))*q8
       dj21=(-dx(1)+dx(2)+dx(3)-dx(4)-dx(5)+dx(6)+dx(7)-dx(8))*q8
       dj22=(-dy(1)+dy(2)+dy(3)-dy(4)-dy(5)+dy(6)+dy(7)-dy(8))*q8
       dj23=(-dz(1)+dz(2)+dz(3)-dz(4)-dz(5)+dz(6)+dz(7)-dz(8))*q8
       dj31=(-dx(1)-dx(2)+dx(3)+dx(4)-dx(5)-dx(6)+dx(7)+dx(8))*q8
       dj32=(-dy(1)-dy(2)+dy(3)+dy(4)-dy(5)-dy(6)+dy(7)+dy(8))*q8
       dj33=(-dz(1)-dz(2)+dz(3)+dz(4)-dz(5)-dz(6)+dz(7)+dz(8))*q8
       dj41=(-dx(1)-dx(2)-dx(3)-dx(4)+dx(5)+dx(6)+dx(7)+dx(8))*q8
       dj42=(-dy(1)-dy(2)-dy(3)-dy(4)+dy(5)+dy(6)+dy(7)+dy(8))*q8
       dj43=(-dz(1)-dz(2)-dz(3)-dz(4)+dz(5)+dz(6)+dz(7)+dz(8))*q8
       dj51=( dx(1)-dx(2)+dx(3)-dx(4)+dx(5)-dx(6)+dx(7)-dx(8))*q8
       dj52=( dy(1)-dy(2)+dy(3)-dy(4)+dy(5)-dy(6)+dy(7)-dy(8))*q8
       dj53=( dz(1)-dz(2)+dz(3)-dz(4)+dz(5)-dz(6)+dz(7)-dz(8))*q8
       dj61=( dx(1)-dx(2)-dx(3)+dx(4)-dx(5)+dx(6)+dx(7)-dx(8))*q8
       dj62=( dy(1)-dy(2)-dy(3)+dy(4)-dy(5)+dy(6)+dy(7)-dy(8))*q8
       dj63=( dz(1)-dz(2)-dz(3)+dz(4)-dz(5)+dz(6)+dz(7)-dz(8))*q8
       dj71=( dx(1)+dx(2)-dx(3)-dx(4)-dx(5)-dx(6)+dx(7)+dx(8))*q8
       dj72=( dy(1)+dy(2)-dy(3)-dy(4)-dy(5)-dy(6)+dy(7)+dy(8))*q8
       dj73=( dz(1)+dz(2)-dz(3)-dz(4)-dz(5)-dz(6)+dz(7)+dz(8))*q8
       dj81=(-dx(1)+dx(2)-dx(3)+dx(4)+dx(5)-dx(6)+dx(7)-dx(8))*q8
       dj82=(-dy(1)+dy(2)-dy(3)+dy(4)+dy(5)-dy(6)+dy(7)-dy(8))*q8
       dj83=(-dz(1)+dz(2)-dz(3)+dz(4)+dz(5)-dz(6)+dz(7)-dz(8))*q8
c       
        call ele(0d0,0d0,0d0,-2)
c       *** loop over all cubature points
        do icubp = 1, ncubp
          
c       *** cubature point on the reference element
          xi1=dxi(icubp,1)
          xi2=dxi(icubp,2)
          xi3=dxi(icubp,3)
          
c       *** jacobian of the bilinear mapping onto the reference element
       djac(1,1)=dj21+dj51*xi2+dj61*xi3+dj81*xi2*xi3
       djac(1,2)=dj31+dj51*xi1+dj71*xi3+dj81*xi1*xi3
       djac(1,3)=dj41+dj61*xi1+dj71*xi2+dj81*xi1*xi2
       djac(2,1)=dj22+dj52*xi2+dj62*xi3+dj82*xi2*xi3
       djac(2,2)=dj32+dj52*xi1+dj72*xi3+dj82*xi1*xi3
       djac(2,3)=dj42+dj62*xi1+dj72*xi2+dj82*xi1*xi2
       djac(3,1)=dj23+dj53*xi2+dj63*xi3+dj83*xi2*xi3
       djac(3,2)=dj33+dj53*xi1+dj73*xi3+dj83*xi1*xi3
       djac(3,3)=dj43+dj63*xi1+dj73*xi2+dj83*xi1*xi2
       detj= djac(1,1)*(djac(2,2)*djac(3,3)-djac(3,2)*djac(2,3))
     *      -djac(2,1)*(djac(1,2)*djac(3,3)-djac(3,2)*djac(1,3))
     *      +djac(3,1)*(djac(1,2)*djac(2,3)-djac(2,2)*djac(1,3))
       om=domega(icubp)*abs(detj)
          
c       *** cubature point on the real element
       xx=dj11+djac(1,1)*xi1+dj31*xi2+dj41*xi3+dj71*xi2*xi3
       yy=dj12+dj22*xi1+djac(2,2)*xi2+dj42*xi3+dj62*xi1*xi3
       zz=dj13+dj23*xi1+dj33*xi2+djac(3,3)*xi3+dj53*xi1*xi2
          
c       *** evaluate the basis functions in the cubature point
          if((ielt.eq.2).or.(ielt.eq.3)) then
            call ele(xx,yy,zz,-3)
          else
            call ele(xi1,xi2,xi3,-3)
          endif
          if (ier.lt.0) goto 99999
          
c       evaluate the solution values and derivatives in the cubature point     
          du1v=0d0 ! value
          du1x=0d0 ! x dreiv.
          du1y=0d0 ! y deriv
          du1z=0d0 ! z deriv
          do i=1,idfl
            ig=kdfg(i)
            du1v=du1v+du1(ig)*dbas(1,kdfl(i),1)
            du1x=du1x+du1(ig)*dbas(1,kdfl(i),2)
            du1y=du1y+du1(ig)*dbas(1,kdfl(i),3)
            du1z=du1z+du1(ig)*dbas(1,kdfl(i),4)
          enddo
          
          du2v=0d0 ! value
          du2x=0d0 ! x dreiv.
          du2y=0d0 ! y deriv
          du2z=0d0 ! y deriv
          do i=1,idfl
            ig=kdfg(i)
            du2v=du2v+du2(ig)*dbas(1,kdfl(i),1)
            du2x=du2x+du2(ig)*dbas(1,kdfl(i),2)
            du2y=du2y+du2(ig)*dbas(1,kdfl(i),3)
            du2z=du2z+du2(ig)*dbas(1,kdfl(i),4)
          enddo
          
          du3v=0d0 ! value
          du3x=0d0 ! x dreiv.
          du3y=0d0 ! y deriv
          du3z=0d0 ! y deriv
          do i=1,idfl
            ig=kdfg(i)
            du3v=du3v+du3(ig)*dbas(1,kdfl(i),1)
            du3x=du3x+du3(ig)*dbas(1,kdfl(i),2)
            du3y=du3y+du3(ig)*dbas(1,kdfl(i),3)
            du3z=du3z+du3(ig)*dbas(1,kdfl(i),4)
          enddo
          
          dav=0d0 ! value
          dax=0d0 ! x dreiv.
          day=0d0 ! y deriv
          daz=0d0 ! y deriv
          do i=1,idfl
            ig=kdfg(i)
            dav=dav+da(ig)*dbas(1,kdfl(i),1)
            dax=dax+da(ig)*dbas(1,kdfl(i),2)
            day=day+da(ig)*dbas(1,kdfl(i),3)
            daz=daz+da(ig)*dbas(1,kdfl(i),4)
          enddo
c       form the integrand

          dnx=-dax
          dny=-day
          dnz=-daz
c
          ah1=-dp(iel)*dnx+dpf1*(du1x*dnx+du1y*dny+du1z*dnz)
          ah2=-dp(iel)*dny+dpf1*(du2x*dnx+du2y*dny+du2z*dnz)
          ah3=-dp(iel)*dnz+dpf1*(du3x*dnx+du3y*dny+du3z*dnz)

          dfwx=dfwx+ah1*om
          dfwy=dfwy+ah2*om
          dfwz=dfwz+ah3*om
          
        enddo
      enddo
c       
c       
99999 end

c=============================================================================
c  the volume integration with linear pressure dp(1..nvt)

      subroutine lstgforce(du1,du2,du3,dp,kvert,karea,kedge,kabd,kelbd,
     *    dcorvg,ele,dfwx,dfwy,dfwz,da,ielt,icub,dpf1,dpf2)
c       
      implicit double precision (a,c-h,o-u,w-z),logical(b)
      character sub*6,fmt*15,cparam*120
c       
      parameter (nnbas=27,nnder=10,nncubp=36)
      parameter (nnve=8,nnee=12,nnae=6,nndim=3)
      dimension dcorvg(3,*),kvert(nnve,*),karea(nnae,*),kedge(nnee,*)
      dimension kabd(*),kelbd(*)
c
      dimension du1(*),du2(*),du3(*),dp(*),da(*)

      PARAMETER (Q2=0.5D0,Q8=0.125D0)
c       
      common /output/ m,mt,mkeyb,mterm,merr,mprot,msys,mtrc,irecl8
      common /errctl/ ier,icheck
      common /char/   sub,fmt(3),cparam
      common /elem/   dx(nnve),dy(nnve),dz(nnve),djac(3,3),detj,
     *    dbas(nndim,nnbas,nnder),bder(nnder),kve(nnve),
     *    iel,ndim
      common /triad/  nel,nvt,net,nat,nve,nee,nae,nvel,neel,nved,
     *    nvar,near,nbct,nvbd,nebd,nabd
      common /cub/    dxi(nncubp,3),domega(nncubp),ncubp,icubp
      common /coaux1/ kdfg(nnbas),kdfl(nnbas),idfl
      common /coaux2/ kdfg1(nnbas),kdfl1(nnbas),idfl1

      save 
      external ele,e011
c       

c       *** preparation - evaluation of parameters
      ier=0

c       *** which derivatives of basis functions are needed?
      do  i = 1,nnder
        bder(i)=.false.
      enddo
      
      bder(1)=.true.
      bder(2)=.true.
      bder(3)=.true.
      bder(4)=.true.
      
c       *** dummy call of ele sets number of element
      ityp=-1
      call ele(0d0,0d0,0d0,ityp)
      idfl=ndfl(ityp)

      ityp1=-1
      call e011(0d0,0d0,0d0,ityp1)
      idfl1=ndfl(ityp1)

      call cb3h(icub)
      if (ier.ne.0) goto 99999
c       
c       
c       *** loop over all elements
      do iel=1,nel

        call ndfgl(iel,1,ityp,kvert,kedge,karea,kdfg,kdfl)
        if (ier.lt.0) goto 99999

        call ndfgl(iel,1,ityp1,kvert,kedge,karea,kdfg1,kdfl1)
        if (ier.lt.0) goto 99999

c       
c       *** evaluation of coordinates of the vertices
        do ive = 1, nve
          jp=kvert(ive,iel)
          kve(ive)=jp
          dx(ive)=dcorvg(1,jp)
          dy(ive)=dcorvg(2,jp)
          dz(ive)=dcorvg(3,jp)
        enddo
c       
       dj11=( dx(1)+dx(2)+dx(3)+dx(4)+dx(5)+dx(6)+dx(7)+dx(8))*q8
       dj12=( dy(1)+dy(2)+dy(3)+dy(4)+dy(5)+dy(6)+dy(7)+dy(8))*q8
       dj13=( dz(1)+dz(2)+dz(3)+dz(4)+dz(5)+dz(6)+dz(7)+dz(8))*q8
       dj21=(-dx(1)+dx(2)+dx(3)-dx(4)-dx(5)+dx(6)+dx(7)-dx(8))*q8
       dj22=(-dy(1)+dy(2)+dy(3)-dy(4)-dy(5)+dy(6)+dy(7)-dy(8))*q8
       dj23=(-dz(1)+dz(2)+dz(3)-dz(4)-dz(5)+dz(6)+dz(7)-dz(8))*q8
       dj31=(-dx(1)-dx(2)+dx(3)+dx(4)-dx(5)-dx(6)+dx(7)+dx(8))*q8
       dj32=(-dy(1)-dy(2)+dy(3)+dy(4)-dy(5)-dy(6)+dy(7)+dy(8))*q8
       dj33=(-dz(1)-dz(2)+dz(3)+dz(4)-dz(5)-dz(6)+dz(7)+dz(8))*q8
       dj41=(-dx(1)-dx(2)-dx(3)-dx(4)+dx(5)+dx(6)+dx(7)+dx(8))*q8
       dj42=(-dy(1)-dy(2)-dy(3)-dy(4)+dy(5)+dy(6)+dy(7)+dy(8))*q8
       dj43=(-dz(1)-dz(2)-dz(3)-dz(4)+dz(5)+dz(6)+dz(7)+dz(8))*q8
       dj51=( dx(1)-dx(2)+dx(3)-dx(4)+dx(5)-dx(6)+dx(7)-dx(8))*q8
       dj52=( dy(1)-dy(2)+dy(3)-dy(4)+dy(5)-dy(6)+dy(7)-dy(8))*q8
       dj53=( dz(1)-dz(2)+dz(3)-dz(4)+dz(5)-dz(6)+dz(7)-dz(8))*q8
       dj61=( dx(1)-dx(2)-dx(3)+dx(4)-dx(5)+dx(6)+dx(7)-dx(8))*q8
       dj62=( dy(1)-dy(2)-dy(3)+dy(4)-dy(5)+dy(6)+dy(7)-dy(8))*q8
       dj63=( dz(1)-dz(2)-dz(3)+dz(4)-dz(5)+dz(6)+dz(7)-dz(8))*q8
       dj71=( dx(1)+dx(2)-dx(3)-dx(4)-dx(5)-dx(6)+dx(7)+dx(8))*q8
       dj72=( dy(1)+dy(2)-dy(3)-dy(4)-dy(5)-dy(6)+dy(7)+dy(8))*q8
       dj73=( dz(1)+dz(2)-dz(3)-dz(4)-dz(5)-dz(6)+dz(7)+dz(8))*q8
       dj81=(-dx(1)+dx(2)-dx(3)+dx(4)+dx(5)-dx(6)+dx(7)-dx(8))*q8
       dj82=(-dy(1)+dy(2)-dy(3)+dy(4)+dy(5)-dy(6)+dy(7)-dy(8))*q8
       dj83=(-dz(1)+dz(2)-dz(3)+dz(4)+dz(5)-dz(6)+dz(7)-dz(8))*q8
c       
        call ele(0d0,0d0,0d0,-2)
        call e011(0d0,0d0,0d0,-2)

c       *** loop over all cubature points
        do icubp = 1, ncubp
          
c       *** cubature point on the reference element
          xi1=dxi(icubp,1)
          xi2=dxi(icubp,2)
          xi3=dxi(icubp,3)
          
c       *** jacobian of the bilinear mapping onto the reference element
       djac(1,1)=dj21+dj51*xi2+dj61*xi3+dj81*xi2*xi3
       djac(1,2)=dj31+dj51*xi1+dj71*xi3+dj81*xi1*xi3
       djac(1,3)=dj41+dj61*xi1+dj71*xi2+dj81*xi1*xi2
       djac(2,1)=dj22+dj52*xi2+dj62*xi3+dj82*xi2*xi3
       djac(2,2)=dj32+dj52*xi1+dj72*xi3+dj82*xi1*xi3
       djac(2,3)=dj42+dj62*xi1+dj72*xi2+dj82*xi1*xi2
       djac(3,1)=dj23+dj53*xi2+dj63*xi3+dj83*xi2*xi3
       djac(3,2)=dj33+dj53*xi1+dj73*xi3+dj83*xi1*xi3
       djac(3,3)=dj43+dj63*xi1+dj73*xi2+dj83*xi1*xi2
       detj= djac(1,1)*(djac(2,2)*djac(3,3)-djac(3,2)*djac(2,3))
     *      -djac(2,1)*(djac(1,2)*djac(3,3)-djac(3,2)*djac(1,3))
     *      +djac(3,1)*(djac(1,2)*djac(2,3)-djac(2,2)*djac(1,3))
       om=domega(icubp)*abs(detj)
          
c       *** cubature point on the real element
       xx=dj11+djac(1,1)*xi1+dj31*xi2+dj41*xi3+dj71*xi2*xi3
       yy=dj12+dj22*xi1+djac(2,2)*xi2+dj42*xi3+dj62*xi1*xi3
       zz=dj13+dj23*xi1+dj33*xi2+djac(3,3)*xi3+dj53*xi1*xi2
          
c       *** evaluate the basis functions in the cubature point
c       for the velocities
          if((ielt.eq.2).or.(ielt.eq.3)) then
            call ele(xx,yy,zz,-3)
          else
            call ele(xi1,xi2,xi3,-3)
          endif
          if (ier.lt.0) goto 99999
 
c       evaluate the solution values and derivatives in the cubature point     
          du1v=0d0 ! value
          du1x=0d0 ! x dreiv.
          du1y=0d0 ! y deriv
          du1z=0d0 ! z deriv
          do i=1,idfl
            ig=kdfg(i)
            du1v=du1v+du1(ig)*dbas(1,kdfl(i),1)
            du1x=du1x+du1(ig)*dbas(1,kdfl(i),2)
            du1y=du1y+du1(ig)*dbas(1,kdfl(i),3)
            du1z=du1z+du1(ig)*dbas(1,kdfl(i),4)
          enddo
          
          du2v=0d0 ! value
          du2x=0d0 ! x dreiv.
          du2y=0d0 ! y deriv
          du2z=0d0 ! y deriv
          do i=1,idfl
            ig=kdfg(i)
            du2v=du2v+du2(ig)*dbas(1,kdfl(i),1)
            du2x=du2x+du2(ig)*dbas(1,kdfl(i),2)
            du2y=du2y+du2(ig)*dbas(1,kdfl(i),3)
            du2z=du2z+du2(ig)*dbas(1,kdfl(i),4)
          enddo
          
          du3v=0d0 ! value
          du3x=0d0 ! x dreiv.
          du3y=0d0 ! y deriv
          du3z=0d0 ! y deriv
          do i=1,idfl
            ig=kdfg(i)
            du3v=du3v+du3(ig)*dbas(1,kdfl(i),1)
            du3x=du3x+du3(ig)*dbas(1,kdfl(i),2)
            du3y=du3y+du3(ig)*dbas(1,kdfl(i),3)
            du3z=du3z+du3(ig)*dbas(1,kdfl(i),4)
          enddo
          
          dav=0d0 ! value
          dax=0d0 ! x dreiv.
          day=0d0 ! y deriv
          daz=0d0 ! y deriv
          do i=1,idfl
            ig=kdfg(i)
            dav=dav+da(ig)*dbas(1,kdfl(i),1)
            dax=dax+da(ig)*dbas(1,kdfl(i),2)
            day=day+da(ig)*dbas(1,kdfl(i),3)
            daz=daz+da(ig)*dbas(1,kdfl(i),4)
          enddo


c       *** evaluate the basis functions in the cubature point
c       for the pressure
          call e011(xi1,xi2,xi3,-3)
          if (ier.lt.0) goto 99999
          
c       evaluate the solution values and derivatives in the cubature point     
          dpv=0d0 ! value
          do i=1,idfl1
            ig=kdfg1(i)
            dpv=dpv+dp(ig)*dbas(1,kdfl1(i),1)
          enddo
          


c       form the integrand

          dnx=-dax
          dny=-day
          dnz=-daz
c
          ah1=-dpv*dnx+dpf1*(du1x*dnx+du1y*dny+du1z*dnz)
          ah2=-dpv*dny+dpf1*(du2x*dnx+du2y*dny+du2z*dnz)
          ah3=-dpv*dnz+dpf1*(du3x*dnx+du3y*dny+du3z*dnz)

          dfwx=dfwx+ah1*om
          dfwy=dfwy+ah2*om
          dfwz=dfwz+ah3*om
          
        enddo
      enddo
c       
c       
99999 end

