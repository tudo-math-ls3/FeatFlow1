************************************************************************
      SUBROUTINE BDRNEU(KABD,KNPR,DCORVG,INEUM,KAREA,KADJ,KVERT,KELBD)
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
       PX=(XV2+XV3+XV6+XV7)*0.25D0
       PY=(YV2+YV3+YV6+YV7)*0.25D0
       PZ=(ZV2+ZV3+ZV6+ZV7)*0.25D0
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
*    Purpose:  updates the right hand with pressure drop data
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
        PX=(XV2+XV3+XV6+XV7)*0.25D0
        PY=(YV2+YV3+YV6+YV7)*0.25D0
        PZ=(ZV2+ZV3+ZV6+ZV7)*0.25D0
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
          IVT3=KVERT(6,-IELBD)
          IVT4=KVERT(7,-IELBD)
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
        PX=(XV2+XV3+XV6+XV7)*0.25D0
        PY=(YV2+YV3+YV6+YV7)*0.25D0
        PZ=(ZV2+ZV3+ZV6+ZV7)*0.25D0
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
          IVT3=KVERT(6,IELBD)
          IVT4=KVERT(7,IELBD)
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
        PX=(XV2+XV3+XV6+XV7)*0.25D0
        PY=(YV2+YV3+YV6+YV7)*0.25D0
        PZ=(ZV2+ZV3+ZV6+ZV7)*0.25D0
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
          IVT3=KVERT(6,IELBD)
          IVT4=KVERT(7,IELBD)
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
         DNZ=0D0
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
