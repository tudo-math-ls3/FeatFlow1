************************************************************************
      SUBROUTINE GUPWD (U1L1,U1L2,U1L3,U2L1,U2L2,U2L3,A1L,A2L,U1,U2,U3,
     *                  D1,D2,D3,A,KCOLA,KLDA,KVERT,KAREA,DCORVG,IDEF)
************************************************************************
*    Purpose: -  Adds the upwind-part on matrix block A after
*                it was initialized by the Stokes matrix
*             -  The input vector (U1,U2,U3) is the old velocity field
*             -  The input vector (UT1,UT2,UT3) is the transport direct.
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      REAL A
C
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION U1L1(*),U1L2(*),U1L3(*),U2L1(*),U2L2(*),U2L3(*)
      DIMENSION U1(*),U2(*),U3(*),D1(*),D2(*),D3(*)
      DIMENSION A(*),KCOLA(*),KLDA(*)
C
C *** Local arrays for informations about one element
      DIMENSION IAREA(6),ISTORE(6,6),FLUX(4),ELMA(6,6),IM(4),IV(4)
      DIMENSION UU1(6),UU2(6),UU3(6),DLAM(4),DLEN(4)
      DIMENSION XV(8),YV(8),ZV(8),XA(6),YA(6),ZA(6)
C
C *** Usual data for mesh management
      DIMENSION KVERT(NNVE,*),KAREA(NNAE,*),DCORVG(3,*)
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
C
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
C
      SAVE 
C
*********************************************************************
*     weighted upwind nach Samarski
      PHIP(X)=(0.5D0+X)/(1D0+X)
      PHIM(X)= 0.5D0   /(1D0-X)
*********************************************************************
C
C *** Loop over all finite elements IEL=1,...,NEL
C
      DO 1  IEL=1,NEL
C
C *** XE,YE,ZE will be coordinates of the center of the element
      XE=0D0
      YE=0D0
      ZE=0D0
C-----------------------------------------------------------------------
C *** 1. Loop over all 8 vertices
C
      DO 10  II=1,8
      IVT   =KVERT(II,IEL)
      XV(II)=DCORVG(1,IVT)
      YV(II)=DCORVG(2,IVT)
      ZV(II)=DCORVG(3,IVT)
      XE=XE+XV(II)
      YE=YE+YV(II)
      ZE=ZE+ZV(II)
10    CONTINUE
C
      XE=0.125D0*XE
      YE=0.125D0*YE
      ZE=0.125D0*ZE
C
C-----------------------------------------------------------------------
C *** 2. Loop over all 6 U-nodes
C
      DO 11  II=1,6
      IAR      =KAREA(II,IEL)
      IAREA(II)=IAR
C
      IF (II.EQ.1) THEN
       IV1=1
       IV2=2
       IV3=3
       IV4=4
      ENDIF
C
      IF (II.EQ.2) THEN
       IV1=1
       IV2=2
       IV3=6
       IV4=5
      ENDIF
C
      IF (II.EQ.3) THEN
       IV1=2
       IV2=3
       IV3=7
       IV4=6
      ENDIF
C
      IF (II.EQ.4) THEN
       IV1=3
       IV2=4
       IV3=8
       IV4=7
      ENDIF
C
      IF (II.EQ.5) THEN
       IV1=4
       IV2=1
       IV3=5
       IV4=8
      ENDIF
C
      IF (II.EQ.6) THEN
       IV1=6
       IV2=5
       IV3=8
       IV4=7
      ENDIF
C
      XA(II)=(XV(IV1)+XV(IV2)+XV(IV3)+XV(IV4))*0.25D0
      YA(II)=(YV(IV1)+YV(IV2)+YV(IV3)+YV(IV4))*0.25D0
      ZA(II)=(ZV(IV1)+ZV(IV2)+ZV(IV3)+ZV(IV4))*0.25D0
C
      UU1(II)=A1L*U1L1(IAR)+A2L*U2L1(IAR)
      UU2(II)=A1L*U1L2(IAR)+A2L*U2L2(IAR)
      UU3(II)=A1L*U1L3(IAR)+A2L*U2L3(IAR)
C
11    CONTINUE
C
C-----------------------------------------------------------------------
C *** 3. Loop over all 6 U-nodes
C
      DO 12  II=1,6
C
C *** Determine the indices ISTORE(II,JJ) to store the element matrix
C     entry ELMA(II,JJ) on array A
      IAR=IAREA(II)
      IA1=KLDA(IAR)
      IA2=KLDA(IAR+1)-1
      DO 120  JJ=1,6
      JAR=IAREA(JJ)
C *** Searching loop
      DO 1201  IA=IA1,IA2
      IF (KCOLA(IA).EQ.JAR)  GOTO 1202
1201  CONTINUE
C *** Error case
      WRITE(MTERM,*) 'ERROR in GUPWD: entry index IA not found'
      RETURN
1202  ISTORE(II,JJ)=IA
C *** Initialization of element matrix ELMA(.,.)
      ELMA(II,JJ)=0D0
120   CONTINUE
C
C
      IM0=II
C
      IF (II.EQ.1) THEN
       IM(1)=2
       IM(2)=3
       IM(3)=4
       IM(4)=5
C
       IV(1)=1
       IV(2)=2
       IV(3)=3
       IV(4)=4
      ENDIF
C
      IF (II.EQ.2) THEN
       IM(1)=1
       IM(2)=3
       IM(3)=6
       IM(4)=5
C
       IV(1)=1
       IV(2)=2
       IV(3)=6
       IV(4)=5
      ENDIF
C
      IF (II.EQ.3) THEN
       IM(1)=1
       IM(2)=4
       IM(3)=6
       IM(4)=2
C
       IV(1)=2
       IV(2)=3
       IV(3)=7
       IV(4)=6
      ENDIF
C
      IF (II.EQ.4) THEN
       IM(1)=1
       IM(2)=5
       IM(3)=6
       IM(4)=3
C
       IV(1)=3
       IV(2)=4
       IV(3)=8
       IV(4)=7
      ENDIF
C
      IF (II.EQ.5) THEN
       IM(1)=1
       IM(2)=2
       IM(3)=6
       IM(4)=4
C
       IV(1)=4
       IV(2)=1
       IV(3)=5
       IV(4)=8
      ENDIF
C
      IF (II.EQ.6) THEN
       IM(1)=2
       IM(2)=5
       IM(3)=4
       IM(4)=3
C
       IV(1)=6
       IV(2)=5
       IV(3)=8
       IV(4)=7
      ENDIF
C
      DO 121 I=1,4
       I1=I
       I2=I+1
       IF (I2.EQ.5) I2=1
       XNH1=XV(IV(I1))-XE
       YNH1=YV(IV(I1))-YE
       ZNH1=ZV(IV(I1))-ZE
       XNH2=XV(IV(I2))-XE
       YNH2=YV(IV(I2))-YE
       ZNH2=ZV(IV(I2))-ZE
       XNH =YNH1*ZNH2-ZNH1*YNH2
       YNH =ZNH1*XNH2-XNH1*ZNH2
       ZNH =XNH1*YNH2-YNH1*XNH2
       DLEN(I)=SQRT(XNH**2+YNH**2+ZNH**2)
       G1=0.5D0*(UU1(IM(I))+UU1(IM0))
       G2=0.5D0*(UU2(IM(I))+UU2(IM0))
       G3=0.5D0*(UU3(IM(I))+UU3(IM0))
       XNHH=XA(IM(I))-XA(IM0)
       YNHH=YA(IM(I))-YA(IM0)
       ZNHH=ZA(IM(I))-ZA(IM0)
       DINN=XNH*XNHH+YNH*YNHH+ZNH*ZNHH
       IF (DINN.LT.0D0) THEN
        XNH=-XNH
        YNH=-YNH
        ZNH=-ZNH
       ENDIF
       FLUX(I)=0.5D0*(XNH*G1+YNH*G2+ZNH*G3)
121   CONTINUE
C
C
C
      IF (UPSAM.GE.0) THEN
C
       DO 128  I=1,4
        IF (FLUX(I).GE.0D0) THEN
         DLAM(I)=PHIP(UPSAM*RE*FLUX(I)/SQRT(DLEN(I)))
        ELSE
         DLAM(I)=PHIM(UPSAM*RE*FLUX(I)/SQRT(DLEN(I)))
        ENDIF
128    CONTINUE
C
      ELSE
C
       DO 129  I=1,4
        DLAM(I)=0D0
        IF (FLUX(I).GE.0.D0) DLAM(I)=1D0
129    CONTINUE
C
      ENDIF
C
C
      H1=FLUX(1)*(1D0-DLAM(1))
      H2=FLUX(2)*(1D0-DLAM(2))
      H3=FLUX(3)*(1D0-DLAM(3))
      H4=FLUX(4)*(1D0-DLAM(4))
      ELMA(IM0,IM0)  =-H1-H2-H3-H4
      ELMA(IM0,IM(1))= H1
      ELMA(IM0,IM(2))= H2
      ELMA(IM0,IM(3))= H3
      ELMA(IM0,IM(4))= H4
C
12    CONTINUE
C-----------------------------------------------------------------------
C *** 4. Loop over all 4 U-nodes:  Addding ELMA(.,.) to matrix A
C
      DO 13   II=1,6
      DO 130  JJ=1,6
      ELMH=THSTEP*ELMA(II,JJ)
      IF (ABS(ELMH).LT.1D-15) GOTO 130
C
      IF (IDEF.LT.2) THEN
       IA=ISTORE(II,JJ)
       A(IA)=A(IA)+REAL(ELMH)
      ENDIF
C
      IF (IDEF.GT.0) THEN
       IAREAI=IAREA(II)
       IAREAJ=IAREA(JJ)
       D1(IAREAI)=D1(IAREAI)-ELMH*U1(IAREAJ)
       D2(IAREAI)=D2(IAREAI)-ELMH*U2(IAREAJ)
       D3(IAREAI)=D3(IAREAI)-ELMH*U3(IAREAJ)
      ENDIF 
C
130   CONTINUE
13    CONTINUE
C-----------------------------------------------------------------------
C
1     CONTINUE
C
      END
      
