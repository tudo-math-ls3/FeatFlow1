************************************************************************
      SUBROUTINE GUPWD (U1L1,U1L2,U1L3,U2L1,U2L2,U2L3,A1L,A2L,U1,U2,U3,
     *                  D1,D2,D3,A,KCOLA,KLDA,KVERT,KAREA,DCORVG,IDEF)
************************************************************************
*    Purpose: -  Adds the upwind-part on matrix block A after
*                it was initialized by the Stokes matrix
*             -  The input vector (U1,U2,U3) is the old velocity field
*             -  The input vectors (UiLj,..) are the transport direct.
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
      DIMENSION IAREA(6),ISTORE(6,6),FLUX(4),ELMA(6,6),IM(4)
      DIMENSION UU1(6),UU2(6),UU3(6),DLAM(4),DLEN(4,6)
      DIMENSION XV(8),YV(8),ZV(8),XVE(8),YVE(8),ZVE(8),
     *          XA(6),YA(6),ZA(6),XN(4,6),YN(4,6),ZN(4,6),
     *          UMEAN1(4,6),UMEAN2(4,6),UMEAN3(4,6)
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
      DO 100  II=1,8
      IVT   =KVERT(II,IEL)
      XV(II)=DCORVG(1,IVT)
      YV(II)=DCORVG(2,IVT)
      ZV(II)=DCORVG(3,IVT)
      XE=XE+XV(II)
      YE=YE+YV(II)
      ZE=ZE+ZV(II)
100   CONTINUE
C
      XE=0.125D0*XE
      YE=0.125D0*YE
      ZE=0.125D0*ZE
C
      DO 110  II=1,8
      XVE(II)=XV(II)-XE
      YVE(II)=YV(II)-YE
      ZVE(II)=ZV(II)-ZE
110   CONTINUE
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
       XA(II)=(XV(IV1)+XV(IV2)+XV(IV3)+XV(IV4))*0.25D0
       YA(II)=(YV(IV1)+YV(IV2)+YV(IV3)+YV(IV4))*0.25D0
       ZA(II)=(ZV(IV1)+ZV(IV2)+ZV(IV3)+ZV(IV4))*0.25D0
      ENDIF
C
      IF (II.EQ.2) THEN
       IV1=1
       IV2=2
       IV3=6
       IV4=5
       XA(II)=(XV(IV1)+XV(IV2)+XV(IV3)+XV(IV4))*0.25D0
       YA(II)=(YV(IV1)+YV(IV2)+YV(IV3)+YV(IV4))*0.25D0
       ZA(II)=(ZV(IV1)+ZV(IV2)+ZV(IV3)+ZV(IV4))*0.25D0
      ENDIF
C
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
      IAR=IAREA(II)
      IA1=KLDA(IAR)
      IA2=KLDA(IAR+1)-1
C
      DO 120  JJ=1,6
      JAR=IAREA(JJ)
C
      DO 1201  IA=IA1,IA2
      IF (KCOLA(IA).EQ.JAR)  GOTO 1202
1201  CONTINUE
      RETURN
C
1202  ISTORE(II,JJ)=IA
      ELMA(II,JJ)=0D0
120   CONTINUE
12    CONTINUE
C
C
C
      XN  (1,1)=YVE(1)*ZVE(2)-ZVE(1)*YVE(2)
      YN  (1,1)=ZVE(1)*XVE(2)-XVE(1)*ZVE(2)
      ZN  (1,1)=XVE(1)*YVE(2)-YVE(1)*XVE(2)
C
      XNHH=XA(2)-XA(1)
      YNHH=YA(2)-YA(1)
      ZNHH=ZA(2)-ZA(1)
      DINN=XN(1,1)*XNHH+YN(1,1)*YNHH+ZN(1,1)*ZNHH
      IF (DINN.LT.0D0) THEN
       DNFACT=-1D0
      ELSE
       DNFACT= 1D0
      ENDIF
C
      IF (DNFACT.LT.0D0) THEN
       XN  (1,1)=-XN(1,1)
       YN  (1,1)=-YN(1,1)
       ZN  (1,1)=-ZN(1,1)
      ENDIF
C
C
      XN  (2,1)= DNFACT*(YVE(2)*ZVE(3)-ZVE(2)*YVE(3))
      YN  (2,1)= DNFACT*(ZVE(2)*XVE(3)-XVE(2)*ZVE(3))
      ZN  (2,1)= DNFACT*(XVE(2)*YVE(3)-YVE(2)*XVE(3))
C
      XN  (3,1)= DNFACT*(YVE(3)*ZVE(4)-ZVE(3)*YVE(4))
      YN  (3,1)= DNFACT*(ZVE(3)*XVE(4)-XVE(3)*ZVE(4))
      ZN  (3,1)= DNFACT*(XVE(3)*YVE(4)-YVE(3)*XVE(4))
C
      XN  (4,1)= DNFACT*(YVE(4)*ZVE(1)-ZVE(4)*YVE(1))
      YN  (4,1)= DNFACT*(ZVE(4)*XVE(1)-XVE(4)*ZVE(1))
      ZN  (4,1)= DNFACT*(XVE(4)*YVE(1)-YVE(4)*XVE(1))
C
C
C
      XN  (1,2)=-XN(1,1)
      YN  (1,2)=-YN(1,1)
      ZN  (1,2)=-ZN(1,1)
C
      XN  (2,2)=-DNFACT*(YVE(2)*ZVE(6)-ZVE(2)*YVE(6))
      YN  (2,2)=-DNFACT*(ZVE(2)*XVE(6)-XVE(2)*ZVE(6))
      ZN  (2,2)=-DNFACT*(XVE(2)*YVE(6)-YVE(2)*XVE(6))
C
      XN  (3,2)=-DNFACT*(YVE(6)*ZVE(5)-ZVE(6)*YVE(5))
      YN  (3,2)=-DNFACT*(ZVE(6)*XVE(5)-XVE(6)*ZVE(5))
      ZN  (3,2)=-DNFACT*(XVE(6)*YVE(5)-YVE(6)*XVE(5))
C
      XN  (4,2)=-DNFACT*(YVE(5)*ZVE(1)-ZVE(5)*YVE(1))
      YN  (4,2)=-DNFACT*(ZVE(5)*XVE(1)-XVE(5)*ZVE(1))
      ZN  (4,2)=-DNFACT*(XVE(5)*YVE(1)-YVE(5)*XVE(1))
C
C
      XN  (1,3)=-XN(2,1)
      YN  (1,3)=-YN(2,1)
      ZN  (1,3)=-ZN(2,1)
C
      XN  (2,3)=-DNFACT*(YVE(3)*ZVE(7)-ZVE(3)*YVE(7))
      YN  (2,3)=-DNFACT*(ZVE(3)*XVE(7)-XVE(3)*ZVE(7))
      ZN  (2,3)=-DNFACT*(XVE(3)*YVE(7)-YVE(3)*XVE(7))
C
      XN  (3,3)=-DNFACT*(YVE(7)*ZVE(6)-ZVE(7)*YVE(6))
      YN  (3,3)=-DNFACT*(ZVE(7)*XVE(6)-XVE(7)*ZVE(6))
      ZN  (3,3)=-DNFACT*(XVE(7)*YVE(6)-YVE(7)*XVE(6))
C
      XN  (4,3)=-XN(2,2)
      YN  (4,3)=-YN(2,2)
      ZN  (4,3)=-ZN(2,2)
C
C
      XN  (1,4)=-XN(3,1)
      YN  (1,4)=-YN(3,1)
      ZN  (1,4)=-ZN(3,1)
C
      XN  (2,4)=-DNFACT*(YVE(4)*ZVE(8)-ZVE(4)*YVE(8))
      YN  (2,4)=-DNFACT*(ZVE(4)*XVE(8)-XVE(4)*ZVE(8))
      ZN  (2,4)=-DNFACT*(XVE(4)*YVE(8)-YVE(4)*XVE(8))
C
      XN  (3,4)=-DNFACT*(YVE(8)*ZVE(7)-ZVE(8)*YVE(7))
      YN  (3,4)=-DNFACT*(ZVE(8)*XVE(7)-XVE(8)*ZVE(7))
      ZN  (3,4)=-DNFACT*(XVE(8)*YVE(7)-YVE(8)*XVE(7))
C
      XN  (4,4)=-XN(2,3)
      YN  (4,4)=-YN(2,3)
      ZN  (4,4)=-ZN(2,3)
C
C
      XN  (1,5)=-XN(4,1)
      YN  (1,5)=-YN(4,1)
      ZN  (1,5)=-ZN(4,1)
C
      XN  (2,5)=-XN(4,2)
      YN  (2,5)=-YN(4,2)
      ZN  (2,5)=-ZN(4,2)
C
      XN  (3,5)=-DNFACT*(YVE(5)*ZVE(8)-ZVE(5)*YVE(8))
      YN  (3,5)=-DNFACT*(ZVE(5)*XVE(8)-XVE(5)*ZVE(8))
      ZN  (3,5)=-DNFACT*(XVE(5)*YVE(8)-YVE(5)*XVE(8))
C
      XN  (4,5)=-XN(2,4)
      YN  (4,5)=-YN(2,4)
      ZN  (4,5)=-ZN(2,4)
C
C
      XN  (1,6)=-XN(3,2)
      YN  (1,6)=-YN(3,2)
      ZN  (1,6)=-ZN(3,2)
C
      XN  (2,6)=-XN(3,5)
      YN  (2,6)=-YN(3,5)
      ZN  (2,6)=-ZN(3,5)
C
      XN  (3,6)=-XN(3,4)
      YN  (3,6)=-YN(3,4)
      ZN  (3,6)=-ZN(3,4)
C
      XN  (4,6)=-XN(3,3)
      YN  (4,6)=-YN(3,3)
      ZN  (4,6)=-ZN(3,3)
C
C
      DLEN(1,1)=SQRT(XN(1,1)**2+YN(1,1)**2+ZN(1,1)**2)
      DLEN(2,1)=SQRT(XN(2,1)**2+YN(2,1)**2+ZN(2,1)**2)
      DLEN(3,1)=SQRT(XN(3,1)**2+YN(3,1)**2+ZN(3,1)**2)
      DLEN(4,1)=SQRT(XN(4,1)**2+YN(4,1)**2+ZN(4,1)**2)
C
      DLEN(1,2)=DLEN(1,1)
      DLEN(2,2)=SQRT(XN(2,2)**2+YN(2,2)**2+ZN(2,2)**2)
      DLEN(3,2)=SQRT(XN(3,2)**2+YN(3,2)**2+ZN(3,2)**2)
      DLEN(4,2)=SQRT(XN(4,2)**2+YN(4,2)**2+ZN(4,2)**2)
C
      DLEN(1,3)=DLEN(2,1)
      DLEN(2,3)=SQRT(XN(2,3)**2+YN(2,3)**2+ZN(2,3)**2)
      DLEN(3,3)=SQRT(XN(3,3)**2+YN(3,3)**2+ZN(3,3)**2)
      DLEN(4,3)=DLEN(2,2)
C
      DLEN(1,4)=DLEN(3,1)
      DLEN(2,4)=SQRT(XN(2,4)**2+YN(2,4)**2+ZN(2,4)**2)
      DLEN(3,4)=SQRT(XN(3,4)**2+YN(3,4)**2+ZN(3,4)**2)
      DLEN(4,4)=DLEN(2,3)
C
      DLEN(1,5)=DLEN(4,1)
      DLEN(2,5)=DLEN(4,2)
      DLEN(3,5)=SQRT(XN(3,5)**2+YN(3,5)**2+ZN(3,5)**2)
      DLEN(4,5)=DLEN(4,4)
C
      DLEN(1,6)=DLEN(3,2)
      DLEN(2,6)=DLEN(3,5)
      DLEN(3,6)=DLEN(3,4)
      DLEN(4,6)=DLEN(3,3)
C
C
      DMEAN=0.5D0
C
      UMEAN1(1,1)=DMEAN*(UU1(2)+UU1(1))
      UMEAN2(1,1)=DMEAN*(UU2(2)+UU2(1))
      UMEAN3(1,1)=DMEAN*(UU3(2)+UU3(1))
C
      UMEAN1(2,1)=DMEAN*(UU1(3)+UU1(1))
      UMEAN2(2,1)=DMEAN*(UU2(3)+UU2(1))
      UMEAN3(2,1)=DMEAN*(UU3(3)+UU3(1))
C
      UMEAN1(3,1)=DMEAN*(UU1(4)+UU1(1))
      UMEAN2(3,1)=DMEAN*(UU2(4)+UU2(1))
      UMEAN3(3,1)=DMEAN*(UU3(4)+UU3(1))
C
      UMEAN1(4,1)=DMEAN*(UU1(5)+UU1(1))
      UMEAN2(4,1)=DMEAN*(UU2(5)+UU2(1))
      UMEAN3(4,1)=DMEAN*(UU3(5)+UU3(1))
C
C
      UMEAN1(1,2)=UMEAN1(1,1)
      UMEAN2(1,2)=UMEAN2(1,1)
      UMEAN3(1,2)=UMEAN3(1,1)
C
      UMEAN1(2,2)=DMEAN*(UU1(3)+UU1(2))
      UMEAN2(2,2)=DMEAN*(UU2(3)+UU2(2))
      UMEAN3(2,2)=DMEAN*(UU3(3)+UU3(2))
C
      UMEAN1(3,2)=DMEAN*(UU1(6)+UU1(2))
      UMEAN2(3,2)=DMEAN*(UU2(6)+UU2(2))
      UMEAN3(3,2)=DMEAN*(UU3(6)+UU3(2))
C
      UMEAN1(4,2)=DMEAN*(UU1(5)+UU1(2))
      UMEAN2(4,2)=DMEAN*(UU2(5)+UU2(2))
      UMEAN3(4,2)=DMEAN*(UU3(5)+UU3(2))
C
C
      UMEAN1(1,3)=UMEAN1(2,1)
      UMEAN2(1,3)=UMEAN2(2,1)
      UMEAN3(1,3)=UMEAN3(2,1)
C
      UMEAN1(2,3)=DMEAN*(UU1(4)+UU1(3))
      UMEAN2(2,3)=DMEAN*(UU2(4)+UU2(3))
      UMEAN3(2,3)=DMEAN*(UU3(4)+UU3(3))
C
      UMEAN1(3,3)=DMEAN*(UU1(6)+UU1(3))
      UMEAN2(3,3)=DMEAN*(UU2(6)+UU2(3))
      UMEAN3(3,3)=DMEAN*(UU3(6)+UU3(3))
C
      UMEAN1(4,3)=UMEAN1(2,2)
      UMEAN2(4,3)=UMEAN2(2,2)
      UMEAN3(4,3)=UMEAN3(2,2)
C
C
      UMEAN1(1,4)=UMEAN1(3,1)
      UMEAN2(1,4)=UMEAN2(3,1)
      UMEAN3(1,4)=UMEAN3(3,1)
C
      UMEAN1(2,4)=DMEAN*(UU1(5)+UU1(4))
      UMEAN2(2,4)=DMEAN*(UU2(5)+UU2(4))
      UMEAN3(2,4)=DMEAN*(UU3(5)+UU3(4))
C
      UMEAN1(3,4)=DMEAN*(UU1(6)+UU1(4))
      UMEAN2(3,4)=DMEAN*(UU2(6)+UU2(4))
      UMEAN3(3,4)=DMEAN*(UU3(6)+UU3(4))
C
      UMEAN1(4,4)=UMEAN1(2,3)
      UMEAN2(4,4)=UMEAN2(2,3)
      UMEAN3(4,4)=UMEAN3(2,3)
C
C
      UMEAN1(1,5)=UMEAN1(4,1)
      UMEAN2(1,5)=UMEAN2(4,1)
      UMEAN3(1,5)=UMEAN3(4,1)
C
      UMEAN1(2,5)=UMEAN1(4,2)
      UMEAN2(2,5)=UMEAN2(4,2)
      UMEAN3(2,5)=UMEAN3(4,2)
C
      UMEAN1(3,5)=DMEAN*(UU1(6)+UU1(5))
      UMEAN2(3,5)=DMEAN*(UU2(6)+UU2(5))
      UMEAN3(3,5)=DMEAN*(UU3(6)+UU3(5))
C
      UMEAN1(4,5)=UMEAN1(2,4)
      UMEAN2(4,5)=UMEAN2(2,4)
      UMEAN3(4,5)=UMEAN3(2,4)
C
C
      UMEAN1(1,6)=UMEAN1(3,2)
      UMEAN2(1,6)=UMEAN2(3,2)
      UMEAN3(1,6)=UMEAN3(3,2)
C
      UMEAN1(2,6)=UMEAN1(3,5)
      UMEAN2(2,6)=UMEAN2(3,5)
      UMEAN3(2,6)=UMEAN3(3,5)
C
      UMEAN1(3,6)=UMEAN1(3,4)
      UMEAN2(3,6)=UMEAN2(3,4)
      UMEAN3(3,6)=UMEAN3(3,4)
C
      UMEAN1(4,6)=UMEAN1(3,3)
      UMEAN2(4,6)=UMEAN2(3,3)
      UMEAN3(4,6)=UMEAN3(3,3)
C
C
C
      DO 13  II=1,6
C
      IM0=II
C
      IF (II.EQ.1) THEN
       IM(1)=2
       IM(2)=3
       IM(3)=4
       IM(4)=5
      ENDIF
C
      IF (II.EQ.2) THEN
       IM(1)=1
       IM(2)=3
       IM(3)=6
       IM(4)=5
      ENDIF
C
      IF (II.EQ.3) THEN
       IM(1)=1
       IM(2)=4
       IM(3)=6
       IM(4)=2
      ENDIF
C
      IF (II.EQ.4) THEN
       IM(1)=1
       IM(2)=5
       IM(3)=6
       IM(4)=3
      ENDIF
C
      IF (II.EQ.5) THEN
       IM(1)=1
       IM(2)=2
       IM(3)=6
       IM(4)=4
      ENDIF
C
      IF (II.EQ.6) THEN
       IM(1)=2
       IM(2)=5
       IM(3)=4
       IM(4)=3
      ENDIF
C
      DO 130 I=1,4
      FLUX(I)=0.5D0*( XN(I,II)*UMEAN1(I,II)+YN(I,II)*UMEAN2(I,II)
     *               +ZN(I,II)*UMEAN3(I,II))
C
      IF (UPSAM.GE.0) THEN
C
       IF (FLUX(I).GE.0D0) THEN
        DLAM(I)=PHIP(UPSAM*RE*FLUX(I)/SQRT(DLEN(I,II)))
       ELSE
        DLAM(I)=PHIM(UPSAM*RE*FLUX(I)/SQRT(DLEN(I,II)))
       ENDIF
C
      ELSE
C
       DLAM(I)=0D0
       IF (FLUX(I).GE.0.D0) DLAM(I)=1D0
C
      ENDIF
C
130   CONTINUE
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
13    CONTINUE
C-----------------------------------------------------------------------
C *** 4. Loop over all 4 U-nodes:  Addding ELMA(.,.) to matrix A
C
      DO 14   II=1,6
      DO 140  JJ=1,6
      ELMH=THSTEP*ELMA(II,JJ)
      IF (ABS(ELMH).LT.1D-15) GOTO 140
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
140   CONTINUE
14    CONTINUE
C-----------------------------------------------------------------------
C
1     CONTINUE
C
      END
      
