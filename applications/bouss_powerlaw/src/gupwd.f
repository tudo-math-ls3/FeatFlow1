************************************************************************
      SUBROUTINE GUPWD (U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,
     *                  A,NA,KCOLA,KLDA,KVERT,KMID,DCORVG,IDEF,
     *                  ICOMP,locny)
************************************************************************
*    Purpose: -  Adds the upwind-part on matrix block A after
*                it was initialized by the Stokes matrix
*             -  The input vector (U1 ,U2 ) is the old velocity field
*             -  The input vectors (UiL1,UiL2) in linear combination  
*                Ail are the transport direction
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL A
C
      DOUBLE PRECISION locny
      dimension locny(*)
      PARAMETER (NNVE=4)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION A(*),KCOLA(*),KLDA(*)
C
C *** Local arrays for informations about one element
      DIMENSION IMID(4),ISTORE(4,4),FLUX(4),UU1(4),UU2(4),XV(4),YV(4)
      DIMENSION ELMA(4,4)
C
C *** Usual data for mesh management
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
C
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
C
      SAVE 
C
C
*********************************************************************
*   weighted upwind nach Samarski
*********************************************************************
c
      PHIP(X)=(0.5D0+X)/(1D0+X)
      PHIM(X)=    0.5D0/(1D0-X)
*********************************************************************
C
C
C *** Loop over all finite elements IEL=1,...,NEL
C
c      write (*,*) 'upwind'
c	return
C
      DO 1  IEL=1,NEL
C
C *** XE,YE will be coordinates of the center of the element
      XE=0.D0
      YE=0.D0
C-----------------------------------------------------------------------
C *** 1. Loop over all 4 U-nodes: IMID,XV,YV,XE,YE,UU1,UU2
C
      DO 11  II=1,4
      I=KMID(II,IEL)-NVT
      IMID(II)=I
      IV=KVERT(II,IEL)
      XV(II)=DCORVG(1,IV)
      YV(II)=DCORVG(2,IV)
      XE=XE+XV(II)
      YE=YE+YV(II)
      UU1(II)=A1L*U1L1(I)+A2L*U2L1(I)
      UU2(II)=A1L*U1L2(I)+A2L*U2L2(I)
  11  CONTINUE
C
      XE=0.25D0*XE
      YE=0.25D0*YE
C-----------------------------------------------------------------------
C *** 2. Loop over all 4 U-nodes:  FLUX(.), ISTORE(.,.), ELMA(.,.)
C
      DO 12  II=1,4
C
C *** Setting II-1 modulo 4 on IM1
      IM1=II-1
      IF (IM1.LT.1) IM1=4
C *** Calculation of the flux  FLUX(II)
      XN=-YV(II)+YE
      YN= XV(II)-XE
      G1=0.5D0*(UU1(IM1)+UU1(II))
      G2=0.5D0*(UU2(IM1)+UU2(II))
      FLUX(II)=XN*G1+YN*G2
C *** Determine the indices ISTORE(II,JJ) to store the element matrix
C     entry ELMA(II,JJ) on array A
      I=IMID(II)
      IA1=KLDA(I)
      IA2=KLDA(I+1)-1
      DO 120  JJ=1,4
      J=IMID(JJ)
C *** Searching loop
      DO 1201  IA=IA1,IA2
      IF (KCOLA(IA).EQ.J)  GOTO 121
 1201 CONTINUE
C
C *** Error case
      WRITE(MTERM,*) 'ERROR in GUPWD: entry index IA not found'
      RETURN
C
  121 ISTORE(II,JJ)=IA
C
C *** Initialization of element matrix ELMA(.,.)
      ELMA(II,JJ)=0D0
C
  120 CONTINUE
  12  CONTINUE
C-----------------------------------------------------------------------
C *** 3. Loop over all 4 U-nodes:  Calculation of ELMA(.,.)
C
      DO 13  II=1,4
C *** Setting II-1 modulo 4 on IM1
      IM0=II
      IM1=II-1
      IF (IM1.LT.1) IM1=4
      IM2=II+1
      IF (IM2.GT.4) IM2=1
C *** Calculate the part corresponding to GAMMA(IM0) and GAMMA(IM2)
      FL0=FLUX(IM0)
      FL2=FLUX(IM2)
C
C
      IF (UPSAM.GE.0) THEN
C
c      UPSRE=UPSAM*RE
c         print *,'upw',locny(iel),1d0/locny(iel)
         UPSRE=UPSAM*RE/locny(iel)
C
       IF (FL0.GE.0D0) THEN
        DL0=PHIP(UPSRE*FL0)
       ELSE
        DL0=PHIM(UPSRE*FL0)
       ENDIF
C
       IF (FL2.GE.0D0) THEN
        DL2=PHIP(UPSRE*FL2)
       ELSE
        DL2=PHIM(UPSRE*FL2)
       ENDIF
C
      ELSE
C
       DL0=0D0
       DL2=0D0
       IF (FL0.GE.0.D0) DL0=1D0
       IF (FL2.GE.0.D0) DL2=1D0
C
      ENDIF
C
C
      H00=DL0*FL0
      H22=(1D0-DL2)*FL2
      ELMA(IM0,IM0) = H00-H22
      ELMA(IM0,IM2) =     H22
      ELMA(IM0,IM1) =-H00
  13  CONTINUE
C-----------------------------------------------------------------------
C *** 4. Loop over all 4 U-nodes:  Addding ELMA(.,.) to matrix A
C
      DO 14   II=1,4
      DO 140  JJ=1,4
      ELMH=THSTEP*ELMA(II,JJ)
C
      IF (IDEF.LT.2) THEN
       IA   =ISTORE(II,JJ)
       A(IA)=A(IA)+REAL(ELMH)
       A(3*NA+IA)=A(3*NA+IA)+REAL(ELMH)
      ENDIF
C
      IF (IDEF.GT.0) THEN 
          D1(IMID(II))=D1(IMID(II))-ELMH*U1(IMID(JJ))
          D2(IMID(II))=D2(IMID(II))-ELMH*U2(IMID(JJ))
      ENDIF 
C
 140  CONTINUE
  14  CONTINUE
C-----------------------------------------------------------------------
C
  1   CONTINUE
C
      END
      
************************************************************************
      SUBROUTINE TGUPWD (U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,
     *                  A,NA,KCOLA,KLDA,KVERT,KMID,DCORVG,IDEF,ICOMP)
************************************************************************
*    Purpose: -  Adds the upwind-part on matrix block A after
*                it was initialized by the Stokes matrix
*             -  The input vector (U1 ,U2 ) is the old velocity field
*             -  The input vectors (UiL1,UiL2) in linear combination  
*                Ail are the transport direction
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL A
C
      PARAMETER (NNVE=4)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION A(*),KCOLA(*),KLDA(*)
C
C *** Local arrays for informations about one element
      DIMENSION IMID(4),ISTORE(4,4),FLUX(4),UU1(4),UU2(4),XV(4),YV(4)
      DIMENSION ELMA(4,4)
C
C *** Usual data for mesh management
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
C
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
C
      SAVE 
C
C
*********************************************************************
*   weighted upwind nach Samarski
*********************************************************************
c
      PHIP(X)=(0.5D0+X)/(1D0+X)
      PHIM(X)=    0.5D0/(1D0-X)
*********************************************************************
C
C
C *** Loop over all finite elements IEL=1,...,NEL
C
c      write (*,*) 'upwind'
c	return
C
      DO 1  IEL=1,NEL
C
C *** XE,YE will be coordinates of the center of the element
      XE=0.D0
      YE=0.D0
C-----------------------------------------------------------------------
C *** 1. Loop over all 4 U-nodes: IMID,XV,YV,XE,YE,UU1,UU2
C
      DO 11  II=1,4
      I=KMID(II,IEL)-NVT
      IMID(II)=I
      IV=KVERT(II,IEL)
      XV(II)=DCORVG(1,IV)
      YV(II)=DCORVG(2,IV)
      XE=XE+XV(II)
      YE=YE+YV(II)
      UU1(II)=A1L*U1L1(I)+A2L*U2L1(I)
      UU2(II)=A1L*U1L2(I)+A2L*U2L2(I)
  11  CONTINUE
C
      XE=0.25D0*XE
      YE=0.25D0*YE
C-----------------------------------------------------------------------
C *** 2. Loop over all 4 U-nodes:  FLUX(.), ISTORE(.,.), ELMA(.,.)
C
      DO 12  II=1,4
C
C *** Setting II-1 modulo 4 on IM1
      IM1=II-1
      IF (IM1.LT.1) IM1=4
C *** Calculation of the flux  FLUX(II)
      XN=-YV(II)+YE
      YN= XV(II)-XE
      G1=0.5D0*(UU1(IM1)+UU1(II))
      G2=0.5D0*(UU2(IM1)+UU2(II))
      FLUX(II)=XN*G1+YN*G2
C *** Determine the indices ISTORE(II,JJ) to store the element matrix
C     entry ELMA(II,JJ) on array A
      I=IMID(II)
      IA1=KLDA(I)
      IA2=KLDA(I+1)-1
      DO 120  JJ=1,4
      J=IMID(JJ)
C *** Searching loop
      DO 1201  IA=IA1,IA2
      IF (KCOLA(IA).EQ.J)  GOTO 121
 1201 CONTINUE
C
C *** Error case
      WRITE(MTERM,*) 'ERROR in GUPWD: entry index IA not found'
      RETURN
C
  121 ISTORE(II,JJ)=IA
C
C *** Initialization of element matrix ELMA(.,.)
      ELMA(II,JJ)=0D0
C
  120 CONTINUE
  12  CONTINUE
C-----------------------------------------------------------------------
C *** 3. Loop over all 4 U-nodes:  Calculation of ELMA(.,.)
C
      DO 13  II=1,4
C *** Setting II-1 modulo 4 on IM1
      IM0=II
      IM1=II-1
      IF (IM1.LT.1) IM1=4
      IM2=II+1
      IF (IM2.GT.4) IM2=1
C *** Calculate the part corresponding to GAMMA(IM0) and GAMMA(IM2)
      FL0=FLUX(IM0)
      FL2=FLUX(IM2)
C
C
      IF (UPSAM.GE.0) THEN
C
       UPSRE=UPSAM*RE
C
       IF (FL0.GE.0D0) THEN
        DL0=PHIP(UPSRE*FL0)
       ELSE
        DL0=PHIM(UPSRE*FL0)
       ENDIF
C
       IF (FL2.GE.0D0) THEN
        DL2=PHIP(UPSRE*FL2)
       ELSE
        DL2=PHIM(UPSRE*FL2)
       ENDIF
C
      ELSE
C
       DL0=0D0
       DL2=0D0
       IF (FL0.GE.0.D0) DL0=1D0
       IF (FL2.GE.0.D0) DL2=1D0
C
      ENDIF
C
C
      H00=DL0*FL0
      H22=(1D0-DL2)*FL2
      ELMA(IM0,IM0) = H00-H22
      ELMA(IM0,IM2) =     H22
      ELMA(IM0,IM1) =-H00
  13  CONTINUE
C-----------------------------------------------------------------------
C *** 4. Loop over all 4 U-nodes:  Addding ELMA(.,.) to matrix A
C
      DO 14   II=1,4
      DO 140  JJ=1,4
      ELMH=THSTEP*ELMA(II,JJ)
C
      IF (IDEF.LT.2) THEN
       IA   =ISTORE(II,JJ)
       A(IA)=A(IA)+REAL(ELMH)
       A(3*NA+IA)=A(3*NA+IA)+REAL(ELMH)
      ENDIF
C
      IF (IDEF.GT.0) THEN 
          D1(IMID(II))=D1(IMID(II))-ELMH*U1(IMID(JJ))
      ENDIF 
C
 140  CONTINUE
  14  CONTINUE
C-----------------------------------------------------------------------
C
  1   CONTINUE
C
      END
      
