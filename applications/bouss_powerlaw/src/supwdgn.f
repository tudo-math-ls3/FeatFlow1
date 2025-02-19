************************************************************************
      DOUBLE PRECISION FUNCTION DVISCO(D)
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
      common /fluidp/ pla,ple
      SAVE 
C     
C     D..square of abs of symm part of grad
C
C     DVISCO=NY*(1d0+(PLE+D)**(-PLA/2d0))
c     DVISCO=NY*(PLE+D)**(-PLA/2d0)
      DVISCO=NY*(PLE+sqrt(D))**(-PLA)
C     
      END
c
c
************************************************************************
      SUBROUTINE SUPWDG(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                  KCOLA,KLDA,KVERT,KMID,DCORVG,ELE,COEF,
     *                  IDEF,ICOMP,DCMASS,locny)
************************************************************************
*     Purpose: -  Adds the SUPG-part on matrix block A after
*                 it was initialized by the linear part
*              -  The input vector Ui is the old velocity field
*              -  The input vectors UjLi are the transport directions
*     PARAMETRIC VERSION
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21)
C
      DOUBLE PRECISION locny
      dimension locny(*)
      DIMENSION A(*),KCOLA(*),KLDA(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRY(NNBAS,NNBAS,4)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
C
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
      common /fluidp/ pla,ple
      SAVE 
C
C
      IF (IPRECA.EQ.4) THEN
         DNY=NY
      ELSE
         DNY=0D0
      ENDIF
C     
      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
         CT0=DCMASS/THSTEP
      ELSE
         CT0=0D0
      ENDIF
C
      ICUB=ICUBN
C
      DO 1 I = 1,NNDER
1     BDER(I)=.FALSE.
C
      DO 2 I=1,3
2     BDER(I)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      ICUBP=ICUB
      CALL ELE(0D0,0D0,-2)
C
      DUMAX=0D0
      IF (A2L.EQ.0D0) THEN
       DO 10 IEQ=1,NMT
       DU1=U1L1(IEQ)
       DU2=U1L2(IEQ)
       DUNORM=SQRT(DU1**2+DU2**2)
       DUMAX=MAX(DUMAX,DUNORM)
10     CONTINUE
      ELSE       
       DO 20 IEQ=1,NMT
       DU1=A1L*U1L1(IEQ)+A2L*U2L1(IEQ)
       DU2=A1L*U1L2(IEQ)+A2L*U2L2(IEQ)
       DUNORM=SQRT(DU1**2+DU2**2)
       DUMAX=MAX(DUMAX,DUNORM)
20     CONTINUE
      ENDIF       
C
      IF (DUMAX.LT.1D-8) DUMAX=1D-8
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine local DELTA for streamline-diffusion
      IF (ISTOK.EQ.0) then
         CALL DELTSD(U1L1,U1L2,U2L1,U2L2,A1L,A2L,IEL,
     *               DUMAX,DELTA,KVERT,KMID,DCORVG)
         else
            DELTA=0D0
         endif
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
         DENTRY(JDOFE,JDOFE,1)=0D0
         DENTRY(JDOFE,JDOFE,2)=0D0
         DENTRY(JDOFE,JDOFE,3)=0D0
         DENTRY(JDOFE,JDOFE,4)=0D0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
          DENTRY(JDOFE,IDOFE,1)=0D0
          DENTRY(JDOFE,IDOFE,2)=0D0
          DENTRY(JDOFE,IDOFE,3)=0D0
          DENTRY(JDOFE,IDOFE,4)=0D0
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices
      DO 120 IVE = 1, NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
120   CONTINUE
C
      DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
      DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
      DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
      DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
C
      avgny=0d0
C *** Loop over all cubature points
      DO 200 ICUBP = 1, NCUBP
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
C *** Jacobian of the bilinear mapping onto the reference element
      DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
      DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
      DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
      DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
      DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
      OM=DOMEGA(ICUBP)*DETJ
C
      CALL ELE(XI1,XI2,-3)
      IF (IER.LT.0) GOTO 99999
C
C *** Evaluation of velocity in cubature points
      IF (DCMASS.NE.0D0) THEN
C
       DU1=0D0
       DU2=0D0
       IF (A2L.EQ.0D0) THEN
        DO 210 JDFL=1,IDFL
        HBAS=DBAS(KDFL(JDFL),1)
        IF (ABS(HBAS).GE.1D-8) THEN
         JDFG=KDFG(JDFL)
         DU1=DU1+U1L1(JDFG)*HBAS
         DU2=DU2+U1L2(JDFG)*HBAS
        ENDIF
210     CONTINUE
       ELSE
        DO 220 JDFL=1,IDFL
        HBAS=DBAS(KDFL(JDFL),1)
        IF (ABS(HBAS).GE.1D-8) THEN
         JDFG=KDFG(JDFL)
         DU1=DU1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
         DU2=DU2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
         ENDIF
220     CONTINUE
       ENDIF
C     and derivatives
          DUX1=0D0
          DUX2=0D0
          IF (A2L.EQ.0D0) THEN
           DO 211 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),2)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUX1=DUX1+U1L1(JDFG)*HBAS
             DUX2=DUX2+U1L2(JDFG)*HBAS
            ENDIF
 211       CONTINUE
          ELSE
           DO 221 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),2)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUX1=DUX1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
             DUX2=DUX2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
            ENDIF
 221       CONTINUE
          ENDIF
          
          DUY1=0D0
          DUY2=0D0
          IF (A2L.EQ.0D0) THEN
           DO 212 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),3)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUY1=DUY1+U1L1(JDFG)*HBAS
             DUY2=DUY2+U1L2(JDFG)*HBAS
            ENDIF
 212       CONTINUE
          ELSE
           DO 222 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),3)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUY1=DUY1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
             DUY2=DUY2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
            ENDIF
 222       CONTINUE
          ENDIF
c
C     
C     norm of symmetric part of velocity gradient
          DUGSQ=DUX1**2+0.5d0*(DUX2+DUY1)**2+DUY2**2

          if (PLA.eq.0d0) then 
           DNY=NY
          else
           DNY=DVISCO(DUGSQ)
          endif
          avgny=avgny+dny/NY
C     
C *** Summing up over all pairs of multiindices
       if (IGRAD.eq.1) then
c     use just grad
C ***  Summing up over all pairs of multiindices
       DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2
C
       DO 235 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH4= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
        HSUMI=HBASI2*DU1+HBASI3*DU2
             AH1= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +   DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH4= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +   DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            
 235        CONTINUE
230    CONTINUE
          else
c     use symmetric part of grad
C ***  Summing up over all pairs of multiindices
       DO 240 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2
C
       DO 245 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(2d0*HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH2= DNY*HBASJ2*HBASJ3
             AH3= DNY*HBASJ3*HBASJ2
             AH4= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(HBASJ2**2+2d0*HBASJ3**2)+CT0*HBASJ1**2
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
        HSUMI=HBASI2*DU1+HBASI3*DU2
             AH1= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +   DNY*(2d0*HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH2= DNY*HBASI3*HBASJ2
             AH3= DNY*HBASI2*HBASJ3
             AH4= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +   DNY*(HBASJ2*HBASI2+2d0*HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,2)=DENTRY(JDOFE,IDOFE,2)+OM*AH2
            DENTRY(JDOFE,IDOFE,3)=DENTRY(JDOFE,IDOFE,3)+OM*AH3
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            
 245        CONTINUE
240    CONTINUE
       endif
C
      ELSE
C
C ***  Summing up over 1 pair of multiindices
       DO 250 JDOFE=1,IDFL
       HBASJ1=DBAS(KDFL(JDOFE),1)
C
       DO 260 IDOFE=1,IDFL
       HBASI1=DBAS(KDFL(IDOFE),1)
C
       AH=-1D0/THSTEP*HBASJ1*HBASI1
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH
260    CONTINUE
250    CONTINUE
C
      ENDIF
C
200   CONTINUE
C     
      locny(iel)=avgny/ncubp
C     
      DO 400 JDOFE=1,IDFL
         DO 410 IDOFE=1,IDFL
            DENTH1=THSTEP*DENTRY(JDOFE,IDOFE,1)
            DENTH2=THSTEP*DENTRY(JDOFE,IDOFE,2)
            DENTH3=THSTEP*DENTRY(JDOFE,IDOFE,3)
            DENTH4=THSTEP*DENTRY(JDOFE,IDOFE,4)
C     
            IF (IDEF.LT.2) THEN
               IA   =KENTRY(JDOFE,IDOFE)
               A(IA)=A(IA)+REAL(DENTH1)
               A(NA+IA)=A(NA+IA)+REAL(DENTH2)
               A(2*NA+IA)=A(2*NA+IA)+REAL(DENTH3)
               A(3*NA+IA)=A(3*NA+IA)+REAL(DENTH4)
            ENDIF
C     
            IF (IDEF.GT.0) THEN 
               IDFG=KDFG(IDOFE)
               JDFG=KDFG(JDOFE)
               D1(JDFG)=D1(JDFG)-DENTH1*U1(IDFG)-DENTH2*U2(IDFG)
               D2(JDFG)=D2(JDFG)-DENTH3*U1(IDFG)-DENTH4*U2(IDFG)
            ENDIF 
C     
 410     CONTINUE
 400  continue
C     
 100  CONTINUE
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE SUPWNP(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                  KCOLA,KLDA,KVERT,KMID,DCORVG,ELE,COEF,
     *                  IDEF,ICOMP,DCMASS,locny)
************************************************************************
*     Purpose: -  Adds the SUPG-part on matrix block A after
*                 it was initialized by the linear part
*              -  The input vector Ui is the old velocity field
*              -  The input vectors UjLi are the transport directions
*     NONPARAMETRIC VERSION
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21,NNCOF=6)
C
      DOUBLE PRECISION locny
      dimension locny(*)
      DIMENSION A(*),KCOLA(*),KLDA(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRY(NNBAS,NNBAS,4)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      COMMON /COFBAS/ COB(NNBAS,NNCOF)
C
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
      common /fluidp/ pla,ple
      SAVE 
C
      IF (IPRECA.EQ.4) THEN
         DNY=NY
      ELSE
         DNY=0D0
      ENDIF
C     
      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
         CT0=DCMASS/THSTEP
      ELSE
         CT0=0D0
      ENDIF
C
C
      ICUB=ICUBN
C
      DO 1 I = 1,NNDER
1     BDER(I)=.FALSE.
C
      DO 2 I=1,3
2     BDER(I)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
C
C
      DUMAX=0D0
      IF (A2L.EQ.0D0) THEN
       DO 10 IEQ=1,NMT
       DU1=U1L1(IEQ)
       DU2=U1L2(IEQ)
       DUNORM=SQRT(DU1**2+DU2**2)
       DUMAX=MAX(DUMAX,DUNORM)
10     CONTINUE
      ELSE       
       DO 20 IEQ=1,NMT
       DU1=A1L*U1L1(IEQ)+A2L*U2L1(IEQ)
       DU2=A1L*U1L2(IEQ)+A2L*U2L2(IEQ)
       DUNORM=SQRT(DU1**2+DU2**2)
       DUMAX=MAX(DUMAX,DUNORM)
20     CONTINUE
      ENDIF       
C
      IF (DUMAX.LT.1D-8) DUMAX=1D-8
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine local DELTA for streamline-diffusion
      IF (ISTOK.EQ.0) then
         CALL DELTSD(U1L1,U1L2,U2L1,U2L2,A1L,A2L,IEL,
     *               DUMAX,DELTA,KVERT,KMID,DCORVG)
         else
            DELTA=0D0
         endif
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
         DENTRY(JDOFE,JDOFE,1)=0D0
         DENTRY(JDOFE,JDOFE,2)=0D0
         DENTRY(JDOFE,JDOFE,3)=0D0
         DENTRY(JDOFE,JDOFE,4)=0D0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
          DENTRY(JDOFE,IDOFE,1)=0D0
          DENTRY(JDOFE,IDOFE,2)=0D0
          DENTRY(JDOFE,IDOFE,3)=0D0
          DENTRY(JDOFE,IDOFE,4)=0D0
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices
      DO 120 IVE = 1, NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
120   CONTINUE
C
      DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
      DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
      DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
      DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
C
C *** Dummy call - ELE may save arithmetic operations
      CALL ELE(0D0,0D0,-2)
      IF (IER.LT.0) GOTO 99999
C
      avgny=0d0
C *** Loop over all cubature points
      DO 200 ICUBP = 1, NCUBP
C
C *** Cubature points on the reference element
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
C
C *** Jacobian of the bilinear mapping onto the reference element
      DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
      DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
      DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
      DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
      DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
C
C *** Cubature points on actual Element + weights
      XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *  +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
      YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1
     *  +0.5D0*(DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2
      OM=DOMEGA(ICUBP)*DETJ
C
      CALL ELE(XX,YY,-3)
      IF (IER.LT.0) GOTO 99999
C
C *** Evaluation of velocity in cubature points
      IF (DCMASS.NE.0D0) THEN
C
       DU1=0D0
       DU2=0D0
       IF (A2L.EQ.0D0) THEN
        DO 210 JDFL=1,IDFL
        HBAS=DBAS(KDFL(JDFL),1)
        IF (ABS(HBAS).GE.1D-8) THEN
         JDFG=KDFG(JDFL)
         DU1=DU1+U1L1(JDFG)*HBAS
         DU2=DU2+U1L2(JDFG)*HBAS
        ENDIF
210     CONTINUE
       ELSE
        DO 220 JDFL=1,IDFL
        HBAS=DBAS(KDFL(JDFL),1)
        IF (ABS(HBAS).GE.1D-8) THEN
         JDFG=KDFG(JDFL)
         DU1=DU1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
         DU2=DU2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
         ENDIF
220     CONTINUE
       ENDIF
C     and derivatives
          DUX1=0D0
          DUX2=0D0
          IF (A2L.EQ.0D0) THEN
           DO 211 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),2)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUX1=DUX1+U1L1(JDFG)*HBAS
             DUX2=DUX2+U1L2(JDFG)*HBAS
            ENDIF
 211       CONTINUE
          ELSE
           DO 221 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),2)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUX1=DUX1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
             DUX2=DUX2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
            ENDIF
 221       CONTINUE
          ENDIF
          
          DUY1=0D0
          DUY2=0D0
          IF (A2L.EQ.0D0) THEN
           DO 212 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),3)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUY1=DUY1+U1L1(JDFG)*HBAS
             DUY2=DUY2+U1L2(JDFG)*HBAS
            ENDIF
 212       CONTINUE
          ELSE
           DO 222 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),3)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUY1=DUY1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
             DUY2=DUY2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
            ENDIF
 222       CONTINUE
          ENDIF
c
C     
C     norm of symmetric part of velocity gradient
          DUGSQ=DUX1**2+0.5d0*(DUX2+DUY1)**2+DUY2**2

          if (PLA.eq.0d0) then 
           DNY=NY
          else
           DNY=DVISCO(DUGSQ)
          endif
          avgny=avgny+dny/NY
C     
C
C *** Summing up over all pairs of multiindices
       if (IGRAD.eq.1) then
c     use just grad
C ***  Summing up over all pairs of multiindices
       DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2
C
       DO 235 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH4= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
        HSUMI=HBASI2*DU1+HBASI3*DU2
             AH1= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +   DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH4= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +   DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            
 235        CONTINUE
230    CONTINUE
          else
c     use symmetric part of grad
C ***  Summing up over all pairs of multiindices
       DO 240 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2
C
       DO 245 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(2d0*HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH2= DNY*HBASJ2*HBASJ3
             AH3= DNY*HBASJ3*HBASJ2
             AH4= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(HBASJ2**2+2d0*HBASJ3**2)+CT0*HBASJ1**2
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
        HSUMI=HBASI2*DU1+HBASI3*DU2
             AH1= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +   DNY*(2d0*HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH2= DNY*HBASI3*HBASJ2
             AH3= DNY*HBASI2*HBASJ3
             AH4= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +   DNY*(HBASJ2*HBASI2+2d0*HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,2)=DENTRY(JDOFE,IDOFE,2)+OM*AH2
            DENTRY(JDOFE,IDOFE,3)=DENTRY(JDOFE,IDOFE,3)+OM*AH3
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            
 245        CONTINUE
240    CONTINUE
       endif
C
      ELSE
C
C ***  Summing up over 1 pair of multiindices
       DO 250 JDOFE=1,IDFL
       HBASJ1=DBAS(KDFL(JDOFE),1)
C
       DO 260 IDOFE=1,IDFL
       HBASI1=DBAS(KDFL(IDOFE),1)
C
       AH=-1D0/THSTEP*HBASJ1*HBASI1
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH
260    CONTINUE
250    CONTINUE
C
      ENDIF
C
 200  CONTINUE
C     
      locny(iel)=avgny/ncubp
C     
      DO 400 JDOFE=1,IDFL
         DO 410 IDOFE=1,IDFL
            DENTH1=THSTEP*DENTRY(JDOFE,IDOFE,1)
            DENTH2=THSTEP*DENTRY(JDOFE,IDOFE,2)
            DENTH3=THSTEP*DENTRY(JDOFE,IDOFE,3)
            DENTH4=THSTEP*DENTRY(JDOFE,IDOFE,4)
C     
            IF (IDEF.LT.2) THEN
               IA   =KENTRY(JDOFE,IDOFE)
               A(IA)=A(IA)+REAL(DENTH1)
               A(NA+IA)=A(NA+IA)+REAL(DENTH2)
               A(2*NA+IA)=A(2*NA+IA)+REAL(DENTH3)
               A(3*NA+IA)=A(3*NA+IA)+REAL(DENTH4)
            ENDIF
C     
            IF (IDEF.GT.0) THEN 
               IDFG=KDFG(IDOFE)
               JDFG=KDFG(JDOFE)
               D1(JDFG)=D1(JDFG)-DENTH1*U1(IDFG)-DENTH2*U2(IDFG)
               D2(JDFG)=D2(JDFG)-DENTH3*U1(IDFG)-DENTH4*U2(IDFG)
            ENDIF 
C     
 410     CONTINUE
 400  continue
C     
 100  CONTINUE
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE  DELTSD  (U1L1,U1L2,U2L1,U2L2,A1L,A2L,IELE,DUMAX,
     *                    DELTA,KVERT,KMID,DCORVG)
************************************************************************
*     Calculates coefficient Delta for SD
*     RELOC = local RE
*     UNORM=local velocity norm
*     HLOCAL=local mesh width
************************************************************************
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNBAS=21,NNLEV=9,NNARR=299,NNWORK=1)
      PARAMETER (NNDER=6,NNCUBP=36,NNVE=4,NNAB=21)
C
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
C
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /COAUX1/ KDFG,KDFL,IDFL
C

C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
      SAVE
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C
C
C
       DU1=0D0
       DU2=0D0
       dux1=0d0
       dux2=0d0
       duy1=0d0
       duy2=0d0
       DO 10 IDOF=1,IDFL
          DU1=DU1+(A1L*U1L1(KDFG(IDOF))+A2L*U2L1(KDFG(IDOF)))
          DU2=DU2+(A1L*U1L2(KDFG(IDOF))+A2L*U2L2(KDFG(IDOF)))
 10   continue
C
       UNORM=0.25D0*SQRT(DU1**2+DU2**2)
C
       IF (UNORM.LE.1D-8) THEN
        DELTA=0D0
       ELSE

         IF (UPSAM.LT.0D0) THEN
C
C          Jetzt wird das lokale H in Element IEL bestimmt
C
           CALL HWAHL (HLOCAL,UNORM, DU1, DU2, IELE, KVERT,KIMD,DCORVG)
c          AREA=DBLE(VWORK(L(KLAREA(ILEV))+IEL-1))
c          HLOCAL=SQRT(AREA)
c
          DELTA=ABS(UPSAM)*HLOCAL
C
         ELSE
C
C
C        Jetzt wird das lokale H in Element IEL bestimmt
C
c                 AREA=DBLE(VWORK(L(KLAREA(ILEV))+IEL-1))
c                 HLOCAL=SQRT(AREA)
                CALL HWAHL (HLOCAL,UNORM,DU1,DU2,IELE,KVERT,KIMD,DCORVG)
c
            RELOC=UNORM*HLOCAL/NY
            DELTA=UPSAM*HLOCAL/DUMAX*2D0*(RELOC/(1D0+RELOC))       
           ENDIF      
C
        ENDIF
C
      END
c
************************************************************************
      SUBROUTINE addstp(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                  KCOLA,KLDA,KVERT,KMID,DCORVG,ELE,COEF,
     *                  IDEF,ICOMP,DCMASS,locny)
************************************************************************
*     Purpose: -  Adds the SUPG-part on matrix block A after
*                 it was initialized by the linear part
*              -  The input vector Ui is the old velocity field
*              -  The input vectors UjLi are the transport directions
*     PARAMETRIC VERSION
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21)
C
      DOUBLE PRECISION locny
      dimension locny(*)
      DIMENSION A(*),KCOLA(*),KLDA(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRY(NNBAS,NNBAS,4)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
C
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
      common /fluidp/ pla,ple
      SAVE 
C
C
      IF (IPRECA.EQ.4) THEN
         DNY=NY
      ELSE
         DNY=0D0
      ENDIF
C     
      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
         CT0=DCMASS/THSTEP
      ELSE
         CT0=0D0
      ENDIF
C
      ICUB=ICUBN
C
      DO 1 I = 1,NNDER
1     BDER(I)=.FALSE.
C
      DO 2 I=1,3
2     BDER(I)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      ICUBP=ICUB
      CALL ELE(0D0,0D0,-2)
C
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
         DENTRY(JDOFE,JDOFE,1)=0D0
         DENTRY(JDOFE,JDOFE,2)=0D0
         DENTRY(JDOFE,JDOFE,3)=0D0
         DENTRY(JDOFE,JDOFE,4)=0D0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
          DENTRY(JDOFE,IDOFE,1)=0D0
          DENTRY(JDOFE,IDOFE,2)=0D0
          DENTRY(JDOFE,IDOFE,3)=0D0
          DENTRY(JDOFE,IDOFE,4)=0D0
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices
      DO 120 IVE = 1, NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
120   CONTINUE
C
      DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
      DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
      DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
      DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
C
      avgny=0d0
C *** Loop over all cubature points
      DO 200 ICUBP = 1, NCUBP
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
C *** Jacobian of the bilinear mapping onto the reference element
      DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
      DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
      DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
      DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
      DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
      OM=DOMEGA(ICUBP)*DETJ
C
      CALL ELE(XI1,XI2,-3)
      IF (IER.LT.0) GOTO 99999
C
C *** Evaluation of velocity in cubature points
      IF (DCMASS.NE.0D0) THEN
C
       DU1=0D0
       DU2=0D0
       IF (A2L.EQ.0D0) THEN
        DO 210 JDFL=1,IDFL
        HBAS=DBAS(KDFL(JDFL),1)
        IF (ABS(HBAS).GE.1D-8) THEN
         JDFG=KDFG(JDFL)
         DU1=DU1+U1L1(JDFG)*HBAS
         DU2=DU2+U1L2(JDFG)*HBAS
        ENDIF
210     CONTINUE
       ELSE
        DO 220 JDFL=1,IDFL
        HBAS=DBAS(KDFL(JDFL),1)
        IF (ABS(HBAS).GE.1D-8) THEN
         JDFG=KDFG(JDFL)
         DU1=DU1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
         DU2=DU2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
         ENDIF
220     CONTINUE
       ENDIF
C     and derivatives
          DUX1=0D0
          DUX2=0D0
          IF (A2L.EQ.0D0) THEN
           DO 211 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),2)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUX1=DUX1+U1L1(JDFG)*HBAS
             DUX2=DUX2+U1L2(JDFG)*HBAS
            ENDIF
 211       CONTINUE
          ELSE
           DO 221 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),2)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUX1=DUX1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
             DUX2=DUX2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
            ENDIF
 221       CONTINUE
          ENDIF
          
          DUY1=0D0
          DUY2=0D0
          IF (A2L.EQ.0D0) THEN
           DO 212 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),3)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUY1=DUY1+U1L1(JDFG)*HBAS
             DUY2=DUY2+U1L2(JDFG)*HBAS
            ENDIF
 212       CONTINUE
          ELSE
           DO 222 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),3)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUY1=DUY1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
             DUY2=DUY2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
            ENDIF
 222       CONTINUE
          ENDIF
c
C     
C     norm of symmetric part of velocity gradient
          DUGSQ=DUX1**2+0.5d0*(DUX2+DUY1)**2+DUY2**2

          if (PLA.eq.0d0) then 
           DNY=NY
          else
           DNY=DVISCO(DUGSQ)
          endif
          avgny=avgny+dny/NY
C     
C *** Summing up over all pairs of multiindices
       if (IGRAD.eq.1) then
c     use just grad
C ***  Summing up over all pairs of multiindices
       DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
C
       DO 235 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH4= AH1
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
             AH1= DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH4= AH1
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            
 235        CONTINUE
230    CONTINUE
          else
c     use symmetric part of grad
C ***  Summing up over all pairs of multiindices
       DO 240 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
C
       DO 245 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= DNY*(2d0*HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH2= DNY*HBASJ2*HBASJ3
             AH3= DNY*HBASJ3*HBASJ2
             AH4= DNY*(HBASJ2**2+2d0*HBASJ3**2)+CT0*HBASJ1**2
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
             AH1= DNY*(2d0*HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH2= DNY*HBASI3*HBASJ2
             AH3= DNY*HBASI2*HBASJ3
             AH4= DNY*(HBASJ2*HBASI2+2d0*HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,2)=DENTRY(JDOFE,IDOFE,2)+OM*AH2
            DENTRY(JDOFE,IDOFE,3)=DENTRY(JDOFE,IDOFE,3)+OM*AH3
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            
 245        CONTINUE
240    CONTINUE
       endif
C
      ELSE
C
C ***  Summing up over 1 pair of multiindices
       DO 250 JDOFE=1,IDFL
       HBASJ1=DBAS(KDFL(JDOFE),1)
C
       DO 260 IDOFE=1,IDFL
       HBASI1=DBAS(KDFL(IDOFE),1)
C
       AH=-1D0/THSTEP*HBASJ1*HBASI1
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH
260    CONTINUE
250    CONTINUE
C
      ENDIF
C
200   CONTINUE
C     
      locny(iel)=avgny/ncubp
C     
      DO 400 JDOFE=1,IDFL
         DO 410 IDOFE=1,IDFL
            DENTH1=THSTEP*DENTRY(JDOFE,IDOFE,1)
            DENTH2=THSTEP*DENTRY(JDOFE,IDOFE,2)
            DENTH3=THSTEP*DENTRY(JDOFE,IDOFE,3)
            DENTH4=THSTEP*DENTRY(JDOFE,IDOFE,4)
C     
            IF (IDEF.LT.2) THEN
               IA   =KENTRY(JDOFE,IDOFE)
               A(IA)=A(IA)+REAL(DENTH1)
               A(NA+IA)=A(NA+IA)+REAL(DENTH2)
               A(2*NA+IA)=A(2*NA+IA)+REAL(DENTH3)
               A(3*NA+IA)=A(3*NA+IA)+REAL(DENTH4)
            ENDIF
C     
            IF (IDEF.GT.0) THEN 
               IDFG=KDFG(IDOFE)
               JDFG=KDFG(JDOFE)
               D1(JDFG)=D1(JDFG)-DENTH1*U1(IDFG)-DENTH2*U2(IDFG)
               D2(JDFG)=D2(JDFG)-DENTH3*U1(IDFG)-DENTH4*U2(IDFG)
            ENDIF 
C     
 410     CONTINUE
 400  continue
C     
 100  CONTINUE
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE addstn(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                  KCOLA,KLDA,KVERT,KMID,DCORVG,ELE,COEF,
     *                  IDEF,ICOMP,DCMASS,locny)
************************************************************************
*     Purpose: -  Adds the SUPG-part on matrix block A after
*                 it was initialized by the linear part
*              -  The input vector Ui is the old velocity field
*              -  The input vectors UjLi are the transport directions
*     NONPARAMETRIC VERSION
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21,NNCOF=6)
C
      DOUBLE PRECISION locny
      dimension locny(*)
      DIMENSION A(*),KCOLA(*),KLDA(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRY(NNBAS,NNBAS,4)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      COMMON /COFBAS/ COB(NNBAS,NNCOF)
C
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
      common /fluidp/ pla,ple
C
      SAVE 
c
      external dvisco
C
      IF (IPRECA.EQ.4) THEN
         DNY=NY
      ELSE
         DNY=0D0
      ENDIF
C     
      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
         CT0=DCMASS/THSTEP
      ELSE
         CT0=0D0
      ENDIF
C
C
      ICUB=ICUBN
C
      DO 1 I = 1,NNDER
1     BDER(I)=.FALSE.
C
      DO 2 I=1,3
2     BDER(I)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
C
C
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
         DENTRY(JDOFE,JDOFE,1)=0D0
         DENTRY(JDOFE,JDOFE,2)=0D0
         DENTRY(JDOFE,JDOFE,3)=0D0
         DENTRY(JDOFE,JDOFE,4)=0D0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
          DENTRY(JDOFE,IDOFE,1)=0D0
          DENTRY(JDOFE,IDOFE,2)=0D0
          DENTRY(JDOFE,IDOFE,3)=0D0
          DENTRY(JDOFE,IDOFE,4)=0D0
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices
      DO 120 IVE = 1, NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
120   CONTINUE
C
      DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
      DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
      DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
      DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
C
C *** Dummy call - ELE may save arithmetic operations
      CALL ELE(0D0,0D0,-2)
      IF (IER.LT.0) GOTO 99999
C
      avgny=0d0
C *** Loop over all cubature points
      DO 200 ICUBP = 1, NCUBP
C
C *** Cubature points on the reference element
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
C
C *** Jacobian of the bilinear mapping onto the reference element
      DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
      DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
      DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
      DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
      DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
C
C *** Cubature points on actual Element + weights
      XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *  +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
      YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1
     *  +0.5D0*(DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2
      OM=DOMEGA(ICUBP)*DETJ
C
      CALL ELE(XX,YY,-3)
      IF (IER.LT.0) GOTO 99999
C
C *** Evaluation of velocity in cubature points
      IF (DCMASS.NE.0D0) THEN
C
       DU1=0D0
       DU2=0D0
       IF (A2L.EQ.0D0) THEN
        DO 210 JDFL=1,IDFL
        HBAS=DBAS(KDFL(JDFL),1)
        IF (ABS(HBAS).GE.1D-8) THEN
         JDFG=KDFG(JDFL)
         DU1=DU1+U1L1(JDFG)*HBAS
         DU2=DU2+U1L2(JDFG)*HBAS
        ENDIF
210     CONTINUE
       ELSE
        DO 220 JDFL=1,IDFL
        HBAS=DBAS(KDFL(JDFL),1)
        IF (ABS(HBAS).GE.1D-8) THEN
         JDFG=KDFG(JDFL)
         DU1=DU1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
         DU2=DU2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
         ENDIF
220     CONTINUE
       ENDIF
C     and derivatives
          DUX1=0D0
          DUX2=0D0
          IF (A2L.EQ.0D0) THEN
           DO 211 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),2)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUX1=DUX1+U1L1(JDFG)*HBAS
             DUX2=DUX2+U1L2(JDFG)*HBAS
            ENDIF
 211       CONTINUE
          ELSE
           DO 221 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),2)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUX1=DUX1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
             DUX2=DUX2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
            ENDIF
 221       CONTINUE
          ENDIF
          
          DUY1=0D0
          DUY2=0D0
          IF (A2L.EQ.0D0) THEN
           DO 212 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),3)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUY1=DUY1+U1L1(JDFG)*HBAS
             DUY2=DUY2+U1L2(JDFG)*HBAS
            ENDIF
 212       CONTINUE
          ELSE
           DO 222 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),3)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUY1=DUY1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
             DUY2=DUY2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
            ENDIF
 222       CONTINUE
          ENDIF
c
C     
C     norm of symmetric part of velocity gradient
          DUGSQ=DUX1**2+0.5d0*(DUX2+DUY1)**2+DUY2**2

          if (PLA.eq.0d0) then 
           DNY=NY
          else
           DNY=DVISCO(DUGSQ)
          endif
          avgny=avgny+dny/NY
C
C
C *** Summing up over all pairs of multiindices
       if (IGRAD.eq.1) then
c     use just grad
C ***  Summing up over all pairs of multiindices
       DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
C
       DO 235 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH4= AH1
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
             AH1= DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH4= AH1
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            
 235        CONTINUE
230    CONTINUE
          else
c     use symmetric part of grad
C ***  Summing up over all pairs of multiindices
       DO 240 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
C
       DO 245 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= DNY*(2d0*HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH2= DNY*HBASJ2*HBASJ3
             AH3= DNY*HBASJ3*HBASJ2
             AH4= DNY*(HBASJ2**2+2d0*HBASJ3**2)+CT0*HBASJ1**2
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
             AH1= DNY*(2d0*HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH2= DNY*HBASI3*HBASJ2
             AH3= DNY*HBASI2*HBASJ3
             AH4= DNY*(HBASJ2*HBASI2+2d0*HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,2)=DENTRY(JDOFE,IDOFE,2)+OM*AH2
            DENTRY(JDOFE,IDOFE,3)=DENTRY(JDOFE,IDOFE,3)+OM*AH3
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            
 245        CONTINUE
240    CONTINUE
       endif
C
      ELSE
C
C ***  Summing up over 1 pair of multiindices
       DO 250 JDOFE=1,IDFL
       HBASJ1=DBAS(KDFL(JDOFE),1)
C
       DO 260 IDOFE=1,IDFL
       HBASI1=DBAS(KDFL(IDOFE),1)
C
       AH=-1D0/THSTEP*HBASJ1*HBASI1
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH
260    CONTINUE
250    CONTINUE
C
      ENDIF
C     
 200  CONTINUE
C     
      locny(iel)=avgny/ncubp
C     
      DO 400 JDOFE=1,IDFL
         DO 410 IDOFE=1,IDFL
            DENTH1=THSTEP*DENTRY(JDOFE,IDOFE,1)
            DENTH2=THSTEP*DENTRY(JDOFE,IDOFE,2)
            DENTH3=THSTEP*DENTRY(JDOFE,IDOFE,3)
            DENTH4=THSTEP*DENTRY(JDOFE,IDOFE,4)
C     
            IF (IDEF.LT.2) THEN
               IA   =KENTRY(JDOFE,IDOFE)
               A(IA)=A(IA)+REAL(DENTH1)
               A(NA+IA)=A(NA+IA)+REAL(DENTH2)
               A(2*NA+IA)=A(2*NA+IA)+REAL(DENTH3)
               A(3*NA+IA)=A(3*NA+IA)+REAL(DENTH4)
            ENDIF
C     
            IF (IDEF.GT.0) THEN 
               IDFG=KDFG(IDOFE)
               JDFG=KDFG(JDOFE)
               D1(JDFG)=D1(JDFG)-DENTH1*U1(IDFG)-DENTH2*U2(IDFG)
               D2(JDFG)=D2(JDFG)-DENTH3*U1(IDFG)-DENTH4*U2(IDFG)
            ENDIF 
C     
 410     CONTINUE
 400  continue
C     
 100  CONTINUE
C
99999 END
C
C
C
