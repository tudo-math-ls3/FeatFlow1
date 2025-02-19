************************************************************************
      DOUBLE PRECISION FUNCTION   COEFST (X,Y,IA,IB,IBLOC,BFIRST)
*
*     Coefficient for the Stokes-block
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
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
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      SAVE
C
C
C
      IF ((IA.EQ.1).AND.(IB.EQ.1)) THEN
       COEFST=1D0
      ELSE
       COEFST=NY
      ENDIF
C
      END
C
************************************************************************
      DOUBLE PRECISION FUNCTION COEFFN(IA,IB,IBLOC,U1L1,U1L2,U2L1,U2L2,
     *                                 A1L,A2L,DELTA,DCMASS)
*
*     Coefficient for the convective-block
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS),U1L1(*),U1L2(*),U2L1(*),U2L2(*)
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
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      SAVE
C
C
      DU1=0D0
      DU2=0D0
C 
      IF (DCMASS.NE.0D0) THEN
C
C
         DO 10 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),1)
            DU1=DU1+(A1L*U1L1(KDFG(JDFL))+A2L*U2L1(KDFG(JDFL)))*HBAS
            DU2=DU2+(A1L*U1L2(KDFG(JDFL))+A2L*U2L2(KDFG(JDFL)))*HBAS
 10      CONTINUE
C
       IF ((IA.EQ.1).AND.(IB.EQ.2)) COEFFN=DU1
       IF ((IA.EQ.1).AND.(IB.EQ.3)) COEFFN=DU2
       IF ((IA.EQ.2).AND.(IB.EQ.3)) COEFFN=DELTA*DU1*DU2
       IF ((IA.EQ.3).AND.(IB.EQ.2)) COEFFN=DELTA*DU1*DU2
C
       IF ((IA.EQ.1).AND.(IB.EQ.1)) COEFFN=DCMASS/THSTEP
C
       IF (IPRECA.EQ.4) THEN      
        IF ((IA.EQ.2).AND.(IB.EQ.2)) COEFFN=DELTA*DU1**2+NY
        IF ((IA.EQ.3).AND.(IB.EQ.3)) COEFFN=DELTA*DU2**2+NY
       ELSE
        IF ((IA.EQ.2).AND.(IB.EQ.2)) COEFFN=DELTA*DU1**2
        IF ((IA.EQ.3).AND.(IB.EQ.3)) COEFFN=DELTA*DU2**2
       ENDIF
C
      ELSE
C
       COEFFN=0D0
       IF ((IA.EQ.1).AND.(IB.EQ.1)) COEFFN=-1D0/THSTEP
C
      ENDIF
C
c      aux=COEFNN(IA,IB,IBLOC,U1L1,U1L2,U2L1,U2L2,
c     *                                 A1L,A2L,DELTA,DCMASS)
c      print *,'original',COEFFN,aux
C
      END
C
C
************************************************************************
      DOUBLE PRECISION FUNCTION COEFNN(IA,IB,IBLOC,U1L1,U1L2,U2L1,U2L2,
     *                                 A1L,A2L,DELTA,DCMASS)
*
*     Coefficient for the convective-block, nonconstant viscosity
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS),U1L1(*),U1L2(*),U2L1(*),U2L2(*)
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
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
c
      common /fluidp/ PLA,PLE
      SAVE
C
C
      DU1=0D0
      DU2=0D0
      dux1=0d0
      dux2=0d0
      duy1=0d0
      duy2=0d0
C
      COEFNN=0d0

      IF (DCMASS.NE.0D0) THEN
C     
         DO 10 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),1)
            DU1=DU1+(A1L*U1L1(KDFG(JDFL))+A2L*U2L1(KDFG(JDFL)))*HBAS
            DU2=DU2+(A1L*U1L2(KDFG(JDFL))+A2L*U2L2(KDFG(JDFL)))*HBAS
            HBAS=DBAS(KDFL(JDFL),2)
            DUX1=DUX1+(A1L*U1L1(KDFG(JDFL))+A2L*U2L1(KDFG(JDFL)))*HBAS
            DUX2=DUX2+(A1L*U1L2(KDFG(JDFL))+A2L*U2L2(KDFG(JDFL)))*HBAS
            HBAS=DBAS(KDFL(JDFL),3)
            DUY1=DUY1+(A1L*U1L1(KDFG(JDFL))+A2L*U2L1(KDFG(JDFL)))*HBAS
            DUY2=DUY2+(A1L*U1L2(KDFG(JDFL))+A2L*U2L2(KDFG(JDFL)))*HBAS
 10      CONTINUE
C
         if(PLA.eq.0d0) then 
            DNY=NY
         else
            DUGSQ=DUX1**2+0.5d0*(DUX2+DUY1)**2+DUY2**2
            DNY=NY*(PLE+sqrt(DUGSQ))**(-PLA)
         endif

         if(IBLOC.eq.1) then
            IF ((IA.EQ.1).AND.(IB.EQ.2)) COEFNN=DU1
            IF ((IA.EQ.1).AND.(IB.EQ.3)) COEFNN=DU2
            IF ((IA.EQ.2).AND.(IB.EQ.3)) COEFNN=DELTA*DU1*DU2
            IF ((IA.EQ.3).AND.(IB.EQ.2)) COEFNN=DELTA*DU1*DU2
C     
            IF ((IA.EQ.1).AND.(IB.EQ.1)) COEFNN=DCMASS/THSTEP

            if(igrad.eq.0) then
            IF ((IA.EQ.2).AND.(IB.EQ.2)) COEFNN=DELTA*DU1**2+2d0*DNY
            IF ((IA.EQ.3).AND.(IB.EQ.3)) COEFNN=DELTA*DU2**2+DNY
            else
            IF ((IA.EQ.2).AND.(IB.EQ.2)) COEFNN=DELTA*DU1**2+DNY
            IF ((IA.EQ.3).AND.(IB.EQ.3)) COEFNN=DELTA*DU2**2+DNY
            endif
         endif
         
         if(IBLOC.eq.2) then
            if(igrad.eq.0) then
            IF ((IA.EQ.2).AND.(IB.EQ.2)) COEFNN=0d0
            IF ((IA.EQ.3).AND.(IB.EQ.3)) COEFNN=DNY
            else
            IF ((IA.EQ.2).AND.(IB.EQ.2)) COEFNN=0d0
            IF ((IA.EQ.3).AND.(IB.EQ.3)) COEFNN=0d0
            endif
         endif
         
         if(IBLOC.eq.3) then
            if(igrad.eq.0) then
            IF ((IA.EQ.2).AND.(IB.EQ.2)) COEFNN=DNY
            IF ((IA.EQ.3).AND.(IB.EQ.3)) COEFNN=0d0
            else
            IF ((IA.EQ.2).AND.(IB.EQ.2)) COEFNN=0d0
            IF ((IA.EQ.3).AND.(IB.EQ.3)) COEFNN=0d0
            endif
         endif
         
         if(IBLOC.eq.4) then
            IF ((IA.EQ.1).AND.(IB.EQ.2)) COEFNN=DU1
            IF ((IA.EQ.1).AND.(IB.EQ.3)) COEFNN=DU2
            IF ((IA.EQ.2).AND.(IB.EQ.3)) COEFNN=DELTA*DU1*DU2
            IF ((IA.EQ.3).AND.(IB.EQ.2)) COEFNN=DELTA*DU1*DU2
C     
            IF ((IA.EQ.1).AND.(IB.EQ.1)) COEFNN=DCMASS/THSTEP
            
            if(igrad.eq.0) then
            IF ((IA.EQ.2).AND.(IB.EQ.2)) COEFNN=DELTA*DU1**2+DNY
            IF ((IA.EQ.3).AND.(IB.EQ.3)) COEFNN=DELTA*DU2**2+2d0*DNY
            else
            IF ((IA.EQ.2).AND.(IB.EQ.2)) COEFNN=DELTA*DU1**2+DNY
            IF ((IA.EQ.3).AND.(IB.EQ.3)) COEFNN=DELTA*DU2**2+DNY
            endif
         endif
C     
         if(IBLOC.eq.5) then
            IF ((IA.EQ.1).AND.(IB.EQ.2)) COEFNN=DU1
            IF ((IA.EQ.1).AND.(IB.EQ.3)) COEFNN=DU2
            IF ((IA.EQ.2).AND.(IB.EQ.3)) COEFNN=DELTA*DU1*DU2
            IF ((IA.EQ.3).AND.(IB.EQ.2)) COEFNN=DELTA*DU1*DU2
C     
            IF ((IA.EQ.1).AND.(IB.EQ.1)) COEFNN=DCMASS/THSTEP
            
            IF ((IA.EQ.2).AND.(IB.EQ.2)) COEFNN=DELTA*DU1**2+DNY
            IF ((IA.EQ.3).AND.(IB.EQ.3)) COEFNN=DELTA*DU2**2+DNY
         endif
      ELSE
C     
         COEFNN=0D0
         if ((ibloc.eq.1).or.(ibloc.eq.4)) then
            IF ((IA.EQ.1).AND.(IB.EQ.1)) COEFNN=-1D0/THSTEP
         endif
C     
      ENDIF
C
C     
99999 END
C
************************************************************************
      DOUBLE PRECISION FUNCTION   COEFFB (X,Y,IA,IB,IBLOC,BFIRST)
*
*     Coefficient for the B1/B2-blocks
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
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
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      SAVE
C
      COEFFB= -1.D0
C
      END
C
************************************************************************
      DOUBLE PRECISION FUNCTION    RHS  (X,Y,IA,IBLOC,BFIRST)
*
*     Right hand side yielding the exact solution UE
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
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
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      SAVE
C
      RHS=FDATIN(5,IBLOC,X,Y,TIMENS,RE)
C
      END
C
C
************************************************************************
      DOUBLE PRECISION FUNCTION    UE  (X,Y,IBLOC)
*
*     Exact solution - velocity,  also for boundary conditions
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
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
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      SAVE
C
      UE=FDATIN(1,IBLOC,X,Y,TIMENS,RE)
C
      END
C
************************************************************************
      DOUBLE PRECISION FUNCTION    PE  (X,Y)
*
*     Exact solution - pressure, only for error analysis
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
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
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      SAVE
C     
      PE=FDATIN(4,IBLOC,X,Y,TIMENS,RE)
C
      END
C
*************************************************************************
      DOUBLE PRECISION FUNCTION    UEX(X,Y,IBLOC)
C
C     x-Ableitung der exakten Loesung, only for error analysis
*************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
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
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      SAVE
C
      UEX=FDATIN(2,IBLOC,X,Y,TIMENS,RE)
C
      END
C
*************************************************************************
      DOUBLE PRECISION FUNCTION    UEY(X,Y,IBLOC)
C
C     y-Ableitung der exakten Loesung, only for error analysis
*************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
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
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      SAVE
C
      UEY=FDATIN(3,IBLOC,X,Y,TIMENS,RE)
C
      END
