************************************************************************
      SUBROUTINE SUPWDG(U1L1,U1L2,U1L3,U2L1,U2L2,U2L3,A1L,A2L,U1,U2,U3,
     *                  D1,D2,D3,A,NA,KCOLA,KLDA,KVERT,KAREA,KEDGE,
     *                  KINT,DCORVG,ELE,COEFFN,IDEF,DCMASS)
************************************************************************
*     Purpose: -  Adds the SUPG-part on matrix block A after
*                 it was initialized by the linear part
*              -  The input vector Ui is the old velocity field
*              -  The input vectors UjLi are the transport directions
*     PARAMETRIC VERSION
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION A(*),KCOLA(*),KLDA(*)
      DIMENSION U1L1(*),U1L2(*),U1L3(*),U2L1(*),U2L2(*),U2L3(*)
      DIMENSION U1(*),U2(*),U3(*),D1(*),D2(*),D3(*),KINT(*)
      DIMENSION DCORVG(3,*),KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRY(NNBAS,NNBAS)
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
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC
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
      ICUB=ICUBN
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
      DO 1 I = 1,NNDER
1     BDER(I)=.FALSE.
C
      DO 2 I=1,4
2     BDER(I)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
ccc      CALL ZTIME(TTT1)
ccc      TTU0=TTU0+TTT1-TTT0
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
ccc      CALL ZTIME(TTT0)
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
ccc      CALL ZTIME(TTT1)
ccc      TTU3=TTU3+TTT1-TTT0
C
ccc      CALL ZTIME(TTT0)
      DUMAX=0D0
      IF (A2L.EQ.0D0) THEN
       DO 10 IEQ=1,NAT
       DU1=U1L1(IEQ)
       DU2=U1L2(IEQ)
       DU3=U1L3(IEQ)
       DUNORM=SQRT(DU1**2+DU2**2+DU3**2)
       DUMAX=MAX(DUMAX,DUNORM)
10     CONTINUE
      ELSE       
       DO 20 IEQ=1,NAT
       DU1=A1L*U1L1(IEQ)+A2L*U2L1(IEQ)
       DU2=A1L*U1L2(IEQ)+A2L*U2L2(IEQ)
       DU3=A1L*U1L3(IEQ)+A2L*U2L3(IEQ)
       DUNORM=SQRT(DU1**2+DU2**2+DU3**2)
       DUMAX=MAX(DUMAX,DUNORM)
20     CONTINUE
      ENDIF       
C
      IF (DUMAX.LT.1D-8) DUMAX=1D-8
C
C
ccc      CALL ZTIME(TTT1)
ccc      TTU0=TTU0+TTT1-TTT0
C *** Loop over all elements
      DO 100 IEL=1,NEL
      ILINT=KINT(IEL)
C      ILINT=0
C
ccc      CALL ZTIME(TTT0)
C *** Set zero elements of Jacobian for axiparallel grid
      IF (ILINT.EQ.2) THEN
       DJAC(1,3)=0D0
       DJAC(2,3)=0D0
       DJAC(3,1)=0D0
       DJAC(3,2)=0D0
      ENDIF
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine local DELTA for streamline-diffusion
      CALL DELTSD(U1L1,U1L2,U1L3,U2L1,U2L2,U2L3,A1L,A2L,IEL,DUMAX,DELTA)
      IF (ISTOK.EQ.1) DELTA=0D0
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
      DENTRY(JDOFE,JDOFE)=0D0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
      DENTRY(JDOFE,IDOFE)=0D0
111   CONTINUE
110   CONTINUE
c
C
C *** Evaluation of coordinates of the vertices
      DO 120 IVE = 1, NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
120   CONTINUE
ccc      CALL ZTIME(TTT1)
ccc      TTU11=TTU11+TTT1-TTT0
ccc      CALL ZTIME(TTT0)
C
      IF (ILINT.EQ.2) THEN
       DJ11=(DX(2)+DX(4))*Q2
       DJ12=(DY(2)+DY(4))*Q2
       DJ13=(DZ(1)+DZ(5))*Q2
       DJAC(1,1)=(-DX(1)+DX(2))*Q2
       DJAC(2,1)=(-DY(1)+DY(2))*Q2
       DJAC(1,2)=(-DX(1)+DX(4))*Q2
       DJAC(2,2)=(-DY(1)+DY(4))*Q2
       DJAC(3,3)=(-DZ(1)+DZ(5))*Q2
       DETJ=DJAC(3,3)*(DJAC(1,1)*DJAC(2,2)-DJAC(2,1)*DJAC(1,2))
      ELSE IF (ILINT.EQ.1) THEN
       DJ11=(DX(1)+DX(2)+DX(3)+DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJ12=(DY(1)+DY(2)+DY(3)+DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJ13=(DZ(1)+DZ(2)+DZ(3)+DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DJAC(1,1)=(-DX(1)+DX(2)+DX(3)-DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
       DJAC(2,1)=(-DY(1)+DY(2)+DY(3)-DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
       DJAC(3,1)=(-DZ(1)+DZ(2)+DZ(3)-DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
       DJAC(1,2)=(-DX(1)-DX(2)+DX(3)+DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
       DJAC(2,2)=(-DY(1)-DY(2)+DY(3)+DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
       DJAC(3,2)=(-DZ(1)-DZ(2)+DZ(3)+DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
       DJAC(1,3)=(-DX(1)-DX(2)-DX(3)-DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJAC(2,3)=(-DY(1)-DY(2)-DY(3)-DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJAC(3,3)=(-DZ(1)-DZ(2)-DZ(3)-DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      ELSE
       DJ11=( DX(1)+DX(2)+DX(3)+DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJ12=( DY(1)+DY(2)+DY(3)+DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJ13=( DZ(1)+DZ(2)+DZ(3)+DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DJ21=(-DX(1)+DX(2)+DX(3)-DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
       DJ22=(-DY(1)+DY(2)+DY(3)-DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
       DJ23=(-DZ(1)+DZ(2)+DZ(3)-DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
       DJ31=(-DX(1)-DX(2)+DX(3)+DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
       DJ32=(-DY(1)-DY(2)+DY(3)+DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
       DJ33=(-DZ(1)-DZ(2)+DZ(3)+DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
       DJ41=(-DX(1)-DX(2)-DX(3)-DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJ42=(-DY(1)-DY(2)-DY(3)-DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJ43=(-DZ(1)-DZ(2)-DZ(3)-DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DJ51=( DX(1)-DX(2)+DX(3)-DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
       DJ52=( DY(1)-DY(2)+DY(3)-DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
       DJ53=( DZ(1)-DZ(2)+DZ(3)-DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
       DJ61=( DX(1)-DX(2)-DX(3)+DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
       DJ62=( DY(1)-DY(2)-DY(3)+DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
       DJ63=( DZ(1)-DZ(2)-DZ(3)+DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
       DJ71=( DX(1)+DX(2)-DX(3)-DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
       DJ72=( DY(1)+DY(2)-DY(3)-DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
       DJ73=( DZ(1)+DZ(2)-DZ(3)-DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
       DJ81=(-DX(1)+DX(2)-DX(3)+DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
       DJ82=(-DY(1)+DY(2)-DY(3)+DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
       DJ83=(-DZ(1)+DZ(2)-DZ(3)+DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
      ENDIF
ccc      CALL ZTIME(TTT1)
ccc      TTU2=TTU2+TTT1-TTT0
C
C *** Loop over all cubature points
      DO 200 ICUBP = 1, NCUBP
ccc      CALL ZTIME(TTT0)
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the bilinear mapping onto the reference element
      IF (ILINT.EQ.0) THEN
       DJAC(1,1)=DJ21+DJ51*XI2+DJ61*XI3+DJ81*XI2*XI3
       DJAC(1,2)=DJ31+DJ51*XI1+DJ71*XI3+DJ81*XI1*XI3
       DJAC(1,3)=DJ41+DJ61*XI1+DJ71*XI2+DJ81*XI1*XI2
       DJAC(2,1)=DJ22+DJ52*XI2+DJ62*XI3+DJ82*XI2*XI3
       DJAC(2,2)=DJ32+DJ52*XI1+DJ72*XI3+DJ82*XI1*XI3
       DJAC(2,3)=DJ42+DJ62*XI1+DJ72*XI2+DJ82*XI1*XI2
       DJAC(3,1)=DJ23+DJ53*XI2+DJ63*XI3+DJ83*XI2*XI3
       DJAC(3,2)=DJ33+DJ53*XI1+DJ73*XI3+DJ83*XI1*XI3
       DJAC(3,3)=DJ43+DJ63*XI1+DJ73*XI2+DJ83*XI1*XI2
       DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      ENDIF
      OM=DOMEGA(ICUBP)*ABS(DETJ)
ccc      CALL ZTIME(TTT1)
ccc      TTU2=TTU2+TTT1-TTT0
C
ccc      CALL ZTIME(TTT0)
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
ccc      CALL ZTIME(TTT1)
ccc      TTU4=TTU4+TTT1-TTT0
C
ccc      CALL ZTIME(TTT0)
C *** Evaluation of velocity in cubature points
      IF (DCMASS.NE.0D0) THEN
C
       DU1=0D0
       DU2=0D0
       DU3=0D0
       IF (A2L.EQ.0D0) THEN
        DO 210 JDFL=1,IDFL
        HBAS=DBAS(1,KDFL(JDFL),1)
        IF (ABS(HBAS).GE.1D-8) THEN
         JDFG=KDFG(JDFL)
         DU1=DU1+U1L1(JDFG)*HBAS
         DU2=DU2+U1L2(JDFG)*HBAS
         DU3=DU3+U1L3(JDFG)*HBAS
        ENDIF
210     CONTINUE
       ELSE
        DO 220 JDFL=1,IDFL
        HBAS=DBAS(1,KDFL(JDFL),1)
        IF (ABS(HBAS).GE.1D-8) THEN
         JDFG=KDFG(JDFL)
         DU1=DU1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
         DU2=DU2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
         DU3=DU3+(A1L*U1L3(JDFG)+A2L*U2L3(JDFG))*HBAS
         ENDIF
220     CONTINUE
       ENDIF
C
C ***  Summing up over all pairs of multiindices
       DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(1,JDOFEH,1)
       HBASJ2=DBAS(1,JDOFEH,2)
       HBASJ3=DBAS(1,JDOFEH,3)
       HBASJ4=DBAS(1,JDOFEH,4)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2+HBASJ4*DU3
C
       DO 240 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
        AH= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *     +  DNY*(HBASJ2**2+HBASJ3**2+HBASJ4**2)+CT0*HBASJ1**2
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(1,IDOFEH,1)
        HBASI2=DBAS(1,IDOFEH,2)
        HBASI3=DBAS(1,IDOFEH,3)
        HBASI4=DBAS(1,IDOFEH,4)
        HSUMI=HBASI2*DU1+HBASI3*DU2+HBASI4*DU3
        AH= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *    +   DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3+HBASJ4*HBASI4)
     *    +   CT0*HBASJ1*HBASI1
       ENDIF
C
       DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
240    CONTINUE
230    CONTINUE
ccc       CALL ZTIME(TTT1)
ccc       TTU12=TTU12+TTT1-TTT0
C
      ELSE
C
C ***  Summing up over 1 pair of multiindices
       DO 250 JDOFE=1,IDFL
       HBASJ1=DBAS(1,KDFL(JDOFE),1)
C
       DO 260 IDOFE=1,IDFL
       HBASI1=DBAS(1,KDFL(IDOFE),1)
C
       AH=-1D0/THSTEP*HBASJ1*HBASI1
C
       DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
260    CONTINUE
250    CONTINUE
ccc       CALL ZTIME(TTT1)
ccc       TTU12=TTU12+TTT1-TTT0
C
      ENDIF
C
200   CONTINUE
C
C
ccc      CALL ZTIME(TTT0)
      DO 400 JDOFE=1,IDFL
      DO 400 IDOFE=1,IDFL
      DENTH=THSTEP*DENTRY(JDOFE,IDOFE)
C
      IF (IDEF.LT.2) THEN
       IA   =KENTRY(JDOFE,IDOFE)
       A(IA)=A(IA)+REAL(DENTH)
      ENDIF
C
      IF (IDEF.GT.0) THEN 
       IDFG=KDFG(IDOFE)
       JDFG=KDFG(JDOFE)
       D1(JDFG)= D1(JDFG)-DENTH*U1(IDFG)
       D2(JDFG)= D2(JDFG)-DENTH*U2(IDFG)
       D3(JDFG)= D3(JDFG)-DENTH*U3(IDFG)
      ENDIF 
C
400   CONTINUE
C
ccc      CALL ZTIME(TTT1)
ccc      TTU13=TTU13+TTT1-TTT0
100   CONTINUE
C
ccc      write(6,*) TTU0,TTU11,TTU12,TTU13,TTU2,TTU3,TTU4,NCUBP
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE SUPWNP(U1L1,U1L2,U1L3,U2L1,U2L2,U2L3,A1L,A2L,U1,U2,U3,
     *                  D1,D2,D3,A,NA,KCOLA,KLDA,KVERT,KAREA,KEDGE,
     *                  KINT,DCORVG,ELE,COEFFN,IDEF,DCMASS)
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
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
C
      DIMENSION A(*),KCOLA(*),KLDA(*)
      DIMENSION U1L1(*),U1L2(*),U1L3(*),U2L1(*),U2L2(*),U2L3(*)
      DIMENSION U1(*),U2(*),U3(*),D1(*),D2(*),D3(*),KINT(*)
      DIMENSION DCORVG(3,*),KVERT(NNVE,*),KAREA(NNAE,*),KEDGE(NNEE,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS),DENTRY(NNBAS,NNBAS)
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
      COMMON /COFBAS/ COB(NNBAS,NNCOF)
C
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC
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
ccc     CALL ZTIME(TTT0)
      ICUB=ICUBN
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
      DO 1 I = 1,NNDER
1     BDER(I)=.FALSE.
C
      DO 2 I=1,4
2     BDER(I)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
C
      DUMAX=0D0
      IF (A2L.EQ.0D0) THEN
       DO 10 IEQ=1,NAT
       DU1=U1L1(IEQ)
       DU2=U1L2(IEQ)
       DU3=U1L3(IEQ)
       DUNORM=SQRT(DU1**2+DU2**2+DU3**2)
       DUMAX=MAX(DUMAX,DUNORM)
10     CONTINUE
      ELSE       
       DO 20 IEQ=1,NAT
       DU1=A1L*U1L1(IEQ)+A2L*U2L1(IEQ)
       DU2=A1L*U1L2(IEQ)+A2L*U2L2(IEQ)
       DU3=A1L*U1L3(IEQ)+A2L*U2L3(IEQ)
       DUNORM=SQRT(DU1**2+DU2**2+DU3**2)
       DUMAX=MAX(DUMAX,DUNORM)
20     CONTINUE
      ENDIF       
C
      IF (DUMAX.LT.1D-8) DUMAX=1D-8
C
C
ccc      CALL ZTIME(TTT1)
ccc      TTU0=TTU0+TTT1-TTT0
C *** Loop over all elements
      DO 100 IEL=1,NEL
      ILINT=KINT(IEL)
C      ILINT=0
C
ccc      CALL ZTIME(TTT0)
C *** Set zero elements of Jacobian for axiparallel grid
      IF (ILINT.EQ.2) THEN
       DJAC(1,3)=0D0
       DJAC(2,3)=0D0
       DJAC(3,1)=0D0
       DJAC(3,2)=0D0
      ENDIF
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine local DELTA for streamline-diffusion
      CALL DELTSD(U1L1,U1L2,U1L3,U2L1,U2L2,U2L3,A1L,A2L,IEL,DUMAX,DELTA)
      IF (ISTOK.EQ.1) DELTA=0D0
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
      DENTRY(JDOFE,JDOFE)=0D0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
      DENTRY(JDOFE,IDOFE)=0D0
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices
      DO 120 IVE = 1, NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
120   CONTINUE
ccc      CALL ZTIME(TTT1)
ccc      TTU11=TTU11+TTT1-TTT0
ccc      CALL ZTIME(TTT0)
C
      IF (ILINT.EQ.2) THEN
       DJ11=(DX(2)+DX(4))*Q2
       DJ12=(DY(2)+DY(4))*Q2
       DJ13=(DZ(1)+DZ(5))*Q2
       DJAC(1,1)=(-DX(1)+DX(2))*Q2
       DJAC(2,1)=(-DY(1)+DY(2))*Q2
       DJAC(1,2)=(-DX(1)+DX(4))*Q2
       DJAC(2,2)=(-DY(1)+DY(4))*Q2
       DJAC(3,3)=(-DZ(1)+DZ(5))*Q2
       DETJ=DJAC(3,3)*(DJAC(1,1)*DJAC(2,2)-DJAC(2,1)*DJAC(1,2))
      ELSE IF (ILINT.EQ.1) THEN
       DJ11=(DX(1)+DX(2)+DX(3)+DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJ12=(DY(1)+DY(2)+DY(3)+DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJ13=(DZ(1)+DZ(2)+DZ(3)+DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DJAC(1,1)=(-DX(1)+DX(2)+DX(3)-DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
       DJAC(2,1)=(-DY(1)+DY(2)+DY(3)-DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
       DJAC(3,1)=(-DZ(1)+DZ(2)+DZ(3)-DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
       DJAC(1,2)=(-DX(1)-DX(2)+DX(3)+DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
       DJAC(2,2)=(-DY(1)-DY(2)+DY(3)+DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
       DJAC(3,2)=(-DZ(1)-DZ(2)+DZ(3)+DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
       DJAC(1,3)=(-DX(1)-DX(2)-DX(3)-DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJAC(2,3)=(-DY(1)-DY(2)-DY(3)-DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJAC(3,3)=(-DZ(1)-DZ(2)-DZ(3)-DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      ELSE
       DJ11=( DX(1)+DX(2)+DX(3)+DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJ12=( DY(1)+DY(2)+DY(3)+DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJ13=( DZ(1)+DZ(2)+DZ(3)+DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DJ21=(-DX(1)+DX(2)+DX(3)-DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
       DJ22=(-DY(1)+DY(2)+DY(3)-DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
       DJ23=(-DZ(1)+DZ(2)+DZ(3)-DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
       DJ31=(-DX(1)-DX(2)+DX(3)+DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
       DJ32=(-DY(1)-DY(2)+DY(3)+DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
       DJ33=(-DZ(1)-DZ(2)+DZ(3)+DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
       DJ41=(-DX(1)-DX(2)-DX(3)-DX(4)+DX(5)+DX(6)+DX(7)+DX(8))*Q8
       DJ42=(-DY(1)-DY(2)-DY(3)-DY(4)+DY(5)+DY(6)+DY(7)+DY(8))*Q8
       DJ43=(-DZ(1)-DZ(2)-DZ(3)-DZ(4)+DZ(5)+DZ(6)+DZ(7)+DZ(8))*Q8
       DJ51=( DX(1)-DX(2)+DX(3)-DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
       DJ52=( DY(1)-DY(2)+DY(3)-DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
       DJ53=( DZ(1)-DZ(2)+DZ(3)-DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
       DJ61=( DX(1)-DX(2)-DX(3)+DX(4)-DX(5)+DX(6)+DX(7)-DX(8))*Q8
       DJ62=( DY(1)-DY(2)-DY(3)+DY(4)-DY(5)+DY(6)+DY(7)-DY(8))*Q8
       DJ63=( DZ(1)-DZ(2)-DZ(3)+DZ(4)-DZ(5)+DZ(6)+DZ(7)-DZ(8))*Q8
       DJ71=( DX(1)+DX(2)-DX(3)-DX(4)-DX(5)-DX(6)+DX(7)+DX(8))*Q8
       DJ72=( DY(1)+DY(2)-DY(3)-DY(4)-DY(5)-DY(6)+DY(7)+DY(8))*Q8
       DJ73=( DZ(1)+DZ(2)-DZ(3)-DZ(4)-DZ(5)-DZ(6)+DZ(7)+DZ(8))*Q8
       DJ81=(-DX(1)+DX(2)-DX(3)+DX(4)+DX(5)-DX(6)+DX(7)-DX(8))*Q8
       DJ82=(-DY(1)+DY(2)-DY(3)+DY(4)+DY(5)-DY(6)+DY(7)-DY(8))*Q8
       DJ83=(-DZ(1)+DZ(2)-DZ(3)+DZ(4)+DZ(5)-DZ(6)+DZ(7)-DZ(8))*Q8
      ENDIF
C
ccc      CALL ZTIME(TTT1)
ccc      TTU2=TTU2+TTT1-TTT0
C *** Dummy call - ELE may save arithmetic operations
ccc      CALL ZTIME(TTT0)
      CALL ELE(0D0,0D0,0D0,-2)
      IF (IER.LT.0) GOTO 99999
ccc      CALL ZTIME(TTT1)
ccc      TTU3=TTU3+TTT1-TTT0
C
C *** Cubature points on the reference element
      DO 200 ICUBP = 1, NCUBP
ccc      CALL ZTIME(TTT0)
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the bilinear mapping onto the reference element
      IF (ILINT.EQ.0) THEN
       DJAC(1,1)=DJ21+DJ51*XI2+DJ61*XI3+DJ81*XI2*XI3
       DJAC(1,2)=DJ31+DJ51*XI1+DJ71*XI3+DJ81*XI1*XI3
       DJAC(1,3)=DJ41+DJ61*XI1+DJ71*XI2+DJ81*XI1*XI2
       DJAC(2,1)=DJ22+DJ52*XI2+DJ62*XI3+DJ82*XI2*XI3
       DJAC(2,2)=DJ32+DJ52*XI1+DJ72*XI3+DJ82*XI1*XI3
       DJAC(2,3)=DJ42+DJ62*XI1+DJ72*XI2+DJ82*XI1*XI2
       DJAC(3,1)=DJ23+DJ53*XI2+DJ63*XI3+DJ83*XI2*XI3
       DJAC(3,2)=DJ33+DJ53*XI1+DJ73*XI3+DJ83*XI1*XI3
       DJAC(3,3)=DJ43+DJ63*XI1+DJ73*XI2+DJ83*XI1*XI2
       DETJ= DJAC(1,1)*(DJAC(2,2)*DJAC(3,3)-DJAC(3,2)*DJAC(2,3))
     *      -DJAC(2,1)*(DJAC(1,2)*DJAC(3,3)-DJAC(3,2)*DJAC(1,3))
     *      +DJAC(3,1)*(DJAC(1,2)*DJAC(2,3)-DJAC(2,2)*DJAC(1,3))
      ENDIF
C
C *** Cubature points on the actual element + weights
      IF (ILINT.EQ.2) THEN
       XX=DJ11+DJAC(1,1)*XI1+DJAC(1,2)*XI2
       YY=DJ12+DJAC(2,1)*XI1+DJAC(2,2)*XI2
       ZZ=DJ13+DJAC(3,3)*XI3
      ELSE IF (ILINT.EQ.1) THEN
       XX=DJ11+DJAC(1,1)*XI1+DJAC(1,2)*XI2+DJAC(1,3)*XI3
       YY=DJ12+DJAC(2,1)*XI1+DJAC(2,2)*XI2+DJAC(2,3)*XI3
       ZZ=DJ13+DJAC(3,1)*XI1+DJAC(3,2)*XI2+DJAC(3,3)*XI3
      ELSE
       XX=DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
       YY=DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
       ZZ=DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2
      ENDIF
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
ccc      CALL ZTIME(TTT1)
ccc      TTU2=TTU2+TTT1-TTT0
ccc      CALL ZTIME(TTT0)
      CALL ELE(XX,YY,ZZ,-3)
      IF (IER.LT.0) GOTO 99999
C
ccc      CALL ZTIME(TTT1)
ccc      TTU4=TTU4+TTT1-TTT0
ccc      CALL ZTIME(TTT0)
C *** Evaluation of velocity in cubature points
      IF (DCMASS.NE.0D0) THEN
C
       DU1=0D0
       DU2=0D0
       DU3=0D0
       IF (A2L.EQ.0D0) THEN
        DO 210 JDFL=1,IDFL
        HBAS=DBAS(1,KDFL(JDFL),1)
        IF (ABS(HBAS).GE.1D-8) THEN
         JDFG=KDFG(JDFL)
         DU1=DU1+U1L1(JDFG)*HBAS
         DU2=DU2+U1L2(JDFG)*HBAS
         DU3=DU3+U1L3(JDFG)*HBAS
        ENDIF
210     CONTINUE
       ELSE
        DO 220 JDFL=1,IDFL
        HBAS=DBAS(1,KDFL(JDFL),1)
        IF (ABS(HBAS).GE.1D-8) THEN
         JDFG=KDFG(JDFL)
         DU1=DU1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
         DU2=DU2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
         DU3=DU3+(A1L*U1L3(JDFG)+A2L*U2L3(JDFG))*HBAS
         ENDIF
220     CONTINUE
       ENDIF
C
C ***  Summing up over all pairs of multiindices
       DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(1,JDOFEH,1)
       HBASJ2=DBAS(1,JDOFEH,2)
       HBASJ3=DBAS(1,JDOFEH,3)
       HBASJ4=DBAS(1,JDOFEH,4)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2+HBASJ4*DU3
C
       DO 240 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
        AH= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *     +  DNY*(HBASJ2**2+HBASJ3**2+HBASJ4**2)+CT0*HBASJ1**2
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(1,IDOFEH,1)
        HBASI2=DBAS(1,IDOFEH,2)
        HBASI3=DBAS(1,IDOFEH,3)
        HBASI4=DBAS(1,IDOFEH,4)
        HSUMI=HBASI2*DU1+HBASI3*DU2+HBASI4*DU3
        AH= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *    +   DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3+HBASJ4*HBASI4)
     *    +   CT0*HBASJ1*HBASI1
       ENDIF
C
       DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
240    CONTINUE
230    CONTINUE
ccc       CALL ZTIME(TTT1)
ccc       TTU12=TTU12+TTT1-TTT0
C
      ELSE
C
C ***  Summing up over 1 pair of multiindices
       DO 250 JDOFE=1,IDFL
       HBASJ1=DBAS(1,KDFL(JDOFE),1)
C
       DO 260 IDOFE=1,IDFL
       HBASI1=DBAS(1,KDFL(IDOFE),1)
C
       AH=-1D0/THSTEP*HBASJ1*HBASI1
C
       DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
260    CONTINUE
250    CONTINUE
ccc       CALL ZTIME(TTT1)
ccc       TTU12=TTU12+TTT1-TTT0
C
      ENDIF
C
200   CONTINUE
C
C
ccc      CALL ZTIME(TTT0)
      DO 400 JDOFE=1,IDFL
      DO 400 IDOFE=1,IDFL
      DENTH=THSTEP*DENTRY(JDOFE,IDOFE)
C
      IF (IDEF.LT.2) THEN
       IA   =KENTRY(JDOFE,IDOFE)
       A(IA)=A(IA)+REAL(DENTH)
      ENDIF
C
      IF (IDEF.GT.0) THEN 
       IDFG=KDFG(IDOFE)
       JDFG=KDFG(JDOFE)
       D1(JDFG)= D1(JDFG)-DENTH*U1(IDFG)
       D2(JDFG)= D2(JDFG)-DENTH*U2(IDFG)
       D3(JDFG)= D3(JDFG)-DENTH*U3(IDFG)
      ENDIF 
C
400   CONTINUE
C
ccc      CALL ZTIME(TTT1)
ccc      TTU13=TTU13+TTT1-TTT0
100   CONTINUE
C
ccc      write(6,*) TTU0,TTU11,TTU12,TTU13,TTU2,TTU3,TTU4,NCUBP
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE  DELTSD  (U1L1,U1L2,U1L3,U2L1,U2L2,U2L3,A1L,A2L,IEL,
     *                     DUMAX,DELTA)
************************************************************************
*     Calculates coefficient Delta for SD
*     RELOC = local RE
*     UNORM=local velocity norm
*     HLOCAL=local mesh width
************************************************************************
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNBAS=27,NNLEV=9,NNARR=299,NNWORK=1)
C
      DIMENSION U1L1(*),U1L2(*),U1L3(*),U2L1(*),U2L2(*),U2L3(*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),
     *                KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),KLB3(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLVOL(NNLEV),LU1OLD,LU2OLD,LU3OLD,
     *                LPOLD,LD1,LD2,LD3,LDP
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC
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
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C
C
C
      AVOL=DBLE(VWORK(L(KLVOL(ILEV))+IEL-1))
      HLOCAL=AVOL**(1D0/3D0)
C
      IF (UPSAM.LT.0D0) THEN
C
       DELTA=ABS(UPSAM)*HLOCAL
C
      ELSE
C
       DU1=0D0
       DU2=0D0
       DU3=0D0
       DO 10 IDOF=1,IDFL
       DU1=DU1+(A1L*U1L1(KDFG(IDOF))+A2L*U2L1(KDFG(IDOF)))
       DU2=DU2+(A1L*U1L2(KDFG(IDOF))+A2L*U2L2(KDFG(IDOF)))
10     DU3=DU3+(A1L*U1L3(KDFG(IDOF))+A2L*U2L3(KDFG(IDOF)))
C
       UNORM=1D0/6D0*SQRT(DU1**2+DU2**2+DU3**2)
C
C       IF (UNORM.LE.1D-8) THEN
C        DELTA=0D0
C       ELSE
        RELOC=UNORM*HLOCAL/NY
        DELTA=UPSAM*HLOCAL/DUMAX*2D0*(RELOC/(1D0+RELOC))
C       ENDIF
C
      ENDIF
C
c      WRITE(6,*) DELTA,DELTA/(UPSAM*HLOCAL),RELOC,UNORM,DUMAX,HLOCAL
C
C
      END
