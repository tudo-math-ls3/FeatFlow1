      SUBROUTINE ELPQU(DU1,DU2,DU3,KVERT,KEDGE,KAREA,DCORVG,ELE,
     *                 ICUB,ILINT,U,UX,UY,UZ,MFILE)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6)
      PARAMETER (NNDIM=3,Q2=0.5D0,Q8=0.125D0)
      CHARACTER FMT*15,SUB*6,CPARAM*120
      DIMENSION DU1(*),DU2(*),DU3(*)
      DIMENSION KVERT(NNVE,*),KEDGE(NNEE,*),KAREA(NNAE,*),DCORVG(3,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
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
      SAVE
C
      SUB='ELPQU'
      IER=0
C
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IF (IER.NE.0) GOTO 99999
      IDFL=NDFL(IELTYP)
      IF (IER.LT.0) GOTO 99999
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
      DO 1 IDER=1,NNDER
1     BDER(IDER)=.FALSE.
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
      BDER(4)=.TRUE.
C
      ERRL2= 0D0
      ERRH1= 0D0
      DIVL2= 0D0
      DIVLI=-1D0
      DNL2=  0D0
      DNH1=  0D0
C
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
      IF (ILINT.EQ.2) THEN
       DJAC(1,3)=0D0
       DJAC(2,3)=0D0
       DJAC(3,1)=0D0
       DJAC(3,2)=0D0
      ENDIF
C
      DO 100 IEL=1,NEL
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
      DO 110 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
110   CONTINUE
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
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
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
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
      IF (ILINT.EQ.2) THEN
       XX= DJ11+DJAC(1,1)*XI1+DJAC(1,2)*XI2
       YY= DJ12+DJAC(2,1)*XI1+DJAC(2,2)*XI2
       ZZ= DJ13+DJAC(3,3)*XI3
      ELSE IF (ILINT.EQ.1) THEN
       XX= DJ11+DJAC(1,1)*XI1+DJAC(1,2)*XI2+DJAC(1,3)*XI3
       YY= DJ12+DJAC(2,1)*XI1+DJAC(2,2)*XI2+DJAC(2,3)*XI3
       ZZ= DJ13+DJAC(3,1)*XI1+DJAC(3,2)*XI2+DJAC(3,3)*XI3
      ELSE
       XX= DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
       YY= DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
       ZZ= DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2
      ENDIF
C
      UH1 =0D0
      UH2 =0D0
      UH3 =0D0
      UH1X=0D0
      UH2X=0D0
      UH3X=0D0
      UH1Y=0D0
      UH2Y=0D0
      UH3Y=0D0
      UH1Z=0D0
      UH2Z=0D0
      UH3Z=0D0
C
      DO 210 JDOFE=1,IDFL
      IEQ=KDFG(JDOFE)
      ILO=KDFL(JDOFE)
      UH1 =UH1 +DU1(IEQ)*DBAS(1,ILO,1)
      UH2 =UH2 +DU2(IEQ)*DBAS(1,ILO,1)
      UH3 =UH3 +DU3(IEQ)*DBAS(1,ILO,1)
      UH1X=UH1X+DU1(IEQ)*DBAS(1,ILO,2)
      UH2X=UH2X+DU2(IEQ)*DBAS(1,ILO,2)
      UH3X=UH3X+DU3(IEQ)*DBAS(1,ILO,2)
      UH1Y=UH1Y+DU1(IEQ)*DBAS(1,ILO,3)
      UH2Y=UH2Y+DU2(IEQ)*DBAS(1,ILO,3)
      UH3Y=UH3Y+DU3(IEQ)*DBAS(1,ILO,3)
      UH1Z=UH1Z+DU1(IEQ)*DBAS(1,ILO,4)
      UH2Z=UH2Z+DU2(IEQ)*DBAS(1,ILO,4)
      UH3Z=UH3Z+DU3(IEQ)*DBAS(1,ILO,4)
210   CONTINUE
C
      ERRL2=ERRL2+OM*(U(XX,YY,ZZ,1)-UH1)**2+OM*(U(XX,YY,ZZ,2)-UH2)**2+
     *            OM*(U(XX,YY,ZZ,3)-UH3)**2
      ERRH1=ERRH1+OM*((UX(XX,YY,ZZ,1)-UH1X)**2+
     *                (UY(XX,YY,ZZ,1)-UH1Y)**2+
     *                (UZ(XX,YY,ZZ,1)-UH1Z)**2)
      ERRH1=ERRH1+OM*((UX(XX,YY,ZZ,2)-UH2X)**2+
     *                (UY(XX,YY,ZZ,2)-UH2Y)**2+
     *                (UZ(XX,YY,ZZ,2)-UH2Z)**2)
      ERRH1=ERRH1+OM*((UX(XX,YY,ZZ,3)-UH3X)**2+
     *                (UY(XX,YY,ZZ,3)-UH3Y)**2+
     *                (UZ(XX,YY,ZZ,3)-UH3Z)**2)
      DIVL2=DIVL2+OM*(UH1X+UH2Y+UH3Z)**2
      DIVLI=MAX(DIVLI,ABS(UH1X+UH2Y+UH3Z))
      DNL2=DNL2+OM*(U(XX,YY,ZZ,1))**2+OM*(U(XX,YY,ZZ,2))**2+
     *          OM*(U(XX,YY,ZZ,3))**2
      DNH1=DNH1+OM*(UX(XX,YY,ZZ,1)**2+UY(XX,YY,ZZ,1)**2+
     *              UZ(XX,YY,ZZ,1)**2)
      DNH1=DNH1+OM*(UX(XX,YY,ZZ,2)**2+UY(XX,YY,ZZ,2)**2+
     *              UZ(XX,YY,ZZ,2)**2)
      DNH1=DNH1+OM*(UX(XX,YY,ZZ,3)**2+UY(XX,YY,ZZ,3)**2+
     *              UZ(XX,YY,ZZ,3)**2)
C
200   CONTINUE
100   CONTINUE
C
      IF (DNL2.LT.1D-15) THEN
       WRITE(MTERM,*)
       WRITE(MTERM,*) '* ELPQU * EXACT SOLUTION ZERO !!!'
       WRITE(MFILE,*)
       WRITE(MFILE,*) '* ELPQU * EXACT SOLUTION ZERO !!!'
       DNL2=1D0
      ENDIF
C
      IF (DNH1.LT.1D-15) THEN
       WRITE(MTERM,*)
       WRITE(MTERM,*) '* ELPQU * EXACT DERIVATIVE ZERO !!!'
       WRITE(MFILE,*)
       WRITE(MFILE,*) '* ELPQU * EXACT DERIVATIVE ZERO !!!'
       DNH1=1D0
      ENDIF
C
c     IF (MT.GT.0) THEN
       WRITE(MTERM,*)
       WRITE(MTERM,*) '*ELPQU*REL. L2-H1-ERROR',
     *                 SQRT(ERRL2/DNL2),SQRT(ERRH1/DNH1),IELTYP,ICUB
       WRITE(MTERM,*) '*ELPQU*DIVERGENZ       ',SQRT(DIVL2),DIVLI
c      ENDIF
C
       IF (MTERM.NE.MFILE) THEN
        WRITE(MFILE,*)
        WRITE(MFILE,*) '*ELPQU*REL. L2-H1-ERROR',
     *                 SQRT(ERRL2/DNL2),SQRT(ERRH1/DNH1),IELTYP,ICUB
        WRITE(MFILE,*) '*ELPQU*DIVERGENZ       ',SQRT(DIVL2),DIVLI
       ENDIF
C
99999 END
C
C=======================================================================
C
      SUBROUTINE ELPQN(DU1,DU2,DU3,KVERT,KEDGE,KAREA,DCORVG,ELE,
     *                 ICUB,ILINT,U,UX,UY,UZ,MFILE)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6,
     *           NNDIM=3,NNCOF=10)
      PARAMETER (Q2=0.5D0,Q8=0.125D0)
      DIMENSION DU1(*),DU2(*),DU3(*)
      DIMENSION KVERT(NNVE,*),KEDGE(NNEE,*),KAREA(NNAE,*),DCORVG(3,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
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
      SAVE
C
      SUB='ELPQN'
      IER=0
C
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      DO 1 IDER=1,NNDER
1     BDER(IDER)=.FALSE.
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
      BDER(4)=.TRUE.
C
      ERRL2= 0D0
      ERRH1= 0D0
      DIVL2= 0D0
      DIVLI=-1D0
      DNL2=  0D0
      DNH1=  0D0
C
C *** Set zero elements of Jacobian for axiparallel grid
      IF (ILINT.EQ.2) THEN
       DJAC(1,3)=0D0
       DJAC(2,3)=0D0
       DJAC(3,1)=0D0
       DJAC(3,2)=0D0
      ENDIF
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Evaluation of coordinates of the vertices
      DO 110 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
110   CONTINUE
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
C *** Dummy call - ELE may save arithmetic operations
      CALL ELE(0D0,0D0,0D0,-2)
      IF (IER.LT.0) GOTO 99999
C
C *** Cubature points on the reference element
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
C *** Jacobian of the trilinear mapping onto the reference element
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
       XX= DJ11+DJAC(1,1)*XI1+DJAC(1,2)*XI2
       YY= DJ12+DJAC(2,1)*XI1+DJAC(2,2)*XI2
       ZZ= DJ13+DJAC(3,3)*XI3
      ELSE IF (ILINT.EQ.1) THEN
       XX= DJ11+DJAC(1,1)*XI1+DJAC(1,2)*XI2+DJAC(1,3)*XI3
       YY= DJ12+DJAC(2,1)*XI1+DJAC(2,2)*XI2+DJAC(2,3)*XI3
       ZZ= DJ13+DJAC(3,1)*XI1+DJAC(3,2)*XI2+DJAC(3,3)*XI3
      ELSE
       XX= DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
       YY= DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
       ZZ= DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2
      ENDIF
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XX,YY,ZZ,-3)
      IF (IER.LT.0) GOTO 99999
C
      UH1 =0D0
      UH2 =0D0
      UH3 =0D0
      UH1X=0D0
      UH2X=0D0
      UH3X=0D0
      UH1Y=0D0
      UH2Y=0D0
      UH3Y=0D0
      UH1Z=0D0
      UH2Z=0D0
      UH3Z=0D0
C
      DO 210 JDOFE=1,IDFL
      IEQ=KDFG(JDOFE)
      ILO=KDFL(JDOFE)
      UH1 =UH1 +DU1(IEQ)*DBAS(1,ILO,1)
      UH2 =UH2 +DU2(IEQ)*DBAS(1,ILO,1)
      UH3 =UH3 +DU3(IEQ)*DBAS(1,ILO,1)
      UH1X=UH1X+DU1(IEQ)*DBAS(1,ILO,2)
      UH2X=UH2X+DU2(IEQ)*DBAS(1,ILO,2)
      UH3X=UH3X+DU3(IEQ)*DBAS(1,ILO,2)
      UH1Y=UH1Y+DU1(IEQ)*DBAS(1,ILO,3)
      UH2Y=UH2Y+DU2(IEQ)*DBAS(1,ILO,3)
      UH3Y=UH3Y+DU3(IEQ)*DBAS(1,ILO,3)
      UH1Z=UH1Z+DU1(IEQ)*DBAS(1,ILO,4)
      UH2Z=UH2Z+DU2(IEQ)*DBAS(1,ILO,4)
      UH3Z=UH3Z+DU3(IEQ)*DBAS(1,ILO,4)
210   CONTINUE
C
      ERRL2=ERRL2+OM*(U(XX,YY,ZZ,1)-UH1)**2+OM*(U(XX,YY,ZZ,2)-UH2)**2+
     *            OM*(U(XX,YY,ZZ,3)-UH3)**2
      ERRH1=ERRH1+OM*((UX(XX,YY,ZZ,1)-UH1X)**2+
     *                (UY(XX,YY,ZZ,1)-UH1Y)**2+
     *                (UZ(XX,YY,ZZ,1)-UH1Z)**2)
      ERRH1=ERRH1+OM*((UX(XX,YY,ZZ,2)-UH2X)**2+
     *                (UY(XX,YY,ZZ,2)-UH2Y)**2+
     *                (UZ(XX,YY,ZZ,2)-UH2Z)**2)
      ERRH1=ERRH1+OM*((UX(XX,YY,ZZ,3)-UH3X)**2+
     *                (UY(XX,YY,ZZ,3)-UH3Y)**2+
     *                (UZ(XX,YY,ZZ,3)-UH3Z)**2)
      DIVL2=DIVL2+OM*(UH1X+UH2Y+UH3Z)**2
      DIVLI=MAX(DIVLI,ABS(UH1X+UH2Y+UH3Z))
      DNL2=DNL2+OM*(U(XX,YY,ZZ,1))**2+OM*(U(XX,YY,ZZ,2))**2+
     *          OM*(U(XX,YY,ZZ,3))**2
      DNH1=DNH1+OM*(UX(XX,YY,ZZ,1)**2+UY(XX,YY,ZZ,1)**2+
     *              UZ(XX,YY,ZZ,1)**2)
      DNH1=DNH1+OM*(UX(XX,YY,ZZ,2)**2+UY(XX,YY,ZZ,2)**2+
     *              UZ(XX,YY,ZZ,2)**2)
      DNH1=DNH1+OM*(UX(XX,YY,ZZ,3)**2+UY(XX,YY,ZZ,3)**2+
     *              UZ(XX,YY,ZZ,3)**2)
C
200   CONTINUE
100   CONTINUE
C
      IF (DNL2.LT.1D-15) THEN
       WRITE(MTERM,*)
       WRITE(MTERM,*) '* ELPQU * EXACT SOLUTION ZERO !!!'
       WRITE(MFILE,*)
       WRITE(MFILE,*) '* ELPQU * EXACT SOLUTION ZERO !!!'
       DNL2=1D0
      ENDIF
C
      IF (DNH1.LT.1D-15) THEN
       WRITE(MTERM,*)
       WRITE(MTERM,*) '* ELPQU * EXACT DERIVATIVE ZERO !!!'
       WRITE(MFILE,*)
       WRITE(MFILE,*) '* ELPQU * EXACT DERIVATIVE ZERO !!!'
       DNH1=1D0
      ENDIF
C
c     IF (MT.GT.0) THEN
       WRITE(MTERM,*)
       WRITE(MTERM,*) '*ELPQU*REL. L2-H1-ERROR',
     *                 SQRT(ERRL2/DNL2),SQRT(ERRH1/DNH1),IELTYP,ICUB
       WRITE(MTERM,*) '*ELPQU*DIVERGENZ       ',SQRT(DIVL2),DIVLI
c      ENDIF
C
       IF (MTERM.NE.MFILE) THEN
        WRITE(MFILE,*)
        WRITE(MFILE,*) '*ELPQU*REL. L2-H1-ERROR',
     *                 SQRT(ERRL2/DNL2),SQRT(ERRH1/DNH1),IELTYP,ICUB
        WRITE(MFILE,*) '*ELPQU*DIVERGENZ       ',SQRT(DIVL2),DIVLI
       ENDIF
C
99999 END
C
C=======================================================================
C
      SUBROUTINE ELPQP(P,PERR,KVERT,KEDGE,KAREA,DCORVG,ELE,ICUB,EP,
     *                 MFILE,ILINT)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNBAS=27,NNDER=10,NNCUBP=36,NNVE=8,NNEE=12,NNAE=6)
      PARAMETER (NNDIM=3,Q2=0.5D0,Q8=0.125D0)
      CHARACTER FMT*15,SUB*6,CPARAM*120
      DIMENSION P(*),PERR(*)
      DIMENSION KVERT(NNVE,*),KEDGE(NNEE,*),KAREA(NNAE,*),DCORVG(3,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
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
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/ELEM/,/TRIAD/,/CUB/,/COAUX1/
C
      SUB='ELPQP'
      IER=0
C
      IELTYP=-1
      CALL ELE(0D0,0D0,0D0,IELTYP)
      IF (IER.NE.0) GOTO 99999
      IDFL=NDFL(IELTYP)
      IF (IER.LT.0) GOTO 99999
      CALL CB3H(ICUB)
      IF (IER.NE.0) GOTO 99999
      DO 1 IDER=1,NNDER
1     BDER(IDER)=.FALSE.
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
      BDER(4)=.TRUE.
C
      ERRL2=0D0
      DNL2=0D0
C
      ICUBP=ICUB
      CALL ELE(0D0,0D0,0D0,-2)
      IF (ILINT.EQ.2) THEN
       DJAC(1,3)=0D0
       DJAC(2,3)=0D0
       DJAC(3,1)=0D0
       DJAC(3,2)=0D0
      ENDIF
C
      DO 100 IEL=1,NEL
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KEDGE,KAREA,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
      DO 110 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
      DZ(IVE)=DCORVG(3,JP)
110   CONTINUE
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
      DO 200 ICUBP=1,NCUBP
C
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
      XI3=DXI(ICUBP,3)
C
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
C      DETJ=ABS(DETJ)
      OM=DOMEGA(ICUBP)*ABS(DETJ)
C
      CALL ELE(XI1,XI2,XI3,-3)
      IF (IER.LT.0) GOTO 99999
C
      IF (ILINT.EQ.2) THEN
       XX= DJ11+DJAC(1,1)*XI1+DJAC(1,2)*XI2
       YY= DJ12+DJAC(2,1)*XI1+DJAC(2,2)*XI2
       ZZ= DJ13+DJAC(3,3)*XI3
      ELSE IF (ILINT.EQ.1) THEN
       XX= DJ11+DJAC(1,1)*XI1+DJAC(1,2)*XI2+DJAC(1,3)*XI3
       YY= DJ12+DJAC(2,1)*XI1+DJAC(2,2)*XI2+DJAC(2,3)*XI3
       ZZ= DJ13+DJAC(3,1)*XI1+DJAC(3,2)*XI2+DJAC(3,3)*XI3
      ELSE
       XX= DJ11+DJAC(1,1)*XI1+DJ31*XI2+DJ41*XI3+DJ71*XI2*XI3
       YY= DJ12+DJ22*XI1+DJAC(2,2)*XI2+DJ42*XI3+DJ62*XI1*XI3
       ZZ= DJ13+DJ23*XI1+DJ33*XI2+DJAC(3,3)*XI3+DJ53*XI1*XI2
      ENDIF
C
      PH =0D0
C
      DO 210 JDOFE=1,IDFL
      IEQ=KDFG(JDOFE)
      ILO=KDFL(JDOFE)
      PH =PH +P(IEQ)*DBAS(1,ILO,1)
210   CONTINUE
C
      ERRL2=ERRL2+OM*(EP(XX,YY,ZZ)-PH)**2
      DNL2=DNL2+OM*EP(XX,YY,ZZ)**2
      PERR(IEQ)=ABS(EP(XX,YY,ZZ)-PH)
C
200   CONTINUE
100   CONTINUE
C
      IF (DNL2.LT.1D-15) THEN
       WRITE(MTERM,*)
       WRITE(MTERM,*) '* ELPQP * EXACT SOLUTION ZERO !!!'
       WRITE(MFILE,*)
       WRITE(MFILE,*) '* ELPQP * EXACT SOLUTION ZERO !!!'
       DNL2=1D0
       DNH1=1D0
      ENDIF

c     IF (MT.GT.0) THEN
       WRITE(MTERM,*) '*ELPQP*REL. L2-ERROR   '
     *                 ,SQRT(ERRL2/DNL2),IELTYP,ICUB
       WRITE(MTERM,*)
c      ENDIF
C
       IF (MTERM.NE.MFILE) THEN
        WRITE(MFILE,*) '*ELPQP*REL. L2-ERROR   '
     *                  ,SQRT(ERRL2/DNL2),IELTYP,ICUB
        WRITE(MFILE,*)
       ENDIF
C
99999 END
