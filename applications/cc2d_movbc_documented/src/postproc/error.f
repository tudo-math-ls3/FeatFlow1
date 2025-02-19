************************************************************************
* The routines in this file realize error analyzing routines
* for velocity and pressure components.
************************************************************************

************************************************************************
* Error analysis for velocity components
*
* Parametric version; for element E030,E031
*
* Computes L2-error, H1-error, divergence,... and prints the results
* onto screen.
*
* In:
*   DU1    - array [1..*] of double
*            Velocity X-component
*   DU2    - array [1..*] of double
*            Velocity Y-component
*   KVERT,
*   KMID,
*   DCORVG - Usual geometry information
*   ELE    - SUBROUTINE; Element routine
*   ICUB   - Cubature formula for error analysis
*   U      - SUBROUTINE; Reference solution
*   UX     - SUBROUTINE; Reference X-derivative
*   UY     - SUBROUTINE; Reference Y-derivative
*   MFILE  - Handle of file where to write the output to, additionally
*            to the screen
************************************************************************

      SUBROUTINE ELPQU(DU1,DU2,KVERT,KMID,DCORVG,ELE,ICUB,U,UX,UY,MFILE)

      IMPLICIT NONE

C main COMMON blocks

      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'

      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'ccub.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      INCLUDE 'casmbly.inc'

C parameters
      
      DOUBLE PRECISION DU1(*), DU2(*),DCORVG(2,*),U,UX,UY
      INTEGER KVERT(NNVE,*),KMID(NNVE,*)
      INTEGER ICUB,MFILE

      EXTERNAL ELE
      
C externals

      INTEGER NDFL
      EXTERNAL NDFL
      
C local variables      

      DOUBLE PRECISION ERRL2,ERRH1,DIVL2,DIVLI,DNL2,DNH1
      INTEGER IVE,JP
      DOUBLE PRECISION XX,YY,OM
      DOUBLE PRECISION DJ1,DJ2,DJ3,DJ4,XI1,XI2
      DOUBLE PRECISION UH1,UH2,UH1X,UH2X,UH1Y,UH2Y
      INTEGER IELTYP,IDER,JDOFE,IEQ,ILO
      
      
      SUB='ELPQU'
      IER=0

C     Ask the element about its type:

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IF (IER.NE.0) GOTO 99999
      
C     Get local degrees of freedom
      
      IDFL=NDFL(IELTYP)
      IF (IER.LT.0) GOTO 99999
      
C     Initialize cubature formula
      
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
      
C     Initialize derivative flags; we want to have
C     function value and 1st derivatives:
      
      DO IDER=1,NNDER
        BDER(IDER)=.FALSE.
      END DO
      
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.

      ERRL2= 0D0
      ERRH1= 0D0
      DIVL2= 0D0
      DIVLI=-1D0
      DNL2=  0D0
      DNH1=  0D0
      
C     Call the element top precalculate information in the cubature 
C     points:
      
      CALL ELE(0D0,0D0,-2)

C     Loop over the elements to calculate the information:

      DO IEL=1,NEL

C       Get the global DOF's on our element

        CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999

C       Initialize the element parameters to match to the current element

        DO IVE=1,NVE
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
        END DO

C       Initialize auxiliary Jacobian factors for the transformaion

        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))

C       Loop over the cubature points to calculate the values there:

        DO ICUBP=1,NCUBP

C         Coordinates of the cubature point on the reference element

          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)
      
C         Calculate Jacobian of the mapping and Jacobian determinant
            	
          DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
          DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
          DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
          DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
          DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)

C         And the weight in the cuibature formula
	
          OM=DOMEGA(ICUBP)*DETJ

C         Get the value of the element in the cubature point:

          CALL ELE(XI1,XI2,-3)
          IF (IER.LT.0) GOTO 99999

C         Map the cubature point to the real element
  
          XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *      +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
          YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1+0.5D0*
     *      (DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2   

C         Calculate value, X-derivative and Y-derivative of
C         DU1 and DU2 in our cubature point by a loop over the
C         DOF's in our element

          UH1 =0D0
          UH2 =0D0
          UH1X=0D0
          UH2X=0D0
          UH1Y=0D0
          UH2Y=0D0
          
          DO JDOFE=1,IDFL
            IEQ=KDFG(JDOFE)
            ILO=KDFL(JDOFE)
            UH1 =UH1 +DU1(IEQ)*DBAS(ILO,1)
            UH2 =UH2 +DU2(IEQ)*DBAS(ILO,1)
            UH1X=UH1X+DU1(IEQ)*DBAS(ILO,2)
            UH2X=UH2X+DU2(IEQ)*DBAS(ILO,2)
            UH1Y=UH1Y+DU1(IEQ)*DBAS(ILO,3)
            UH2Y=UH2Y+DU2(IEQ)*DBAS(ILO,3)
          END DO

C         Subtract calculated value from real solution.
C         Compute L2-error,

          ERRL2=ERRL2+OM*((U (XX,YY,1)-UH1 )**2+(U (XX,YY,2)-UH2 )**2)
          
C         H1-error,
          
          ERRH1=ERRH1+OM*((UX(XX,YY,1)-UH1X)**2+(UY(XX,YY,1)-UH1Y)**2)
          ERRH1=ERRH1+OM*((UX(XX,YY,2)-UH2X)**2+(UY(XX,YY,2)-UH2Y)**2)
          
C         divergence in L2-norm
          
          DIVL2=DIVL2+OM*(UH1X+UH2Y)**2
          
C         divergence in SUP-norm,          
          
          DIVLI=MAX(DIVLI,ABS(UH1X+UH2Y))
          
C         as well as L2-norm and H1-norm of the reference solution.

          DNL2=DNL2+OM*(U (XX,YY,1)**2+U (XX,YY,2)**2)
          DNH1=DNH1+OM*(UX(XX,YY,1)**2+UY(XX,YY,1)**2)
          DNH1=DNH1+OM*(UX(XX,YY,2)**2+UY(XX,YY,2)**2)

        END DO ! ICUBP
        
      END DO ! IEL

C     Print the results onto screen:

      IF (DNL2.LT.1D-15) THEN
        WRITE(MTERM,*)
        WRITE(MTERM,*) '* ELPQU * EXACT SOLUTION ZERO !!!'
        WRITE(MFILE,*)
        WRITE(MFILE,*) '* ELPQU * EXACT SOLUTION ZERO !!!'
        DNL2=1D0
      ENDIF

      IF (DNH1.LT.1D-15) THEN
        WRITE(MTERM,*)
        WRITE(MTERM,*) '* ELPQU * EXACT DERIVATIVE ZERO !!!'
        WRITE(MFILE,*)
        WRITE(MFILE,*) '* ELPQU * EXACT DERIVATIVE ZERO !!!'
        DNH1=1D0
      ENDIF

      WRITE(MTERM,*)
      WRITE(MTERM,*) '*ELPQU*REL. L2-H1-ERROR',
     *                 SQRT(ERRL2/DNL2),SQRT(ERRH1/DNH1),IELTYP,ICUB
      WRITE(MTERM,*) '*ELPQU*DIVERGENZ       ',SQRT(DIVL2),DIVLI

      IF (MTERM.NE.MFILE) THEN
        WRITE(MFILE,*)
        WRITE(MFILE,*) '*ELPQU*REL. L2-H1-ERROR',
     *                 SQRT(ERRL2/DNL2),SQRT(ERRH1/DNH1),IELTYP,ICUB
        WRITE(MFILE,*) '*ELPQU*DIVERGENZ       ',SQRT(DIVL2),DIVLI
      ENDIF

99999 END

************************************************************************
* Error analysis for pressure components
*
* Parametric version; for E010
*
* Computes L2-error, and prints the results onto screen.
*
* In:
*   P      - array [1..*] of double
*            Pressure solution
*   KVERT,
*   KMID,
*   DCORVG - Usual geometry information
*   ELE    - SUBROUTINE; Element routine
*   ICUB   - Cubature formula for error analysis
*   EP     - SUBROUTINE; Reference solution
*   MFILE  - Handle of file where to write the output to, additionally
*            to the screen
*
* Out:
*   PERR   - array [1..*] of double
*     This vector corresponds to the pressure vector P. For each DOF
*     in P, the corresponding entry in this vector receives the
*     absolute error to the reference pressure in this DOF.
*     
*     Remark: For the moment, this works only correctly if E010 is used.
************************************************************************

      SUBROUTINE ELPQP(P,PERR,KVERT,KMID,DCORVG,ELE,ICUB,EP,MFILE)

      IMPLICIT NONE

C main COMMON blocks

      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'

      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'ccub.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      INCLUDE 'casmbly.inc'

C parameters
      
      DOUBLE PRECISION P(*),PERR(*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*)
      DOUBLE PRECISION DCORVG(2,*)
      INTEGER ICUB,MFILE
      DOUBLE PRECISION EP

      EXTERNAL ELE
      
C externals

      INTEGER NDFL
      EXTERNAL NDFL

C local variables      

      DOUBLE PRECISION ERRL2,DNL2,DNH1
      INTEGER IVE,JP
      DOUBLE PRECISION XX,YY,OM,PH
      DOUBLE PRECISION DJ1,DJ2,DJ3,DJ4,XI1,XI2
      INTEGER IELTYP,IDER,JDOFE,IEQ,ILO

      SUB='ELPQP'
      IER=0

C     Ask the element about its type:

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IF (IER.NE.0) GOTO 99999
      
C     Get local degrees of freedom
      
      IDFL=NDFL(IELTYP)
      IF (IER.LT.0) GOTO 99999
      
C     Initialize cubature formula
      
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
      
C     Initialize derivative flags; we only need function values.
      
      DO IDER=1,NNDER
        BDER(IDER)=.FALSE.
      END DO
      
      BDER(1)=.TRUE.

      ERRL2=0D0
      DNL2=0D0
      
C     Call the element top precalculate information in the cubature 
C     points:

      CALL ELE(0D0,0D0,-2)

C     Loop over the elements to calculate the information:

      DO IEL=1,NEL

C       Get the global DOF's on our element

        CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999

C       Initialize the element parameters to match to the current element

        DO IVE=1,NVE
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
        END DO

C       Initialize auxiliary Jacobian factors for the transformaion

        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))

C       Loop over the cubature points to calculate the values there:

        DO ICUBP=1,NCUBP

C         Coordinates of the cubature point on the reference element

          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)
      
C         Calculate Jacobian of the mapping and Jacobian determinant
            	
          DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
          DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
          DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
          DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
          DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)

C         And the weight in the cuibature formula
	
          OM=DOMEGA(ICUBP)*DETJ

C         Get the value of the element in the cubature point:

          CALL ELE(XI1,XI2,-3)
          IF (IER.LT.0) GOTO 99999

C         Map the cubature point to the real element
  
          XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *      +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
          YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1+0.5D0*
     *      (DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2   

C         Calculate the pressure in the cubature point:

          PH=0D0
          DO JDOFE=1,IDFL
            IEQ=KDFG(JDOFE)
            ILO=KDFL(JDOFE)
            PH =PH +P(IEQ)*DBAS(ILO,1)
          END DO

C         Subtract calculated value from real solution.
C         Compute L2-error,

          ERRL2=ERRL2+OM*(EP(XX,YY)-PH)**2
          
C         as well as L2-norm of the reference solution.

          DNL2 =DNL2 +OM* EP(XX,YY)**2
          
C         Calculate the pointwise error of each pressure DOF:
          
          PERR(IEQ)=ABS(EP(XX,YY)-PH)
        
        END DO  
        
      END DO
 
C     Write the results to screen
          
      IF (DNL2.LT.1D-15) THEN
        WRITE(MTERM,*)
        WRITE(MTERM,*) '* ELPQP * EXACT SOLUTION ZERO !!!'
        WRITE(MFILE,*)
        WRITE(MFILE,*) '* ELPQP * EXACT SOLUTION ZERO !!!'
        DNL2=1D0
        DNH1=1D0
      ENDIF

      WRITE(MTERM,*) '*ELPQP*REL. L2-ERROR   '
     *                 ,SQRT(ERRL2/DNL2),IELTYP,ICUB
      WRITE(MTERM,*)

      IF (MTERM.NE.MFILE) THEN
        WRITE(MFILE,*) '*ELPQP*REL. L2-ERROR   '
     *                  ,SQRT(ERRL2/DNL2),IELTYP,ICUB
        WRITE(MFILE,*)
      ENDIF

99999 END

************************************************************************
* Error analysis for velocity components
*
* Nonparametric version; for element EM30,EM31
*
* Computes L2-error, H1-error, divergence,... and prints the results
* onto screen.
*
* In:
*   DU1    - array [1..NMT] of double
*            Velocity X-component
*   DU2    - array [1..NMT] of double
*            Velocity Y-component
*   KVERT,
*   KMID,
*   DCORVG - Usual geometry information
*   ELE    - SUBROUTINE; Element routine
*   ICUB   - Cubature formula for error analysis
*   U      - SUBROUTINE; Reference solution
*   UX     - SUBROUTINE; Reference X-derivative
*   UY     - SUBROUTINE; Reference Y-derivative
*   MFILE  - Handle of file where to write the output to, additionally
*            to the screen
************************************************************************

      SUBROUTINE ELPQN(DU1,DU2,KVERT,KMID,DCORVG,ELE,ICUB,U,UX,UY,MFILE)

      IMPLICIT NONE

C main COMMON blocks

      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'

      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'ccub.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      INCLUDE 'casmbly.inc'

C parameters
      
      DOUBLE PRECISION DU1(*), DU2(*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*)
      DOUBLE PRECISION DCORVG(2,*)
      INTEGER ICUB,MFILE
      DOUBLE PRECISION U,UX,UY

      EXTERNAL ELE
      
C externals

      INTEGER NDFL
      EXTERNAL NDFL

C local variables      

      DOUBLE PRECISION ERRL2,ERRH1,DIVL2,DIVLI,DNL2,DNH1
      DOUBLE PRECISION UH1,UH2,UH1X,UH2X,UH1Y,UH2Y
      INTEGER IVE,JP
      DOUBLE PRECISION XX,YY,OM
      DOUBLE PRECISION DJ1,DJ2,DJ3,DJ4,XI1,XI2
      INTEGER IELTYP,IDER,JDOFE,IEQ,ILO

      SUB='ELPQN'
      IER=0
C     Ask the element about its type:

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IF (IER.NE.0) GOTO 99999
      
C     Get local degrees of freedom
      
      IDFL=NDFL(IELTYP)
      IF (IER.LT.0) GOTO 99999
      
C     Initialize cubature formula
      
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
      
C     Initialize derivative flags; we want to have
C     function value and 1st derivatives:
      
      DO IDER=1,NNDER
        BDER(IDER)=.FALSE.
      END DO
      
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.

      ERRL2= 0D0
      ERRH1= 0D0
      DIVL2= 0D0
      DIVLI=-1D0
      DNL2=  0D0
      DNH1=  0D0
      
C     Loop over the elements to calculate the information:

      DO IEL=1,NEL

C       Get the global DOF's on our element

        CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999

C       Initialize the element parameters to match to the current element

        DO IVE=1,NVE
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
        END DO

C       Initialize auxiliary Jacobian factors for the transformaion

        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))

C       Call the element to precalculate values in the cubature
C       points of our element:
      
        CALL ELE(0D0,0D0,-2)

C       Loop over the cubature points to calculate the values there:

        DO ICUBP=1,NCUBP

C         Coordinates of the cubature point on the reference element

          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)
      
C         Calculate Jacobian of the mapping and Jacobian determinant
            	
          DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
          DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
          DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
          DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
          DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)

C         And the weight in the cuibature formula
	
          OM=DOMEGA(ICUBP)*DETJ

C         Map the cubature point to the real element
  
          XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *      +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
          YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1+0.5D0*
     *      (DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2   

C         Get the value of the element in the cubature point:

          CALL ELE(XX,YY,-3)
          IF (IER.LT.0) GOTO 99999

C         Calculate value, X-derivative and Y-derivative of
C         DU1 and DU2 in our cubature point by a loop over the
C         DOF's in our element

          UH1 =0D0
          UH2 =0D0
          UH1X=0D0
          UH2X=0D0
          UH1Y=0D0
          UH2Y=0D0
          
          DO JDOFE=1,IDFL
            IEQ=KDFG(JDOFE)
            ILO=KDFL(JDOFE)
            UH1 =UH1 +DU1(IEQ)*DBAS(ILO,1)
            UH2 =UH2 +DU2(IEQ)*DBAS(ILO,1)
            UH1X=UH1X+DU1(IEQ)*DBAS(ILO,2)
            UH2X=UH2X+DU2(IEQ)*DBAS(ILO,2)
            UH1Y=UH1Y+DU1(IEQ)*DBAS(ILO,3)
            UH2Y=UH2Y+DU2(IEQ)*DBAS(ILO,3)
          END DO

C         Subtract calculated value from real solution.
C         Compute L2-error,

          ERRL2=ERRL2+OM*((U (XX,YY,1)-UH1 )**2+(U (XX,YY,2)-UH2 )**2)
          
C         H1-error,
          
          ERRH1=ERRH1+OM*((UX(XX,YY,1)-UH1X)**2+(UY(XX,YY,1)-UH1Y)**2)
          ERRH1=ERRH1+OM*((UX(XX,YY,2)-UH2X)**2+(UY(XX,YY,2)-UH2Y)**2)
          
C         divergence in L2-norm
          
          DIVL2=DIVL2+OM*(UH1X+UH2Y)**2
          
C         divergence in SUP-norm,          
          
          DIVLI=MAX(DIVLI,ABS(UH1X+UH2Y))
          
C         as well as L2-norm and H1-norm of the reference solution.

          DNL2=DNL2+OM*(U (XX,YY,1)**2+U (XX,YY,2)**2)
          DNH1=DNH1+OM*(UX(XX,YY,1)**2+UY(XX,YY,1)**2)
          DNH1=DNH1+OM*(UX(XX,YY,2)**2+UY(XX,YY,2)**2)

        END DO ! ICUBP
        
      END DO ! IEL

C     Print the results onto screen:

      IF (DNL2.LT.1D-15) THEN
        WRITE(MTERM,*)
        WRITE(MTERM,*) '* ELPQU * EXACT SOLUTION ZERO !!!'
        WRITE(MFILE,*)
        WRITE(MFILE,*) '* ELPQU * EXACT SOLUTION ZERO !!!'
        DNL2=1D0
      ENDIF

      IF (DNH1.LT.1D-15) THEN
        WRITE(MTERM,*)
        WRITE(MTERM,*) '* ELPQU * EXACT DERIVATIVE ZERO !!!'
        WRITE(MFILE,*)
        WRITE(MFILE,*) '* ELPQU * EXACT DERIVATIVE ZERO !!!'
        DNH1=1D0
      ENDIF

      WRITE(MTERM,*)
      WRITE(MTERM,*) '*ELPQU*REL. L2-H1-ERROR',
     *                 SQRT(ERRL2/DNL2),SQRT(ERRH1/DNH1),IELTYP,ICUB
      WRITE(MTERM,*) '*ELPQU*DIVERGENZ       ',SQRT(DIVL2),DIVLI

      IF (MTERM.NE.MFILE) THEN
        WRITE(MFILE,*)
        WRITE(MFILE,*) '*ELPQU*REL. L2-H1-ERROR',
     *                 SQRT(ERRL2/DNL2),SQRT(ERRH1/DNH1),IELTYP,ICUB
        WRITE(MFILE,*) '*ELPQU*DIVERGENZ       ',SQRT(DIVL2),DIVLI
      ENDIF

99999 END

