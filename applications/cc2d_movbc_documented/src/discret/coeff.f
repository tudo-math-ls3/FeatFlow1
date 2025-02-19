************************************************************************
      DOUBLE PRECISION FUNCTION   COEFST (X,Y,IA,IB,IBLOC,BFIRST)
*
*     Coefficient for the Stokes-block
************************************************************************
      
      IMPLICIT NONE

C main COMMON blocks
      
      INCLUDE 'cinidat.inc'
      INCLUDE 'cnsparfrac.inc'

C parameters

      INTEGER IA,IB,IBLOC
      DOUBLE PRECISION X,Y
      LOGICAL BFIRST
      
C local variables

      IF ((IA.EQ.1).AND.(IB.EQ.1)) THEN
       COEFST=1D0
      ELSE
       COEFST=NY
      ENDIF
C
      END
C
************************************************************************
      DOUBLE PRECISION FUNCTION COEFFN(IA,IB,U1L1,U1L2,U2L1,U2L2,
     *                                 A1L,A2L,DELTA,DCMASS)
*
*     Coefficient for the convective-block
************************************************************************
      IMPLICIT NONE

C main COMMON blocks
      
      INCLUDE 'cinidat.inc'
      INCLUDE 'cnsparfrac.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      INCLUDE 'casmbly.inc'

C parameters

      DOUBLE PRECISION U1L1(*),U1L2(*),U2L1(*),U2L2(*)
      DOUBLE PRECISION A1L,A2L,DELTA,DCMASS
      INTEGER IA, IB

C local variables

      DOUBLE PRECISION DU1, DU2, HBAS
      INTEGER JDFL

      DU1=0D0
      DU2=0D0
C 
      IF (DCMASS.NE.0D0) THEN
C
       DO 10 JDFL=1,IDFL
       HBAS=DBAS(KDFL(JDFL),1)
       DU1=DU1+(A1L*U1L1(KDFG(JDFL))+A2L*U2L1(KDFG(JDFL)))*HBAS
       DU2=DU2+(A1L*U1L2(KDFG(JDFL))+A2L*U2L2(KDFG(JDFL)))*HBAS
10     CONTINUE
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
C
      END
C
************************************************************************
      DOUBLE PRECISION FUNCTION   COEFFB (X,Y,IA,IB,IBLOC,BFIRST)
*
*     Coefficient for the B1/B2-blocks
************************************************************************
      
      IMPLICIT NONE

C main COMMON blocks
      
      INCLUDE 'cinidat.inc'
      INCLUDE 'cnsparfrac.inc'

C parameters

      DOUBLE PRECISION X,Y
      INTEGER IA, IB, IBLOC
      LOGICAL BFIRST

      COEFFB= -1.D0
C
      END
C
************************************************************************
      DOUBLE PRECISION FUNCTION    RHS  (X,Y,IA,IBLOC,BFIRST)
*
*     Right hand side yielding the exact solution UE
************************************************************************
      
      IMPLICIT NONE

C main COMMON blocks
      
      INCLUDE 'cinidat.inc'
      INCLUDE 'cnsparfrac.inc'

C parameters

      DOUBLE PRECISION X,Y
      INTEGER IA, IB, IBLOC
      LOGICAL BFIRST
      
C externals
      DOUBLE PRECISION FDATIN
      EXTERNAL FDATIN
      
C local variables

      INTEGER TRIA
      
C
      RHS=FDATIN(5,IBLOC,X,Y,TIMENS,RE)
C
      END
C
C
************************************************************************
* Exact solution - velocity,  also for boundary conditions
*
* This routine determines - depending on a point - the exact
* solution value in that point.
*
* In:
*  X,Y    - coordinates of the point
*  IBLOC  - matrix block that specifies the type of the solution to
*           return.
*           = 0: matrix block for x-velocity
*           = 1: matrix block for y-velocity
*
************************************************************************
      DOUBLE PRECISION FUNCTION UE  (X,Y,IBLOC)
      
      IMPLICIT NONE

C main COMMON blocks
      
      INCLUDE 'cinidat.inc'
      INCLUDE 'cnsparfrac.inc'

C parameters

      DOUBLE PRECISION X,Y
      INTEGER IA, IB, IBLOC
      LOGICAL BFIRST
      
C externals
      DOUBLE PRECISION FDATIN
      EXTERNAL FDATIN

      UE=FDATIN(1,IBLOC,X,Y,TIMENS,RE)
C
      END
C
************************************************************************
*     Exact solution - pressure, only for error analysis
************************************************************************
      DOUBLE PRECISION FUNCTION    PE  (X,Y, IBLOC)
      
      IMPLICIT NONE

C main COMMON blocks
      
      INCLUDE 'cinidat.inc'
      INCLUDE 'cnsparfrac.inc'

C parameters

      DOUBLE PRECISION X,Y
      INTEGER IBLOC

C externals
      DOUBLE PRECISION FDATIN
      EXTERNAL FDATIN
      
C local variables    
      INTEGER TRIA  

      PE=FDATIN(4,IBLOC,X,Y,TIMENS,RE)
C
      END
C
*************************************************************************
      DOUBLE PRECISION FUNCTION    UEX(X,Y,IBLOC)
C
C     x-derivative of exact solution, only for error analysis
*************************************************************************
      
      IMPLICIT NONE

C main COMMON blocks
      
      INCLUDE 'cinidat.inc'
      INCLUDE 'cnsparfrac.inc'

C parameters

      DOUBLE PRECISION X,Y
      INTEGER IBLOC

C externals
      DOUBLE PRECISION FDATIN
      EXTERNAL FDATIN
      
C local variables
      INTEGER TRIA

      UEX=FDATIN(2,IBLOC,X,Y,TIMENS,RE)
C
      END
C
*************************************************************************
      DOUBLE PRECISION FUNCTION    UEY(X,Y,IBLOC)
C
C     y-derivative of exact solution,, only for error analysis
*************************************************************************
      
      IMPLICIT NONE

C main COMMON blocks
      
      INCLUDE 'cinidat.inc'
      INCLUDE 'cnsparfrac.inc'

C parameters

      DOUBLE PRECISION X,Y
      INTEGER IBLOC

C externals
      DOUBLE PRECISION FDATIN
      EXTERNAL FDATIN
      
C local variables
      INTEGER TRIA

      UEY=FDATIN(3,IBLOC,X,Y,TIMENS,RE)   
C
      END
