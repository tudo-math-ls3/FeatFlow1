************************************************************************
      SUBROUTINE  MATML (D1,D2,U1,U2,C1,C2,A,KCOLA,KLDA,NU)
************************************************************************
*
*   Purpose:  performs the matrix-vector-operation
* 
*                   D:= C1*(A*U) + C2*D
*
*             for the vectors  D=(D1,D2)  and  U=(U1,U2)  with the 
*             given scalar variables C1,C2
*  
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL A
C
      PARAMETER (NNVE=4)
      DIMENSION D1(*),D2(*),U1(*),U2(*)
      DIMENSION A(*),KCOLA(*),KLDA(*)
      
      SAVE
C
      NA=KLDA(NU+1)-1
C
C-----------------------------------------------------------------------
C     Calculation of  D1=C2*D1+C1*(A1*U1+A2*U2)
C-----------------------------------------------------------------------
      CALL  LAX37 (A(1),KCOLA,KLDA,NU,U1,D1,C1,C2)
      CALL  LAX37 (A(NA+1),KCOLA,KLDA,NU,U2,D1,C1,1d0)
C
C-----------------------------------------------------------------------
C     Calculation of  D2=C2*D2+C1*(A3*U1+A4*U2)
C-----------------------------------------------------------------------
      CALL  LAX37 (A(2*NA+1),KCOLA,KLDA,NU,U1,D2,C1,C2)
      CALL  LAX37 (A(3*NA+1),KCOLA,KLDA,NU,U2,D2,C1,1d0)
C
C *** Set the Dirichlet-components of (D1,D2) to zero
c      CALL  BDRY0 (D1,D2,KMBD,NMBD)
C
      END
      



