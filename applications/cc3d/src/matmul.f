************************************************************************
      SUBROUTINE  MATML (D1,D2,D3,DP,U1,U2,U3,P,C1,C2,A,KCOLA,KLDA,
     *                    B1,B2,B3,KCOLB,KLDB,NU,NP,KABD,NABD,INEUM)
************************************************************************
*
*   Purpose:  performs the matrix-vector-operation
* 
*                   D:= C1*(A*U) + C2*D
*
*           for the vectors  D=(D1,D2,D3,DP) and U=(U1,U2,U3,P) with 
*           the given scalar variables C1,C2
*  
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL A,B1,B2,B3
C
      PARAMETER (NNVE=8)
      DIMENSION D1(*),D2(*),D3(*),DP(*),U1(*),U2(*),U3(*),P(*)
      DIMENSION A(*),KCOLA(*),KLDA(*),B1(*),B2(*),B3(*),KCOLB(*),KLDB(*)
      DIMENSION KABD(*)
C
      SAVE
C
C
C-----------------------------------------------------------------------
C     Calculation of  D1=C2*D1+C1*(A*U1+B1*P)
C-----------------------------------------------------------------------
      CALL  LAX37 (A,KCOLA,KLDA,NU,U1,D1,C1,C2)
      IF (INEUM.EQ.1) CALL  BDRDEF(D1,KABD,NABD,0.5D0)
      CALL  LAX39 (B1,KCOLB,KLDB,NU,P,D1,C1,1D0)
      IF (INEUM.EQ.1) CALL  BDRDEF(D1,KABD,NABD,2D0)
C
C-----------------------------------------------------------------------
C     Calculation of  D2=C2*D2+C1*(A*U2+B2*P)
C-----------------------------------------------------------------------
      CALL  LAX37 (A,KCOLA,KLDA,NU,U2,D2,C1,C2)
C *** Compute  D:=D-B*P
      IF (INEUM.EQ.1) CALL  BDRDEF(D2,KABD,NABD,0.5D0)
      CALL  LAX39 (B2,KCOLB,KLDB,NU,P,D2,C1,1D0)
      IF (INEUM.EQ.1) CALL  BDRDEF(D2,KABD,NABD,2D0)
C
C-----------------------------------------------------------------------
C     Calculation of  D3=C2*D3+C1*(A*U3+B3*P)
C-----------------------------------------------------------------------
      CALL  LAX37 (A,KCOLA,KLDA,NU,U3,D3,C1,C2)
C *** Compute  D:=D-B*P
      IF (INEUM.EQ.1) CALL  BDRDEF(D3,KABD,NABD,0.5D0)
      CALL  LAX39 (B3,KCOLB,KLDB,NU,P,D3,C1,1D0)
      IF (INEUM.EQ.1) CALL  BDRDEF(D3,KABD,NABD,2D0)
C
C *** Set the Dirichlet-components of (D1,D2,D3) to zero
      CALL  BDRY0 (D1,D2,D3,KABD,NABD)
C
C-----------------------------------------------------------------------
C     Calculation of  DP=C2*DP+C1*(B1T*U1+B2T*U2+B3T*U3)
C-----------------------------------------------------------------------
      CALL  LTX39 (B1,KCOLB,KLDB,NU,NP,U1,DP,C1,C2)
      CALL  LTX39 (B2,KCOLB,KLDB,NU,NP,U2,DP,C1,1D0)
      CALL  LTX39 (B3,KCOLB,KLDB,NU,NP,U3,DP,C1,1D0)
C
C
C
      END
      
