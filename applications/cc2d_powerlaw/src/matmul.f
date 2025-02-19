************************************************************************
      SUBROUTINE  MATML (D1,D2,DP,U1,U2,P,C1,C2,A,KCOLA,KLDA,
     *                    B1,B2,KCOLB,KLDB,NU,NP,KMBD,NMBD,INEUM)
************************************************************************
*
*   Purpose:  performs the matrix-vector-operation
* 
*                   D:= C1*(A*U) + C2*D
*
*             for the vectors  D=(D1,D2,DP)  and  U=(U1,U2,P)  with the 
*             given scalar variables C1,C2
*  
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL A,B1,B2
C
      PARAMETER (NNVE=4)
      DIMENSION D1(*),D2(*),DP(*),U1(*),U2(*),P(*)
      DIMENSION A(*),KCOLA(*),KLDA(*),B1(*),B2(*),KCOLB(*),KLDB(*)
      DIMENSION KMBD(*)
C
      
      SAVE
C
      NA=KLDA(NU+1)-1
C
C-----------------------------------------------------------------------
C     Calculation of  D1=C2*D1+C1*(A*U12+B1*P)
C-----------------------------------------------------------------------
      CALL  LAX37 (A(1),KCOLA,KLDA,NU,U1,D1,C1,C2)
      CALL  LAX37 (A(NA+1),KCOLA,KLDA,NU,U2,D1,C1,1d0)
      IF (INEUM.EQ.1) CALL  BDRDEF(D1,KMBD,NMBD,0.5D0)
      CALL  LAX39 (B1,KCOLB,KLDB,NU,P,D1,C1,1D0)
      IF (INEUM.EQ.1) CALL  BDRDEF(D1,KMBD,NMBD,2D0)
C
C-----------------------------------------------------------------------
C     Calculation of  D2=C2*D2+C1*(A*U12+B2*P)
C-----------------------------------------------------------------------
      CALL  LAX37 (A(2*NA+1),KCOLA,KLDA,NU,U1,D2,C1,C2)
      CALL  LAX37 (A(3*NA+1),KCOLA,KLDA,NU,U2,D2,C1,1d0)
C *** Compute  D:=D-B*P
      IF (INEUM.EQ.1) CALL  BDRDEF(D2,KMBD,NMBD,0.5D0)
      CALL  LAX39 (B2,KCOLB,KLDB,NU,P,D2,C1,1D0)
      IF (INEUM.EQ.1) CALL  BDRDEF(D2,KMBD,NMBD,2D0)
C
C *** Set the Dirichlet-components of (D1,D2) to zero
      CALL  BDRY0 (D1,D2,KMBD,NMBD)
C
C-----------------------------------------------------------------------
C     Calculation of  DP=C2*DP+C1*(B1T*U1+B2T*U2)
C-----------------------------------------------------------------------
      CALL  LTX39 (B1,KCOLB,KLDB,NU,NP,U1,DP,C1,C2)
      CALL  LTX39 (B2,KCOLB,KLDB,NU,NP,U2,DP,C1,1D0)
C
C
      END
************************************************************************
      SUBROUTINE  MATML1 (D1,D2,DP,U1,U2,P,C1,C2,A,KCOLA,KLDA,
     *                    B1,B2,KCOLB,KLDB,NU,NP,KMBD,NMBD,INEUM)
************************************************************************
*
*   Purpose:  performs the matrix-vector-operation
* 
*                   D:= C1*(A*U) + C2*D
*
*             for the vectors  D=(D1,D2,DP)  and  U=(U1,U2,P)  with the 
*             given scalar variables C1,C2
*  
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL A,B1,B2
C
      PARAMETER (NNVE=4)
      DIMENSION D1(*),D2(*),DP(*),U1(*),U2(*),P(*)
      DIMENSION A(*),KCOLA(*),KLDA(*),B1(*),B2(*),KCOLB(*),KLDB(*)
      DIMENSION KMBD(*)
C
      
      SAVE
C
      NA=KLDA(NU+1)-1
C
C-----------------------------------------------------------------------
C     Calculation of  D1=C2*D1+C1*(A*U12+B1*P)
C-----------------------------------------------------------------------
      CALL  LAX37 (A(1),KCOLA,KLDA,NU,U1,D1,C1,C2)
      CALL  LAX37 (A(NA+1),KCOLA,KLDA,NU,U2,D1,C1,1d0)
      IF (INEUM.EQ.1) CALL  BDRDEF(D1,KMBD,NMBD,0.5D0)
      CALL  LAX39 (B1,KCOLB,KLDB,NU,P,D1,C1,1D0)
      IF (INEUM.EQ.1) CALL  BDRDEF(D1,KMBD,NMBD,2D0)
C
C-----------------------------------------------------------------------
C     Calculation of  D2=C2*D2+C1*(A*U12+B2*P)
C-----------------------------------------------------------------------
      CALL  LAX37 (A(2*NA+1),KCOLA,KLDA,NU,U1,D2,C1,C2)
      CALL  LAX37 (A(3*NA+1),KCOLA,KLDA,NU,U2,D2,C1,1d0)
C *** Compute  D:=D-B*P
      IF (INEUM.EQ.1) CALL  BDRDEF(D2,KMBD,NMBD,0.5D0)
      CALL  LAX39 (B2,KCOLB,KLDB,NU,P,D2,C1,1D0)
      IF (INEUM.EQ.1) CALL  BDRDEF(D2,KMBD,NMBD,2D0)
C
C *** Set the Dirichlet-components of (D1,D2) to zero
      CALL  BDRY0 (D1,D2,KMBD,NMBD)
C
C-----------------------------------------------------------------------
C     Calculation of  DP=C2*DP+C1*(B1T*U1+B2T*U2)
C-----------------------------------------------------------------------
c      CALL  LTX39 (B1,KCOLB,KLDB,NU,NP,U1,DP,C1,C2)
c      CALL  LTX39 (B2,KCOLB,KLDB,NU,NP,U2,DP,C1,1D0)
C
C
      END
************************************************************************
      SUBROUTINE  MATML2 (D1,D2,DP,U1,U2,P,C1,C2,A,KCOLA,KLDA,
     *                    B1,B2,KCOLB,KLDB,NU,NP,KMBD,NMBD,INEUM)
************************************************************************
*
*   Purpose:  performs the matrix-vector-operation
* 
*                   D:= C1*(A*U) + C2*D
*
*             for the vectors  D=(D1,D2,DP)  and  U=(U1,U2,P)  with the 
*             given scalar variables C1,C2
*  
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL A,B1,B2
C
      PARAMETER (NNVE=4)
      DIMENSION D1(*),D2(*),DP(*),U1(*),U2(*),P(*)
      DIMENSION A(*),KCOLA(*),KLDA(*),B1(*),B2(*),KCOLB(*),KLDB(*)
      DIMENSION KMBD(*)
C
      
      SAVE
C
      NA=KLDA(NU+1)-1
C
C
C-----------------------------------------------------------------------
C     Calculation of  DP=C2*DP+C1*(B1T*U1+B2T*U2)
C-----------------------------------------------------------------------
      CALL  LTX39 (B1,KCOLB,KLDB,NU,NP,U1,DP,C1,C2)
      CALL  LTX39 (B2,KCOLB,KLDB,NU,NP,U2,DP,C1,1D0)
C
C
      END
