************************************************************************
* Matrix-vector-multiplication
*
* This performs the following matrix vector multiplication:
* 
*          D:= C1*(A*U) + C2*D
*
* for the vectors  D=(D1,D2,DP)  and  U=(U1,U2,P)  with the 
* given scalar variables C1,C2. The matrix A is the linearised
* nonlinear system matrix for velocity and pressure.
*
* D is assumed to be a defect vector. All entries corresponding to
* Dirichlet nodes are set to 0.0.
*
* In:
*   D1     - array [1..NU] of double
*   D2     - array [1..NU] of double
*   DP     - array [1..NP] of double
*            D = (D1,D2,DP) is the first vector
*   U1     - array [1..NU] of double
*   U2     - array [1..NU] of double
*   P      - array [1..NP] of double
*            U = (U1,U2,P) is the second vector
*   C1,
*   C2     - The constants
*   A,
*   KCOLA,
*   KLDA   - Velocity matrix
*   B1     - First B-matrix
*   B2     - Second B-matrix
*   KCOLB,
*   KLDB   - Structure of the two B-matrices, matrix-structure 9
*   NU     - Number of velocity components in Ux and Dx
*   NP     - Number of pressure components in DP and P 
*   NMBD   - Number of boundary vertices
*   KMBD   - array [1..NMBD] of integer
*            "Shortcut nodal property" array. Array with the vertex
*            numbers of all boundary nodes. KMBD(.) > 0 for Dirichlet
*            nodes, KMBD(.) < 0 for Neumann type nodes
*   INEUM  - =0, if there are no Neumann nodes
*            =1, if there are Neumann nodes
*
* Out:
*   D      = C1*(A*U) + C2*D
*
* Warning:
*   Due to a bug in VANCA, this routine actually performs
*            D:= 0.5 * C1*(A*U)  +  0.5 * C2*D
*   for all Neumann boundary nodes! This compensates the accidently 
*   introduced factor 2.0 in the VANCA smoothing/solving!
************************************************************************

      SUBROUTINE  MTMUL (D1,D2,DP,U1,U2,P,C1,C2,A,KCOLA,KLDA,
     *                    B1,B2,KCOLB,KLDB,NU,NP,KMBD,NMBD,INEUM)
      IMPLICIT NONE

C common blocks

      INCLUDE 'cbasictria.inc'

C parameters

      DOUBLE PRECISION A(*),B1(*),B2(*)

      DOUBLE PRECISION D1(*),D2(*),DP(*),U1(*),U2(*),P(*)
      INTEGER KCOLA(*),KLDA(*),KCOLB(*),KLDB(*),KMBD(*)
      
      DOUBLE PRECISION C1,C2
      INTEGER NU,NP,NMBD,INEUM
      
C-----------------------------------------------------------------------
C     Calculation of  D1=C2*D1+C1*(A*U1+B1*P)
C-----------------------------------------------------------------------

C     Multiply with A:

      CALL  LAX17 (A,KCOLA,KLDA,NU,U1,D1,C1,C2)
      
C     Halfen the value of Neumann boundary nodes before multiplying
C     with B to compensate the VANCA bug:
      
      IF (INEUM.EQ.1) CALL  BDRDEF(D1,KMBD,NMBD,0.5D0)
      
C     Multiply with B:
      
      CALL  LAX19 (B1,KCOLB,KLDB,NU,P,D1,C1,1D0)
      
C     Restore the old value in D for the Neumann boundary nodes by
C     multiplying with 2.0 again...
      
      IF (INEUM.EQ.1) CALL  BDRDEF(D1,KMBD,NMBD,2D0)

C-----------------------------------------------------------------------
C     Calculation of  D2=C2*D2+C1*(A*U2+B2*P)
C-----------------------------------------------------------------------

C     The same handling as above for the 2nd velocity component...

      CALL  LAX17 (A,KCOLA,KLDA,NU,U2,D2,C1,C2)
      IF (INEUM.EQ.1) CALL  BDRDEF(D2,KMBD,NMBD,0.5D0)
      CALL  LAX19 (B2,KCOLB,KLDB,NU,P,D2,C1,1D0)
      IF (INEUM.EQ.1) CALL  BDRDEF(D2,KMBD,NMBD,2D0)

C     Set the Dirichlet-components of the defect vector
C     (D1,D2) to zero

      CALL  BDRY0 (D1,D2,KMBD,NMBD)

C-----------------------------------------------------------------------
C     Calculation of  DP=C2*DP+C1*(B1T*U1+B2T*U2)
C-----------------------------------------------------------------------

      CALL  LTX19 (B1,KCOLB,KLDB,NU,NP,U1,DP,C1,C2)
      CALL  LTX19 (B2,KCOLB,KLDB,NU,NP,U2,DP,C1,1D0)

      END
      
