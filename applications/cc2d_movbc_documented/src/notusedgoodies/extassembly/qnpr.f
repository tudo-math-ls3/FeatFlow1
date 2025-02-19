***********************************************************************
* Creation of the prolongation- and restriction-matrix
*
* This routine creates the prolongation/restriction matrix that is
* needed for the prolongation/restriction of Qn.
*
* In:
*   ELE   - The element
*   NDFC  - Number of local degrees of freedom of the element used in
*           the coarse-grid triangulation
*   NDFF  - Number of local degrees of freedom of the element used in
*           the fine-grid triangulation; usually = NDFC
*   COORC - Coordinates of the DOF's on the reference element
*           in the coarse-grid FE-space
*   COORF - Coordinates of the DOF's on the reference element
*           in the fine-grid FE-space; usually = COORC except for
*           when different FE-spaces are used when switching levels.
*
* Out:
*   PMAT  - array [1..NDFC,1..NDFF,1..4] of double
*           prolongation/restriction matrix
***********************************************************************

      SUBROUTINE PRCMAT (ELE,NDFC,NDFF,COORC,COORF,PMAT)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      
      EXTERNAL ELE
      
      INTEGER NDFC, NDFF
      DOUBLE PRECISION COORC(2,NDFC), COORF(2,NDFF)
      DOUBLE PRECISION PMAT(NDFF,NDFC,4)
      
      INTEGER I,J
      DOUBLE PRECISION X,Y

C Initialise the element-COMMON blocks to directly evaluate
C on the reference element:

      DJAC(1,1) = 1D0
      DJAC(1,2) = 0D0
      DJAC(1,1) = 0D0
      DJAC(2,2) = 1D0
      DETJ = 1D0
      
C Only the function value is necessary
      
      DO I=1,NNDER
        BDER(I) = .FALSE.
      END DO
      BDER (1) = .TRUE.
      
C The weights=entries in the prolongation matrix are created
C by simple evaluation of the element routine in the degrees of
C freedom. We have one coarse-grid element containing 4 fine-grid
C sub-elements:
C
C 4-----+-----3
C |     |     |
C |     |     |
C |     |     |
C +-----+-----+
C |     |     |
C |     |     |
C |     |     |
C 1-----+-----2
C
C The DOF's 1..4 are the same for all Qn-elements, they mark the corners.
C We perform a loop through the smaller sub-elements. On each sub-element
C we evaluate the element in all fine-grid DOF's, e.g. for the first
C sub-element with Q1:
C
C +-----+-----+
C |           |
C |           |
C |           |
C 4-----3     +
C |     |     |
C |     |     |
C |     |     |
C 1-----2-----+
C
C The weights that are calculated by the element are then stored 
C in the first PMAT-submatrix, so this matrix contains all weights
C for the "fully" interpolation - since the prolongation is nothing
C else than an interpolation of the coarse-grid space into
C the fine-grid space.
C PMAT (.,.,1), the values of the other sub-elements in the other
C three PMAT-submatrices PMAT(.,.,2..4).
C
C The order in which the evaluated is a little bit special. We store
C the values in the order of the local degrees of freedom, so we
C have to take care about the rotation of the sub-elements when
C evaluating in the coordinates of the DOF's of the fine-grid.
C The second sub-element for example has the local DOF's
C
C +-----+-----+
C |           |
C |           |
C |           |
C +     3-----2
C |     |     |
C |     |     |
C |     |     |
C +-----4-----1
C
C which is the reference-element rotated by 90°!
      
C Build the (local) prolongation matrix rowwise
      
      DO I=1,NDFF
        
C Project the fine-grid point from the fine-grid reference element
C into the coarse-grid reference element by scaling + shifting such
C that the local dof 1 is the same in both elements.
C
C Evaluate the element in that point on the fine grid.

        X = 0.5D0*COORF(1,I) - 0.5D0
        Y = 0.5D0*COORF(2,I) - 0.5D0
        CALL ELE (X,Y,0)
      
C The element returns the weights of the prolongation in the
C DBAS-array in the /ELEM/-COMMON block (as function values of the
C basis functions). Copy these values to the prolongation matrix:

        DO J=1,NDFC
          PMAT (I,J,1) = DBAS (J,1)
        END DO

C The same for the other fine-grid elements that the coarse grid
C element contain. Take care about the rotation of the fine-grid
C elements because the number of the dof's are not so easy
C to determine.
C
C 2nd smaller element = reference element rotated by 90°:

        X = -0.5D0*COORF(2,I) + 0.5D0
        Y = 0.5D0*COORF(1,I) - 0.5D0
        CALL ELE (X,Y,0)
      
        DO J=1,NDFC
          PMAT (I,J,2) = DBAS (J,1)
        END DO

C 3rd smaller element = reference element rotated by 180°:

        X = -0.5D0*COORF(1,I) + 0.5D0
        Y = -0.5D0*COORF(2,I) + 0.5D0
        CALL ELE (X,Y,0)
      
        DO J=1,NDFC
          PMAT (I,J,3) = DBAS (J,1)
        END DO

C 4th smaller element = reference element rotated by 270°:

        X = 0.5D0*COORF(2,I) - 0.5D0
        Y = -0.5D0*COORF(1,I) + 0.5D0
        CALL ELE (X,Y,0)
      
        DO J=1,NDFC
          PMAT (I,J,4) = DBAS (J,1)
        END DO
      
      END DO
      
      END
      
***********************************************************************
* Prolongation
*
* This implements general prolongation for quadrilateral elements.
*
* In:
*   IELTP  - Type of element (11=Q1, 13=Q2, ...)
*   CVEC   - array [1..NEQC] of double
*            Coarse grid vector - to be prolongated
*   NEQC   - Length of CVEC
*   NEQF   - Length of FVEC
*   NDFC   - Number of local degrees of freedom on the coarse grid
*   NDFF   - Number of local degrees of freedom on the fine grid
*   PMAT   - The prolongation matrix; created by PRCMAT
*   TRIAC  - STRIA-triangulation structure on the coarse grid
*   TRIAF  - STRIA-triangulation structure on the fine grid
*   KVERTC - Vertex array on the coarse grid. Must point to the vertex
*            array identified by the handle in TRIAC, otherwise the
*            result of this routine is undefined. This is for
*            speed reasons.
*   KMIDC  - Midpoint array on the coarse grid. Must point to the 
*            midpoint array identified by the handle in TRIAC, otherwise 
*            the result of this routine is undefined. This is for
*            speed reasons.
*   KVERTF - Vertex array on the fine grid. Must point to the vertex
*            array identified by the handle in TRIAF, otherwise the
*            result of this routine is undefined. This is for
*            speed reasons.
*   KMIDF  - Midpoint array on the fine grid. Must point to the 
*            midpoint array identified by the handle in TRIAF, otherwise 
*            the result of this routine is undefined. This is for
*            speed reasons.
*   KADJF  - Adjacent-information array on the fine grid. Must point to  
*            the KADJ-array identified by the handle in TRIAF, otherwise 
*            the result of this routine is undefined. This is for
*            speed reasons.
*   KNAEF  - array [1..NEQF] of integer
*            Number of elements meeting in a vertex of the fine grid.
*
* Out:
*   FVEC   - array [1..NEQF] of double
*            The prolongated vector.
***********************************************************************
      
      SUBROUTINE PROLQN (IELTP,CVEC, FVEC, NEQC, NEQF, NDFC, NDFF, PMAT,
     *                   TRIAC, TRIAF,
     *                   KVERTC, KMIDC, KVERTF, KMIDF, 
     *                   KADJF, KNAEF)
      
      IMPLICIT NONE

      INCLUDE 'cbasictria.inc'
      INCLUDE 'cbasicelem.inc'
      
      INCLUDE 'stria.inc'
      
      INTEGER IELTP, NEQC, NEQF, NDFC, NDFF
      INTEGER KVERTC(NNVE,*), KVERTF(NNVE,*), KNAEF(NEQF)
      INTEGER KADJF(NNVE,*), KMIDC (NNVE,*), KMIDF (NNVE,*)
      INTEGER TRIAC(SZTRIA), TRIAF(SZTRIA)
      DOUBLE PRECISION CVEC(NEQC), FVEC(NEQF), PMAT(NDFF,NDFC,4)
      
      INTEGER IOFF, IOFC, IELH(4), IELC, IELF, IVC, IVF, I
      INTEGER KDFGC(NNBAS), KDFLC(NNBAS)
      INTEGER KDFGF(NNBAS), KDFLF(NNBAS)
      
C Clear the "target"-vector

      CALL LCL1(FVEC,NEQF)
      
C Loop through the elements of the coarse grid

      DO IELC = 1,TRIAC(ONEL)
      
C Determine the global numbers of the DOF's. We need them sorted in
C the order of the local DOF's, as the prolongation/restriction
C matrix was build according to these. So determine local/global
C DOF's and resort them so that the local DOF's are in consecutive
C order (1,2,3,...):

        CALL NDFGLX(TRIAC,IELC,1,IELTP,KVERTC,KMIDC,KDFGC,KDFLC)
        CALL NGLSD(KDFLC,KDFGC,NDFC)

C Determine the fine-grid elements that reside in our coarse-grid
C element:

        IELH(1) = IELC
        IELH(2) = KADJF(2,IELH(1))
        IELH(3) = KADJF(2,IELH(2))
        IELH(4) = KADJF(2,IELH(3))
      
C Now loop through the elements on the fine grid.
C IELF is something like a "local" element number 1..4
C corresponding to the 4 subelements in the coarse-grid element.
C The "real" number of the subelement is then IELH(IELF).

        DO IELF = 1,4
          
          IF (IELH(IELF).NE.0) THEN
          
C Determine the global numbers of the DOF's. We need them sorted in
C the order of the local DOF's, as the prolongation/restriction
C matrix was build according to these. So determine local/global
C DOF's and resort them so that the local DOF's are in consecutive
C order (1,2,3,...):

            CALL NDFGLX(TRIAF,
     *                  IELH(IELF),1,IELTP,KVERTF,KMIDF,KDFGF,KDFLF)
            CALL NGLSD(KDFLF,KDFGF,NDFF)
          
C Ok, we are in the fine grid element IELF residing in the
C coarse grid element IELC. Loop through all DOF's of this
C fine grid element to calculate the values in these DOF's:

            DO IOFF = 1,NDFF
          
C The current fine-grid vertex corresponding to our current
C DOF is...
      
              IVF = KDFGF(IOFF)
            
C Multiply the values in the DOF's of the coarse grid with the
C weights in the prolongation matrix to obtain the values in the
C fine-grid DOF; sum them up to the full value in each fine-grid DOF:
  
              DO IOFC = 1,NDFC
            
C The vertex number of the current coarse grid vertex is...
 
                IVC = KDFGC(IOFC)
            
C We have 4 submatrices in the prolongation matrix - one for each
C sub-element in the coarse grid element. The current prolongation
C sub-matrix corresponts to the current sub-element IELF
            
                FVEC (IVF) = FVEC (IVF) +
     *             PMAT (IOFF,IOFC,IELF) * CVEC (IVC)
            
C The prolongation therefore implements a simple interpolation
C between the spaces. But from another point of view you can
C even look at the prolongation as a "weighted distribution" of 
C the value in one node to the values in the other nodes.
C E.g. for the Q1-element the coarse-grid node C1 distributes
C its value to the fine grid nodes Fi with the weights like in 
C the following picture:
C
C        C4                   C3
C          +--------+--------+
C          |                 |
C          |                 |
C          |                 |
C   0.5 +> F4-------F3       +
C       |  |    ___^|        |
C       |  |  _/0.25|        |
C       |  | /      |        |
C       \  F1-------F2-------+
C        C1-^ 1.0   ^         C2
C        |----------|0.5
            
              END DO
              
            END DO

          END IF
        
        END DO
        
      END DO

C So far so good. Unfortunately the values are not correct for the
C moment. We must be aware of that we summed up some values more than
C once. To be exact, we summed up every value as often as there
C are elements adjacent to the DOF. So by dividing the value at
C each DOF by the number of adjacent vertices we scale the values
C back to the correct magnitude.
C Furthermore we added the value of every coarse grid DOF 4x (for the
C four passes of the four sub-elements), so we have to divide the 
C value here by 4, too!

      DO I = 1,NEQF 
        FVEC(I) = FVEC(I) / DBLE(KNAEF(I))
      END DO

      END
      
      
***********************************************************************
* Restriction
*
* This implements general restriction for quadrilateral elements.
*
* In:
*   IELTP  - Type of element (11=Q1, 13=Q2, ...)
*   FVE    - array [1..NEQF] of double
*            The fine-grid vector - to be restricted.
*   NEQC   - Length of CVEC
*   NEQF   - Length of FVEC
*   NDFC   - Number of local degrees of freedom on the coarse grid
*   NDFF   - Number of local degrees of freedom on the fine grid
*   PMAT   - The prolongation matrix; created by PRCMAT
*   TRIAC  - STRIA-triangulation structure on the coarse grid
*   TRIAF  - STRIA-triangulation structure on the fine grid
*   KVERTC - Vertex array on the coarse grid. Must point to the vertex
*            array identified by the handle in TRIAC, otherwise the
*            result of this routine is undefined. This is for
*            speed reasons.
*   KMIDC  - Midpoint array on the coarse grid. Must point to the 
*            midpoint array identified by the handle in TRIAC, otherwise 
*            the result of this routine is undefined. This is for
*            speed reasons.
*   KVERTF - Vertex array on the fine grid. Must point to the vertex
*            array identified by the handle in TRIAF, otherwise the
*            result of this routine is undefined. This is for
*            speed reasons.
*   KMIDF  - Midpoint array on the fine grid. Must point to the 
*            midpoint array identified by the handle in TRIAF, otherwise 
*            the result of this routine is undefined. This is for
*            speed reasons.
*   KADJF  - Adjacent-information array on the fine grid. Must point to  
*            the KADJ-array identified by the handle in TRIAF, otherwise 
*            the result of this routine is undefined. This is for
*            speed reasons.
*   KNAEF  - array [1..NEQF] of integer
*            Number of elements meeting in a vertex of the fine grid.
*
* Out:
*   CVEC   - array [1..NEQC] of double
*            The restricted FVEC-vector.
***********************************************************************
      
      SUBROUTINE RESTQN (IELTP,CVEC, FVEC, NEQC, NEQF, NDFC, NDFF, PMAT,
     *                   TRIAC, TRIAF, 
     *                   KVERTC, KMIDC, KVERTF, KMIDF, KADJF, KNAEF)
      
      IMPLICIT NONE

      INCLUDE 'cbasictria.inc'
      INCLUDE 'cbasicelem.inc'
      
      INCLUDE 'stria.inc'
      
      INTEGER IELTP, NEQC, NEQF, NDFC, NDFF
      INTEGER KVERTC(NNVE,*), KVERTF(NNVE,*), KNAEF(NEQF)
      INTEGER KADJF(NNVE,*), KMIDC (NNVE,*), KMIDF (NNVE,*)
      DOUBLE PRECISION CVEC(NEQC), FVEC(NEQF), PMAT(NDFF,NDFC,4)
      
      INTEGER TRIAC(SZTRIA), TRIAF(SZTRIA)
      
      INTEGER IOFF, IOFC, IELH(4), IELC, IELF, IVC, IVF
      DOUBLE PRECISION WEIGHT
      INTEGER KDFGC(NNBAS), KDFLC(NNBAS)
      INTEGER KDFGF(NNBAS), KDFLF(NNBAS)
      
C Clear the "target"-vector

      CALL LCL1(CVEC,NEQC)

C Loop through the elements of the coarse grid

      DO IELC = 1,TRIAC(ONEL)
      
C Determine the global numbers of the DOF's. We need them sorted in
C the order of the local DOF's, as the prolongation/restriction
C matrix was build according to these. So determine local/global
C DOF's and resort them so that the local DOF's are in consecutive
C order (1,2,3,...):

        CALL NDFGLX(TRIAC,IELC,1,IELTP,KVERTC,KMIDC,KDFGC,KDFLC)
        CALL NGLSD(KDFLC,KDFGC,NDFC)

C Determine the fine-grid elements that reside in our coarse-grid
C element:

        IELH(1) = IELC
        IELH(2) = KADJF(2,IELH(1))
        IELH(3) = KADJF(2,IELH(2))
        IELH(4) = KADJF(2,IELH(3))
      
C Now loop through the elements on the fine grid.
C IELF is something like a "local" element number 1..4
C corresponding to the 4 subelements in the coarse-grid element.
C The "real" number of the subelement is then IELH(IELF).

        DO IELF = 1,4
          
          IF (IELH(IELF).NE.0) THEN
          
C Determine the global numbers of the DOF's. We need them sorted in
C the order of the local DOF's, as the prolongation/restriction
C matrix was build according to these. So determine local/global
C DOF's and resort them so that the local DOF's are in consecutive
C order (1,2,3,...):

            CALL NDFGLX(TRIAF,
     *                  IELH(IELF),1,IELTP,KVERTF,KMIDF,KDFGF,KDFLF)
            CALL NGLSD(KDFLF,KDFGF,NDFF)
          
C Ok, we are in the fine grid element IELF residing in the
C coarse grid element IELC. Loop through all DOF's of this
C fine grid element to calculate the values in these DOF's:

            DO IOFF = 1,NDFF
          
C The current fine-grid vertex corresponding to our current
C DOF is...
      
              IVF = KDFGF(IOFF)

C Calculate the additional "weight" from the number of adjacent
C elements to the vertex IVF. In the prolongation we had been 
C lucky that this weight was constant through the whole calculation,
C so we were able to incorporate them quickly afterwards.
C Here we are not so lucky, we have to incorporate the weight
C during the caluculation, so restriction is a little bit more
C expensive.
C This "weight" is again the number of adjacent elements to
C a fine-grid node. This is due tu the fact that we loop
C over the fine grid elements and therefore we "distribute" each 
C fine-grid-node as often to the coarse-grid nodes as there are
C adjacent elements to it.
            
              WEIGHT = 1D0/DBLE(KNAEF(IVF))

C Multiply the values in the DOF's of the fine grid with the
C weights in the restriction matrix to obtain the values in the
C coarse-grid DOF; sum them up to the full value in each coarse-grid
C DOF:

              DO IOFC = 1,NDFC
            
C The vertex number of the current coarse grid vertex is...

                IVC = KDFGC (IOFC)
            
C We have 4 submatrices in the restriction matrix - one for each
C sub-element in the coarse grid element. The current restriction
C sub-matrix corresponts to the current sub-element IELF.
C The restriction matrix is exactly the transposed of the
C prolongation matrix, so we don't have to use a specially
C designed matrix here.
C
C But in contrast to the prolongation we cannot sum the values
C up and obtain the correct value by a postprocessing-like
C multiplication with the number of adjacent elements. This
C time we have to include this number directly in the formula
C as this "weight" is not constant for all values incorporated
C into the coarse-grid value!

                CVEC (IVC) = CVEC (IVC) +
     *             PMAT (IOFF,IOFC,IELF) * FVEC (IVF) * WEIGHT

C The restriction is the "adjoint" operator of the prolongation.   
C Its matrix defining how to "create" the coarse grid right-hand-side
C by weighted summation of the fine-grid right-hand-side 
C (not solution! This is a *restriction*, not an *interpolation*
C to the lower level! We are working in the dual space!)
C is defined by the transposed matrix of the prolongation operator
C (therefore we can use PMAT here, but with transposed indices).
C While the prolongation can be viewed as a "distribution" of the
C weights of one node to the neighbour nodes, the restriction
C is the "collection" of the values in the fine-grid nodes
C to the coarse grid node - with the same weights as the prolongation.
C
C
C        C4                   C3
C          +--------+--------+
C          |                 |
C          |                 |
C          |                 |
C      +-- F4-------F3       +
C      |   |    ___/|        |
C      |   | _/0.25 |        |
C      |   |v       |        |
C  0.5 |   F1-------F2-------+
C      +> C1<| 1.0   |         C2
C         ^----------| 0.5
            
              END DO
          
            END DO
            
          END IF
        
        END DO
        
      END DO

      END 

************************************************************************
* DOF's-per-element counting routine
*
* This routine counts how many elements meet in each global DOF.
*
* In:
*   TRIA  : array [1..SZTRIA] of integer
*           Triangulation structure of the underlying grid
*   ELE   : Element subroutine
* Out:
*   LAE   : Handle to array [1..NEQ] of integer
*           KAE(I) is the number of elements meeting in a global DOF I;
*           this is 1 for points in elements, 2 for points on edges. 
*           For corner vertices this is >= 1.
************************************************************************

      SUBROUTINE GDFPEC (TRIA,ELE,LAE)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'stria.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'cbasicelem.inc'
      
      INTEGER TRIA(SZTRIA),LAE
      EXTERNAL ELE
      
C     local variables

      INTEGER IEL,IELTYP
      INTEGER NEQ,KAE,I,KVERT,KMID
      INTEGER KDFG(NNBAS),KDFL(NNBAS),IDFL
      
C     externals

      INTEGER NDFGX,NDFL
      EXTERNAL NDFGX,NDFL
      
C     Get the type of the element

      IELTYP = -1
      CALL ELE(0D0,0D0,IELTYP)
      
C     Get the number of global DOF's

      NEQ = NDFGX(IELTYP,TRIA)
      
C     Allocate memory

      CALL ZNEW(NEQ,3,LAE,'KAE   ')
      KAE   = L(LAE)
      
C     Get the number of local DOF's

      IDFL = NDFL(IELTYP)
      
C     Loop through the elements     

      KVERT = L(TRIA(OLVERT))
      KMID  = L(TRIA(OLMID)) 
      
      DO IEL = 1,TRIA(ONEL)

C       The global DOF's are e.g. numbered as:
C
C           4 --10---9-- 3
C           |            |
C          11   16  15   8
C           |            |
C          12   13  14   7
C           |            |
C           1 ---5---6-- 2
C
C       Get the local and global DOF's into KDFG on this element.

        CALL NDFGLX(TRIA,IEL,0,IELTYP,KWORK(KVERT),KWORK(KMID),
     *              KDFG,KDFL)
     
C       Increase the KAE-value for all vertices on this element by one.
C       So at the end, KAE is the number of elements meeting in a 
C       DOF.

        DO I=1,IDFL
          KWORK(KAE+KDFG(I)-1) = KWORK(KAE+KDFG(I)-1) + 1
        END DO

      END DO
        
C     That's it.
        
      END

!***********************************************************************
!* Build KAE for corner vertices
!***********************************************************************
!
!      SUBROUTINE BKNAEV (KNAEF, KVEL, NVT, NVEL)
!      
!      IMPLICIT NONE
!      
!      INTEGER NVT,NVEL
!      INTEGER KNAEF(*), KVEL(NVEL,*)
!      
!      INTEGER I,J    
!      
!      DO I=1,NVT
!        KNAEF(I) = 0
!        DO J=1,NVEL
!          IF (KVEL(J,I).NE.0) KNAEF(I)=KNAEF(I)+1
!        END DO
!      END DO
!      
!      END
!      
!***********************************************************************
!* Build KAE for vertices on edges 
!***********************************************************************
!
!      SUBROUTINE BKNAEM (KNAEF, KMEL, NVT, NMT)
!      
!      IMPLICIT NONE
!      
!      INTEGER KNAEF(*),KMEL(2,*)
!      INTEGER NVT,NMT
!      
!      INTEGER I,J
!      
!      DO I=1,NMT
!        KNAEF(NVT+I) = 0
!        DO J=1,2
!          IF (KMEL(J,I).NE.0) KNAEF(NVT+I)=KNAEF(NVT+I)+1
!        END DO
!      END DO
!      
!      END
!      
!***********************************************************************
!***********************************************************************
!           
!      SUBROUTINE BKNAEV (KNAEF, KVEL, NVT, NVEL)
!      
!      IMPLICIT NONE
!      
!      INTEGER NVT,NVEL
!      INTEGER KNAEF(*), KVEL(NVEL,*)
!      
!      INTEGER I,J    
!      
!      DO I=1,NVT
!        KNAEF(I) = 0
!        DO J=1,NVEL
!          IF (KVEL(J,I).NE.0) KNAEF(I)=KNAEF(I)+1
!        END DO
!      END DO
!      
!      END
!      
!
!      SUBROUTINE BKNAEM (KNAEF, KMEL, NVT, NMT)
!      
!      IMPLICIT NONE
!      
!      INTEGER KNAEF(*),KMEL(2,*)
!      INTEGER NVT,NMT
!      
!      INTEGER I,J
!      
!      DO I=1,NMT
!        KNAEF(NVT+I) = 0
!        DO J=1,2
!          IF (KMEL(J,I).NE.0) KNAEF(NVT+I)=KNAEF(NVT+I)+1
!        END DO
!      END DO
!      
!      END
!      
!            
!      SUBROUTINE TESTPR
!      
!      IMPLICIT NONE
!      
!      INCLUDE 'cmem.inc'
!      
!      INCLUDE 'cbasicmg.inc'
!      INCLUDE 'cmgpar.inc'
!      INCLUDE 'cmgtria.inc'
!
!      INCLUDE 'commonsini.for'
!      
!      INCLUDE 'stria.inc'
!
!      INTEGER KLA,KLCOLA,KLLDA,KLU,KLB,KLD,KLAUX
!      COMMON /MGFLD/  KLA(NNLEV),KLCOLA(NNLEV),
!     *                KLLDA(NNLEV),KLU(NNLEV),
!     *                KLB(NNLEV),KLD(NNLEV),KLAUX(NNLEV)
!     
!      INTEGER KNA,KNEQ
!      COMMON /MGDIM/  KNA(NNLEV),KNEQ(NNLEV)
!      
!      INCLUDE 'cbasictria.inc'
!      INCLUDE 'ctria.inc'
!      
!      INTEGER I,KAE
!      INTEGER TRIA1(SZTRIA), TRIA2(SZTRIA)
!      DOUBLE PRECISION COORC(2,9)
!      DATA COORC/-1,-1, 1,-1, 1,1, -1,1,
!     *            0,-1, 1,0,  0,1, -1,0,
!     *            0,0 /
!      
!      DOUBLE PRECISION PMAT(9,9,4)
!      
!      EXTERNAL EXXXQ
!      
!      I=2
!      CALL SETLEV (NLMAX-1,I)
!      CALL C2TRIA(TRIA1)
!
!      CALL SETLEV (NLMAX,I)
!      CALL C2TRIA(TRIA2)
!      
!      CALL ZNEW (NVT+NMT,3,KAE,'KAE   ')
!      DO I=1,NVT+NMT+NEL
!        KWORK(L(KAE)+I-1) = 1
!      END DO
!      
!      CALL BKNAEV (KWORK(L(KAE)), KWORK(L(KLVEL(ILEV))),KNVT(ILEV), 
!     *             NVEL)
!      CALL BKNAEM (KWORK(L(KAE)), KWORK(L(KLMEL(ILEV))),KNVT(ILEV),
!     *             KNMT(ILEV))
!      
!      CALL LCL1(DWORK(L(KLU(ILEV-1))),KNEQ(ILEV-1))
!      
!      DO I=1,KNEQ(ILEV)
!        DWORK(L(KLU(ILEV))+I-1) = 1D0
!      END DO
!      
!      CALL PRCMAT (EXXXQ,9,9,COORC,COORC,PMAT)
!      
!      CALL MR013(DWORK(L(KLU(ILEV))),DWORK(L(KLU(ILEV-1))),
!     *           KWORK(L(KLVERT(ILEV))),KWORK(L(KLVERT(ILEV-1))),
!     *           KWORK(L(KLMID(ILEV))),KWORK(L(KLMID(ILEV-1))),
!     *           KWORK(L(KLADJ(ILEV))),KWORK(L(KLADJ(ILEV-1))),
!     *           KNEL(ILEV),KNEL(ILEV-1),
!     *           KNVT(ILEV),KNVT(ILEV-1),
!     *           KNMT(ILEV),KNMT(ILEV-1))
!      
!      CALL RESTQN (IELTYP, DWORK(L(KLU(ILEV-1))), DWORK(L(KLU(ILEV))), 
!     *             KNEQ(ILEV-1), KNEQ(ILEV), 9, 9, PMAT, 
!     *             TRIA1, TRIA2,
!     *             KWORK(L(KLVERT(ILEV-1))), KWORK(L(KLMID(ILEV-1))), 
!     *             KWORK(L(KLVERT(ILEV))), KWORK(L(KLMID(ILEV))), 
!     *             KWORK(L(KLADJ(ILEV))), KWORK(L(KAE)))      
!      
!      END 
      