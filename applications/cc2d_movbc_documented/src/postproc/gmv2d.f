C AREA  - array [1..NEL] of double. Area of each cell.

      SUBROUTINE XGMV2D(MUNIT,CFILE,NEL,NVT,KVERT,DCORVG,KNPR,
     *                  VU,VV,VP,VISO,AREA,TIMENS)
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'

C parameters

      INTEGER MUNIT, NEL, NVT
      CHARACTER CFILE*(*)
      INTEGER KVERT(4,*),KNPR(*)
      DOUBLE PRECISION DCORVG(2,*),TIMENS
      REAL VU(*),VV(*),VP(*),VISO(*)
      DOUBLE PRECISION AREA(NEL)
      
C externals

      INTEGER NFBDYC
      EXTERNAL NFBDYC

C local variables

      INTEGER IFIRST, IVT, IEL, FLAG

      DO 1 IFIRST = 1,60
      IF (CFILE(IFIRST:IFIRST).EQ.' ') GOTO 2
1     CONTINUE
2     WRITE(CFILE(IFIRST:IFIRST+4),'(A4)') '.gmv'
C
      OPEN (UNIT=MUNIT,FILE=CFILE)
C
C
      WRITE(MUNIT,'(A)')'gmvinput ascii'
      WRITE(MUNIT,*)'nodes ',NVT
C
      DO 100 IVT=1,NVT
100   WRITE(MUNIT,1000) REAL(DCORVG(1,IVT))
      DO 101 IVT=1,NVT
101   WRITE(MUNIT,1000) REAL(DCORVG(2,IVT))
      DO 102 IVT=1,NVT
102   WRITE(MUNIT,1000) 0.E0
C
      WRITE(MUNIT,*)'cells ',NEL
      DO  IEL=1,NEL
        WRITE(MUNIT,*)'quad 4'
        WRITE(MUNIT,1100) KVERT(1,IEL),KVERT(2,IEL),
     *                    KVERT(3,IEL),KVERT(4,IEL)
      END DO
C
C      WRITE(MUNIT,*)'materials 1 0'
C      WRITE(MUNIT,*)'mat1'
C      DO 120 IEL=1,NEL
C 120  WRITE(MUNIT,*)'1'
C
      WRITE(MUNIT,*)  'velocity 1'
      DO 130 IVT=1,NVT
130   WRITE(MUNIT,1000) VU(IVT)
      DO 140 IVT=1,NVT
140   WRITE(MUNIT,1000) VV(IVT)
      DO 150 IVT=1,NVT
150   WRITE(MUNIT,1000) 0E0

      WRITE(MUNIT,*)  'variable'

      WRITE(MUNIT,*)  'pressure 1'
      DO 160 IVT=1,NVT
160   WRITE(MUNIT,1000) VP(IVT)

      WRITE(MUNIT,*)  'streamfunction 1'
      DO 170 IVT=1,NVT
170   WRITE(MUNIT,1000) VISO(IVT)

      WRITE(MUNIT,*)  'cell_size 0'
      DO IEL=1,NEL
        WRITE(MUNIT,1000) AREA(IEL)
      END DO

      WRITE(MUNIT,*)  'endvars'

C Export all nodes/cells of moving boundary components with different
C node/cell flags. Using this technique allows the user in the
C postprocessing to filter out all cells/nodes virtually not belonging to
C the geometry.

      WRITE (MUNIT,*) 'flags'

C First export all nodes... (there are inner nodes, boundary nodes and
C nodes belonging to a fict. boundary component)
      
      WRITE (MUNIT,*) 'nodetype 3 1'
      WRITE (MUNIT,*) 'inner boundary movbc'
      DO IVT=1,NVT
        FLAG=1
        IF(KNPR(IVT).GT.0) FLAG=2
        IF (KNPR(IVT).LT.-NEL) FLAG=3
        WRITE (MUNIT,'(I1)') FLAG
      END DO
      
C ... then all cells. All cells with all 4 vertices belonging to a
C fictitious boundary component are treated as inner cells of this
C component.
C Normal cells receive the status "inner", cells of a fict. boundary
C component receive the status "movbc"
      
      WRITE (MUNIT,*) 'celltype 2 0'
      WRITE (MUNIT,*) 'inner movbc'
      DO IEL=1,NEL
        FLAG=2
        DO IVT=1,4
          IF (KNPR(KVERT(IVT,IEL)).GE.-NEL) FLAG=1
        END DO
        WRITE (MUNIT,'(I1)') FLAG
      END DO

C Flag-section ended...
      WRITE (MUNIT,*) 'endflags'
      
C Write out the moving boundary data

      WRITE(MUNIT,*)  'polygons'

C Write analytical moving boundaries

      CALL FBDGMV (MUNIT,3)

      WRITE(MUNIT,*)  'endpoly'

C Write out statistical data      
      
      WRITE(MUNIT,*)  'probtime ',REAL(TIMENS)
      
C That's all...
      
      WRITE(MUNIT,*)  'endgmv'
      REWIND (MUNIT)
      CLOSE  (MUNIT)

1000  FORMAT(E15.8)
1100  FORMAT(8I8)

      END

************************************************************************** 
* Write scalar solution to GMV-file
*
* This routine writes a scalar solution DU with gradient vectors
* DUX/DUY into a GMV-file with filename CFILE. The unit number MUNIT
* will be used for output.
*
* In:
*   MUNIT   - integer
*             Unit number to use
*   CFILE   - string character array
*             Filename to use
*   NEL,
*   NVT,
*   KVERT,
*   DCORVG  - usual geometry informaion
*   DU      - array [1..NVT] of double
*             Solution vector in every node
*   DUX     - array [1..NVT] of double
*             X-dericative in every node
*   DUY     - array [1..NVT] of double
*             Y-derivative in every node
*   TAG     - double precision
*             A user-defineable, double precision tag that is written
*             to the GMV-file. The tag arises in the upper-right corner
*             of the GMV-file as time-indicator.
************************************************************************** 

      SUBROUTINE XGMVSC(MUNIT,CFILE,NEL,NVT,KVERT,DCORVG,
     *                  DU,DUX,DUY,TAG)
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cbasictria.inc'

C parameters

      INTEGER MUNIT, NEL, NVT
      CHARACTER CFILE*(*)
      INTEGER KVERT(NNVE,*)
      DOUBLE PRECISION DCORVG(2,*),TAG
      DOUBLE PRECISION DU(*),DUX(*),DUY(*)
      
C externals

      INTEGER NFBDYC
      EXTERNAL NFBDYC

C local variables

      INTEGER IFIRST, IVT, IEL, FLAG

      OPEN (UNIT=MUNIT,FILE=CFILE)
 
      WRITE(MUNIT,'(A)')'gmvinput ascii'

C Write out the nodes
      
      WRITE(MUNIT,*)'nodes ',NVT

      DO 100 IVT=1,NVT
100   WRITE(MUNIT,1000) REAL(DCORVG(1,IVT))
      DO 101 IVT=1,NVT
101   WRITE(MUNIT,1000) REAL(DCORVG(2,IVT))
      DO 102 IVT=1,NVT
102   WRITE(MUNIT,1000) 0.E0

C ...the cells...

      WRITE(MUNIT,*)'cells ',NEL
      DO  IEL=1,NEL
        WRITE(MUNIT,*)'quad 4'
        WRITE(MUNIT,1100) KVERT(1,IEL),KVERT(2,IEL),
     *                    KVERT(3,IEL),KVERT(4,IEL)
      END DO
C
C      WRITE(MUNIT,*)'materials 1 0'
C      WRITE(MUNIT,*)'mat1'
C      DO 120 IEL=1,NEL
C 120  WRITE(MUNIT,*)'1'
C

C ... the velocity vectors ...

      WRITE(MUNIT,*)  'velocity 1'
      DO 130 IVT=1,NVT
130   WRITE(MUNIT,1000) DUX(IVT)
      DO 140 IVT=1,NVT
140   WRITE(MUNIT,1000) DUY(IVT)
      DO 150 IVT=1,NVT
150   WRITE(MUNIT,1000) 0E0

C ... the solution ...

      WRITE(MUNIT,*)  'variable'

      WRITE(MUNIT,*)  'solution 1'
      DO 160 IVT=1,NVT
160   WRITE(MUNIT,1000) DU(IVT)

      WRITE(MUNIT,*)  'endvars'

C Write out the tag
      
      WRITE(MUNIT,*)  'probtime ',REAL(TAG)

C That's all...
      
      WRITE(MUNIT,*)  'endgmv'
      REWIND (MUNIT)
      CLOSE  (MUNIT)

1000  FORMAT(E15.8)
1100  FORMAT(8I8)

      END
