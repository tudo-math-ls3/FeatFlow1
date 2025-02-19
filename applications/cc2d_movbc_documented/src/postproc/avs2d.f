      SUBROUTINE XAVS2D(MUNIT,CFILE,NEL,NVT,KVERT,DCORVG,VU,VV,VP,VISO)

      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'

C parameters

      INTEGER MUNIT, NEL, NVT
      CHARACTER CFILE*(*)
      INTEGER KVERT(4,*)
      DOUBLE PRECISION DCORVG(2,*)
      REAL VU(*),VV(*),VP(*),VISO(*)

C local variables

      INTEGER IFIRST, NVEC, IVT, IEL

      DO 1 IFIRST = 1,60
      IF (CFILE(IFIRST:IFIRST).EQ.' ') GOTO 2
1     CONTINUE
2     WRITE(CFILE(IFIRST:IFIRST+4),'(A4)') '.inp'
C
C=======================================================================
C     NVEC = FUNKT.WERTE PRO KNOTEN
C=======================================================================
      NVEC = 5
C
      OPEN (UNIT=MUNIT,FILE=CFILE)
C
C
      WRITE(MUNIT,*) NVT,NEL,NVEC,0,0
C
      DO 100 IVT=1,NVT
100   WRITE(MUNIT,1000) IVT,REAL(DCORVG(1,IVT)),REAL(DCORVG(2,IVT)),0E0
C
      DO 110 IEL=1,NEL
110   WRITE(MUNIT,1100) IEL,1,'quad',KVERT(1,IEL),KVERT(2,IEL),
     *                               KVERT(3,IEL),KVERT(4,IEL)
C
C
      WRITE(MUNIT,*) 3,3,1,1
      WRITE(MUNIT,*)  'vel'//',','m/s'
      WRITE(MUNIT,*)  'p'  //',','m/s'
      WRITE(MUNIT,*)  'stream'//',','m/s'      
C
      DO 200 IVT=1,NVT
200   WRITE(MUNIT,2000) IVT,VU(IVT),VV(IVT),0E0,VP(IVT),VISO(IVT)
C
C
      REWIND (MUNIT)
      CLOSE  (MUNIT)
C
C
1000  FORMAT(I8,3E15.8)
1100  FORMAT(2I8,A7,8I8)
2000  FORMAT(I8,5E15.8)
C
      END
