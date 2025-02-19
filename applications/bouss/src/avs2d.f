      SUBROUTINE XGMV2D(MUNIT,CFILE,NEL,NVT,KVERT,DCORVG,VU,VV,VP,
     *                  VISO,VT,TIMENS)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CFILE*(*)
      DIMENSION KVERT(4,*),DCORVG(2,*)
      DIMENSION VU(*),VV(*),VP(*),VISO(*),VT(*)
C
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      SAVE
C
C
C
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
      DO 110 IEL=1,NEL
      WRITE(MUNIT,*)'quad 4'
110   WRITE(MUNIT,1100) KVERT(1,IEL),KVERT(2,IEL),
     *                  KVERT(3,IEL),KVERT(4,IEL)
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

      WRITE(MUNIT,*)  'temperature 1'
      DO 180 IVT=1,NVT
180   WRITE(MUNIT,1000) VT(IVT)

      WRITE(MUNIT,*)  'endvars'
      WRITE(MUNIT,*)  'probtime ',REAL(TIMENS)
      WRITE(MUNIT,*)  'endgmv'
C
C
C
      REWIND (MUNIT)
      CLOSE  (MUNIT)
C
C
1000  FORMAT(E15.8)
1100  FORMAT(8I8)
C
      END
C
C
C
      SUBROUTINE XAVS2D(MUNIT,CFILE,NEL,NVT,KVERT,DCORVG,
     *                   VU,VV,VP,VISO,VT)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CFILE*(*)
      DIMENSION KVERT(4,*),DCORVG(2,*)
      DIMENSION VU(*),VV(*),VP(*),VISO(*),VT(*)
C
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      SAVE
C
C
C
      DO 1 IFIRST = 1,60
      IF (CFILE(IFIRST:IFIRST).EQ.' ') GOTO 2
1     CONTINUE
2     WRITE(CFILE(IFIRST:IFIRST+4),'(A4)') '.inp'
C
C=======================================================================
C     NVEC = FUNKT.WERTE PRO KNOTEN
C=======================================================================
      NVEC = 6
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
      WRITE(MUNIT,*) 4,3,1,1,1
      WRITE(MUNIT,*)  'vel'//',','m/s'
      WRITE(MUNIT,*)  'p'  //',','m/s'
      WRITE(MUNIT,*)  'stream'//',','m/s'      
      WRITE(MUNIT,*)  'Temp'//',','K'      
C
      DO 200 IVT=1,NVT
200   WRITE(MUNIT,2000) IVT,VU(IVT),VV(IVT),0E0,VP(IVT),
     *                  VISO(IVT),VT(IVT)
C
C
      REWIND (MUNIT)
      CLOSE  (MUNIT)
C
C
1000  FORMAT(I8,3E15.8)
1100  FORMAT(2I8,A7,8I8)
2000  FORMAT(I8,6E15.8)
      END
