      PROGRAM AVS2GMV2D
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVT=200000)
      CHARACTER CFILE*15,CFILEH*15
      DIMENSION KVERT(4,NNVT),VCORVG(3,NNVT)
      DIMENSION VU(NNVT),VV(NNVT),VW(NNVT),VP(NNVT),VISO(NNVT),VT(NNVT)
C
C
C
      WRITE(6,*) 'ISTART,ISTOP,IADD : '
      READ(5,*) ISTART,ISTOP,IADD
C
C
      DO 999 ITNS=ISTART,ISTOP,IADD
C
      CFILE='u.        '
      IF ((ITNS.GE.0).AND.(ITNS.LT.10)) 
     *     WRITE(CFILE(3:3),'(I1.1)') ITNS
      IF ((ITNS.GE.10).AND.(ITNS.LT.100)) 
     *     WRITE(CFILE(3:4),'(I2.2)') ITNS
      IF ((ITNS.GE.100).AND.(ITNS.LT.1000)) 
     *     WRITE(CFILE(3:5),'(I3.3)') ITNS
      IF ((ITNS.GE.1000).AND.(ITNS.LT.10000)) 
     *     WRITE(CFILE(3:6),'(I4.4)') ITNS
      IF (ITNS.GE.10000) STOP
C
      DO 1 IFIRST = 1,60
      IF (CFILE(IFIRST:IFIRST).EQ.' ') GOTO 2
1     CONTINUE
2     WRITE(CFILE(IFIRST:IFIRST+4),'(A4)') '.inp'
C
C=======================================================================
C     read AVS file
C=======================================================================
      NVEC = 6
      MUNIT=50
C
      OPEN (UNIT=MUNIT,FILE=CFILE)
C
C
      READ(MUNIT,*) NVT,NEL,NVEC,IH,IH
C
      DO 100 IVT=1,NVT
100   READ(MUNIT,1000) IVTH,VCORVG(1,IVT),VCORVG(2,IVT),VCORVG(3,IVT)
C
      DO 110 IEL=1,NEL
110   READ(MUNIT,1100) IELH,IH,CFILEH,KVERT(1,IEL),KVERT(2,IEL),
     *                                KVERT(3,IEL),KVERT(4,IEL)
C
C
      READ(MUNIT,*) IH1,IH2,IH3,IH4,IH5
      READ(MUNIT,*) CFILEH
      READ(MUNIT,*) CFILEH
      READ(MUNIT,*) CFILEH
      READ(MUNIT,*) CFILEH
C
      DO 120 IVT=1,NVT
120   READ(MUNIT,1200) IVTH,VU(IVT),VV(IVT),VW(IVT),VP(IVT),
     *                      VISO(IVT),VT(IVT)
C
      REWIND (MUNIT)
      CLOSE  (MUNIT)
C
C=======================================================================
C     finish AVS file + prepare for GMV file
C=======================================================================
C
      CFILE='u.        '
      IF ((ITNS.GE.0).AND.(ITNS.LT.10)) 
     *     WRITE(CFILE(3:3),'(I1.1)') ITNS
      IF ((ITNS.GE.10).AND.(ITNS.LT.100)) 
     *     WRITE(CFILE(3:4),'(I2.2)') ITNS
      IF ((ITNS.GE.100).AND.(ITNS.LT.1000)) 
     *     WRITE(CFILE(3:5),'(I3.3)') ITNS
      IF ((ITNS.GE.1000).AND.(ITNS.LT.10000)) 
     *     WRITE(CFILE(3:6),'(I4.4)') ITNS
      IF (ITNS.GE.10000) STOP
C
      DO 198 IFIRST = 1,60
      IF (CFILE(IFIRST:IFIRST).EQ.' ') GOTO 199
198   CONTINUE
199   WRITE(CFILE(IFIRST:IFIRST+4),'(A4)') '.gmv'
C
C=======================================================================
C     write GMV file
C=======================================================================
      NVEC = 5
      MUNIT=51
C
      OPEN (UNIT=MUNIT,FILE=CFILE)
C
C
      WRITE(MUNIT,'(A)')'gmvinput ascii'
      WRITE(MUNIT,*)'nodes ',NVT
C
      DO 200 IVT=1,NVT
200   WRITE(MUNIT,2000) VCORVG(1,IVT)
      DO 201 IVT=1,NVT
201   WRITE(MUNIT,2000) VCORVG(2,IVT)
      DO 202 IVT=1,NVT
202   WRITE(MUNIT,2000) VCORVG(3,IVT)
C
      WRITE(MUNIT,*)'cells ',NEL
      DO 210 IEL=1,NEL
      WRITE(MUNIT,*)'quad 4'
210   WRITE(MUNIT,2100) KVERT(1,IEL),KVERT(2,IEL),
     *                  KVERT(3,IEL),KVERT(4,IEL)
C
C
      WRITE(MUNIT,*)  'velocity 1'
      DO 230 IVT=1,NVT
230   WRITE(MUNIT,2000) VU(IVT)
      DO 240 IVT=1,NVT
240   WRITE(MUNIT,2000) VV(IVT)
      DO 250 IVT=1,NVT
250   WRITE(MUNIT,2000) VW(IVT)

      WRITE(MUNIT,*)  'variable'

      WRITE(MUNIT,*)  'pressure 1'
      DO 260 IVT=1,NVT
260   WRITE(MUNIT,2000) VP(IVT)

      WRITE(MUNIT,*)  'streamfunction 1'
      DO 270 IVT=1,NVT
270   WRITE(MUNIT,2000) VISO(IVT)

      WRITE(MUNIT,*)  'temperature 1'
      DO 280 IVT=1,NVT
280   WRITE(MUNIT,2000) VT(IVT)

      WRITE(MUNIT,*)  'endvars'
      WRITE(MUNIT,*)  'endgmv'
C
C
      REWIND (MUNIT)
      CLOSE  (MUNIT)
C
C
999   CONTINUE
C
C
1000  FORMAT(I8,3E12.5)
1100  FORMAT(2I8,A7,8I8)
1200  FORMAT(I8,6E12.5)
C
2000  FORMAT(E12.5)
2100  FORMAT(8I8)
C
C
C
      END
