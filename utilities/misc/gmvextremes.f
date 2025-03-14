      PROGRAM GMVEXTREME
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVT=100000)
      CHARACTER CFILE*15,CFILEH*15
      DIMENSION KVERT(4,NNVT),VCORVG(3,NNVT)
      DIMENSION VU(NNVT),VV(NNVT),VW(NNVT),VP(NNVT),VISO(NNVT),VT(NNVT)
C
      VMAXU=-1E10
      VMAXV=-1E10
      VMAXW=-1E10
      VMAXP=-1E10
      VMAXI=-1E10
      VMAXT=-1E10
C
      VMINU=1E10
      VMINV=1E10
      VMINW=1E10
      VMINP=1E10
      VMINI=1E10
      VMINT=1E10
C
C
C
      WRITE(6,*) 'ISTART,ISTOP,IADD : '
      READ(5,*) ISTART,ISTOP,IADD
C
C
      DO 999 ITNS=ISTART,ISTOP,IADD
C
      VMAXUH=-1E10
      VMAXVH=-1E10
      VMAXWH=-1E10
      VMAXPH=-1E10
      VMAXIH=-1E10
      VMAXTH=-1E10
C
      VMINUH=1E10
      VMINVH=1E10
      VMINWH=1E10
      VMINPH=1E10
      VMINIH=1E10
      VMINTH=1E10
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
2     WRITE(CFILE(IFIRST:IFIRST+4),'(A4)') '.gmv'
C
C=======================================================================
C     NVEC = FUNKT.WERTE PRO KNOTEN
C=======================================================================
      NVEC = 5
      MUNIT=50
C
      OPEN (UNIT=MUNIT,FILE=CFILE)
C
C
      READ(MUNIT,*) CFILEH
      READ(MUNIT,*) CFILEH,NVT
C
      DO 100 IVT=1,NVT
100   READ(MUNIT,1000) VCORVG(1,IVT)
      DO 101 IVT=1,NVT
101   READ(MUNIT,1000) VCORVG(2,IVT)
      DO 102 IVT=1,NVT
102   READ(MUNIT,1000) VCORVG(3,IVT)
C
      READ(MUNIT,*) CFILEH,NEL
      DO 110 IEL=1,NEL
      READ(MUNIT,*) CFILEH
110   READ(MUNIT,1100) KVERT(1,IEL),KVERT(2,IEL),
     *                 KVERT(3,IEL),KVERT(4,IEL)
C
C
      READ(MUNIT,*)  CFILEH
      DO 130 IVT=1,NVT
130   READ(MUNIT,1000) VU(IVT)
      DO 140 IVT=1,NVT
140   READ(MUNIT,1000) VV(IVT)
      DO 150 IVT=1,NVT
150   READ(MUNIT,1000) VW(IVT)

      READ(MUNIT,*)  CFILEH

      READ(MUNIT,*)  CFILEH
      DO 160 IVT=1,NVT
160   READ(MUNIT,1000) VP(IVT)

      READ(MUNIT,*)  CFILEH
      DO 170 IVT=1,NVT
170   READ(MUNIT,1000) VISO(IVT)

      READ(MUNIT,*)  CFILEH
      DO 180 IVT=1,NVT
180   READ(MUNIT,1000) VT(IVT)
C
C
      REWIND (MUNIT)
      CLOSE  (MUNIT)
C
      DO 10 IVT=1,NVT
      VMAXUH=MAX(VMAXUH,VU(IVT))
      VMAXVH=MAX(VMAXVH,VV(IVT))
      VMAXWH=MAX(VMAXWH,VW(IVT))
      VMAXPH=MAX(VMAXPH,VP(IVT))
      VMAXIH=MAX(VMAXIH,VISO(IVT))
      VMAXTH=MAX(VMAXTH,VT(IVT))
C
      VMINUH=MIN(VMINUH,VU(IVT))
      VMINVH=MIN(VMINVH,VV(IVT))
      VMINWH=MIN(VMINWH,VW(IVT))
      VMINPH=MIN(VMINPH,VP(IVT))
      VMINIH=MIN(VMINIH,VISO(IVT))
      VMINTH=MIN(VMINTH,VT(IVT))
C
      VMAXU=MAX(VMAXUH,VMAXU)
      VMAXV=MAX(VMAXVH,VMAXV)
      VMAXW=MAX(VMAXWH,VMAXW)
      VMAXP=MAX(VMAXPH,VMAXP)
      VMAXI=MAX(VMAXIH,VMAXI)
      VMAXT=MAX(VMAXTH,VMAXT)
C
      VMINU=MIN(VMINUH,VMINU)
      VMINV=MIN(VMINVH,VMINV)
      VMINW=MIN(VMINWH,VMINW)
      VMINP=MIN(VMINPH,VMINP)
      VMINI=MIN(VMINIH,VMINI)
      VMINT=MIN(VMINTH,VMINT)
C
10    CONTINUE
C
      write(6,*) 
      write(6,*) ITNS,NVT,CFILE
      WRITE(6,*) ' U :',VMINUH,VMAXUH
      WRITE(6,*) ' V :',VMINVH,VMAXVH
      WRITE(6,*) ' W :',VMINWH,VMAXWH
      WRITE(6,*) ' P :',VMINPH,VMAXPH
      WRITE(6,*) ' I :',VMINIH,VMAXIH
      WRITE(6,*) ' T :',VMINTH,VMAXTH
C
      write(6,*) 'TOTAL'
      WRITE(6,*) ' U :',VMINU,VMAXU
      WRITE(6,*) ' V :',VMINV,VMAXV
      WRITE(6,*) ' W :',VMINW,VMAXW
      WRITE(6,*) ' P :',VMINP,VMAXP
      WRITE(6,*) ' I :',VMINI,VMAXI
      WRITE(6,*) ' T :',VMINT,VMAXT
C
999   CONTINUE
C
C
1000  FORMAT(E12.5)
1100  FORMAT(8I8)
C
C
C
      END
