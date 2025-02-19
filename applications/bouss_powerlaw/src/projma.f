************************************************************************
      SUBROUTINE   PROJST  (KCOLC,KLDC,KADJ,NEL,NC)
************************************************************************
*    Purpose:  calculates the matrix positions for projection of C
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=4)
      DIMENSION KCOLC(*),KLDC(*),KADJ(NNVE,*)
C
C
C
      NC=0
      KLDC(1)=1
C
      DO 10 IEL=1,NEL
      NC=NC+1
      KCOLC(NC)=IEL
C
      DO 20 IVE=1,4
      IADJ=KADJ(IVE,IEL)
      IF (IADJ.EQ.0) GOTO 20
C
      NC=NC+1
      KCOLC(NC)=IADJ
20    CONTINUE
C
      KLDC(IEL+1)=NC+1
C
10    CONTINUE
C
C
      DO 30 IEQ=1,NEL
C
31    BSORT=.TRUE.
      DO 32 ICOL=KLDC(IEQ)+1,KLDC(IEQ+1)-2
      IF (KCOLC(ICOL).GT.KCOLC(ICOL+1)) THEN
       IHELP=KCOLC(ICOL)
       KCOLC(ICOL)=KCOLC(ICOL+1)
       KCOLC(ICOL+1)=IHELP
       BSORT=.FALSE.
      ENDIF
32    CONTINUE
      IF (.NOT.BSORT) GOTO 31
C
30    CONTINUE
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE   PROJMA  (DC,KCOLC,KLDC,KNPR,KMID,KADJ,
     *                      DM,B1,B2,KCOLB,KLDB,NEL,NVT,NMT)
************************************************************************
*    Purpose:  calculates the matrix entries for projection of C
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DOUBLE PRECISION B1,B2
C
      PARAMETER (NNVE=4)
      DIMENSION DC(*),KCOLC(*),KLDC(*)
      DIMENSION KNPR(*),KMID(NNVE,*),KADJ(NNVE,*)
      DIMENSION DM(*),B1(*),B2(*),KCOLB(*),KLDB(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
C
C
C
      DO 10 IEL=1,NEL
C
      DO 20 IVE=1,4
      IMID=KMID(IVE,IEL)-NVT
      IADJ=KADJ(IVE,IEL)
      INPR=KNPR(IMID+NVT)
      IF (INPR.NE.0) GOTO 20
C
      ILD1=KLDC(IEL)
      ILD2=KLDC(IEL+1)-1
C
      DO 30 ILDB1=KLDB(IMID),KLDB(IMID+1)-1
      IF (KCOLB(ILDB1).EQ.IEL) GOTO 32
30    CONTINUE
C
32    CONTINUE 
C
      DH1= B1(ILDB1)**2+B2(ILDB1)**2
      DC(ILD1)=DC(ILD1)+DH1/DM(IMID)
C
      IF (IADJ.EQ.0) GOTO 20
C
      DO 40 ILDB2=KLDB(IMID),KLDB(IMID+1)-1
      IF (KCOLB(ILDB2).EQ.IADJ) GOTO 42
40    CONTINUE
C
42    CONTINUE 
C
      DO 50 ILDC=ILD1+1,ILD2
      IF (KCOLC(ILDC).EQ.IADJ) GOTO 52
50    CONTINUE
C
52    CONTINUE
C
C
      DH2= B1(ILDB1)*B1(ILDB2)+B2(ILDB1)*B2(ILDB2)
      DC(ILDC)=DH2/DM(IMID)
C
C
20    CONTINUE
C
10    CONTINUE
C
C
      IF (MT.GT.9) THEN
      DO 200 IEQ=1,NEL
      DO 210 ILD=KLDC(IEQ),KLDC(IEQ+1)-1
210   WRITE(6,*) IEQ,ILD,KCOLC(ILD),DC(ILD)
      WRITE(6,*) '***'
200   CONTINUE
      ENDIF
C
      END
