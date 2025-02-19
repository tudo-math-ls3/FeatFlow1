C *******************************************************************
      SUBROUTINE TRCORE(KVERT,DCORVG,DCORE,NEL)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KVERT(4,*),DCORVG(2,*),DCORE(2,*)
C
      DO 10 IEL=1,NEL
      DCORE(1,IEL)=0D0
      DCORE(2,IEL)=0D0
C
      DO 20 IVE=1,4
      IVT=KVERT(IVE,IEL)
      PX=DCORVG(1,IVT)
      PY=DCORVG(2,IVT)
C
      DCORE(1,IEL)=DCORE(1,IEL)+0.25D0*PX
      DCORE(2,IEL)=DCORE(2,IEL)+0.25D0*PY
C
20    CONTINUE
10    CONTINUE
C
      END
C
C *******************************************************************
      SUBROUTINE TRCORM(KVERT,KMID,DCORVG,DCORM,NEL,NVT)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KVERT(4,*),KMID(4,*),DCORVG(2,*),DCORM(2,*)
C
      NVE=4
      NMT=NVT
C
      DO 10 IEL=1,NEL
      DO 20 IVE=1,NVE
C
      IMID=KMID(IVE,IEL)
      IF (IMID.LE.NMT) GOTO 20
C
      NMT =NMT+1
      IVEH=IVE+1
      IF (IVEH.EQ.NVE+1) IVEH=1
C
      IVT1=KVERT(IVE,IEL)
      PX1 =DCORVG(1,IVT1)
      PY1 =DCORVG(2,IVT1)
C
      IVT2=KVERT(IVEH,IEL)
      PX2 =DCORVG(1,IVT2)
      PY2 =DCORVG(2,IVT2)
C
      PX=0.5D0*(PX1+PX2)
      PY=0.5D0*(PY1+PY2)
      DCORM(1,IMID-NVT)=PX
      DCORM(2,IMID-NVT)=PY
C
20    CONTINUE
10    CONTINUE
C
      END
C
C *******************************************************************
C
      SUBROUTINE TRSRT(KV1,KV2,INUM,N,DCORVG)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (EPS=1D-12)
      DIMENSION KV1(*),KV2(*),DCORVG(2,*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      SAVE /OUTPUT/
C
C *** Initialization
      DO 1 I=1,N
1     KV1(I)=I
C
      IF (INUM.EQ.1) INU=2
      IF (INUM.EQ.2) INU=1
C
C *** Sorting by INUM-components first (Heap-Sort)
      CALL SORT(KV1,N,DCORVG,INUM)
C
C *** Partially sorting by INU-components then (Heap-Sort)
      I1=1
      Y1=DCORVG(INUM,KV1(1))
C
      DO 2 I=2,N
      I2=I
      Y2=DCORVG(INUM,KV1(I))
      IF (ABS(Y1-Y2).GT.EPS) THEN
       N1=I2-I1
       IF (N1.NE.1) THEN
       CALL SORT(KV1(I1),N1,DCORVG,INU)
       ENDIF
       I1=I2
       Y1=Y2
      ELSEIF (I2.EQ.N) THEN
       N1=I2-I1+1
       IF (N1.NE.1) THEN
       CALL SORT(KV1(I1),N1,DCORVG,INU)
       ENDIF
      ENDIF
2     CONTINUE
C
      CALL RANK(KV1,KV2,N)
C
      END
C
C *******************************************************************
C	
      SUBROUTINE SORT(KV1,N,DCORVG,IDIM)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KV1(*),DCORVG(2,*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      SAVE /OUTPUT/
C
      L=N/2+1
      IR=N
C
10    CONTINUE
      IF (L.GT.1) THEN
       L=L-1
       KKV1=KV1(L)
      ELSE
       KKV1=KV1(IR)
       KV1(IR)=KV1(1)
       IR=IR-1
       IF (IR.EQ.1) THEN
        KV1(1)=KKV1
        RETURN
       ENDIF
      ENDIF
      I=L
      J=L+L
20    IF (J.LE.IR) THEN
       IF (J.LT.IR) THEN
        IF (DCORVG(IDIM,KV1(J)).LT.DCORVG(IDIM,KV1(J+1))) J=J+1
       ENDIF
       IF (DCORVG(IDIM,KKV1).LT.DCORVG(IDIM,KV1(J))) THEN
        KV1(I)=KV1(J)
        I=J
        J=J+J
       ELSE
        J=IR+1
       ENDIF
       GOTO 20
      ENDIF
      KV1(I)=KKV1
      GOTO 10
C
      END
C	
C *******************************************************************
C
      SUBROUTINE RANK(KV1,KV2,N)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KV1(*),KV2(*)
C
      DO 1 J=1,N
1     KV2(KV1(J))=J
C
      END
C	
C *******************************************************************
C
      SUBROUTINE MTSRTD(DA,DAH,KCOL,KCOLH,KLD,KLDH,KTR1,KTR2,NEQ)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),DAH(*),KCOL(*),KCOLH(*),KLD(*),KLDH(*),
     *          KTR1(*),KTR2(*)
      DIMENSION KH1(100),KH2(100)
C
C
C
      ILD=1
C
      DO 10 I=1,NEQ
C
      I1       =KTR1(I)
      DA(ILD)  =DAH(KLDH(I1))
      KLD(I)   =ILD
      KCOL(ILD)=I
      ILD      =ILD+1
C
      IH1=KLDH(I1)+1
      IH2=KLDH(I1+1)-1
      ID1=IH2-IH1+1
C
      DO 11 JH=IH1,IH2
      KH1(JH-IH1+1)=KTR2(KCOLH(JH))
      KH2(JH-IH1+1)=JH
11    CONTINUE
C
      CALL KHSORT(KH1,KH2,ID1)
C
      DO 20 J=1,ID1
      DH       =DAH(KH2(J))
      ICOL     =KH1(J)
      DA(ILD)  =DH
      KCOL(ILD)=ICOL
      ILD      =ILD+1
20    CONTINUE
C
10    CONTINUE
C
      KLD(NEQ+1)=KLDH(NEQ+1)
C
      END
C
C *******************************************************************
C
      SUBROUTINE MTSRTV(VA,VAH,KCOL,KCOLH,KLD,KLDH,KTR1,KTR2,NEQ)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),VAH(*),KCOL(*),KCOLH(*),KLD(*),KLDH(*),
     *          KTR1(*),KTR2(*)
      DIMENSION KH1(100),KH2(100)
C
C
C
      ILD=1
C
      DO 10 I=1,NEQ
C
      I1       =KTR1(I)
      VA(ILD)  =VAH(KLDH(I1))
      KLD(I)   =ILD
      KCOL(ILD)=I
      ILD      =ILD+1
C
      IH1=KLDH(I1)+1
      IH2=KLDH(I1+1)-1
      ID1=IH2-IH1+1
C
      DO 11 JH=IH1,IH2
      KH1(JH-IH1+1)=KTR2(KCOLH(JH))
      KH2(JH-IH1+1)=JH
11    CONTINUE
C
      CALL KHSORT(KH1,KH2,ID1)
C
      DO 20 J=1,ID1
      VH       =VAH(KH2(J))
      ICOL     =KH1(J)
      VA(ILD)  =VH
      KCOL(ILD)=ICOL
      ILD      =ILD+1
20    CONTINUE
C
10    CONTINUE
C
      KLD(NEQ+1)=KLDH(NEQ+1)
C
      END
C
C *******************************************************************
C
      SUBROUTINE MTSRTR(KCOL,KCOLH,KLD,KLDH,KTR1,KTR2,NEQ)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KCOL(*),KCOLH(*),KLD(*),KLDH(*),KTR1(*),KTR2(*)
      DIMENSION KH1(100),KH2(100)
C
C
C
      ILD=1
C
      DO 10 I=1,NEQ
C
      I1       =KTR2(I)
      KLD(I)   =ILD
      KCOL(ILD)=I
      ILD      =ILD+1
C
      IH1=KLDH(I1)+1
      IH2=KLDH(I1+1)-1
      ID1=IH2-IH1+1
C
      DO 11 JH=IH1,IH2
      KH1(JH-IH1+1)=KTR1(KCOLH(JH))
      KH2(JH-IH1+1)=JH
11    CONTINUE
C
      CALL KHSORT(KH1,KH2,ID1)
C
      DO 20 J=1,ID1
      ICOL     =KH1(J)
      KCOL(ILD)=ICOL
      ILD      =ILD+1
20    CONTINUE
C
10    CONTINUE
C
      KLD(NEQ+1)=KLDH(NEQ+1)
C
      END
C
C *******************************************************************
C
      SUBROUTINE KHSORT(KH1,KH2,IDIM)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KH1(*),KH2(*)
C
      BMORE=.TRUE.
5     IF (.NOT.BMORE) GOTO 99999
      BMORE=.FALSE.
C
      DO 10 ICOMP=1,IDIM-1
C
      IF (KH1(ICOMP).GT.KH1(ICOMP+1)) THEN
       IAUX1=KH1(ICOMP)
       IAUX2=KH2(ICOMP)
       KH1(ICOMP)=KH1(ICOMP+1)
       KH2(ICOMP)=KH2(ICOMP+1)
       KH1(ICOMP+1)=IAUX1
       KH2(ICOMP+1)=IAUX2
       BMORE=.TRUE.
      ENDIF
C
10    CONTINUE
C
      GOTO 5
C
99999 END
C
C *******************************************************************
C
      SUBROUTINE VECSRT(DX,DD,KTR1,KTR2,NEQ,IPAR)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DX(*),DD(*),KTR1(*),KTR2(*)
C
      IF (IPAR.EQ.1) THEN
       DO 10 IEQ=1,NEQ
10     DD(IEQ)=DX(KTR1(IEQ))
      ELSE
       DO 20 IEQ=1,NEQ
20     DD(IEQ)=DX(KTR2(IEQ))
      ENDIF
C
C
C
      END
C
C *******************************************************************
C
      SUBROUTINE CUTCE0(KLD,KCOL,KCON,KDEG,NEQ,NDEG)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KLD(*),KCOL(*),KCON(*),KDEG(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      SAVE
C
C
C
       DO 10 IEQ=1,NEQ
C       
       DO 20 IDEG=1,NDEG
20     KDEG(IDEG)=0
C
       DO 22 ILD=KLD(IEQ)+1,KLD(IEQ+1)-1
22     KDEG(ILD-KLD(IEQ))=KCOL(ILD)
C
       DO 30 IDEG1=1,KLD(IEQ+1)-KLD(IEQ)-1
       INDMIN=1
C
       IF (KDEG(INDMIN).EQ.0) THEN
        MINDEG=NEQ
       ELSE
        MINDEG=KLD(KDEG(INDMIN)+1)-KLD(KDEG(INDMIN))-1
       ENDIF
C       
       DO 40 IDEG2=1,KLD(IEQ+1)-KLD(IEQ)-1
       IF (KDEG(IDEG2).EQ.0) THEN
        INDDEG=NEQ
       ELSE
        INDDEG=KLD(KDEG(IDEG2)+1)-KLD(KDEG(IDEG2))-1
       ENDIF
C
       IF (INDDEG.LT.MINDEG) THEN
        INDMIN=IDEG2
        MINDEG=INDDEG
       ENDIF
C       
40     CONTINUE
C
       KCON(KLD(IEQ)+IDEG1)=KDEG(INDMIN)
       KDEG(INDMIN)        =0
C
30     CONTINUE
C
10     CONTINUE
C
C
C
99999 END      
C
C *******************************************************************
C
      SUBROUTINE CUTCE1(KLD,KCON,NEQ,KTR1,KTR2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KLD(*),KCON(*),KTR1(*),KTR2(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      SAVE
C
      ICOUNT=1
      INDNEQ=1
      KTR1(1)=1
      KTR2(1)=1
C
1     DO 10 ILD=KLD(KTR1(ICOUNT)),KLD(KTR1(ICOUNT)+1)-1
      ICOL=KCON(ILD)
      IF (KTR2(ICOL).NE.0) GOTO 10
C
      INDNEQ      =INDNEQ+1
      KTR2(ICOL  )=INDNEQ
      KTR1(INDNEQ)=ICOL
10    CONTINUE
C
      ICOUNT=ICOUNT+1
      IF (ICOUNT.NE.NEQ) GOTO 1
C
99999 END      
