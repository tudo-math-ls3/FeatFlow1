      SUBROUTINE XOWSC1(MFILE,CFILE)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CFILE*(*)
C
      PARAMETER (NNVE=4,NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/
C
      SUB='XOWSC1'
      IF (ICHECK.GE.997) CALL OTRC('XOWSC1 ','05/15/95')
      IER=0
C
C *** Open I/O-file
C
      CALL OF0(MFILE,CFILE,1)
      IF (IER.NE.0) GOTO 99999
C
      WRITE (MFILE,*) 'Coarse mesh'
      WRITE (MFILE,*) 'modified by trigen2d'
      WRITE (MFILE,*) NEL,NVT,NMT,NVE,NBCT
C
      CALL OWSC1(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *           KWORK(L(LVERT)),KWORK(L(LNPR)),KWORK(L(LMM)),MFILE)
      IF (IER.NE.0) GOTO 99999
C
      IF (LADJ.GT.0) THEN
       CALL ZDISP(0,LADJ,'KADJ  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
C
99999 END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE OWSC1(DCORVG,DVBDP,KVBD,KVERT,KNPR,KMM,MFILE)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=4)
      DIMENSION DCORVG(2,*),DVBDP(*)
      DIMENSION KVERT(NNVE,*),KNPR(*),KMM(2,*),KVBD(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='OWSC1'
      IF (ICHECK.GE.997) CALL OTRC('OWSC1  ','05/15/95')
C
      DO 5 IVBD=1,NVBD
      DCORVG(1,KVBD(IVBD))=DVBDP(IVBD)
      DCORVG(2,KVBD(IVBD))=0D0
5     CONTINUE
C
      WRITE (MFILE,*) 'DCORVG'
      DO 10 IVT=1,NVT
10    WRITE (MFILE,*) DCORVG(1,IVT),DCORVG(2,IVT)
      WRITE (MFILE,*) 'KVERT'
      DO 20 IEL=1,NEL
20    WRITE (MFILE,*) KVERT(1,IEL),KVERT(2,IEL),
     *                KVERT(3,IEL),KVERT(4,IEL)
      WRITE (MFILE,*) 'KNPR'
      DO 30 IVT=1,NVT
30    WRITE (MFILE,*) KNPR(IVT)
      WRITE (MFILE,*) 'KMM'
      DO 40 IBCT=1,NBCT
40    WRITE (MFILE,*) KMM(1,IBCT),KMM(2,IBCT)
C
C
C
99999 END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE SBC2P(DCORVG,DVBDP,KVBD,NVBD)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DCORVG(2,*),DVBDP(*),KVBD(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='SBC2P'
      IF (ICHECK.GE.997) CALL OTRC('SBC2P ','05/15/95')
C
      DO 1 IVBD=1,NVBD
      DCORVG(1,KVBD(IVBD))=DVBDP(IVBD)
      DCORVG(2,KVBD(IVBD))=0D0
1     CONTINUE
C
99999 END
