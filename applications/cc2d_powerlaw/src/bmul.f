************************************************************************
      SUBROUTINE BBUILD(KVERT,KMID,KADJ,DCORVG,B1,B2,KCOLB,KLDB,
     *                  NB,NEL,NVT,NMT)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DOUBLE PRECISION B1,B2
C
      PARAMETER (NNVE=4)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*),DCORVG(2,*)
      DIMENSION B1(*),B2(*),KCOLB(*),KLDB(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      SAVE
C
C
C
      ILD=0
C
C   
      DO 10 IEL=1,NEL
C
      DO 20 IVE=1,4
      IADJ=KADJ(IVE,IEL)
      IF ((IADJ.GT.0).AND.(IADJ.LT.IEL)) GOTO 20
C
      IVEN=IVE+1
      IF (IVEN.EQ.5) IVEN=1
C
      IMID=KMID (IVE ,IEL)-NVT
      IVT =KVERT(IVE ,IEL)
      IVTN=KVERT(IVEN,IEL)
C
      PX =DCORVG(1,IVT)
      PY =DCORVG(2,IVT)
      PXN=DCORVG(1,IVTN)
      PYN=DCORVG(2,IVTN)
C
      DN1=-PYN+PY
      DN2= PXN-PX
C
      ILD=ILD+1
      KLDB (IMID)=ILD
      KCOLB(ILD )=IEL
      B1   (ILD )=DN1
      B2   (ILD )=DN2
C
      IF (IADJ.GT.0) THEN
       ILD=ILD+1
       KCOLB(ILD)= IADJ
       B1   (ILD)=-DN1
       B2   (ILD)=-DN2
      ENDIF
C
20    CONTINUE
C
10    CONTINUE
C
      NB         =ILD
      KLDB(NMT+1)=ILD+1
C
      END
C
C
C
************************************************************************
      SUBROUTINE BMUL1 (KVERT,KMID,KADJ,DCORVG,DP,DF1,DF2,
     *                  NEL,NVT,NMT,A1,A2)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*),DCORVG(2,*)
      DIMENSION DP(*),DF1(*),DF2(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      SAVE
C
C
C
      DO 10 IEL=1,NEL
C
      DO 20 IVE=1,4
      IADJ=KADJ(IVE,IEL)
      IF ((IADJ.GT.0).AND.(IADJ.LT.IEL)) GOTO 20
C
      IVEN=IVE+1
      IF (IVEN.EQ.5) IVEN=1
C
      IMID=KMID (IVE ,IEL)-NVT
      IVT =KVERT(IVE ,IEL)
      IVTN=KVERT(IVEN,IEL)
C
      PX =DCORVG(1,IVT)
      PY =DCORVG(2,IVT)
      PXN=DCORVG(1,IVTN)
      PYN=DCORVG(2,IVTN)
C
      DN1=-PYN+PY
      DN2= PXN-PX
C
      IF (IADJ.GT.0) THEN
       DF1(IMID)=A2*DF1(IMID)+A1*DN1*(DP(IEL)-DP(IADJ))
       DF2(IMID)=A2*DF2(IMID)+A1*DN2*(DP(IEL)-DP(IADJ))
      ELSE
       DF1(IMID)=A2*DF1(IMID)+A1*DN1*DP(IEL)
       DF2(IMID)=A2*DF2(IMID)+A1*DN2*DP(IEL)
      ENDIF
C
20    CONTINUE
C
10    CONTINUE
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE BTMUL1(KVERT,KMID,KADJ,DCORVG,DU1,DU2,DFP,
     *                  NEL,NVT,NMT,A1)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*),DCORVG(2,*)
      DIMENSION DFP(*),DU1(*),DU2(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      SAVE
C
C
C
      DO 10 IEL=1,NEL
      DFP(IEL)=0D0
C
      DO 20 IVE=1,4
      IVEN=IVE+1
      IF (IVEN.EQ.5) IVEN=1
C
      IMID=KMID (IVE ,IEL)-NVT
      IVT =KVERT(IVE ,IEL)
      IVTN=KVERT(IVEN,IEL)
C
      PX =DCORVG(1,IVT)
      PY =DCORVG(2,IVT)
      PXN=DCORVG(1,IVTN)
      PYN=DCORVG(2,IVTN)
C
      DN1=-PYN+PY
      DN2= PXN-PX
C
      DFP(IEL)=DFP(IEL)+A1*(DN1*DU1(IMID)+DN2*DU2(IMID))
C
20    CONTINUE
C
10    CONTINUE
C
C
C
      END
