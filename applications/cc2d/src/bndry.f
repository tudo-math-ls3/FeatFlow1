************************************************************************
      SUBROUTINE    BDRNEU  (KMBD,KVBD,KEBD,KVERT,KMID,KNPR,DDBD,DMBDP,
     *                       NMBD,INEUM)
************************************************************************
*    Purpose:  sets the DIRICHLET- and NEUMANN components
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=4)
      DIMENSION KMBD(*),KVBD(*),KEBD(*),KVERT(NNVE,*),KMID(NNVE,*)
      DIMENSION KNPR(*),DDBD(*),DMBDP(*)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      SAVE
C
      NMBD =0
      INEUM=0
C
      DO 1 IVBD=1,NVBD
      CALL GETMBD(IMID,IV1,IV2,IVBD,KVBD,KEBD,KVERT,KMID,KNPR,INPR)
      IF (IMID.EQ.0) GOTO 1
C
      DPAR=DMBDP(IVBD)
      NMBD=NMBD+1
      KMBD(NMBD)=IMID
      DDBD(NMBD)=DPAR
C
      INPART=0
      CALL NEUDAT(INPART,INPRN,DPARN1,DPARN2,TIMENS)
      NPART=INPART
C      
      DO 10 INPART=1,NPART
      CALL NEUDAT(INPART,INPRN,DPARN1,DPARN2,TIMENS)
      IF ((DPAR.GT.DPARN1).AND.(DPAR.LT.DPARN2)
     *                    .AND.(INPR.EQ.INPRN)) THEN
       INEUM=1
       KNPR(NVT+IMID)=0
       KMBD(NMBD)=-KMBD(NMBD)
      ENDIF
10    CONTINUE
C
1     CONTINUE
C
      END
C
************************************************************************
      SUBROUTINE PDSET (KVBD,KEBD,KVERT,KMID,KNPR,DCORVG,DMBDP,DF1,DF2,
     *                  TSTEPB)
************************************************************************
*    Purpose:  sets the PRESSURE DROP values
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      DIMENSION KVBD(*),KEBD(*),KVERT(NNVE,*),KMID(NNVE,*),KNPR(*)
      DIMENSION DCORVG(2,*),DMBDP(*),DF1(*),DF2(*)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
      SAVE
C
C
      DO 10 IVBD=1,NVBD
      CALL GETMBD(IMID,IV1,IV2,IVBD,KVBD,KEBD,KVERT,KMID,KNPR,INPR)
      IF (IMID.EQ.0) GOTO 10
C
      DPAR=DMBDP(IVBD)
C
      INPART=0
      CALL NEUDAT(INPART,INPRN,DPARN1,DPARN2,TIMENS)
      NPART=INPART
C      
      DO 20 INPART=1,NPART
      CALL NEUDAT(INPART,INPRN,DPARN1,DPARN2,TIMENS)
      IF ((DPAR.GT.DPARN1).AND.(DPAR.LT.DPARN2)
     *                    .AND.(INPR.EQ.INPRN)) THEN
       PX1  =DCORVG(1,IV1)
       PY1  =DCORVG(2,IV1)
       PX2  =DCORVG(1,IV2)
       PY2  =DCORVG(2,IV2)
       DN1  =-PY2+PY1
       DN2  = PX2-PX1
       PMEAN= FDATIN(7,INPR,DPAR,DPAR,TIMENS,RE)
       DF1(IMID)=DF1(IMID)+PMEAN*DN1*TSTEPB
       DF2(IMID)=DF2(IMID)+PMEAN*DN2*TSTEPB
      ENDIF
20    CONTINUE
C
10    CONTINUE
      END
C
************************************************************************
      SUBROUTINE    BDRSET  (DU1,DU2,DF1,DF2,KMBD,DDBD,KNPR,NMBD,NVT,
     *                       PARX,PARY,UE)
************************************************************************
*    Purpose:  updates the solution vector (DU1,DU2) and the right hand 
*              side (DF1,DF2) for all DIRICHLET boundary nodes
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      DIMENSION DU1(*),DU2(*),DF1(*),DF2(*),DDBD(*)
      DIMENSION KMBD(*),KNPR(*)
      EXTERNAL PARX,PARY,UE
C
C
      DO 1 IMBD=1,NMBD
      IMID=KMBD(IMBD)
      IF (IMID.LT.0) GOTO 1
      INPR=KNPR(KNPR(NVT+IMID))
      DPAR=DDBD(IMBD)
C
      PX=PARX(DPAR,INPR)
      PY=PARY(DPAR,INPR)
C
      U1=UE(PX,PY,1)
      U2=UE(PX,PY,2)
      DU1(IMID)=U1
      DU2(IMID)=U2
      DF1(IMID)=U1
      DF2(IMID)=U2
1     CONTINUE
      END
C
************************************************************************
      SUBROUTINE    BDRDEF  (DX,KMBD,NMBD,A1)
************************************************************************
*    Purpose:  sets the NEUMANN-components of the vector DX to A1
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DX(*),KMBD(*)
C
C *** loop over all boundary vertices
C
      DO 1 IMBD=1,NMBD
      IMID=KMBD(IMBD)
      IF (IMID.GT.0) GOTO 1
      DX(-IMID)=A1*DX(-IMID)
1     CONTINUE
      END
C
************************************************************************
      SUBROUTINE   BDRYA  (VA,KCOL,KLD,KMBD,NMBD)
************************************************************************
*    Purpose:  updates the matrix entries for all DIRICHLET boundary
*              nodes
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL VA
C
      PARAMETER (NNVE=4)
      DIMENSION VA(*),KCOL(*),KLD(*),KMBD(*)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      SAVE 
C
      DO 1 IMBD=1,NMBD
      IMID=KMBD(IMBD)
      IF (IMID.LT.0) GOTO 1
C
C *** The diagonal element is set to 1
      VA(KLD(IMID))=1E0
C
      DO 2 ICOL=KLD(IMID)+1,KLD(IMID+1)-1
2     VA(ICOL)=0E0
C
1     CONTINUE
      END
C
************************************************************************
      SUBROUTINE    BDRY0  (D1,D2,KMBD,NMBD)
************************************************************************
*    Purpose:  sets the DIRICHLET-components of the vector (D1,D2) to
*              zero
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      DIMENSION D1(*),D2(*)
      DIMENSION KMBD(*)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      SAVE
C
C
      DO 1 IMBD=1,NMBD
      IMID=KMBD(IMBD)
      IF (IMID.LT.0) GOTO 1
      D1(IMID)=0D0
      D2(IMID)=0D0
1     CONTINUE
C
      END
C
************************************************************************
      SUBROUTINE GETMBD (IMID,IVT1,IVT2,IVBD,KVBD,KEBD,KVERT,KMID,KNPR,
     *                   INPR)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      DIMENSION KVBD(*),KEBD(*),KVERT(NNVE,*),KMID(NNVE,*),KNPR(*)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      SAVE
C
      IVT1=KVBD(IVBD)
      INPR=KNPR(IVT1)
      IEL=KEBD(IVBD)
      IF (IEL.EQ.0) THEN
        IMID=0
        GOTO 99999
      ENDIF
C
      DO 1  I=1,4
      IVERT=I
      IF (KVERT(I,IEL).EQ.IVT1) GOTO 2
   1  CONTINUE
C
C *** Error
      WRITE(MTERM,*) 'ERROR in GETMBD: vertice not found'
      RETURN
C
   2  IMID=KMID(IVERT,IEL)-NVT
      IVERT2=IVERT+1
      IF (IVERT2.GT.NNVE) IVERT2=1
      IVT2=KVERT(IVERT2,IEL)
C
99999 END
C
************************************************************************
      SUBROUTINE BDPRES(DP,KVERT,KNPR,KVBD,KMM,DCORVG,DVBDP,NVBD,P1,P2)
************************************************************************
*    Purpose:  Calculates integral boundary pressure
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      DIMENSION DP(*),KVERT(NNVE,*),KNPR(*),KVBD(*),KMM(2,*)
      DIMENSION DCORVG(2,*),DVBDP(*)
C
      COMMON /NSPTS/  KPU(2),KPP(4),KPX(4),KPI(2),DPI(2,2),DPF(2)
      SAVE
C
C
C
      P1=0D0
      DLEN=0D0
C
C
      IF (KPI(1).EQ.0) GOTO 19
C
      DO 10 IVBD=1,NVBD
C      
      DPAR=DVBDP(IVBD)
      IF ((DPAR.GE.DPI(1,1)).AND.(DPAR.LT.DPI(2,1))
     *                      .AND.(KNPR(KVBD(IVBD)).EQ.KPI(1))) THEN
       IVT1=KVBD(IVBD)
       PX1=DCORVG(1,IVT1)
       PY1=DCORVG(2,IVT1)
      ELSE
       GOTO 10
      ENDIF
C
      IF (IVT1.EQ.KMM(2,KPI(1))) THEN
       IVT2=KMM(1,KPI(1))
       PX2=DCORVG(1,IVT2)
       PY2=DCORVG(2,IVT2)
      ELSE
       IVBDH=IVBD+1
       DPARH=DVBDP(IVBDH)
       IVT2=KVBD(IVBDH)
       PX2=DCORVG(1,IVT2)
       PY2=DCORVG(2,IVT2)
      ENDIF
C
      DL=SQRT((PX2-PX1)**2+(PY2-PY1)**2)
      DLEN=DLEN+DL
      P1=P1+0.5D0*DL*(DP(IVT1)+DP(IVT2))
C
10    CONTINUE
C
      P1=P1/DLEN
C
C
C
19    P2=0D0
      DLEN=0D0
C
C
      IF (KPI(2).EQ.0) RETURN
C
      DO 20 IVBD=1,NVBD
C      
      DPAR=DVBDP(IVBD)
      IF ((DPAR.GE.DPI(1,2)).AND.(DPAR.LT.DPI(2,2))
     *                      .AND.(KNPR(KVBD(IVBD)).EQ.KPI(2))) THEN
       IVT1=KVBD(IVBD)
       PX1=DCORVG(1,IVT1)
       PY1=DCORVG(2,IVT1)
      ELSE
       GOTO 20
      ENDIF
C
      IF (IVT1.EQ.KMM(2,KPI(2))) THEN
       IVT2=KMM(1,KPI(2))
       PX2=DCORVG(1,IVT2)
       PY2=DCORVG(2,IVT2)
      ELSE
       IVBDH=IVBD+1
       DPARH=DVBDP(IVBDH)
       IVT2=KVBD(IVBDH)
       PX2=DCORVG(1,IVT2)
       PY2=DCORVG(2,IVT2)
      ENDIF
C
      DL=SQRT((PX2-PX1)**2+(PY2-PY1)**2)
      DLEN=DLEN+DL
      P2=P2+0.5D0*DL*(DP(IVT1)+DP(IVT2))
C
20    CONTINUE
C
      P2=P2/DLEN
C
C
C
      END
C
************************************************************************
      SUBROUTINE BDFORC(DU1,DU2,DP,KVERT,KMID,KVBD,KEBD,KMM,DCORVG,ELE,
     *                  DFW,DAW)
************************************************************************
*    Purpose:  Calculates lift (DFW) and drag (DAW)
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21,NNCOF=6)
      DIMENSION DU1(*),DU2(*),DP(*),DCORVG(2,*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KVBD(*),KEBD(*),KMM(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      COMMON /COFBAS/ COB(NNBAS,NNCOF)
C
      COMMON /NSPTS/  KPU(2),KPP(4),KPX(4),KPI(2),DPI(2,2),DPF(2)
      SAVE
C
C
C
      DFW=0D0
      DAW=0D0
      DLEN=0D0
C
C
      IF ((DPF(1).EQ.0D0).OR.(DPF(2).EQ.0D0)) RETURN
C
C
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
C
      NCUBP=1
      ICUBP=1
C
      IW1=1
      IW2=1
C
      DO 10 IVBD=1,NVBD
      IVT=KVBD(IVBD)
      IF (IVT.EQ.KMM(1,2)) GOTO 15
10    CONTINUE
C
      WRITE(6,*) 'ERROR 1 IN BDFORC'
      RETURN
C
15    IVBD1=IVBD
      ISTOP=0
C
      DO 100 IVBD=IVBD1,NVBD
      IF (ISTOP.EQ.1) GOTO 1000
      IEL =KEBD(IVBD)
      IVT1=KVBD(IVBD)
      IF (IVT1.EQ.KMM(2,2)) THEN
       IVT2=KVBD(IVBD1)
       ISTOP=1
      ELSE
       IVT2=KVBD(IVBD+1)
      ENDIF
C
      DO 101 II=1,4
101   IF (KVERT(II,IEL).EQ.IVT1) GOTO 102
      WRITE(6,*) 'WRONG 1'
      IW1=IW1+1
102   DO 103 II=1,4
103   IF (KVERT(II,IEL).EQ.IVT2) GOTO 104
      WRITE(6,*) 'WRONG 2'
      IW2=IW2+1
104   CONTINUE
C
      PX1=DCORVG(1,IVT1)
      PX2=DCORVG(1,IVT2)
      PY1=DCORVG(2,IVT1)
      PY2=DCORVG(2,IVT2)
      PXM=0.5D0*(PX1+PX2)
      PYM=0.5D0*(PY1+PY2)
      DLH=SQRT((PX2-PX1)**2+(PY2-PY1)**2)
C
      DTX= (PX2-PX1)/DLH
      DTY= (PY2-PY1)/DLH
      DNX=-(PY2-PY1)/DLH
      DNY= (PX2-PX1)/DLH
      DPCONT=DP(IEL)
      DLEN=DLEN+DLH
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
C
      DO 120 IVE = 1,4
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
120   CONTINUE
C
      DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
      DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
      DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
      DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
C
      XX=PXM
      YY=PYM
C
      CALL ELE(0D0,0D0,-2)
      CALL ELE(XX,YY,-3)
C
      DUT=0
      DO 130 JDFL=1,IDFL
      DUT=DUT+DU1(KDFG(JDFL))*DBAS(KDFL(JDFL),2)*DTX*DNX
     *       +DU2(KDFG(JDFL))*DBAS(KDFL(JDFL),2)*DTY*DNX
     *       +DU1(KDFG(JDFL))*DBAS(KDFL(JDFL),3)*DTX*DNY
     *       +DU2(KDFG(JDFL))*DBAS(KDFL(JDFL),3)*DTY*DNY
130   CONTINUE
C
      DFW=DFW+DLH*(DPF(1)*DUT*DNY-DPCONT*DNX)
      DAW=DAW-DLH*(DPF(1)*DUT*DNX+DPCONT*DNY)
C
100   CONTINUE
C
1000  DFW=2D0*DFW/DPF(2)      
      DAW=2D0*DAW/DPF(2) 
C
C
C
99999 END
C
