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
      CHARACTER SUB*6,FMT*15,CPARAM*120
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
      INCLUDE 'jump.inc'
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
      DUTX1=0
      DUTX2=0
      DUTY1=0
      DUTY2=0    
      DO 130 JDFL=1,IDFL
C
      DUX1=DU1(KDFG(JDFL))*DBAS(KDFL(JDFL),2)
      DUX2=DU2(KDFG(JDFL))*DBAS(KDFL(JDFL),2)
      DUY1=DU1(KDFG(JDFL))*DBAS(KDFL(JDFL),3)
      DUY2=DU2(KDFG(JDFL))*DBAS(KDFL(JDFL),3)
C
      DUT=DUT+DUX1*DTX*DNX
     *       +DUX2*DTY*DNX
     *       +DUY1*DTX*DNY
     *       +DUY2*DTY*DNY
c
      DUTX1=DUTX1+DUX1
      DUTX2=DUTX2+DUX2
      DUTY1=DUTY1+DUY1
      DUTY2=DUTY2+DUY2
C
130   CONTINUE
C
C ***  norm of symmetric part of velocity gradient
       IF (ISPGRAD .EQ. 1)DUGSQ=DUTX1**2+0.5d0*(DUTX2+DUTY1)**2+DUTY2**2
c
C ***  norm of gradient
       IF (ISPGRAD .EQ. 0)DUGSQ=DUTX1**2+DUTX2**2+DUTY1**2+DUTY2**2
C
       DDNY=DVISCO(DUGSQ)
C
c
      DFW=DFW+DLH*(DDNY*DUT*DNY-DPCONT*DNX)
      DAW=DAW-DLH*(DDNY*DUT*DNX+DPCONT*DNY)
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
************************************************************************
      SUBROUTINE VIBDFORC1(DU1,DU2,DP,KVERT,BST,NA1,KCOLA,KLDA,
     *                    B1,B2,KCOLB,KLDB,KNPR,KMID,KVBD,
     *                    KEBD,KMM,DCORVG,ELE,DFWVI,DAWVI,
     *                    dviappp,dvifppp,dviauuu,dvifuuu,
     *                    VD1,VD2,VL1,VL2,DAU1,DAU2,DVDBU,DVLBU)

************************************************************************
* Purpose: Calculates lift (DFW) and drag (DAW) by using volume integration
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-V,W-Z)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL B1,B2,BST
c      DOUBLE PRECISION BST

C
      LOGICAL  BDER
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21,NNCOF=6)
      DIMENSION DU1(*),DU2(*),DP(*),DCORVG(2,*)
      DIMENSION B1(*),B2(*),KCOLB(*),KLDB(*),KNPR(*)     
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KVBD(*),KEBD(*),KMM(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION BST(*),KCOLA(*),KLDA(*)
C
      DIMENSION VD1(*),VD2(*),VL1(*),VL2(*)
      DIMENSION DAU1(*),DAU2(*),DVDBU(*),DVLBU(*)      
c                                  
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      COMMON /COFBAS/ COB(NNBAS,NNCOF)
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /NSPTS/  KPU(2),KPP(4),KPX(4),KPI(2),DPI(2,2),DPF(2)
      
      SAVE
C
C
      CALL  LCL1( VD1,NU)
      CALL  LCL1( VD2,NU) 
      CALL  LCL1( VL1,NU)     
      CALL  LCL1( VL2,NU)   
      
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.      
C
      IF ((DPF(1).EQ.0D0).OR.(DPF(2).EQ.0D0)) RETURN
C
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
      IF (IVT.EQ.KMM(1,2)) GOTO 15!erster Knoten auf Komponente 2!
10    CONTINUE
C
      WRITE(6,*) 'ERROR 1 IN BDFORC'
      RETURN
C
15    IVBD1=IVBD
      ISTOP=0       
      DO 100 IVBD=IVBD1,NVBD!schleife über rankkomponente 2  
      IF (ISTOP.EQ.1) GOTO 1001
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
      CALL GETMBD(IMID,IV1,IV2,IVBD,KVBD,KEBD,KVERT,KMID,KNPR,INPR)
C       
       VD1(IMID) = 1.0D0
       VL2(IMID) = 1.0D0
C        
100   CONTINUE
C
1001   CONTINUE
C
      CALL LAX37(BST(0*NA+1),KCOLA,KLDA,NU,DU1,DAU1,1.0D0,0.0D0)
      CALL LAX37(BST(1*NA+1),KCOLA,KLDA,NU,DU2,DAU1,1.0D0,1.0D0)
      CALL LAX37(BST(2*NA+1),KCOLA,KLDA,NU,DU1,DAU2,1.0D0,0.0D0)
      CALL LAX37(BST(3*NA+1),KCOLA,KLDA,NU,DU2,DAU2,1.0D0,1.0D0)
C
       CALL   LSP1 (VD1,DAU1,NU,VD1AU1)
       CALL   LSP1 (VD2,DAU2,NU,VD2AU2)    
c 
C----------   ----------------
      CALL LTX39(B1,KCOLB,KLDB,NU,NP,VD1,DVDBU,1.0D0,0.0D0)
      CALL LTX39(B2,KCOLB,KLDB,NU,NP,VD2,DVDBU,1.0D0,1.0D0) 
C           
       CALL   LSP1 (DP,DVDBU,NP,DPVDBU)
C
      DAWVI= -(VD1AU1+VD2AU2)-(DPVDBU)
C  
c----------lift------------
       CALL   LSP1 (VL1,DAU1,NU,VL1AU1)
       CALL   LSP1 (VL2,DAU2,NU,VL2AU2)

      CALL LTX39(B1,KCOLB,KLDB,NU,NP,VL1,DVLBU,1.0D0,0.0D0)
      CALL LTX39(B2,KCOLB,KLDB,NU,NP,VL2,DVLBU,1.0D0,1.0D0)            
       CALL   LSP1 (DP,DVLBU,NP,DPVLBU)
    
      DFWVI= -(VL1AU1+VL2AU2)-(DPVLBU)

c-----------------------------------------------------------
c
1000     dviappp=-DPVDBU
         dvifppp=-DPVLBU
         dviauuu=-(VD1AU1+VD2AU2)
         dvifuuu=-(VL1AU1+VL2AU2)
C
         dviappp=2D0*dviappp/DPF(2)
         dvifppp=2D0*dvifppp/DPF(2) 
         dviauuu=2D0*dviauuu/DPF(2)
         dvifuuu=2D0*dvifuuu/DPF(2)           

      DFWVI=2D0*DFWVI/DPF(2)      
      DAWVI=2D0*DAWVI/DPF(2) 



99999 END
C
************************************************************************
      SUBROUTINE VIBDFORC(DU1,DU2,DP,KVERT,BST,NA1,KCOLA,KLDA,
     *                    B1,B2,KCOLB,KLDB,KNPR,KMID,KVBD,
     *                    KEBD,KMM,DCORVG,ELE,DFWVI,DAWVI,
     *                    dviappp,dvifppp,dviauuu,dvifuuu,
     *                    VD1,VD2,VL1,VL2,DAU1,DAU2,DVDBU,DVLBU)

************************************************************************
* Purpose: Calculates lift (DFW) and drag (DAW) by using volume integration
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-V,W-Z)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL B1,B2,BST
c      DOUBLE PRECISION BST

C
      LOGICAL  BDER
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21,NNCOF=6)
      DIMENSION DU1(*),DU2(*),DP(*),DCORVG(2,*)
      DIMENSION B1(*),B2(*),KCOLB(*),KLDB(*),KNPR(*)     
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KVBD(*),KEBD(*),KMM(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION BST(*),KCOLA(*),KLDA(*)
C
      DIMENSION VD1(*),VD2(*),VL1(*),VL2(*)
      DIMENSION DAU1(*),DAU2(*),DVDBU(*),DVLBU(*)      
c                                  
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      COMMON /COFBAS/ COB(NNBAS,NNCOF)
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /NSPTS/  KPU(2),KPP(4),KPX(4),KPI(2),DPI(2,2),DPF(2)
      
      SAVE
C
C
      CALL  LCL1( VD1,NU)
      CALL  LCL1( VD2,NU) 
      CALL  LCL1( VL1,NU)     
      CALL  LCL1( VL2,NU)   
      
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.      
C
      IF ((DPF(1).EQ.0D0).OR.(DPF(2).EQ.0D0)) RETURN
C
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
      IF (IVT.EQ.KMM(1,2)) GOTO 15!erster Knoten auf Komponente 2!
10    CONTINUE
C
      WRITE(6,*) 'ERROR 1 IN BDFORC'
      RETURN
C
15    IVBD1=IVBD
      ISTOP=0       
      DO 100 IVBD=IVBD1,NVBD!schleife über rankkomponente 2  
      IF (ISTOP.EQ.1) GOTO 1001
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
      CALL GETMBD(IMID,IV1,IV2,IVBD,KVBD,KEBD,KVERT,KMID,KNPR,INPR)
C       
       VD1(IMID) = 1.0D0
       VL2(IMID) = 1.0D0
C        
100   CONTINUE
C
1001   CONTINUE
C
      CALL LAX37(BST(0*NA+1),KCOLA,KLDA,NU,DU1,DAU1,1.0D0,0.0D0)
      CALL LAX37(BST(1*NA+1),KCOLA,KLDA,NU,DU2,DAU1,1.0D0,1.0D0)
      CALL LAX37(BST(2*NA+1),KCOLA,KLDA,NU,DU1,DAU2,1.0D0,0.0D0)
      CALL LAX37(BST(3*NA+1),KCOLA,KLDA,NU,DU2,DAU2,1.0D0,1.0D0)
C
       CALL   LSP1 (VD1,DAU1,NU,VD1AU1)
       CALL   LSP1 (VD2,DAU2,NU,VD2AU2)      
C
      CALL LTX39(B1,KCOLB,KLDB,NU,NP,VD1,DVDBU,1.0D0,0.0D0)
      CALL LTX39(B2,KCOLB,KLDB,NU,NP,VD2,DVDBU,1.0D0,1.0D0) 
C           
       CALL   LSP1 (DP,DVDBU,NP,DPVDBU)
C
      DAWVI= -(VD1AU1+VD2AU2)-(DPVDBU)
      
c----------lift------------
       CALL   LSP1 (VL1,DAU1,NU,VL1AU1)
       CALL   LSP1 (VL2,DAU2,NU,VL2AU2)

      CALL LTX39(B1,KCOLB,KLDB,NU,NP,VL1,DVLBU,1.0D0,0.0D0)
      CALL LTX39(B2,KCOLB,KLDB,NU,NP,VL2,DVLBU,1.0D0,1.0D0)            
       CALL   LSP1 (DP,DVLBU,NP,DPVLBU)
    
      DFWVI= -(VL1AU1+VL2AU2)-(DPVLBU)

c-----------------------------------------------------------
1000     dviappp=-DPVDBU
         dvifppp=-DPVLBU
         dviauuu=-(VD1AU1+VD2AU2)
         dvifuuu=-(VL1AU1+VL2AU2)
C
         dviappp=2D0*dviappp/DPF(2)
         dvifppp=2D0*dvifppp/DPF(2) 
         dviauuu=2D0*dviauuu/DPF(2)
         dvifuuu=2D0*dvifuuu/DPF(2)           

      DFWVI=2D0*DFWVI/DPF(2)      
      DAWVI=2D0*DAWVI/DPF(2) 



99999 END
c
************************************************************************
      SUBROUTINE VIBF_STG(DU1,DU2,DP,KVERT,KNPR,KMID,
     *                                 DCORVG,ELE,DFWVI,DAWVI)

************************************************************************
*       Purpose: Calculates lift (DFW) and drag (DAW) by STRONG integration
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)  
      
C       
      LOGICAL  BDER
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21,NNCOF=6)
      DIMENSION DU1(*),DU2(*),DP(*),DCORVG(2,*)
      DIMENSION KNPR(*)   
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
C       
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *    DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      COMMON /COFBAS/ COB(NNBAS,NNCOF)
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /NSPTS/  KPU(2),KPP(4),KPX(4),KPI(2),DPI(2,2),DPF(2)
      
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *    EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *    RLXSM,RLXSL,AMINMG,AMAXMG
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *    IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *    INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *    NSM,NSL,NSMFAC 

      DIMENSION DVD2(NU)     
      
C       SAVE
      EXTERNAL ELE
C       

      CALL  LCL1( DVD2, NU) 
      
C-------------------------------------------------------------------------
      
      ICHOICE=1!  (=1: KNPR; =2: ALFA)            
      
c-----------------------------------------------------------------

      DO 2 IEL0=1,NEL
        DO 20 IVE=1,4

          IM1=KMID(IVE,IEL0)
          if(knpr(im1).gt.0) then
            ivt=KNPR(IM1)
              IF(KNPR(ivt).EQ.2) THEN
                
                DVD2(IM1-NVT) = 1.0D0

              ENDIF
            endif
C       
  20      CONTINUE
  2     CONTINUE

      CALL STGFORCE(DFWVI,DAWVI,DU1,DU2,DP,NY,DVD2,
     *    KVERT,KMID,DCORVG,ELE,IELT,8)

      DFWVI=2D0*DFWVI/DPF(2)      
      DAWVI=2D0*DAWVI/DPF(2) 


99999 END





C=============================================================================
      subroutine stgforce(df1,df2,du1,du2,dp,dny,da,
     *    kvert,kmid,dcorvg,ele,ielt,icub)
c       
      implicit double precision (a,c-h,o-u,w-z),logical(b)
      character sub*6,fmt*15,cparam*120
c       
      parameter (nnbas=21,nnder=6,nncubp=36,nnve=4,nnab=21)
      dimension kvert(nnve,*),kmid(nnve,*),du1(*),du2(*),dp(*)
      dimension kdfg(nnbas),kdfl(nnbas),dcorvg(2,*),DA(*)
      dimension kentry(nnbas,nnbas),dentry(nnbas,nnbas)
c       
      common /output/ m,mt,mkeyb,mterm,merr,mprot,msys,mtrc,irecl8
      common /errctl/ ier,icheck
      common /char/   sub,fmt(3),cparam
      common /elem/   dx(nnve),dy(nnve),djac(2,2),detj,
     *    dbas(nnbas,nnder),bder(nnder),kve(nnve),iel
      common /triad/  nel,nvt,nmt,nve,nvel,nbct,nvbd
      common /cub/    dxi(nncubp,3),domega(nncubp),ncubp,icubp
      common /coaux1/ kdfg,kdfl,idfl
      save /output/,/errctl/,/char/,/elem/,/triad/,/cub/,/coaux1/
      COMMON /NEWTON/ DNNS,DNPW,DNPWSTART,DNP,DTOLERNEW,
     *                INEWTON,IELTNEWTON,IFIXMIN,IFIXMAX,
     *                ICUBNEWTON,ISPGRAD,INEWPRE  
      EXTERNAL ELE, DVISCO,DVISCOPRIM
c       
c       *** preparation - evaluation of parameters
      ier=0
c       *** which derivatives of basis functions are needed?
      do  i = 1,nnder
        bder(i)=.false.
      enddo
      
      bder(1)=.true.
      bder(2)=.true.
      bder(3)=.true.
      
c       *** dummy call of ele sets number of element
      ieltyp=-1
      call ele(0d0,0d0,ieltyp)

      idfl=ndfl(ieltyp)
      call cb2q(icub)
      if (ier.ne.0) goto 99999
c       
c       *** dummy call - ele may save arithmetic operations
      icubp=icub
      call ele(0d0,0d0,-2)
      
      
      df1=0.0
      df2=0.0
c       
c       *** loop over all elements
      do iel=1,nel
        call ndfgl(iel,1,ieltyp,kvert,kmid,kdfg,kdfl)
        if (ier.lt.0) goto 99999
c       
c       *** evaluation of coordinates of the vertices
        do ive = 1, nve
          jp=kvert(ive,iel)
          kve(ive)=jp
          dx(ive)=dcorvg(1,jp)
          dy(ive)=dcorvg(2,jp)
        enddo
c       
        dj1=0.5d0*(-dx(1)-dx(2)+dx(3)+dx(4))
        dj2=0.5d0*( dx(1)-dx(2)+dx(3)-dx(4))
        dj3=0.5d0*(-dy(1)+dy(2)-dy(3)+dy(4))
        dj4=0.5d0*(-dy(1)+dy(2)+dy(3)-dy(4))
c       
        call ele(0d0,0d0,-2)
c       *** loop over all cubature points
        do icubp = 1, ncubp
          
c       *** cubature point on the reference element
          xi1=dxi(icubp,1)
          xi2=dxi(icubp,2)
          
c       *** jacobian of the bilinear mapping onto the reference element
          djac(1,1)=0.5d0*(dx(2)-dx(1)+dj2)+0.5d0*dj2*xi2
          djac(1,2)=0.5d0*dj1+0.5d0*dj2*xi1
          djac(2,1)=0.5d0*dj4-0.5d0*dj3*xi2
          djac(2,2)=0.5d0*(dy(3)-dy(1)-dj4)-0.5d0*dj3*xi1
          detj=djac(1,1)*djac(2,2)-djac(1,2)*djac(2,1)
          om=domega(icubp)*detj
          
c       *** cubature point on the real element
          xx=0.5d0*(dx(1)+dx(2)+dj1)+0.5d0*(dx(2)-dx(1)+dj2)*xi1
     *        +0.5d0*dj1*xi2+0.5d0*dj2*xi1*xi2
          yy=0.5d0*(dy(1)+dy(3)+dj3)+0.5d0*dj4*xi1
     *        +0.5d0*(dy(3)-dy(1)-dj4)*xi2-0.5d0*dj3*xi1*xi2
          
c       *** evaluate the basis functions in the cubature point
          if((ielt.eq.2).or.(ielt.eq.3)) then
            call ele(xx,yy,-3)
          else
            call ele(xi1,xi2,-3)
          endif
          if (ier.lt.0) goto 99999
          
c       evaluate the solution values and derivatives in the cubature point     
          du1v=0d0 ! value
          du1x=0d0 ! x dreiv.
          du1y=0d0 ! y deriv
          do i=1,idfl
            ig=kdfg(i)
            du1v=du1v+du1(ig)*dbas(kdfl(i),1)
            du1x=du1x+du1(ig)*dbas(kdfl(i),2)
            du1y=du1y+du1(ig)*dbas(kdfl(i),3)
          enddo
          
          du2v=0d0 ! value
          du2x=0d0 ! x dreiv.
          du2y=0d0 ! y deriv
          do i=1,idfl
            ig=kdfg(i)
            du2v=du2v+du2(ig)*dbas(kdfl(i),1)
            du2x=du2x+du2(ig)*dbas(kdfl(i),2)
            du2y=du2y+du2(ig)*dbas(kdfl(i),3)
          enddo
          
          dav=0d0 ! value
          dax=0d0 ! x dreiv.
          day=0d0 ! y deriv
          do i=1,idfl
            ig=kdfg(i)
            dav=dav+da(ig)*dbas(kdfl(i),1)
            dax=dax+da(ig)*dbas(kdfl(i),2)
            day=day+da(ig)*dbas(kdfl(i),3)
          enddo
c       form the integrand
          
c       dn1=dav*dax
c       dn2=dav*day

          dn1=-dax
          dn2=-day
C ***  norm of symmetric part of velocity gradient
       IF (ISPGRAD .EQ. 1)DUGSQ=DU1X**2+0.5d0*(DU2X+DU1Y)**2+DU2Y**2
c
C ***  norm of gradient
       IF (ISPGRAD .EQ. 0)DUGSQ=DU1X**2+DU2X**2+DU1Y**2+DU2Y**2
C
       DDNY=DVISCO(DUGSQ)
C         dd=du1x**2+0.5*(du2x+du1y)**2+du2y**2

c$$$          if((dax.ne.0d0).or.(day.ne.0d0)) 
c$$$     *         print *,sqrt((xx-0.2)**2+(yy-0.2)**2)
          IF (ISPGRAD .EQ. 1)THEN
             ah1=-dp(iel)*dn1+DDNY*(2.0*du1x*dn1+(du2x+du1y)*dn2)
             ah2=-dp(iel)*dn2+DDNY*((du2x+du1y)*dn1+2.0*du2y*dn2)
          ENDIF
          IF (ISPGRAD .EQ. 0)THEN
             ah1=-dp(iel)*dn1+DDNY*(du1x*dn1+du1y*dn2)
             ah2=-dp(iel)*dn2+DDNY*(du2x*dn1+du2y*dn2)
          ENDIF
          df1=df1+ah1*om
          df2=df2+ah2*om
          
        enddo
      enddo
c       
c       
99999 end



