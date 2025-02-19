************************************************************************
      SUBROUTINE   INTUVD (DU1,DU2,DL1,DL2,DAUX,NVT,NEL,NVBD,
     *                     KMID,KVERT,KVBD,KMBD,DCORVG,UE)
************************************************************************
*
*    Purpose:  - Interpolates the solution vector (DU1,DU2) to
*                the vector (DL1,DL2) of dimension NVT with
*                values in the vertices
*              - the values of vertices at the boundary are computed
*                via the exact function UE
*
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=4)
C
      DIMENSION DU1(*),DU2(*),DL1(*),DL2(*),DAUX(*),
     *          KVERT(NNVE,*),KMID(NNVE,*),KVBD(*),KMBD(*),DCORVG(2,*)
C
      SAVE
C-----------------------------------------------------------------------
      DO 1 IVT=1,NVT
      DL1 (IVT)=0D0
      DL2 (IVT)=0D0
1     DAUX(IVT)=0D0
C
C
      DO 10 IEL=1,NEL
C
      IM1=KMID(1,IEL)-NVT
      IM2=KMID(2,IEL)-NVT
      IM3=KMID(3,IEL)-NVT
      IM4=KMID(4,IEL)-NVT
C
      IV1=KVERT(1,IEL)
      IV2=KVERT(2,IEL)
      IV3=KVERT(3,IEL)
      IV4=KVERT(4,IEL)
C
      DUH1=DU1(IM1)
      DUH2=DU1(IM2)
      DUH3=DU1(IM3)
      DUH4=DU1(IM4)
C
      DVH1=DU2(IM1)
      DVH2=DU2(IM2)
      DVH3=DU2(IM3)
      DVH4=DU2(IM4)
C
      DAUX(IV1)=DAUX(IV1)+1D0
      DAUX(IV2)=DAUX(IV2)+1D0
      DAUX(IV3)=DAUX(IV3)+1D0
      DAUX(IV4)=DAUX(IV4)+1D0
C
      DL1(IV1)=DL1(IV1) + 0.75D0*(DUH1+DUH4) - 0.25D0*(DUH2+DUH3)
      DL1(IV2)=DL1(IV2) + 0.75D0*(DUH2+DUH1) - 0.25D0*(DUH3+DUH4)
      DL1(IV3)=DL1(IV3) + 0.75D0*(DUH3+DUH2) - 0.25D0*(DUH4+DUH1)
      DL1(IV4)=DL1(IV4) + 0.75D0*(DUH4+DUH3) - 0.25D0*(DUH1+DUH2)
C
      DL2(IV1)=DL2(IV1) + 0.75D0*(DVH1+DVH4) - 0.25D0*(DVH2+DVH3)
      DL2(IV2)=DL2(IV2) + 0.75D0*(DVH2+DVH1) - 0.25D0*(DVH3+DVH4)
      DL2(IV3)=DL2(IV3) + 0.75D0*(DVH3+DVH2) - 0.25D0*(DVH4+DVH1)
      DL2(IV4)=DL2(IV4) + 0.75D0*(DVH4+DVH3) - 0.25D0*(DVH1+DVH2)
C
10    CONTINUE
C
      DO 20  IV=1,NVT
      DL1(IV)=DL1(IV)/DAUX(IV)
      DL2(IV)=DL2(IV)/DAUX(IV)
 20   CONTINUE
C
C=======================================================================
C *** boundary values
C=======================================================================
C
      DO 30  I=1,NVBD
      I1=I
      IF (I.EQ.1) THEN 
       I0=NVBD
      ELSE
       I0=I-1
      ENDIF
C
      IV =KVBD(I)
      IM1=KMBD(I1)
      IM0=KMBD(I0)
      IF ((IM1.LT.0).OR.(IM0.LT.0)) GOTO 30
C
      X=DCORVG(1,IV)
      Y=DCORVG(2,IV)
      DL1(IV)=UE(X,Y,1)
      DL2(IV)=UE(X,Y,2)
30    CONTINUE
C
      END
c
c
************************************************************************
      SUBROUTINE   INTPV (DP,DPL,DAUX,AREA,KVERT)
************************************************************************
*    Purpose:  - Interpolates the solution pressure DP to
*                the (REAL) vector VPL of dimension NVT with
*                values in the vertices
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL AREA
      PARAMETER (NNVE=4)
C
      DIMENSION DP(*),DPL(*),DAUX(*),KVERT(NNVE,*),AREA(*)
C
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      SAVE
C-----------------------------------------------------------------------
c
      DO 1 IVT=1,NVT
      DPL (IVT)=0D0
1     DAUX(IVT)=0D0
C
C
      DO 10 IEL=1,NEL
C
      DPIEL=DP(IEL)
      DAREA=DBLE(AREA(IEL))
C
      IV1=KVERT(1,IEL)
      IV2=KVERT(2,IEL)
      IV3=KVERT(3,IEL)
      IV4=KVERT(4,IEL)
C
      DPL(IV1)=DPL(IV1)+0.25D0*DAREA*DPIEL
      DPL(IV2)=DPL(IV2)+0.25D0*DAREA*DPIEL
      DPL(IV3)=DPL(IV3)+0.25D0*DAREA*DPIEL
      DPL(IV4)=DPL(IV4)+0.25D0*DAREA*DPIEL
C
      DAUX(IV1)=DAUX(IV1)+0.25D0*DAREA
      DAUX(IV2)=DAUX(IV2)+0.25D0*DAREA
      DAUX(IV3)=DAUX(IV3)+0.25D0*DAREA
      DAUX(IV4)=DAUX(IV4)+0.25D0*DAREA
C
10    CONTINUE
C
C
      DO 20 IVT=1,NVT
20    DPL(IVT)=DPL(IVT)/DAUX(IVT)
C
C
C      VPH=0E0
C      DO 30 IVT=1,NVT
C30    VPH=VPH+VPL(IVT)
C      VMWP=VPH/REAL(NVT)
C
C
C      DO 40 IVT=1,NVT
C40    VPL(IVT)=VPL(IVT)-VMWP
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE   SETARE  (AREA,NEL,KVERT,DCORVG)
************************************************************************
*
*   Purpose: - writes on  AREA(IEL)  the area of the element IEL,
*              IEL=1,...,NEL
*            - writes on  AREA(NEL+1) the sum of all  AREA(IEL)
*            - KVERT,DCORVG are the usual FEAT arrays
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL AREA
C
      PARAMETER (NNVE=4)
C
      DIMENSION  AREA(*),KVERT(NNVE,*),DCORVG(2,*)
C=======================================================================
      SUM=0.D0
      DO  11  IEL=1,NEL
C
      I1=KVERT(1,IEL)
      I2=KVERT(2,IEL)
      I3=KVERT(3,IEL)
      I4=KVERT(4,IEL)
C
      X1=DCORVG(1,I1)
      X2=DCORVG(1,I2)
      X3=DCORVG(1,I3)
      X4=DCORVG(1,I4)
C
      Y1=DCORVG(2,I1)
      Y2=DCORVG(2,I2)
      Y3=DCORVG(2,I3)
      Y4=DCORVG(2,I4)
C
      AAA=0.5D0*(  DABS((X1-X2)*(Y3-Y2)-(Y1-Y2)*(X3-X2))
     *            +DABS((X1-X4)*(Y3-Y4)-(Y1-Y4)*(X3-X4)) )
      AREA(IEL)=REAL(AAA)
      SUM=SUM+AAA
  11  CONTINUE
C
      AREA(NEL+1)=REAL(SUM)
C
      END
c
c
c
************************************************************************
      SUBROUTINE    TOL20A  (P,AREA,NEL,INEUM)
************************************************************************
*
*    Purpose: - Transforms the vector P into the space L2_0
*             - uses the vector AREA with the areas of all elements
*               and on AREA(NEL+1) the sum of all AREA(IEL)
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL AREA
C
      DIMENSION P(*),AREA(*)
C
      IF (INEUM.EQ.1) RETURN
C
      MIN=0.0D0
      PINT=0.0d0
      DO 2  IEL=1,NEL
      PINT=PINT+P(IEL)*DBLE(AREA(IEL))
 2    CONTINUE
C
        C=(PINT)/DBLE(AREA(NEL+1))

C
      DO 3  IEL=1,NEL
      P(IEL)=P(IEL)-C
 3    CONTINUE
c

c
      END
c
************************************************************************
      SUBROUTINE    MEAN  (P,AREA,NEL)!,DMEAN)
************************************************************************
*
*    Purpose: - Transforms the vector P into the space L2_0
*             - uses the vector AREA with the areas of all elements
*               and on AREA(NEL+1) the sum of all AREA(IEL)
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL AREA
      COMMON /BEARING/ DCU1,DCU2,DMEAN
      DIMENSION P(*),AREA(*)
C
      RETURN
C
      MIN=0.0D0
      PINT=0.0d0
      DO 2  IEL=1,NEL
      PINT=PINT+P(IEL)*DBLE(AREA(IEL))
 2    CONTINUE
C
        C=(PINT)/DBLE(AREA(NEL+1))

C
      DO 3  IEL=1,NEL
      P(IEL)=P(IEL)-C+DMEAN
 3    CONTINUE
c

c
      END
c
************************************************************************
      SUBROUTINE    prtavg  (P,AREA,NEL,INEUM)
************************************************************************
*
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL AREA
C
      DIMENSION P(*),AREA(*)
C
      return
c
      MIN=0.0D0
      PINT=0.D0
      DO 2  IEL=1,NEL
      PINT=PINT+P(IEL)*DBLE(AREA(IEL))
 2    CONTINUE
C
      C=PINT/DBLE(AREA(NEL+1))
C
      print *,'pressure avg',c
c

c
      END
c
c
c
************************************************************************
      SUBROUTINE   SETLEV (ISETLV)
************************************************************************
*
*   Purpose:  sets all data for current level ILEV (from /MGPAR/)
*            
*   Input:
*   -------
*     ILEV        - current level number (from /MGPAR/)
*     ISETLV >=1  - update of /TRIAA/,TRIAD/ 
*            >=2  - update of /LEVDIM/,/ADRFLD/
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299,NNLEV=9,  NNWORK=1)
C
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C *** Standard COMMON blocks
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
C
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
C
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
C
      common /locvis/ klny(NNLEV),kny
C
      SAVE 
C-----------------------------------------------------------------------
C *** elementary check
      IF (ILEV.LT.NLMIN .OR. ILEV.GT.NLMAX .OR. ISETLV.LT.1  .OR.
     *    ISETLV.GT.2 )  THEN
        WRITE(MTERM,*) 'ERROR in SETLEV: ILEV or ISETLV is wrong'
        STOP
      ENDIF
C
C *** update of /TRIAD/,/TRIAA/
C
      NEL =KNEL  (ILEV)
      NVT =KNVT  (ILEV)
      NMT =KNMT  (ILEV)
      NVEL=KNVEL (ILEV)
      NVBD=KNVBD (ILEV)
C
      LCORVG=KLCVG  (ILEV)
      LCORMG=KLCMG  (ILEV)
      LVERT =KLVERT (ILEV)
      LMID  =KLMID  (ILEV)
      LADJ  =KLADJ  (ILEV)
      LVEL  =KLVEL  (ILEV)
      LMEL  =KLMEL  (ILEV)
      LNPR  =KLNPR  (ILEV)
      LMM   =KLMM   (ILEV)
      LVBD  =KLVBD  (ILEV)
      LEBD  =KLEBD  (ILEV)
      LBCT  =KLBCT  (ILEV)
      LVBDP =KLVBDP (ILEV)
      LMBDP =KLMBDP (ILEV)
C
C *** update of /LEVDIM/,/ADRFLD/ if  ISETLV=2
C
      IF (ISETLV.EQ.2)  THEN
C
         NA =KNA  (ILEV)
         NB =KNB  (ILEV)
         NU =KNU  (ILEV)
         NP =KNP  (ILEV)
         NUP=KNUP (ILEV)
C
         KA1  =L(KLA    (ILEV))
         KST1 =L(KLST   (ILEV))
         KM1  =L(KLM    (ILEV))
         KCOLA=L(KLCOLA (ILEV))
         KLDA =L(KLLDA  (ILEV))
         KB1  =L(KLB1   (ILEV))
         KB2  =L(KLB2   (ILEV))
         KCOLB=L(KLCOLB (ILEV))
         KLDB =L(KLLDB  (ILEV))
         KU1  =L(KLUP   (ILEV))
         KU2  =KU1+NU
         KP   =KU2+NU
         KF1  =L(KLF12P (ILEV))
         KF2  =KF1+NU
         KFP  =KF2+NU
         KAUX1=L(KLAUX  (ILEV))
         KAUX2=KAUX1+NU
         KAUXP=KAUX2+NU
         kny=L(klny(ILEV))
C
      ENDIF
C
      END
c
c
c
************************************************************************
      SUBROUTINE   U2ISO (DCORVG,KVERT,KMID,KADJ,DVIND,DX,DU1,DU2)
************************************************************************
*
*   Purpose: - calculates streamlines
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNVE=4,NNARR=299,NNWORK=1)
      DIMENSION DCORVG(2,*),KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*),
     *          DVIND(*),DX(*),DU1(*),DU2(*)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
C
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      SAVE
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      REAL  VWORK
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C=======================================================================
      SUB='U2ISO '
C
C
C
C-----------------------------------------------------------------------
      DO 1 IVT=1,NVT
      DX   (IVT)=0D0
1     DVIND(IVT)=0D0
C
C
C *** start with element 1 and its first vertex
      DX(KVERT(1,1))   =0D0
      DVIND(KVERT(1,1))=1D0
C
C
      DO 10 IEH=1,NEL
C
      IEL=IEH
      IH=0
C
      DO 12 IVE= 1,NVE
      JVE=KVERT(IVE,IEL)
      IHV=INT(DVIND(JVE))
      IH=IH+IHV
      IF (IHV.GE.1) IND=IVE
12    CONTINUE
C
C
      IF ((IH.GE.NVE).OR.(IH.EQ.0)) GOTO 20
C
C
13    DO 14 IVE=1,NVE-1
C
      INDH=MOD(IND,NVE)+1
      IVTH=KVERT(INDH,IEL)
C
      IF (DVIND(IVTH).GE.1D0) GOTO 14
      DVIND(IVTH)=1D0
C
      IVT =KVERT(IND,IEL)
      IMID=KMID (IND,IEL)-NVT
C
      PX1=DCORVG(1,IVT)
      PY1=DCORVG(2,IVT)
      PX2=DCORVG(1,IVTH)
      PY2=DCORVG(2,IVTH)
      DN1 = PY2-PY1
      DN2 =-PX2+PX1
C
      DX(IVTH)=DX(IVT)+(DU1(IMID)*DN1+DU2(IMID)*DN2)
14    IND=INDH
C
C
20    DO 22 IME=1,NVE
C
      IELH=KADJ(IME,IEL)
      IF (IELH.EQ.0) GOTO 22
      IH=0
C
      DO 24 IVE=1,NVE
      JVE=KVERT(IVE,IELH)
      IHV=INT(DVIND(JVE))
      IH=IH+IHV
      IF (IHV.GE.1) INDH=IVE
24    CONTINUE
C
C
      IF ((IH.LT.NVE).AND.(IH.GT.0)) GOTO 30
C
C
22    CONTINUE
      GOTO 10
C
C
30    IEL=IELH
      IND=INDH
      GOTO 13
C
10    CONTINUE
C
      DXH=DX(1)
      DO 40 IVT=1,NVT
      DX(IVT)=DX(IVT)-DXH
40    CONTINUE
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE C2N2DM (DPC,DPL,KMID,KADJ,NEL,NMT,NVT,IPAR)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=4)
      DIMENSION DPC(*),DPL(*),KMID(NNVE,*),KADJ(NNVE,*)
C
C-----------------------------------------------------------------------
C
C *** constant to nonconforming = 0
      IF (IPAR.EQ.0) THEN
C
       DO 10 IEL=1,NEL
       DPH  =DPC(IEL)
C
       DO 20 IVE=1,4
       IADJ=KADJ(IVE,IEL)
       IMID=KMID(IVE,IEL)-NVT
C
       IF (IADJ.EQ.0)   DPL(IMID)=DPH
       IF (IADJ.GT.IEL) DPL(IMID)=0.5D0*(DPH+DPC(IADJ))
C
20     CONTINUE
10     CONTINUE
C
      ELSE
C
       DO 110 IEL=1,NEL
       DPC(IEL)=0.25D0*( DPL(KMID(1,IEL)-NVT)+DPL(KMID(2,IEL)-NVT)
     *                  +DPL(KMID(3,IEL)-NVT)+DPL(KMID(4,IEL)-NVT))
110    CONTINUE
C
      ENDIF      
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE GRDIST(DCORVG,KNPR,NVT,DIEPS)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DCORVG(2,*),KNPR(*)
C
      H=1D0/(SQRT(DBLE(NVT))-1)
      HDIST=DIEPS*H
C
      DO 1 IVT=1,NVT
      IF (KNPR(IVT).EQ.0) THEN
       DCORVG(1,IVT)=DCORVG(1,IVT)+DBLE((-1)**MOD(IVT,17))*HDIST
       DCORVG(2,IVT)=DCORVG(2,IVT)+DBLE((-1)**IVT)*HDIST
      ENDIF
1     CONTINUE
C
      END
C
C
C
************************************************************************
      SUBROUTINE CHCOOR(DCORVG,KVERT,KMID,KADJ,KNPR,NEL,NVT)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DCORVG(2,*),KVERT(4,*),KMID(4,*),KADJ(4,*),KNPR(*)
C
      DO 10 IEL=1,NEL/4
C
      IADJ3=KADJ(2,IEL)
      IADJ4=KADJ(3,IEL)
      IVT1=KVERT(2,IEL)
      IVT2=KVERT(4,IEL)
      IVT3=KVERT(2,IADJ3)
      IVT4=KVERT(4,IADJ4)
      IVTM=KVERT(3,IEL)
C
      PX1=DCORVG(1,IVT1)
      PX2=DCORVG(1,IVT2)
      PX3=DCORVG(1,IVT3)
      PX4=DCORVG(1,IVT4)
C
      PY1=DCORVG(2,IVT1)
      PY2=DCORVG(2,IVT2)
      PY3=DCORVG(2,IVT3)
      PY4=DCORVG(2,IVT4)
C
      PXM=DCORVG(1,IVTM)
      PYM=DCORVG(2,IVTM)
C
      PX=0.25D0*(PX1+PX2+PX3+PX4)
      PY=0.25D0*(PY1+PY2+PY3+PY4)
C
      DCORVG(1,IVTM)=PX
      DCORVG(2,IVTM)=PY
C
10    CONTINUE
C
      END
C
C
C
C
************************************************************************
      SUBROUTINE EM30(XI1,XI2,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNCOF=6)
      DIMENSION DXM(4),DYM(4),DLX(4),DLY(4),A(4,4),F(4),CKH(4),CK(4,4)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COFBAS/ COB(NNBAS,NNCOF)
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/,/CHAR/,/CUB/
C
      F1(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)=1D0
      F2(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)=CA1*X  +CB1*Y
      F3(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)=CA2*X  +CB2*Y
      F4(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)=CA3*X*X+CB3*X*Y+CC3*Y*Y
C
C
      SUB='EM30'
      IF (ICHECK.GE.998) CALL OTRC('EM30  ','01/08/94')
C
C
C *** Dummy call
      IF (IPAR.EQ.-1) THEN
       IER=0
       IPAR=30
       GOTO 99999
      ENDIF
C
C
      IF (IPAR.EQ.-2) THEN
       DO 20 IVE=1,NVE
       DXM(IVE)=0.5D0*(DX(IVE)+DX(MOD(IVE,4)+1))
       DYM(IVE)=0.5D0*(DY(IVE)+DY(MOD(IVE,4)+1))
       DLX(IVE)=0.5D0*(DX(MOD(IVE,4)+1)-DX(IVE))
       DLY(IVE)=0.5D0*(DY(MOD(IVE,4)+1)-DY(IVE))
20     CONTINUE
C
       CA1=(DXM(2)-DXM(4))/SQRT((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
       CB1=(DYM(2)-DYM(4))/SQRT((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
       CA2=(DXM(3)-DXM(1))/SQRT((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)
       CB2=(DYM(3)-DYM(1))/SQRT((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)
       CA3=CA1**2-CA2**2
       CB3=2D0*(CA1*CB1-CA2*CB2)
       CC3=CB1**2-CB2**2
C
       DO 22 IA=1,4
       PXL=DXM(IA)-SQRT(1D0/3D0)*DLX(IA)
       PYL=DYM(IA)-SQRT(1D0/3D0)*DLY(IA)
       PXU=DXM(IA)+SQRT(1D0/3D0)*DLX(IA)
       PYU=DYM(IA)+SQRT(1D0/3D0)*DLY(IA)
       A(IA,1)=0.5D0*( F1(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
     *                +F1(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
       A(IA,2)=0.5D0*( F2(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
     *                +F2(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
       A(IA,3)=0.5D0*( F3(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
     *                +F3(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
       A(IA,4)=0.5D0*( F4(PXL,PYL,CA1,CB1,CA2,CB2,CA3,CB3,CC3)
     *                +F4(PXU,PYU,CA1,CB1,CA2,CB2,CA3,CB3,CC3))
22     CONTINUE
C
       CALL INVERT(A,F,CKH,0)       
C
       DO 24 IK1=1,4
       DO 24 IK2=1,4
24     CK(IK1,IK2)=A(IK2,IK1)
C
       DO 26 IK=1,4
       COB(IK,1)=CK(IK,4)*CA3
       COB(IK,2)=CK(IK,4)*CC3
       COB(IK,3)=CK(IK,4)*CB3
       COB(IK,4)=CK(IK,2)*CA1+CK(IK,3)*CA2
       COB(IK,5)=CK(IK,2)*CB1+CK(IK,3)*CB2
       COB(IK,6)=CK(IK,1)
26     CONTINUE
      ENDIF
C
C

      IER=-1
C *** No second order derivatives available
C *** Used for second order problems only
      IF (BDER(4).OR.BDER(5).OR.BDER(6)) GOTO 99999
C
      IER=0
C
      IF (.NOT.BDER(1)) GOTO 101
C
C *** Function values
      DBAS(1,1)= COB(1,1)*XI1**2+COB(1,2)*XI2**2+COB(1,3)*XI1*XI2
     *          +COB(1,4)*XI1   +COB(1,5)*XI2   +COB(1,6)
      DBAS(2,1)= COB(2,1)*XI1**2+COB(2,2)*XI2**2+COB(2,3)*XI1*XI2
     *          +COB(2,4)*XI1   +COB(2,5)*XI2   +COB(2,6)
      DBAS(3,1)= COB(3,1)*XI1**2+COB(3,2)*XI2**2+COB(3,3)*XI1*XI2
     *          +COB(3,4)*XI1   +COB(3,5)*XI2   +COB(3,6)
      DBAS(4,1)= COB(4,1)*XI1**2+COB(4,2)*XI2**2+COB(4,3)*XI1*XI2
     *          +COB(4,4)*XI1   +COB(4,5)*XI2   +COB(4,6)
101   IF (.NOT.(BDER(2).OR.BDER(3))) GOTO 99999
C
C *** First order derivatives
      IF (.NOT.BDER(2)) GOTO 102
      DBAS(1,2)= 2D0*COB(1,1)*XI1+COB(1,3)*XI2+COB(1,4)
      DBAS(2,2)= 2D0*COB(2,1)*XI1+COB(2,3)*XI2+COB(2,4)
      DBAS(3,2)= 2D0*COB(3,1)*XI1+COB(3,3)*XI2+COB(3,4)
      DBAS(4,2)= 2D0*COB(4,1)*XI1+COB(4,3)*XI2+COB(4,4)
C
102   IF (.NOT.BDER(3)) GOTO 99999
      DBAS(1,3)= 2D0*COB(1,2)*XI2+COB(1,3)*XI1+COB(1,5)
      DBAS(2,3)= 2D0*COB(2,2)*XI2+COB(2,3)*XI1+COB(2,5)
      DBAS(3,3)= 2D0*COB(3,2)*XI2+COB(3,3)*XI1+COB(3,5)
      DBAS(4,3)= 2D0*COB(4,2)*XI2+COB(4,3)*XI1+COB(4,5)
C
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE EM31(XI1,XI2,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNCOF=6)
      DIMENSION DXM(4),DYM(4),A(4,4),F(4),CKH(4),CK(4,4)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COFBAS/ COB(NNBAS,NNCOF)
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/,/CHAR/,/CUB/
C
      F1(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)=1D0
      F2(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)=CA1*X  +CB1*Y
      F3(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)=CA2*X  +CB2*Y
      F4(X,Y,CA1,CB1,CA2,CB2,CA3,CB3,CC3)=CA3*X*X+CB3*X*Y+CC3*Y*Y
C
C
      SUB='EM31'
      IF (ICHECK.GE.998) CALL OTRC('EM31  ','01/08/94')
C
C
C *** Dummy call
      IF (IPAR.EQ.-1) THEN
       IER=0
       IPAR=31
       GOTO 99999
      ENDIF
C
C
      IF (IPAR.EQ.-2) THEN
       DO 20 IVE=1,NVE
       DXM(IVE)=0.5D0*(DX(IVE)+DX(MOD(IVE,4)+1))
       DYM(IVE)=0.5D0*(DY(IVE)+DY(MOD(IVE,4)+1))
20     CONTINUE
C
       CA1=(DXM(2)-DXM(4))/SQRT((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
       CB1=(DYM(2)-DYM(4))/SQRT((DXM(2)-DXM(4))**2+(DYM(2)-DYM(4))**2)
       CA2=(DXM(3)-DXM(1))/SQRT((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)
       CB2=(DYM(3)-DYM(1))/SQRT((DXM(1)-DXM(3))**2+(DYM(1)-DYM(3))**2)
       CA3=CA1**2-CA2**2
       CB3=2D0*(CA1*CB1-CA2*CB2)
       CC3=CB1**2-CB2**2
C
       DO 22 IA=1,4
       A(IA,1)=F1(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)       
       A(IA,2)=F2(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)       
       A(IA,3)=F3(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)       
       A(IA,4)=F4(DXM(IA),DYM(IA),CA1,CB1,CA2,CB2,CA3,CB3,CC3)
22     CONTINUE
C
       CALL INVERT(A,F,CKH,0)       
C
       DO 24 IK1=1,4
       DO 24 IK2=1,4
24     CK(IK1,IK2)=A(IK2,IK1)
C
       DO 26 IK=1,4
       COB(IK,1)=CK(IK,4)*CA3
       COB(IK,2)=CK(IK,4)*CC3
       COB(IK,3)=CK(IK,4)*CB3
       COB(IK,4)=CK(IK,2)*CA1+CK(IK,3)*CA2
       COB(IK,5)=CK(IK,2)*CB1+CK(IK,3)*CB2
       COB(IK,6)=CK(IK,1)
26     CONTINUE
      ENDIF
C
C

      IER=-1
C *** No second order derivatives available
C *** Used for second order problems only
      IF (BDER(4).OR.BDER(5).OR.BDER(6)) GOTO 99999
C
      IER=0
C
      IF (.NOT.BDER(1)) GOTO 101
C
C *** Function values
      DBAS(1,1)= COB(1,1)*XI1**2+COB(1,2)*XI2**2+COB(1,3)*XI1*XI2
     *          +COB(1,4)*XI1   +COB(1,5)*XI2   +COB(1,6)
      DBAS(2,1)= COB(2,1)*XI1**2+COB(2,2)*XI2**2+COB(2,3)*XI1*XI2
     *          +COB(2,4)*XI1   +COB(2,5)*XI2   +COB(2,6)
      DBAS(3,1)= COB(3,1)*XI1**2+COB(3,2)*XI2**2+COB(3,3)*XI1*XI2
     *          +COB(3,4)*XI1   +COB(3,5)*XI2   +COB(3,6)
      DBAS(4,1)= COB(4,1)*XI1**2+COB(4,2)*XI2**2+COB(4,3)*XI1*XI2
     *          +COB(4,4)*XI1   +COB(4,5)*XI2   +COB(4,6)
101   IF (.NOT.(BDER(2).OR.BDER(3))) GOTO 99999
C
C *** First order derivatives
      IF (.NOT.BDER(2)) GOTO 102
      DBAS(1,2)= 2D0*COB(1,1)*XI1+COB(1,3)*XI2+COB(1,4)
      DBAS(2,2)= 2D0*COB(2,1)*XI1+COB(2,3)*XI2+COB(2,4)
      DBAS(3,2)= 2D0*COB(3,1)*XI1+COB(3,3)*XI2+COB(3,4)
      DBAS(4,2)= 2D0*COB(4,1)*XI1+COB(4,3)*XI2+COB(4,4)
C
102   IF (.NOT.BDER(3)) GOTO 99999
      DBAS(1,3)= 2D0*COB(1,2)*XI2+COB(1,3)*XI1+COB(1,5)
      DBAS(2,3)= 2D0*COB(2,2)*XI2+COB(2,3)*XI1+COB(2,5)
      DBAS(3,3)= 2D0*COB(3,2)*XI2+COB(3,3)*XI1+COB(3,5)
      DBAS(4,3)= 2D0*COB(4,2)*XI2+COB(4,3)*XI1+COB(4,5)
C
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE INVERT(A,F,X,IPAR)
      DOUBLE PRECISION A,B,F,X
      PARAMETER (NDIM=4)
      DIMENSION A(NDIM,NDIM),B(NDIM,NDIM),F(NDIM),X(NDIM),
     *          MERKX(NDIM),MERKY(NDIM)
C
C
      IF (IPAR.EQ.0) THEN
       CALL AUSTAU(NDIM,NDIM,A,B,MERKX,MERKY,IFEHL)
       DO 10 IA=1,NDIM
       DO 10 IB=1,NDIM
10     A(IA,IB)=B(IA,IB)
      ENDIF
C
      IF (IPAR.EQ.1) THEN
       DO 20 IA=1,NDIM
       X(IA)=0D0
       DO 22 IB=1,NDIM
       X(IA)=X(IA)+A(IA,IB)*F(IB)
22     CONTINUE
20     CONTINUE
      ENDIF
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE AUSTAU(NDIM,N,A,B,MERKX,MERKY,IFEHL)
      DOUBLE PRECISION A,B,HILF,PIVOT
      DIMENSION A(NDIM,N),B(NDIM,N),MERKX(N),MERKY(N)
C
      IFEHL=1
      DO 100 I=1,N
      MERKX(I)=0
      MERKY(I)=0
      DO 100 L=1,N
100   B(I,L)=A(I,L)
C
      DO 400 I=1,N
      PIVOT=0D0
      DO 200 IX=1,N
      IF (MERKX(IX).NE.0) GOTO 200
      DO 180 IY=1,N
      IF (MERKY(IY).NE.0) GOTO 180
      IF (ABS(B(IX,IY)).LE.ABS(PIVOT)) GOTO 180
      PIVOT=B(IX,IY)
      INDX=IX
      INDY=IY
180   CONTINUE
200   CONTINUE
C
      IF (ABS(PIVOT).LE.0.0) GOTO 770
      MERKX(INDX)=INDY
      MERKY(INDY)=INDX
      B(INDX,INDY)=1D0/PIVOT
      DO 300 L=1,N
      IF (L.EQ.INDX) GOTO 300
      DO 280 M=1,N
      IF (M.EQ.INDY) GOTO 280
      B(L,M)=B(L,M)-B(L,INDY)*B(INDX,M)/PIVOT
280   CONTINUE
300   CONTINUE
C
      DO 390 IX=1,N
      IF (IX.NE.INDX) B(IX,INDY)=B(IX,INDY)/PIVOT
390   CONTINUE
      DO 400 IY=1,N
      IF (IY.NE.INDY) B(INDX,IY)=-B(INDX,IY)/PIVOT
400   CONTINUE
C
C
      DO 500 I=2,N
      IX=I-1
      IF (MERKX(IX).EQ.IX) GOTO 500
      DO 450 J=1,N
      IY=J
      IF (MERKX(IY).EQ.IX) GOTO 460
450   CONTINUE
460   DO 490 K=1,N
      HILF=B(IX,K)
      B(IX,K)=B(IY,K)
490   B(IY,K)=HILF
      MERKX(IY)=MERKX(IX)
500   MERKX(IX)=IX
      DO 600 I=2,N
      IX=I-1
      IF (MERKY(IX).EQ.IX) GOTO 600
      DO 550 J=1,N
      IY=J
      IF (MERKY(IY).EQ.IX) GOTO 560
550   CONTINUE
560   DO 590 K=1,N
      HILF=B(K,IX)
      B(K,IX)=B(K,IY)
590   B(K,IY)=HILF
      MERKY(IY)=MERKY(IX)
600   MERKY(IX)=IX
C
      IFEHL=0
770   RETURN
      END
C
C
C
************************************************************************
      SUBROUTINE CRITAD(TIMEIN,TIMENS,TIMEST,EPSADI,EPSADL,EPSAD,IADIN)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (PI=3.1415926535897931D0)
C
      EPSAD=EPSADL
C
      IF (TIMEIN.GT.0) THEN
       TDIFF=TIMENS-TIMEST
C
       IF (IADIN.EQ.0) THEN
        IF (TDIFF.LE.TIMEIN) THEN
         EPSAD=EPSADI
        ELSE
         EPSAD=EPSADL
        ENDIF
       ENDIF
C
       IF (IADIN.EQ.1) THEN
        IF (TDIFF.LE.TIMEIN) THEN
         EPSAD=EPSADI+TDIFF/TIMEIN*(EPSADL-EPSADI)
        ELSE
         EPSAD=EPSADL
        ENDIF
       ENDIF
C
       IF (IADIN.EQ.2) THEN
        IF (TDIFF.LE.TIMEIN) THEN
         EPSAD=EPSADI**(1D0-TDIFF/TIMEIN)*EPSADL**(TDIFF/TIMEIN)
        ELSE
         EPSAD=EPSADL
        ENDIF
       ENDIF
C
      ENDIF
C
C
C     
      END
C
c================================================================
      SUBROUTINE VORT(DU1,DU2,DVORT,KVERT,KMID,DCORVG,ELE)
c================================================================
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21)
      CHARACTER FMT*15,SUB*6,CPARAM*120
      DIMENSION DU1(*), DU2(*),DVORT(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      SAVE
C
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
C
      DO 100 IEL=1,NEL
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
      DO 120 IVE = 1,4
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
120   CONTINUE
      PXM=0.25D0*(DX(1)+DX(2)+DX(3)+DX(4))
      PYM=0.25D0*(DY(1)+DY(2)+DY(3)+DY(4))
C
      CALL ELE(0D0,0D0,-2)
      CALL ELE(PXM,PYM,-3)
C
      UH1 =0D0
      UH2 =0D0
      UH1X=0D0
      UH2X=0D0
      UH1Y=0D0
      UH2Y=0D0
C
      DO 210 JDOFE=1,IDFL
      IEQ=KDFG(JDOFE)
      ILO=KDFL(JDOFE)
      UH1 =UH1 +DU1(IEQ)*DBAS(ILO,1)
      UH2 =UH2 +DU2(IEQ)*DBAS(ILO,1)
      UH1X=UH1X+DU1(IEQ)*DBAS(ILO,2)
      UH2X=UH2X+DU2(IEQ)*DBAS(ILO,2)
      UH1Y=UH1Y+DU1(IEQ)*DBAS(ILO,3)
      UH2Y=UH2Y+DU2(IEQ)*DBAS(ILO,3)
C
210   CONTINUE
c
      DVORT(IEL)=UH1Y-UH2X
C
100   CONTINUE      
c
99999 END
C
c================================================================
      SUBROUTINE NGRAD(DU1,DU2,DNG,KVERT,KMID,DCORVG,ELE)
c================================================================
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21)
      CHARACTER FMT*15,SUB*6,CPARAM*120
      DIMENSION DU1(*), DU2(*),DNG(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      common /fluidp/ pla,ple
      SAVE
C
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
C
      DO 100 IEL=1,NEL
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
      DO 120 IVE = 1,4
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
120   CONTINUE
      PXM=0.25D0*(DX(1)+DX(2)+DX(3)+DX(4))
      PYM=0.25D0*(DY(1)+DY(2)+DY(3)+DY(4))
C
      CALL ELE(0D0,0D0,-2)
      CALL ELE(PXM,PYM,-3)
C
      UH1 =0D0
      UH2 =0D0
      UH1X=0D0
      UH2X=0D0
      UH1Y=0D0
      UH2Y=0D0
C
      DO 210 JDOFE=1,IDFL
      IEQ=KDFG(JDOFE)
      ILO=KDFL(JDOFE)
      UH1 =UH1 +DU1(IEQ)*DBAS(ILO,1)
      UH2 =UH2 +DU2(IEQ)*DBAS(ILO,1)
      UH1X=UH1X+DU1(IEQ)*DBAS(ILO,2)
      UH2X=UH2X+DU2(IEQ)*DBAS(ILO,2)
      UH1Y=UH1Y+DU1(IEQ)*DBAS(ILO,3)
      UH2Y=UH2Y+DU2(IEQ)*DBAS(ILO,3)
C
210   CONTINUE
c
      DNG(IEL)=(UH1X**2
     *          +0.5d0*(UH2X+UH1Y)**2+UH2Y**2)**(0.5d0) 
C
100   CONTINUE
C
99999 END
c================================================================
      SUBROUTINE TENSOR(DU1,DU2,DP,DNY,DT,KVERT,KMID,DCORVG,ELE)
c================================================================
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21)
      CHARACTER FMT*15,SUB*6,CPARAM*120
      DIMENSION DU1(*), DU2(*),DP(*),DNY(*),DT(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      common /fluidp/ pla,ple
      EXTERNAL DVISCO
      SAVE
C
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
C
      DO 100 IEL=1,NEL
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
      DO 120 IVE = 1,4
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
120   CONTINUE
      PXM=0.25D0*(DX(1)+DX(2)+DX(3)+DX(4))
      PYM=0.25D0*(DY(1)+DY(2)+DY(3)+DY(4))
C
      CALL ELE(0D0,0D0,-2)
      CALL ELE(PXM,PYM,-3)
C
      UH1 =0D0
      UH2 =0D0
      UH1X=0D0
      UH2X=0D0
      UH1Y=0D0
      UH2Y=0D0
C
      DO 210 JDOFE=1,IDFL
      IEQ=KDFG(JDOFE)
      ILO=KDFL(JDOFE)
      UH1 =UH1 +DU1(IEQ)*DBAS(ILO,1)
      UH2 =UH2 +DU2(IEQ)*DBAS(ILO,1)
      UH1X=UH1X+DU1(IEQ)*DBAS(ILO,2)
      UH2X=UH2X+DU2(IEQ)*DBAS(ILO,2)
      UH1Y=UH1Y+DU1(IEQ)*DBAS(ILO,3)
      UH2Y=UH2Y+DU2(IEQ)*DBAS(ILO,3)
C
210   CONTINUE
c
      DNORM= (UH1X**2
     *          +0.5d0*(UH2X+UH1Y)**2+UH2Y**2)**(0.5d0) 
      DVISCOSITY=DVISCO(DNORM)
      DT(IEL)=-2.0D0*DVISCOSITY*UH1X+DP(IEL)
C
100   CONTINUE
C
99999 END
c================================================================
      SUBROUTINE TENSOR1(DU1,DU2,DP,DNY,DT1,KVERT,KMID,DCORVG,ELE)
c================================================================
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21)
      CHARACTER FMT*15,SUB*6,CPARAM*120
      DIMENSION DU1(*), DU2(*),DP(*),DNY(*),DT1(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      common /fluidp/ pla,ple
      EXTERNAL DVISCO
      SAVE
C
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
C
      DO 100 IEL=1,NEL
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
      DO 120 IVE = 1,4
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
120   CONTINUE
      PXM=0.25D0*(DX(1)+DX(2)+DX(3)+DX(4))
      PYM=0.25D0*(DY(1)+DY(2)+DY(3)+DY(4))
C
      CALL ELE(0D0,0D0,-2)
      CALL ELE(PXM,PYM,-3)
C
      UH1 =0D0
      UH2 =0D0
      UH1X=0D0
      UH2X=0D0
      UH1Y=0D0
      UH2Y=0D0
C
      DO 210 JDOFE=1,IDFL
      IEQ=KDFG(JDOFE)
      ILO=KDFL(JDOFE)
      UH1 =UH1 +DU1(IEQ)*DBAS(ILO,1)
      UH2 =UH2 +DU2(IEQ)*DBAS(ILO,1)
      UH1X=UH1X+DU1(IEQ)*DBAS(ILO,2)
      UH2X=UH2X+DU2(IEQ)*DBAS(ILO,2)
      UH1Y=UH1Y+DU1(IEQ)*DBAS(ILO,3)
      UH2Y=UH2Y+DU2(IEQ)*DBAS(ILO,3)
C
210   CONTINUE
c
      DNORM= (UH1X**2
     *          +0.5d0*(UH2X+UH1Y)**2+UH2Y**2)**(0.5d0) 
      DVISCOSITY=DVISCO(DNORM)
      DT1(IEL)=2.0D0*DVISCOSITY*(0.5d0*(UH2X+UH1Y))

100   CONTINUE
C
99999 END
c================================================================
      SUBROUTINE TENSOR2(DU1,DU2,DP,DNY,DT2,KVERT,KMID,DCORVG,ELE)
c================================================================
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21)
      CHARACTER FMT*15,SUB*6,CPARAM*120
      DIMENSION DU1(*), DU2(*),DP(*),DNY(*),DT2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      common /fluidp/ pla,ple
      EXTERNAL DVISCO
      SAVE
C
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
C
      DO 100 IEL=1,NEL
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
      DO 120 IVE = 1,4
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
120   CONTINUE
      PXM=0.25D0*(DX(1)+DX(2)+DX(3)+DX(4))
      PYM=0.25D0*(DY(1)+DY(2)+DY(3)+DY(4))
C
      CALL ELE(0D0,0D0,-2)
      CALL ELE(PXM,PYM,-3)
C
      UH1 =0D0
      UH2 =0D0
      UH1X=0D0
      UH2X=0D0
      UH1Y=0D0
      UH2Y=0D0
C
      DO 210 JDOFE=1,IDFL
      IEQ=KDFG(JDOFE)
      ILO=KDFL(JDOFE)
      UH1 =UH1 +DU1(IEQ)*DBAS(ILO,1)
      UH2 =UH2 +DU2(IEQ)*DBAS(ILO,1)
      UH1X=UH1X+DU1(IEQ)*DBAS(ILO,2)
      UH2X=UH2X+DU2(IEQ)*DBAS(ILO,2)
      UH1Y=UH1Y+DU1(IEQ)*DBAS(ILO,3)
      UH2Y=UH2Y+DU2(IEQ)*DBAS(ILO,3)
C
210   CONTINUE
c
      DNORM= (UH1X**2
     *          +0.5d0*(UH2X+UH1Y)**2+UH2Y**2)**(0.5d0) 
      DVISCOSITY=DVISCO(DNORM)
      DT2(IEL)=-2.0D0*DVISCOSITY*UH2Y+DP(IEL)     
C
100   CONTINUE
C
99999 END
c================================================================
      SUBROUTINE BNBUILD(U1L1,U1L2,U2L1,U2L2,A1L,A2L,
     *                VB1,VB2,VBN1,VBN2,
     *                KCOLB,KLDB,NP,KVERT,KMID,KADJ,NB,DCORVG,ELE)
c================================================================
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21)
      CHARACTER FMT*15,SUB*6,CPARAM*120
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*)
      DIMENSION VB1(*),VB2(*),KCOLB(*),KLDB(*)
      DIMENSION VBN1(*),VBN2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS),KENETRY(NNVE)
      DIMENSION DB1ENTRYP(NNBAS),DBN1ENTRYP(NNBAS)
      DIMENSION KBENTRY(NNBAS) ,DB2ENTRYP(NNBAS),DBN2ENTRYP(NNBAS)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL

C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
      INCLUDE 'jump.inc'
      EXTERNAL E031,E010,COEFFB,ELE
      EXTERNAL DVISPO
      SAVE
C
      CALL LCL2(VB1,NB )! CLEAR 
      CALL LCL2(VB2,NB )
c
      CALL LCL2(VBN1,NB )! CLEAR 
      CALL LCL2(VBN2,NB )
C
      ICUB=ICUBN
C *** 
c
      DCONST=DNP
         IF(INEWPRE .EQ. 1 .AND. IUSENEWT .EQ. 1) DCONST=1.0
         IF(INEWPRE .EQ. 0) DCONST=0.0D0
C
      BDER(1)=.TRUE.
      BDER(2)=.TRUE.
      BDER(3)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
C
      DO 100 IEL=1,NEL
C
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
      DO 110 JDOFE=1,IDFL
      JCOL0=KLDB(KDFG(JDOFE))
      IDFG=IEL
      DO 112 JCOL=JCOL0,NB
      IF (KCOLB(JCOL).EQ.IDFG) GOTO 111
112   CONTINUE
111   JCOL0=JCOL+1
          KBENTRY(JDOFE)=JCOL 
          DB1ENTRYP(JDOFE)=0D0
          DB2ENTRYP(JDOFE)=0D0 
          DBN1ENTRYP(JDOFE)=0D0
          DBN2ENTRYP(JDOFE)=0D0 
110   CONTINUE
C
      DO 120 IVE = 1,4
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
120   CONTINUE
      DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
      DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
      DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
      DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
c
      PXM=0.25D0*(DX(1)+DX(2)+DX(3)+DX(4))
      PYM=0.25D0*(DY(1)+DY(2)+DY(3)+DY(4))
C
      CALL ELE(0D0,0D0,-2)

C
C *** Loop over all cubature points
      DO 200 ICUBP = 1, NCUBP
c
C *** Cubature points on the reference element
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
c
C *** Jacobian of the bilinear mapping onto the reference element
      DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
      DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
      DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
      DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
      DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
C
C *** Cubature points on actual Element + weights
      XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *  +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
      YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1
     *  +0.5D0*(DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2
      OM=DOMEGA(ICUBP)*DETJ
C
      CALL ELE(XX,YY,-3)
      IF (IER.LT.0) GOTO 99999
c
C *** Evaluation of velocity in cubature points
C
       DU1=0D0
       DU2=0D0
       IF (A2L.EQ.0D0) THEN
        DO 210 JDFL=1,IDFL
        HBAS=DBAS(KDFL(JDFL),1)
        IF (ABS(HBAS).GE.1D-8) THEN
         JDFG=KDFG(JDFL)
         DU1=DU1+U1L1(JDFG)*HBAS
         DU2=DU2+U1L2(JDFG)*HBAS
        ENDIF
210     CONTINUE
       ELSE
        DO 220 JDFL=1,IDFL
        HBAS=DBAS(KDFL(JDFL),1)
        IF (ABS(HBAS).GE.1D-8) THEN
         JDFG=KDFG(JDFL)
         DU1=DU1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
         DU2=DU2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
         ENDIF
220     CONTINUE
       ENDIF
C
C *** Evaluation of derivatives of velocity in cubature points
          DUX1=0D0
          DUX2=0D0
          IF (A2L.EQ.0D0) THEN
           DO 211 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),2)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUX1=DUX1+U1L1(JDFG)*HBAS
             DUX2=DUX2+U1L2(JDFG)*HBAS
            ENDIF
 211       CONTINUE
          ELSE
           DO 221 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),2)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUX1=DUX1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
             DUX2=DUX2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
            ENDIF
 221       CONTINUE
          ENDIF
          
          DUY1=0D0
          DUY2=0D0
          IF (A2L.EQ.0D0) THEN
           DO 212 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),3)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUY1=DUY1+U1L1(JDFG)*HBAS
             DUY2=DUY2+U1L2(JDFG)*HBAS
            ENDIF
 212       CONTINUE
          ELSE
           DO 222 JDFL=1,IDFL
            HBAS=DBAS(KDFL(JDFL),3)
            IF (ABS(HBAS).GE.1D-8) THEN
             JDFG=KDFG(JDFL)
             DUY1=DUY1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
             DUY2=DUY2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
            ENDIF
 222       CONTINUE
          ENDIF
C
C ***  norm of symmetric part of velocity gradient
       IF (ISPGRAD .EQ. 1)DUGSQ=DUX1**2+0.5d0*(DUX2+DUY1)**2+DUY2**2
c
C ***  norm of gradient
       IF (ISPGRAD .EQ. 0)DUGSQ=DUX1**2+DUX2**2+DUY1**2+DUY2**2
C
       IF (ISPGRAD .EQ. 1)DNY=2.0D0*DVISPO(DUGSQ)
       IF (ISPGRAD .EQ. 0)DNY=DVISPO(DUGSQ)
C
          IF(ISPGRAD .EQ. 1)THEN 
             HBN11=DUX1
             HBN12=0.5D0*(DUX2+DUY1)
             HBN21=0.5D0*(DUX2+DUY1)
             HBN22=DUY2
          ELSE
             HBN11=DUX1
             HBN12=DUY1
             HBN21=DUX2
             HBN22=DUY2
          ENDIF
C     
C ***  Summing up over all pairs of multiindices
       DO 230 JDOFE=1,IDFL
          JDOFEH=KDFL(JDOFE)
          HBASJ1=DBAS(JDOFEH,1)
          HBASJ2=DBAS(JDOFEH,2)
          HBASJ3=DBAS(JDOFEH,3)
C
          DH1=DNY*(HBASJ2*HBN11+HBASJ3*HBN12)*DCONST
          DH2=DNY*(HBASJ2*HBN21+HBASJ3*HBN22)*DCONST    
C     
          DB1ENTRYP(JDOFE)=DB1ENTRYP(JDOFE)-OM*HBASJ2
          DB2ENTRYP(JDOFE)=DB2ENTRYP(JDOFE)-OM*HBASJ3
          DBN1ENTRYP(JDOFE)=DBN1ENTRYP(JDOFE)-OM*HBASJ2+OM*DH1
          DBN2ENTRYP(JDOFE)=DBN2ENTRYP(JDOFE)-OM*HBASJ3+OM*DH2
230    CONTINUE
c

200   CONTINUE
C
      DO 400 JDOFE=1,IDFL
         IB=KBENTRY(JDOFE)
         JDFG=KDFG(JDOFE)
         JDL=KLDB(KDFG(JDOFE))
         JCOL=KCOLB(KDFG(JDOFE))
         VB1(IB)=VB1(IB)+REAL(DB1ENTRYP(JDOFE))
         VB2(IB)=VB2(IB)+REAL(DB2ENTRYP(JDOFE))
         VBN1(IB)=VBN1(IB)+REAL(DBN1ENTRYP(JDOFE))
         VBN2(IB)=VBN2(IB)+REAL(DBN2ENTRYP(JDOFE))
400   CONTINUE
C
100   CONTINUE
c
C
99999 END
C
************************************************************************
      SUBROUTINE  VANCR2 (U1,U2,P,F1,F2,FP,
     *                    A,KACOL,KALD,B1,B2,BN1,BN2,KBCOL,KBLD,
     *                    KMBD,KVERT,KMID,DCORVG,KNPR,NMBD,
     *                    Inum,AKTELE,ite)
************************************************************************
C
C----------------------------------------------------------------------
C Purpose:  Vanca-type Glaetter. 
C     Auf jedem Element wird die lokale 
c     8x8-Matrix (AA(i,j)) EXAKT invertiert
c     
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'
      INCLUDE 'bouss.inc'
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION BB1,BB2,BBN1,BBN2
      REAL  A,B1,B2,BN1,BN2

      INTEGER AKTELE
C
c      PARAMETER (NNVE=4)
      DIMENSION U1(*),U2(*),P(*)
      DIMENSION F1(*),F2(*),FP(*),BN1(*),BN2(*)
      DIMENSION A(*),KACOL(*),KALD(*),B1(*),B2(*),KBCOL(*),KBLD(*)
      DIMENSION KMBD(*),KVERT(NNVE,*),KMID(NNVE,*),KNPR(*),DCORVG(2,*)
C
C *** Local arrays for informations about one element
      DIMENSION AA(8,8),BB1(4),BB2(4),FF(8),BBN1(4),BBN2(4)
      DIMENSION IU(4),BDBC(4),UV(8)
C
      NA=KALD(NU+1)-1
C
C=======================================================================
C     Block Gauss-Seidel on Schur Complement
C=======================================================================
c$$$C=======================================================================
c$$$      IF ((ILEV.NE.NLMAX)) THEN
c$$$          CALL GARLI3 (ILEV,NEL,NVT,NMT,
c$$$     *             DWORK(L(LCORVG)), KWORK(L(LVERT)),KWORK(L(LMID)),
c$$$     *            f1,f2,20+ite)
c$$$          ENDIF
C=======================================================================
c
      IF (INUM.eq.0) THEN
         Nloop=1
      ELSE
         NLOOP=NEL
      ENDIF

      CALL VKADD (AVK1,AVK2)
      ALPHM=0d0

      qumax=0d0
      qumin=1d10
      qumax2=0d0
      qumin2=1d10
      DO 1  IELN=1,NLOOP!NEL
C
      IF (INUM.eq.0) THEN
       IEL=AKTELE
      ELSE
       IEL=IELN
      ENDIF
C
      FFP=FP(IEL)
C
C
      DO 10  II=1,4
      IMID=KMID(II,IEL)
      I=IMID-NVT
      IU(II)=I
10    CONTINUE
C
C
c

C-----------------------------------------------------------------------
C *** Loop over all 4 U-nodes of that element
      DO 11  II=1,4
C
      IMID=KMID(II,IEL)
      I=IMID-NVT
      IF (KNPR(IMID).EQ.0)  THEN
       BDBC(II)=.FALSE.
      ELSE
       BDBC(II)=.TRUE.
      ENDIF
C
C *** Put on AA(.) the diagonal entry of matrix A
      IA1=KALD(I)
      AA(II,II)=AVK1*DBLE(A(IA1))+AVK2*DBLE(A(4*NA+ia1))
      AA(II,II+4)=AVK1*DBLE(A(NA+IA1))!*0d0
      AA(II+4,II)=AVK1*DBLE(A(2*NA+IA1))!*0d0
      AA(II+4,II+4)=AVK1*DBLE(A(3*NA+IA1))+AVK2*DBLE(A(4*NA+ia1))
C

C *** Initial setting of FF (formerly known as FF1,FF2)
      FF(II)=F1(I)
      FF(4+II)=F2(I)
C
      IF (BDBC(II)) GOTO 13
C
C *** Put all off-diagonal entries of A to the right hand side
ccc      GOTO 111
      IA2=KALD(I+1)-1
      DO 110  IA=IA1+1,IA2
      J=KACOL(IA)
      IF ((J.NE.IU(1)).AND.(J.NE.IU(2)).AND.
     *    (J.NE.IU(3)).AND.(J.NE.IU(4))) THEN!(behaelt 4 Komponenten)
      AOFF1=AVK1*DBLE(A(IA))+AVK2*DBLE(A(4*NA+ia))
      AOFF2=AVK1*DBLE(A(NA+IA))
      AOFF3=AVK1*DBLE(A(2*NA+IA))
      AOFF4=AVK1*DBLE(A(3*NA+IA))+AVK2*DBLE(A(4*NA+ia))
C
c      write (*,*) ia,'aoff', aoff1,aoff2
      FF(II)  =FF(II)  -AOFF1*U1(J)-AOFF2*U2(j)
      FF(4+II)=FF(4+II)-AOFF3*U1(J)-AOFF4*U2(J)
      ENDIF
110   CONTINUE
C
C-----------------------------------------------------------------------
C *** Get BB1,BB2 and modify FF1,FF2
      IB1=KBLD(I)
      IB2=KBLD(I+1)-1
      JP1=KBCOL(IB1)
      JP2=KBCOL(IB2)
C
      IF (JP1.EQ.IEL)  THEN
          PJP=P(JP2)
          BB1(II)=DBLE(B1(IB1))
          BB2(II)=DBLE(B2(IB1))
          BBN1(II)=DBLE(BN1(IB1))
          BBN2(II)=DBLE(BN2(IB1))
          FF(II)  =FF(II)  -DBLE(BN1(IB2))*PJP
          FF(4+II)=FF(4+II)-DBLE(BN2(IB2))*PJP
      ELSE IF (JP2.EQ.IEL)  THEN
          PJP=P(JP1)
          BB1(II)=DBLE(B1(IB2))
          BB2(II)=DBLE(B2(IB2))
          BBN1(II)=DBLE(BN1(IB2))
          BBN2(II)=DBLE(BN2(IB2))
          FF(II)  =FF(II)  -DBLE(BN1(IB1))*PJP
          FF(4+II)=FF(4+II)-DBLE(BN2(IB1))*PJP
      ELSE
          WRITE(MTERM,*) 'ERROR in SMOOTH: IEL entry in B not found'
          STOP
      ENDIF
C
      GOTO 11
C
C-----------------------------------------------------------------------
C *** The case that II is a Dirichlet boundary node
13    CONTINUE
      IB=KBLD(I)
      JP=KBCOL(IB)
      BB1(II)=DBLE(B1(IB))
      BB2(II)=DBLE(B2(IB))
      BBN1(II)=DBLE(BN1(IB))
      BBN2(II)=DBLE(BN2(IB))
C
C
11    CONTINUE
C
C-----------------------------------------------------------------------
      DO 21  II=1,4 
      I=IU(II)
C
      DO 22  JJ=1,4
      IF (JJ.EQ.II) GOTO 22
C
      AA(II,JJ)=0D0
      AA(II+4,JJ)=0D0
      AA(II,JJ+4)=0D0
      AA(4+II,4+JJ)=0D0
      IF (BDBC(II)) GOTO 22
C
      J=IU(JJ)
      IA1=KALD(I)
      IA2=KALD(I+1)-1
C
      DO 220  IA=IA1+1,IA2
      JH=KACOL(IA)
      IF (J.EQ.JH) THEN
      AA(II,JJ)=AVK1*DBLE(A(IA))+AVK2*DBLE(A(4*NA+ia))
      AA(II,JJ+4)=AVK1*DBLE(A(NA+IA))!*0d0
      AA(II+4,JJ)=AVK1*DBLE(A(2*NA+IA))!*0d0
      AA(II+4,JJ+4)=AVK1*DBLE(A(3*NA+IA))+AVK2*DBLE(A(4*NA+ia))
C
c       AA(jj,ii)=DBLE(A(IA))
c       AA(II,JJ)    =DBLE(A(     IA))
c       AA(II,4+JJ)  =DBLE(A(  NA+IA))
c       AA(4+II,JJ)  =DBLE(A(2*NA+IA))
c       AA(4+II,4+JJ)=DBLE(A(3*NA+IA))
       GOTO 22
      ENDIF
 220  CONTINUE
C
 22   CONTINUE
C
 21   CONTINUE
C


C-----------------------------------------------------------------------
C *** Calling the subroutine ELUPDN for element update 
      CALL  ELUPN2 (U1,U2,P,IEL,IU,BDBC,AA,BB1,BB2,FF,FFP)
c$$$      CALL  ELUPN2 (UV,P,IEL,IU,BDBC,AA,BB1,BB2,FF,FFP)
C

c$$$c      IF (ilev.ne.nlmax) then
c$$$      CALL GARLIC (IEL,ILEV,NVT,DCORVG, KVERT,
c$$$     *             KMID, U1,U2, Dul2,dum,dalfa,dquot, 
c$$$     *             dwork(l(lfil1)),dwork(l(lfil2)),0)
c$$$
c      endif
1     CONTINUE
c$$$C=======================================================================
c$$$          CALL GARLI3 (ILEV,NEL,NVT,NMT,
c$$$     *             DWORK(L(LCORVG)), KWORK(L(LVERT)),KWORK(L(LMID)),
c$$$     *            U1,U2,50+ite)
C=======================================================================

C=======================================================================
C
      END
c
************************************************************************
      SUBROUTINE  shift (DX,NEQ)  
************************************************************************
*
*   Purpose: - shift the pressure 
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
C-----------------------------------------------------------------------
C
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C=======================================================================
C     Getting all parameters for DX
C=======================================================================
C *** addresses for the current level ILEV
C
      ISETLV=2
      CALL  SETLEV (ISETLV)
      I2=1+NU
      IP=I2+NU
C
C
      IF(ILEV .EQ. NLEV)THEN
         Do I=1,NP
            DX(IP-1+I)=DX(IP-1+I)+0.0D0
         ENDDO
      ENDIF 
C
C
      END
C
c$$$c==================================================================
c$$$      subroutine mext7p1(neq,kcol,kld,na,na1,iaux)
c$$$c==================================================================
c$$$c computes the new/extended number of nozeros (=NA1)
c$$$c needs working integer vector iaux with dimension neq
c$$$      
c$$$      integer neq,na,na1,iaux(*),kcol(*),kld(*)
c$$$      integer i,j,k,i1,j1,k1,k2,isum,ia
c$$$
c$$$      isum=1
c$$$      do i=1,neq
c$$$        ia=1
c$$$        
c$$$        do k=kld(i),kld(i+1)-1
c$$$          j=kcol(k)
c$$$          iaux(ia)=j
c$$$          ia=ia+1
c$$$          isum=isum+1
c$$$        enddo
c$$$        
c$$$        do k=kld(i),kld(i+1)-1
c$$$          j=kcol(k)
c$$$
c$$$          do k1=kld(j),kld(j+1)-1
c$$$            j1=kcol(k1)
c$$$            iflag=0
c$$$
c$$$            do k2=1,ia-1
c$$$              j2=iaux(k2)
c$$$              if(j1.eq.j2) then
c$$$                iflag=iflag+1
c$$$              endif
c$$$            enddo
c$$$            
c$$$            if(.not.(iflag.eq.0).or.(iflag.eq.1)) then
c$$$              print *,'error'
c$$$            endif
c$$$
c$$$            if(iflag.eq.0) then
c$$$              iaux(ia)=j1
c$$$              ia=ia+1
c$$$              isum=isum+1
c$$$            endif
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$
c$$$      na1=isum-1
c$$$
c$$$      end
c$$$
c$$$      subroutine mext7p2(neq,kcol,kld,na,kcol1,kld1,na1)
c$$$c create the extended kld/kcol vectors
c$$$c needs na1, kcol1(na1), kld(neq)
c$$$
c$$$      
c$$$      integer neq,na,na1,kcol(*),kld(*),kcol1(*),kld1(*)
c$$$      integer i,j,k,i1,j1,k1,k2,isum,ia
c$$$
c$$$      isum=1
c$$$      do i=1,neq
c$$$        ia=1
c$$$        kld1(i)=isum
c$$$
c$$$        do k=kld(i),kld(i+1)-1
c$$$          j=kcol(k)
c$$$          kcol1(isum)=j
c$$$          ia=ia+1
c$$$          isum=isum+1
c$$$        enddo
c$$$          
c$$$        do k=kld(i),kld(i+1)-1
c$$$          j=kcol(k)
c$$$
c$$$          do k1=kld(j),kld(j+1)-1
c$$$            j1=kcol(k1)
c$$$            iflag=0
c$$$
c$$$            do k2=kld1(i),ia-1
c$$$              j2=kcol1(k2)
c$$$              if(j1.eq.j2) then
c$$$                iflag=iflag+1
c$$$              endif
c$$$            enddo
c$$$            
c$$$            if(.not.(iflag.eq.0).or.(iflag.eq.1)) then
c$$$              print *,'error'
c$$$            endif
c$$$
c$$$            if(iflag.eq.0) then
c$$$              kcol1(isum)=j1
c$$$              ia=ia+1
c$$$              isum=isum+1
c$$$            endif
c$$$          enddo
c$$$        enddo
c$$$      enddo
c$$$      
c$$$      na1=isum-1
c$$$      lda(neq+1)=isum
c$$$
c$$$      end

************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* XAP7                                                                 *
*                                                                      *
* Purpose  Call AP7                                                    *
*          Allocate KLD and KCOL on DWORK                              *
*                                                                      *
* Subroutines/functions called  AP7, EA00, ZNEW, ZDISP, ZFREE          *
*                                                                      *
* Version from  12/11/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* ELE      SUBR   EXTERNAL Subroutine - evaluation of basis functions  *
*                 for dummy call only                                  *
* For the description of the remaining input parameter see AP7         *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* LCOL     I*4    Number of KCOL                                       *
* LLD      I*4    Number of KLD                                        *
* NA       I*4    Length of KCOL (and VA (DA))                         *
* NEQ      I*4    Number of equations                                  *
*                                                                      *
************************************************************************
*                                                                      *
* AP7                                                                  *
*                                                                      *
* Purpose  Calculation of the pointer vectors KLD and KCOL             *
*          for a matrix corresponding to a given element type          *
*          Storage technique 7/8                                       *
*                                                                      *
* Subroutines/functions called  NDFL, NDFGL, EA00                      *
*                                                                      *
* Version from  12/11/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* KCOL1    I*4    Auxiliary array                                      *
* KIND     I*4    Auxiliary array                                      *
* NA       I*4    Maximal length of KCOL                               *
* NEQ      I*4    Number of equations                                  *
* ELE      SUBR   EXTERNAL Subroutine - evaluation of basis functions- *
*                 for dummy call only                                  *
* ISYMM    I*4    >=1 matrix symmetric - storage technique 8           *
* KVERT    I*4    Arrays describing the triangulation                  *
* KMID     I*4                                                         *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KCOL     I*4    Pointer vector containing column indices             *
* KLD      I*4    Pointer vector to diagonal elements                  *
* NA       I*4    Effective number of elements in KCOL                 *
* IER      I*4    Error indicator                                      *
*                 -118  Not enough space for KCOL                      *
*                       Error occured on element IEL                   *
*                                                                      *
************************************************************************
C
      SUBROUTINE XAPP7(LCOL,LLD,NA,NEQ,ELE,ISYMM)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL ELE
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/
C
      SUB='XAP7'
      IF (ICHECK.GE.997) CALL OTRC('XAP7  ','12/11/89')
C
C *** Determine total number of degrees of freedom
      CALL EA00(ELE,ELE,NVE,IELTYP)
      NEQ=NDFG(IELTYP)
      IF (IER.NE.0) GOTO 99999
C
C *** Allocate KLD on DWORK ***
      CALL ZNEW(NEQ+1,-3,LLD,'KLD   ')
      IF (IER.NE.0) GOTO 99999
C *** Allocate KLD1 on DWORK ***
      CALL ZNEW(NEQ+1,-3,LLD1,'KLD1   ')
      IF (IER.NE.0) GOTO 99999
C
C *** Determine free space on DWORK 
      IWMAX0=IWMAX
      CALL ZFREE(3,IFREE)
      IF (IER.NE.0) GOTO 99999
C
      NA=IFREE/4-4
      CALL ZNEW(NA,-3,LCOL,'KCOL  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NA,-3,LCOL1,'KCOL1 ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NA,-3,LIND,'KIND  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NEQ,-3,laux,'iaux  ')
      IF (IER.NE.0) GOTO 99999

C
      CALL APP7(KWORK(L(LCOL)),KWORK(L(LCOL1)),
     *         KWORK(L(LIND)),KWORK(L(LLD)),
     *         NA,NEQ,ELE,ISYMM, 
     *         KWORK(L(LVERT)),KWORK(L(LMID)))
      IF (IER.NE.0) GOTO 99999
C
      call apex71(neq,kwork(l(lcol)),kwork(l(lld)),na,
     *            na1,
     *            kwork(l(laux)))
      PRINT*,'NA,NA1',NA,NA1
c
C *** Release space on DWORK not used for KCOL ***
      CALL ZDISP(NA,LIND,'KIND  ')
      CALL ZDISP(NA,LCOL1,'KCOL1 ')
      CALL ZDISP(NA,LCOL,'KCOL  ')
C
      IWMAX=MAX(IWORK,IWMAX0)
C 
      CALL ZDISP(0,LIND,'KIND  ')
      CALL ZDISP(0,LCOL1,'KCOL1 ')
      CALL ZNEW(na1,-3,LCOL1,'KCOL1 ')
      IF (IER.NE.0) GOTO 99999
      
C      CALL ZNEW(NEQ+1,-3,LLD1,'KLD1  ')
C      IF (IER.NE.0) GOTO 99999

      call apex72(neq,kwork(l(lcol)),kwork(l(lld)),na,
     *            kwork(l(lcol1)),kwork(l(lld1)),na1,
     *            kwork(l(laux)))
      CALL ZNEW(na1,-3,LCOL,'KCOL ')
      na=na1
       CALL XLCP3(lcol1,lcol,Na1)
       CALL XLCP3(lld1,lld,Neq+1)

      CALL ZDISP(0,laux,'iaux  ')
99999 END
C
C
C
      SUBROUTINE APP7(KCOL,KCOL1,KIND,KLD,NA,NEQ,ELE,ISYMM,
     *                                              KVERT,KMID)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=21,NNVE=4)
      DIMENSION KLD(*),KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KDFG1(NNBAS),KDFL1(NNBAS)
      DIMENSION KCOL(*),KCOL1(*),KIND(*),KVERT(NNVE,*),KMID(NNVE,*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      EXTERNAL ELE
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='AP7'
      IF (ICHECK.GE.997) CALL OTRC('AP7   ','12/11/89')
C
      IER=0
C
      IF (NA.LT.NEQ) THEN
C *** DWORK exhausted
       WRITE (CPARAM,'(I15)') 1
       CALL WERR(-118,'AP7   ')
       GOTO 99999
      ENDIF
      IFREE=NA-NEQ
      NA=NEQ
C *** Initialization of KIND and KCOL1
      DO 10 IEQ=1,NEQ
      KIND(IEQ)=0
10    KCOL1(IEQ)=IEQ
C
C *** Set element number by dummy call ***
      CALL EA00(ELE,ELE,NVE,IELTYP)
C *** Determine number of degrees of freedom per element ***
      IDFL=NDFL(IELTYP)
      IF (IER.NE.0) GOTO 99999
C
      BSYMM=ISYMM.EQ.1
      IDOFE1=1
C 
C *** Loop over elements
      DO 100 IEL=1,NEL
C
C *** KDFG returns the global degrees of freedom in increasing order
      CALL NDFGL(IEL,0,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.NE.0) GOTO 99999
C
C *** Loop over local number of degrees of freedom
      DO 110 IDOFE=1,IDFL
      IROW=KDFG(IDOFE)
      IF (BSYMM) THEN
       IF (IROW.EQ.NEQ) GOTO 110
       IDOFE1=IDOFE
      ENDIF
C
C *** Loop over off-diagonal elements
C *** only upper triangular part in symmetric case
      DO 120 JDOFE=IDOFE1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 120
      JCOL=KDFG(JDOFE)
C *** JCOL is the global number of d.o.f. to be inserted into row IROW
C
C *** Look whether entry (IROW,JCOL) is already provided
      IPOS=IROW
121   IF (KIND(IPOS).EQ.0) THEN
C *** final entry in row IROW up to now
C
      IF (IFREE.LE.0) THEN
C *** DWORK exhausted
        WRITE (CPARAM,'(I15)') IEL
        CALL WERR(-118,'AP7   ')
        GOTO 99999
       ENDIF
C
C *** insert new entry
       NA=NA+1
       IFREE=IFREE-1
       KCOL1(NA)=JCOL
       KIND(IPOS)=NA
       KIND(NA)=0
       GOTO 120
C
      ELSE
C
       IPOS=KIND(IPOS)
       IF (KCOL1(IPOS).NE.JCOL) GOTO 121
C
      ENDIF
C
120   CONTINUE
c
c$$$c***  Extended of the matrix
c$$$c
c$$$      IELADJ=KMEL(2,IROW)
c$$$      IF(IELADJ .NE. 0)THEN
c$$$      CALL NDFGL(IELADJ,0,IELTYP,KVERT,KMID,KDFG1,KDFL)
c$$$      IF (IER.NE.0) GOTO 99999
c$$$      DO 1201 JDOFE=IDOFE1,IDFL
c$$$      IF (IDOFE.EQ.JDOFE) GOTO 1201
c$$$      JCOL=KDFG1(JDOFE)
c$$$C *** JCOL is the global number of d.o.f. to be inserted into row IROW
c$$$C
c$$$C *** Look whether entry (IROW,JCOL) is already provided
c$$$      IPOS=IROW
c$$$ 1211 IF (KIND(IPOS).EQ.0) THEN
c$$$C *** final entry in row IROW up to now
c$$$C
c$$$      IF (IFREE.LE.0) THEN
c$$$C *** DWORK exhausted
c$$$        WRITE (CPARAM,'(I15)') IELADJ
c$$$        CALL WERR(-118,'AP7   ')
c$$$        GOTO 99999
c$$$       ENDIF
c$$$C
c$$$C *** insert new entry
c$$$       NA=NA+1
c$$$       IFREE=IFREE-1
c$$$       KCOL1(NA)=JCOL
c$$$       KIND(IPOS)=NA
c$$$       KIND(NA)=0
c$$$       GOTO 1201
c$$$C
c$$$      ELSE
c$$$C
c$$$       IPOS=KIND(IPOS)
c$$$       IF (KCOL1(IPOS).NE.JCOL) GOTO 1211
c$$$C
c$$$      ENDIF
c$$$C
c$$$ 1201 CONTINUE
c$$$      ENDIF
110   CONTINUE
100   CONTINUE
C
C
C *** Collect entries on KCOL1 separately for each row
C
      NA=0
      DO 200 IEQ=1,NEQ
      NA=NA+1
      KLD(IEQ)=NA
      KCOL(NA)=IEQ
      IPOS=IEQ
C
201   IF (KIND(IPOS).NE.0) THEN
       NA=NA+1
       IPOS=KIND(IPOS)
       KCOL(NA)=KCOL1(IPOS)
       GOTO 201
      ENDIF
C
200   CONTINUE
      KLD(NEQ+1)=NA+1
C
C
C *** Sort off-diagonal entries on KCOL separately for each row
C
      DO 300 IEQ=1,NEQ
C
301   BSORT=.TRUE.
      DO 302 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-2
      IF (KCOL(ICOL).GT.KCOL(ICOL+1)) THEN
       IHELP=KCOL(ICOL)
       KCOL(ICOL)=KCOL(ICOL+1)
       KCOL(ICOL+1)=IHELP
       BSORT=.FALSE.
      ENDIF
302   CONTINUE
      IF (.NOT.BSORT) GOTO 301
C
300   CONTINUE
C
99999 END
      SUBROUTINE XMAPP7(KLCOL,KLLD,KNA,KNEQ,KMEL,ELE,ISYMM)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNARR=299,NNLEV=9)
      DIMENSION KLCOL(*),KLLD(*),KNA(*),KNEQ(*),KMEL(2,*)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL ELE
      SAVE /TRIAA/,/TRIAD/,/MGTRD/,/MGTRA/,/MGPAR/,/ERRCTL/,/CHAR/
C
      SUB='XMAPP7 '
      IF (ICHECK.GE.997) CALL OTRC('XMAP7 ','08/25/90')
      IER=0
C
      DO 10 ILEV=NLMIN,NLMAX
C
      NEL =KNEL(ILEV)
      NVT =KNVT(ILEV)
      NMT =KNMT(ILEV)
      NVEL=KNVEL(ILEV)
      NVBD=KNVBD(ILEV)
C
      LCORVG=KLCVG(ILEV)
      LCORMG=KLCMG(ILEV)
      LVERT =KLVERT(ILEV)
      LMID  =KLMID(ILEV)
      LADJ  =KLADJ(ILEV)
      LVEL  =KLVEL(ILEV)
      LMEL  =KLMEL(ILEV)
      LNPR  =KLNPR(ILEV)
      LMM   =KLMM(ILEV)
      LVBD  =KLVBD(ILEV)
      LEBD  =KLEBD(ILEV)
      LBCT  =KLBCT(ILEV)
      LVBDP =KLVBDP(ILEV)
      LMBDP =KLMBDP(ILEV)
C
      CALL XAPP7(KLCOL(ILEV),KLLD(ILEV),KNA(ILEV),KNEQ(ILEV),ELE,ISYMM)
      IF (IER.NE.0) GOTO 99999
C     
10    CONTINUE
C     
99999 END
      
************************************************************************
*                                                                      *
* XAPEX7                                                               *
*                                                                      *
* Purpose  Allocate and find the extended LCOL1, LLD1 arrays           *
*          Call of APEX71 APEX72                                       *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LCOL     I*4    Numbers of the pointer arrays KCOL and KLD           *
* LLD      I*4    calculated by AP7                                    *
* NA       I*4    Number of nonzeros                                   *
* NEQ      I*4    Number of equations (=number of matrix rows)         *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* LCOL1    I*4    Numbers of the pointer arrays KCOL1 and KLD1         *
* LLD1     I*4    of the matrix structure extended by 1 level          *
* NA1      I*4    Number of nonzeros in the new structure              *
*                                                                      *
************************************************************************

      subroutine xapex7(LCOL,LLD,NA,NEQ,lcol1,lld1,na1)

      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))

      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      
      SUB='XAPEX7'
      IF (ICHECK.GE.997) CALL OTRC('XAPEX7','12/01/05')

      ! allocate the working vector iaux()
      CALL ZNEW(NEQ,-3,laux,'iaux  ')
      IF (IER.NE.0) GOTO 99999
      
      ! find the new number of nonzeros na1
      call apex71(neq,kwork(l(lcol)),kwork(l(lld)),na,
     *            na1,
     *            kwork(l(laux)))

      ! allocate the arrays for the new structure
      CALL ZNEW(na1,-3,LCOL1,'KCOL1 ')
      IF (IER.NE.0) GOTO 99999
      
      CALL ZNEW(NEQ+1,-3,LLD1,'KLD1  ')
      IF (IER.NE.0) GOTO 99999

      ! now fill in the new structure
      call apex72(neq,kwork(l(lcol)),kwork(l(lld)),na,
     *            kwork(l(lcol1)),kwork(l(lld1)),na1,
     *            kwork(l(laux)))

      CALL ZDISP(0,laux,'iaux  ')
99999 END      


c======================================================================


      subroutine apex71(neq,kcol,kld,na,na1,iaux)
c
c Matrix Pointer EXtend format 7 pass 1
c
c returns the new/extended number of nozeros (=NA1)
c needs working integer vector iaux with dimension neq
c      
      integer neq,na,na1,iaux(*),kcol(*),kld(*)
      integer i,j,k,i1,j1,k1,k2,isum,ia

      isum=0
      
      do i=1,neq ! loop over all eq (=mat. rows)

        ! clear the vector iaux(.)
        call LCL3(iaux,neq)

        ! create the matrix row image in the vector iaux(.)
        ! iaux(i)=0/1 if zero/nonzero
        do k=kld(i),kld(i+1)-1
          j=kcol(k)
          iaux(j)=1
        enddo
        
        ! now extend the the row image by including coupling 
        ! from the other rows like
        ! for given row i, all j and all j1 do
        ! if a(i,j) is nonzero and a(j,j1) is nonzero 
        ! then set a(i,j1) to nonzero too
        do k=kld(i)+1,kld(i+1)-1
          j=kcol(k)

          do k1=kld(j),kld(j+1)-1
            j1=kcol(k1)
            iaux(j1)=1
          enddo
        enddo

        ! count the nonzeros in the extended row
        ia=0
        do j=1,neq
          if(iaux(j).eq.1) ia=ia+1
        enddo

        isum=isum+ia
      enddo ! loop over all eq (=mat. rows)

      na1=isum
      end
c======================================================================

      subroutine apex72(neq,kcol,kld,na,kcol1,kld1,na1,iaux)
c
c Matrix Pointer EXtend format 7 pass 2
c
c create the extended kld/kcol vectors
c needs na1 (computed by a previous call to mext7p1), 
c kcol(na), kld(neq) ... description of the old matrix
c kcol1(na1), kld1(neq) ... allocated arrays for the new structure
c needs working integer vector iaux with dimension neq
c

      integer neq,na,na1,kcol(*),kld(*),kcol1(*),kld1(*)
      integer i,j,k,i1,j1,k1,k2,isum,iaux(*)

      icol1=1

      do i=1,neq ! loop over all eq (=mat. rows)

        ! clear the vector iaux(.)
        call LCL3(iaux,neq)

        ! create the matrix row image in the vector iaux(.)
        ! iaux(i)=0/1 if zero/nonzero
        do k=kld(i),kld(i+1)-1
          j=kcol(k)
          iaux(j)=1
        enddo
        
        ! now extend the the row image by including coupling 
        ! from the other rows like
        ! a(i,j) is nonzero and a(j,j1) is nonzero 
        ! then set a(i,j1) to nonzero too
        do k=kld(i)+1,kld(i+1)-1
          j=kcol(k)

          do k1=kld(j),kld(j+1)-1
            j1=kcol(k1)
            iaux(j1)=1
          enddo
        enddo

        kld1(i)=icol1

        ! do the diagonal first
        if(iaux(i).eq.0) print *,'error: no diagonal!'
        kcol1(icol1)=i
        iaux(i)=0
        icol1=icol1+1

        ! extract the remaining nonzeros to the new kcol1 
        ! automaticaly in the correct ordering
        do j=1,neq
          if(iaux(j).eq.1) then
            kcol1(icol1)=j
            icol1=icol1+1
          endif
        enddo

      enddo ! loop over all eq (=mat. rows)

      kld1(neq+1)=icol1
      if(na1.ne.icol1-1) print *,'error: wrong count of nnz'

      end
      
