************************************************************************
      SUBROUTINE   INTUVD (DU1,DU2,DT,DL1,DL2,DL3,DAUX,NVT,NEL,NVBD,
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
      PARAMETER (NNVE=4,NNLEV=9)
C
      DIMENSION DU1(*),DU2(*),DT(*),DL1(*),DL2(*),DL3(*),DAUX(*),
     *          KVERT(NNVE,*),KMID(NNVE,*),KVBD(*),KMBD(*),DCORVG(2,*)
C
c      INCLUDE 'bouss.inc'
      SAVE
C-----------------------------------------------------------------------
      DO 1 IVT=1,NVT
      DL1 (IVT)=0D0
      DL2 (IVT)=0D0
      DL3 (IVT)=0D0
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
      DTH1=DT(IM1)
      DTH2=DT(IM2)
      DTH3=DT(IM3)
      DTH4=DT(IM4)
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
      DL3(IV1)=DL3(IV1) + 0.75D0*(DTH1+DTH4) - 0.25D0*(DTH2+DTH3)
      DL3(IV2)=DL3(IV2) + 0.75D0*(DTH2+DTH1) - 0.25D0*(DTH3+DTH4)
      DL3(IV3)=DL3(IV3) + 0.75D0*(DTH3+DTH2) - 0.25D0*(DTH4+DTH1)
      DL3(IV4)=DL3(IV4) + 0.75D0*(DTH4+DTH3) - 0.25D0*(DTH1+DTH2)
C
10    CONTINUE
C
      DO 20  IV=1,NVT
      DL1(IV)=DL1(IV)/DAUX(IV)
      DL2(IV)=DL2(IV)/DAUX(IV)
      DL3(IV)=DL3(IV)/DAUX(IV)
20    CONTINUE
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
C
      BDFLAG=.TRUE.
      CALL TRAND (KNPR,DDBD,NVT,NMBD,BDFLAG,IMID,IMBD,0,X,Y)
C
      IF (BDFLAG) THEN
         DL3(IV)=UE(X,Y,3)
      ENDIF
30    CONTINUE
C
      END
c
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
      DIMENSION  AREA(*),KVERT(NNVE,*),DCORVG(2,*)
c
C=======================================================================
c
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
*    Purpose: - Transforms the vector P into the space L2_0
*             - uses the vector AREA with the areas of all elements
*               and on AREA(NEL+1) the sum of all AREA(IEL)
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL AREA
      DIMENSION P(*),AREA(*)
C
      IF (INEUM.EQ.1) RETURN
C
      PINT=0.D0
      DO 2  IEL=1,NEL
  2   PINT=PINT+P(IEL)*DBLE(AREA(IEL))
C
      C=PINT/DBLE(AREA(NEL+1))
      DO 3  IEL=1,NEL
  3   P(IEL)=P(IEL)-C
C
      END
c
c
c
************************************************************************
      SUBROUTINE   S E T L E V   (ISETLV)
************************************************************************
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
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
C *** global constants
      PARAMETER (NNARR=299,NNLEV=9,  NNWORK=1)
C
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
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLMASS(NNLEV),KLM(NNLEV),
     *                KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGDIM/  KNA(NNLEV),KNB(NNLEV),KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
C
C
      INCLUDE 'bouss.inc'
      SAVE 
C
C-----------------------------------------------------------------------
C
      IF (ILEV.LT.NLMIN .OR. ILEV.GT.NLMAX .OR. ISETLV.LT.1  .OR.
     *    ISETLV.GT.2 )  THEN
        WRITE(MTERM,*) 'ERROR in SETLEV: ILEV or ISETLV is wrong'
        STOP
      ENDIF
C
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
C
      IF (ISETLV.EQ.2)  THEN
C
         NA  =KNA  (ILEV)
         NB  =KNB  (ILEV)
         NU  =KNU  (ILEV)
         NP  =KNP  (ILEV)
         NUP =KNUP (ILEV)
         NUPT=KNUPT(ILEV)
C
         KA1   =L(KLA    (ILEV))
         KST1  =L(KLST   (ILEV))
         KMASS1=L(KLMASS (ILEV))
         KM1   =L(KLM    (ILEV))
         KCOLA =L(KLCOLA (ILEV))
         KLDA  =L(KLLDA  (ILEV))
         KB1   =L(KLB1   (ILEV))
         KB2   =L(KLB2   (ILEV))
         KCOLB =L(KLCOLB (ILEV))
         KLDB  =L(KLLDB  (ILEV))
         KU1   =L(KLUP   (ILEV))
         KU2   =KU1+NU
         KP    =KU2+NU
         KT    =KP+NP
         KF1   =L(KLF12P (ILEV))
         KF2   =KF1+NU
         KFP   =KF2+NU
         KFT   =KFP+NP
         KAUX1 =L(KLAUX  (ILEV))
         KAUX2 =KAUX1+NU
         KAUXP =KAUX2+NU
         KAUXT =KAUXP+NP
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
      SUBROUTINE C2L2D (DPC,DPL,AREA,VAUX,KVERT,NEL,NVT,IPAR)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL AREA,VAUX
C
      PARAMETER (NNVE=4)
      DIMENSION DPC(*),DPL(*),VAUX(*),AREA(*),KVERT(NNVE,*)
C
C-----------------------------------------------------------------------
C
C *** constant to linear = 0
      IF (IPAR.EQ.0) THEN
C
       DO 10 IEL=1,NEL
C
       DPH  =DPC(IEL)
       DAREA=DBLE(AREA(IEL))
C
       IV1=KVERT(1,IEL)
       IV2=KVERT(2,IEL)
       IV3=KVERT(3,IEL)
       IV4=KVERT(4,IEL)
C
       DPL(IV1)=DPL(IV1)+DAREA*DPH
       DPL(IV2)=DPL(IV2)+DAREA*DPH
       DPL(IV3)=DPL(IV3)+DAREA*DPH
       DPL(IV4)=DPL(IV4)+DAREA*DPH
C
       VAUX(IV1)=VAUX(IV1)+REAL(DAREA)
       VAUX(IV2)=VAUX(IV2)+REAL(DAREA)
       VAUX(IV3)=VAUX(IV3)+REAL(DAREA)
       VAUX(IV4)=VAUX(IV4)+REAL(DAREA)
C
10     CONTINUE
C
       DO 20 IVT=1,NVT
20     DPL(IVT)=DPL(IVT)/DBLE(VAUX(IVT))
C
      ELSE
C
       DO 110 IEL=1,NEL
       DPC(IEL)=0.25D0*( DPL(KVERT(1,IEL))+DPL(KVERT(2,IEL))
     *                  +DPL(KVERT(3,IEL))+DPL(KVERT(4,IEL)))
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
      SUBROUTINE C2N2D (DPC,DPL,AREA,VAUX,KMID,NEL,NMT,NVT,IPAR)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      REAL AREA,VAUX
C
      PARAMETER (NNVE=4)
      DIMENSION DPC(*),DPL(*),VAUX(*),AREA(*),KMID(NNVE,*)
C
C-----------------------------------------------------------------------
C
C *** constant to nonconforming = 0
      IF (IPAR.EQ.0) THEN
C
       DO 10 IEL=1,NEL
C
       DPH  =DPC(IEL)
       DAREA=DBLE(AREA(IEL))
C
       IV1=KMID(1,IEL)-NVT
       IV2=KMID(2,IEL)-NVT
       IV3=KMID(3,IEL)-NVT
       IV4=KMID(4,IEL)-NVT
C
       DPL(IV1)=DPL(IV1)+DAREA*DPH
       DPL(IV2)=DPL(IV2)+DAREA*DPH
       DPL(IV3)=DPL(IV3)+DAREA*DPH
       DPL(IV4)=DPL(IV4)+DAREA*DPH
C
       VAUX(IV1)=VAUX(IV1)+REAL(DAREA)
       VAUX(IV2)=VAUX(IV2)+REAL(DAREA)
       VAUX(IV3)=VAUX(IV3)+REAL(DAREA)
       VAUX(IV4)=VAUX(IV4)+REAL(DAREA)
C
10     CONTINUE
C
       DO 20 IVT=1,NMT
20     DPL(IVT)=DPL(IVT)/DBLE(VAUX(IVT))
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
      SUBROUTINE IF217 (DA,DC,KCOL,KLD,DX,DB,DD,NEQ,NIT,OMEGA)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),DC(*),KCOL(*),KLD(*),DX(*),DB(*),DD(*)
C
      DO 10 ITE=1,NIT
      CALL LCP1 (DB,DD,NEQ)
      CALL LAX17(DA,KCOL,KLD,NEQ,DX,DD,-1D0,1D0)
      CALL IF117(DC,KCOL,KLD,DD,NEQ)
      CALL LLC1 (DD,DX,NEQ,OMEGA,1D0)
10    CONTINUE
C
      END
C
C
C
************************************************************************
      SUBROUTINE IF227 (VA,VC,KCOL,KLD,DX,DB,DD,NEQ,NIT,OMEGA)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),VC(*),KCOL(*),KLD(*),DX(*),DB(*),DD(*)
C
      DO 10 ITE=1,NIT
      CALL LCP1 (DB,DD,NEQ)
      CALL LAX37(VA,KCOL,KLD,NEQ,DX,DD,-1D0,1D0)
      CALL IF137(VC,KCOL,KLD,DD,NEQ)
      CALL LLC1 (DD,DX,NEQ,OMEGA,1D0)
10    CONTINUE
C
      END
C
C
C
************************************************************************
      SUBROUTINE IF237 (DA,VC,KCOL,KLD,DX,DB,DD,NEQ,NIT,OMEGA)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),VC(*),KCOL(*),KLD(*),DX(*),DB(*),DD(*)
C
      DO 10 ITE=1,NIT
      CALL LCP1 (DB,DD,NEQ)
      CALL LAX17(DA,KCOL,KLD,NEQ,DX,DD,-1D0,1D0)
      CALL IF137(VC,KCOL,KLD,DD,NEQ)
      CALL LLC1 (DD,DX,NEQ,OMEGA,1D0)
10    CONTINUE
C
      END
C
C
C
************************************************************************
      SUBROUTINE IF027 (VA,VC,KCOL,KLD,DX,DB,DD,NEQ,NIT,ITE,EPS,OMEGA,
     *                  RHO)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),VC(*),KCOL(*),KLD(*),DX(*),DB(*),DD(*)
C
C
C
      CALL LCP1 (DB,DD,NEQ)
      CALL LAX37(VA,KCOL,KLD,NEQ,DX,DD,-1D0,1D0)
      CALL LL21 (DD,NEQ,FD)
      IF (ABS(FD).LT.1D-12) THEN
       ITE=1
       RHO=0D0
       RETURN
      ENDIF       
C
      DO 10 ITE=1,NIT
C
      CALL LCP1 (DB,DD,NEQ)
      CALL LAX37(VA,KCOL,KLD,NEQ,DX,DD,-1D0,1D0)
      CALL LL21 (DD,NEQ,RES)
      IF (ABS(RES).LT.EPS*FD) THEN
       RHO=(RES/FD)**(1D0/DBLE(ITE-1))
       RETURN       
      ENDIF
C
      CALL IF137(VC,KCOL,KLD,DD,NEQ)
C
      CALL LLC1 (DD,DX,NEQ,OMEGA,1D0)
C
10    CONTINUE
C
      RHO=(RES/FD)**(1D0/DBLE(NIT))
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE IF037 (DA,VC,KCOL,KLD,DX,DB,DD,NEQ,NIT,ITE,EPS,OMEGA,
     *                  RHO)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),VC(*),KCOL(*),KLD(*),DX(*),DB(*),DD(*)
C
C
C
      CALL LCP1 (DB,DD,NEQ)
      CALL LAX17(DA,KCOL,KLD,NEQ,DX,DD,-1D0,1D0)
      CALL LL21 (DD,NEQ,FD)
      IF (ABS(FD).LT.1D-12) THEN
       ITE=1
       RHO=0D0
       RETURN
      ENDIF       
C
      DO 10 ITE=1,NIT
C
      CALL LCP1 (DB,DD,NEQ)
      CALL LAX17(DA,KCOL,KLD,NEQ,DX,DD,-1D0,1D0)
      CALL LL21 (DD,NEQ,RES)
      IF (ABS(RES).LT.EPS*FD) THEN
       RHO=(RES/FD)**(1D0/DBLE(ITE-1))
       RETURN       
      ENDIF
C
      CALL IF137(VC,KCOL,KLD,DD,NEQ)
C
      CALL LLC1 (DD,DX,NEQ,OMEGA,1D0)
C
10    CONTINUE
C
      RHO=(RES/FD)**(1D0/DBLE(NIT))
C
C
      END
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
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DCORVG(2,*),KVERT(4,*),KMID(4,*),KADJ(4,*),KNPR(*)
C
C
C
      DO 100 IEL=1,NEL/4
C
      IEL1=IEL
      IEL2=KADJ(2,IEL1)
      IEL3=KADJ(2,IEL2)
      IEL4=KADJ(2,IEL3)
C
      IADJ1=KADJ(1,IEL1)
      IADJ2=KADJ(1,IEL2)
      IADJ3=KADJ(1,IEL3)
      IADJ4=KADJ(1,IEL4)
C
      IVT1=KVERT(2,IEL1)
      IVT2=KVERT(2,IEL2)
      IVT3=KVERT(2,IEL3)
      IVT4=KVERT(2,IEL4)
C
      NADJ0=0
      IF (IADJ1.EQ.0) NADJ0=NADJ0+1
      IF (IADJ2.EQ.0) NADJ0=NADJ0+1
      IF (IADJ3.EQ.0) NADJ0=NADJ0+1
      IF (IADJ4.EQ.0) NADJ0=NADJ0+1
C
C
      IF (NADJ0.EQ.0) GOTO 100
C
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
      IVTM=KVERT(3,IEL)
      PXM=DCORVG(1,IVTM)
      PYM=DCORVG(2,IVTM)
C
C
      IF (NADJ0.EQ.1) THEN
       IF ((IADJ1.EQ.0).OR.(IADJ3.EQ.0)) THEN
        PX=0.5D0*(PX1+PX3)
        PY=0.5D0*(PY1+PY3)
       ELSE
        PX=0.5D0*(PX2+PX4)
        PY=0.5D0*(PY2+PY4)
       ENDIF
       DCORVG(1,IVTM)=PX
       DCORVG(2,IVTM)=PY
      ENDIF
C
C      
      IF (NADJ0.GT.1) THEN
       PX=0.25D0*(PX1+PX2+PX3+PX4)
       PY=0.25D0*(PY1+PY2+PY3+PY4)
       DCORVG(1,IVTM)=PX
       DCORVG(2,IVTM)=PY
      ENDIF
C     
100   CONTINUE
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE DYAX7(DA,KCOL,KLD,NEQ,DY,DX1,DX2,DYAX1,DYAX2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),DX1(*),DX2(*),DY(*)
C
      DYAX1=0D0
      DYAX2=0D0
C
      DO 1 IROW=1,NEQ
      DYIROW=DY(IROW)
      IF (DYIROW.EQ.0D0) GOTO 1
      AUX1=0D0
      AUX2=0D0
C
      DO 2 ICOL=KLD(IROW),KLD(IROW+1)-1
      JCOL=KCOL(ICOL)
      DAICOL=DA(ICOL)
      IF (DAICOL.EQ.0D0) GOTO 2
      AUX1=AUX1+DAICOL*DX1(JCOL)
      AUX2=AUX2+DAICOL*DX2(JCOL)
2     CONTINUE
C
      DYAX1=DYAX1+AUX1*DYIROW
      DYAX2=DYAX2+AUX2*DYIROW
1     CONTINUE
C
      END
C
C
C
************************************************************************
      SUBROUTINE VYAX7(VA,KCOL,KLD,NEQ,DY,DX1,DX2,DYAX1,DYAX2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VA(*),KCOL(*),KLD(*),DX1(*),DX2(*),DY(*)
C
      DYAX1=0D0
      DYAX2=0D0
C
      DO 1 IROW=1,NEQ
      DYIROW=DY(IROW)
      IF (DYIROW.EQ.0D0) GOTO 1
      AUX1=0D0
      AUX2=0D0
C
      DO 2 ICOL=KLD(IROW),KLD(IROW+1)-1
      JCOL=KCOL(ICOL)
      DAICOL=DBLE(VA(ICOL))
      IF (DAICOL.EQ.0D0) GOTO 2
      AUX1=AUX1+DAICOL*DX1(JCOL)
      AUX2=AUX2+DAICOL*DX2(JCOL)
2     CONTINUE
C
      DYAX1=DYAX1+AUX1*DYIROW
      DYAX2=DYAX2+AUX2*DYIROW
1     CONTINUE
C
      END
C
C
C
************************************************************************
      SUBROUTINE CRITAD(TIMEIN,TIMENS,TIMEST,EPSADI,EPSADL,EPSAD,IADIN)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
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
C
C
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
* II010                                                                *
*                                                                      *
* Purpose  Solution of a linear system  A*X = B  using                 *
*          a preconditioned BI-CGSTAB method                           *
*          Double precision version                                    *
*                                                                      *
* Subroutines/functions called  LSP1 , LLC1 , LL21 , LLI1, LCL1        *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX       R*8    Starting vector                                      *
* DB       R*8    Right hand side                                      *
* NEQ      I*4    Number of equations                                  *
* ITE      I*4    Minimum number of iterations (used for IREL=1 only)  *
* NIT      I*4    Maximum number of iterations                         *
* EPS      R*8    Desired precision                                    *
*                 IREL=0: Stop if !!RES!! < EPS                        *
*                 IREL=1: Stop if !!RES!!/!!RES0!! < EPS               *
*                         and a minimum of ITE iterations is performed *
* DAX0     SUBR   EXTERNAL Subroutine DAX0(DX,DAX,NEQ,A1,A2)           *
*                 Results  DAX := A1 * A * DX + A2 * DAX               *
* DCG0C    SUBR   EXTERNAL Subroutine DCG0C(DG,NEQ)                    *
*                 Results  DG := C**(-1) * DG  for the precondioning   *
*                 matrix C                                             *
* BNOCON   L*4    .TRUE.   No Preconditioning                          *
* DR,DD    R*8    Workspace vectors of length NEQ                      *
* DD1,DG   R*8    For BNOCON , DG must be replaced by DR               *
* IREL     I*4    See above                                            *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Solution vector                                      *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
*                 +1  Precision EPS not achieved after NIT iterations  *
*                                                                      *
************************************************************************
C
      SUBROUTINE II010(DX,DB,NEQ,NIT,ITE,EPS,DAX0,DCG0C,BNOCON,
     *                 DR,DR0,DP,DPA,DSA,IREL,CAPPA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DX(*),DB(*),DR(*),DR0(*),DP(*),DPA(*),DSA(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/
C
      SUB='II010'
      IF (ICHECK.GE.997) CALL OTRC('II010 ','04/14/93')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      NIT0=MAX(ITE,0)
      BREL=IREL.EQ.1
C
      CALL LCP1(DB,DR,NEQ)
      IF (BNOCON) CALL DCG0C(DR,NEQ)
C
      IF (ICHECK.GT.0) THEN
       CALL LLI1(DR,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL1(DX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IE010 ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      RHO0  =1D0
      DALPHA=1D0
      OMEGA0=1D0
C
      CALL LCP1(DB,DR,NEQ)
      CALL DAX0(DX,DR,NEQ,-1D0,1D0)
      IF (BNOCON) CALL DCG0C(DR,NEQ)
      CALL LL21(DR,NEQ,RES)
      IF (RES.LE.EPS.AND..NOT.BREL) THEN
       ITE=0
       FR=RES
       GOTO 200
      ENDIF
C
      CALL LCP1(DR,DR0,NEQ)
C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      CALL LSP1(DR0,DR,NEQ,RHO1)
      DBETA=(RHO1*DALPHA)/(RHO0*OMEGA0)
      RHO0 =RHO1
C
      CALL LLC1(DR ,DP,NEQ,1D0,DBETA)
      CALL LLC1(DPA,DP,NEQ,-DBETA*OMEGA0,1D0)
C
      CALL DAX0(DP,DPA,NEQ,1D0,0D0)
      IF (BNOCON) CALL DCG0C(DPA,NEQ)
C
      CALL LSP1(DR0,DPA,NEQ,DALPHA)
      DALPHA=RHO1/DALPHA
C
      CALL LLC1(DPA,DR,NEQ,-DALPHA,1D0)
C
      CALL DAX0(DR,DSA,NEQ,1D0,0D0)
      IF (BNOCON) CALL DCG0C(DSA,NEQ)
C
      CALL LSP1(DSA,DR ,NEQ,OMEGA1)
      CALL LSP1(DSA,DSA,NEQ,OMEGA2)
      OMEGA0=OMEGA1/OMEGA2
C
      CALL LLC1(DP ,DX ,NEQ,DALPHA,1D0)
      CALL LLC1(DR ,DX ,NEQ,OMEGA0,1D0)
C
      CALL LLC1(DSA,DR,NEQ,-OMEGA0,1D0)
C
C
      CALL LL21(DR,NEQ,FR)
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'IE010 ')
      ENDIF
      IF (BREL) THEN
       IF (FR.LE.RES*EPS.AND.ITE.GE.NIT0) GOTO 200
      ELSE
       IF (FR.LE.EPS) GOTO 200
      ENDIF
C
100   CONTINUE
C
      WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
c      CALL OMSG(71,'II010 ')
c      CALL OMSG(72,'II010 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/DBLE(NIT))
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'II010 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'II010 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/DBLE(ITE))
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'II010 ')
C
99999 END
C
C
C
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
* IE010                                                                *
*                                                                      *
* Purpose  Solution of a linear system  A*X = B  using                 *
*          a preconditioned conjugate gradient method                  *
*          Double precision version                                    *
*                                                                      *
* Subroutines/functions called  LSP1 , LLC1 , LL21 , LLI1, LCL1        *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX       R*8    Starting vector                                      *
* DB       R*8    Right hand side                                      *
* NEQ      I*4    Number of equations                                  *
* ITE      I*4    Minimum number of iterations (used for IREL=1 only)  *
* NIT      I*4    Maximum number of iterations                         *
* EPS      R*8    Desired precision                                    *
*                 IREL=0: Stop if !!RES!! < EPS                        *
*                 IREL=1: Stop if !!RES!!/!!RES0!! < EPS               *
*                         and a minimum of ITE iterations is performed *
* DAX0     SUBR   EXTERNAL Subroutine DAX0(DX,DAX,NEQ,A1,A2)           *
*                 Results  DAX := A1 * A * DX + A2 * DAX               *
* DCG0C    SUBR   EXTERNAL Subroutine DCG0C(DG,NEQ)                    *
*                 Results  DG := C**(-1) * DG  for the precondioning   *
*                 matrix C                                             *
* BNOCON   L*4    .TRUE.   No Preconditioning                          *
* DR,DD    R*8    Workspace vectors of length NEQ                      *
* DD1,DG   R*8    For BNOCON , DG must be replaced by DR               *
* IREL     I*4    See above                                            *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Solution vector                                      *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
*                 +1  Precision EPS not achieved after NIT iterations  *
*                                                                      *
************************************************************************
C
      SUBROUTINE IE010(DX,DB,NEQ,NIT,ITE,EPS,DAX0,DCG0C,BNOCON,
     *                 DR,DD,DD1,DG,IREL,CAPPA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DX(*),DB(*),DR(*),DG(*),DD(*),DD1(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
      DATA FR/0D0/
C
      SUB='IE010'
      IF (ICHECK.GE.997) CALL OTRC('IE010 ','01/02/89')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      NIT0=MAX(ITE,0)
      BREL=IREL.EQ.1
C
      IF (ICHECK.GT.0) THEN
       CALL LLI1(DB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL1(DX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IE010 ')
        GOTO 99999
       ENDIF
      ENDIF
C
C *** Initialization
      CALL DAX0(DX,DR,NEQ,1D0,0D0)
      CALL LLC1(DB,DR,NEQ,-1D0,1D0)
      CALL LL21(DR,NEQ,RES)
      IF (RES.LE.EPS.AND..NOT.BREL) THEN
       ITE=0
       FR=RES
       GOTO 200
      ENDIF
C
      IF (BNOCON) THEN
       SIGMA0=RES*RES
      ELSE
       CALL LCP1(DR,DG,NEQ)
       CALL DCG0C(DG,NEQ)
       CALL LSP1(DR,DG,NEQ,SIGMA0)
      ENDIF
C
      CALL LLC1(DG,DD,NEQ,-1D0,0D0)
C
C *** Iterative correction
      DO 100 ITE=1,NIT
C
      CALL DAX0(DD,DD1,NEQ,1D0,0D0)
      CALL LSP1(DD,DD1,NEQ,ALPHA)
      ALPHA=SIGMA0/ALPHA
      CALL LLC1(DD,DX,NEQ,ALPHA,1D0)
      CALL LLC1(DD1,DR,NEQ,ALPHA,1D0)
C
      CALL LL21(DR,NEQ,FR)
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,FR
       CALL OMSG(73,'IE010 ')
      ENDIF
      IF (BREL) THEN
       IF (FR.LE.RES*EPS.AND.ITE.GE.NIT0) GOTO 200
      ELSE
       IF (FR.LE.EPS) GOTO 200
      ENDIF
C
      IF (BNOCON) THEN
       SIGMA1=FR*FR
      ELSE
       CALL LCP1(DR,DG,NEQ)
       CALL DCG0C(DG,NEQ)
       CALL LSP1(DR,DG,NEQ,SIGMA1)
      ENDIF
C
      GAMMA=SIGMA1/SIGMA0
      SIGMA0=SIGMA1
      CALL LLC1(DG,DD,NEQ,-1D0,GAMMA)
100   CONTINUE
C
      WRITE (CPARAM,'(I15,2D25.16)') NIT,FR,RES
      CALL OMSG(71,'IE010 ')
      CALL OMSG(72,'IE010 ')
C
      IF (RES.GE.1D-70) THEN
       CAPPA=(FR/RES)**(1D0/NIT)
      ELSE
       CAPPA=0D0
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE010 ')
C
      IER=1
      GOTO 99999
C
200   IER=0
      IF (RES.GE.1D-70) RES=FR/RES
      WRITE (CPARAM,'(I15,2D25.16)') ITE,FR,RES
      CALL OMSG(72,'IE010 ')
C
      IF (ITE.EQ.0) THEN
       CAPPA=0D0
      ELSE
       CAPPA=RES**(1D0/ITE)
      ENDIF
      WRITE(CPARAM,'(D25.16)') CAPPA
      CALL OMSG(76,'IE010 ')
C
99999 END
