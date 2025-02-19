************************************************************************
      DOUBLE PRECISION FUNCTION DVISPO(D)
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (PI=3.1415926535897931D0)
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4)
      DIMENSION VWORK(1),KWORK(1)
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM) 
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
C

      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
c
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
C
c
      INCLUDE 'jump.inc'
      COMMON /NSCOUN/ NNONL,NMG    
C   
      common /fluidp/ pla,ple
      SAVE 
C
       IF(ISA .EQ. 1)DP =DWORK((KP)+IEL-1)/TSTEPH 
       IF(ISA .EQ. 0)DP =DWORK((KP)+IEL-1)
c
      IF(IMODEL .EQ. 3) THEN
         IF (DP.GT. DEPSI1)DPPRIM=DBETA
         IF (DP.LT. DEPSI1)DPPRIM=0.0d0
      ENDIF   
C
      IF(IMODEL .EQ. 4) DPPPRIM=(-DBETA/2)*(PLE+exp(DBETA*DP))**(-1/2-1) 
C
      IF(IMODEL .EQ. 5) DPPRIM=ALPHA*DBETA*exp(DBETA*DP) 
C
      IF(IMODEL .EQ. 6) DPPRIM=ALPHA*DBETA*exp(DBETA*DP)   
C
      IF(IMODEL .EQ. 7) DPP=DBETA
c

        IF (PLA .EQ. 0D0) THEN
C
           IF (IMODEL .LE. 2)THEN
                   DVISPO=0.0D0
           ENDIF
           IF (IMODEL .EQ. 3)THEN
                   DVISPO=NY*DPPRIM
           ENDIF
C
           IF (IMODEL .EQ. 4)THEN
                   DVISPO=NY*DPPRIM
           ENDIF
C
           IF (IMODEL .GE. 5)THEN
                   DVISPO=NY*DPPRIM
           ENDIF
C
        ELSE 
C
        IF(IMODEL .LT. 3)DVISPO=0.0D0
C
        IF (IMODEL .EQ. 3)THEN
           DVISPO=NY*(PLE+D)**(-PLA/2)*DPPRIM !POWER LAW X PRESSURE
        ENDIF
C
        IF (IMODEL .EQ. 4)THEN
           DVISPO=NY*(PLE+DPPRIM+D)**(-PLA/2) 
        ENDIF
C
        IF (IMODEL .GE. 5)THEN
           DVISPO=NY*(PLE+D)**(-PLA/2)*DPPRIM 
        ENDIF
C
        ENDIF ! (PLA .EQ. 0D0)
       END
C
************************************************************************
      DOUBLE PRECISION FUNCTION DVISCO(D)
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (PI=3.1415926535897931D0)
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4)
      DIMENSION VWORK(1),KWORK(1)
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM) 
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC
c
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
c
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
C
      INCLUDE 'jump.inc'
c
      COMMON /NSCOUN/ NNONL,NMG    
C   
      common /fluidp/ pla,ple
      SAVE 
C
C       DP =DWORK((KP)+IEL-1)/TSTEPH  
       IF(ISA .EQ. 1)DP =DWORK((KP)+IEL-1)/TSTEPH 
       IF(ISA .EQ. 0)DP =DWORK((KP)+IEL-1)
c
      IF(IMODEL .EQ. 3) THEN
         IF (DP.GT. DEPSI1)DPP=DBETA*DP
         IF (DP.LT. DEPSI1)DPP=1.0D-1
      ENDIF 
C
      IF(IMODEL .EQ. 4)DPP=(PLE+exp(DBETA*DP))**(-1/2) 
C
      IF(IMODEL .EQ. 5) DPP=ALPHA*(1.0D-3+exp(DBETA*DP))
C
      IF(IMODEL .EQ. 6) DPP=ALPHA*exp(DBETA*DP)   
C
      IF(IMODEL .EQ. 7) DPP=1.0D-3+DBETA*DP
C
      IF(IMODEL .EQ. 9) THEN
         IF(DP .LT. 1.0D0)DPP=DBETA*exp(DP-1.0d0)
         if(DP .GE. 1.0D0)DPP=DBETA*DP
      ENDIF   
c
      IF(IMODEL .EQ. 11)THEN
         IF(DP .LE. 0.1D0)DPP=0.1D0
         IF(DP .GT. 0.1D0)DPP=DP
      ENDIF
c 
        IF (PLA .EQ. 0D0) THEN
           DVISCO=NY
           IF (IMODEL .EQ. 3)THEN
                   DVISCO=NY*DPP
           ENDIF
C
           IF (IMODEL .EQ. 4)THEN
                   DVISCO=NY*DPP
           ENDIF
C
           IF (IMODEL .GE. 5)THEN
                   DVISCO=NY*DPP
           ENDIF
C
        ELSE 
C
        IF(IMODEL .EQ. 0)THEN ! PLASTICITY
           IF (D .LE. 1.0d0)DVISCO=NY
           IF (D .GT. 1.0d0)DVISCO=NY*(PLE+D)**(-PLA/2)
        ENDIF   
C
        IF (IMODEL .EQ. 1)
     *  DVISCO=NY*(PLE+D)**(-PLA/2) !POWER LAW
C
        IF(IMODEL .EQ. 2)
     *  DVISCO=NY*(1.0D0+PLE*D)**(-PLA/2) !CARREAU LAW 
C
        IF (IMODEL .EQ. 3)THEN
           DVISCO=NY*(PLE+D)**(-PLA/2)*DPP !POWER LAW X PRESSURE
        ENDIF
C
        IF (IMODEL .EQ. 4)THEN
           DVISCO=NY*(PLE+DPP+D)**(-PLA/2) 
        ENDIF
C
        IF (IMODEL .GE. 5)THEN
           DVISCO=NY*(PLE+D)**(-PLA/2)*DPP 
        ENDIF
C
        ENDIF ! (PLA .EQ. 0D0)
       END
C
************************************************************************
      DOUBLE PRECISION FUNCTION DVISCOPRIM(D)
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (PI=3.1415926535897931D0)
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4)
      DIMENSION VWORK(1),KWORK(1)
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
C
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
c
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
C
      INCLUDE 'jump.inc'

      COMMON /NSCOUN/ NNONL,NMG     
C
      common /fluidp/ pla,ple
      SAVE 
C
C       DP =DWORK((KP)+IEL-1)/TSTEPH  
       IF(ISA .EQ. 1)DP =DWORK((KP)+IEL-1)/TSTEPH 
       IF(ISA .EQ. 0)DP =DWORK((KP)+IEL-1)
C
      IF(IMODEL .EQ. 3) THEN
         IF (DP.GT. DEPSI1)DPP=DBETA*DP
         IF (DP.LT. DEPSI1)DPP=1.0D-1
      ENDIF 
C
      IF(IMODEL .EQ. 4) DPP=(PLE+exp(DBETA*DP))**(-1/2) 
C
      IF(IMODEL .EQ. 5) DPP=ALPHA*(1.0D-3+exp(DBETA*DP))
C
      IF(IMODEL .EQ. 6) DPP=ALPHA*exp(DBETA*DP)
C
      IF(IMODEL .EQ. 7) DPP=1.0D-3+DBETA*DP
C
      IF(IMODEL .EQ. 9) THEN
         IF(DP .LT. 1.0D0)DPP=DBETA*exp(DP-1.0d0)
         if(DP .GE. 1.0D0)DPP=DBETA*DP
      ENDIF  
C
      IF(IMODEL .EQ. 11)THEN
         IF(DP .LE. 0.1D0)DPP=0.1D0
         IF(DP .GT. 0.1D0)DPP=DP
      ENDIF
C
      IF(PLA .EQ. 0D0)THEN
         DVISCOPRIM=0D0
         ELSE
C
         IF(IMODEL .EQ. 0)THEN ! PLASTICITY
           IF (D .LE. 1.0d0)DVISCOPRIM=0D0
           IF (D .GT. 1.0d0)
     *       DVISCOPRIM=NY*(-PLA/2)*(PLE+D)**(-PLA/2-1)
         ENDIF   
C
         IF (IMODEL .EQ. 1)THEN
               DVISCOPRIM=NY*(-PLA/2)*(PLE+D)**(-PLA/2-1)
         ENDIF
C
         IF (IMODEL .EQ. 2)
     *       DVISCOPRIM=NY*PLE*(-PLA/2)*
     *                         (1.0D0+PLE*D)**(-PLA/2-1)
C
         IF (IMODEL .EQ. 3)THEN ! RELATED TO POWER LAW X PRESSURE
           DVISCOPRIM=NY*(-PLA/2)*(PLE+D)**(-PLA/2-1)*DPP
         ENDIF   
C
         IF (IMODEL .EQ. 4)THEN 
           DVISCOPRIM=NY*(-PLA/2)*(PLE+DPP+D)**(-PLA/2-1)
         ENDIF 
C
         IF (IMODEL .GE. 5)THEN 
           DVISCOPRIM=NY*(-PLA/2)*(PLE+D)**(-PLA/2-1)*DPP
         ENDIF
c
      ENDIF ! (PLA .EQ. 0D0)  
        END
C
************************************************************************
      SUBROUTINE SUPWDG(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                  KCOLA,KLDA,KVERT,KMID,DCORVG,ELE,COEFFN,
     *                  IDEF,DCMASS,locny)
************************************************************************
*     Purpose: -  Adds the SUPG-part on matrix block A after
*                 it was initialized by the linear part
*              -  The input vector Ui is the old velocity field
*              -  The input vectors UjLi are the transport directions
*     PARAMETRIC VERSION
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4)
      PARAMETER (NNLEV=9)
C
      DOUBLE PRECISION locny
      dimension locny(*)
      DIMENSION A(*),KCOLA(*),KLDA(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS)
      DIMENSION DENTRY(NNBAS,NNBAS,4),DHILF(NNBAS,NNBAS)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
C
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
C
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
      INCLUDE 'bouss.inc'
      INCLUDE 'jump.inc'
      EXTERNAL DVISCO,DVISCOPRIM
      SAVE 
C
C
      ICUB=ICUBN
C
      IF (IPRECA.EQ.4) THEN
       DNY=NY
      ELSE
       DNY=0D0
      ENDIF
C
      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
       CT0=DCMASS/THSTEP
      ELSE
       CT0=0D0
      ENDIF
C
C
      DO 1 I = 1,NNDER
1     BDER(I)=.FALSE.
C
      DO 2 I=1,3
2     BDER(I)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      ICUBP=ICUB
      CALL ELE(0D0,0D0,-2)
C
      DUMAX=0D0
      IF (A2L.EQ.0D0) THEN
       DO 10 IEQ=1,NMT
       DU1=U1L1(IEQ)
       DU2=U1L2(IEQ)
       DUNORM=SQRT(DU1**2+DU2**2)
       DUMAX=MAX(DUMAX,DUNORM)
10     CONTINUE
      ELSE       
       DO 20 IEQ=1,NMT
       DU1=A1L*U1L1(IEQ)+A2L*U2L1(IEQ)
       DU2=A1L*U1L2(IEQ)+A2L*U2L2(IEQ)
       DUNORM=SQRT(DU1**2+DU2**2)
       DUMAX=MAX(DUMAX,DUNORM)
20     CONTINUE
      ENDIF       
C
      IF (DUMAX.LT.1D-8) DUMAX=1D-8
C
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine local DELTA for streamline-diffusion
      CALL DELTSD(U1L1,U1L2,U2L1,U2L2,A1L,A2L,IEL,DUMAX,DELTA)
      IF (ISTOK.EQ.1) DELTA=0D0
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
         DENTRY(JDOFE,JDOFE,1)=0D0
         DENTRY(JDOFE,JDOFE,2)=0D0
         DENTRY(JDOFE,JDOFE,3)=0D0
         DENTRY(JDOFE,JDOFE,4)=0D0
            dhilf(jdofe,jdofe)=0d0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
          DENTRY(JDOFE,IDOFE,1)=0D0
          DENTRY(JDOFE,IDOFE,2)=0D0
          DENTRY(JDOFE,IDOFE,3)=0D0
          DENTRY(JDOFE,IDOFE,4)=0D0
            dhilf(jdofe,idofe)=0d0
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices
      DO 120 IVE = 1, NVE
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
      avgny=0d0
C *** Loop over all cubature points
      DO 200 ICUBP = 1, NCUBP
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
C *** Jacobian of the bilinear mapping onto the reference element
      DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
      DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
      DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
      DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
      DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
      OM=DOMEGA(ICUBP)*DETJ
C
      CALL ELE(XI1,XI2,-3)
      IF (IER.LT.0) GOTO 99999
C
C *** Evaluation of velocity in cubature points
      IF (DCMASS.NE.0D0) THEN
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
C     and derivatives
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
c
c
C ***  norm of symmetric part of velocity gradient
       IF (ISPGRAD .EQ. 1)DUGSQ=DUX1**2+0.5d0*(DUX2+DUY1)**2+DUY2**2
c
C ***  norm of gradient
       IF (ISPGRAD .EQ. 0)DUGSQ=DUX1**2+DUX2**2+DUY1**2+DUY2**2
       DNY=DVISCO(DUGSQ)
       avgny=avgny+dny/NY
c
        IF (ISPGRAD .EQ. 0)THEN
C ***  use just grad
C ***  Summing up over all pairs of multiindices
       DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2
C
       DO 235 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
        AH1= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *     + DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
        AH4= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *     + DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
        HSUMI=HBASI2*DU1+HBASI3*DU2
             AH1= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +  CT0*HBASJ1*HBASI1
             AH4= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +  CT0*HBASJ1*HBASI1
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            DHILF (JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE,4)+OM*AH1           
235    CONTINUE
230    CONTINUE
C
       ELSE
C
C***   use symmetric part of grad  
C   
C ***  Summing up over all pairs of multiindices
       DO 240 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2
C
       DO 245 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *          + DNY*(2d0*HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH2= DNY*HBASJ2*HBASJ3
             AH3= DNY*HBASJ3*HBASJ2
             AH4= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *          + DNY*(HBASJ2**2+2d0*HBASJ3**2)+CT0*HBASJ1**2
             ahilf=DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
        HSUMI=HBASI2*DU1+HBASI3*DU2
             AH1= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(2d0*HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +  CT0*HBASJ1*HBASI1
             AH2= DNY*HBASI2*HBASJ3!HBASI3*HBASJ2
             AH3= DNY*HBASI3*HBASJ2!HBASI2*HBASJ3
             AH4= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(HBASJ2*HBASI2+2d0*HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             ahilf=DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +  CT0*HBASJ1*HBASI1
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,2)=DENTRY(JDOFE,IDOFE,2)+OM*AH2
            DENTRY(JDOFE,IDOFE,3)=DENTRY(JDOFE,IDOFE,3)+OM*AH3
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            DHILF (JDOFE,IDOFE)=DHILF (JDOFE,IDOFE)+OM*AHilf            
            
 245        CONTINUE
240    CONTINUE
       ENDIF ! ISPGRAD
      ELSE
C
C ***  Summing up over 1 pair of multiindices
       DO 250 JDOFE=1,IDFL
       HBASJ1=DBAS(KDFL(JDOFE),1)
C
       DO 260 IDOFE=1,IDFL
       HBASI1=DBAS(KDFL(IDOFE),1)
C
       AH=-1D0/THSTEP*HBASJ1*HBASI1
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH
260    CONTINUE
250    CONTINUE
C
      ENDIF
C
200   CONTINUE
C
      locny(iel)=avgny/ncubp
C
      DO 400 JDOFE=1,IDFL
      DO 400 IDOFE=1,IDFL
          DENTH1=THSTEP*DENTRY(JDOFE,IDOFE,1)
          DENTH2=THSTEP*DENTRY(JDOFE,IDOFE,2)
          DENTH3=THSTEP*DENTRY(JDOFE,IDOFE,3)
          DENTH4=THSTEP*DENTRY(JDOFE,IDOFE,4)
          DHHILF=THSTEP*DHILF (JDOFE,IDOFE)
C
      IF (IDEF.LT.2) THEN
       IA   =KENTRY(JDOFE,IDOFE)
         IF (IPRECO.EQ.0) THEN !Def-Tensor as Precond.
           A(IA)=A(IA)+REAL(DENTH1)
           A(NA+IA)=A(NA+IA)+REAL(DENTH2)
           A(2*NA+IA)=A(2*NA+IA)+REAL(DENTH3)
           A(3*NA+IA)=A(3*NA+IA)+REAL(DENTH4)
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
c
          ELSE !GRAD as precond.
          A(IA)=A(IA)+REAL(DHHILF)
           A(NA+IA)=0d0
           A(2*NA+IA)=0d0
           A(3*NA+IA)=A(3*NA+IA)+REAL(DHHILF)
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
           ENDIf
      ENDIF
C
      IF (IDEF.GT.0) THEN 
       IDFG=KDFG(IDOFE)
       JDFG=KDFG(JDOFE)
           D1(JDFG)= D1(JDFG)-DENTH1*U1(IDFG)-DENTH2*U2(IDFG)
           D2(JDFG)= D2(JDFG)-DENTH3*U1(IDFG)-DENTH4*U2(IDFG)
      ENDIF 
C
400   CONTINUE
C
100   CONTINUE
C
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE SUPWNP(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                  KCOLA,KLDA,KVERT,KMID,DCORVG,ELE,COEFFN,
     *                  IDEF,DCMASS,locny)
************************************************************************
*     Purpose: -  Adds the SUPG-part on matrix block A after
*                 it was initialized by the linear part
*              -  The input vector Ui is the old velocity field
*              -  The input vectors UjLi are the transport directions
*     NONPARAMETRIC VERSION
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNCOF=6)
      PARAMETER (NNLEV=9)
C
      DOUBLE PRECISION locny
      dimension locny(*)
      DIMENSION A(*),KCOLA(*),KLDA(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS)
      DIMENSION DENTRY(NNBAS,NNBAS,4),DHILF(NNBAS,NNBAS)
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
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
C
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
ccc   COMMON /UPTIME/ TTU0,TTU11,TTU12,TTU13,TTU2,TTU3,TTU4
C
      INCLUDE 'jump.inc'
      INCLUDE 'bouss.inc'
      EXTERNAL DVISCO
      SAVE 
C
c      print *,'A1sup=',(A(i),i=1,25)

C
ccc   CALL ZTIME(TTT0)
      ICUB=ICUBN
C
c$$$      IF (IPRECA.EQ.4) THEN
c$$$       DNY=NY
c$$$      ELSE
c$$$       DNY=0D0
c$$$      ENDIF
C
      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
       CT0=DCMASS/THSTEP
      ELSE
       CT0=0D0
      ENDIF
C
C
      DO 1 I = 1,NNDER
1     BDER(I)=.FALSE.
C
      DO 2 I=1,3
2     BDER(I)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      DUMAX=0D0
      IF (A2L.EQ.0D0) THEN
       DO 10 IEQ=1,NMT
       DU1=U1L1(IEQ)
       DU2=U1L2(IEQ)
       DUNORM=SQRT(DU1**2+DU2**2)
       DUMAX=MAX(DUMAX,DUNORM)
10     CONTINUE
      ELSE       
       DO 20 IEQ=1,NMT
       DU1=A1L*U1L1(IEQ)+A2L*U2L1(IEQ)
       DU2=A1L*U1L2(IEQ)+A2L*U2L2(IEQ)
       DUNORM=SQRT(DU1**2+DU2**2)
       DUMAX=MAX(DUMAX,DUNORM)
20     CONTINUE
      ENDIF       
C
      IF (DUMAX.LT.1D-8) DUMAX=1D-8
C
C
ccc   CALL ZTIME(TTT1)
ccc   TTU0=TTU0+TTT1-TTT0
C *** Loop over all elements
      DO 100 IEL=1,NEL
ccc   CALL ZTIME(TTT0)
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine local DELTA for streamline-diffusion
      CALL DELTSD(U1L1,U1L2,U2L1,U2L2,A1L,A2L,IEL,DUMAX,DELTA)
      IF (ISTOK.EQ.1) DELTA=0D0
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
         DENTRY(JDOFE,JDOFE,1)=0D0
         DENTRY(JDOFE,JDOFE,2)=0D0
         DENTRY(JDOFE,JDOFE,3)=0D0
         DENTRY(JDOFE,JDOFE,4)=0D0
            dhilf(jdofe,jdofe)=0d0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
          DENTRY(JDOFE,IDOFE,1)=0D0
          DENTRY(JDOFE,IDOFE,2)=0D0
          DENTRY(JDOFE,IDOFE,3)=0D0
          DENTRY(JDOFE,IDOFE,4)=0D0
            dhilf(jdofe,idofe)=0d0
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices
      DO 120 IVE = 1, NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
120   CONTINUE
ccc   CALL ZTIME(TTT1)
ccc   TTU11=TTU11+TTT1-TTT0
ccc   CALL ZTIME(TTT0)
C
      DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
      DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
      DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
      DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
C
ccc   CALL ZTIME(TTT1)
ccc   TTU2=TTU2+TTT1-TTT0
C *** Dummy call - ELE may save arithmetic operations
ccc   CALL ZTIME(TTT0)
      CALL ELE(0D0,0D0,-2)
      IF (IER.LT.0) GOTO 99999
ccc   CALL ZTIME(TTT1)
ccc   TTU3=TTU3+TTT1-TTT0
C
      avgny=0d0
C *** Loop over all cubature points
      DO 200 ICUBP = 1, NCUBP
ccc   CALL ZTIME(TTT0)
C
C *** Cubature points on the reference element
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
C
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
ccc   CALL ZTIME(TTT1)
ccc   TTU2=TTU2+TTT1-TTT0
ccc   CALL ZTIME(TTT0)
      CALL ELE(XX,YY,-3)
      IF (IER.LT.0) GOTO 99999
C
ccc   CALL ZTIME(TTT1)
ccc   TTU4=TTU4+TTT1-TTT0
ccc   CALL ZTIME(TTT0)
C *** Evaluation of velocity in cubature points
      IF (DCMASS.NE.0D0) THEN
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
C     and derivatives
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
c
C ***  norm of symmetric part of velocity gradient
       IF (ISPGRAD .EQ. 1)DUGSQ=DUX1**2+0.5d0*(DUX2+DUY1)**2+DUY2**2
c
C ***  norm of gradient
       IF (ISPGRAD .EQ. 0)DUGSQ=DUX1**2+DUX2**2+DUY1**2+DUY2**2
C
          DNY=DVISCO(DUGSQ)
          avgny=avgny+dny/NY
C   
        IF (ISPGRAD .EQ. 0)THEN
C ***  use just grad
C ***  Summing up over all pairs of multiindices
       DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
       HSUMJ=(HBASJ2*DU1+HBASJ3*DU2)!*0d0!!!WEG zu testzwecken
C
       DO 235 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *          + DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH4= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *          + DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
        HSUMI=(HBASI2*DU1+HBASI3*DU2)!*0d0!!!WEG zu testzwecken
             AH1= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH4= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            DHILF (JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE,4)+OM*AH1           
c      IF (KENTRY(JDOFE,IDOFE).eq.5) print *,'su[', DENTRY(JDOFE,IDOFE,1)
            
 235        CONTINUE
230    CONTINUE
c
ccc   CALL ZTIME(TTT1)
ccc   TTU12=TTU12+TTT1-TTT0
C
       ELSE
C
C***   use symmetric part of grad
c       DNY=2D0*DNY
C    
C ***  Summing up over all pairs of multiindices
       DO 240 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2
C
       DO 245 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *          + DNY*(2d0*HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH2= DNY*HBASJ2*HBASJ3
             AH3= DNY*HBASJ3*HBASJ2
             AH4= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *          + DNY*(HBASJ2**2+2d0*HBASJ3**2)+CT0*HBASJ1**2
             ahilf=DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
        HSUMI=HBASI2*DU1+HBASI3*DU2
             AH1= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(2d0*HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH2= DNY*HBASI2*HBASJ3!HBASI3*HBASJ2
             AH3= DNY*HBASI3*HBASJ2!HBASI2*HBASJ3
             AH4= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(HBASJ2*HBASI2+2d0*HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             ahilf=DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +  CT0*HBASJ1*HBASI1
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,2)=DENTRY(JDOFE,IDOFE,2)+OM*AH2
            DENTRY(JDOFE,IDOFE,3)=DENTRY(JDOFE,IDOFE,3)+OM*AH3
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            DHILF (JDOFE,IDOFE)=DHILF (JDOFE,IDOFE)+OM*AHilf            
            
 245        CONTINUE
240    CONTINUE
       ENDIF ! ISPGRAD
C
      ELSE
C
C ***  Summing up over 1 pair of multiindices
       DO 250 JDOFE=1,IDFL
       HBASJ1=DBAS(KDFL(JDOFE),1)
C
       DO 260 IDOFE=1,IDFL
       HBASI1=DBAS(KDFL(IDOFE),1)
C
       AH=-1D0/THSTEP*HBASJ1*HBASI1
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH
 
260    CONTINUE
250    CONTINUE
ccc   CALL ZTIME(TTT1)
ccc   TTU12=TTU12+TTT1-TTT0
C
      ENDIF
C
200   CONTINUE
C
      locny(iel)=avgny/ncubp
C
ccc   CALL ZTIME(TTT0)
      DO 400 JDOFE=1,IDFL
      DO 400 IDOFE=1,IDFL
          DENTH1=THSTEP*DENTRY(JDOFE,IDOFE,1)
          DENTH2=THSTEP*DENTRY(JDOFE,IDOFE,2)
          DENTH3=THSTEP*DENTRY(JDOFE,IDOFE,3)
          DENTH4=THSTEP*DENTRY(JDOFE,IDOFE,4)
          DHHILF=THSTEP*DHILF (JDOFE,IDOFE)
C
          BSONST=.false.
      IF (IDEF.LT.2) THEN
       IA   =KENTRY(JDOFE,IDOFE)
         IF ((IPRECO.EQ.0).or.(BSONST)) THEN !Def-Tensor as Precond.
c            IF (ia.eq.1) write (*,*) ilev,IPRECO, MOD(IZBV2,IPRECO)
            A(IA)=A(IA)+REAL(DENTH1)
           A(NA+IA)  =A(NA+IA)+REAL(DENTH2)
           A(2*NA+IA)=A(2*NA+IA)+REAL(DENTH3)
           A(3*NA+IA)=A(3*NA+IA)+REAL(DENTH4)
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
c
          ELSEIF (IPRECO.EQ.1) THEN !GRAD as precond.
          A(IA)=A(IA)+REAL(DHHILF)
           A(NA+IA)=0d0
           A(2*NA+IA)=0d0
           A(3*NA+IA)=A(3*NA+IA)+REAL(DHHILF)
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
          ELSE
           A(IA)     =A(IA)     +0.5d0*(REAL(DHHILF)+REAL(DENTH1))
           A(NA+IA)  =A(  NA+IA)+0.5d0*(0d0+REAL(DENTH2))
           A(2*NA+IA)=A(2*NA+IA)+0.5d0*(0d0+REAL(DENTH3))
           A(3*NA+IA)=A(3*NA+IA)+0.5d0*(REAL(DHHILF)+REAL(DENTH4))
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
           ENDIf
      ENDIF

c$$$
c$$$         IF (IPRECO.EQ.0) THEN !Def-Tensor as Precond.
c$$$           A(IA)=A(IA)+REAL(DENTH1)
c$$$           A(NA+IA)=A(NA+IA)+REAL(DENTH2)
c$$$           A(2*NA+IA)=A(2*NA+IA)+REAL(DENTH3)
c$$$           A(3*NA+IA)=A(3*NA+IA)+REAL(DENTH4)
c$$$           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
c$$$c
c$$$          ELSE !GRAD as precond.
c$$$          A(IA)=A(IA)+REAL(DHHILF)
c$$$           A(NA+IA)=0d0
c$$$           A(2*NA+IA)=0d0
c$$$           A(3*NA+IA)=A(3*NA+IA)+REAL(DHHILF)
c$$$           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
c$$$           ENDIf
c$$$      ENDIF
C
      IF (IDEF.GT.0) THEN 
       IDFG=KDFG(IDOFE)
       JDFG=KDFG(JDOFE)
           D1(JDFG)= D1(JDFG)-DENTH1*U1(IDFG)-DENTH2*U2(IDFG)
           D2(JDFG)= D2(JDFG)-DENTH3*U1(IDFG)-DENTH4*U2(IDFG)
      ENDIF 
C
C
400   CONTINUE
C
C
ccc   CALL ZTIME(TTT1)
ccc   TTU13=TTU13+TTT1-TTT0
100   CONTINUE
c
C
ccc   write(6,*) TTU0,TTU11,TTU12,TTU13,TTU2,TTU3,TTU4,NCUBP
C
c      print *,ilev'Asup=',(A(i),i=1,25)

99999 END
C
C
C
************************************************************************
      SUBROUTINE  DELTSD  (U1L1,U1L2,U2L1,U2L2,A1L,A2L,IEL,DUMAX,DELTA)
************************************************************************
*     Calculates coefficient Delta for SD
*     RELOC = local RE
*     UNORM=local velocity norm
*     HLOCAL=local mesh width
************************************************************************
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNBAS=21,NNLEV=9,NNARR=299,NNWORK=1)
C
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
C
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
      SAVE
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C
C
C
      AREA=DBLE(VWORK(L(KLAREA(ILEV))+IEL-1))
      HLOCAL=SQRT(AREA)
C
      IF (UPSAM.LT.0D0) THEN
C
       DELTA=ABS(UPSAM)*HLOCAL
C
      ELSE
C
       DU1=0D0
       DU2=0D0
       DO 10 IDOF=1,IDFL
       DU1=DU1+(A1L*U1L1(KDFG(IDOF))+A2L*U2L1(KDFG(IDOF)))
10     DU2=DU2+(A1L*U1L2(KDFG(IDOF))+A2L*U2L2(KDFG(IDOF)))
C
       UNORM=0.25D0*SQRT(DU1**2+DU2**2)
C
C       IF (UNORM.LE.1D-8) THEN
C        DELTA=0D0
C       ELSE
        RELOC=UNORM*HLOCAL/NY
        DELTA=UPSAM*HLOCAL/DUMAX*2D0*(RELOC/(1D0+RELOC))
C       ENDIF
C
      ENDIF
C
C      WRITE(6,*) DELTA,DELTA/(UPSAM*HLOCAL),RELOC,UNORM,HLOCAL
C
C
      END
************************************************************************
      SUBROUTINE ADDSTP(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *  KCOLA,KLDA,KVERT,KMID,DCORVG,ELE,COEFFN,
     *  IDEF,DCMASS,locny)
************************************************************************
*     Purpose: -  Adds the SUPG-part on matrix block A after
*                 it was initialized by the linear part
*              -  The input vector Ui is the old velocity field
*              -  The input vectors UjLi are the transport directions
*     PARAMETRIC VERSION
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4)
      PARAMETER (NNLEV=9)
C
      DOUBLE PRECISION locny
      dimension locny(*)
      DIMENSION A(*),KCOLA(*),KLDA(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS)
      DIMENSION DENTRY(NNBAS,NNBAS,4),DHILF(NNBAS,NNBAS)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
C
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
C
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
      INCLUDE 'bouss.inc'
      EXTERNAL DVISCO
      SAVE 
C
C
      ICUB=ICUBN
C
       DNY=NY
C
      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
       CT0=DCMASS/THSTEP
      ELSE
       CT0=0D0
      ENDIF
C
C
      DO 1 I = 1,NNDER
1     BDER(I)=.FALSE.
C
      DO 2 I=1,3
2     BDER(I)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      ICUBP=ICUB
      CALL ELE(0D0,0D0,-2)
C
C
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
         DENTRY(JDOFE,JDOFE,1)=0D0
         DENTRY(JDOFE,JDOFE,2)=0D0
         DENTRY(JDOFE,JDOFE,3)=0D0
         DENTRY(JDOFE,JDOFE,4)=0D0
            dhilf(jdofe,jdofe)=0d0
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
          DENTRY(JDOFE,IDOFE,1)=0D0
          DENTRY(JDOFE,IDOFE,2)=0D0
          DENTRY(JDOFE,IDOFE,3)=0D0
          DENTRY(JDOFE,IDOFE,4)=0D0
            dhilf(jdofe,idofe)=0d0
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices
      DO 120 IVE = 1, NVE
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
      avgny=0d0
C *** Loop over all cubature points
      DO 200 ICUBP = 1, NCUBP
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
C *** Jacobian of the bilinear mapping onto the reference element
      DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
      DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
      DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
      DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
      DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
      OM=DOMEGA(ICUBP)*DETJ
C
      CALL ELE(XI1,XI2,-3)
      IF (IER.LT.0) GOTO 99999
C
C *** Evaluation of velocity in cubature points
      IF (DCMASS.NE.0D0) THEN
C
C     and derivatives
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
c
C ***  norm of symmetric part of velocity gradient
       IF (ISPGRAD .EQ. 1)DUGSQ=DUX1**2+0.5d0*(DUX2+DUY1)**2+DUY2**2
c
C ***  norm of gradient
       IF (ISPGRAD .EQ. 0)DUGSQ=DUX1**2+DUX2**2+DUY1**2+DUY2**2
C
          DNY=DVISCO(DUGSQ)
          avgny=avgny+dny/NY

C
        IF (ISPGRAD .EQ. 0)THEN
C ***  use just grad    
C ***  Summing up over all pairs of multiindices
       DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
C
       DO 235 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
        HSUMI=HBASI2*DU1+HBASI3*DU2
             AH1= DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH1
            DHILF (JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE,4)+OM*AH1           
            
 235        CONTINUE
230    CONTINUE
C     
       ELSE
C
C***  use symmetric part of grad
c
C
C ***  Summing up over all pairs of multiindices
       DO 240 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2
C
       DO 245 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= DNY*(2d0*HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH2= DNY*HBASJ2*HBASJ3
             AH3= DNY*HBASJ3*HBASJ2
             AH4= DNY*(HBASJ2**2+2d0*HBASJ3**2)+CT0*HBASJ1**2
             ahilf=DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
        HSUMI=HBASI2*DU1+HBASI3*DU2
             AH1= DNY*(2d0*HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH2= DNY*HBASI2*HBASJ3!HBASI3*HBASJ2
             AH3= DNY*HBASI3*HBASJ2!HBASI2*HBASJ3
             AH4= DNY*(HBASJ2*HBASI2+2d0*HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             ahilf=DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +  CT0*HBASJ1*HBASI1
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,2)=DENTRY(JDOFE,IDOFE,2)+OM*AH2
            DENTRY(JDOFE,IDOFE,3)=DENTRY(JDOFE,IDOFE,3)+OM*AH3
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            DHILF (JDOFE,IDOFE)=DHILF (JDOFE,IDOFE)+OM*AHilf            
            
 245        CONTINUE
240    CONTINUE
       ENDIF ! ISPGRAD       
      ELSE
C
C ***  Summing up over 1 pair of multiindices
       DO 250 JDOFE=1,IDFL
       HBASJ1=DBAS(KDFL(JDOFE),1)
C
       DO 260 IDOFE=1,IDFL
       HBASI1=DBAS(KDFL(IDOFE),1)
C
       AH=-1D0/THSTEP*HBASJ1*HBASI1
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH
260    CONTINUE
250    CONTINUE
C
      ENDIF
C
200   CONTINUE
C
      locny(iel)=avgny/ncubp
C
      DO 400 JDOFE=1,IDFL
      DO 400 IDOFE=1,IDFL
          DENTH1=THSTEP*DENTRY(JDOFE,IDOFE,1)
          DENTH2=THSTEP*DENTRY(JDOFE,IDOFE,2)
          DENTH3=THSTEP*DENTRY(JDOFE,IDOFE,3)
          DENTH4=THSTEP*DENTRY(JDOFE,IDOFE,4)
          DHHILF=THSTEP*DHILF (JDOFE,IDOFE)
C
      IF (IDEF.LT.2) THEN
       IA   =KENTRY(JDOFE,IDOFE)
         IF (IPRECO.EQ.0) THEN !Def-Tensor as Precond.
           A(IA)=A(IA)+REAL(DENTH1)
           A(NA+IA)=A(NA+IA)+REAL(DENTH2)
           A(2*NA+IA)=A(2*NA+IA)+REAL(DENTH3)
           A(3*NA+IA)=A(3*NA+IA)+REAL(DENTH4)
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
c
          ELSE !GRAD as precond.
          A(IA)=A(IA)+REAL(DHHILF)
           A(NA+IA)=0d0
           A(2*NA+IA)=0d0
           A(3*NA+IA)=A(3*NA+IA)+REAL(DHHILF)
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
           ENDIf
      ENDIF
C
      IF (IDEF.GT.0) THEN 
       IDFG=KDFG(IDOFE)
       JDFG=KDFG(JDOFE)
           D1(JDFG)= D1(JDFG)-DENTH1*U1(IDFG)-DENTH2*U2(IDFG)
           D2(JDFG)= D2(JDFG)-DENTH3*U1(IDFG)-DENTH4*U2(IDFG)
      ENDIF 
C
400   CONTINUE
C
100   CONTINUE
C
C     
99999  END
C
C
C     
************************************************************************
      SUBROUTINE ADDSTN(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *  KCOLA,KLDA,KVERT,KMID,DCORVG,ELE,COEFFN,
     *  IDEF,DCMASS,locny)
************************************************************************
*     Purpose: -  Adds the M/ST-part on matrix block A 
*     -  The input vector Ui is the old velocity field
*     NONPARAMETRIC VERSION
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C     
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNCOF=6)
      PARAMETER (NNLEV=9)
C     
      DOUBLE PRECISION locny
      dimension locny(*)
      DIMENSION A(*),KCOLA(*),KLDA(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS)
      DIMENSION DENTRY(NNBAS,NNBAS,4),DHILF(NNBAS,NNBAS)
C     
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *  DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      COMMON /COFBAS/ COB(NNBAS,NNCOF)
C     
C***  COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
C     
C***  user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *  IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *  INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *  NSM,NSL,NSMFAC
C     
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *  EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *  RLXSM,RLXSL,AMINMG,AMAXMG
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
C
      INCLUDE 'jump.inc'
      INCLUDE 'bouss.inc'
      EXTERNAL DVISCO
      SAVE 
C    
      CALL VKADD(AVK1,AVK2)
C
ccc   CALL ZTIME(TTT0)
      ICUB=ICUBN
C     
       DNY=NY
C 
c

c
c      print *,ilev,'A1=',(A(i),i=1,25)
      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
       CT0=DCMASS/THSTEP
      ELSE
       CT0=0D0
      ENDIF
C     
C     
      DO 1 I = 1,NNDER
       BDER(I)=.FALSE.
 1    continue
C     
      DO 2 I=1,3
       BDER(I)=.TRUE.
 2    continue
C     
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C     
************************************************************************
C***  Calculation of the matrix - storage technique 7 or 8
************************************************************************
C     
C     
ccc   CALL ZTIME(TTT1)
ccc   TTU0=TTU0+TTT1-TTT0
C***  Loop over all elements
       DO 100 IEL=1,NEL
ccc   CALL ZTIME(TTT0)
        CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999
C     
C     
C***  Determine entry positions in matrix
        DO 110 JDOFE=1,IDFL
         ILD=KLDA(KDFG(JDOFE))
         KENTRY(JDOFE,JDOFE)=ILD
         DENTRY(JDOFE,JDOFE,1)=0D0
         DENTRY(JDOFE,JDOFE,2)=0D0
         DENTRY(JDOFE,JDOFE,3)=0D0
         DENTRY(JDOFE,JDOFE,4)=0D0
            dhilf(jdofe,jdofe)=0d0
         JCOL0=ILD
         DO 111 IDOFE=1,IDFL
          IF (IDOFE.EQ.JDOFE) GOTO 111
          IDFG=KDFG(IDOFE)
          DO 112 JCOL=JCOL0,NA
           IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
 112      CONTINUE
 113      JCOL0=JCOL+1
          KENTRY(JDOFE,IDOFE)=JCOL
          DENTRY(JDOFE,IDOFE,1)=0D0
          DENTRY(JDOFE,IDOFE,2)=0D0
          DENTRY(JDOFE,IDOFE,3)=0D0
          DENTRY(JDOFE,IDOFE,4)=0D0
            dhilf(jdofe,idofe)=0d0
 111     CONTINUE
 110    CONTINUE
C     
C***  Evaluation of coordinates of the vertices
        BRANDE=.FALSE.
        DO 120 IVE = 1, NVE
         JP=KVERT(IVE,IEL)
         KVE(IVE)=JP
         DX(IVE)=DCORVG(1,JP)
         DY(IVE)=DCORVG(2,JP)
c =========================================
c     Hier mache ich den Randcheck!
c =========================================
c        IF ((DX(IVE)*DY(IVE).eq.1d0).or.(DX(IVE)*DY(IVE).eq.0d0))
        IF ((DX(IVE)).ge.0.5d0)!.or.(DX(IVE)*DY(IVE).eq.0d0))
     *  BRANDE=.TRUE.! DANN MACHE ICH, BEI IGRAD/IPRECO=2 DEF, Sonst GRAD
c        BRANDE=.TRUE.
c
c =========================================
 120    CONTINUE
ccc   CALL ZTIME(TTT1)
ccc   TTU11=TTU11+TTT1-TTT0
ccc   CALL ZTIME(TTT0)
C     

        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
C     
ccc   CALL ZTIME(TTT1)
ccc   TTU2=TTU2+TTT1-TTT0
C***  Dummy call - ELE may save arithmetic operations
ccc   CALL ZTIME(TTT0)
        CALL ELE(0D0,0D0,-2)
        IF (IER.LT.0) GOTO 99999
ccc   CALL ZTIME(TTT1)
ccc   TTU3=TTU3+TTT1-TTT0
C     
      avgny=0d0
C***  Loop over all cubature points
        DO 200 ICUBP = 1, NCUBP
ccc   CALL ZTIME(TTT0)
C     
C***  Cubature points on the reference element
         XI1=DXI(ICUBP,1)
         XI2=DXI(ICUBP,2)
C     
C***  Jacobian of the bilinear mapping onto the reference element
         DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
         DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
         DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
         DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
         DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
C     
C***  Cubature points on actual Element + weights
         XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *     +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
         YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1
     *     +0.5D0*(DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2
         OM=DOMEGA(ICUBP)*DETJ
C  
ccc   CALL ZTIME(TTT1)
ccc   TTU2=TTU2+TTT1-TTT0
ccc   CALL ZTIME(TTT0)
         CALL ELE(XX,YY,-3)
         IF (IER.LT.0) GOTO 99999
C     
ccc   CALL ZTIME(TTT1)
ccc   TTU4=TTU4+TTT1-TTT0
ccc   CALL ZTIME(TTT0)
C***  Evaluation of velocity in cubature points
         IF (DCMASS.NE.0D0) THEN
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
c
C     and derivatives
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
          DNY=DVISCO(DUGSQ)
          avgny=avgny+dny/NY
C 
        IF (ISPGRAD .EQ. 0)THEN
c     use grad
C***  Summing up over all pairs of multiindices
          DO 230 JDOFE=1,IDFL
           JDOFEH=KDFL(JDOFE)
           HBASJ1=DBAS(JDOFEH,1)
           HBASJ2=DBAS(JDOFEH,2)
           HBASJ3=DBAS(JDOFEH,3)
C     
C     
           DO 235 IDOFE=1,IDFL
            IF (IDOFE.EQ.JDOFE) THEN
             AH1= DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
            ELSE
             IDOFEH=KDFL(IDOFE)
             HBASI1=DBAS(IDOFEH,1)
             HBASI2=DBAS(IDOFEH,2)
             HBASI3=DBAS(IDOFEH,3)
c     
c     J -test function   I -trial function
c     
             AH1= DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +  CT0*HBASJ1*HBASI1
            ENDIF
C     
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH1
            DHILF (JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE,4)+OM*AH1           
c       IF (KENTRY(JDOFE,IDOFE).eq.5) print *,iel,'up',
c     * DENTRY(JDOFE,IDOFE,1)
           
 235     CONTINUE
 230      CONTINUE
C     
       ELSE
C
C***  use symmetric part of grad
c          
C
C ***  Summing up over all pairs of multiindices
       DO 240 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2
C
       DO 245 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= DNY*(2d0*HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH2= DNY*HBASJ2*HBASJ3
             AH3= DNY*HBASJ3*HBASJ2
             AH4= DNY*(HBASJ2**2+2d0*HBASJ3**2)+CT0*HBASJ1**2
             ahilf=DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
        HSUMI=HBASI2*DU1+HBASI3*DU2
             AH1= DNY*(2d0*HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH2= DNY*HBASI2*HBASJ3!HBASI3*HBASJ2
             AH3= DNY*HBASI3*HBASJ2!HBASI2*HBASJ3
             AH4= DNY*(HBASJ2*HBASI2+2d0*HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             ahilf=DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +  CT0*HBASJ1*HBASI1
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,2)=DENTRY(JDOFE,IDOFE,2)+OM*AH2
            DENTRY(JDOFE,IDOFE,3)=DENTRY(JDOFE,IDOFE,3)+OM*AH3
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            DHILF (JDOFE,IDOFE)=DHILF (JDOFE,IDOFE)+OM*AHilf            
            
 245        CONTINUE
240    CONTINUE
       ENDIF ! ISPGRAD
         ELSE
c     dcmass.eq.0d0
C     
C***  Summing up over 1 pair of multiindices
          DO 250 JDOFE=1,IDFL
           HBASJ1=DBAS(KDFL(JDOFE),1)
C     
           DO 260 IDOFE=1,IDFL
            HBASI1=DBAS(KDFL(IDOFE),1)
C     
            AH=-1D0/THSTEP*HBASJ1*HBASI1
C     
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH
 260       CONTINUE
 250      CONTINUE

C     
         ENDIF
C     
 200    CONTINUE
C
      locny(iel)=avgny/ncubp
C     
ccc   CALL ZTIME(TTT0)
        DO 400 JDOFE=1,IDFL
         DO 400 IDOFE=1,IDFL
          DENTH1=THSTEP*DENTRY(JDOFE,IDOFE,1)
          DENTH2=THSTEP*DENTRY(JDOFE,IDOFE,2)
          DENTH3=THSTEP*DENTRY(JDOFE,IDOFE,3)
          DENTH4=THSTEP*DENTRY(JDOFE,IDOFE,4)
          DHHILF=THSTEP*DHILF (JDOFE,IDOFE)
C     
          BSONST=.false.
          goto 9911
c          IF ((IPRECO.EQ.0).or.(IPRECO.eq.1)) GOTO 9911
c          IF (ILEV.eq.NLMAX) THEN
c             IF ((MOD(IZBV2+1,IPRECO).eq.0)) BSONST=.true.
c          ELSE
c             IF ((MOD(IZBV2,IPRECO).eq.0)) BSONST=.true.
c          ENDIF
c          IF (IZBV2.le.5) BSONST=.true.
 9911     CONTINUE
      IF (IDEF.LT.2) THEN
       IA   =KENTRY(JDOFE,IDOFE)
         IF ((IPRECO.EQ.0).or.(BSONST)) THEN !Def-Tensor as Precond.
c            IF (ia.eq.1) write (*,*) ilev,IPRECO, MOD(IZBV2,IPRECO)
            A(IA)=A(IA)+REAL(DENTH1)
           A(NA+IA)  =A(NA+IA)+REAL(DENTH2)
           A(2*NA+IA)=A(2*NA+IA)+REAL(DENTH3)
           A(3*NA+IA)=A(3*NA+IA)+REAL(DENTH4)
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
c
          ELSEIF (IPRECO.EQ.1) THEN !GRAD as precond.
           A(IA)=A(IA)+REAL(DHHILF)
           A(NA+IA)=0d0
           A(2*NA+IA)=0d0
           A(3*NA+IA)=A(3*NA+IA)+REAL(DHHILF)
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
          ELSE
c        IF (brande) print *,'up final', denth2,denth1,denth4
          A(IA)=A(IA)+REAL(DENTH1)
           A(NA+IA)  =A(NA+IA)+REAL(DENTH2)
           A(2*NA+IA)=A(2*NA+IA)+REAL(DENTH3)
           A(3*NA+IA)=A(3*NA+IA)+REAL(DENTH4)
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
           ENDIf
      ENDIF
C     
          IF (IDEF.GT.0) THEN 
           IDFG=KDFG(IDOFE)
           JDFG=KDFG(JDOFE)
           D1(JDFG)= D1(JDFG)-DENTH1*U1(IDFG)-DENTH2*U2(IDFG)
           D2(JDFG)= D2(JDFG)-DENTH3*U1(IDFG)-DENTH4*U2(IDFG)
c           write (*,*) 'in addstn', denth1
          ENDIF 
C     
 400     CONTINUE
C     
ccc   CALL ZTIME(TTT1)
ccc   TTU13=TTU13+TTT1-TTT0
 100    CONTINUE
C     
ccc   write(6,*) TTU0,TTU11,TTU12,TTU13,TTU2,TTU3,TTU4,NCUBP
C     
cc       if (IGRAD.eq.2) THEN
cc          IGRAD=ig1
cc          IPRECO=ig2
cc       ENDIF
c  
c      print *,ilev,'A=',(A(i),i=1,25)
C     
99999  END
C
C************************************************************************
      SUBROUTINE ADDSN1(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *  KCOLA,KLDA,KVERT,KMID,DCORVG,ELE,COEFFN,
     *  IDEF,DCMASS,locny)
************************************************************************
*     Purpose: -  Adds the M/ST-part on matrix block A 
*     -  The input vector Ui is the old velocity field
*     NONPARAMETRIC VERSION. Hier: GRAD-GRAD TEIL
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C     
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNCOF=6)
      PARAMETER (NNLEV=9)
C     
      DOUBLE PRECISION locny
      dimension locny(*)
      DIMENSION A(*),KCOLA(*),KLDA(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS)
      DIMENSION DENTRY(NNBAS,NNBAS,4),DHILF(NNBAS,NNBAS)
C     
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *  DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      COMMON /COFBAS/ COB(NNBAS,NNCOF)
C     
C***  COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
C     
C***  user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *  IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *  INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *  NSM,NSL,NSMFAC
C     
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *  EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *  RLXSM,RLXSL,AMINMG,AMAXMG
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
C
      INCLUDE 'jump.inc'
      INCLUDE 'bouss.inc'
      EXTERNAL DVISCO
      SAVE 
C    
      CALL VKADD(AVK1,AVK2)
C
      ICUB=ICUBN
C     
       DNY=NY
C     
      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
       CT0=DCMASS/THSTEP
      ELSE
       CT0=0D0
      ENDIF
C  

C     
      DO 1 I = 1,NNDER
       BDER(I)=.FALSE.
 1    continue
C     
      DO 2 I=1,3
       BDER(I)=.TRUE.
 2    continue
C     
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C     
************************************************************************
C***  Calculation of the matrix - storage technique 7 or 8
************************************************************************
C     
C     
       DO 100 IEL=1,NEL
c
        CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999
C     
C     
C***  Determine entry positions in matrix
        DO 110 JDOFE=1,IDFL
         ILD=KLDA(KDFG(JDOFE))
         KENTRY(JDOFE,JDOFE)=ILD
         DENTRY(JDOFE,JDOFE,1)=0D0
         DENTRY(JDOFE,JDOFE,2)=0D0
         DENTRY(JDOFE,JDOFE,3)=0D0
         DENTRY(JDOFE,JDOFE,4)=0D0
            dhilf(jdofe,jdofe)=0d0
         JCOL0=ILD
         DO 111 IDOFE=1,IDFL
          IF (IDOFE.EQ.JDOFE) GOTO 111
          IDFG=KDFG(IDOFE)
          DO 112 JCOL=JCOL0,NA
           IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
 112      CONTINUE
 113      JCOL0=JCOL+1
          KENTRY(JDOFE,IDOFE)=JCOL
          DENTRY(JDOFE,IDOFE,1)=0D0
          DENTRY(JDOFE,IDOFE,2)=0D0
          DENTRY(JDOFE,IDOFE,3)=0D0
          DENTRY(JDOFE,IDOFE,4)=0D0
            dhilf(jdofe,idofe)=0d0
 111     CONTINUE
 110    CONTINUE
C     
C***  Evaluation of coordinates of the vertices
        DO 120 IVE = 1, NVE
         JP=KVERT(IVE,IEL)
         KVE(IVE)=JP
         DX(IVE)=DCORVG(1,JP)
         DY(IVE)=DCORVG(2,JP)
 120    CONTINUE

C     
        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
C     
C***  Dummy call - ELE may save arithmetic operations
c
        CALL ELE(0D0,0D0,-2)
        IF (IER.LT.0) GOTO 99999
C     
      avgny=0d0
C***  Loop over all cubature points
        DO 200 ICUBP = 1, NCUBP
c
C     
C***  Cubature points on the reference element
         XI1=DXI(ICUBP,1)
         XI2=DXI(ICUBP,2)
C     
C***  Jacobian of the bilinear mapping onto the reference element
         DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
         DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
         DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
         DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
         DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
C     
C***  Cubature points on actual Element + weights
         XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *     +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
         YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1
     *     +0.5D0*(DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2
         OM=DOMEGA(ICUBP)*DETJ
C     
         CALL ELE(XX,YY,-3)
         IF (IER.LT.0) GOTO 99999
C     
C***  Evaluation of velocity in cubature points
         IF (DCMASS.NE.0D0) THEN
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
c
C     and derivatives
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
          DNY=DVISCO(DUGSQ)
          avgny=avgny+dny/NY
C       
          IF (ISPGRAD .EQ. 0)THEN

C *** use grad
C***  Summing up over all pairs of multiindices

          DO 230 JDOFE=1,IDFL
           JDOFEH=KDFL(JDOFE)
           HBASJ1=DBAS(JDOFEH,1)
           HBASJ2=DBAS(JDOFEH,2)
           HBASJ3=DBAS(JDOFEH,3)
C     
C     
           DO 235 IDOFE=1,IDFL
            IF (IDOFE.EQ.JDOFE) THEN
             AH1= DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
            ELSE
             IDOFEH=KDFL(IDOFE)
             HBASI1=DBAS(IDOFEH,1)
             HBASI2=DBAS(IDOFEH,2)
             HBASI3=DBAS(IDOFEH,3)
c     
c     J -test function   I -trial function
c     
             AH1= DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +  CT0*HBASJ1*HBASI1
            ENDIF
C     
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH1
            DHILF (JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE,4)+OM*AH1           
c           
 235     CONTINUE
 230      CONTINUE
C
       ELSE
C
C***   use symmetric part of grad
c        
C
C ***  Summing up over all pairs of multiindices
       DO 240 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2
C
       DO 245 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= DNY*(2d0*HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH2= DNY*HBASJ2*HBASJ3
             AH3= DNY*HBASJ3*HBASJ2
             AH4= DNY*(HBASJ2**2+2d0*HBASJ3**2)+CT0*HBASJ1**2
             ahilf=DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
        HSUMI=HBASI2*DU1+HBASI3*DU2
             AH1= DNY*(2d0*HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH2= DNY*HBASI2*HBASJ3!HBASI3*HBASJ2
             AH3= DNY*HBASI3*HBASJ2!HBASI2*HBASJ3
             AH4= DNY*(HBASJ2*HBASI2+2d0*HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             ahilf=DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +  CT0*HBASJ1*HBASI1
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,2)=DENTRY(JDOFE,IDOFE,2)+OM*AH2
            DENTRY(JDOFE,IDOFE,3)=DENTRY(JDOFE,IDOFE,3)+OM*AH3
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            DHILF (JDOFE,IDOFE)=DHILF (JDOFE,IDOFE)+OM*AHilf            
            
 245        CONTINUE
240    CONTINUE
       ENDIF ! ISPGRAD
c
         ELSE!     dcmass.eq.0d0
C     
C***  Summing up over 1 pair of multiindices
          DO 250 JDOFE=1,IDFL
           HBASJ1=DBAS(KDFL(JDOFE),1)
C     
           DO 260 IDOFE=1,IDFL
            HBASI1=DBAS(KDFL(IDOFE),1)
C     
            AH=-1D0/THSTEP*HBASJ1*HBASI1
C     
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH
 260       CONTINUE
 250      CONTINUE
C     
         ENDIF
C     
 200    CONTINUE
C
      locny(iel)=avgny/ncubp
C     
c
        DO 400 JDOFE=1,IDFL
         DO 400 IDOFE=1,IDFL
          DENTH1=THSTEP*DENTRY(JDOFE,IDOFE,1)
          DENTH2=THSTEP*DENTRY(JDOFE,IDOFE,2)
          DENTH3=THSTEP*DENTRY(JDOFE,IDOFE,3)
          DENTH4=THSTEP*DENTRY(JDOFE,IDOFE,4)
          DHHILF=THSTEP*DHILF (JDOFE,IDOFE)
C     
          BSONST=.false.
          IF ((IPRECO.EQ.0).or.(IPRECO.eq.1)) GOTO 9911
          IF (ILEV.eq.NLMAX) THEN
             IF ((MOD(IZBV2+1,IPRECO).eq.0)) BSONST=.true.
          ELSE
             IF ((MOD(IZBV2,IPRECO).eq.0)) BSONST=.true.
          ENDIF
 9911     CONTINUE
c
      IF (IDEF.LT.2) THEN
       IA   =KENTRY(JDOFE,IDOFE)
       A(IA)=A(IA)+REAL(DHHILF)
       A(3*NA+IA)=A(3*NA+IA)+REAL(DHHILF)
      ENDIF! Auf die Matrix kommt nur der Grad Teil, Rest in ADDSN2
C     
          IF (IDEF.GT.0) THEN 
           IDFG=KDFG(IDOFE)
           JDFG=KDFG(JDOFE)
           D1(JDFG)= D1(JDFG)-DENTH1*U1(IDFG)-DENTH2*U2(IDFG)
           D2(JDFG)= D2(JDFG)-DENTH3*U1(IDFG)-DENTH4*U2(IDFG)
          ENDIF !Den Defekt berechne ich aber voll hier!
C     
 400     CONTINUE
c      
C     
 100    CONTINUE
C     
C     
c 
C     
99999  END
*##########################            #######################################
C
************************************************************************
      SUBROUTINE SUPWDGNEW(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                  KCOLA,KLDA,KVERT,KMID,DCORVG,ELE,COEFFN,
     *                  IDEF,DCMASS,locny)
************************************************************************
*     Purpose: -  Adds the SUPG-part on matrix block A after
*                 it was initialized by the linear part
*              -  The input vector Ui is the old velocity field
*              -  The input vectors UjLi are the transport directions
*     PARAMETRIC VERSION
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4)
      PARAMETER (NNLEV=9)
C
      DOUBLE PRECISION locny
      dimension locny(*)
      DIMENSION A(*),KCOLA(*),KLDA(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS)
      DIMENSION DENTRY(NNBAS,NNBAS,4),DHILF(NNBAS,NNBAS)
      DIMENSION DEJACOBNS(NNBAS,NNBAS,4),DEJACOBPW(NNBAS,NNBAS,4)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
C
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
C
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
      INCLUDE 'bouss.inc'
      EXTERNAL DVISCO,DVISCOPRIM
      SAVE 
C
      ICUB=ICUBNEWTON
C
      IF (IPRECA.EQ.4) THEN
       DNY=NY
      ELSE
       DNY=0D0
      ENDIF
C
      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
       CT0=DCMASS/THSTEP
      ELSE
       CT0=0D0
      ENDIF
C
C
      DO 1 I = 1,NNDER
1     BDER(I)=.FALSE.
C
      DO 2 I=1,3
2     BDER(I)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      ICUBP=ICUB
      CALL ELE(0D0,0D0,-2)
C
      DUMAX=0D0
      IF (A2L.EQ.0D0) THEN
       DO 10 IEQ=1,NMT
       DU1=U1L1(IEQ)
       DU2=U1L2(IEQ)
       DUNORM=SQRT(DU1**2+DU2**2)
       DUMAX=MAX(DUMAX,DUNORM)
10     CONTINUE
      ELSE       
       DO 20 IEQ=1,NMT
       DU1=A1L*U1L1(IEQ)+A2L*U2L1(IEQ)
       DU2=A1L*U1L2(IEQ)+A2L*U2L2(IEQ)
       DUNORM=SQRT(DU1**2+DU2**2)
       DUMAX=MAX(DUMAX,DUNORM)
20     CONTINUE
      ENDIF       
C
      IF (DUMAX.LT.1D-8) DUMAX=1D-8
C
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine local DELTA for streamline-diffusion
      CALL DELTSD(U1L1,U1L2,U2L1,U2L2,A1L,A2L,IEL,DUMAX,DELTA)
c      IF (ISTOK.EQ.1 .OR. DNNS .NE. 0D0) DELTA=0D0 
      IF (ISTOK.EQ.1 ) DELTA=0D0
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
         DENTRY(JDOFE,JDOFE,1)=0D0
         DENTRY(JDOFE,JDOFE,2)=0D0
         DENTRY(JDOFE,JDOFE,3)=0D0
         DENTRY(JDOFE,JDOFE,4)=0D0
            dhilf(jdofe,jdofe)=0d0
C
         DEJACOBNS(JDOFE,JDOFE,1)=0D0
         DEJACOBNS(JDOFE,JDOFE,2)=0D0
         DEJACOBNS(JDOFE,JDOFE,3)=0D0
         DEJACOBNS(JDOFE,JDOFE,4)=0D0
C 
         DEJACOBPW(JDOFE,JDOFE,1)=0D0
         DEJACOBPW(JDOFE,JDOFE,2)=0D0
         DEJACOBPW(JDOFE,JDOFE,3)=0D0
         DEJACOBPW(JDOFE,JDOFE,4)=0D0
C
C
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
          DENTRY(JDOFE,IDOFE,1)=0D0
          DENTRY(JDOFE,IDOFE,2)=0D0
          DENTRY(JDOFE,IDOFE,3)=0D0
          DENTRY(JDOFE,IDOFE,4)=0D0
            dhilf(jdofe,idofe)=0d0
C
         DEJACOBNS(JDOFE,IDOFE,1)=0D0
         DEJACOBNS(JDOFE,IDOFE,2)=0D0
         DEJACOBNS(JDOFE,IDOFE,3)=0D0
         DEJACOBNS(JDOFE,IDOFE,4)=0D0
C 
         DEJACOBPW(JDOFE,IDOFE,1)=0D0
         DEJACOBPW(JDOFE,IDOFE,2)=0D0
         DEJACOBPW(JDOFE,IDOFE,3)=0D0
         DEJACOBPW(JDOFE,IDOFE,4)=0D0
C
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices
      DO 120 IVE = 1, NVE
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
      avgny=0d0
C *** Loop over all cubature points
      DO 200 ICUBP = 1, NCUBP
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
C *** Jacobian of the bilinear mapping onto the reference element
      DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
      DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
      DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
      DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
      DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
      OM=DOMEGA(ICUBP)*DETJ
C
      CALL ELE(XI1,XI2,-3)
      IF (IER.LT.0) GOTO 99999
C
C *** Evaluation of velocity in cubature points
      IF (DCMASS.NE.0D0) THEN
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
C     and derivatives
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
c
C
C ***  norm of symmetric part of velocity gradient
       IF (ISPGRAD .EQ. 1)DUGSQ=DUX1**2+0.5d0*(DUX2+DUY1)**2+DUY2**2
c
C ***  norm of gradient
       IF (ISPGRAD .EQ. 0)DUGSQ=DUX1**2+DUX2**2+DUY1**2+DUY2**2
c
          DNY    =DVISCO(DUGSQ)
          DNYPRIM=2.0D0*DVISCOPRIM(DUGSQ)
          avgny=avgny+dny/NY
C
          DUDOTGRADU1= DU1*DUX1+DU2*DUY1
          DUDOTGRADU2= DU1*DUX2+DU2*DUY2  
C  
C
        IF (ISPGRAD .EQ. 0)THEN
C
C ***  use just grad
C ***  Summing up over all pairs of multiindices
       DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2
C
       HYUMJ=DUX1*HBASJ2+DUY1*HBASJ3
       HZUMJ=DUX2*HBASJ2+DUY2*HBASJ3
C
       DO 235 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
        AH1= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *     +  DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
        AH4= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *     +  DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
c
             AJ1=  DUX1*HBASJ1**2
             AJ2=  DUY1*HBASJ1**2
             AJ3=  DUX2*HBASJ1**2
             AJ4=  DUY2*HBASJ1**2
C
             AK1=  DNYPRIM*HYUMJ**2
             AK2=  DNYPRIM*HZUMJ*HYUMJ
             AK3=  DNYPRIM*HYUMJ*HZUMJ
             AK4=  DNYPRIM*HZUMJ**2
C
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
        HSUMI=HBASI2*DU1+HBASI3*DU2
C
       HYUMI=DUX1*HBASI2+DUY1*HBASI3
       HZUMI=DUX2*HBASI2+DUY2*HBASI3
C
             AH1= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +   DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH4= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +   DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
C
             AJ1=  DUX1*HBASJ1*HBASI1
             AJ2=  DUY1*HBASJ1*HBASI1
             AJ3=  DUX2*HBASJ1*HBASI1
             AJ4=  DUY2*HBASJ1*HBASI1
C
             AK1=  DNYPRIM*HYUMI*HYUMJ
             AK2=  DNYPRIM*HZUMI*HYUMJ
             AK3=  DNYPRIM*HYUMI*HZUMJ
             AK4=  DNYPRIM*HZUMI*HZUMJ
C
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            DHILF (JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE,4)+OM*AH1    
C
         DEJACOBNS(JDOFE,IDOFE,1)=DEJACOBNS(JDOFE,IDOFE,1)+OM*AJ1
         DEJACOBNS(JDOFE,IDOFE,2)=DEJACOBNS(JDOFE,IDOFE,2)+OM*AJ2
         DEJACOBNS(JDOFE,IDOFE,3)=DEJACOBNS(JDOFE,IDOFE,3)+OM*AJ3
         DEJACOBNS(JDOFE,IDOFE,4)=DEJACOBNS(JDOFE,IDOFE,4)+OM*AJ4
C
         DEJACOBPW(JDOFE,IDOFE,1)=DEJACOBPW(JDOFE,IDOFE,1)+OM*AK1
         DEJACOBPW(JDOFE,IDOFE,2)=DEJACOBPW(JDOFE,IDOFE,2)+OM*AK2
         DEJACOBPW(JDOFE,IDOFE,3)=DEJACOBPW(JDOFE,IDOFE,3)+OM*AK3
         DEJACOBPW(JDOFE,IDOFE,4)=DEJACOBPW(JDOFE,IDOFE,4)+OM*AK4 
C     
235    CONTINUE
230    CONTINUE
c
       ELSE
C
C***  use symmetric part of grad
c      
       DNYPRIM=2D0*DNYPRIM
C
C ***  Summing up over all pairs of multiindices
       DO 240 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2
C
       HYUMJ=DUX1*HBASJ2+0.5D0*(DUY1+DUX2)*HBASJ3
       HZUMJ=0.5D0*(DUY1+DUX2)*HBASJ2+DUY2*HBASJ3
c
       DO 245 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(2d0*HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH2= DNY*HBASJ2*HBASJ3
             AH3= DNY*HBASJ3*HBASJ2
             AH4= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(HBASJ2**2+2d0*HBASJ3**2)+CT0*HBASJ1**2
             ahilf=DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
c
             AJ1=  DUX1*HBASJ1**2
             AJ2=  DUY1*HBASJ1**2
             AJ3=  DUX2*HBASJ1**2
             AJ4=  DUY2*HBASJ1**2
C
             AK1=  DNYPRIM*HYUMJ**2
             AK2=  DNYPRIM*HZUMJ*HYUMJ
             AK3=  DNYPRIM*HYUMJ*HZUMJ
             AK4=  DNYPRIM*HZUMJ**2
c
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
        HSUMI=HBASI2*DU1+HBASI3*DU2
C
       HYUMI=DUX1*HBASI2+0.5D0*(DUY1+DUX2)*HBASI3
       HZUMI=0.5D0*(DUY1+DUX2)*HBASI2+DUY2*HBASI3
C
             AH1= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +   DNY*(2d0*HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH2= DNY*HBASI2*HBASJ3!HBASI3*HBASJ2
             AH3= DNY*HBASI3*HBASJ2!HBASI2*HBASJ3
             AH4= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +   DNY*(HBASJ2*HBASI2+2d0*HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             ahilf=DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +  CT0*HBASJ1*HBASI1
c
             AJ1=  DUX1*HBASJ1*HBASI1
             AJ2=  DUY1*HBASJ1*HBASI1
             AJ3=  DUX2*HBASJ1*HBASI1
             AJ4=  DUY2*HBASJ1*HBASI1
C
             AK1=  DNYPRIM*HYUMI*HYUMJ
             AK2=  DNYPRIM*HZUMI*HYUMJ
             AK3=  DNYPRIM*HYUMI*HZUMJ
             AK4=  DNYPRIM*HZUMI*HZUMJ
c
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,2)=DENTRY(JDOFE,IDOFE,2)+OM*AH2
            DENTRY(JDOFE,IDOFE,3)=DENTRY(JDOFE,IDOFE,3)+OM*AH3
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            DHILF (JDOFE,IDOFE)=DHILF (JDOFE,IDOFE)+OM*AHilf 
c
         DEJACOBNS(JDOFE,IDOFE,1)=DEJACOBNS(JDOFE,IDOFE,1)+OM*AJ1
         DEJACOBNS(JDOFE,IDOFE,2)=DEJACOBNS(JDOFE,IDOFE,2)+OM*AJ2
         DEJACOBNS(JDOFE,IDOFE,3)=DEJACOBNS(JDOFE,IDOFE,3)+OM*AJ3
         DEJACOBNS(JDOFE,IDOFE,4)=DEJACOBNS(JDOFE,IDOFE,4)+OM*AJ4
C
         DEJACOBPW(JDOFE,IDOFE,1)=DEJACOBPW(JDOFE,IDOFE,1)+OM*AK1
         DEJACOBPW(JDOFE,IDOFE,2)=DEJACOBPW(JDOFE,IDOFE,2)+OM*AK2
         DEJACOBPW(JDOFE,IDOFE,3)=DEJACOBPW(JDOFE,IDOFE,3)+OM*AK3
         DEJACOBPW(JDOFE,IDOFE,4)=DEJACOBPW(JDOFE,IDOFE,4)+OM*AK4 
c           
            
 245        CONTINUE
240    CONTINUE
       ENDIF ! ISPGRAD
C
      ELSE
C
C ***  Summing up over 1 pair of multiindices
       DO 250 JDOFE=1,IDFL
       HBASJ1=DBAS(KDFL(JDOFE),1)
C
       DO 260 IDOFE=1,IDFL
       HBASI1=DBAS(KDFL(IDOFE),1)
C
       AH=-1D0/THSTEP*HBASJ1*HBASI1
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH
260    CONTINUE
250    CONTINUE
C
      ENDIF
C
200   CONTINUE
c

      DCONST=0D0
      locny(iel)=avgny/ncubp
C
      DO 400 JDOFE=1,IDFL
      DO 400 IDOFE=1,IDFL

          DENTH1=THSTEP*DENTRY(JDOFE,IDOFE,1)*DCONST
          DENTH2=THSTEP*DENTRY(JDOFE,IDOFE,2)*DCONST
          DENTH3=THSTEP*DENTRY(JDOFE,IDOFE,3)*DCONST
          DENTH4=THSTEP*DENTRY(JDOFE,IDOFE,4)*DCONST
          DHHILF=THSTEP*DHILF (JDOFE,IDOFE)*DCONST
C
          DENTJ1=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,1)
          DENTJ2=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,2)
          DENTJ3=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,3)
          DENTJ4=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,4)
C
          DENTK1=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,1)
          DENTK2=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,2)
          DENTK3=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,3)
          DENTK4=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,4)
C
C
      IF (IDEF.LT.2) THEN
       IA   =KENTRY(JDOFE,IDOFE)
         IF (IPRECO.EQ.0) THEN !Def-Tensor as Precond.
           A(IA)=A(IA)+REAL(DENTH1)+REAL(DENTJ1)+REAL(DENTK1)
           A(NA+IA)=A(NA+IA)+REAL(DENTH2)+REAL(DENTJ2)+REAL(DENTK2)
           A(2*NA+IA)=A(2*NA+IA)+REAL(DENTH3)+REAL(DENTJ3)+REAL(DENTK3)
           A(3*NA+IA)=A(3*NA+IA)+REAL(DENTH4)+REAL(DENTJ4)+REAL(DENTK4)
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
c
          ELSE !GRAD as precond.
          A(IA)=A(IA)+REAL(DHHILF)
           A(NA+IA)=0d0
           A(2*NA+IA)=0d0
           A(3*NA+IA)=A(3*NA+IA)+REAL(DHHILF)
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
           ENDIf
      ENDIF

      IF (IDEF.GT.0) THEN 
       IDFG=KDFG(IDOFE)
       JDFG=KDFG(JDOFE)
           D1(JDFG)= D1(JDFG)-DENTH1*U1(IDFG)-DENTH2*U2(IDFG)
           D2(JDFG)= D2(JDFG)-DENTH3*U1(IDFG)-DENTH4*U2(IDFG)
      ENDIF 
C
400   CONTINUE
C
100   CONTINUE
C
99999 END
C
C
************************************************************************
      SUBROUTINE SUPWNPNEW(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                  KCOLA,KLDA,KVERT,KMID,DCORVG,ELE,COEFFN,
     *                  IDEF,DCMASS,locny)
************************************************************************
*     Purpose: -  Adds the SUPG-part on matrix block A after
*                 it was initialized by the linear part
*              -  The input vector Ui is the old velocity field
*              -  The input vectors UjLi are the transport directions
*     NONPARAMETRIC VERSION
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNCOF=6)
      PARAMETER (NNLEV=9)
C
      DOUBLE PRECISION locny
      dimension locny(*)
      DIMENSION A(*),KCOLA(*),KLDA(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS)
      DIMENSION DENTRY(NNBAS,NNBAS,4),DHILF(NNBAS,NNBAS)
      DIMENSION DEJACOBNS(NNBAS,NNBAS,4),DEJACOBPW(NNBAS,NNBAS,4)
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
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
C
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
ccc   COMMON /UPTIME/ TTU0,TTU11,TTU12,TTU13,TTU2,TTU3,TTU4
C
      INCLUDE 'jump.inc'
      COMMON /NSCOUN/ NNONL,NMG
C
      INCLUDE 'bouss.inc'
      EXTERNAL DVISCO,DVISCOPRIM
      SAVE 
C
c      print *,'A1sup=',(A(i),i=1,25)
  
C
ccc   CALL ZTIME(TTT0)
c
       ICUB=ICUBNEWTON
C
      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
       CT0=DCMASS/THSTEP
      ELSE
       CT0=0D0
      ENDIF
C
      DO 1 I = 1,NNDER
1     BDER(I)=.FALSE.
C
      DO 2 I=1,3
2     BDER(I)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      DUMAX=0D0
      IF (A2L.EQ.0D0) THEN
       DO 10 IEQ=1,NMT
       DU1=U1L1(IEQ)
       DU2=U1L2(IEQ)
       DUNORM=SQRT(DU1**2+DU2**2)
       DUMAX=MAX(DUMAX,DUNORM)
10     CONTINUE
      ELSE       
       DO 20 IEQ=1,NMT
       DU1=A1L*U1L1(IEQ)+A2L*U2L1(IEQ)
       DU2=A1L*U1L2(IEQ)+A2L*U2L2(IEQ)
       DUNORM=SQRT(DU1**2+DU2**2)
       DUMAX=MAX(DUMAX,DUNORM)
20     CONTINUE
      ENDIF       
C
      IF (DUMAX.LT.1D-8) DUMAX=1D-8
C
C
ccc   CALL ZTIME(TTT1)
ccc   TTU0=TTU0+TTT1-TTT0
C *** Loop over all elements
      DO 100 IEL=1,NEL
ccc   CALL ZTIME(TTT0)
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine local DELTA for streamline-diffusion
      CALL DELTSD(U1L1,U1L2,U2L1,U2L2,A1L,A2L,IEL,DUMAX,DELTA)
c       IF (ISTOK.EQ.1 .OR. DNNS .NE. 0D0) DELTA=0D0
        IF (ISTOK.EQ.1) DELTA=0D0
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
         DENTRY(JDOFE,JDOFE,1)=0D0
         DENTRY(JDOFE,JDOFE,2)=0D0
         DENTRY(JDOFE,JDOFE,3)=0D0
         DENTRY(JDOFE,JDOFE,4)=0D0
            dhilf(jdofe,jdofe)=0d0
C
         DEJACOBNS(JDOFE,JDOFE,1)=0D0
         DEJACOBNS(JDOFE,JDOFE,2)=0D0
         DEJACOBNS(JDOFE,JDOFE,3)=0D0
         DEJACOBNS(JDOFE,JDOFE,4)=0D0
C 
         DEJACOBPW(JDOFE,JDOFE,1)=0D0
         DEJACOBPW(JDOFE,JDOFE,2)=0D0
         DEJACOBPW(JDOFE,JDOFE,3)=0D0
         DEJACOBPW(JDOFE,JDOFE,4)=0D0
C
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
          DENTRY(JDOFE,IDOFE,1)=0D0
          DENTRY(JDOFE,IDOFE,2)=0D0
          DENTRY(JDOFE,IDOFE,3)=0D0
          DENTRY(JDOFE,IDOFE,4)=0D0
            dhilf(jdofe,idofe)=0d0
C
         DEJACOBNS(JDOFE,IDOFE,1)=0D0
         DEJACOBNS(JDOFE,IDOFE,2)=0D0
         DEJACOBNS(JDOFE,IDOFE,3)=0D0
         DEJACOBNS(JDOFE,IDOFE,4)=0D0
C
         DEJACOBPW(JDOFE,IDOFE,1)=0D0
         DEJACOBPW(JDOFE,IDOFE,2)=0D0
         DEJACOBPW(JDOFE,IDOFE,3)=0D0
         DEJACOBPW(JDOFE,IDOFE,4)=0D0
C
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices
      DO 120 IVE = 1, NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
120   CONTINUE
ccc   CALL ZTIME(TTT1)
ccc   TTU11=TTU11+TTT1-TTT0
ccc   CALL ZTIME(TTT0)
C
      DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
      DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
      DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
      DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
C
ccc   CALL ZTIME(TTT1)
ccc   TTU2=TTU2+TTT1-TTT0
C *** Dummy call - ELE may save arithmetic operations
ccc   CALL ZTIME(TTT0)
      CALL ELE(0D0,0D0,-2)
      IF (IER.LT.0) GOTO 99999
ccc   CALL ZTIME(TTT1)
ccc   TTU3=TTU3+TTT1-TTT0
C
      avgny=0d0
C *** Loop over all cubature points
      DO 200 ICUBP = 1, NCUBP
ccc   CALL ZTIME(TTT0)
C
C *** Cubature points on the reference element
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
C
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
ccc   CALL ZTIME(TTT1)
ccc   TTU2=TTU2+TTT1-TTT0
ccc   CALL ZTIME(TTT0)
      CALL ELE(XX,YY,-3)
      IF (IER.LT.0) GOTO 99999
C
ccc   CALL ZTIME(TTT1)
ccc   TTU4=TTU4+TTT1-TTT0
ccc   CALL ZTIME(TTT0)
C *** Evaluation of velocity in cubature points
      IF (DCMASS.NE.0D0) THEN
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
C     and derivatives
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
c
C ***  norm of symmetric part of velocity gradient
       IF (ISPGRAD .EQ. 1)DUGSQ=DUX1**2+0.5d0*(DUX2+DUY1)**2+DUY2**2
c
C ***  norm of gradient
       IF (ISPGRAD .EQ. 0)DUGSQ=DUX1**2+DUX2**2+DUY1**2+DUY2**2
       DNY=DVISCO(DUGSQ)
       DNYPRIM=2D0*DVISCOPRIM(DUGSQ)
c
       avgny=avgny+dny/NY
C
        IF (ISPGRAD .EQ. 0)THEN
C ***  use just grad
C ***  Summing up over all pairs of multiindices
       DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
       HSUMJ=(HBASJ2*DU1+HBASJ3*DU2)!*0d0!!!WEG zu testzwecken
C
       HYUMJ=DUX1*HBASJ2+DUY1*HBASJ3
       HZUMJ=DUX2*HBASJ2+DUY2*HBASJ3
C
       DO 235 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH4= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
C
C
             AJ1= DUX1*HBASJ1**2
             AJ2= DUY1*HBASJ1**2
             AJ3= DUX2*HBASJ1**2
             AJ4= DUY2*HBASJ1**2 
c
             AK1=  DNYPRIM*HYUMJ**2
             AK2=  DNYPRIM*HZUMJ*HYUMJ
             AK3=  DNYPRIM*HYUMJ*HZUMJ
             AK4=  DNYPRIM*HZUMJ**2
c
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
        HSUMI=(HBASI2*DU1+HBASI3*DU2)!*0d0!!!WEG zu testzwecken
C
       HYUMI=DUX1*HBASI2+DUY1*HBASI3
       HZUMI=DUX2*HBASI2+DUY2*HBASI3
C
             AH1= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +   DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH4= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +   DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
C
             AJ1=  DUX1*HBASJ1*HBASI1
             AJ2=  DUY1*HBASJ1*HBASI1
             AJ3=  DUX2*HBASJ1*HBASI1
             AJ4=  DUY2*HBASJ1*HBASI1
c
             AK1=  DNYPRIM*HYUMI*HYUMJ
             AK2=  DNYPRIM*HZUMI*HYUMJ
             AK3=  DNYPRIM*HYUMI*HZUMJ
             AK4=  DNYPRIM*HZUMI*HZUMJ
C
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            DHILF (JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE,4)+OM*AH1           
c      IF (KENTRY(JDOFE,IDOFE).eq.5) print *,'su[', DENTRY(JDOFE,IDOFE,1)
C
         DEJACOBNS(JDOFE,IDOFE,1)=DEJACOBNS(JDOFE,IDOFE,1)+OM*AJ1
         DEJACOBNS(JDOFE,IDOFE,2)=DEJACOBNS(JDOFE,IDOFE,2)+OM*AJ2
         DEJACOBNS(JDOFE,IDOFE,3)=DEJACOBNS(JDOFE,IDOFE,3)+OM*AJ3
         DEJACOBNS(JDOFE,IDOFE,4)=DEJACOBNS(JDOFE,IDOFE,4)+OM*AJ4
C
         DEJACOBPW(JDOFE,IDOFE,1)=DEJACOBPW(JDOFE,IDOFE,1)+OM*AK1
         DEJACOBPW(JDOFE,IDOFE,2)=DEJACOBPW(JDOFE,IDOFE,2)+OM*AK2
         DEJACOBPW(JDOFE,IDOFE,3)=DEJACOBPW(JDOFE,IDOFE,3)+OM*AK3
         DEJACOBPW(JDOFE,IDOFE,4)=DEJACOBPW(JDOFE,IDOFE,4)+OM*AK4 
c            
235        CONTINUE
230    CONTINUE
C
       ELSE
C
C***  use symmetric part of grad
c
       DNYPRIM=2D0*DNYPRIM
C
C ***  Summing up over all pairs of multiindices
       DO 240 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2
C
       HYUMJ=DUX1*HBASJ2+0.5D0*(DUY1+DUX2)*HBASJ3
       HZUMJ=0.5D0*(DUY1+DUX2)*HBASJ2+DUY2*HBASJ3
C
       DO 245 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(2d0*HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH2= DNY*HBASJ2*HBASJ3
             AH3= DNY*HBASJ3*HBASJ2
             AH4= HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *         +  DNY*(HBASJ2**2+2d0*HBASJ3**2)+CT0*HBASJ1**2
             ahilf=DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
c
             AJ1=  DUX1*HBASJ1**2
             AJ2=  DUY1*HBASJ1**2
             AJ3=  DUX2*HBASJ1**2
             AJ4=  DUY2*HBASJ1**2
C
             AK1=  DNYPRIM*HYUMJ**2
             AK2=  DNYPRIM*HZUMJ*HYUMJ
             AK3=  DNYPRIM*HYUMJ*HZUMJ
             AK4=  DNYPRIM*HZUMJ**2
c
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
        HSUMI=HBASI2*DU1+HBASI3*DU2
C
       HYUMI=DUX1*HBASI2+0.5D0*(DUY1+DUX2)*HBASI3
       HZUMI=0.5D0*(DUY1+DUX2)*HBASI2+DUY2*HBASI3
C
             AH1= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +   DNY*(2d0*HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH2= DNY*HBASI2*HBASJ3!HBASI3*HBASJ2
             AH3= DNY*HBASI3*HBASJ2!HBASI2*HBASJ3
             AH4= HSUMI*(DELTA*HSUMJ+HBASJ1)
     *         +   DNY*(HBASJ2*HBASI2+2d0*HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             ahilf=DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +  CT0*HBASJ1*HBASI1
c
             AJ1=  DUX1*HBASJ1*HBASI1
             AJ2=  DUY1*HBASJ1*HBASI1
             AJ3=  DUX2*HBASJ1*HBASI1
             AJ4=  DUY2*HBASJ1*HBASI1
C
             AK1=  DNYPRIM*HYUMI*HYUMJ
             AK2=  DNYPRIM*HZUMI*HYUMJ
             AK3=  DNYPRIM*HYUMI*HZUMJ
             AK4=  DNYPRIM*HZUMI*HZUMJ
c
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,2)=DENTRY(JDOFE,IDOFE,2)+OM*AH2
            DENTRY(JDOFE,IDOFE,3)=DENTRY(JDOFE,IDOFE,3)+OM*AH3
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            DHILF (JDOFE,IDOFE)=DHILF (JDOFE,IDOFE)+OM*AHilf    
c
         DEJACOBNS(JDOFE,IDOFE,1)=DEJACOBNS(JDOFE,IDOFE,1)+OM*AJ1
         DEJACOBNS(JDOFE,IDOFE,2)=DEJACOBNS(JDOFE,IDOFE,2)+OM*AJ2
         DEJACOBNS(JDOFE,IDOFE,3)=DEJACOBNS(JDOFE,IDOFE,3)+OM*AJ3
         DEJACOBNS(JDOFE,IDOFE,4)=DEJACOBNS(JDOFE,IDOFE,4)+OM*AJ4
C
         DEJACOBPW(JDOFE,IDOFE,1)=DEJACOBPW(JDOFE,IDOFE,1)+OM*AK1
         DEJACOBPW(JDOFE,IDOFE,2)=DEJACOBPW(JDOFE,IDOFE,2)+OM*AK2
         DEJACOBPW(JDOFE,IDOFE,3)=DEJACOBPW(JDOFE,IDOFE,3)+OM*AK3
         DEJACOBPW(JDOFE,IDOFE,4)=DEJACOBPW(JDOFE,IDOFE,4)+OM*AK4 
c        
            
 245        CONTINUE
240    CONTINUE
       ENDIF ! ISPGRAD
c
c
ccc   CALL ZTIME(TTT1)
ccc   TTU12=TTU12+TTT1-TTT0
C
      ELSE
C
C ***  Summing up over 1 pair of multiindices
       DO 250 JDOFE=1,IDFL
       HBASJ1=DBAS(KDFL(JDOFE),1)
C
       DO 260 IDOFE=1,IDFL
       HBASI1=DBAS(KDFL(IDOFE),1)
C
       AH=-1D0/THSTEP*HBASJ1*HBASI1
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH
 
260    CONTINUE
250    CONTINUE
ccc   CALL ZTIME(TTT1)
ccc   TTU12=TTU12+TTT1-TTT0
C
      ENDIF
C
200   CONTINUE
C
C
ccc   CALL ZTIME(TTT0)
      DCONST=0D0
      locny(iel)=avgny/ncubp
      DO 400 JDOFE=1,IDFL
      DO 400 IDOFE=1,IDFL
          DENTH1=THSTEP*DENTRY(JDOFE,IDOFE,1)*DCONST
          DENTH2=THSTEP*DENTRY(JDOFE,IDOFE,2)*DCONST
          DENTH3=THSTEP*DENTRY(JDOFE,IDOFE,3)*DCONST
          DENTH4=THSTEP*DENTRY(JDOFE,IDOFE,4)*DCONST
          DHHILF=THSTEP*DHILF (JDOFE,IDOFE)*DCONST
C
          DENTJ1=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,1)
          DENTJ2=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,2)
          DENTJ3=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,3)
          DENTJ4=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,4)
C
          DENTK1=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,1)
          DENTK2=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,2)
          DENTK3=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,3)
          DENTK4=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,4)
C
          BSONST=.false.
C
      IF (IDEF.LT.2) THEN
       IA   =KENTRY(JDOFE,IDOFE)
         IF ((IPRECO.EQ.0).or.(BSONST)) THEN !Def-Tensor as Precond.
c            IF (ia.eq.1) write (*,*) ilev,IPRECO, MOD(IZBV2,IPRECO)
            A(IA)=A(IA)+REAL(DENTH1)+REAL(DENTJ1)+REAL(DENTK1)
           A(NA+IA)  =A(NA+IA)+REAL(DENTH2)+REAL(DENTJ2)+REAL(DENTK2)
           A(2*NA+IA)=A(2*NA+IA)+REAL(DENTH3)+REAL(DENTJ3)+REAL(DENTK3)
           A(3*NA+IA)=A(3*NA+IA)+REAL(DENTH4)+REAL(DENTJ4)+REAL(DENTK4)
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
c
          ELSEIF (IPRECO.EQ.1) THEN !GRAD as precond.
           A(IA)=A(IA)+REAL(DHHILF)
           A(NA+IA)=0d0
           A(2*NA+IA)=0d0
           A(3*NA+IA)=A(3*NA+IA)+REAL(DHHILF)
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
          ELSE
           A(IA)     =A(IA)     +0.5d0*(REAL(DHHILF)+REAL(DENTH1)
     *                          +REAL(DENTJ2))
           A(NA+IA)  =A(  NA+IA)+0.5d0*(0d0+REAL(DENTH2)
     *                          +REAL(DENTJ2))
           A(2*NA+IA)=A(2*NA+IA)+0.5d0*(0d0+REAL(DENTH3)
     *                          +REAL(DENTJ3))
           A(3*NA+IA)=A(3*NA+IA)+0.5d0*(REAL(DHHILF)+REAL(DENTH4)
     *                          +REAL(DENTJ4))
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
           ENDIf
      ENDIF

c$$$
c$$$         IF (IPRECO.EQ.0) THEN !Def-Tensor as Precond.
c$$$           A(IA)=A(IA)+REAL(DENTH1)
c$$$           A(NA+IA)=A(NA+IA)+REAL(DENTH2)
c$$$           A(2*NA+IA)=A(2*NA+IA)+REAL(DENTH3)
c$$$           A(3*NA+IA)=A(3*NA+IA)+REAL(DENTH4)
c$$$           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
c$$$c
c$$$          ELSE !GRAD as precond.
c$$$          A(IA)=A(IA)+REAL(DHHILF)
c$$$           A(NA+IA)=0d0
c$$$           A(2*NA+IA)=0d0
c$$$           A(3*NA+IA)=A(3*NA+IA)+REAL(DHHILF)
c$$$           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
c$$$           ENDIf
c$$$      ENDIF
C
      IF (IDEF.GT.0) THEN 
       IDFG=KDFG(IDOFE)
       JDFG=KDFG(JDOFE)
           D1(JDFG)= D1(JDFG)-DENTH1*U1(IDFG)-DENTH2*U2(IDFG)
           D2(JDFG)= D2(JDFG)-DENTH3*U1(IDFG)-DENTH4*U2(IDFG)
      ENDIF 
C
400   CONTINUE
C
ccc   CALL ZTIME(TTT1)
ccc   TTU13=TTU13+TTT1-TTT0
100   CONTINUE
C
ccc   write(6,*) TTU0,TTU11,TTU12,TTU13,TTU2,TTU3,TTU4,NCUBP
C
c      print *,ilev'Asup=',(A(i),i=1,25)

99999 END
C
C************************************************************************
      SUBROUTINE ADDSTPNEW(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *  KCOLA,KLDA,KVERT,KMID,DCORVG,ELE,COEFFN,
     *  IDEF,DCMASS,locny)
************************************************************************
*     Purpose: -  Adds the SUPG-part on matrix block A after
*                 it was initialized by the linear part
*              -  The input vector Ui is the old velocity field
*              -  The input vectors UjLi are the transport directions
*     PARAMETRIC VERSION
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4)
      PARAMETER (NNLEV=9)
C
      DOUBLE PRECISION locny
      dimension locny(*)
      DIMENSION A(*),KCOLA(*),KLDA(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS)
      DIMENSION DENTRY(NNBAS,NNBAS,4),DHILF(NNBAS,NNBAS)
      DIMENSION DEJACOBNS(NNBAS,NNBAS,4),DEJACOBPW(NNBAS,NNBAS,4)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
C
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
C
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
      INCLUDE 'bouss.inc'
      EXTERNAL DVISCO,DVISCOPRIM
      SAVE 
C
      ICUB=ICUBNEWTON
C
      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
       CT0=DCMASS/THSTEP
      ELSE
       CT0=0D0
      ENDIF
C
      DO 1 I = 1,NNDER
1     BDER(I)=.FALSE.
C
      DO 2 I=1,3
2     BDER(I)=.TRUE.
C
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************
      ICUBP=ICUB
      CALL ELE(0D0,0D0,-2)
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLDA(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
         DENTRY(JDOFE,JDOFE,1)=0D0
         DENTRY(JDOFE,JDOFE,2)=0D0
         DENTRY(JDOFE,JDOFE,3)=0D0
         DENTRY(JDOFE,JDOFE,4)=0D0
            dhilf(jdofe,jdofe)=0d0
C
         DEJACOBNS(JDOFE,JDOFE,1)=0D0
         DEJACOBNS(JDOFE,JDOFE,2)=0D0
         DEJACOBNS(JDOFE,JDOFE,3)=0D0
         DEJACOBNS(JDOFE,JDOFE,4)=0D0
C
         DEJACOBPW(JDOFE,JDOFE,1)=0D0
         DEJACOBPW(JDOFE,JDOFE,2)=0D0
         DEJACOBPW(JDOFE,JDOFE,3)=0D0
         DEJACOBPW(JDOFE,JDOFE,4)=0D0
C
      JCOL0=ILD
      DO 111 IDOFE=1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
          DENTRY(JDOFE,IDOFE,1)=0D0
          DENTRY(JDOFE,IDOFE,2)=0D0
          DENTRY(JDOFE,IDOFE,3)=0D0
          DENTRY(JDOFE,IDOFE,4)=0D0
            dhilf(jdofe,idofe)=0d0
C
         DEJACOBNS(JDOFE,IDOFE,1)=0D0
         DEJACOBNS(JDOFE,IDOFE,2)=0D0
         DEJACOBNS(JDOFE,IDOFE,3)=0D0
         DEJACOBNS(JDOFE,IDOFE,4)=0D0
C
         DEJACOBPW(JDOFE,IDOFE,1)=0D0
         DEJACOBPW(JDOFE,IDOFE,2)=0D0
         DEJACOBPW(JDOFE,IDOFE,3)=0D0
         DEJACOBPW(JDOFE,IDOFE,4)=0D0
C
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices
      DO 120 IVE = 1, NVE
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
      avgny=0d0
c
C *** Loop over all cubature points
      DO 200 ICUBP = 1, NCUBP
      XI1=DXI(ICUBP,1)
      XI2=DXI(ICUBP,2)
C *** Jacobian of the bilinear mapping onto the reference element
      DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
      DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
      DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
      DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
      DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
      OM=DOMEGA(ICUBP)*DETJ
C
      CALL ELE(XI1,XI2,-3)
      IF (IER.LT.0) GOTO 99999
C
C *** Evaluation of velocity in cubature points
      IF (DCMASS.NE.0D0) THEN
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
C     and derivatives
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
c
C ***  norm of symmetric part of velocity gradient
       IF (ISPGRAD .EQ. 1)DUGSQ=DUX1**2+0.5d0*(DUX2+DUY1)**2+DUY2**2
c
C ***  norm of gradient
       IF (ISPGRAD .EQ. 0)DUGSQ=DUX1**2+DUX2**2+DUY1**2+DUY2**2
       DNY=DVISCO(DUGSQ)
       DNYPRIM=2D0*DVISCOPRIM(DUGSQ)
C       DNY    =DVISCO(DUGSQ)
C       DNYPRIM=2.0D0*DVISCOPRIM(DUGSQ)
       avgny=avgny+dny/NY
C 
        IF (ISPGRAD .EQ. 0)THEN
C ***  use just grad
C ***  Summing up over all pairs of multiindices
       DO 230 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
C
       HYUMJ=DUX1*HBASJ2+DUY1*HBASJ3
       HZUMJ=DUX2*HBASJ2+DUY2*HBASJ3
c
       DO 235 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
c
             AH1= DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
C
             AJ1= DUX1*HBASJ1**2
             AJ2= DUY1*HBASJ1**2
             AJ3= DUX2*HBASJ1**2
             AJ4= DUY2*HBASJ1**2
c
             AK1=  DNYPRIM*HYUMJ**2
             AK2=  DNYPRIM*HZUMJ*HYUMJ
             AK3=  DNYPRIM*HYUMJ*HZUMJ
             AK4=  DNYPRIM*HZUMJ**2
c 
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
c
        HSUMI=HBASI2*DU1+HBASI3*DU2
             AH1= DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
C
       HYUMI=DUX1*HBASI2+DUY1*HBASI3
       HZUMI=DUX2*HBASI2+DUY2*HBASI3
c
             AJ1=  DUX1*HBASJ1*HBASI1
             AJ2=  DUY1*HBASJ1*HBASI1
             AJ3=  DUX2*HBASJ1*HBASI1
             AJ4=  DUY2*HBASJ1*HBASI1
c
             AK1=  DNYPRIM*HYUMI*HYUMJ
             AK2=  DNYPRIM*HZUMI*HYUMJ
             AK3=  DNYPRIM*HYUMI*HZUMJ
             AK4=  DNYPRIM*HZUMI*HZUMJ
C
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH1
            DHILF (JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE,4)+OM*AH1     
C
         DEJACOBNS(JDOFE,IDOFE,1)=DEJACOBNS(JDOFE,IDOFE,1)+OM*AJ1
         DEJACOBNS(JDOFE,IDOFE,2)=DEJACOBNS(JDOFE,IDOFE,2)+OM*AJ2
         DEJACOBNS(JDOFE,IDOFE,3)=DEJACOBNS(JDOFE,IDOFE,3)+OM*AJ3
         DEJACOBNS(JDOFE,IDOFE,4)=DEJACOBNS(JDOFE,IDOFE,4)+OM*AJ4
C
         DEJACOBPW(JDOFE,IDOFE,1)=DEJACOBPW(JDOFE,IDOFE,1)+OM*AK1
         DEJACOBPW(JDOFE,IDOFE,2)=DEJACOBPW(JDOFE,IDOFE,2)+OM*AK2
         DEJACOBPW(JDOFE,IDOFE,3)=DEJACOBPW(JDOFE,IDOFE,3)+OM*AK3
         DEJACOBPW(JDOFE,IDOFE,4)=DEJACOBPW(JDOFE,IDOFE,4)+OM*AK4  
c     
            
 235        CONTINUE
230    CONTINUE
C
       ELSE
C
C***  use symmetric part of grad
c
      DNYPRIM=2D0*DNYPRIM
C
C ***  Summing up over all pairs of multiindices
       DO 240 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2
C
       HYUMJ=DUX1*HBASJ2+0.5D0*(DUY1+DUX2)*HBASJ3
       HZUMJ=0.5D0*(DUY1+DUX2)*HBASJ2+DUY2*HBASJ3
c
       DO 245 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= DNY*(2d0*HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH2= DNY*HBASJ2*HBASJ3
             AH3= DNY*HBASJ3*HBASJ2
             AH4= DNY*(HBASJ2**2+2d0*HBASJ3**2)+CT0*HBASJ1**2
             ahilf=DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
c
c
             AJ1=  DUX1*HBASJ1**2
             AJ2=  DUY1*HBASJ1**2
             AJ3=  DUX2*HBASJ1**2
             AJ4=  DUY2*HBASJ1**2
C
             AK1=  DNYPRIM*HYUMJ**2
             AK2=  DNYPRIM*HZUMJ*HYUMJ
             AK3=  DNYPRIM*HYUMJ*HZUMJ
             AK4=  DNYPRIM*HZUMJ**2
c
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
        HSUMI=HBASI2*DU1+HBASI3*DU2
c
       HYUMI=DUX1*HBASI2+0.5D0*(DUY1+DUX2)*HBASI3
       HZUMI=0.5D0*(DUY1+DUX2)*HBASI2+DUY2*HBASI3
c
             AH1=  DNY*(2d0*HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH2= DNY*HBASI2*HBASJ3!HBASI3*HBASJ2
             AH3= DNY*HBASI3*HBASJ2!HBASI2*HBASJ3
             AH4= DNY*(HBASJ2*HBASI2+2d0*HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             ahilf=DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +  CT0*HBASJ1*HBASI1
c
             AJ1=  DUX1*HBASJ1*HBASI1
             AJ2=  DUY1*HBASJ1*HBASI1
             AJ3=  DUX2*HBASJ1*HBASI1
             AJ4=  DUY2*HBASJ1*HBASI1
C
             AK1=  DNYPRIM*HYUMI*HYUMJ
             AK2=  DNYPRIM*HZUMI*HYUMJ
             AK3=  DNYPRIM*HYUMI*HZUMJ
             AK4=  DNYPRIM*HZUMI*HZUMJ
c
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,2)=DENTRY(JDOFE,IDOFE,2)+OM*AH2
            DENTRY(JDOFE,IDOFE,3)=DENTRY(JDOFE,IDOFE,3)+OM*AH3
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            DHILF (JDOFE,IDOFE)=DHILF (JDOFE,IDOFE)+OM*AHilf  
c
         DEJACOBNS(JDOFE,IDOFE,1)=DEJACOBNS(JDOFE,IDOFE,1)+OM*AJ1
         DEJACOBNS(JDOFE,IDOFE,2)=DEJACOBNS(JDOFE,IDOFE,2)+OM*AJ2
         DEJACOBNS(JDOFE,IDOFE,3)=DEJACOBNS(JDOFE,IDOFE,3)+OM*AJ3
         DEJACOBNS(JDOFE,IDOFE,4)=DEJACOBNS(JDOFE,IDOFE,4)+OM*AJ4
C
         DEJACOBPW(JDOFE,IDOFE,1)=DEJACOBPW(JDOFE,IDOFE,1)+OM*AK1
         DEJACOBPW(JDOFE,IDOFE,2)=DEJACOBPW(JDOFE,IDOFE,2)+OM*AK2
         DEJACOBPW(JDOFE,IDOFE,3)=DEJACOBPW(JDOFE,IDOFE,3)+OM*AK3
         DEJACOBPW(JDOFE,IDOFE,4)=DEJACOBPW(JDOFE,IDOFE,4)+OM*AK4 
c                      
 245        CONTINUE
240    CONTINUE
       ENDIF ! ISPGRAD  
C   
      ELSE
C
C ***  Summing up over 1 pair of multiindices
       DO 250 JDOFE=1,IDFL
       HBASJ1=DBAS(KDFL(JDOFE),1)
C
       DO 260 IDOFE=1,IDFL
       HBASI1=DBAS(KDFL(IDOFE),1)
C
       AH=-1D0/THSTEP*HBASJ1*HBASI1
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH
260    CONTINUE
250    CONTINUE
C
      ENDIF
C
200   CONTINUE
C
      DCONST=0D0
      locny(iel)=avgny/ncubp
c
      DO 400 JDOFE=1,IDFL
      DO 400 IDOFE=1,IDFL
          DENTH1=THSTEP*DENTRY(JDOFE,IDOFE,1)*DCONST
          DENTH2=THSTEP*DENTRY(JDOFE,IDOFE,2)*DCONST
          DENTH3=THSTEP*DENTRY(JDOFE,IDOFE,3)*DCONST
          DENTH4=THSTEP*DENTRY(JDOFE,IDOFE,4)*DCONST
          DHHILF=THSTEP*DHILF (JDOFE,IDOFE)*DCONST
C
          DENTJ1=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,1)
          DENTJ2=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,2)
          DENTJ3=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,3)
          DENTJ4=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,4)
C
          DENTK1=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,1)
          DENTK2=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,2)
          DENTK3=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,3)
          DENTK4=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,4)
c
      IF (IDEF.LT.2) THEN
       IA   =KENTRY(JDOFE,IDOFE)
         IF (IPRECO.EQ.0) THEN !Def-Tensor as Precond.
           A(IA)=A(IA)+REAL(DENTH1)+REAL(DENTJ1)+REAL(DENTK1)
           A(NA+IA)=A(NA+IA)+REAL(DENTH2)+REAL(DENTJ2)+REAL(DENTK2)
           A(2*NA+IA)=A(2*NA+IA)+REAL(DENTH3)+REAL(DENTJ3)+REAL(DENTK3)
           A(3*NA+IA)=A(3*NA+IA)+REAL(DENTH4)+REAL(DENTJ4)+REAL(DENTK4)
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
c
          ELSE !GRAD as precond.
          A(IA)=A(IA)+REAL(DHHILF)
           A(NA+IA)=0d0
           A(2*NA+IA)=0d0
           A(3*NA+IA)=A(3*NA+IA)+REAL(DHHILF)
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
           ENDIf
      ENDIF
C
      IF (IDEF.GT.0) THEN 
       IDFG=KDFG(IDOFE)
       JDFG=KDFG(JDOFE)
           D1(JDFG)= D1(JDFG)-DENTH1*U1(IDFG)-DENTH2*U2(IDFG)
           D2(JDFG)= D2(JDFG)-DENTH3*U1(IDFG)-DENTH4*U2(IDFG)
      ENDIF 
C
400   CONTINUE
C
100   CONTINUE
C   
99999  END
C
C     
************************************************************************
      SUBROUTINE ADDSTNNEW(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *  KCOLA,KLDA,KVERT,KMID,DCORVG,ELE,COEFFN,
     *  IDEF,DCMASS,locny)
************************************************************************
*     Purpose: -  Adds the M/ST-part on matrix block A 
*     -  The input vector Ui is the old velocity field
*     NONPARAMETRIC VERSION
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C     
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNCOF=6)
      PARAMETER (NNLEV=9)
C     
      DOUBLE PRECISION locny
      dimension locny(*)
      DIMENSION A(*),KCOLA(*),KLDA(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS)
      DIMENSION DENTRY(NNBAS,NNBAS,4),DHILF(NNBAS,NNBAS)
      DIMENSION DEJACOBNS(NNBAS,NNBAS,4),DEJACOBPW(NNBAS,NNBAS,4)
C     
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *  DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      COMMON /COFBAS/ COB(NNBAS,NNCOF)
C     
C***  COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
C     
C***  user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *  IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *  INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *  NSM,NSL,NSMFAC
C     
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *  EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *  RLXSM,RLXSL,AMINMG,AMAXMG
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
C
      INCLUDE 'jump.inc'
      INCLUDE 'bouss.inc'
      EXTERNAL DVISCO,DVISCOPRIM
      SAVE 
C    
      CALL VKADD(AVK1,AVK2)
C
ccc   CALL ZTIME(TTT0)
      ICUB=ICUBNEWTON
c
c      print *,ilev,'A1=',(A(i),i=1,25)
      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
       CT0=DCMASS/THSTEP
      ELSE
       CT0=0D0
      ENDIF
C     
C     
      DO 1 I = 1,NNDER
       BDER(I)=.FALSE.
 1    continue
C     
      DO 2 I=1,3
       BDER(I)=.TRUE.
 2    continue
C     
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C     
************************************************************************
C***  Calculation of the matrix - storage technique 7 or 8
************************************************************************
C     
C     
ccc   CALL ZTIME(TTT1)
ccc   TTU0=TTU0+TTT1-TTT0
C***  Loop over all elements
       DO 100 IEL=1,NEL
ccc   CALL ZTIME(TTT0)
        CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999
C     
C     
C***  Determine entry positions in matrix
        DO 110 JDOFE=1,IDFL
         ILD=KLDA(KDFG(JDOFE))
         KENTRY(JDOFE,JDOFE)=ILD
         DENTRY(JDOFE,JDOFE,1)=0D0
         DENTRY(JDOFE,JDOFE,2)=0D0
         DENTRY(JDOFE,JDOFE,3)=0D0
         DENTRY(JDOFE,JDOFE,4)=0D0
            dhilf(jdofe,jdofe)=0d0
C
         DEJACOBNS(JDOFE,JDOFE,1)=0D0
         DEJACOBNS(JDOFE,JDOFE,2)=0D0
         DEJACOBNS(JDOFE,JDOFE,3)=0D0
         DEJACOBNS(JDOFE,JDOFE,4)=0D0
C
         DEJACOBPW(JDOFE,JDOFE,1)=0D0
         DEJACOBPW(JDOFE,JDOFE,2)=0D0
         DEJACOBPW(JDOFE,JDOFE,3)=0D0
         DEJACOBPW(JDOFE,JDOFE,4)=0D0
C
C
         JCOL0=ILD
         DO 111 IDOFE=1,IDFL
          IF (IDOFE.EQ.JDOFE) GOTO 111
          IDFG=KDFG(IDOFE)
          DO 112 JCOL=JCOL0,NA
           IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
 112      CONTINUE
 113      JCOL0=JCOL+1
          KENTRY(JDOFE,IDOFE)=JCOL
          DENTRY(JDOFE,IDOFE,1)=0D0
          DENTRY(JDOFE,IDOFE,2)=0D0
          DENTRY(JDOFE,IDOFE,3)=0D0
          DENTRY(JDOFE,IDOFE,4)=0D0
            dhilf(jdofe,idofe)=0d0
C
         DEJACOBNS(JDOFE,IDOFE,1)=0D0
         DEJACOBNS(JDOFE,IDOFE,2)=0D0
         DEJACOBNS(JDOFE,IDOFE,3)=0D0
         DEJACOBNS(JDOFE,IDOFE,4)=0D0
C 
         DEJACOBPW(JDOFE,IDOFE,1)=0D0
         DEJACOBPW(JDOFE,IDOFE,2)=0D0
         DEJACOBPW(JDOFE,IDOFE,3)=0D0
         DEJACOBPW(JDOFE,IDOFE,4)=0D0
C
 111     CONTINUE
 110    CONTINUE
C     
C***  Evaluation of coordinates of the vertices
        BRANDE=.FALSE.
        DO 120 IVE = 1, NVE
         JP=KVERT(IVE,IEL)
         KVE(IVE)=JP
         DX(IVE)=DCORVG(1,JP)
         DY(IVE)=DCORVG(2,JP)
c =========================================
c     Hier mache ich den Randcheck!
c =========================================
c        IF ((DX(IVE)*DY(IVE).eq.1d0).or.(DX(IVE)*DY(IVE).eq.0d0))
        IF ((DX(IVE)).ge.0.5d0)!.or.(DX(IVE)*DY(IVE).eq.0d0))
     *  BRANDE=.TRUE.! DANN MACHE ICH, BEI IGRAD/IPRECO=2 DEF, Sonst GRAD
c        BRANDE=.TRUE.
c
c =========================================
 120    CONTINUE
ccc   CALL ZTIME(TTT1)
ccc   TTU11=TTU11+TTT1-TTT0
ccc   CALL ZTIME(TTT0)
C     

        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
C     
ccc   CALL ZTIME(TTT1)
ccc   TTU2=TTU2+TTT1-TTT0
C***  Dummy call - ELE may save arithmetic operations
ccc   CALL ZTIME(TTT0)
        CALL ELE(0D0,0D0,-2)
        IF (IER.LT.0) GOTO 99999
ccc   CALL ZTIME(TTT1)
ccc   TTU3=TTU3+TTT1-TTT0
C
      avgny=0d0    
c 
C***  Loop over all cubature points
        DO 200 ICUBP = 1, NCUBP
ccc   CALL ZTIME(TTT0)
C     
C***  Cubature points on the reference element
         XI1=DXI(ICUBP,1)
         XI2=DXI(ICUBP,2)
C     
C***  Jacobian of the bilinear mapping onto the reference element
         DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
         DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
         DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
         DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
         DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
C     
C***  Cubature points on actual Element + weights
         XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *     +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
         YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1
     *     +0.5D0*(DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2
         OM=DOMEGA(ICUBP)*DETJ
C  
ccc   CALL ZTIME(TTT1)
ccc   TTU2=TTU2+TTT1-TTT0
ccc   CALL ZTIME(TTT0)
         CALL ELE(XX,YY,-3)
         IF (IER.LT.0) GOTO 99999
C     
ccc   CALL ZTIME(TTT1)
ccc   TTU4=TTU4+TTT1-TTT0
ccc   CALL ZTIME(TTT0)
C***  Evaluation of velocity in cubature points
         IF (DCMASS.NE.0D0) THEN
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
c
C     and derivatives
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
c
C ***  norm of symmetric part of velocity gradient
       IF (ISPGRAD .EQ. 1)DUGSQ=DUX1**2+0.5d0*(DUX2+DUY1)**2+DUY2**2
c
C ***  norm of gradient
       IF (ISPGRAD .EQ. 0)DUGSQ=DUX1**2+DUX2**2+DUY1**2+DUY2**2
       DNY=DVISCO(DUGSQ)
       DNYPRIM=2D0*DVISCOPRIM(DUGSQ)
       avgny=avgny+dny/NY
C 
        IF (ISPGRAD .EQ. 0)THEN
C *** use grad
C***  Summing up over all pairs of multiindices
          DO 230 JDOFE=1,IDFL
           JDOFEH=KDFL(JDOFE)
           HBASJ1=DBAS(JDOFEH,1)
           HBASJ2=DBAS(JDOFEH,2)
           HBASJ3=DBAS(JDOFEH,3)
C  
           HYUMJ=DUX1*HBASJ2+DUY1*HBASJ3
           HZUMJ=DUX2*HBASJ2+DUY2*HBASJ3   
C     
           DO 235 IDOFE=1,IDFL
            IF (IDOFE.EQ.JDOFE) THEN
             AH1= DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
C
             AJ1= DUX1*HBASJ1**2
             AJ2= DUY1*HBASJ1**2
             AJ3= DUX2*HBASJ1**2
             AJ4= DUY2*HBASJ1**2 
c
             AK1=  DNYPRIM*HYUMJ**2
             AK2=  DNYPRIM*HZUMJ*HYUMJ
             AK3=  DNYPRIM*HYUMJ*HZUMJ
             AK4=  DNYPRIM*HZUMJ**2
c
            ELSE
             IDOFEH=KDFL(IDOFE)
             HBASI1=DBAS(IDOFEH,1)
             HBASI2=DBAS(IDOFEH,2)
             HBASI3=DBAS(IDOFEH,3)
c
             HYUMI=DUX1*HBASI2+DUY1*HBASI3
             HZUMI=DUX2*HBASI2+DUY2*HBASI3
c
c     
c     J -test function   I -trial function
c     
             AH1= DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +  CT0*HBASJ1*HBASI1
C
             AJ1=  DUX1*HBASJ1*HBASI1
             AJ2=  DUY1*HBASJ1*HBASI1
             AJ3=  DUX2*HBASJ1*HBASI1
             AJ4=  DUY2*HBASJ1*HBASI1
c
             AK1=  DNYPRIM*HYUMI*HYUMJ
             AK2=  DNYPRIM*HZUMI*HYUMJ
             AK3=  DNYPRIM*HYUMI*HZUMJ
             AK4=  DNYPRIM*HZUMI*HZUMJ
C
            ENDIF
C     
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,2)=DENTRY(JDOFE,IDOFE,2)+0d0
            DENTRY(JDOFE,IDOFE,3)=DENTRY(JDOFE,IDOFE,3)+0d0
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH1
C
         DEJACOBNS(JDOFE,IDOFE,1)=DEJACOBNS(JDOFE,IDOFE,1)+OM*AJ1
         DEJACOBNS(JDOFE,IDOFE,2)=DEJACOBNS(JDOFE,IDOFE,2)+OM*AJ2
         DEJACOBNS(JDOFE,IDOFE,3)=DEJACOBNS(JDOFE,IDOFE,3)+OM*AJ3
         DEJACOBNS(JDOFE,IDOFE,4)=DEJACOBNS(JDOFE,IDOFE,4)+OM*AJ4 
c
         DEJACOBPW(JDOFE,IDOFE,1)=DEJACOBPW(JDOFE,IDOFE,1)+OM*AK1
         DEJACOBPW(JDOFE,IDOFE,2)=DEJACOBPW(JDOFE,IDOFE,2)+OM*AK2
         DEJACOBPW(JDOFE,IDOFE,3)=DEJACOBPW(JDOFE,IDOFE,3)+OM*AK3
         DEJACOBPW(JDOFE,IDOFE,4)=DEJACOBPW(JDOFE,IDOFE,4)+OM*AK4 
C
            DHILF (JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE,4)+OM*AH1           
c       IF (KENTRY(JDOFE,IDOFE).eq.5) print *,iel,'up',
c     * DENTRY(JDOFE,IDOFE,1)
           
 235     CONTINUE
 230      CONTINUE
C
       ELSE
C
C***  use symmetric part of grad
c     
       DNYPRIM=2D0*DNYPRIM
C
C ***  Summing up over all pairs of multiindices
       DO 240 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2
C
       HYUMJ=DUX1*HBASJ2+0.5D0*(DUY1+DUX2)*HBASJ3
       HZUMJ=0.5D0*(DUY1+DUX2)*HBASJ2+DUY2*HBASJ3
c
       DO 245 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= DNY*(2d0*HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH2= DNY*HBASJ2*HBASJ3
             AH3= DNY*HBASJ3*HBASJ2
             AH4= DNY*(HBASJ2**2+2d0*HBASJ3**2)+CT0*HBASJ1**2
             ahilf=DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
c
             AJ1=  DUX1*HBASJ1**2
             AJ2=  DUY1*HBASJ1**2
             AJ3=  DUX2*HBASJ1**2
             AJ4=  DUY2*HBASJ1**2
C
             AK1=  DNYPRIM*HYUMJ**2
             AK2=  DNYPRIM*HZUMJ*HYUMJ
             AK3=  DNYPRIM*HYUMJ*HZUMJ
             AK4=  DNYPRIM*HZUMJ**2
c
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
        HSUMI=HBASI2*DU1+HBASI3*DU2
c
       HYUMI=DUX1*HBASI2+0.5D0*(DUY1+DUX2)*HBASI3
       HZUMI=0.5D0*(DUY1+DUX2)*HBASI2+DUY2*HBASI3
c
             AH1= DNY*(2d0*HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH2= DNY*HBASI2*HBASJ3!HBASI3*HBASJ2
             AH3= DNY*HBASI3*HBASJ2!HBASI2*HBASJ3
             AH4= DNY*(HBASJ2*HBASI2+2d0*HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             ahilf=DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +  CT0*HBASJ1*HBASI1
c
             AJ1=  DUX1*HBASJ1*HBASI1
             AJ2=  DUY1*HBASJ1*HBASI1
             AJ3=  DUX2*HBASJ1*HBASI1
             AJ4=  DUY2*HBASJ1*HBASI1
C
             AK1=  DNYPRIM*HYUMI*HYUMJ
             AK2=  DNYPRIM*HZUMI*HYUMJ
             AK3=  DNYPRIM*HYUMI*HZUMJ
             AK4=  DNYPRIM*HZUMI*HZUMJ
c
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,2)=DENTRY(JDOFE,IDOFE,2)+OM*AH2
            DENTRY(JDOFE,IDOFE,3)=DENTRY(JDOFE,IDOFE,3)+OM*AH3
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            DHILF (JDOFE,IDOFE)=DHILF (JDOFE,IDOFE)+OM*AHilf  
c
         DEJACOBNS(JDOFE,IDOFE,1)=DEJACOBNS(JDOFE,IDOFE,1)+OM*AJ1
         DEJACOBNS(JDOFE,IDOFE,2)=DEJACOBNS(JDOFE,IDOFE,2)+OM*AJ2
         DEJACOBNS(JDOFE,IDOFE,3)=DEJACOBNS(JDOFE,IDOFE,3)+OM*AJ3
         DEJACOBNS(JDOFE,IDOFE,4)=DEJACOBNS(JDOFE,IDOFE,4)+OM*AJ4
C
         DEJACOBPW(JDOFE,IDOFE,1)=DEJACOBPW(JDOFE,IDOFE,1)+OM*AK1
         DEJACOBPW(JDOFE,IDOFE,2)=DEJACOBPW(JDOFE,IDOFE,2)+OM*AK2
         DEJACOBPW(JDOFE,IDOFE,3)=DEJACOBPW(JDOFE,IDOFE,3)+OM*AK3
         DEJACOBPW(JDOFE,IDOFE,4)=DEJACOBPW(JDOFE,IDOFE,4)+OM*AK4 
c          
            
 245        CONTINUE
240    CONTINUE
       ENDIF ! ISPGRAD
C     
         ELSE
c     dcmass.eq.0d0
C     
C***  Summing up over 1 pair of multiindices
          DO 250 JDOFE=1,IDFL
           HBASJ1=DBAS(KDFL(JDOFE),1)
C     
           DO 260 IDOFE=1,IDFL
            HBASI1=DBAS(KDFL(IDOFE),1)
C     
            AH=-1D0/THSTEP*HBASJ1*HBASI1
C     
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH
 260       CONTINUE
 250      CONTINUE

C     
         ENDIF
C     
 200    CONTINUE
C
      locny(iel)=avgny/ncubp
      DCONST=0D0
C     
ccc   CALL ZTIME(TTT0)
        DO 400 JDOFE=1,IDFL
         DO 400 IDOFE=1,IDFL
          DENTH1=THSTEP*DENTRY(JDOFE,IDOFE,1)*DCONST
          DENTH2=THSTEP*DENTRY(JDOFE,IDOFE,2)*DCONST
          DENTH3=THSTEP*DENTRY(JDOFE,IDOFE,3)*DCONST
          DENTH4=THSTEP*DENTRY(JDOFE,IDOFE,4)*DCONST
          DHHILF=THSTEP*DHILF (JDOFE,IDOFE)*DCONST
C
          DENTJ1=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,1)
          DENTJ2=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,2)
          DENTJ3=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,3)
          DENTJ4=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,4)
c
          DENTK1=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,1)
          DENTK2=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,2)
          DENTK3=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,3)
          DENTK4=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,4)
C     
          BSONST=.false.
          goto 9911
c          IF ((IPRECO.EQ.0).or.(IPRECO.eq.1)) GOTO 9911
c          IF (ILEV.eq.NLMAX) THEN
c             IF ((MOD(IZBV2+1,IPRECO).eq.0)) BSONST=.true.
c          ELSE
c             IF ((MOD(IZBV2,IPRECO).eq.0)) BSONST=.true.
c          ENDIF
c          IF (IZBV2.le.5) BSONST=.true.
 9911     CONTINUE
      IF (IDEF.LT.2) THEN
       IA   =KENTRY(JDOFE,IDOFE)
         IF ((IPRECO.EQ.0).or.(BSONST)) THEN !Def-Tensor as Precond.
           A(IA)=A(IA)+REAL(DENTH1)+REAL(DENTJ1)+REAL(DENTK1)
           A(NA+IA)=A(NA+IA)+REAL(DENTH2)+REAL(DENTJ2)+REAL(DENTK2)
           A(2*NA+IA)=A(2*NA+IA)+REAL(DENTH3)+REAL(DENTJ3)+REAL(DENTK3)
           A(3*NA+IA)=A(3*NA+IA)+REAL(DENTH4)+REAL(DENTJ4)+REAL(DENTK4)
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
c
          ELSEIF (IPRECO.EQ.1) THEN !GRAD as precond.
           A(IA)=A(IA)+REAL(DHHILF)
           A(NA+IA)=0d0
           A(2*NA+IA)=0d0
           A(3*NA+IA)=A(3*NA+IA)+REAL(DHHILF)
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
          ELSE
c        IF (brande) print *,'up final', denth2,denth1,denth4
          A(IA)=A(IA)+REAL(DENTH1)+REAL(DENTJ1)
           A(NA+IA)  =A(NA+IA)+REAL(DENTH2)+REAL(DENTJ2)
           A(2*NA+IA)=A(2*NA+IA)+REAL(DENTH3)+REAL(DENTJ3)
           A(3*NA+IA)=A(3*NA+IA)+REAL(DENTH4)+REAL(DENTJ4)
           A(4*NA+IA)=A(4*NA+IA)+REAL(DHHILF)!Vorkonditionierer
           ENDIf
      ENDIF
C     
          IF (IDEF.GT.0) THEN 
           IDFG=KDFG(IDOFE)
           JDFG=KDFG(JDOFE)
           D1(JDFG)= D1(JDFG)-DENTH1*U1(IDFG)-DENTH2*U2(IDFG)
           D2(JDFG)= D2(JDFG)-DENTH3*U1(IDFG)-DENTH4*U2(IDFG)
c           write (*,*) 'in addstn', denth1
          ENDIF 
C     
 400     CONTINUE
C     
ccc   CALL ZTIME(TTT1)
ccc   TTU13=TTU13+TTT1-TTT0
 100    CONTINUE
C     
ccc   write(6,*) TTU0,TTU11,TTU12,TTU13,TTU2,TTU3,TTU4,NCUBP
C     
cc       if (IGRAD.eq.2) THEN
cc          IGRAD=ig1
cc          IPRECO=ig2
cc       ENDIF
c  
c      print *,ilev,'A=',(A(i),i=1,25)
C     
        
99999  END
C
C************************************************************************
      SUBROUTINE ADDSN1NEW(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *  KCOLA,KLDA,KVERT,KMID,DCORVG,ELE,COEFFN,
     *  IDEF,DCMASS,locny)
************************************************************************
*     Purpose: -  Adds the M/ST-part on matrix block A 
*     -  The input vector Ui is the old velocity field
*     NONPARAMETRIC VERSION. Hier: GRAD-GRAD TEIL
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C     
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNCOF=6)
      PARAMETER (NNLEV=9)
C     
      DOUBLE PRECISION locny
      dimension locny(*)
      DIMENSION A(*),KCOLA(*),KLDA(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KENTRY(NNBAS,NNBAS)
      DIMENSION DENTRY(NNBAS,NNBAS,4),DHILF(NNBAS,NNBAS)
      DIMENSION DEJACOBNS(NNBAS,NNBAS,4),DEJACOBPW(NNBAS,NNBAS,4)
C     
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *  DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      COMMON /COFBAS/ COB(NNBAS,NNCOF)
C     
C***  COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
C     
C***  user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *  IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *  INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *  NSM,NSL,NSMFAC
C     
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *  EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *  RLXSM,RLXSL,AMINMG,AMAXMG
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
C
      INCLUDE 'jump.inc'
      INCLUDE 'bouss.inc'
      EXTERNAL DVISCO,DVISCOPRIM
      SAVE 
C 
      CALL VKADD(AVK1,AVK2)
C
      ICUB=ICUBNEWTON
C     
      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
       CT0=DCMASS/THSTEP
      ELSE
       CT0=0D0
      ENDIF
C     
C     
      DO 1 I = 1,NNDER
       BDER(I)=.FALSE.
 1    continue
C     
      DO 2 I=1,3
       BDER(I)=.TRUE.
 2    continue
C     
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C     
************************************************************************
C***  Calculation of the matrix - storage technique 7 or 8
************************************************************************
C     
C     
       DO 100 IEL=1,NEL
c
        CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999
C     
C     
C***  Determine entry positions in matrix
        DO 110 JDOFE=1,IDFL
         ILD=KLDA(KDFG(JDOFE))
         KENTRY(JDOFE,JDOFE)=ILD
         DENTRY(JDOFE,JDOFE,1)=0D0
         DENTRY(JDOFE,JDOFE,2)=0D0
         DENTRY(JDOFE,JDOFE,3)=0D0
         DENTRY(JDOFE,JDOFE,4)=0D0
            dhilf(jdofe,jdofe)=0d0
C
         DEJACOBNS(JDOFE,JDOFE,1)=0D0
         DEJACOBNS(JDOFE,JDOFE,2)=0D0
         DEJACOBNS(JDOFE,JDOFE,3)=0D0
         DEJACOBNS(JDOFE,JDOFE,4)=0D0
C
         DEJACOBPW(JDOFE,JDOFE,1)=0D0
         DEJACOBPW(JDOFE,JDOFE,2)=0D0
         DEJACOBPW(JDOFE,JDOFE,3)=0D0
         DEJACOBPW(JDOFE,JDOFE,4)=0D0
C
         JCOL0=ILD
         DO 111 IDOFE=1,IDFL
          IF (IDOFE.EQ.JDOFE) GOTO 111
          IDFG=KDFG(IDOFE)
          DO 112 JCOL=JCOL0,NA
           IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
 112      CONTINUE
 113      JCOL0=JCOL+1
          KENTRY(JDOFE,IDOFE)=JCOL
          DENTRY(JDOFE,IDOFE,1)=0D0
          DENTRY(JDOFE,IDOFE,2)=0D0
          DENTRY(JDOFE,IDOFE,3)=0D0
          DENTRY(JDOFE,IDOFE,4)=0D0
            dhilf(jdofe,idofe)=0d0
C
         DEJACOBNS(JDOFE,IDOFE,1)=0D0
         DEJACOBNS(JDOFE,IDOFE,2)=0D0
         DEJACOBNS(JDOFE,IDOFE,3)=0D0
         DEJACOBNS(JDOFE,IDOFE,4)=0D0
C 
         DEJACOBPW(JDOFE,IDOFE,1)=0D0
         DEJACOBPW(JDOFE,IDOFE,2)=0D0
         DEJACOBPW(JDOFE,IDOFE,3)=0D0
         DEJACOBPW(JDOFE,IDOFE,4)=0D0
C
 111     CONTINUE
 110    CONTINUE
C     
C***  Evaluation of coordinates of the vertices
        DO 120 IVE = 1, NVE
         JP=KVERT(IVE,IEL)
         KVE(IVE)=JP
         DX(IVE)=DCORVG(1,JP)
         DY(IVE)=DCORVG(2,JP)
 120    CONTINUE

C     
        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
C     
C***  Dummy call - ELE may save arithmetic operations
c
        CALL ELE(0D0,0D0,-2)
        IF (IER.LT.0) GOTO 99999
C
      avgny=0d0 
c    
C***  Loop over all cubature points
        DO 200 ICUBP = 1, NCUBP
c
C     
C***  Cubature points on the reference element
         XI1=DXI(ICUBP,1)
         XI2=DXI(ICUBP,2)
C     
C***  Jacobian of the bilinear mapping onto the reference element
         DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
         DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
         DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
         DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
         DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
C     
C***  Cubature points on actual Element + weights
         XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *     +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
         YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1
     *     +0.5D0*(DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2
         OM=DOMEGA(ICUBP)*DETJ
C     
         CALL ELE(XX,YY,-3)
         IF (IER.LT.0) GOTO 99999
C     
C***  Evaluation of velocity in cubature points
         IF (DCMASS.NE.0D0) THEN
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
c
C     and derivatives
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
c
C ***  norm of symmetric part of velocity gradient
       IF (ISPGRAD .EQ. 1)DUGSQ=DUX1**2+0.5d0*(DUX2+DUY1)**2+DUY2**2
c
C ***  norm of gradient
       IF (ISPGRAD .EQ. 0)DUGSQ=DUX1**2+DUX2**2+DUY1**2+DUY2**2
       DNY=DVISCO(DUGSQ)
       DNYPRIM=2D0*DVISCOPRIM(DUGSQ)
       avgny=avgny+dny/NY
C
        IF (ISPGRAD .EQ. 0)THEN
c     use grad
C***  Summing up over all pairs of multiindices

          DO 230 JDOFE=1,IDFL
           JDOFEH=KDFL(JDOFE)
           HBASJ1=DBAS(JDOFEH,1)
           HBASJ2=DBAS(JDOFEH,2)
           HBASJ3=DBAS(JDOFEH,3)
C     
           HYUMJ=DUX1*HBASJ2+DUY1*HBASJ3
           HZUMJ=DUX2*HBASJ2+DUY2*HBASJ3
C     
           DO 235 IDOFE=1,IDFL
            IF (IDOFE.EQ.JDOFE) THEN
             AH1= DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
C
             AJ1= DUX1*HBASJ1**2
             AJ2= DUY1*HBASJ1**2
             AJ3= DUX2*HBASJ1**2
             AJ4= DUY2*HBASJ1**2 
c
             AK1=  DNYPRIM*HYUMJ**2
             AK2=  DNYPRIM*HZUMJ*HYUMJ
             AK3=  DNYPRIM*HYUMJ*HZUMJ
             AK4=  DNYPRIM*HZUMJ**2
C
            ELSE
             IDOFEH=KDFL(IDOFE)
             HBASI1=DBAS(IDOFEH,1)
             HBASI2=DBAS(IDOFEH,2)
             HBASI3=DBAS(IDOFEH,3)
c
             HYUMI=DUX1*HBASI2+DUY1*HBASI3
             HZUMI=DUX2*HBASI2+DUY2*HBASI3
c     
c     J -test function   I -trial function
c     
             AH1= DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +  CT0*HBASJ1*HBASI1
C
             AJ1=  DUX1*HBASJ1*HBASI1
             AJ2=  DUY1*HBASJ1*HBASI1
             AJ3=  DUX2*HBASJ1*HBASI1
             AJ4=  DUY2*HBASJ1*HBASI1
c
             AK1=  DNYPRIM*HYUMI*HYUMJ
             AK2=  DNYPRIM*HZUMI*HYUMJ
             AK3=  DNYPRIM*HYUMI*HZUMJ
             AK4=  DNYPRIM*HZUMI*HZUMJ
C
            ENDIF
C     
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH1
            DHILF (JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE,4)+OM*AH1      
C
         DEJACOBNS(JDOFE,IDOFE,1)=DEJACOBNS(JDOFE,IDOFE,1)+OM*AJ1
         DEJACOBNS(JDOFE,IDOFE,1)=DEJACOBNS(JDOFE,IDOFE,1)+OM*AJ2
         DEJACOBNS(JDOFE,IDOFE,1)=DEJACOBNS(JDOFE,IDOFE,1)+OM*AJ3
         DEJACOBNS(JDOFE,IDOFE,4)=DEJACOBNS(JDOFE,IDOFE,4)+OM*AJ4  
C     
         DEJACOBPW(JDOFE,IDOFE,1)=DEJACOBPW(JDOFE,IDOFE,1)+OM*AK1
         DEJACOBPW(JDOFE,IDOFE,2)=DEJACOBPW(JDOFE,IDOFE,2)+OM*AK2
         DEJACOBPW(JDOFE,IDOFE,3)=DEJACOBPW(JDOFE,IDOFE,3)+OM*AK3
         DEJACOBPW(JDOFE,IDOFE,4)=DEJACOBPW(JDOFE,IDOFE,4)+OM*AK4 
c           
235      CONTINUE
230      CONTINUE
c
       ELSE
C
C***  use symmetric part of grad
c
       DNYPRIM=2D0*DNYPRIM
C
C ***  Summing up over all pairs of multiindices
       DO 240 JDOFE=1,IDFL
       JDOFEH=KDFL(JDOFE)
       HBASJ1=DBAS(JDOFEH,1)
       HBASJ2=DBAS(JDOFEH,2)
       HBASJ3=DBAS(JDOFEH,3)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2
C
       HYUMJ=DUX1*HBASJ2+0.5D0*(DUY1+DUX2)*HBASJ3
       HZUMJ=0.5D0*(DUY1+DUX2)*HBASJ2+DUY2*HBASJ3
c
       DO 245 IDOFE=1,IDFL
       IF (IDOFE.EQ.JDOFE) THEN
             AH1= DNY*(2d0*HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
             AH2= DNY*HBASJ2*HBASJ3
             AH3= DNY*HBASJ3*HBASJ2
             AH4= DNY*(HBASJ2**2+2d0*HBASJ3**2)+CT0*HBASJ1**2
             ahilf=DNY*(HBASJ2**2+HBASJ3**2)+CT0*HBASJ1**2
c
             AJ1=  DUX1*HBASJ1**2
             AJ2=  DUY1*HBASJ1**2
             AJ3=  DUX2*HBASJ1**2
             AJ4=  DUY2*HBASJ1**2
C
             AK1=  DNYPRIM*HYUMJ**2
             AK2=  DNYPRIM*HZUMJ*HYUMJ
             AK3=  DNYPRIM*HYUMJ*HZUMJ
             AK4=  DNYPRIM*HZUMJ**2
c
       ELSE
        IDOFEH=KDFL(IDOFE)
        HBASI1=DBAS(IDOFEH,1)
        HBASI2=DBAS(IDOFEH,2)
        HBASI3=DBAS(IDOFEH,3)
        HSUMI=HBASI2*DU1+HBASI3*DU2
c
       HYUMI=DUX1*HBASI2+0.5D0*(DUY1+DUX2)*HBASI3
       HZUMI=0.5D0*(DUY1+DUX2)*HBASI2+DUY2*HBASI3
c
             AH1= DNY*(2d0*HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             AH2= DNY*HBASI2*HBASJ3!HBASI3*HBASJ2
             AH3= DNY*HBASI3*HBASJ2!HBASI2*HBASJ3
             AH4= DNY*(HBASJ2*HBASI2+2d0*HBASJ3*HBASI3)
     *         +   CT0*HBASJ1*HBASI1
             ahilf=DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *         +  CT0*HBASJ1*HBASI1
c
             AJ1=  DUX1*HBASJ1*HBASI1
             AJ2=  DUY1*HBASJ1*HBASI1
             AJ3=  DUX2*HBASJ1*HBASI1
             AJ4=  DUY2*HBASJ1*HBASI1
C
             AK1=  DNYPRIM*HYUMI*HYUMJ
             AK2=  DNYPRIM*HZUMI*HYUMJ
             AK3=  DNYPRIM*HYUMI*HZUMJ
             AK4=  DNYPRIM*HZUMI*HZUMJ
c
       ENDIF
C
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH1
            DENTRY(JDOFE,IDOFE,2)=DENTRY(JDOFE,IDOFE,2)+OM*AH2
            DENTRY(JDOFE,IDOFE,3)=DENTRY(JDOFE,IDOFE,3)+OM*AH3
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH4
            DHILF (JDOFE,IDOFE)=DHILF (JDOFE,IDOFE)+OM*AHilf    
c
         DEJACOBNS(JDOFE,IDOFE,1)=DEJACOBNS(JDOFE,IDOFE,1)+OM*AJ1
         DEJACOBNS(JDOFE,IDOFE,2)=DEJACOBNS(JDOFE,IDOFE,2)+OM*AJ2
         DEJACOBNS(JDOFE,IDOFE,3)=DEJACOBNS(JDOFE,IDOFE,3)+OM*AJ3
         DEJACOBNS(JDOFE,IDOFE,4)=DEJACOBNS(JDOFE,IDOFE,4)+OM*AJ4
C
         DEJACOBPW(JDOFE,IDOFE,1)=DEJACOBPW(JDOFE,IDOFE,1)+OM*AK1
         DEJACOBPW(JDOFE,IDOFE,2)=DEJACOBPW(JDOFE,IDOFE,2)+OM*AK2
         DEJACOBPW(JDOFE,IDOFE,3)=DEJACOBPW(JDOFE,IDOFE,3)+OM*AK3
         DEJACOBPW(JDOFE,IDOFE,4)=DEJACOBPW(JDOFE,IDOFE,4)+OM*AK4 
c        
            
 245        CONTINUE
240    CONTINUE
       ENDIF ! ISPGRAD
         ELSE!     dcmass.eq.0d0
C     
C***  Summing up over 1 pair of multiindices
          DO 250 JDOFE=1,IDFL
           HBASJ1=DBAS(KDFL(JDOFE),1)
C     
           DO 260 IDOFE=1,IDFL
            HBASI1=DBAS(KDFL(IDOFE),1)
C     
            AH=-1D0/THSTEP*HBASJ1*HBASI1
C     
            DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM*AH
            DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM*AH
 260       CONTINUE
 250      CONTINUE
C     
         ENDIF
C     
 200    CONTINUE
C
      locny(iel)=avgny/ncubp
C    
      DCONST=0D0
c
        DO 400 JDOFE=1,IDFL
         DO 400 IDOFE=1,IDFL
          DENTH1=THSTEP*DENTRY(JDOFE,IDOFE,1)*DCONST
          DENTH2=THSTEP*DENTRY(JDOFE,IDOFE,2)*DCONST
          DENTH3=THSTEP*DENTRY(JDOFE,IDOFE,3)*DCONST
          DENTH4=THSTEP*DENTRY(JDOFE,IDOFE,4)*DCONST
          DHHILF=THSTEP*DHILF (JDOFE,IDOFE)*DCONST
C
          DENTJ1=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,1)
          DENTJ2=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,2)
          DENTJ3=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,3)
          DENTJ4=DNNS*THSTEP*DEJACOBNS(JDOFE,IDOFE,4)
C
          DENTK1=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,1)
          DENTK2=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,2)
          DENTK3=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,3)
          DENTK4=DNPW*THSTEP*DEJACOBPW(JDOFE,IDOFE,4) 
c    
          BSONST=.false.
          IF ((IPRECO.EQ.0).or.(IPRECO.eq.1)) GOTO 9911
          IF (ILEV.eq.NLMAX) THEN
             IF ((MOD(IZBV2+1,IPRECO).eq.0)) BSONST=.true.
          ELSE
             IF ((MOD(IZBV2,IPRECO).eq.0)) BSONST=.true.
          ENDIF
 9911     CONTINUE
c
      IF (IDEF.LT.2) THEN
       IA   =KENTRY(JDOFE,IDOFE)
       A(IA)=A(IA)+REAL(DHHILF)
       A(3*NA+IA)=A(3*NA+IA)+REAL(DHHILF)
      ENDIF! Auf die Matrix kommt nur der Grad Teil, Rest in ADDSN2
C     
          IF (IDEF.GT.0) THEN 
           IDFG=KDFG(IDOFE)
           JDFG=KDFG(JDOFE)
           D1(JDFG)= D1(JDFG)-DENTH1*U1(IDFG)-DENTH2*U2(IDFG)
           D2(JDFG)= D2(JDFG)-DENTH3*U1(IDFG)-DENTH4*U2(IDFG)
          ENDIF !Den Defekt berechne ich aber voll hier!
C     
 400     CONTINUE
c      
C     
 100    CONTINUE
C     
C     
99999  END
