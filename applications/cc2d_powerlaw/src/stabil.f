************************************************************************
      SUBROUTINE STABIL (UTR1,VTR1,UTR2,VTR2,A1L,A2L,U,V,
     *                  D1,D2,A,NA1,KCOLA1,KLDA1,
     *                  VB1,VB2,NBB,VBN1,VBN2,KCOLBB,KLDBB,
     *                  KVERT,KMID,KADJ,KMEL,
     *                  DCORVG,ELE,COEFFN,IDEF,DCMASS)
************************************************************************
c
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'bouss.inc'
      INCLUDE 'block.inc'
c
      INCLUDE 'jump.inc'
C
      REAL *4 VB1,VB2,VBN1,VBN2
C
      DIMENSION UTR1(*),UTR2(*),VTR1(*),VTR2(*) !Transportrichtungen
      DIMENSION U(*),V(*) !Geschwindigkeiten
      DIMENSION D1(*),D2(*) ! Defekte oder rechte Seiten
      DIMENSION A(*),KCOLA1(*),KLDA1(*)
      DIMENSION VB1(*),VB2(*),VBN1(*),VBN2(*),KCOLBB(*),KLDBB(*)
      DIMENSION KADJ(NNVE,*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*),KMEL(2,*)
c
C *** definition of finite elements
      EXTERNAL E030,E031,EM31,EM30,ELE
c
      DNPW=DNPWSTART
      IF(IUSENEWT .EQ. 1) DNPW=1.0D0
      IF(IUSENEWT .EQ. 1)  DNNS=1.0D0     
c      PRINT*,'IUSENEWT,DNPW ',IUSENEWT,DNPW
c      PRINT*,'IUSENEWT,DNNS ',IUSENEWT,DNNS
      KNY=L(KLNY(NLEV))
C
c1.0D-4   DMPD    (limit for defect improvement)
c1.0D-4   DMPMG   (damping of residuals for mg-it.)
c1.0D-2   EPSMG   (limit for residuals for mg-it.)
c1.0D-2   DMPSL   (damping of residuals for solving)
c1.0D-4   EPSSL   (limit for residuals for solving)
c$$$      IF(DNNS .GE. 0.99D0)THEN
c$$$         DMPD=1.0D-4
c$$$         DMPMG=1.0D-4
c$$$         EPSMG=1.0D-2
c$$$         DMPSL=1.0D-2
c$$$         EPSSL=1.0D-4
c$$$      ENDIF   
c      print*,DMPD,DMPMG,EPSMG,DMPSL,EPSSL
c
      CALL BNBUILD(UTR1,VTR1,UTR2,VTR2,A1L,A2L,VB1,VB2,VBN1,VBN2,
     *                KCOLBB,KLDBB,NP,KVERT,KMID,KADJ,NBB,DCORVG,EM30)
c
                  CALL  STABILCNV (UTR1,VTR1,UTR2,VTR2,A1L,A2L,U,V,
     *                  D1,D2,A,NA1,KCOLA1,KLDA1,KVERT,KMID,DCORVG,
     *                  ELE,COEFFN,IDEF,DCMASS)
c
           IF(IJUMP .GE. 1)
     *             CALL JUMP (UTR1,VTR1,UTR2,VTR2,A1L,A2L,U,V,
     *                  D1,D2,A,NA1,KCOLA1,KLDA1,KVERT,KMID,KADJ,KMEL,
     *                  DCORVG,ELE,COEFFN,IDEF,DCMASS,DWORK(KNY))
C
                 CALL  STABILNEWT (UTR1,VTR1,UTR2,VTR2,A1L,A2L,U,V,
     *                  D1,D2,A,NA1,KCOLA1,KLDA1,KVERT,KMID,DCORVG,
     *                  ELE,COEFFN,IDEF,DCMASS)
C
      END
C
************************************************************************
      SUBROUTINE STABILCNV (UTR1,VTR1,UTR2,VTR2,A1L,A2L,U,V,
     *                  D1,D2,A,NA1,KCOLA1,KLDA1,KVERT,KMID,DCORVG,
     *                  ELE,COEFFN,IDEF,DCMASS)
************************************************************************
c
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'bouss.inc' 
      INCLUDE 'block.inc'
c
      real A
      DOUBLE PRECISION UTR1(*),UTR2(*),VTR1(*),VTR2(*) !Transportrichtungen
      DOUBLE PRECISION U(*),V(*) !Geschwindigkeiten
      DOUBLE PRECISION D1(*),D2(*) ! Defekte oder rechte Seiten
      DIMENSION A(*),KCOLA1(*),KLDA1(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
c
C *** definition of finite elements
      EXTERNAL E030,E031,EM31,EM30,ELE
c
c      write (*,*) 'in stabil', na,na1
        IF (IUPW.EQ.1) THEN
           IF (IPRECA.EQ.4) THEN
              IF (IELT.EQ.0) 
     *          CALL ADDSTP(UTR1,VTR1,UTR2,VTR2,A1L,A2L,
     *                     U,V,D1,D2,A,NA1,KCOLA1,KLDA1,
     *                     KVERT,KMID,DCORVG,
     *                     E031,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))
              IF (IELT.EQ.1) 
     *          CALL ADDSTP(UTR1,VTR1,UTR2,VTR2,A1L,A2L,
     *                     U,V,D1,D2,A,NA1,KCOLA1,KLDA1,
     *                     KVERT,KMID,DCORVG,
     *                     E031,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))
              IF (IELT.EQ.2) 
     *          CALL ADDSTN(UTR1,VTR1,UTR2,VTR2,A1L,A2L,
     *                     U,V,D1,D2,A,NA1,KCOLA1,KLDA1,
     *                     KVERT,KMID,DCORVG,
     *                     EM31,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))
              IF (IELT.EQ.3) THEN
                IF (IZBV1.GE.0) THEN
                   CALL ADDSTN(UTR1,VTR1,UTR2,VTR2,A1L,A2L,
     *                     U,V,D1,D2,A,NA1,KCOLA1,KLDA1,
     *                     KVERT,KMID,DCORVG,
     *                     EM30,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))
                 ELSE
                   CALL ADDSN1(UTR1,VTR1,UTR2,VTR2,A1L,A2L,
     *                     U,V,D1,D2,A,NA1,KCOLA1,KLDA1,
     *                     KVERT,KMID,DCORVG,
     *                     EM30,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))
                 ENDIF!izbv1
              ENDIF!ielt=3
           ENDIF!ipreca=4
          if (ISTOK.NE.1)
     *    CALL GUPWD (UTR1,VTR1,UTR2,VTR2,
     *               A1L,A2L,U,V,D1,D2,
     *               A,NA1,KCOLA1,KLDA1,KVERT,KMID,DCORVG,IDEF,
     *               DWORK(L(KLNY(ILEV))))
        ELSE! (IUPWD.NE.1)
         IF (IELT.EQ.0) 
     *    CALL SUPWDG(UTR1,VTR1,UTR2,VTR2,
     *               A1L,A2L,U,V,D1,D2,
     *               A,NA1,KCOLA1,KLDA1,KVERT,KMID,DCORVG,
     *                E031,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))
         IF (IELT.EQ.1) 
     *    CALL SUPWDG(UTR1,VTR1,UTR2,VTR2,
     *               A1L,A2L,U,V,D1,D2,
     *               A,NA1,KCOLA1,KLDA1,KVERT,KMID,DCORVG,
     *               E030,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))

         IF (IELT.EQ.2) 
     *    CALL SUPWNP(UTR1,VTR1,UTR2,VTR2,
     *               A1L,A2L,U,V,D1,D2,
     *               A,NA1,KCOLA1,KLDA1,KVERT,KMID,DCORVG,
     *               EM31,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))
         IF (IELT.EQ.3)   
     *    CALL SUPWNP(UTR1,VTR1,UTR2,VTR2,
     *               A1L,A2L,U,V,D1,D2,
     *               A,NA1,KCOLA1,KLDA1,KVERT,KMID,DCORVG,
     *               EM30,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))
         ENDIF
c
c
         
c
      END
c
c
************************************************************************
      SUBROUTINE STABL1 (UTR1,VTR1,UTR2,VTR2,A1L,A2L,U,V,
     *                  D1,D2,A,NA1,KCOLA1,KLDA1,KVERT,KMID,DCORVG,
     *                  ELE,COEFFN,IDEF,DCMASS)
************************************************************************
c
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'bouss.inc'
      INCLUDE 'block.inc'
c
      REAL A
      DOUBLE PRECISION UTR1(*),UTR2(*),VTR1(*),VTR2(*) !Transportrichtungen
      DOUBLE PRECISION U(*),V(*) !Geschwindigkeiten
      DOUBLE PRECISION D1(*),D2(*) ! Defekte oder rechte Seiten
      DIMENSION A(*),KCOLA1(*),KLDA1(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
c
C *** definition of finite elements
      EXTERNAL E030,E031,EM31,EM30,ELE
c
      write (*,*) 'simma irgendwann mal in stabl1?'
         IF (IELT.EQ.0) 
     *    CALL SUPWDG(UTR1,VTR1,UTR2,VTR2,
     *               A1L,A2L,U,V,D1,D2,
     *               A,NA1,KCOLA1,KLDA1,KVERT,KMID,DCORVG,
     *                E031,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))
         IF (IELT.EQ.1) 
     *    CALL SUPWDG(UTR1,VTR1,UTR2,VTR2,
     *               A1L,A2L,U,V,D1,D2,
     *               A,NA1,KCOLA1,KLDA1,KVERT,KMID,DCORVG,
     *                E030,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))

         IF (IELT.EQ.2) 
     *    CALL SUPWNP(UTR1,VTR1,UTR2,VTR2,
     *               A1L,A2L,U,V,D1,D2,
     *               A,NA1,KCOLA1,KLDA1,KVERT,KMID,DCORVG,
     *                EM31,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))
         IF (IELT.EQ.3) 
     *    CALL SUPWNP(UTR1,VTR1,UTR2,VTR2,
     *               A1L,A2L,U,V,D1,D2,
     *               A,NA1,KCOLA1,KLDA1,KVERT,KMID,DCORVG,
     *                EM30,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))
      END

************************************************************************
      SUBROUTINE STABILNEWT (UTR1,VTR1,UTR2,VTR2,A1L,A2L,U,V,
     *                  D1,D2,A,NA1,KCOLA1,KLDA1,KVERT,KMID,DCORVG,
     *                  ELE,COEFFN,IDEF,DCMASS)
************************************************************************
c
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'bouss.inc'
      INCLUDE 'block.inc'
      INCLUDE 'jump.inc'

      REAL A
      DOUBLE PRECISION UTR1(*),UTR2(*),VTR1(*),VTR2(*) !Transportrichtungen
      DOUBLE PRECISION U(*),V(*) !Geschwindigkeiten
      DOUBLE PRECISION D1(*),D2(*) ! Defekte oder rechte Seiten
      DIMENSION A(*),KCOLA1(*),KLDA1(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*)
c
C *** definition of finite elements
      EXTERNAL E030,E031,EM31,EM30,ELE
      
C
        IF (IUPW.EQ.1) THEN
           IF (IPRECA.EQ.4) THEN
              IF (IELTNEWTON.EQ.0) 
     *          CALL ADDSTPNEW(UTR1,VTR1,UTR2,VTR2,A1L,A2L,
     *                     U,V,D1,D2,A,NA1,KCOLA1,KLDA1,
     *                     KVERT,KMID,DCORVG,
     *                     E031,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))
              IF (IELTNEWTON.EQ.1) 
     *          CALL ADDSTPNEW(UTR1,VTR1,UTR2,VTR2,A1L,A2L,
     *                     U,V,D1,D2,A,NA1,KCOLA1,KLDA1,
     *                     KVERT,KMID,DCORVG,
     *                     E031,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))
              IF (IELTNEWTON.EQ.2) 
     *          CALL ADDSTNNEW(UTR1,VTR1,UTR2,VTR2,A1L,A2L,
     *                     U,V,D1,D2,A,NA1,KCOLA1,KLDA1,
     *                     KVERT,KMID,DCORVG,
     *                     EM31,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))
              IF (IELTNEWTON.EQ.3) THEN
                IF (IZBV1.GE.0) THEN
                   CALL ADDSTNNEW(UTR1,VTR1,UTR2,VTR2,A1L,A2L,
     *                     U,V,D1,D2,A,NA1,KCOLA1,KLDA1,
     *                     KVERT,KMID,DCORVG,
     *                     EM30,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))
                 ELSE
                   CALL ADDSN1NEW(UTR1,VTR1,UTR2,VTR2,A1L,A2L,
     *                     U,V,D1,D2,A,NA1,KCOLA1,KLDA1,
     *                     KVERT,KMID,DCORVG,
     *                     EM30,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))
                 ENDIF!izbv1
              ENDIF!ielt=3
           ENDIF!ipreca=4
c          if (ISTOK.NE.1)
c     *    CALL GUPWD (UTR1,VTR1,UTR2,VTR2,
c     *               A1L,A2L,U,V,D1,D2,
c     *               A,NA1,KCOLA1,KLDA1,KVERT,KMID,DCORVG,IDEF,
c     *               DWORK(L(KLNY(ILEV))))
        ELSE! (IUPWD.NE.1)
         IF (IELTNEWTON.EQ.0) 
     *     CALL SUPWDGNEW(UTR1,VTR1,UTR2,VTR2,
     *               A1L,A2L,U,V,D1,D2,
     *               A,NA1,KCOLA1,KLDA1,KVERT,KMID,DCORVG,
     *                E031,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))
         IF (IELTNEWTON.EQ.1) 
     *    CALL SUPWDGNEW(UTR1,VTR1,UTR2,VTR2,
     *               A1L,A2L,U,V,D1,D2,
     *               A,NA1,KCOLA1,KLDA1,KVERT,KMID,DCORVG,
     *                E030,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))

         IF (IELTNEWTON.EQ.2) 
     *    CALL SUPWNPNEW(UTR1,VTR1,UTR2,VTR2,
     *               A1L,A2L,U,V,D1,D2,
     *               A,NA1,KCOLA1,KLDA1,KVERT,KMID,DCORVG,
     *                EM31,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))
         IF (IELTNEWTON.EQ.3)   
     *    CALL SUPWNPNEW(UTR1,VTR1,UTR2,VTR2,
     *               A1L,A2L,U,V,D1,D2,
     *               A,NA1,KCOLA1,KLDA1,KVERT,KMID,DCORVG,
     *                EM30,COEFF,IDEF,DCMASS,DWORK(L(KLNY(ILEV))))
         ENDIF
c
c
         
c
      END
c
************************************************************************
      SUBROUTINE JUMP(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                KCOLA,KLDA,KVERT,KMID,KADJ,KMEL,DCORVG,ELE,
     *                COEFFN,IDEF,DCMASS,DNNY)
************************************************************************
*     Purpose: -  Adds the JUMP to the  matrix block 
*
*
*     
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4)
      PARAMETER (NNLEV=9)
C
      DIMENSION A(*),KCOLA(*),KLDA(*),DNNY(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KADJ(NNVE,*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*),KMEL(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
c      DIMENSION KENTRY(NNBAS,NNBAS)
c      DIMENSION DENTRY(NNBAS,NNBAS,4),DHILF(NNBAS,NNBAS)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
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
      EXTERNAL E030,E031,EM31,EM30,ELE
      
      IF (IELTJUMP .EQ. 0)
     *   CALL JUMPCP(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                KCOLA,KLDA,KVERT,KMID,KADJ,KMEL,DCORVG,E031,
     *                COEFFN,IDEF,DCMASS,DNNY)
C
      IF (IELTJUMP .EQ. 1)
     *   CALL JUMPCP(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                KCOLA,KLDA,KVERT,KMID,KADJ,KMEL,DCORVG,E030,
     *                COEFFN,IDEF,DCMASS,DNNY)
C
      IF (IELTJUMP .EQ. 2)
     *   CALL JUMPC(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                KCOLA,KLDA,KVERT,KMID,KADJ,KMEL,DCORVG,EM31,
     *                COEFFN,IDEF,DCMASS,DNNY)
C
      IF (IELTJUMP .EQ. 3)
     *   CALL JUMPC(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                KCOLA,KLDA,KVERT,KMID,KADJ,KMEL,DCORVG,EM30,
     *                COEFFN,IDEF,DCMASS,DNNY)
C
      END
c     
c
************************************************************************
      SUBROUTINE JUMPL(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                KCOLA,KLDA,KVERT,KMID,KADJ,KMEL,DCORVG,ELE,
     *                COEFFN,IDEF,DCMASS,DNNY)
************************************************************************
*     Purpose: -  Adds the JUMP to the  matrix block 
*
*
*     
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4)
      PARAMETER (NNLEV=9)
C
      DIMENSION A(*),KCOLA(*),KLDA(*),DNNY(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*),KMEL(2,*)
      DIMENSION KADJ(NNVE,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
c      DIMENSION KENTRY(NNBAS,NNBAS)
c      DIMENSION DENTRY(NNBAS,NNBAS,4),DHILF(NNBAS,NNBAS)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
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
      EXTERNAL DVISCO,E030,E031,EM31,EM30,ELE
      SAVE
c      EXTERNAL E030,E031,EM31,EM30 
c$$$      IF (IELTJUMP .EQ. 0)
c$$$     *   CALL JUMPLCP(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
c$$$     *                KCOLA,KLDA,KVERT,KMID,KADJ,KMEL,DCORVG,E031,
c$$$     *                COEFFN,IDEF,DCMASS,DNNY)
c$$$C
c$$$      IF (IELTJUMP .EQ. 1)
c$$$     *   CALL JUMPLCP(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
c$$$     *                KCOLA,KLDA,KVERT,KMID,KADJ,KMEL,DCORVG,E030,
c$$$     *                COEFFN,IDEF,DCMASS,DNNY)
c$$$C
c$$$      IF (IELTJUMP .EQ. 2)
c$$$     *   CALL JUMPLC(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
c$$$     *                KCOLA,KLDA,KVERT,KMID,KADJ,KMEL,DCORVG,EM31,
c$$$     *                COEFFN,IDEF,DCMASS,DNNY)
c$$$C
c$$$      IF (IELTJUMP .EQ. 3)
c$$$     *   CALL JUMPLC(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
c$$$     *                KCOLA,KLDA,KVERT,KMID,KADJ,KMEL,DCORVG,EM30,
c$$$     *                COEFFN,IDEF,DCMASS,DNNY)
c$$$C
      END
c
C
************************************************************************
      SUBROUTINE JUMPC(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                KCOLA,KLDA,KVERT,KMID,KADJ,KMEL,DCORVG,ELE,
     *                COEFFN,IDEF,DCMASS,DNNY)
************************************************************************
*     Purpose: -  Adds the JUMP to the  matrix block 
*
*
*     
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4)
      PARAMETER (NNLEV=9)
C
      DIMENSION A(*),KCOLA(*),KLDA(*),DNNY(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*),KMEL(2,*)
      DIMENSION KADJ(NNVE,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION DENTRY(NNBAS,NNBAS,4)
      DIMENSION DXG(3),DYG(3),DALPHA(3),ILEADJ(2),INUG(8),INUG7(7)
      DIMENSION DBASSF(9,7),DBASEJUMP(9,8)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
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
      EXTERNAL DVISCO,ELE
      SAVE 
C
      ICUB=3
      IDFL7=7
C
      CALL ZTIME(T1)
C
      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
       CT0=DCMASS/THSTEP
      ELSE
       CT0=0D0
      ENDIF
C
      DNY=NY
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
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************

      CALL ELE(0D0,0D0,-2)
c
C
C *** Loop over all midpoints 
C
      DO 100  IMT=1,NMT
C
         IF (A2L.EQ.0D0) THEN
            DU1=U1L1(IMT)
            DU2=U1L2(IMT)        
         ELSE       
            DU1=A1L*U1L1(IMT)+A2L*U2L1(IMT)
            DU2=A1L*U1L2(IMT)+A2L*U2L2(IMT)
         ENDIF   
C
C *** deter. of the contributed element and the 7x7 numer.
      CALL ADJELEM(KVERT,KMID,KADJ,KMEL,DCORVG,IMT,NVT,NMT,
     *             	ILEADJ,IVT1,IVT2,INUG,INUG7)
c 
      
c      DNNY1=DNNY(ILEADJ(1))
c      IF(ILEADJ(2) .NE. 0) THEN
c         DNNY2=DNNY(ILEADJ(2))
c         DNNY3=0.5D0*(DNNY1+DNNY2)
c        ELSE
c            DNNY3=DNNY1
c      ENDIF  
c      DNY=NY*DNNY3 
c
C *** Determine entry positions in matrix
      DO 110 JDOFE=1, IDFL7
         DENTRY(JDOFE,JDOFE,1)=0D0
         DENTRY(JDOFE,JDOFE,2)=0D0
         DENTRY(JDOFE,JDOFE,3)=0D0
         DENTRY(JDOFE,JDOFE,4)=0D0
      DO 111 IDOFE=1, IDFL7
          DENTRY(JDOFE,IDOFE,1)=0D0
          DENTRY(JDOFE,IDOFE,2)=0D0
          DENTRY(JDOFE,IDOFE,3)=0D0
          DENTRY(JDOFE,IDOFE,4)=0D0
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices on the side IMT
C
      PX1=DCORVG(1,IVT1)
      PY1=DCORVG(2,IVT1)	
      PX2=DCORVG(1,IVT2)
      PY2=DCORVG(2,IVT2)
C
      HLOCAL=SQRT((PX2-PX1)**2+(PY2-PY1)**2)	
C
      DT1=(PX2-PX1)/HLOCAL
      DT2=(PY2-PY1)/HLOCAL
C
      DUDOTN=ABS(-DU1*DT2+DU2*DT1)
C
C
C*** CALCUL OF THE THREE GAUSS POINTS OF NUMERIQUE INTEGRATION 
C***               ON THE SIDE INT 
C                 
       DCONST1= 0D0
       DCONST2= SQRT(3D0/5D0)
       DCONST3= -DCONST2
C
       DXG(1)=0.5D0*(DCONST1*(PX2-PX1)+(PX1+PX2))
       DXG(2)=0.5D0*(DCONST2*(PX2-PX1)+(PX1+PX2))
       DXG(3)=0.5D0*(DCONST3*(PX2-PX1)+(PX1+PX2))
       DYG(1)=0.5D0*(DCONST1*(PY2-PY1)+(PY1+PY2))
       DYG(2)=0.5D0*(DCONST2*(PY2-PY1)+(PY1+PY2)) 
       DYG(3)=0.5D0*(DCONST3*(PY2-PY1)+(PY1+PY2))
C  
      DALPHA(1)=8D0/9D0
      DALPHA(2)=5D0/9D0
      DALPHA(3)=DALPHA(2)
C
C *** Loop over all cubature points
C
      DO 200 ICUBP = 1, ICUB
      XI1=DXG(ICUBP)
      XI2=DYG(ICUBP)
c	
         OM1 =DALPHA(ICUBP)
         OM2 =DALPHA(ICUBP)*HLOCAL**2	
         OM3 =DALPHA(ICUBP)*HLOCAL*HLOCAL**2	
C
C
      CALL BASEJUMP(KVERT,KMID,KADJ,KMEL,DCORVG,IMT,NVT,NMT,
     *             ILEADJ,DT1,DT2,INUG,ELE,XI1,XI2,DBASEJUMP)
c
C *** CALCUL OF THE APPROXIMATE VELOCITY 
      DU1=0.0D0
      DU2=0.0D0
      DO I=1,4
         HBAS=DBASEJUMP(1,I)
         J=INUG(I)
         IF (A2L.EQ.0D0) THEN
            DU1=DU1+U1L1(J)*HBAS
            DU2=DU1+U1L2(J)*HBAS          
         ELSE       
            DU1=DU1+A1L*U1L1(J)+A2L*U2L1(J)
            DU2=DU1+A1L*U1L2(J)+A2L*U2L2(J)  
         ENDIF  
       ENDDO  
C
      CALL BASESF(KVERT,KMID,KADJ,KMEL,DCORVG,IMT,NVT,NMT,
     *             ILEADJ,INUG,XI1,XI2,DBASEJUMP,DBASSF)
c
      IF (DCMASS.NE.0D0) THEN
C

        IF (ISPGRAD .EQ. 1 .OR. ISPGRAD .EQ. 0)THEN
C   
C ***  Summing up over all pairs of multiindices
       DO 240 JDOFE=1, IDFL7
       HBASJ1=DBASSF(1,JDOFE)
       HBASJ2=DBASSF(2,JDOFE)
       HBASJ3=DBASSF(3,JDOFE)
       HBASJN=DBASSF(4,JDOFE)
       HBASJT=DBASSF(5,JDOFE)
       HBASJTN1=DBASSF(6,JDOFE)
       HBASJTN2=DBASSF(7,JDOFE)
       HBASJNT1=DBASSF(8,JDOFE)
       HBASJNT2=DBASSF(9,JDOFE)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2
C       
       DO 245 IDOFE=1, IDFL7
       IF (IDOFE.EQ.JDOFE) THEN
             IF(IJUMP .EQ. 1)AH1= DNY*HBASJ1**2
             IF(IJUMP .EQ. 2)AH1= (HBASJ2**2+HBASJ3**2)
             IF(IJUMP .EQ. 3)AH1=  HBASJN**2
             IF(IJUMP .EQ. 4)AH1=  HBASJT**2
             IF(IJUMP .EQ. 5)THEN
                AH1=HBASJTN1**2
                AH2=0.0D0
                AH3=0.0D0
                AH4=HBASJTN2**2
             ENDIF   
             IF(IJUMP .EQ. 6)THEN
                AH1=HBASJNT1**2
                AH2=0.0D0
                AH3=0.0D0
                AH4=HBASJNT2**2
             ENDIF 
              IF(IJUMP .EQ. 7)THEN
                AAH1=HSUMJ**2
                ABH1=DNY*HBASJT**2
                ACH1=DUDOTN*HBASJT**2
                ADH1=HBASJTN1**2
                ADH2=0.0D0
                ADH3=0.0D0
                ADH4=HBASJTN2**2
             ENDIF  
C
       ELSE      
       HBASI1=DBASSF(1,IDOFE)
       HBASI2=DBASSF(2,IDOFE)
       HBASI3=DBASSF(3,IDOFE)
       HBASIN=DBASSF(4,IDOFE)
       HBASIT=DBASSF(5,IDOFE)
       HBASITN1=DBASSF(6,IDOFE)
       HBASITN2=DBASSF(7,IDOFE)
       HBASINT1=DBASSF(8,IDOFE)
       HBASINT2=DBASSF(9,IDOFE)
       HSUMI=HBASI2*DU1+HBASI3*DU2       
c
             IF(IJUMP .EQ. 1)AH1= DNY*HBASJ1*HBASI1
             IF(IJUMP .EQ. 2)AH1= (HBASJ2*HBASI2+HBASJ3*HBASI3)
             IF(IJUMP .EQ. 3)AH1= HBASIN*HBASJN
             IF(IJUMP .EQ. 4)AH1= HBASIT*HBASJT
             IF(IJUMP .EQ. 5)THEN
                AH1=HBASITN1*HBASJTN1
                AH2=0.0D0
                AH3=0.0D0
                AH4=HBASITN2*HBASJTN2
             ENDIF
             IF(IJUMP .EQ. 6)THEN
                AH1=HBASINT1*HBASJNT1
                AH2=0.0D0
                AH3=0.0D0
                AH4=HBASINT2*HBASJNT2
             ENDIF
              IF(IJUMP .EQ. 7)THEN
                AAH1=HSUMJ*HSUMI
                ABH1=DNY*HBASJT*HBASIT
                ACH1=DUDOTN*HBASJT*HBASIT
                ADH1=HBASJTN1*HBASITN1
                ADH2=0.0D0
                ADH3=0.0D0
                ADH4=HBASJTN2*HBASITN2
             ENDIF  
C
      ENDIF        
C    
       IF(IJUMP .EQ. 1)THEN
               DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM1*AH1
               DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM1*AH1
       ENDIF   
       IF(IJUMP .EQ. 2)THEN
               DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM3*AH1
               DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM3*AH1
       ENDIF     
       IF(IJUMP .EQ. 3)THEN
               DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM3*AH1
               DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM3*AH1
       ENDIF      
       IF(IJUMP .EQ. 4)THEN
               DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM3*AH1
               DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM3*AH1
       ENDIF   
        IF(IJUMP .EQ. 5 .OR. IJUMP .EQ. 6)THEN
           DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM3*AH1
           DENTRY(JDOFE,IDOFE,2)=DENTRY(JDOFE,IDOFE,2)+OM3*AH2
           DENTRY(JDOFE,IDOFE,3)=DENTRY(JDOFE,IDOFE,3)+OM3*AH3
           DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM3*AH4
        ENDIF  
        IF(IJUMP .EQ. 7)THEN
           AH1=AAH1*OM3+ABH1*OM2+ACH1*OM3+ADH1*OM2
           AH2=0.0D0
           AH3=0.0D0
           AH4=AAH1*OM3+ABH1*OM2+ACH1*OM3+ADH4*OM2
C
           DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+AH1
           DENTRY(JDOFE,IDOFE,2)=DENTRY(JDOFE,IDOFE,2)+AH2
           DENTRY(JDOFE,IDOFE,3)=DENTRY(JDOFE,IDOFE,3)+AH3
           DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+AH4
        ENDIF  
C            
 245        CONTINUE
240    CONTINUE
       ENDIF 
      ELSE
C
C ***  Summing up over 1 pair of multiindices
       DO 250 JDOFE=1, IDFL7
       HBASJ1=DBASSF(1,JDOFE)
C
       DO 260 IDOFE=1, IDFL7
       HBASI1=DBASSF(1,IDOFE)
C
       AH=-1D0/THSTEP*HBASJ1*HBASI1
C
c$$$            DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
260    CONTINUE
250    CONTINUE
C
      ENDIF
C
200   CONTINUE
C
      DO 400 JDOFE=1, IDFL7
      DO 400 IDOFE=1, IDFL7
          DENTH1=THSTEP*DENTRY(JDOFE,IDOFE,1)*DJUMP
          DENTH2=THSTEP*DENTRY(JDOFE,IDOFE,2)*DJUMP
          DENTH3=THSTEP*DENTRY(JDOFE,IDOFE,3)*DJUMP
          DENTH4=THSTEP*DENTRY(JDOFE,IDOFE,4)*DJUMP
C
      IF (IDEF.LT.2 .AND. ILEADJ(2) .NE.0) THEN
       IDFG=INUG7(IDOFE)
       JDFG=INUG7(JDOFE)
           DO KCOL=KLDA(JDFG),KLDA(JDFG+1)-1
             IF(KCOLA(KCOL) .EQ. IDFG)THEN
                 A(0*NA+KCOL)=A(0*NA+KCOL)+REAL(DENTH1)
                 A(1*NA+KCOL)=A(1*NA+KCOL)+REAL(DENTH2)
                 A(2*NA+KCOL)=A(2*NA+KCOL)+REAL(DENTH3)
                 A(3*NA+KCOL)=A(3*NA+KCOL)+REAL(DENTH4)
             ENDIF
           ENDDO 
      ENDIF
C
      IF (IDEF.GT.0 .AND. ILEADJ(2) .NE.0) THEN 
       IDFG=INUG7(IDOFE)
       JDFG=INUG7(JDOFE)
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
************************************************************************
      SUBROUTINE ADJELEM(KVERT,KMID,KADJ,KMEL,DCORVG,IMT,NVT,NMT,
     *             	ILEADJ,IVT1,IVT2,INUG,INUG7)
************************************************************************
*     determination of elements adjacent to the curent midpoint
*     and the thre others contributed in numerical integration   
*     on the curenet side IMT
*     output: 1. Array ILEADJ 
*             2. IVT1, IVT2 : vertices containing curent midpoint
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4,NNBAS=21,NNDER=6)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION ILEADJ(2),KMEL(2,*),INUG(8),INUG7(7)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /COAUX1/ KDFG,KDFL,IDFL

      SAVE
C
C
C ***    
C
        DO I=1,2
	   ILEADJ(I)=0
	ENDDO 
C
C ***   DETRMINATION OF THE TWO ADJ ELEMENT TO THE MID IMT   
         IEL1=KMEL(1,IMT) 
         IEL2=KMEL(2,IMT) 
         ILEADJ(1)=KMEL(1,IMT) 
         ILEADJ(2)=KMEL(2,IMT) 
         IEL=IEL1

C ***   Determine of the other adjacent elenents to IEL   
C
C ********* Determine the coresponding vertices         
         DO 10 IVE =1,NNVE
           IVEN=IVE+1
           IF (IVEN.EQ.5) IVEN=1
         	  IMID1=KMID(IVE,IEL)-NVT
          	  IF(IMID1 .EQ. IMT)THEN
         		IVT1=KVERT(IVE ,IEL)
               		IVT2=KVERT(IVEN,IEL)
	      	 ENDIF
10       CONTINUE
C	
        DO IN=1,8
           INUG(IN)=0
        ENDDO
C
         DO  IVE=1,NNVE
             IVEN=IVE+1
             IF(IVEN .EQ. 5)IVEN=1
             IMID=KMID(IVE,ILEADJ(1))-NVT
             INUG(IVE)=IMID
         ENDDO
C
         IF(ILEADJ(2) .NE.0)THEN
            DO IVE=1,NNVE
               IVEN=IVE+1
               IF(IVEN .EQ. 5)IVEN=1
               IMID=KMID(IVE,ILEADJ(2))-NVT
               INUG(NNVE+IVE)=IMID
            ENDDO   
         ENDIF
C
         DO IN=1,7
            INUG7(IN)=0
         ENDDO
         DO IVE=1,NNVE
             IVEN=IVE+1
             IF(IVEN .EQ. 5)IVEN=1
             IMID=KMID(IVE,ILEADJ(1))-NVT
             INUG7(IVE)=IMID
         ENDDO
         IF(ILEADJ(2) .NE.0)THEN
         ICONT=0
         DO IVE=1,NNVE
            IVEN=IVE+1
             IF(IVEN .EQ. 5)IVEN=1
             IMID=KMID(IVE,ILEADJ(2))-NVT
             IF(IMID .NE. IMT)THEN
                ICONT=ICONT+1
                INUG7(NNVE+ICONT)=IMID
             ENDIF
         ENDDO
         ENDIF
c         PRINT*,'INUG7',(INUG7(I),I=1,7)
C
        END
C
************************************************************************
      SUBROUTINE BASEJUMP(KVERT,KMID,KADJ,KMEL,DCORVG,IMT,NVT,NMT,
     *             	ILEADJ,DT1,DT2,INUG,ELE,XI1,XI2,DBASEJUMP)
************************************************************************
*      
*     
*        
*     
*    
*             
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4,NNBAS=21,NNDER=6)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*),DCORVG(2,*)
      DIMENSION ILEADJ(2),KMEL(2,*),DBASEJUMP(9,8),INUG(8)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      SAVE

C
C ***    
C 
      IELTYP=-1
      DO I=1,8
         DBASEJUMP(1,I)=0D0
         DBASEJUMP(2,I)=0D0
         DBASEJUMP(3,I)=0D0
         DBASEJUMP(4,I)=0D0
         DBASEJUMP(5,I)=0D0
         DBASEJUMP(6,I)=0D0
         DBASEJUMP(7,I)=0D0
         DBASEJUMP(8,I)=0D0
         DBASEJUMP(9,I)=0D0
      ENDDO
C
      CALL ELE(0D0,0D0,-2)
C
      IEL=ILEADJ(1)
      DN1=-DT2
      DN2=DT1
C
C *** Evaluation of coordinates of the vertices
         DO  IVE = 1, NNVE
            JP=KVERT(IVE,IEL)
            KVE(IVE)=JP
            DX(IVE)=DCORVG(1,JP)
            DY(IVE)=DCORVG(2,JP)
          ENDDO 
      CALL ELE(0D0,0D0,-2)
      CALL ELE(XI1,XI2,-3)
         DO  IVE=1,NNVE
            DBASEJUMP(1,IVE)=DBAS(IVE,1)
            DBASEJUMP(2,IVE)=DBAS(IVE,2)
            DBASEJUMP(3,IVE)=DBAS(IVE,3)
            DBASEJUMP(4,IVE)= DN1*DBAS(IVE,2)+DN2*DBAS(IVE,3)
            DBASEJUMP(5,IVE)= DT1*DBAS(IVE,2)+DT2*DBAS(IVE,3)
            DBASEJUMP(6,IVE)= DBASEJUMP(5,IVE)*DN1
            DBASEJUMP(7,IVE)= DBASEJUMP(5,IVE)*DN2
            DBASEJUMP(8,IVE)= DBASEJUMP(4,IVE)*DT1
            DBASEJUMP(9,IVE)= DBASEJUMP(4,IVE)*DT2
         ENDDO 
C
         CALL ELE(0D0,0D0,-2)
         IF(ILEADJ(2) .NE.0)THEN
         IEL=ILEADJ(2) 
C
C *** Evaluation of coordinates of the vertices
         DO  IVE = 1, NNVE
            JP=KVERT(IVE,IEL)
            KVE(IVE)=JP
            DX(IVE)=DCORVG(1,JP)
            DY(IVE)=DCORVG(2,JP)
          ENDDO 
            CALL ELE(0d0,0d0,-2)
            CALL ELE(XI1,XI2,-3) 
            DO IVE=1,NNVE
            DBASEJUMP(1,NNVE+IVE)=DBAS(IVE,1)
            DBASEJUMP(2,NNVE+IVE)=DBAS(IVE,2)
            DBASEJUMP(3,NNVE+IVE)=DBAS(IVE,3)
            DBASEJUMP(4,NNVE+IVE)= DN1*DBAS(IVE,2)+DN2*DBAS(IVE,3)
            DBASEJUMP(5,NNVE+IVE)= DT1*DBAS(IVE,2)+DT2*DBAS(IVE,3)
            DBASEJUMP(6,NNVE+IVE)= DBASEJUMP(5,NNVE+IVE)*DN1
            DBASEJUMP(7,NNVE+IVE)= DBASEJUMP(5,NNVE+IVE)*DN2
            DBASEJUMP(8,NNVE+IVE)= DBASEJUMP(4,NNVE+IVE)*DT1
            DBASEJUMP(9,NNVE+IVE)= DBASEJUMP(4,NNVE+IVE)*DT2
            ENDDO 
            ENDIF      

c          PRINT*,'DBASJUMP',(DBASEJUMP(I),I=1,8)
C
99999   END
C
************************************************************************
      SUBROUTINE BASESF(KVERT,KMID,KADJ,KMEL,DCORVG,IMT,NVT,NMT,
     *             ILEADJ,INUG,XI1,XI2,DBASEJUMP,DBASSF)
************************************************************************
*      
*     
*        
*     
*    
*             
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4,NNBAS=21,NNDER=6)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*),DCORVG(2,*)
      DIMENSION ILEADJ(2),KMEL(2,*),DBASEJUMP(9,8),INUG(8),DBASSF(9,7)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      SAVE
C
C ***    
C 
      IELTYP=-1
      DO I=1,7
         DBASSF(1,I)=0D0
         DBASSF(2,I)=0D0
         DBASSF(3,I)=0D0
         DBASSF(4,I)=0D0
         DBASSF(5,I)=0D0
         DBASSF(6,I)=0D0
         DBASSF(7,I)=0D0
         DBASSF(8,I)=0D0
         DBASSF(9,I)=0D0
      ENDDO

           DO I=1,4
              IMID=INUG(I)
              IF(IMID .EQ. IMT)IND=I
              DBASSF(1,I)=DBASEJUMP(1,I)
              DBASSF(2,I)=DBASEJUMP(2,I)
              DBASSF(3,I)=DBASEJUMP(3,I)
              DBASSF(4,I)=DBASEJUMP(4,I)
              DBASSF(5,I)=DBASEJUMP(5,I)
              DBASSF(6,I)=DBASEJUMP(6,I)
              DBASSF(7,I)=DBASEJUMP(7,I)
              DBASSF(8,I)=DBASEJUMP(8,I)
              DBASSF(9,I)=DBASEJUMP(9,I)
           ENDDO
C
        ICONT=0
        DO I=5,8
           IMID=INUG(I)
           IF(IMID .EQ. IMT)THEN
              DBASSF(1,IND)=DBASSF(1,IND)-DBASEJUMP(1,I)
              DBASSF(2,IND)=DBASSF(2,IND)-DBASEJUMP(2,I)
              DBASSF(3,IND)=DBASSF(3,IND)-DBASEJUMP(3,I)
              DBASSF(4,IND)=DBASSF(4,IND)-DBASEJUMP(4,I)
              DBASSF(5,IND)=DBASSF(5,IND)-DBASEJUMP(5,I)
              DBASSF(6,IND)=DBASSF(6,IND)-DBASEJUMP(6,I)
              DBASSF(7,IND)=DBASSF(7,IND)-DBASEJUMP(7,I)
              DBASSF(8,IND)=DBASSF(8,IND)-DBASEJUMP(8,I)
              DBASSF(9,IND)=DBASSF(9,IND)-DBASEJUMP(9,I)
              GOTO 1
           ENDIF 
           ICONT=ICONT+1 
           DBASSF(1,NNVE+ICONT)=-DBASEJUMP(1,I)
           DBASSF(2,NNVE+ICONT)=-DBASEJUMP(2,I)
           DBASSF(3,NNVE+ICONT)=-DBASEJUMP(3,I)
           DBASSF(4,NNVE+ICONT)=-DBASEJUMP(4,I)
           DBASSF(5,NNVE+ICONT)=-DBASEJUMP(5,I)
           DBASSF(6,NNVE+ICONT)=-DBASEJUMP(6,I)
           DBASSF(7,NNVE+ICONT)=-DBASEJUMP(7,I)
           DBASSF(8,NNVE+ICONT)=-DBASEJUMP(8,I)
           DBASSF(9,NNVE+ICONT)=-DBASEJUMP(9,I)
 1         CONTINUE
        ENDDO   
c        PRINT*,'DBASSF(I)',(DBASSF(I),I=1,7)
c        READ(*,*)
           
C
99999   END

************************************************************************
      SUBROUTINE JUMPLC(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                KCOLA,KLDA,KVERT,KMID,KADJ,KMEL,DCORVG,ELE,
     *                COEFFN,IDEF,DCMASS,DNNY)
************************************************************************
*     Purpose: -  Adds the JUMP to the  matrix block 
*
*
*     
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4)
      PARAMETER (NNLEV=9)
C
      DIMENSION A(*),KCOLA(*),KLDA(*),DNNY(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*),KMEL(2,*)
      DIMENSION KADJ(NNVE,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION DENTRY(NNBAS,NNBAS,4)
      DIMENSION DXG(3),DYG(3),DALPHA(3),ILEADJ(2),INUG(8),INUG7(7)
      DIMENSION DBASSF(9,7),DBASEJUMP(9,8)	
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
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
      EXTERNAL DVISCO,ELE
      SAVE 
C
      ICUB=3
      IDFL7=7
C
      CALL ZTIME(T1)
C
      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
       CT0=DCMASS/THSTEP
      ELSE
       CT0=0D0
      ENDIF
C
      DNY=NY
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
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************

      CALL ELE(0D0,0D0,-2)

C
C *** Loop over all midpoints 
C
      DO 100  IMT=1,NMT
C
         IF (A2L.EQ.0D0) THEN
            DU1=U1L1(IMT)
            DU2=U1L2(IMT)        
         ELSE       
            DU1=A1L*U1L1(IMT)+A2L*U2L1(IMT)
            DU2=A1L*U1L2(IMT)+A2L*U2L2(IMT)
         ENDIF 
C
C *** deter. of the contributed element and the 7x7 numer.
      CALL ADJELEM(KVERT,KMID,KADJ,KMEL,DCORVG,IMT,NVT,NMT,
     *             	ILEADJ,IVT1,IVT2,INUG,INUG7)
C
C
c      DNNY1=DNNY(ILEADJ(1))
c      IF(ILEADJ(2) .NE. 0) THEN
c         DNNY2=DNNY(ILEADJ(2))
c         DNNY3=0.5D0*(DNNY1+DNNY2)
c         ELSE
c             DNNY3=DNNY1
c      ENDIF  
c      DNY=NY*DNNY3
c
C 	
C *** Determine entry positions in matrix
      DO 110 JDOFE=1, IDFL7
         DENTRY(JDOFE,JDOFE,1)=0D0
         DENTRY(JDOFE,JDOFE,2)=0D0
         DENTRY(JDOFE,JDOFE,3)=0D0
         DENTRY(JDOFE,JDOFE,4)=0D0
      DO 111 IDOFE=1, IDFL7
          DENTRY(JDOFE,IDOFE,1)=0D0
          DENTRY(JDOFE,IDOFE,2)=0D0
          DENTRY(JDOFE,IDOFE,3)=0D0
          DENTRY(JDOFE,IDOFE,4)=0D0
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices on the side IMT
C
      PX1=DCORVG(1,IVT1)
      PY1=DCORVG(2,IVT1)	
      PX2=DCORVG(1,IVT2)
      PY2=DCORVG(2,IVT2)
C
      HLOCAL=SQRT((PX2-PX1)**2+(PY2-PY1)**2)	
C
      DT1=(PX2-PX1)/HLOCAL
      DT2=(PY2-PY1)/HLOCAL
C
      DUDOTN=ABS(-DU1*DT2+DU2*DT1)
C
C
C*** CALCUL OF THE THREE GAUSS POINTS OF NUMERIQUE INTEGRATION 
C***               ON THE SIDE INT 
C                 
       DCONST1= 0D0
       DCONST2= SQRT(3D0/5D0)
       DCONST3= -DCONST2
C
       DXG(1)=0.5D0*(DCONST1*(PX2-PX1)+(PX1+PX2))
       DXG(2)=0.5D0*(DCONST2*(PX2-PX1)+(PX1+PX2))
       DXG(3)=0.5D0*(DCONST3*(PX2-PX1)+(PX1+PX2))
       DYG(1)=0.5D0*(DCONST1*(PY2-PY1)+(PY1+PY2))
       DYG(2)=0.5D0*(DCONST2*(PY2-PY1)+(PY1+PY2)) 
       DYG(3)=0.5D0*(DCONST3*(PY2-PY1)+(PY1+PY2))
C  
      DALPHA(1)=8D0/9D0
      DALPHA(2)=5D0/9D0
      DALPHA(3)=DALPHA(2)
C
C *** Loop over all cubature points
C
      DO 200 ICUBP = 1, ICUB
      XI1=DXG(ICUBP)
      XI2=DYG(ICUBP)	
c
         OM1 =DALPHA(ICUBP)	
         OM2 =DALPHA(ICUBP)*HLOCAL**2	
         OM3 =DALPHA(ICUBP)*HLOCAL*HLOCAL**2	
C
      CALL BASEJUMP(KVERT,KMID,KADJ,KMEL,DCORVG,IMT,NVT,NMT,
     *             ILEADJ,DT1,DT2,INUG,ELE,XI1,XI2,DBASEJUMP)
C
C *** CALCUL OF THE APPROXIMATE VELOCITY 
      DU1=0.0D0
      DU2=0.0D0
      DO I=1,4
         HBAS=DBASEJUMP(1,I)
         J=INUG(I)
         IF (A2L.EQ.0D0) THEN
            DU1=DU1+U1L1(J)*HBAS
            DU2=DU1+U1L2(J)*HBAS          
         ELSE       
            DU1=DU1+A1L*U1L1(J)+A2L*U2L1(J)
            DU2=DU1+A1L*U1L2(J)+A2L*U2L2(J)  
         ENDIF  
       ENDDO 

         
      CALL BASESF(KVERT,KMID,KADJ,KMEL,DCORVG,IMT,NVT,NMT,
     *             ILEADJ,INUG,XI1,XI2,DBASEJUMP,DBASSF)
c
      IF (DCMASS.NE.0D0) THEN
C
        IF (ISPGRAD .EQ. 1 .OR. ISPGRAD .EQ. 0)THEN
C
C***   use symmetric part of grad
C   
C ***  Summing up over all pairs of multiindices
       DO 240 JDOFE=1, IDFL7
       HBASJ1 =DBASSF(1,JDOFE)
       HBASJ2 =DBASSF(2,JDOFE)
       HBASJ3 =DBASSF(3,JDOFE)
       HBASJN =DBASSF(4,JDOFE)
       HBASJT =DBASSF(5,JDOFE)
       HBASJTN1=DBASSF(6,JDOFE)
       HBASJTN2=DBASSF(7,JDOFE)
       HBASJNT1=DBASSF(8,JDOFE)
       HBASJNT2=DBASSF(9,JDOFE)
       HSUMJ=HBASJ2*DU1+HBASJ3*DU2
C
       DO 245 IDOFE=1, IDFL7
       IF (IDOFE.EQ.JDOFE) THEN
             IF(IJUMP .EQ. 1)AH1= DNY*HBASJ1**2
             IF(IJUMP .EQ. 2)AH1= (HBASJ2**2+HBASJ3**2)
             IF(IJUMP .EQ. 3)AH1= HBASJN**2
             IF(IJUMP .EQ. 4)AH1= HBASJT**2
             IF(IJUMP .EQ. 5)THEN
                AH1=HBASJTN1**2
                AH2=0.0D0
                AH3=0.0D0
                AH4=HBASJTN2**2
             ENDIF 
             IF(IJUMP .EQ. 6)THEN
                AH1=HBASJNT1**2
                AH2=0.0D0
                AH3=0.0D0
                AH4=HBASJNT2**2
             ENDIF 
              IF(IJUMP .EQ. 7)THEN
                AAH1=HSUMJ**2
                ABH1=DNY*HBASJT**2
                ACH1=DUDOTN*HBASJT**2
                ADH1=HBASJTN1**2
                ADH2=0.0D0
                ADH3=0.0D0
                ADH4=HBASJTN2**2
             ENDIF  
C 
       ELSE
       HBASI1 =DBASSF(1,IDOFE)
       HBASI2 =DBASSF(2,IDOFE)
       HBASI3 =DBASSF(3,IDOFE)
       HBASIN =DBASSF(4,IDOFE)
       HBASIT =DBASSF(5,IDOFE)
       HBASITN1=DBASSF(6,IDOFE)
       HBASITN2=DBASSF(7,IDOFE)
       HBASINT1=DBASSF(8,IDOFE)
       HBASINT2=DBASSF(9,IDOFE)
       HSUMI=HBASI2*DU1+HBASI3*DU2
c
             IF(IJUMP .EQ. 1)AH1= DNY*HBASJ1*HBASI1
             IF(IJUMP .EQ. 2)AH1= (HBASJ2*HBASI2+HBASJ3*HBASI3)
             IF(IJUMP .EQ. 3)AH1= HBASIN*HBASJN
             IF(IJUMP .EQ. 4)AH1= HBASIT*HBASJT
             IF(IJUMP .EQ. 5)THEN
                AH1=HBASITN1*HBASJTN1
                AH2=0.0D0
                AH3=0.0D0
                AH4=HBASITN2*HBASJTN2
             ENDIF
             IF(IJUMP .EQ. 6)THEN
                AH1=HBASINT1*HBASJNT1
                AH2=0.0D0
                AH3=0.0D0
                AH4=HBASINT2*HBASJNT2
             ENDIF
              IF(IJUMP .EQ. 7)THEN
                AAH1=HSUMJ*HSUMI
                ABH1=DNY*HBASJT*HBASIT
                ACH1=DUDOTN*HBASJT*HBASIT
                ADH1=HBASJTN1*HBASITN1
                ADH2=0.0D0
                ADH3=0.0D0
                ADH4=HBASJTN2*HBASITN2
             ENDIF  
C
       ENDIF
C
       IF(IJUMP .EQ. 1)THEN
               DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM1*AH1
               DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM1*AH1
       ENDIF   
       IF(IJUMP .EQ. 2)THEN
               DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM3*AH1
               DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM3*AH1
       ENDIF     
       IF(IJUMP .EQ. 3)THEN
               DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM3*AH1
               DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM3*AH1
       ENDIF      
       IF(IJUMP .EQ. 4)THEN
               DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM3*AH1
               DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM3*AH1
       ENDIF   
        IF(IJUMP .EQ. 5 .OR. IJUMP .EQ. 6)THEN
                DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+OM3*AH1
                DENTRY(JDOFE,IDOFE,2)=DENTRY(JDOFE,IDOFE,2)+OM3*AH2
                DENTRY(JDOFE,IDOFE,3)=DENTRY(JDOFE,IDOFE,3)+OM3*AH3
                DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+OM3*AH4
        ENDIF   
        IF(IJUMP .EQ. 7)THEN
           AH1=AAH1*OM3+ABH1*OM2+ACH1*OM3+ADH1*OM2
           AH2=0.0D0
           AH3=0.0D0
           AH4=AAH1*OM3+ABH1*OM2+ACH1*OM3+ADH4*OM2
C
           DENTRY(JDOFE,IDOFE,1)=DENTRY(JDOFE,IDOFE,1)+AH1
           DENTRY(JDOFE,IDOFE,2)=DENTRY(JDOFE,IDOFE,2)+AH2
           DENTRY(JDOFE,IDOFE,3)=DENTRY(JDOFE,IDOFE,3)+AH3
           DENTRY(JDOFE,IDOFE,4)=DENTRY(JDOFE,IDOFE,4)+AH4
        ENDIF 
C           
 245        CONTINUE
240    CONTINUE
       ENDIF ! ISPGRAD
      ELSE
C
C ***  Summing up over 1 pair of multiindices
       DO 250 JDOFE=1, IDFL7
       HBASJ1=DBASSF(1,JDOFE)
C
       DO 260 IDOFE=1, IDFL7
       HBASI1=DBASSF(1,IDOFE)
C
       AH=-1D0/THSTEP*HBASJ1*HBASI1
C
c$$$            DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
260    CONTINUE
250    CONTINUE
C
      ENDIF
C
200   CONTINUE
C
      DO 400 JDOFE=1, IDFL7
      DO 400 IDOFE=1, IDFL7
          DENTH1=THSTEP*DENTRY(JDOFE,IDOFE,1)*DJUMP
          DENTH2=THSTEP*DENTRY(JDOFE,IDOFE,2)*DJUMP
          DENTH3=THSTEP*DENTRY(JDOFE,IDOFE,3)*DJUMP
          DENTH4=THSTEP*DENTRY(JDOFE,IDOFE,4)*DJUMP
C 
c      IF (IDEF.GT.0) THEN 
c      IF (IDEF.GT.0 .AND. ILEADJ(2) .NE.0) THEN
       IF (ILEADJ(2) .NE.0) THEN 
       IDFG=INUG7(IDOFE)
       JDFG=INUG7(JDOFE)
      DO  KCOL=KLDA(JDFG),KLDA(JDFG+1)-1
             IF(KCOLA(KCOL) .EQ. IDFG) GOTO 777 
         ENDDO
           D1(JDFG)= D1(JDFG)-DENTH1*U1(IDFG)-DENTH2*U2(IDFG)
           D2(JDFG)= D2(JDFG)-DENTH3*U1(IDFG)-DENTH4*U2(IDFG)
 777       CONTINUE
C
         ENDIF 
C
400   CONTINUE
C
100   CONTINUE
C
C
99999 END
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
************************************************************************
      SUBROUTINE JUMPCP(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                KCOLA,KLDA,KVERT,KMID,KADJ,KMEL,DCORVG,ELE,
     *                COEFFN,IDEF,DCMASS,DNNY)
************************************************************************
*     Purpose: -  Adds the JUMP to the  matrix block 
*
*
*     
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4)
      PARAMETER (NNLEV=9)
C
      DIMENSION A(*),KCOLA(*),KLDA(*),DNNY(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*),KMEL(2,*)
      DIMENSION KADJ(NNVE,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION DENTRY(NNBAS,NNBAS)
      DIMENSION DXG(3),DYG(3),DALPHA(3),ILEADJ(4),INUG(8),INUG7(7)
      DIMENSION DBASSF(3,7),DBASEJUMP(3,8),DXXG(3),DYYG(3)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
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
      EXTERNAL DVISCO,ELE
      SAVE 
C
      ICUB=3
      IDFL7=7
C
C
      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
       CT0=DCMASS/THSTEP
      ELSE
       CT0=0D0
      ENDIF
C
      DNY=NY
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
      CALL ELE(0D0,0D0,-2)
c
C
C *** Loop over all midpoints 
C
      DO 100  IMT=1,NMT
c         print*,'IMT,nmt',IMT,nmt
c      WRITE(*,*),(I,I=1,NMT)
C
C
C *** deter. of the contributed element and the 7x7 numer.
      CALL ADJELEMP(KVERT,KMID,KADJ,KMEL,DCORVG,IMT,NVT,NMT,
     *             	ILEADJ,IVT1,IVT2,INUG,INUG7)
c
c      print*,'ileadj',(ILEADJ(I),I=1,4)
c      DNNY1=DNNY(ILEADJ(1))
c      IF(ILEADJ(2) .NE. 0) THEN
c         DNNY2=DNNY(ILEADJ(2))
c         DNNY3=0.5D0*(DNNY1+DNNY2)
c        ELSE
c            DNNY3=DNNY1
c      ENDIF  
c      DNY=NY*DNNY3 
c
C *** Determine entry positions in matrix
      DO 110 JDOFE=1, IDFL7
         DENTRY(JDOFE,JDOFE)=0D0
      DO 111 IDOFE=1, IDFL7
          DENTRY(JDOFE,IDOFE)=0D0
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices on the side IMT
C
      PPX1=DCORVG(1,IVT1)
      PPY1=DCORVG(2,IVT1)	
      PPX2=DCORVG(1,IVT2)
      PPY2=DCORVG(2,IVT2)
      IF(ILEADJ(3) .EQ. 1)THEN
         PX1=-1.0D0
         PY1=-1.0D0	
         PX2=+1.0D0
         PY2=-1.0D0
      ENDIF 
      IF(ILEADJ(3) .EQ. 2)THEN
         PX1=+1.0D0
         PY1=-1.0D0	
         PX2=+1.0D0
         PY2=+1.0D0
      ENDIF 
      IF(ILEADJ(3) .EQ. 3)THEN
         PX1=+1.0D0
         PY1=+1.0D0	
         PX2=-1.0D0
         PY2=+1.0D0
      ENDIF 
      IF(ILEADJ(3) .EQ. 4)THEN
         PX1=-1.0D0
         PY1=+1.0D0	
         PX2=-1.0D0
         PY2=-1.0D0
      ENDIF    
C
      IF(ILEADJ(4) .EQ. 1)THEN
         PXX1=+1.0D0
         PYY1=-1.0D0	
         PXX2=-1.0D0
         PYY2=-1.0D0
      ENDIF 
      IF(ILEADJ(4) .EQ. 2)THEN
         PXX1=+1.0D0
         PYY1=+1.0D0	
         PXX2=+1.0D0
         PYY2=-1.0D0
      ENDIF 
      IF(ILEADJ(4) .EQ. 3)THEN
         PXX1=-1.0D0
         PYY1=+1.0D0	
         PXX2=+1.0D0
         PYY2=+1.0D0
      ENDIF 
      IF(ILEADJ(4) .EQ. 4)THEN
         PXX1=-1.0D0
         PYY1=-1.0D0	
         PXX2=-1.0D0
         PYY2=+1.0D0
      ENDIF    
C
      HLOCAL=SQRT((PPX2-PPX1)**2+(PPY2-PPY1)**2)	
C
C*** CALCUL OF THE THREE GAUSS POINTS OF NUMERIQUE INTEGRATION 
C***               ON THE SIDE INT 
C                 
       DCONST1= 0D0
       DCONST2= SQRT(3D0/5D0)
       DCONST3= -DCONST2
C
       DXG(1)=0.5D0*(DCONST1*(PX2-PX1)+(PX1+PX2))
       DXG(2)=0.5D0*(DCONST2*(PX2-PX1)+(PX1+PX2))
       DXG(3)=0.5D0*(DCONST3*(PX2-PX1)+(PX1+PX2))
       DYG(1)=0.5D0*(DCONST1*(PY2-PY1)+(PY1+PY2))
       DYG(2)=0.5D0*(DCONST2*(PY2-PY1)+(PY1+PY2)) 
       DYG(3)=0.5D0*(DCONST3*(PY2-PY1)+(PY1+PY2))
C  
       DXXG(1)=0.5D0*(DCONST1*(PXX2-PXX1)+(PXX1+PXX2))
       DXXG(2)=0.5D0*(DCONST2*(PXX2-PXX1)+(PXX1+PXX2))
       DXXG(3)=0.5D0*(DCONST3*(PXX2-PXX1)+(PXX1+PXX2))
       DYYG(1)=0.5D0*(DCONST1*(PYY2-PYY1)+(PYY1+PYY2))
       DYYG(2)=0.5D0*(DCONST2*(PYY2-PYY1)+(PYY1+PYY2)) 
       DYYG(3)=0.5D0*(DCONST3*(PYY2-PYY1)+(PYY1+PYY2))
C  
      DALPHA(1)=8D0/9D0
      DALPHA(2)=5D0/9D0
      DALPHA(3)=DALPHA(2)
C
C *** Loop over all cubature points
C
      DO 200 ICUBP = 1, ICUB
      XI1=DXG(ICUBP)
      XI2=DYG(ICUBP)
      XXI1=DXXG(ICUBP)
      XXI2=DYYG(ICUBP)	
c	
             IF(IJUMP .EQ. 2)OM =DALPHA(ICUBP)*HLOCAL*(HLOCAL**2)
             IF(IJUMP .EQ. 1)OM=DALPHA(ICUBP)	
C
      CALL BASEJUMPP(KVERT,KMID,KADJ,KMEL,DCORVG,IMT,NVT,NMT,
     *             	ILEADJ,INUG,ELE,XI1,XI2,XXI1,XXI2,DBASEJUMP)
C
      CALL BASESFP(KVERT,KMID,KADJ,KMEL,DCORVG,IMT,NVT,NMT,
     *             ILEADJ,INUG,XI1,XI2,XXI1,XXI2,DBASEJUMP,DBASSF)
c
      IF (DCMASS.NE.0D0) THEN
C
        IF (ISPGRAD .EQ. 1 .OR. ISPGRAD .EQ. 0)THEN
C   
C ***  Summing up over all pairs of multiindices
       DO 240 JDOFE=1, IDFL7
       HBASJ1=DBASSF(1,JDOFE)
       HBASJ2=DBASSF(2,JDOFE)
       HBASJ3=DBASSF(3,JDOFE)
C
       DO 245 IDOFE=1, IDFL7
       IF (IDOFE.EQ.JDOFE) THEN
              IF(IJUMP .EQ. 1)AH1= DNY*HBASJ1**2
              IF(IJUMP .EQ. 2)AH1= (HBASJ2**2+HBASJ3**2)
       ELSE
       HBASI1=DBASSF(1,IDOFE)
       HBASI2=DBASSF(2,IDOFE)
       HBASI3=DBASSF(3,IDOFE)
             IF(IJUMP .EQ. 1)AH1= DNY*HBASJ1*HBASI1
             IF(IJUMP .EQ. 2)AH1= (HBASJ2*HBASI2+HBASJ3*HBASI3)
       ENDIF
C
            DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH1
C            
 245        CONTINUE
240    CONTINUE
       ENDIF 
      ELSE
C
C ***  Summing up over 1 pair of multiindices
       DO 250 JDOFE=1, IDFL7
       HBASJ1=DBASSF(1,JDOFE)
C
       DO 260 IDOFE=1, IDFL7
       HBASI1=DBASSF(1,IDOFE)
C
       AH=-1D0/THSTEP*HBASJ1*HBASI1
C
c$$$            DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
260    CONTINUE
250    CONTINUE
C
      ENDIF
C
200   CONTINUE
C
      DO 400 JDOFE=1, IDFL7
      DO 400 IDOFE=1, IDFL7
          DENTH1=THSTEP*DENTRY(JDOFE,IDOFE)*DJUMP
C
      IF (IDEF.LT.2 .AND. ILEADJ(2) .NE.0) THEN
       IDFG=INUG7(IDOFE)
       JDFG=INUG7(JDOFE)
           DO KCOL=KLDA(JDFG),KLDA(JDFG+1)-1
             IF(KCOLA(KCOL) .EQ. IDFG)THEN
                 A(KCOL)=A(KCOL)+REAL(DENTH1)
                 A(3*NA+KCOL)=A(3*NA+KCOL)+REAL(DENTH1)
             ENDIF
           ENDDO 
      ENDIF
C
c         print*,'IMT,nmt2',IMT,nmt
      IF (IDEF.GT.0 .AND. ILEADJ(2) .NE.0) THEN 
       IDFG=INUG7(IDOFE)
       JDFG=INUG7(JDOFE)
           D1(JDFG)= D1(JDFG)-DENTH1*U1(IDFG)
           D2(JDFG)= D2(JDFG)-DENTH1*U2(IDFG)
      ENDIF 
C
400   CONTINUE
C
100   CONTINUE
C
99999 END
C
************************************************************************
      SUBROUTINE ADJELEMP(KVERT,KMID,KADJ,KMEL,DCORVG,IMT,NVT,NMT,
     *             	ILEADJ,IVT1,IVT2,INUG,INUG7)
************************************************************************
*     determination of elements adjacent to the curent midpoint
*     and the thre others contributed in numerical integration   
*     on the curenet side IMT
*     output: 1. Array ILEADJ 
*             2. IVT1, IVT2 : vertices containing curent midpoint
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4,NNBAS=21,NNDER=6)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*),DCORVG(2,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION ILEADJ(4),KMEL(2,*),INUG(8),INUG7(7)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /COAUX1/ KDFG,KDFL,IDFL

      SAVE
C
C
C ***    
C
        DO I=1,4
	   ILEADJ(I)=0
	ENDDO 
C
C ***   DETRMINATION OF THE TWO ADJ ELEMENT TO THE MID IMT   
c
         IEL1=KMEL(1,IMT) 
         IEL2=KMEL(2,IMT) 
         ILEADJ(1)=KMEL(1,IMT) 
         ILEADJ(2)=KMEL(2,IMT) 
         IEL=IEL1

C ***   Determine of the other adjacent elenents to IEL   
C
C ********* Determine the coresponding vertices         
         DO 10 IVE =1,NNVE
           IVEN=IVE+1
           IF (IVEN.EQ.5) IVEN=1
         	  IMID1=KMID(IVE,IEL)-NVT
          	  IF(IMID1 .EQ. IMT)THEN
         		IVT1=KVERT(IVE ,IEL)
               		IVT2=KVERT(IVEN,IEL)
	      	 ENDIF
10       CONTINUE
C	
        DO IN=1,8
           INUG(IN)=0
        ENDDO
C
         DO  IVE=1,NNVE
             IVEN=IVE+1
             IF(IVEN .EQ. 5)IVEN=1
             IMID=KMID(IVE,ILEADJ(1))-NVT
             IF(IMID .EQ. IMT)ILEADJ(3)=IVE
             INUG(IVE)=IMID
         ENDDO
C
         IF(ILEADJ(2) .NE.0)THEN
            DO IVE=1,NNVE
               IVEN=IVE+1
               IF(IVEN .EQ. 5)IVEN=1
               IMID=KMID(IVE,ILEADJ(2))-NVT
               IF(IMID .EQ. IMT)ILEADJ(4)=IVE
               INUG(NNVE+IVE)=IMID
            ENDDO   
         ENDIF
C
         DO IN=1,7
            INUG7(IN)=0
         ENDDO
         DO IVE=1,NNVE
             IVEN=IVE+1
             IF(IVEN .EQ. 5)IVEN=1
             IMID=KMID(IVE,ILEADJ(1))-NVT
             INUG7(IVE)=IMID
         ENDDO
         IF(ILEADJ(2) .NE.0)THEN
         ICONT=0
         DO IVE=1,NNVE
            IVEN=IVE+1
             IF(IVEN .EQ. 5)IVEN=1
             IMID=KMID(IVE,ILEADJ(2))-NVT
             IF(IMID .NE. IMT)THEN
                ICONT=ICONT+1
                INUG7(NNVE+ICONT)=IMID
             ENDIF
         ENDDO
         ENDIF
c         PRINT*,'INUG7',(INUG7(I),I=1,7)
C
        END
C

************************************************************************
      SUBROUTINE BASEJUMPP(KVERT,KMID,KADJ,KMEL,DCORVG,IMT,NVT,NMT,
     *             	ILEADJ,INUG,ELE,XI1,XI2,XXI1,XXI2,DBASEJUMP)
************************************************************************
*      
*     
*        
*     
*    
*             
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4,NNBAS=21,NNDER=6)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*),DCORVG(2,*)
      DIMENSION ILEADJ(4),KMEL(2,*),DBASEJUMP(3,8),INUG(8)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      SAVE
      EXTERNAL ELE
C
C ***    
C 
      IELTYP=-1
      DO I=1,8
         DBASEJUMP(1,I)=0D0
         DBASEJUMP(2,I)=0D0
         DBASEJUMP(3,I)=0D0
      ENDDO
C
      CALL ELE(0D0,0D0,-2)
C
      IEL=ILEADJ(1)
C
C *** Evaluation of coordinates of the vertices
      DO  IVE = 1, NNVE
         JP=KVERT(IVE,IEL)
         KVE(IVE)=JP
         DX(IVE)=DCORVG(1,JP)
         DY(IVE)=DCORVG(2,JP)
      ENDDO
      CALL ELE(0D0,0D0,-2)
      CALL ELE(XI1,XI2,-3)
         DO  IVE=1,NNVE
            DBASEJUMP(1,IVE)=DBAS(IVE,1)
            DBASEJUMP(2,IVE)=DBAS(IVE,2)
            DBASEJUMP(3,IVE)=DBAS(IVE,3)
         ENDDO 
C
         CALL ELE(0D0,0D0,-2)
         IF(ILEADJ(2) .NE.0)THEN
         IEL=ILEADJ(2)   
C
C *** Evaluation of coordinates of the vertices
         DO  IVE = 1, NNVE
            JP=KVERT(IVE,IEL)
            KVE(IVE)=JP
            DX(IVE)=DCORVG(1,JP)
            DY(IVE)=DCORVG(2,JP)
          ENDDO 
            CALL ELE(0d0,0d0,-2)
            CALL ELE(XXI1,XXI2,-3) 
            DO IVE=1,NNVE
                  DBASEJUMP(1,NNVE+IVE)=DBAS(IVE,1)
                  DBASEJUMP(2,NNVE+IVE)=DBAS(IVE,2)
                  DBASEJUMP(3,NNVE+IVE)=DBAS(IVE,3)
            ENDDO 
      ENDIF      

c          PRINT*,'DBASJUMP',(DBASEJUMP(I),I=1,8)
C
99999   END
C
************************************************************************
      SUBROUTINE BASESFP(KVERT,KMID,KADJ,KMEL,DCORVG,IMT,NVT,NMT,
     *             ILEADJ,INUG,XI1,XI2,XXI1,XXI2,DBASEJUMP,DBASSF)
************************************************************************
*      
*     
*        
*     
*    
*             
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4,NNBAS=21,NNDER=6)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*),DCORVG(2,*)
      DIMENSION ILEADJ(4),KMEL(2,*),DBASEJUMP(3,8),INUG(8),DBASSF(3,7)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      SAVE

C
C ***    
C 
      IELTYP=-1
      DO I=1,7
         DBASSF(1,I)=0D0
         DBASSF(2,I)=0D0
         DBASSF(3,I)=0D0
      ENDDO

           DO I=1,4
              IMID=INUG(I)
              IF(IMID .EQ. IMT)IND=I
              DBASSF(1,I)=DBASEJUMP(1,I)
              DBASSF(2,I)=DBASEJUMP(2,I)
              DBASSF(3,I)=DBASEJUMP(3,I)
           ENDDO
C
        ICONT=0
        DO I=5,8
           IMID=INUG(I)
           IF(IMID .EQ. IMT)THEN
              DBASSF(1,IND)=DBASSF(1,IND)-DBASEJUMP(1,I)
              DBASSF(2,IND)=DBASSF(2,IND)-DBASEJUMP(2,I)
              DBASSF(3,IND)=DBASSF(3,IND)-DBASEJUMP(3,I)
              GOTO 1
           ENDIF 
           ICONT=ICONT+1 
           DBASSF(1,NNVE+ICONT)=-DBASEJUMP(1,I)
           DBASSF(2,NNVE+ICONT)=-DBASEJUMP(2,I)
           DBASSF(3,NNVE+ICONT)=-DBASEJUMP(3,I)
 1         CONTINUE
        ENDDO   
c        PRINT*,'DBASSF(I)',(DBASSF(I),I=1,7)
c        READ(*,*)
           
C
99999   END
C     
************************************************************************
      SUBROUTINE JUMPLCP(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                KCOLA,KLDA,KVERT,KMID,KADJ,KMEL,DCORVG,ELE,
     *                COEFFN,IDEF,DCMASS,DNNY)
************************************************************************
*     Purpose: -  Adds the JUMP to the  matrix block 
*
*
*     
************************************************************************
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      REAL A
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4)
      PARAMETER (NNLEV=9)
C
      DIMENSION A(*),KCOLA(*),KLDA(*),DNNY(*)
      DIMENSION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*),D1(*),D2(*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),DCORVG(2,*),KMEL(2,*)
      DIMENSION KADJ(NNVE,*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION DENTRY(NNBAS,NNBAS)
      DIMENSION DXG(3),DYG(3),DALPHA(3),ILEADJ(4),INUG(8),INUG7(7)
      DIMENSION DBASSF(3,7),DBASEJUMP(3,8),DXXG(3),DYYG(3)	
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
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
      EXTERNAL DVISCO,ELE
      SAVE 
C
      ICUB=3
      IDFL7=7
C
      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
       CT0=DCMASS/THSTEP
      ELSE
       CT0=0D0
      ENDIF
C
      DNY=NY
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
C
************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************

      CALL ELE(0D0,0D0,-2)

C
C *** Loop over all midpoints 
C
      DO 100  IMT=1,NMT
C
C
C *** deter. of the contributed element and the 7x7 numer.
      CALL ADJELEMP(KVERT,KMID,KADJ,KMEL,DCORVG,IMT,NVT,NMT,
     *             	ILEADJ,IVT1,IVT2,INUG,INUG7)
C
c      DNNY1=DNNY(ILEADJ(1))
c      IF(ILEADJ(2) .NE. 0) THEN
c         DNNY2=DNNY(ILEADJ(2))
c         DNNY3=0.5D0*(DNNY1+DNNY2)
c         ELSE
c             DNNY3=DNNY1
c      ENDIF  
c      DNY=NY*DNNY3
c
C 	
C *** Determine entry positions in matrix
      DO 110 JDOFE=1, IDFL7
         DENTRY(JDOFE,JDOFE)=0D0
      DO 111 IDOFE=1, IDFL7
          DENTRY(JDOFE,IDOFE)=0D0
111   CONTINUE
110   CONTINUE
C
C *** Evaluation of coordinates of the vertices on the side IMT
C
      PPX1=DCORVG(1,IVT1)
      PPY1=DCORVG(2,IVT1)	
      PPX2=DCORVG(1,IVT2)
      PPY2=DCORVG(2,IVT2)
C
      IF(ILEADJ(3) .EQ. 1)THEN
         PX1=-1.0D0
         PY1=-1.0D0	
         PX2=+1.0D0
         PY2=-1.0D0
      ENDIF 
      IF(ILEADJ(3) .EQ. 2)THEN
         PX1=+1.0D0
         PY1=-1.0D0	
         PX2=+1.0D0
         PY2=+1.0D0
      ENDIF 
      IF(ILEADJ(3) .EQ. 3)THEN
         PX1=+1.0D0
         PY1=+1.0D0	
         PX2=-1.0D0
         PY2=+1.0D0
      ENDIF 
      IF(ILEADJ(3) .EQ. 4)THEN
         PX1=-1.0D0
         PY1=+1.0D0	
         PX2=-1.0D0
         PY2=-1.0D0
      ENDIF    
C
      IF(ILEADJ(4) .EQ. 1)THEN
         PXX1=+1.0D0
         PYY1=-1.0D0	
         PXX2=-1.0D0
         PYY2=-1.0D0
      ENDIF 
      IF(ILEADJ(4) .EQ. 2)THEN
         PXX1=+1.0D0
         PYY1=+1.0D0	
         PXX2=+1.0D0
         PYY2=-1.0D0
      ENDIF 
      IF(ILEADJ(4) .EQ. 3)THEN
         PXX1=-1.0D0
         PYY1=+1.0D0	
         PXX2=+1.0D0
         PYY2=+1.0D0
      ENDIF 
      IF(ILEADJ(4) .EQ. 4)THEN
         PXX1=-1.0D0
         PYY1=-1.0D0	
         PXX2=-1.0D0
         PYY2=+1.0D0
      ENDIF    
c
      HLOCAL=SQRT((PPX2-PPX1)**2+(PPY2-PPY1)**2)	
C
C*** CALCUL OF THE THREE GAUSS POINTS OF NUMERIQUE INTEGRATION 
C***               ON THE SIDE INT 
C                 
       DCONST1= 0D0
       DCONST2= SQRT(3D0/5D0)
       DCONST3= -DCONST2
C
       DXG(1)=0.5D0*(DCONST1*(PX2-PX1)+(PX1+PX2))
       DXG(2)=0.5D0*(DCONST2*(PX2-PX1)+(PX1+PX2))
       DXG(3)=0.5D0*(DCONST3*(PX2-PX1)+(PX1+PX2))
       DYG(1)=0.5D0*(DCONST1*(PY2-PY1)+(PY1+PY2))
       DYG(2)=0.5D0*(DCONST2*(PY2-PY1)+(PY1+PY2)) 
       DYG(3)=0.5D0*(DCONST3*(PY2-PY1)+(PY1+PY2))
C  
       DXXG(1)=0.5D0*(DCONST1*(PXX2-PXX1)+(PXX1+PXX2))
       DXXG(2)=0.5D0*(DCONST2*(PXX2-PXX1)+(PXX1+PXX2))
       DXXG(3)=0.5D0*(DCONST3*(PXX2-PXX1)+(PXX1+PXX2))
       DYYG(1)=0.5D0*(DCONST1*(PYY2-PYY1)+(PYY1+PYY2))
       DYYG(2)=0.5D0*(DCONST2*(PYY2-PYY1)+(PYY1+PYY2)) 
       DYYG(3)=0.5D0*(DCONST3*(PYY2-PYY1)+(PYY1+PYY2))
C
      DALPHA(1)=8D0/9D0
      DALPHA(2)=5D0/9D0
      DALPHA(3)=DALPHA(2)
C
C *** Loop over all cubature points
C
      DO 200 ICUBP = 1, ICUB
      XI1=DXG(ICUBP)
      XI2=DYG(ICUBP)	
      XXI1=DXXG(ICUBP)
      XXI2=DYYG(ICUBP)	
c
             IF(IJUMP .EQ. 2)OM =DALPHA(ICUBP)*HLOCAL*(HLOCAL**2)	
             IF(IJUMP .EQ. 1)OM=DALPHA(ICUBP)
C
      CALL BASEJUMPP(KVERT,KMID,KADJ,KMEL,DCORVG,IMT,NVT,NMT,
     *             	ILEADJ,INUG,ELE,XI1,XI2,XXI1,XXI2,DBASEJUMP)
C
      CALL BASESFP(KVERT,KMID,KADJ,KMEL,DCORVG,IMT,NVT,NMT,
     *             ILEADJ,INUG,XI1,XI2,XXI1,XXI2,DBASEJUMP,DBASSF)
c
      IF (DCMASS.NE.0D0) THEN
C
        IF (ISPGRAD .EQ. 1 .OR. ISPGRAD .EQ. 0)THEN
C
C***   use symmetric part of grad
C   
C ***  Summing up over all pairs of multiindices
       DO 240 JDOFE=1, IDFL7
       HBASJ1=DBASSF(1,JDOFE)
       HBASJ2=DBASSF(2,JDOFE)
       HBASJ3=DBASSF(3,JDOFE)
C
       DO 245 IDOFE=1, IDFL7
       IF (IDOFE.EQ.JDOFE) THEN
             IF(IJUMP .EQ. 1)AH1= DNY*HBASJ1**2
             IF(IJUMP .EQ. 2)AH1= (HBASJ2**2+HBASJ3**2)
       ELSE
       HBASI1=DBASSF(1,IDOFE)
       HBASI2=DBASSF(2,IDOFE)
       HBASI3=DBASSF(3,IDOFE)
             IF(IJUMP .EQ. 1)AH1= DNY*HBASJ1*HBASI1
             IF(IJUMP .EQ. 2)AH1= (HBASJ2*HBASI2+HBASJ3*HBASI3)
       ENDIF
C
            DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH1
C           
 245        CONTINUE
240    CONTINUE
       ENDIF ! ISPGRAD
      ELSE
C
C ***  Summing up over 1 pair of multiindices
       DO 250 JDOFE=1, IDFL7
       HBASJ1=DBASSF(1,JDOFE)
C
       DO 260 IDOFE=1, IDFL7
       HBASI1=DBASSF(1,IDOFE)
C
       AH=-1D0/THSTEP*HBASJ1*HBASI1
C
c$$$            DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
260    CONTINUE
250    CONTINUE
C
      ENDIF
C
200   CONTINUE
C
      DO 400 JDOFE=1, IDFL7
      DO 400 IDOFE=1, IDFL7
          DENTH1=THSTEP*DENTRY(JDOFE,IDOFE)*DJUMP
C 
c      IF (IDEF.GT.0) THEN 
c      IF (IDEF.GT.0 .AND. ILEADJ(2) .NE.0) THEN
       IF (ILEADJ(2) .NE.0) THEN 
       IDFG=INUG7(IDOFE)
       JDFG=INUG7(JDOFE)
      DO  KCOL=KLDA(JDFG),KLDA(JDFG+1)-1
             IF(KCOLA(KCOL) .EQ. IDFG) GOTO 777 
         ENDDO
           D1(JDFG)= D1(JDFG)-DENTH1*U1(IDFG)
           D2(JDFG)= D2(JDFG)-DENTH1*U2(IDFG)
 777       CONTINUE
C
         ENDIF 
C
400   CONTINUE
C
100   CONTINUE
C
99999 END
C



