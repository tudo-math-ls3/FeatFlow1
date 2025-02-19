************************************************************************
      SUBROUTINE VANCAS (U1,U2,U3,P,U1OLD,U2OLD,U3OLD,POLD,F1,F2,F3,FP,
     *                   A,KCOLA,KLDA,B1,B2,B3,KCOLB,KLDB,NU,NP,
     *                   KABD,KVERT,KAREA,KNPR,DMAXU,DMAXP)
************************************************************************
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DOUBLE PRECISION  BB1,BB2,BB3
      REAL  A,B1,B2,B3
C
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION U1(*),U2(*),U3(*),P(*)
      DIMENSION U1OLD(*),U2OLD(*),U3OLD(*),POLD(*)
      DIMENSION F1(*),F2(*),F3(*),FP(*)
      DIMENSION A(*),KCOLA(*),KLDA(*),B1(*),KCOLB(*),KLDB(*),B2(*),B3(*)
      DIMENSION KABD(*),KVERT(NNVE,*),KAREA(NNAE,*),KNPR(*)
C
C *** Local arrays for informations about one element
      DIMENSION AA(6),BB1(6),BB2(6),BB3(6),FF1(6),FF2(6),FF3(6)
      DIMENSION IU(6),BDBC(6)
C
C *** Usual data for mesh management
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
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
C *** Update  UOLD:=U
      CALL  LCP1 (U1,U1OLD,NU)
      CALL  LCP1 (U2,U2OLD,NU)
      CALL  LCP1 (U3,U3OLD,NU)
      CALL  LCP1 (P ,POLD ,NP)
C
C=======================================================================
C     Block Gauss-Seidel on Schur Complement
C=======================================================================
C *** Loop over all elements
C
      DO 1  IEL=1,NEL
      FFP =FP(IEL)
C
C-----------------------------------------------------------------------
C *** Loop over all 6 U-nodes of that element
      DO 11  II=1,6
C
      IAT=KAREA(II,IEL)
      I=IAT
      IU(II)=I
      IF (KNPR(IAT+NVT).LE.0)  THEN
         BDBC(II)=.FALSE.
      ELSE
         BDBC(II)=.TRUE.
      ENDIF
C
C *** Put on AA(.) the diagonal entry of matrix A
      IA1=KLDA(I)
      AA(II)=DBLE(A(IA1))
C
C *** Initial setting of FF1,FF2,FF3
      FF1(II)=F1(I)
      FF2(II)=F2(I)
      FF3(II)=F3(I)
C
      IF (BDBC(II)) GOTO 13

C *** Put all off-diagonal entries of A to the right hand side
      IA2=KLDA(I+1)-1
      DO 110  IA=IA1+1,IA2
      J=KCOLA(IA)
      AOFF=DBLE(A(IA))
      FF1(II)=FF1(II)-AOFF*U1(J)
      FF2(II)=FF2(II)-AOFF*U2(J)
      FF3(II)=FF3(II)-AOFF*U3(J)
 110  CONTINUE
C
C-----------------------------------------------------------------------
C *** Get BB1,BB2,BB3 and modify FF1,FF2,FF3
      IB1=KLDB(I)
      IB2=KLDB(I+1)-1
      JP1=KCOLB(IB1)
      JP2=KCOLB(IB2)
C
      IF (JP1.EQ.IEL)  THEN
          BB1(II)=DBLE(B1(IB1))
          BB2(II)=DBLE(B2(IB1))
          BB3(II)=DBLE(B3(IB1))
          FF1(II)=FF1(II)-DBLE(B1(IB2))*P(JP2)
          FF2(II)=FF2(II)-DBLE(B2(IB2))*P(JP2)
          FF3(II)=FF3(II)-DBLE(B3(IB2))*P(JP2)
          IB3    =IB1
          JP3    =JP1
      ELSE IF (JP2.EQ.IEL)  THEN
          BB1(II)=DBLE(B1(IB2))
          BB2(II)=DBLE(B2(IB2))
          BB3(II)=DBLE(B3(IB2))
          FF1(II)=FF1(II)-DBLE(B1(IB1))*P(JP1)
          FF2(II)=FF2(II)-DBLE(B2(IB1))*P(JP1)
          FF3(II)=FF3(II)-DBLE(B3(IB1))*P(JP1)
          IB3    =IB2
          JP3    =JP2
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
      IB=KLDB(I)
      JP=KCOLB(IB)
      BB1(II)=DBLE(B1(IB))
      BB2(II)=DBLE(B2(IB))
      BB3(II)=DBLE(B3(IB))
C
11    CONTINUE
C
C-----------------------------------------------------------------------
C *** Calling the subroutine ELUPDT for element update 
      CALL  ELUPDT (U1,U2,U3,P,IEL,IU,BDBC,AA,BB1,BB2,BB3,FF1,FF2,
     *              FF3,FFP)
C
1     CONTINUE
C
C=======================================================================
C
C *** Relaxation   U:= RLXSM*U + (1-RLXSM)*UOLD
      IF (RLXSM.NE.1D0) THEN
       A1= 1D0-RLXSM
       A2=     RLXSM
       CALL  LLC1 (U1OLD,U1,NU,A1,A2)
       CALL  LLC1 (U2OLD,U2,NU,A1,A2)
       CALL  LLC1 (U3OLD,U3,NU,A1,A2)
       CALL  LLC1 (POLD ,P ,NP,A1,A2)
      ENDIF
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE VANCAE (U1,U2,U3,P,U1OLD,U2OLD,U3OLD,POLD,F1,F2,F3,FP,
     *                   A,KCOLA,KLDA,B1,B2,B3,KCOLB,KLDB,NU,NP,
     *                   KABD,KVERT,KAREA,KNPR,RLXPAR,DMAXU,DMAXP)
************************************************************************
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DOUBLE PRECISION  BB1,BB2,BB3,FFP
      REAL  A,B1,B2,B3
C
      PARAMETER (NNVE=8,NNAE=6)
      DIMENSION U1(*),U2(*),U3(*),P(*)
      DIMENSION U1OLD(*),U2OLD(*),U3OLD(*),POLD(*)
      DIMENSION F1(*),F2(*),F3(*),FP(*)
      DIMENSION A(*),KCOLA(*),KLDA(*),B1(*),KCOLB(*),KLDB(*),B2(*),B3(*)
      DIMENSION KABD(*),KVERT(NNVE,*),KAREA(NNAE,*),KNPR(*)
C
C *** Local arrays for informations about one element
      DIMENSION AA(6),BB1(6),BB2(6),BB3(6),FF1(6),FF2(6),FF3(6)
      DIMENSION IU(6),BDBC(6)
C
C *** Usual data for mesh management
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
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
      DMAXU=0D0
      DMAXP=0D0
C
C *** Update  UOLD:=U
      CALL  LCP1 (U1,U1OLD,NU)
      CALL  LCP1 (U2,U2OLD,NU)
      CALL  LCP1 (U3,U3OLD,NU)
      CALL  LCP1 (P ,POLD ,NP)
C
C=======================================================================
C     Block Gauss-Seidel on Schur Complement
C=======================================================================
C *** Loop over all elements
C
      DO 1  IEL=1,NEL
      FFP =FP(IEL)
      FFPH=FFP
C
C-----------------------------------------------------------------------
C *** Loop over all 6 U-nodes of that element
      DO 11  II=1,6
C
      IAT=KAREA(II,IEL)
      I=IAT
      IU(II)=I
      IF (KNPR(IAT+NVT).LE.0)  THEN
         BDBC(II)=.FALSE.
      ELSE
         BDBC(II)=.TRUE.
      ENDIF
C
C *** Put on AA(.) the diagonal entry of matrix A
      IA1=KLDA(I)
      AA(II)=DBLE(A(IA1))
C
C *** Initial setting of FF1,FF2,FF3
      FF1(II)=F1(I)
      FF2(II)=F2(I)
      FF3(II)=F3(I)
C
      IF (BDBC(II)) GOTO 13
C
C *** Put all off-diagonal entries of A to the right hand side
      IA2=KLDA(I+1)-1
      DO 110  IA=IA1+1,IA2
      J=KCOLA(IA)
      AOFF=DBLE(A(IA))
      FF1(II)=FF1(II)-AOFF*U1(J)
      FF2(II)=FF2(II)-AOFF*U2(J)
      FF3(II)=FF3(II)-AOFF*U3(J)
 110  CONTINUE
C
C-----------------------------------------------------------------------
C *** Get BB1,BB2,BB3 and modify FF1,FF2,FF3
      IB1=KLDB(I)
      IB2=KLDB(I+1)-1
      JP1=KCOLB(IB1)
      JP2=KCOLB(IB2)
C
      IF (JP1.EQ.IEL)  THEN
          BB1(II)=DBLE(B1(IB1))
          BB2(II)=DBLE(B2(IB1))
          BB3(II)=DBLE(B3(IB1))
          FF1(II)=FF1(II)-DBLE(B1(IB2))*P(JP2)
          FF2(II)=FF2(II)-DBLE(B2(IB2))*P(JP2)
          FF3(II)=FF3(II)-DBLE(B3(IB2))*P(JP2)
          IB3    =IB1
          JP3    =JP1
      ELSE IF (JP2.EQ.IEL)  THEN
          BB1(II)=DBLE(B1(IB2))
          BB2(II)=DBLE(B2(IB2))
          BB3(II)=DBLE(B3(IB2))
          FF1(II)=FF1(II)-DBLE(B1(IB1))*P(JP1)
          FF2(II)=FF2(II)-DBLE(B2(IB1))*P(JP1)
          FF3(II)=FF3(II)-DBLE(B3(IB1))*P(JP1)
          IB3    =IB2
          JP3    =JP2
      ELSE
          WRITE(MTERM,*) 'ERROR in SMOOTH: IEL entry in B not found'
          STOP
      ENDIF
C
      AH=DBLE(A(IA1))
      ICOLA1=KCOLA(IA1)
      PJP=P(JP3)
      DMAXU=MAX(DMAXU,ABS(FF1(II)-AH*U1(ICOLA1)-DBLE(B1(IB3))*PJP))
      DMAXU=MAX(DMAXU,ABS(FF2(II)-AH*U2(ICOLA1)-DBLE(B2(IB3))*PJP))
      DMAXU=MAX(DMAXU,ABS(FF3(II)-AH*U3(ICOLA1)-DBLE(B3(IB3))*PJP))
      FFPH=FFPH-BB1(II)*U1(IU(II))-BB2(II)*U2(IU(II))-BB3(II)*U3(IU(II))
C
      GOTO 11

C-----------------------------------------------------------------------
C *** The case that II is a Dirichlet boundary node
13    CONTINUE
      IB=KLDB(I)
      JP=KCOLB(IB)
      BB1(II)=DBLE(B1(IB))
      BB2(II)=DBLE(B2(IB))
      BB3(II)=DBLE(B3(IB))
C
      FFPH=FFPH-BB1(II)*U1(IU(II))
     *         -BB2(II)*U2(IU(II))-BB3(II)*U3(IU(II))
C
11    CONTINUE
C
      DMAXP=MAX(DMAXP,ABS(FFPH))
C
C-----------------------------------------------------------------------
C *** Calling the subroutine ELUPDT for element update 
      CALL  ELUPDT (U1,U2,U3,P,IEL,IU,BDBC,AA,BB1,BB2,BB3,FF1,FF2,
     *              FF3,FFP)
C
1     CONTINUE
C
C=======================================================================
C
C *** Relaxation   U:= RLXPAR*U + (1-RLXPAR)*UOLD
      IF (RLXPAR.NE.1D0) THEN
       A1= 1D0-RLXPAR
       A2=     RLXPAR
       CALL  LLC1 (U1OLD,U1,NU,A1,A2)
       CALL  LLC1 (U2OLD,U2,NU,A1,A2)
       CALL  LLC1 (U3OLD,U3,NU,A1,A2)
       CALL  LLC1 (POLD ,P ,NP,A1,A2)
      ENDIF
C
C
C
      END
C
************************************************************************
      SUBROUTINE  ELUPDT  (U1,U2,U3,P,IEL,IU,BDBC,AA,BB1,BB2,BB3,
     *                     FF1,FF2,FF3,FFP)
************************************************************************
*    Purpose: - Executes block Gauss-Seidel update on (U1,U2,U3,P)
*               for all unknowns of the element IEL
*-----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DOUBLE PRECISION  BB1,BB2,BB3
C
      DIMENSION U1(*),U2(*),U3(*),P(*)
C
C *** Local arrays for informations about the element
      DIMENSION AA(6),BB1(6),BB2(6),BB3(6),FF1(6),FF2(6),FF3(6),
     *          IU(6),BDBC(6)
      DIMENSION AI(6),DD1(6),DD2(6),DD3(6),UU1(6),UU2(6),UU3(6)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      SAVE 
C
C-----------------------------------------------------------------------
C
      DO 1  II=1,6
      IF (BDBC(II))  THEN
          DD1(II)=0D0
          DD2(II)=0D0
          DD3(II)=0D0
      ELSE
          DD1(II)=BB1(II)
          DD2(II)=BB2(II)
          DD3(II)=BB3(II)
      ENDIF
      IF (DABS(AA(II)).LT.1.D-10)  THEN
         WRITE(MTERM,*)'ERROR in ELUPDT: diagonal entry is nearly zero'
         WRITE(MTERM,*)'II, I, BNDR(II) : ',II,IU(II),BDBC(II)
         STOP
      ENDIF
      AI(II)=1.D0/AA(II)
  1   CONTINUE
C
C *** Factorization loop
C
      DP=0.D0
      DO 2  II=1,6
      DP =DP -AI(II)*(BB1(II)*DD1(II)+BB2(II)*DD2(II)+BB3(II)*DD3(II))
      FFP=FFP-AI(II)*(BB1(II)*FF1(II)+BB2(II)*FF2(II)+BB3(II)*FF3(II))
  2   CONTINUE
C
C *** Solution loop
C
      IF (DABS(DP).LT.1.D-10)  THEN
         WRITE(MTERM,*)'ERROR in ELUPDT: DP is nearly zero'
         STOP
      ENDIF
      PP=FFP/DP
      DO 3  II=1,6
      UU1(II)=AI(II)*(FF1(II)-DD1(II)*PP)
      UU2(II)=AI(II)*(FF2(II)-DD2(II)*PP)
      UU3(II)=AI(II)*(FF3(II)-DD3(II)*PP)
  3   CONTINUE
C
C *** Updating the global solution vector (U1,U2,U3,P)
C
      P(IEL)=PP
      DO 4  II=1,6
      IF (BDBC(II)) GOTO 4
      I=IU(II)
      U1(I)=UU1(II)
      U2(I)=UU2(II)
      U3(I)=UU3(II)
  4   CONTINUE
C
C
C
      END
