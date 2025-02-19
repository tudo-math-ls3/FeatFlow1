************************************************************************
      SUBROUTINE  VANCR2old (U1,U2,P,F1,F2,FP,
     *                    A,KACOL,KALD,B1,B2,KBCOL,KBLD,
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
      DOUBLE PRECISION BB1,BB2
      REAL  A,B1,B2

      INTEGER AKTELE
C
c      PARAMETER (NNVE=4)
      DIMENSION U1(*),U2(*),P(*)
      DIMENSION F1(*),F2(*),FP(*)
      DIMENSION A(*),KACOL(*),KALD(*),B1(*),B2(*),KBCOL(*),KBLD(*)
      DIMENSION KMBD(*),KVERT(NNVE,*),KMID(NNVE,*),KNPR(*),DCORVG(2,*)
C
C *** Local arrays for informations about one element
      DIMENSION AA(8,8),BB1(4),BB2(4),FF(8)
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
          FF(II)  =FF(II)  -DBLE(B1(IB2))*PJP
          FF(4+II)=FF(4+II)-DBLE(B2(IB2))*PJP
      ELSE IF (JP2.EQ.IEL)  THEN
          PJP=P(JP1)
          BB1(II)=DBLE(B1(IB2))
          BB2(II)=DBLE(B2(IB2))
          FF(II)  =FF(II)  -DBLE(B1(IB1))*PJP
          FF(4+II)=FF(4+II)-DBLE(B2(IB1))*PJP
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
      SUBROUTINE  VANCAR (U1,U2,P,F1,F2,FP,
     *                    A,KACOL,KALD,B1,B2,KBCOL,KBLD,
     *                    KMBD,KVERT,KMID,KNPR,NMBD,Inum,AKTELE)
************************************************************************
C
C----------------------------------------------------------------------
C Purpose:  Vanca-type Glaetter. 
C     Auf jedem Element werden zwei lokale 
c     4x4 Matrizen (AA(i,j)) EXAKT invertiert
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'
      INCLUDE 'bouss.inc'
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION BB1,BB2
      REAL  A,B1,B2
      INTEGER AKTELE
C
c      PARAMETER (NNVE=4)
      DIMENSION U1(*),U2(*),P(*)
      DIMENSION F1(*),F2(*),FP(*)
      DIMENSION A(*),KACOL(*),KALD(*),B1(*),B2(*),KBCOL(*),KBLD(*)
      DIMENSION KMBD(*),KVERT(NNVE,*),KMID(NNVE,*),KNPR(*)
C
C *** Local arrays for informations about one element
      DIMENSION AA1(4,4),AA2(4,4),BB1(4),BB2(4),FF1(4),FF2(4)
      DIMENSION IU(4),BDBC(4)
C
      NA=KALD(NU+1)-1
C
C=======================================================================
C     Block Gauss-Seidel on Schur Complement
C=======================================================================
C *** Loop over all elements (or only one)
C
      IF (INUM.eq.0) THEN
         Nloop=1
      ELSE
         NLOOP=NEL
      ENDIF

      CALL VKADD (AVK1,AVK2)
      DO 1  IELN=1,NLOOP!NELC

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
      AA1(II,II)=AVK1*DBLE(A(IA1))+AVK2*DBLE(A(4*NA+ia1))
      AA2(II,II)=AVK1*DBLE(A(3*NA+IA1))+AVK2*DBLE(A(4*NA+ia1))
C

C *** Initial setting of FF1,FF2
      FF1(II)=F1(I)
      FF2(II)=F2(I)
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
      AOFF2=AVK1*DBLE(A(3*NA+IA))+AVK2*DBLE(A(4*NA+ia))
c      write (*,*) ia,'aoff', aoff1,aoff2
      FF1(II)=FF1(II)-AOFF1*U1(J)
      FF2(II)=FF2(II)-AOFF2*U2(J)
      ENDIF
110   CONTINUE
      DO 111  IA=IA1,IA2
      J=KACOL(IA)
      AOFF1=AVK1*DBLE(A(NA+IA))
      AOFF2=AVK1*DBLE(A(2*NA+IA))
c     
      FF1(II)=FF1(II)-AOFF1*U2(J)
      FF2(II)=FF2(II)-AOFF2*U1(J)
111   CONTINUE
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
          FF1(II)=FF1(II)-DBLE(B1(IB2))*PJP
          FF2(II)=FF2(II)-DBLE(B2(IB2))*PJP
      ELSE IF (JP2.EQ.IEL)  THEN
          PJP=P(JP1)
          BB1(II)=DBLE(B1(IB2))
          BB2(II)=DBLE(B2(IB2))
          FF1(II)=FF1(II)-DBLE(B1(IB1))*PJP
          FF2(II)=FF2(II)-DBLE(B2(IB1))*PJP
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
      AA1(II,JJ)=0D0
      AA2(II,JJ)=0D0
      IF (BDBC(II)) GOTO 22
C
      J=IU(JJ)
      IA1=KALD(I)
      IA2=KALD(I+1)-1
C
      DO 220  IA=IA1+1,IA2
      JH=KACOL(IA)
      IF (J.EQ.JH) THEN
c       AA(jj,ii)=DBLE(A(IA))
       AA1(II,JJ)=AVK1*DBLE(A(IA))+AVK2*DBLE(A(4*NA+ia))
       AA2(II,JJ)=AVK1*DBLE(A(3*NA+IA))+AVK2*DBLE(A(4*NA+ia))
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
      CALL  ELUPDN (U1,U2,P,IEL,IU,BDBC,AA1,AA2,BB1,BB2,FF1,FF2,FFP)
C
1     CONTINUE
C
C=======================================================================
C
C
      END
c
************************************************************************
      SUBROUTINE  VANCAN (U1,U2,P,F1,F2,FP,
     *                    A,KACOL,KALD,B1,B2,KBCOL,KBLD,
     *                    KMBD,KVERT,KMID,KNPR,NMBD)
************************************************************************
C
C-----------------------------------------------------------------------
C Purpose:  Vanca-type Glaetter. 
C     Auf jedem Element wird die lokale Matrix
c     (AA(i,j)) EXAKT invertiert.
C     Unterschied zu VANCAR: Es werden die Einträge
C     aus den 'Nachbarzellen nicht nach rechts gebracht
c     Fuer den Deformationstensor wird der Glaetter
c     als nicht so clever erachtet
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'
      INCLUDE 'bouss.inc'
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION BB1,BB2
      REAL  A,B1,B2
C
c      PARAMETER (NNVE=4)
      DIMENSION U1(*),U2(*),P(*)
      DIMENSION F1(*),F2(*),FP(*)
      DIMENSION A(*),KACOL(*),KALD(*),B1(*),B2(*),KBCOL(*),KBLD(*)
      DIMENSION KMBD(*),KVERT(NNVE,*),KMID(NNVE,*),KNPR(*)
C
C *** Local arrays for informations about one element
      DIMENSION AA1(4,4),AA2(4,4),BB1(4),BB2(4),FF1(4),FF2(4)
      DIMENSION IU(4),BDBC(4)
C
      NA=KALD(NU+1)-1
C
C=======================================================================
C     Block Gauss-Seidel on Schur Complement
C=======================================================================
C *** Loop over all elements
C
      DO 1  IELN=1,NEL
C
c      IF (ISORTP.GT.0) THEN
c         write (*,*) 'in VANCAN geht was schief', isortp
c      ELSE
        IEL=IELN
c      ENDIF
C
      FFP=FP(IEL)
C
C-----------------------------------------------------------------------
C *** Loop over all 4 U-nodes of that element
      DO 11  II=1,4
C
      IMID=KMID(II,IEL)
      I=IMID-NVT
      IU(II)=I
      IF (KNPR(IMID).EQ.0)  THEN
       BDBC(II)=.FALSE.
      ELSE
       BDBC(II)=.TRUE.
      ENDIF
C
C *** Put on AA(.) the diagonal entry of matrix A
      IA1=KALD(I)
      AA1(II,II)=DBLE(A(IA1))
      AA2(II,II)=DBLE(A(3*NA+IA1))
C
C *** Initial setting of FF1,FF2
      FF1(II)=F1(I)
      FF2(II)=F2(I)
C
      IF (BDBC(II)) GOTO 13
C
C *** Put all off-diagonal entries of A to the right hand side
ccc      IA2=KALD(I+1)-1
ccc      DO 110  IA=IA1+1,IA2
ccc      J=KACOL(IA)
ccc      AOFF=DBLE(A(IA))
ccc      FF1(II)=FF1(II)-AOFF*U1(J)
ccc      FF2(II)=FF2(II)-AOFF*U2(J)
ccc110   CONTINUE
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
          FF1(II)=FF1(II)-DBLE(B1(IB2))*PJP
          FF2(II)=FF2(II)-DBLE(B2(IB2))*PJP
      ELSE IF (JP2.EQ.IEL)  THEN
          PJP=P(JP1)
          BB1(II)=DBLE(B1(IB2))
          BB2(II)=DBLE(B2(IB2))
          FF1(II)=FF1(II)-DBLE(B1(IB1))*PJP
          FF2(II)=FF2(II)-DBLE(B2(IB1))*PJP
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
      AA1(II,JJ)=0D0
      AA2(II,JJ)=0D0
      IF (BDBC(II)) GOTO 22
C
      J=IU(JJ)
      IA1=KALD(I)
      IA2=KALD(I+1)-1
C
      DO 220  IA=IA1+1,IA2
      JH=KACOL(IA)
      IF (J.EQ.JH) THEN
       AA1(II,JJ)=DBLE(A(IA))
       AA2(II,JJ)=DBLE(A(3*NA+IA))
       GOTO 22
      ENDIF
 220  CONTINUE
C
 22   CONTINUE
C
 21   CONTINUE
C
C-----------------------------------------------------------------------
C *** Calling the subroutine ELUPDT for element update 
      CALL  ELUPDN (U1,U2,P,IEL,IU,BDBC,AA1,AA2,BB1,BB2,FF1,FF2,FFP)
C
1     CONTINUE
C
C=======================================================================
C
C
C
      END
c

************************************************************************
      SUBROUTINE  VANCAS (U1,U2,P,U1OLD,U2OLD,POLD,F1,F2,FP,
     *                    A,KCOLA,KLDA,B1,B2,KCOLB,KLDB,NU,NP,
     *                    KMBD,KVERT,KMID,KNPR,NMBD)
************************************************************************
C Purpose:  Vanca-type Glaetter. 
C     Auf jedem Element wird von der lokalen Matrix
c     (AA(i,j)) nur die Diagonale invertiert
************************************************************************
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DOUBLE PRECISION BB1,BB2
      REAL  A,B1,B2
C
      PARAMETER (NNVE=4)
      PARAMETER (NNLEV=9)
      DIMENSION U1(*),U2(*),P(*),U1OLD(*),U2OLD(*),POLD(*)
      DIMENSION F1(*),F2(*),FP(*)
      DIMENSION A(*),KCOLA(*),KLDA(*),B1(*),B2(*),KCOLB(*),KLDB(*)
      DIMENSION KMBD(*),KVERT(NNVE,*),KMID(NNVE,*),KNPR(*)
C
C *** Local arrays for informations about one element
      DIMENSION AA1(4),AA2(4),BB1(4),BB2(4),FF1(4),FF2(4)
      DIMENSION IU(4),BDBC(4)
C
C *** Usual data for mesh management
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
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
      INCLUDE 'bouss.inc'
      include 'block.inc'
      SAVE 
C
      NA=KLDA(NU+1)-1
C-----------------------------------------------------------------------
C *** Update  UOLD:=U
      CALL  LCP1 (U1,U1OLD,NU)
      CALL  LCP1 (U2,U2OLD,NU)
      CALL  LCP1 (P ,POLD ,NP)
C
C=======================================================================
C     Block Gauss-Seidel on Schur Complement
C=======================================================================
C *** Loop over all elements
C
      DO 1  IEL=1,NEL
      FFP=FP(IEL)
C
C-----------------------------------------------------------------------
C *** Loop over all 4 U-nodes of that element
      DO 11  II=1,4
C
      IMID=KMID(II,IEL)
      I=IMID-NVT
      IU(II)=I
      IF (KNPR(IMID).EQ.0)  THEN
       BDBC(II)=.FALSE.
      ELSE
       BDBC(II)=.TRUE.
      ENDIF
C
C *** Put on AA(.) the diagonal entry of matrix A
      IA1=KLDA(I)
      AA1(II)=DBLE(A(IA1))
      AA2(II)=DBLE(A(3*NA+IA1))
C
C *** Initial setting of FF1,FF2
      FF1(II)=F1(I)
      FF2(II)=F2(I)
C
      IF (BDBC(II)) GOTO 13
C
C *** Put all off-diagonal entries of A to the right hand side
      IA2=KLDA(I+1)-1
      DO 110  IA=IA1+1,IA2
      J=KCOLA(IA)
      AOFF1=DBLE(A(IA))
      AOFF2=DBLE(A(3*NA+IA))
      FF1(II)=FF1(II)-AOFF1*U1(J)
      FF2(II)=FF2(II)-AOFF2*U2(J)
 110  CONTINUE
      DO 111  IA=IA1,IA2
      J=KCOLA(IA)
      AOFF1=DBLE(A(NA+IA))
      AOFF2=DBLE(A(2*NA+IA))
      FF1(II)=FF1(II)-AOFF1*U2(J)
      FF2(II)=FF2(II)-AOFF2*U1(J)
 111  CONTINUE
C
C-----------------------------------------------------------------------
C *** Get BB1,BB2 and modify FF1,FF2
      IB1=KLDB(I)
      IB2=KLDB(I+1)-1
      JP1=KCOLB(IB1)
      JP2=KCOLB(IB2)
C
      IF (JP1.EQ.IEL)  THEN
          PJP=P(JP2)
          BB1(II)=DBLE(B1(IB1))
          BB2(II)=DBLE(B2(IB1))
          FF1(II)=FF1(II)-DBLE(B1(IB2))*PJP
          FF2(II)=FF2(II)-DBLE(B2(IB2))*PJP
      ELSE IF (JP2.EQ.IEL)  THEN
          PJP=P(JP1)
          BB1(II)=DBLE(B1(IB2))
          BB2(II)=DBLE(B2(IB2))
          FF1(II)=FF1(II)-DBLE(B1(IB1))*PJP
          FF2(II)=FF2(II)-DBLE(B2(IB1))*PJP
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
C
11    CONTINUE
C
C-----------------------------------------------------------------------
C *** Calling the subroutine ELUPDT for element update 
      CALL  ELUPDT (U1,U2,P,IEL,IU,BDBC,AA1,AA2,BB1,BB2,FF1,FF2,FFP)
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
       CALL  LLC1 (POLD ,P ,NP,A1,A2)
      ENDIF
C
C
C
      END
C
C
************************************************************************
      SUBROUTINE  VANCAM (U1,U2,P,F1,F2,FP,
     *                    A,KACOL,KALD,B1,B2,KBCOL,KBLD,
     *                    KMBD,KVERT,KMID,KNPR,NMBD,Inum,AKTELE)
************************************************************************
C
C-----------------------------------------------------------------------
C Purpose:   So ähnlich wie VANCAS, nur auf Defektbasis
C
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'bouss.inc'
      INCLUDE 'block.inc'
C-----------------------------------------------------------------------
C
      DOUBLE PRECISION BB1,BB2
      REAL  A,B1,B2
      INTEGER AKTELE
C
c      PARAMETER (NNVE=4)
      DIMENSION U1(*),U2(*),P(*)
      DIMENSION F1(*),F2(*),FP(*)
      DIMENSION A(*),KACOL(*),KALD(*),B1(*),B2(*),KBCOL(*),KBLD(*)
      DIMENSION KMBD(*),KVERT(NNVE,*),KMID(NNVE,*),KNPR(*)
C
C *** Local arrays for informations about one element
      DIMENSION AA1(4),AA2(4),BB1(4),BB2(4),FF1(4),FF2(4)
      DIMENSION IU(4),BDBC(4)
C
C
      NA=KALD(NU+1)-1

C=======================================================================
C     Block Gauss-Seidel on Schur Complement
C=======================================================================
C *** Loop over all elements (or only one)
C
      IF (INUM.eq.0) THEN
         Nloop=1
      ELSE
         NLOOP=NEL
      ENDIF
      CALL VKADD (AVK1,AVK2)

      DO 1  IELN=1,NLOOP!NELC

      IF (INUM.eq.0) THEN
       IEL=AKTELE
      ELSE
       IEL=IELN
      ENDIF
C
      FFP=FP(IEL)
C
C-----------------------------------------------------------------------
C *** Loop over all 4 U-nodes of that element
      DO 11  II=1,4
C
      IMID=KMID(II,IEL)
      I=IMID-NVT
      IU(II)=I
      IF (KNPR(IMID).EQ.0)  THEN
       BDBC(II)=.FALSE.
      ELSE
       BDBC(II)=.TRUE.
      ENDIF
C
C *** Put on AA(.) the diagonal entry of matrix A
      IA1=KALD(I)
      AA1(II)=AVK1*DBLE(A(IA1))+AVK2*DBLE(A(4*NA+ia1))
      AA2(II)=AVK1*DBLE(A(3*NA+IA1))+AVK2*DBLE(A(4*NA+ia1))
C
C *** Initial setting of FF1,FF2
      FF1(II)=F1(I)
      FF2(II)=F2(I)
C
      IF (BDBC(II)) GOTO 13
C
C *** Put all off-diagonal entries of A to the right hand side
      IA2=KaLD(I+1)-1
      DO 110  IA=IA1+1,IA2
      J=KaCOL(IA)
      AOFF1=AVK1*DBLE(A(IA))+AVK2*DBLE(A(4*NA+ia))
      AOFF2=AVK1*DBLE(A(3*NA+IA))+AVK2*DBLE(A(4*NA+ia))
      FF1(II)=FF1(II)-AOFF1*U1(J)
      FF2(II)=FF2(II)-AOFF2*U2(J)
 110  CONTINUE
      DO 111  IA=IA1,IA2
      J=KaCOL(IA)
      AOFF1=AVK1*DBLE(A(NA+IA))
      AOFF2=AVK1*DBLE(A(2*NA+IA))
      FF1(II)=FF1(II)-AOFF1*U2(J)
      FF2(II)=FF2(II)-AOFF2*U1(J)
 111  CONTINUE

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
          FF1(II)=FF1(II)-DBLE(B1(IB2))*PJP
          FF2(II)=FF2(II)-DBLE(B2(IB2))*PJP
      ELSE IF (JP2.EQ.IEL)  THEN
          PJP=P(JP1)
          BB1(II)=DBLE(B1(IB2))
          BB2(II)=DBLE(B2(IB2))
          FF1(II)=FF1(II)-DBLE(B1(IB1))*PJP
          FF2(II)=FF2(II)-DBLE(B2(IB1))*PJP
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
C
11    CONTINUE
C
C-----------------------------------------------------------------------
C *** Calling the subroutine ELUPDT for element update 
      CALL  ELUPDT (U1,U2,P,IEL,IU,BDBC,AA1,AA2,BB1,BB2,
     *  FF1,FF2,FFP)
C
1     CONTINUE
C
C=======================================================================
C
C
C
      END
C
C
************************************************************************
      SUBROUTINE  ELUPDT  (U1,U2,P,IEL,IU,BDBC,AA1,AA2,BB1,BB2,FF1,FF2,
     *                    FFP)
************************************************************************
*    Purpose: - Executes block Gauss-Seidel update on (U1,U2,P)
*               for all unknowns of the element IEL
************************************************************************
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DOUBLE PRECISION  BB1,BB2
C
      DIMENSION U1(*),U2(*),P(*)
C
C *** Local arrays for informations about the element
      DIMENSION AA1(4),AA2(4),BB1(4),BB2(4),FF1(4),FF2(4),IU(4),BDBC(4)
      DIMENSION AI1(4),AI2(4),DD1(4),DD2(4),UU1(4),UU2(4)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      SAVE 
C
C-----------------------------------------------------------------------
C
      DO 1  II=1,4
      IF (BDBC(II))  THEN
          DD1(II)=0D0
          DD2(II)=0D0
      ELSE
          DD1(II)=BB1(II)
          DD2(II)=BB2(II)
      ENDIF
      IF (DABS(AA1(II)).LT.1.D-10.or.DABS(AA2(II)).LT.1.D-10)  THEN
         WRITE(MTERM,*)'ERROR in ELUPDT: diagonal entry is nearly zero'
         WRITE(MTERM,*)'II, I, BNDR(II) : ',II,IU(II),BDBC(II)
         RETURN
      ENDIF
      AI1(II)=1D0/AA1(II)
      AI2(II)=1D0/AA2(II)
1     CONTINUE
C
C *** Factorization loop
C
      DP=0.D0
      DO 2  II=1,4
      DP =DP  - AI1(II)*BB1(II)*DD1(II)-AI2(II)*BB2(II)*DD2(II)
      FFP=FFP - AI1(II)*BB1(II)*FF1(II)-AI2(II)*BB2(II)*FF2(II)
2     CONTINUE
C
C *** Solution loop
C
      IF (DABS(DP).LT.1.D-10)  THEN
         WRITE(MTERM,*)'ERROR in ELUPDT: DP is nearly zero'
         RETURN
      ENDIF
      PP=FFP/DP
      DO 3  II=1,4
      UU1(II)=AI1(II)*(FF1(II)-DD1(II)*PP)
      UU2(II)=AI2(II)*(FF2(II)-DD2(II)*PP)
3     CONTINUE
C
C *** Updating the global solution vector (U1,U2,P)
C
      P(IEL)=PP
      DO 4  II=1,4
      IF (BDBC(II)) GOTO 4
      I=IU(II)
      U1(I)=UU1(II)
      U2(I)=UU2(II)
4     CONTINUE
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE  VANCAE (U1,U2,P,U1OLD,U2OLD,POLD,F1,F2,FP,
     *                    A,KCOLA,KLDA,B1,B2,KCOLB,KLDB,NU,NP,
     *                    KMBD,KVERT,KMID,KNPR,NMBD,RLXPAR,DMAXU,DMAXP)
************************************************************************
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DOUBLE PRECISION BB1,BB2
      REAL  A,B1,B2
C
      PARAMETER (NNVE=4)
      PARAMETER (NNLEV=9)
      DIMENSION U1(*),U2(*),P(*),U1OLD(*),U2OLD(*),POLD(*)
      DIMENSION F1(*),F2(*),FP(*)
      DIMENSION A(*),KCOLA(*),KLDA(*),B1(*),B2(*),KCOLB(*),KLDB(*)
      DIMENSION KMBD(*),KVERT(NNVE,*),KMID(NNVE,*),KNPR(*)
C
C *** Local arrays for informations about one element
      DIMENSION AA1(4),AA2(4),BB1(4),BB2(4),FF1(4),FF2(4),IU(4),BDBC(4)
C
C *** Usual data for mesh management
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
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
      INCLUDE 'bouss.inc'
      INCLUDE 'block.inc'
      SAVE 
C
C-----------------------------------------------------------------------
      NA=KLDA(NU+1)-1
      DMAXU=0D0
      DMAXP=0D0
C
C *** Update  UOLD:=U
      CALL  LCP1 (U1,U1OLD,NU)
      CALL  LCP1 (U2,U2OLD,NU)
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
C *** Loop over all 4 U-nodes of that element
      DO 11  II=1,4
C
      IMID=KMID(II,IEL)
      I=IMID-NVT
      IU(II)=I
      IF (KNPR(IMID).EQ.0)  THEN
         BDBC(II)=.FALSE.
      ELSE
         BDBC(II)=.TRUE.
      ENDIF
C
C *** Put on AA(.) the diagonal entry of matrix A
      IA1=KLDA(I)
      AA1(II)=DBLE(A(IA1))
      AA2(II)=DBLE(A(3*NA+IA1))
C
C *** Initial setting of FF1,FF2
      FF1(II)=F1(I)
      FF2(II)=F2(I)
C
      IF (BDBC(II)) GOTO 13
C
C *** Put all off-diagonal entries of A to the right hand side
      IA2=KLDA(I+1)-1
      DO 110  IA=IA1+1,IA2
      J=KCOLA(IA)
      AOFF1=DBLE(A(IA))
      AOFF2=DBLE(A(3*NA+IA))
      FF1(II)=FF1(II)-AOFF1*U1(J)
      FF2(II)=FF2(II)-AOFF2*U2(J)
 110  CONTINUE
      DO 111  IA=IA1,IA2
      J=KCOLA(IA)
      AOFF1=DBLE(A(NA+IA))
      AOFF2=DBLE(A(2*NA+IA))
      FF1(II)=FF1(II)-AOFF1*U2(J)
      FF2(II)=FF2(II)-AOFF2*U1(J)
 111  CONTINUE
c
C-----------------------------------------------------------------------
C *** Get BB1,BB2 and modify FF1,FF2
      IB1=KLDB(I)
      IB2=KLDB(I+1)-1
      JP1=KCOLB(IB1)
      JP2=KCOLB(IB2)
C
      IF (JP1.EQ.IEL)  THEN
          PJP=P(JP2)
          BB1(II)=DBLE(B1(IB1))
          BB2(II)=DBLE(B2(IB1))
          FF1(II)=FF1(II)-DBLE(B1(IB2))*PJP
          FF2(II)=FF2(II)-DBLE(B2(IB2))*PJP
          IB3    =IB1
          JP3    =JP1
      ELSE IF (JP2.EQ.IEL)  THEN
          PJP=P(JP1)
          BB1(II)=DBLE(B1(IB2))
          BB2(II)=DBLE(B2(IB2))
          FF1(II)=FF1(II)-DBLE(B1(IB1))*PJP
          FF2(II)=FF2(II)-DBLE(B2(IB1))*PJP
          IB3    =IB2
          JP3    =JP2
      ELSE
          WRITE(MTERM,*) 'ERROR in SMOOTH: IEL entry in B not found'
          STOP
      ENDIF
C
      AH1=DBLE(A(IA1))
      AH2=DBLE(A(3*NA+IA1))
      ICOLA1=KCOLA(IA1)
      PJP=P(JP3)
      DMAXU=MAX(DMAXU,ABS(FF1(II)-AH1*U1(ICOLA1)-DBLE(B1(IB3))*PJP))
      DMAXU=MAX(DMAXU,ABS(FF2(II)-AH2*U2(ICOLA1)-DBLE(B2(IB3))*PJP))
      FFPH=FFPH-BB1(II)*U1(IU(II))-BB2(II)*U2(IU(II))
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
C
      FFPH=FFPH-BB1(II)*U1(IU(II))-BB2(II)*U2(IU(II))
C
11    CONTINUE
C
      DMAXP=MAX(DMAXP,FFPH)
C
C-----------------------------------------------------------------------
C *** Calling the subroutine ELUPDT for element update 
      CALL  ELUPDT (U1,U2,P,IEL,IU,BDBC,AA1,AA2,BB1,BB2,FF1,FF2,FFP)
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
      SUBROUTINE  ELUPDN  (U1,U2,P,IEL,IU,BDBC,AA1,AA2,
     *                    BB1,BB2,FF1,FF2,FFP)
************************************************************************
C
C-----------------------------------------------------------------------
*    Purpose: - Executes block Gauss-Seidel update on (U1,U2,P)
*               for all unknowns of the element IEL
C
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'block.inc'
      INCLUDE 'dwork.inc'
C-----------------------------------------------------------------------
C
C
      DOUBLE PRECISION  BB1,BB2
      INTEGER KPIVOT
C
      DIMENSION U1(*),U2(*),P(*)
C
C *** Local arrays for informations about the element
      DIMENSION AA1(4,4),AA2(4,4),BB1(4),BB2(4)
      DIMENSION FF1(4),FF2(4),IU(4),BDBC(4)
      DIMENSION AI1(4),AI2(4),DD1(4),DD2(4),UU1(4),UU2(4)
      DIMENSION DX1(4),DX2(4),DR1(4),DR2(4)
      dimension KPIVO1 (4),KPIVO2 (4)
C
C
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
c$$$      IF (iel.eq.3) THEN
c$$$         DO 79 i=1,4
c$$$               write (*,*) ff1(i), ff2(i)! (aa1(i,j), j=1,4)
c$$$ 79             continue
c$$$       write (*,*) '======'
c$$$      ENDIF
c      return
C-----------------------------------------------------------------------
c$$$      call ztime (t1)! WEGMACHEN!!!
      DO 1  II=1,4
      IF (BDBC(II))  THEN
          DD1(II)=0D0
          DD2(II)=0D0
      ELSE
          DD1(II)=BB1(II)
          DD2(II)=BB2(II)
      ENDIF
      IF (DABS(AA1(II,II)).LT.1.D-10.or.DABS(AA2(II,II)).LT.1.D-10) THEN
         WRITE(MTERM,*)'ERROR in ELUPDN: diagonal entry is nearly zero'
         WRITE(MTERM,*)'II, I, BNDR(II) : ',II,IU(II),BDBC(II)
         RETURN
      ENDIF
1     CONTINUE
C
C *** Factorization loop
C
      CALL INVERT1(AA1,DD1,DD1,KPIVO1,0)
      CALL INVERT1(AA2,DD1,DD1,KPIVO2,0)
C
      DP=0.D0
      CALL INVERT1(AA1,DX1,DD1,KPIVO1,1)
      CALL INVERT1(AA2,DX2,DD2,KPIVO2,1)
      DO 20  II=1,4
      DP =DP -(BB1(II)*DX1(II)+BB2(II)*DX2(II))
 20   CONTINUE
C
      CALL INVERT1(AA1,DX1,FF1,KPIVO1,1)
      CALL INVERT1(AA2,DX2,FF2,KPIVO2,1)
      DO 30  II=1,4
      FFP =FFP -(BB1(II)*DX1(II)+BB2(II)*DX2(II))
 30   CONTINUE
C
C *** Solution loop
C
      IF (DABS(DP).LT.1.D-10)  THEN
         WRITE(MTERM,*)'ERROR in ELUPDN: DP is nearly zero'
         RETURN
      ENDIF
      PP=FFP/DP
      DO 3  II=1,4
      DR1(II)=FF1(II)-DD1(II)*PP
      DR2(II)=FF2(II)-DD2(II)*PP
3     CONTINUE
C
      CALL INVERT1(AA1,UU1,DR1,KPIVO1,1)
      CALL INVERT1(AA2,UU2,DR2,KPIVO2,1)
C
C *** Updating the global solution vector (U1,U2,P)
C
      P(IEL)=PP
      DO 4  II=1,4
      IF (BDBC(II)) GOTO 4
      I=IU(II)
      U1(I)=UU1(II)
      U2(I)=UU2(II)
4     CONTINUE
C
C
C
c$$$      call ztime (t2)! WEGMACHEN!!!
c$$$      write (*,*) 'Zeit für elupdn', t2-t1
      END

************************************************************
************************************************************************
      SUBROUTINE  ELUPN2  (U1,U2,P,IEL,IU,BDBC,AA,
c$$$      SUBROUTINE  ELUPN2  (UV,P,IEL,IU,BDBC,AA,
     *                    BB1,BB2,FF,FFP)
************************************************************************
C
C-----------------------------------------------------------------------
*    Purpose: - Executes block Gauss-Seidel update on (U1,U2,P)
*               for all unknowns of the element IEL
C               U1,U2 are updated SIMULTANUOUSLY
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'block.inc'
      INCLUDE 'dwork.inc'
C-----------------------------------------------------------------------
C
C
      DOUBLE PRECISION  BB1,BB2
      INTEGER KPIVO
C
      DIMENSION U1(*),U2(*),P(*)
c      DIMENSION P(*)

C
C *** Local arrays for informations about the element
      DIMENSION AA(8,8),BB1(4),BB2(4)
      DIMENSION FF(8),IU(4),BDBC(4)
      DIMENSION AI(8),DD(8),UV(8)
      DIMENSION DX(8),DR(8)
      dimension KPIVO (8)
      SAVE
C
C-----------------------------------------------------------------------
C
c$$$C-----------------------------------------------------------------------
c$$$      IF (iel.eq.3) THEN
c$$$         DO 79 i=1,8
c$$$               write (*,*) ff(i)!(real(aa(i,j)), j=1,8), '|'
c$$$ 79             continue
c$$$       write (*,*) '======'
c$$$      ENDIF
c$$$c      return
C-----------------------------------------------------------------------
c
      DO 1  II=1,4
      IF (BDBC(II))  THEN
          DD(II)=0D0
          DD(4+II)=0D0
      ELSE
          DD(II)=BB1(II)
          DD(4+II)=BB2(II)
      ENDIF
      IF (DABS(AA(II,II)).LT.1.D-10.or.DABS(AA(4+II,4+II)).LT.1.D-10)
     *   THEN
         WRITE(MTERM,*)'ERROR in ELUPN2: diagonal entry is nearly zero'
         WRITE(MTERM,*)'II, I, BNDR(II) : ',II,IU(II),BDBC(II)
         RETURN
      ENDIF
1     CONTINUE
C
C *** Factorization loop
C
      CALL INVER8(AA,DD,DD,KPIVO,0)
C
      DP=0.D0
      CALL INVER8(AA,DX,DD,KPIVO,1)
      DO 20  II=1,4
      DP =DP -(BB1(II)*DX(II)+BB2(II)*DX(4+II))
 20   CONTINUE
C
      CALL INVER8(AA,DX,FF,KPIVO,1)
      DO 30  II=1,4
      FFP =FFP -(BB1(II)*DX(II)+BB2(II)*DX(4+II))
 30   CONTINUE
C
C *** Solution loop
C
      IF (DABS(DP).LT.1.D-10)  THEN
         WRITE(MTERM,*)'ERROR in ELUPDN: DP is nearly zero'
         RETURN
      ENDIF
      PP=FFP/DP
      DO 3  II=1,8
      DR(II)=FF(II)-DD(II)*PP
3     CONTINUE
C
      CALL INVER8(AA,UV,DR,KPIVO,1)
C

C *** Updating the global solution vector (U1,U2,P)
C    Lasse ich jetzt mal, uebergebe dafür UV!!
      P(IEL)=PP
      DO 4  II=1,4
      IF (BDBC(II)) GOTO 4
      I=IU(II)
      U1(I)=UV(II)
      U2(I)=UV(4+II)
4     CONTINUE
C
C

c$$$      call ztime (t2)! WEGMACHEN!!!
c$$$      write (*,*) 'Zeit für elupdn', t2-t1
      END

**************************************************************
C
C
C
      SUBROUTINE INVERT1(A,X,F,KPIVOT,IPAR)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DOUBLE PRECISION A,B,F,X,work,det
      INTEGER KPIVOT, lwor, iinver
      PARAMETER (NDIM=4,NNLEV=9)
      COMMON /SBLOCK/ PSBARM,PSBARP,PSBARR,PSBARB,PSBARS,PSBVOL,
     *                KLLDPA(NNLEV),KLPAT(NNLEV),KLPAEL(NNLEV),
     *                KNPTCH(NNLEV),KLELC1(NNLEV),KLELC2(NNLEV),
     *                NSBLSM,NSBLSL,NELMAC,IBLOCK,BLOCSL,
     *                IFACTO,IINVER
c      COMMON /SBLOCK/ PSBARM,PSBARP,PSBARR,PSBARB,PSBVOL,
c     *                KLLDPA(NNLEV),KLPAT(NNLEV),KLPAEL(NNLEV),
c     *                KNPTCH(NNLEV),KLELC1(NNLEV),KLELC2(NNLEV),
c     *                NSBLSM,NSBLSL,NELMAC,IBLOCK,BLOCSL,
c     *                IFACTO,IINVER
      DIMENSION A(NDIM,NDIM),B(NDIM,NDIM),F(NDIM),X(NDIM),
     *          MERKX(NDIM),MERKY(NDIM),KPIVOT(NDIM),
     *          work(NDIM)
C
C
c
      IF (IINVER.GE.0) THEN
         IF (IPAR.EQ.0) THEN
          CALL AUSTAU1(NDIM,NDIM,A,B,MERKX,MERKY,IFEHL)
          DO 10 IA=1,NDIM
          DO 10 IB=1,NDIM
10        A(IA,IB)=B(IA,IB)
         ENDIF
C
         IF (IPAR.EQ.1) THEN
          DO 20 IA=1,NDIM
          X(IA)=0D0
          DO 22 IB=1,NDIM
          X(IA)=X(IA)+A(IA,IB)*F(IB)
22        CONTINUE
20        CONTINUE
         ENDIF
c
      ELSE
         IF (IPAR.EQ.0) THEN        
         CALL DGETRF (Ndim,Ndim,a,Ndim,kPIVot,IFEHL)
         ENDIF
C
         IF (IPAR.EQ.1) THEN
         do 300 i=1,ndim
 300        work(i)=f(i)
          CALL DGETRS ('N',NDIM,1,A,NDIM,KPIVOT,work,NDIM,INFO)
          do 310 i=1,ndim
 310        x(i)=work(i)
         ENDIF
      ENDIF
C

C
      END
C
C
      SUBROUTINE AUSTAU1(NDIM,N,A,B,MERKX,MERKY,IFEHL)
      DOUBLE PRECISION A,B,HILF,PIVOT
      DIMENSION A(NDIM,N),B(NDIM,N),MERKX(N),MERKY(N)
C
      IFEHL=1
      DO 100 I=1,N
      MERKX(I)=0
      MERKY(I)=0!O
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
**************************************************************
C
C
C
      SUBROUTINE INVER8(A,X,F,KPIVOT,IPAR)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DOUBLE PRECISION A,B,F,X,work,det
      INTEGER KPIVOT, lwor, iinver
      PARAMETER (NDIM=8,NNLEV=9)
      COMMON /SBLOCK/ PSBARM,PSBARP,PSBARR,PSBARB,PSBARS,PSBVOL,
     *                KLLDPA(NNLEV),KLPAT(NNLEV),KLPAEL(NNLEV),
     *                KNPTCH(NNLEV),KLELC1(NNLEV),KLELC2(NNLEV),
     *                NSBLSM,NSBLSL,NELMAC,IBLOCK,BLOCSL,
     *                IFACTO,IINVER
c      COMMON /SBLOCK/ PSBARM,PSBARP,PSBARR,PSBARB,PSBVOL,
c     *                KLLDPA(NNLEV),KLPAT(NNLEV),KLPAEL(NNLEV),
c     *                KNPTCH(NNLEV),KLELC1(NNLEV),KLELC2(NNLEV),
c     *                NSBLSM,NSBLSL,NELMAC,IBLOCK,BLOCSL,
c     *                IFACTO,IINVER
      DIMENSION A(NDIM,NDIM),B(NDIM,NDIM),F(NDIM),X(NDIM),
     *          MERKX(NDIM),MERKY(NDIM),KPIVOT(NDIM),
     *          work(NDIM)
C
C
c
      IF (IINVER.GE.0) THEN
         IF (IPAR.EQ.0) THEN
          CALL AUSTAU1(NDIM,NDIM,A,B,MERKX,MERKY,IFEHL)
          DO 10 IA=1,NDIM
          DO 10 IB=1,NDIM
10        A(IA,IB)=B(IA,IB)
         ENDIF
C
         IF (IPAR.EQ.1) THEN
          DO 20 IA=1,NDIM
          X(IA)=0D0
          DO 22 IB=1,NDIM
          X(IA)=X(IA)+A(IA,IB)*F(IB)
22        CONTINUE
20        CONTINUE
         ENDIF
c
      ELSE
         IF (IPAR.EQ.0) THEN        
         CALL DGETRF (Ndim,Ndim,a,Ndim,kPIVot,IFEHL)
         ENDIF
C
         IF (IPAR.EQ.1) THEN
         do 300 i=1,ndim
 300        work(i)=f(i)
          CALL DGETRS ('N',NDIM,1,A,NDIM,KPIVOT,work,NDIM,INFO)
          do 310 i=1,ndim
 310        x(i)=work(i)
         ENDIF
      ENDIF
C

C
      END
C
C

c========================================================================
C      Filterroutine, die tangentiale Kreisgeschwindigkeiten rausfiltert
c========================================================================
c
      SUBROUTINE GARLIC (IEL,ILEV,nvt, DCORVG, KVERT,KMID,
     *                   U1,U2, Dul2,dum,dalfa,dZBV,
     *                   DFILT1,DFILT2,ICHECK)
c
      PARAMETER (NNVE=4)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
c      
      DIMENSION KVERT(NNVE,*), DCORVG(2,*),KMID(NNVE,*)
      DIMENSION XV(4),YV(4), TX(4),TY(4),BRAN(4)
      DIMENSION XE(4),YE(4),UV(8),U1(*),U2(*),IU(4),DFILT1(*),DFILT2(*)
      SAVE
C-----------------------------------------------------------------------
C *** Hier will ich den Filter austesten.
C
c     Erstmal die Tangenten berechnen:
c
      dul2=0d0
      dum=0d0
      dalfa=0d0
      dquot=0d0
      do 77 i=1,8
 77   UV(i)=0d0
        S1=0d0
        S2=0d0
      XM=0d0
      YM=0d0   !Ergibt Elementmittelpunkt
c
      DO 10  II=1,4
      IMID=KMID(II,IEL)
      I=IMID-NVT
      IU(II)=I
c
      IV=KVERT(II,IEL)
      XV(II)=DCORVG(1,IV)
      YV(II)=DCORVG(2,IV)
      XM=XM+XV(II)
      YM=YM+YV(II)
 10   CONTINUE
C
      XM=XM/4D0
      YM=YM/4D0
        XM=2D0**(ILEV)*XM
        YM=2D0**(ILEV)*YM
        IF (MOD(INT(XM+YM),4).eq.0)THEN
           SIGNUM=1d0
        ELSE
           SIGNUM=1d0
        ENDIF

      DO 20 II=1,4
         II1=II
         II2=II+1
         IF (II2.GT.4) II2=1
         TX(II1) =(XV(II1)-XV(II2))
         TY(II1) =(YV(II1)-YV(II2))
         XE(II) =0.5d0*(XV(II)+XV(II2))
         YE(II) =0.5d0*(YV(II)+YV(II2))
         IF ((ABS(MIN(XE(II),YE(II))).LE.1D-9).OR.
     *      ((ABS(1D0-MAX(YE(II),XE(II))).LE.1D-9))) THEN
           BRAN(ii)=.TRUE.
c        if (ABS(1D0-(XE(II))).LE.1D-9) write (*,*) iu(ii)
           ELSE
           BRAN(ii)=.false.
           ENDIF
         TN = SQRT(Tx(ii1)**2+Ty(ii1)**2)
         TX(ii1)=Tx(ii1)/TN*signum
         TY(II1)=TY(II1)/TN*signum
20      CONTINUE

C     Check, welche  Orientierung die Tangenten haben sollen
C     Man geht hier von einem regelmaessig verfeinerten EQ aus
c


        DO 30 J=1,4
           i=iu(j)
           S1=S1+ U1(i)*Tx(j)+ U2(i)*ty(j)
           S2=S2+ tx(j)*tx(j)+ ty(j)*ty(j)
 30     CONTINUE
c
           ALPHA=-S1/S2

           IF (ABS (ALPHA).ge.AMAX) AMAX=ABS(ALPHA)
c
         
        DO 100 I=1,4
           UV(  i)=UV(  i)+ALPHA*tx(i)!*xhelp!*0d0
           UV(4+i)=UV(4+i)+ALPHA*Ty(i)!*xhelp!*0d0
 100       CONTINUE
           ul2=0d0
           um=0d0
           uneul2=0d0
           bgo=.false.
           DO 444 JJ=1,4
              if (bran(jj))  bgo=.true.
 444          continue
c              if (bgo) goto 505
      DO 400  II=1,4
      I=IU(II)
      ul2=ul2+u1(i)**2+u2(i)**2
      um=um+U1(i)+U2(i)
c
      xhelp=0.5d0
      IF (BRAN(II))then 
         xhelp=0d0
         endif
      IF (ICHECK.EQ.0) THEN
      U1(I)=U1(i)+UV(II)
      U2(I)=U2(i)+UV(4+II)
      uneul2=uneul2+u1(i)**2+u2(i)**2
      ELSE

      DFILT1(i)=DFILT1(i)+UV(ii)*xhelp
      DFILT2(i)=DFILT2(i)+UV(II+4)*xhelp
      uneul2=uneul2+DFILT1(i)**2+DFILT2(i)**2
      ENDIF
c

c      if (i.eq.12) write (*,*) u1(i), UV(II)
 400  CONTINUE
 505  dul2=SQRT(UL2)
      dum=um
      dalfa=(alpha)
      dZBV=SQRT(uneul2)
        END


c========================================================================
C      Filterroutine, die tangentiale Kreisgeschwindigkeiten rausfiltert
c========================================================================
c
      SUBROUTINE GARLI2 (ILEV,NEL,NVT, DCORVG,KVERT,KMID,KNPR,U1,U2,
     *    D1,D2,ite,nlmax)
c
      PARAMETER (NNVE=4,NNLEV=9)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299,NNWORK=1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
c      
      DIMENSION KVERT(NNVE,*), DCORVG(2,*),KNPR(*),KMID(NNVE,*)
      DIMENSION XV(4),YV(4), TX(4),TY(4), UU(4),UV(4)
      DIMENSION XE(4),YE(4),U1(*),U2(*),IU(4),BDBC(4)
      DIMENSION D1(*),D2(*),DU(4),DV(4),BRAN(4)
c      EXTERNAL f_u,f_v
      INCLUDE 'block.inc'
      SAVE
C-----------------------------------------------------------------------
C *** Hier will ich den Filter austesten.
C
c     Erstmal die Tangenten berechnen:
c
c      return
c
c      write (*,*) 'ponn'
      S1U=0d0
      S1D=0d0
      S2=0d0
c
      DO 1 IEL=1,NEL
C *** Loop over all 4 U-nodes of that element


      XM=0d0
      YM=0d0   !Ergibt Elementmittelpunkt
c
      DO 10  II=1,4
      IMID=KMID(II,IEL)
      I=IMID-NVT
      IU(II)=I
      UU(II)=U1(i)
      UV(II)=U2(i)
      DU(II)=D1(i)
      DV(II)=D2(i)
c
      IV=KVERT(II,IEL)
      IV=KVERT(II,IEL)
      XV(II)=DCORVG(1,IV)
      YV(II)=DCORVG(2,IV)
      
      XM=XM+XV(II)
      YM=YM+YV(II)
 10   CONTINUE
C
      XM=XM/4D0
      YM=YM/4D0
      XM=2D0**(ILEV)*XM
      YM=2D0**(ILEV)*YM
C     Check, welche  Orientierung die Tangenten haben sollen
C     Man geht hier von einem regelmaessig verfeinerten EQ aus
      IF (MOD(INT(XM+YM),4).eq.0)THEN
         SIGNUM=-1d0
      ELSE
         SIGNUM=1d0
      ENDIF

      DO 20 II=1,4
         II1=II
         II2=II+1
         IF (II2.GT.4) II2=1
         TX(II1) =XV(II1)-XV(II2)
         TY(II1) =YV(II1)-YV(II2)
         TN = SQRT(Tx(ii1)**2+Ty(ii1)**2)
         TX(ii1)=signum*Tx(ii1)/TN
         TY(II1)=signum*TY(II1)/TN
         XE(II1) =0.5d0*(XV(II1)+XV(II2))
         YE(II1) =0.5d0*(YV(II1)+YV(II2))
20      CONTINUE
c
        DO 30 J=1,4
           S1U=S1U+ UU(j)*Tx(j)+ UV(j)*ty(j)
           S1D=S1D+ DU(j)*Tx(j)+ DV(j)*ty(j)
           S2=S2+tx(j)*tx(j)+ty(j)*ty(j)
 30     CONTINUE
c

 1    CONTINUE
c
c-----------------------------------------------------------------------
c     So, jetzt habe ich Alpha. Nun muß ich leider nochmal über alle
c     Zellen gehen :((
c-----------------------------------------------------------------------
c
        ALPHA=-S1u/S2!*0.50d0!  (Haelfte der optimalen korrektur!)
c        ALPHA=-S1d/S2*0.50d0!  (Haelfte der optimalen korrektur!)
c-----------------------------------------------------------------------
        lul2=0
        lum=0
        lalfa=0
        lquot=0
        CALL ZNEW (nel,1,lul2,'lul2  ')
        CALL ZNEW (nel,1,lum,'lum   ')
        CALL ZNEW (nel,1,lquot,'lalfa ')
        CALL ZNEW (nel,1,lalfa,'lalfa ')

        umax=0d0
        ummax=0d0
        qumax=0d0
        qumin=10d10
        umin=1d11
        ummin=1d11
        almaX=0d0
c        write (*,*) 'alphaSalfa', ilev,alpha
        Do 2 IEL =1,NEL
c
c
      XM=0d0
      YM=0d0   !Ergibt Elementmittelpunkt
      UM=0d0
      DM=0d0
      Ul2=0d0
      Dl2=0d0
      S1u=0d0
      S1d=0d0
      S2=0d0
c
      DO 2210  II=1,4
      IMID=KMID(II,IEL)
      I=IMID-NVT
      IU(II)=I
c
      UU(II)=U1(i)
      UV(II)=U2(i)
      DU(II)=D1(i)
      DV(II)=D2(i)
c
c     Berechne Mittelwert von u, sowie l_2 Wert
c

      UM=Um+UU(II)+UV(II)
      DM=Dm+DU(II)+DV(II)
c      IF (IEL.eq.14) write (*,*) ii,um, uu(ii),uv(ii)
      Ul2=Ul2+UU(ii)**2+UV(II)**2
      Dl2=Dl2+DU(ii)**2+DV(II)**2


      IF (KNPR(IMID).EQ.0)  THEN
         BDBC(II)=.FALSE.
      ELSE
         BDBC(II)=.TRUE.
      ENDIF
cc

      IV=KVERT(II,IEL)
      IV=KVERT(II,IEL)
      XV(II)=DCORVG(1,IV)
      YV(II)=DCORVG(2,IV)
      XM=XM+XV(II)
      YM=YM+YV(II)
 2210 continue
      do 2216 II=1,4
      II2=II+1
      IF (ii2.gt.4) II2=1
      XE(II) =0.5d0*(XV(II)+XV(II2))
      YE(II) =0.5d0*(YV(II)+YV(II2))
      IF ((ABS(MIN(XE(II),YE(II))).LE.1D-9).OR.
     *      ((ABS(1D0-MAX(YE(II),XE(II))).LE.1D-9))) THEN
           BRAN(ii)=.TRUE.
           ELSE
           BRAN(ii)=.false.
           ENDIF
 2216    CONTINUE
C
 
        Ul2=SQRT(Ul2)
        IF (ABS(Ul2).ge.UMAX) UMAX=abs(Ul2)
        IF (ABS(Um).ge.UmMAX) UMMAX=ABS(UM)
        IF (ABS(Ul2).le.UMin) UMin=abs(Ul2)
        IF (ABS(Um).le.UmMin) UMMin=ABS(UM)
        IF (ABS(ul2/um).ge.QUMax) QUmax=ABS(ul2/um) 
        IF (ABS(ul2/um).le.QUMin) QUmin=ABS(ul2/um) 
        
        DWORK(l(lul2)-1+iel)=abs(ul2)
cc        DWORK(l(lum)-1+iel)=(um)
c        DWORK(l(lquot)-1+iel)=ABS(dl2/dm)
        DWORK(l(lquot)-1+iel)=ABS(ul2/um)

c
      XM=XM/4D0
      YM=YM/4D0
      XM=2D0**(ILEV)*XM
      YM=2D0**(ILEV)*YM
C     Check, welche  Orientierung die Tangenten haben sollen
C     Man geht hier von einem regelmaessig verfeinerten EQ aus

      IF (MOD(INT(XM+YM),4).eq.0)THEN
         SIGNUM=-1d0
      ELSE
         SIGNUM=1d0
      ENDIF

      DO 2220 II=1,4
         II1=II
         II2=II+1
         IF (II2.GT.4) II2=1
         TX(II1) =XV(II1)-XV(II2)
         TY(II1) =YV(II1)-YV(II2)
         TN = SQRT(Tx(ii1)**2+Ty(ii1)**2)
         TX(ii1)=signum*Tx(ii1)/TN
         TY(II1)=signum*TY(II1)/TN
         XE(II1) =0.5d0*(XV(II1)+XV(II2))
         YE(II1) =0.5d0*(YV(II1)+YV(II2))
2220      CONTINUE
        DO 2030 J=1,4
           S1u=S1u+ UU(j)*Tx(j)+ UV(j)*ty(j)
           S1d=S1d+ DU(j)*Tx(j)+ DV(j)*ty(j)
           S2=S2+tx(j)*tx(j)+ty(j)*ty(j)
 2030     CONTINUE

cc        ALPHA=-S1u/S2!*0.50d0!  (Haelfte der optimalen korrektur!)
cc        ALPHAd=-S1d/S2!*0.50d0!  (Haelfte der optimalen korrektur!)
c
      UNEUl2=0d0
      UMNEU=0d0
c
      DO 2400  II=1,4
         xhelp=0.5d0
      I=IU(II)
c      IF (BDBC(II)) xhelp=1d0
      IF (BRAN(II)) xhelp=1d0
c           if (iel.eq.8) THEN
c              print *,ii ,U1(i),U2(i),'|',TX(ii),Ty(ii), ALPHA
c              ENDIF
c      UNEUL2=UNEUL2+(u1(i)+ALPHA*tx(ii))**2+(U2(i)+ALPHA*Ty(ii))**2
      UMNEU=UMNEU+(u1(i)+ALPHA*tx(ii))+(U2(i)+ALPHA*Ty(ii))
      U1(I)=u1(i)+ALPHA*tx(ii)*xhelp
      U2(I)=U2(i)+ALPHA*Ty(ii)*xhelp
      uneul2=uneul2+u1(i)**2+u2(i)**2
2400  CONTINUE

 5000 CONTINUE
      UNEuL2=SQRT(UNEUL2)
        DWORK(l(lalfa)-1+iel)=ABS(ALPHA/UL2)!alpha!UNEUL2!
        IF (ABS(alpha/UL2).ge.ALMAX) ALMAX=ABS(ALPHA/UL2)
        DWORK(l(lum)-1+iel)=uneul2!ABS(ALPHA/UL2)

 2    CONTINUE
      IF (ILEV.EQ.NLMAX) THEN
c         IF (ite.lt.50)
c     * print *, '111 Max-Qut ', QUMAX, ' Min-Quot: ', QUMIN , ite
         IF (ite.eq.50)
     * print *, '222 Max-Qut ', QUMAX, ' Min-Quot: ', QUMIN
         ENDIF
c

             CFILNM='#mat/alfa.          '
             WRITE(CFILNM(11:11),'(I1.1)') ILEV
             WRITE(CFILNM(12:12),'(A1)') '.'
             IF (ite.le.9) THEN
             WRITE(CFILNM(13:13),'(I1.1)') ite
             ELSEIF (ite.le.99) THEN 
             WRITE(CFILNM(13:14),'(I2.2)') ite
             ENDIF
             CALL XGMVXX (56,CFILNM,KVERT,DCORVG,DWORK(l(lul2))
     *       ,DWORK(l(lum)),DWORK(l(lalfa)),DWORK(l(lquot)))

        CALL ZDISP (0,lul2,'lul2  ')
        CALL ZDISP (0,lum,'lum   ')
        CALL ZDISP (0,lalfa,'lalfa ')
        CALL ZDISP (0,lquot,'lquot ')

        END


c========================================================================
C      Filterroutine, die tangentiale Kreisgeschwindigkeiten rausfiltert
c========================================================================
c
      SUBROUTINE GARLI3 (ILEV,NEL,NVT,NMT,DCORVG,KVERT,KMID,
     *                   U1,U2,ite)
      PARAMETER (NNVE=4,NNLEV=9)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299,NNWORK=1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
c      
      DIMENSION KVERT(NNVE,*), DCORVG(2,*),KMID(NNVE,*)
      DIMENSION U1(*),U2(*)
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
c      EXTERNAL f_u,f_v
      EXTERNAL PARX,PARY,TMAX,ue
      INCLUDE 'block.inc'
      SAVE
C-----------------------------------------------------------------------
C *** 

c      return
      lul2=0
      lum=0
      lalfa=0
      lquot=0
      lfil1=0
      lfil2=0
      CALL ZNEW (nel,1,lul2,'lul2  ')
      CALL ZNEW (nel,1,lum,'lum   ')
      CALL ZNEW (nel,1,lquot,'lalfa ')
      CALL ZNEW (nel,1,lalfa,'lalfa ')
      call ZNEW (nmt,1,lfil1,'lfil1 ')
      call ZNEW (nmt,1,lfil2,'lfil2 ')

      umax=0d0
      ummax=0d0
      qumax=0d0
      qumin=10d10
      umin=1d11
      ummin=1d11
      almaX=0d0

      DO  1 IEL=1,NEL

         CALL GARLIC (IEL,ILEV,NVT,DCORVG, KVERT,
     *               KMID,U1,U2, Dul2,dum,dalfa,dZBV, 
     *               dwork(l(lfil1)),dwork(l(lfil2)),1)

         DWORK(l(lul2)-1+iel) =dul2
         DWORK(l(lum)-1+iel)  =dZBV!im MomentL Ul2neu
         DWORK(l(lalfa)-1+iel)=ABS(dalfa/dul2)
         DWORK(l(lquot)-1+iel)=ABS(DUL2/dum)
 1    CONTINUE
c
      do 2 I=1,nmt
         U1(i)=u1(i)+DWORK(l(lfil1)-1+i)!
         U2(i)=u2(i)+DWORK(l(lfil2)-1+i)!
2     continue
c


c$$$      do 10 IMBD=1,NMBD
c$$$         IMID=KWORK(L(KLMBD(ILEV))-1+imbd)
c$$$         IF (IMID.LT.0) GOTO 10
c$$$         INPR=KWORK(L(KLNPR(iLEV))-1+KWORK(L(KLNPR(iLEV))-1+(NVT+IMID)))
c$$$         DPAR=DWORK(L(KLDBD(ILEV))-1+imbd)
c$$$         u1(imid)=0d0
c$$$         u2(imid)=0d0
c$$$ 10      continue

      CFILNM='#mat/alfa.          '
      WRITE(CFILNM(11:11),'(I1.1)') ILEV
      WRITE(CFILNM(12:12),'(A1)') '.'
      IF (ite.le.9) THEN
      WRITE(CFILNM(13:13),'(I1.1)') ite
      ELSEIF (ite.le.99) THEN 
      WRITE(CFILNM(13:14),'(I2.2)') ite
      ENDIF
      CALL XGMVXX (56,CFILNM,KVERT,DCORVG,DWORK(l(lul2)),
     *       DWORK(l(lum)),DWORK(l(lalfa)),DWORK(l(lquot)))

      CALL ZDISP (0,lul2,'lul2  ')
      CALL ZDISP (0,lum,'lum   ')
      CALL ZDISP (0,lalfa,'lalfa ')
      CALL ZDISP (0,lquot,'lquot ')
      CALL ZDISP (0,lfil1,'lfilt ')   
      CALL ZDISP (0,lfil2,'lfilt ')   

        END
