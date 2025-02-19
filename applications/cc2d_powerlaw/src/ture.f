From rainer@hermes.iwr.uni-heidelberg.de Thu Feb  1 16:40:56 2001
Date: Thu, 1 Feb 2001 16:35:45 +0100 (MET)
From: Rainer Schmachtel <rainer@hermes.iwr.uni-heidelberg.de>
To: Rainer Schmachtel <rainer@Math.Uni-Dortmund.DE>
Subject: /home/user/featflow/data/cc2d_powerlaw (fwd)



Father Antonio
Royal Spanish Administrator of Foreign Affairs


---------- Forwarded message ----------
Date: Thu, 1 Feb 2001 16:28:19 +0100 (MET)
From: Stefan Turek <ture@Math.Uni-Dortmund.DE>
To: rainer@hermes.iwr.uni-heidelberg.de
Subject: /home/user/featflow/data/cc2d_powerlaw

----------
X-Sun-Data-Type: default
X-Sun-Data-Description: default
X-Sun-Data-Name: vanca.f
X-Sun-Charset: us-ascii
X-Sun-Content-Lines: 461

************************************************************************
      SUBROUTINE  VANCAS (U1,U2,P,U1OLD,U2OLD,POLD,F1,F2,FP,
     *                    A,KCOLA,KLDA,B1,B2,KCOLB,KLDB,NU,NP,
     *                    KMBD,KVERT,KMID,KNPR,NMBD)
************************************************************************
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DOUBLE PRECISION BB1,BB2
      REAL  A,B1,B2
C
      PARAMETER (NNVE=4,NNLEV=9)
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
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC,IGRAD
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
      NA=KLDA(NU+1)-1
C-----------------------------------------------------------------------
C *** Update  UOLD:=U
      CALL  LCP1 (U1,U1OLD,NU)
      CALL  LCP1 (U2,U2OLD,NU)
      CALL  LCP1 (P ,POLD ,NP)
C
      distmx=0d0
      distmi=1d99
      deacmx=0d0
      deacmi=1d99
      deanmx=0d0
      deanmi=1d99
c
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
ccc      write(6,*) '14',AOFF1,AOFF2
      FF1(II)=FF1(II)-AOFF1*U1(J)
      FF2(II)=FF2(II)-AOFF2*U2(J)
 110  CONTINUE
      DO 111  IA=IA1,IA2
      J=KCOLA(IA)
      AOFF1=DBLE(A(NA+IA))
      AOFF2=DBLE(A(2*NA+IA))
ccc      write(6,*) '23',AOFF1,AOFF2
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
      CALL  ELUPDT (U1,U2,P,IEL,IU,BDBC,AA1,AA2,BB1,BB2,FF1,FF2,FFP,
     *              udist,umeanc,umean)
C
      distmx=max(distmx,udist)
      distmi=min(distmi,udist)
      deacmx=max(deacmx,umeanc)
      deacmi=min(deacmi,umeanc)
      deanmx=max(deanmx,umean)
      deanmi=min(deanmi,umean)
c
1     CONTINUE
C
      write(6,*) 'ILEV-MAX ',ILEV,distmx,deacmx,deanmx
      write(6,*) 'ILEV-MIN ',ILEV,distmi,deacmi,deanmi
      write(6,*) 
c
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
C
************************************************************************
      SUBROUTINE  ELUPDT  (U1,U2,P,IEL,IU,BDBC,AA1,AA2,BB1,BB2,FF1,FF2,
     *                     FFP,udist,umeanc,umean)
************************************************************************
*    Purpose: - Executes block Gauss-Seidel update on (U1,U2,P)
*               for all unknowns of the element IEL
************************************************************************
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DOUBLE PRECISION  BB1,BB2
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
C
      DIMENSION U1(*),U2(*),P(*)
C
C *** Local arrays for informations about the element
      DIMENSION AA1(4),AA2(4),BB1(4),BB2(4),FF1(4),FF2(4),IU(4),BDBC(4)
      DIMENSION AI1(4),AI2(4),DD1(4),DD2(4),UU1(4),UU2(4)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
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
      umean1=0d0
      umean2=0d0
      DO 10 II=1,4
      umean1=umean1+0.25D0*UU1(II)
      umean2=umean2+0.25D0*UU2(II)
10    CONTINUE
      umeanc=sqrt(umean1**2d0+umean2**2d0)
C
      udist=0d0
      umean=0d0
      DO 20 II=1,4
      udist=udist+0.25D0*(UU1(II)-umean1)**2+(UU2(II)-umean2)**2
      umean=umean+0.25D0*(UU1(II))**2+(UU2(II))**2
20    CONTINUE
      udist=sqrt(udist)
      umean=sqrt(umean)
      if (umeanc.ne.0d0) udist=udist/umeanc
      if (umeanc.ne.0d0) then
       udist=umean/umeanc
      else
       udist=1d0
      endif
c
c      if (ilev.eq.nlev) write(6,*) 'IEL ', iel,udist,umeanc,umean
c
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
     *               NSM,NSL,NSMFAC,IGRAD
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
      CALL  ELUPDT (U1,U2,P,IEL,IU,BDBC,AA1,AA2,BB1,BB2,FF1,FF2,FFP,
     *              udist,umeanc,umean)
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

----------
X-Sun-Data-Type: default
X-Sun-Data-Description: default
X-Sun-Data-Name: util.f
X-Sun-Charset: us-ascii
X-Sun-Content-Lines: 1064

************************************************************************
      SUBROUTINE FILT1C(DU1,DU2,DUN1,DUN2,DUC1,DUC2,
     *                  KVERT,KMID,KADJ,KMBD,NMBD,NMT,NVT,NEL)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KVERT(4,*),KMID(4,*),KADJ(4,*),KMBD(*)
      DIMENSION DU1(*),DU2(*),DUN1(*),DUN2(*),DUC1(*),DUC2(*)
C
C-----------------------------------------------------------------------
C
c      write(6,*) 'DUN'
      do 10 imt=1,nmt
c      write(6,*) imt,dun1(imt),dun2(imt)
10    continue
c
C *** nonconforming->constant:1
      CALL C2N2DM(DUC1,DUN1,KMID,KADJ,NEL,NMT,NVT,1)
      CALL C2N2DM(DUC2,DUN2,KMID,KADJ,NEL,NMT,NVT,1)
C
c      write(6,*) 'DUC'
      do 20 iel=1,nel
c      write(6,*) iel,duc1(iel),duc2(iel)
20    continue
c
C *** constant->nonconforming:0
      CALL C2N2DM(DUC1,DUN1,KMID,KADJ,NEL,NMT,NVT,0)
      CALL C2N2DM(DUC2,DUN2,KMID,KADJ,NEL,NMT,NVT,0)
C
      CALL  BDRY0 (DUN1,DUN2,KMBD,NMBD)
C
      dmaxu1=0d0
      dmaxu2=0d0
      dmaxn1=0d0
      dmaxn2=0d0
      dmax1 =0d0
      dmax2 =0d0
      dl1u  =0d0
      dl2u  =0d0
      dl1n  =0d0
      dl2n  =0d0
      dl1   =0d0
      dl2   =0d0
c      write(6,*) 'DU-DUN'
      do 30 imt=1,nmt
      dmaxu1=max(dmaxu1,abs(du1(imt)))
      dmaxu2=max(dmaxu2,abs(du2(imt)))
      dmaxn1=max(dmaxn1,abs(dun1(imt)))
      dmaxn2=max(dmaxn2,abs(dun2(imt)))
      dmax1 =max(dmax1 ,abs(dun1(imt)-du1(imt)))
      dmax2 =max(dmax2 ,abs(dun2(imt)-du2(imt)))
      dl1u  =dl1u+du1(imt)**2
      dl2u  =dl2u+du2(imt)**2
      dl1n  =dl1n+dun1(imt)**2
      dl2n  =dl2n+dun2(imt)**2
      dl1   =dl1 +(dun1(imt)-du1(imt))**2
      dl2   =dl2 +(dun2(imt)-du2(imt))**2
c      write(6,*) imt,dun1(imt),du1(imt),dun2(imt),du2(imt)
30    continue
ccc      write(6,*) 'MAX',max(dmax1/dmaxn1,dmax2/dmaxn2),
ccc     *                 max(dmaxu1,dmaxu2),max(dmaxn1,dmaxn2)
c      write(6,*) 'MAX1',dmax1,dmax1/dmaxn1,dmaxu1,dmaxn1
c      write(6,*) 'MAX2',dmax2,dmax2/dmaxn2,dmaxu2,dmaxn2
ccc      write(6,*) 'L2U',max(sqrt(dl1)/sqrt(dl1n),sqrt(dl2)/sqrt(dl2n)),
ccc     *                 max(sqrt(dl1u),sqrt(dl2u))/sqrt(dble(nmt)),
ccc     *                 max(sqrt(dl1n),sqrt(dl2n))/sqrt(dble(nmt))
c      write(6,*) 'L2U1',sqrt(dl1)/sqrt(dble(nmt)),sqrt(dl1)/sqrt(dl1n),
c     *            sqrt(dl1u)/sqrt(dble(nmt)),sqrt(dl1n)/sqrt(dble(nmt))
c      write(6,*) 'L2U2',sqrt(dl2)/sqrt(dble(nmt)),sqrt(dl2)/sqrt(dl2n),
c     *            sqrt(dl2u)/sqrt(dble(nmt)),sqrt(dl2n)/sqrt(dble(nmt))
c
      return
c
      omega=max(sqrt(dl1)/sqrt(dl1n),sqrt(dl2)/sqrt(dl2n))
      fomega=0.01d0*omega/(1d0+0.01d0*omega)
      write(6,*) fomega,omega,max(dmaxu1,dmaxu2),max(dmaxn1,dmaxn2)
c
ccc      return
      do 100 imt=1,nmt
      du1(imt)=du1(imt)+fomega*(dun1(imt)-du1(imt))
      du2(imt)=du2(imt)+fomega*(dun2(imt)-du2(imt))
100   continue
C
C
      END
C
C
C
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
      PINT=0.D0
      DO 2  IEL=1,NEL
  2   PINT=PINT+P(IEL)*DBLE(AREA(IEL))
C
      C=PINT/DBLE(AREA(NEL+1))
C
      DO 3  IEL=1,NEL
  3   P(IEL)=P(IEL)-C
C
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
C
C
----------
X-Sun-Data-Type: default
X-Sun-Data-Description: default
X-Sun-Data-Name: mgrout.f
X-Sun-Charset: us-ascii
X-Sun-Content-Lines: 2800

************************************************************************
*    M011   (edited from M010)                                         *
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
* M011                                                                 *
*                                                                      *
* Purpose  Solution of a linear system  A*X = B  using                 *
*          multigrid iteration                                         *
*          Double precision version                                    *
*                                                                      *
* Subroutines/functions called  LL21 , LLC1                            *
*                                                                      *
* Version from  08/25/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX       R*8    Starting address of vectors containing the           *
* DB       R*8    solution and the right hand side, DD is used as      *
* DD       R*8    auxiliary vector only                                *
* KOFFX           The actual starting address of DX on level ILEV      *
* KOFFB           is DX(1+KOFFX(ILEV)) (analogously for DB and DD)     *
* KOFFD           Total space required for all vectors is              *
*                 KNEQ(NLMIN)+...+KNEQ(NLMAX)                          *
*                 DX(1+KOFFX(NLMAX)) contains initial solution         *
*                 DB(1+KOFFB(NLMAX)) contains right hand side          *
* KNEQ     I*4    Number of equations for all levels                   *
* NLMAX    I*4    Iteration uses levels NLMIN to NLMAX,                *
* NLMIN    I*4    NLMAX  is the finest level                           *
* NIT      I*4    Maximum number of iterations                         *
*                 Iteration completed after reaching the finest level  *
* EPS1      R*8    Desired precision                                   *
*                 Stop if !!DEF!! < EPS1 !!DEFOLD!!                    *
* EPS2      R*8    Desired precision                                   *
*                 Stop if !!DEF!! < EPS2                               *
* KPRSM    I*4    Number of pre -smoothing steps for all levels        *
* KPOSM    I*4    Number of post-smoothing steps for all levels        *
* ICYCLE   I*4    <0: special cycle types (not yet implemented)        *
*                 =0: F-Cycle                                          *
*                 =1: V-Cycle                                          *
*                 =2: W-Cycle                                          *
*                 >2: Cycle of higher order                            *
* DAX      SUBR   CALL DAX(DX,DAX,NEQ,A1,A2)                           *
*                 Returns  DAX := A1*A*DX+A2*DAX                       *
* DPROL    SUBR   DPROL(DX1,DX2)                                       *
*                 Returns  DX2 := Prolongation(DX1)  to higher level   *
* DSTEP    SUBR   DPROL(DX,DD,DB,DSTEPP)                               *
*                 Returns  DSTEPP := optimal step size for correction  *
* DREST    SUBR   DREST(DD1,DD2)                                       *
*                 Returns  DD2 := Restriction(DD1)  to lower level     *
* DPRSM    SUBR   DPRSM(DX,DB,DD,NEQ,NPRSM)                            *
*                 Returns  DX  after  NPRSM:=KPRSM(ILEV)               *
*                 pre-smoothing steps, DD  is used as auxiliary vector *
* DPOSM    SUBR   Same as above, used for post-smoothing               *
* DEX      SUBR   DEX(DX,DB,DD,NEQ)                                    *
*                 Returns exact solution                               *
* DBC      SUBR   DBC(DX,NEQ)                                          *
*                 Copies boundary data onto components of  DX          *
* KIT0     I*4    auxiliary vectors of length  NLMAX                   *
* KIT      I*4                                                         *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Solution vector on DX(KOFFX(NLMAX))                  *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
* RHOLMG   R*8    Multigrid convergence rate                           *
*                                                                      *
************************************************************************
C
      SUBROUTINE  M011 (DX,DB,DD,KOFFX,KOFFB,KOFFD,KNEQ,NIT,ITE,
     *                  EPS1,EPS2,DAX,DPROL,DREST,DPRSM,DPOSM,DEX,DEXA,
     *                  DBC,DSTEP,KIT0,KIT,IREL,IDEFMG,RHOLMG,BMGEND)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNARR=299,NNLEV=9)
      DIMENSION DX(*),DB(*),DD(*),KOFFX(*),KOFFB(*),KOFFD(*)
      DIMENSION KNEQ(*),KIT0(*),KIT(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTIME/ TTMG,TTS,TTE,TTD,TTP,TTR,IMTIME
      SAVE 
C
      SUB='M011  '
      IF (ICHECK.GE.997) CALL OTRC('M011  ','05/13/91')
      IER=0
C
C
      BMGEND=.FALSE.
      BREL=IREL.EQ.1
      BMSG2=M.GE.2.OR.MT.GE.2
      MT0=MT
      MT=0
C
      BTIME=IMTIME.GT.0
      IF (BTIME) THEN
       IF (IMTIME.EQ.1) THEN
        TTMG=0D0
        TTS=0D0
        TTE=0D0
        TTD=0D0
        TTP=0D0
        TTR=0D0
       ENDIF
       CALL ZTIME(TTMG0)
      ENDIF
C
C
      NIT0=MAX(ITE,0)
      IF (ICYCLE.LT.0) THEN
       CALL WERR(-181,'M011  ')
       GOTO 99999
      ENDIF
      IF (NLMAX.LT.NLMIN.OR.NLMIN.LE.0.OR.NLMAX.GT.NNLEV) THEN
       CALL WERR(-182,'M011  ')
       GOTO 99999
      ENDIF
      ILEV=NLMAX
C
C *** special case - zero solution
      IF (BTIME) CALL ZTIME(TTD0)
      CALL LLI1(DB(1+KOFFB(NLMAX)),KNEQ(NLMAX),DMAX,INDMAX)
      IF (DMAX.LE.1D-12) THEN
       CALL LCP1(DB(1+KOFFB(NLMAX)),DX(1+KOFFX(NLMAX)),KNEQ(NLMAX))
       RHOLMG=0D0
       GOTO 99999
      ENDIF
C
      IF (BTIME) THEN
       CALL ZTIME(TTD1)
       TTD=TTD+TTD1-TTD0
      ENDIF
C
C *** special case - only one level
      IF (NLMIN.EQ.NLMAX) THEN
       IF (BTIME) CALL ZTIME(TTE0)
       CALL DEXA(DX(1+KOFFX(NLMAX)),DB(1+KOFFB(NLMAX)),
     *           DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),RHOLMG,EPS2,ITECG)
       IF (BTIME) THEN
        CALL ZTIME(TTE1)
        TTE=TTE+TTE1-TTE0
       ENDIF
       ITE=ITECG
       GOTO 99999
      ENDIF
C
C *** level counts
      KIT0(NLMAX)=1
      DO 2 ILEV=NLMIN+1,NLMAX-1
      IF (ICYCLE.EQ.0) THEN
       KIT0(ILEV)=2
      ELSE
       KIT0(ILEV)=ICYCLE
      ENDIF
2     CONTINUE
C
CCC      IF (BTIME) CALL ZTIME(TTS0)
CCC      IF (KPRSM(NLMAX).GT.0) 
CCC     * CALL DPRSM(DX(1+KOFFX(NLMAX)),DB(1+KOFFB(NLMAX)),
CCC     *            DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),KPRSM(NLMAX))
CCC
CCC      IF (BTIME) THEN
CCC       CALL ZTIME(TTS1)
CCC       TTS=TTS+TTS1-TTS0
CCC      ENDIF
C
      IF (IDEFMG.EQ.1) THEN
       IF (BTIME) CALL ZTIME(TTD0)
       CALL LCP1(DB(1+KOFFB(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX))
       CALL DAX(DX(1+KOFFX(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),
     *          -1D0,1D0)
C
       CALL LL21(DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),DEF)
       DEFOLD=DEF
C *** FD  is considered as initial defect
       FD=DEF
C
       IF (BTIME) THEN
        CALL ZTIME(TTD1)
        TTD=TTD+TTD1-TTD0
       ENDIF
C
       IF (DEF.LE.EPS2.AND..NOT.BREL) THEN
        ITE=0
        GOTO 1000
       ENDIF
      ENDIF
C
C
C *** Start multigrid iteration
      DO 100 ITE=1,NIT
C
C *** initialize level counts for all levels
      DO 101 ILEV=NLMIN,NLMAX
101   KIT(ILEV)=KIT0(ILEV)
C
      ILEV=NLMAX
110   IF (ILEV.NE.NLMIN) THEN
C
CCC      IF (ILEV.NE.NLMAX) THEN
CCC ***  defect on finest level already available
C
        IF (BTIME) CALL ZTIME(TTS0)
        IF (KPRSM(ILEV).GT.0) 
     *   CALL DPRSM(DX(1+KOFFX(ILEV)),DB(1+KOFFB(ILEV)),
     *              DD(1+KOFFD(ILEV)),KNEQ(ILEV),KPRSM(ILEV))
C
        IF (BTIME) THEN
         CALL ZTIME(TTS1)
         TTS=TTS+TTS1-TTS0
        ENDIF
C
        IF (BTIME) CALL ZTIME(TTD0)
        CALL LCP1(DB(1+KOFFB(ILEV)),DD(1+KOFFD(ILEV)),KNEQ(ILEV))
        CALL DAX(DX(1+KOFFX(ILEV)),DD(1+KOFFD(ILEV)),KNEQ(ILEV),
     *           -1D0,1D0)
C
        IF (BTIME) THEN
         CALL ZTIME(TTD1)
         TTD=TTD+TTD1-TTD0
        ENDIF
C
CCC       ENDIF
C
       ILEV=ILEV-1
C
C *** restriction of defect
       IF (BTIME) CALL ZTIME(TTR0)
       CALL DREST(DD(1+KOFFD(ILEV+1)),DB(1+KOFFB(ILEV)))
C
       IF (BTIME) THEN
        CALL ZTIME(TTR1)
        TTR=TTR+TTR1-TTR0
       ENDIF
C
C ***  choose zero as initial vector on lower level
       IF (BTIME) CALL ZTIME(TTD0)
       CALL LCL1(DX(1+KOFFX(ILEV)),KNEQ(ILEV))
       CALL DBC(DB(1+KOFFB(ILEV)),KNEQ(ILEV))
       IF (BTIME) THEN
        CALL ZTIME(TTD1)
        TTD=TTD+TTD1-TTD0
       ENDIF
C
       GOTO 110
C
      ENDIF
C
C *** exact solution on lowest level
      IF (BTIME) CALL ZTIME(TTE0)
      CALL DEX(DX(1+KOFFX(NLMIN)),DB(1+KOFFB(NLMIN)),DD(1+KOFFD(NLMIN)),
     *         KNEQ(NLMIN),RHOLMG)
C
      IF (BTIME) THEN
       CALL ZTIME(TTE1)
       TTE=TTE+TTE1-TTE0
      ENDIF
C
C
130   IF (ILEV.NE.NLMAX) THEN
       ILEV=ILEV+1
C
C *** DPROL  returns  DD:=PROL(DX)
      IF (BTIME) CALL ZTIME(TTP0)
      CALL DPROL(DX(1+KOFFX(ILEV-1)),DD(1+KOFFD(ILEV)))
      CALL DBC(DD(1+KOFFD(ILEV)),KNEQ(ILEV))
C
      CALL DSTEP(DX(1+KOFFX(ILEV)),DD(1+KOFFD(ILEV)),DB(1+KOFFB(ILEV)),
     *           DSTEPP)
      CALL LLC1(DD(1+KOFFD(ILEV)),DX(1+KOFFX(ILEV)),KNEQ(ILEV),DSTEPP,
     *          1D0)
C
      IF (BTIME) THEN
       CALL ZTIME(TTP1)
       TTP=TTP+TTP1-TTP0
      ENDIF
C
C
C *** Post-smoothing
      IF (BTIME) CALL ZTIME(TTS0)
      IF (KPOSM(ILEV).GT.0) 
     * CALL DPOSM(DX(1+KOFFX(ILEV)),DB(1+KOFFB(ILEV)),DD(1+KOFFD(ILEV)),
     *            KNEQ(ILEV),KPOSM(ILEV))
C
      IF (BTIME) THEN
       CALL ZTIME(TTS1)
       TTS=TTS+TTS1-TTS0
      ENDIF
C
       KIT(ILEV)=KIT(ILEV)-1
       IF (KIT(ILEV).EQ.0) THEN
        IF (ICYCLE.EQ.0) THEN
         KIT(ILEV)=1
        ELSE
         KIT(ILEV)=KIT0(ILEV)
        ENDIF
        GOTO 130
       ELSE
        GOTO 110
       ENDIF
C
      ENDIF
C
C
CCC      IF (BTIME) CALL ZTIME(TTS0)
CCC      IF (KPRSM(NLMAX).GT.0)  
CCC     * CALL DPRSM(DX(1+KOFFX(NLMAX)),DB(1+KOFFB(NLMAX)),
CCC     *           DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),KPRSM(NLMAX))
CCC
CCC      IF (BTIME) THEN
CCC       CALL ZTIME(TTS1)
CCC       TTS=TTS+TTS1-TTS0
CCC      ENDIF
C
C
      IF (IDEFMG.EQ.1) THEN
       IF (BTIME) CALL ZTIME(TTD0)
       CALL LCP1(DB(1+KOFFB(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX))
       CALL DAX(DX(1+KOFFX(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),
     *          -1D0,1D0)
C
      CALL LL21(DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),DEF)
       IF (MT.GE.9) WRITE (6,10001) ITE,DEF
       IF (BMSG2) THEN
        MT=MT0
        WRITE (CPARAM,'(I15,D25.16)') ITE,DEF
        CALL OMSG(73,'M011  ')
        WRITE (MPROT,10001) ITE,DEF
        WRITE (*,10001) ITE,DEF
        MT=0
       ENDIF
C
       IF (BTIME) THEN
        CALL ZTIME(TTD1)
        TTD=TTD+TTD1-TTD0
       ENDIF
C
C=======================================================================
C ***  Unexpected STOP !!!
C=======================================================================
       IF ((DEF/DEFOLD.GT.1D3).OR.(DEF/FD.GT.5D3)) THEN
        BMGEND=.TRUE.
        GOTO 1000
       ENDIF
C
      DEFOLD=DEF
       IF (BREL) THEN
        IF (DEF.LE.FD*EPS1.AND.DEF.LE.EPS2.AND.ITE.GE.NIT0) GOTO 1000
       ELSE
        IF (DEF.LE.EPS2) GOTO 1000
       ENDIF
      ENDIF
C
100   CONTINUE
C
      ITE=NIT
      IF (IDEFMG.EQ.1) THEN
       MT=MT0
       RHOLMG=DEF/FD
       RHOLMG=RHOLMG**(1D0/DBLE(NIT))
       IF ((RHOLMG).GE.1D0) BMGEND=.TRUE.
       IF (MT.GE.9) WRITE (6,10000) NIT,DEF,FD,DEF/FD,RHOLMG
       IF (BMSG2) THEN
        WRITE (CPARAM,'(I15,4D25.16)') NIT,DEF,FD,DEF/FD,RHOLMG
        CALL OMSG(72,'M011  ')
        WRITE (MPROT,10000) NIT,DEF,FD,DEF/FD,RHOLMG
        WRITE (MPROT,*) 'IER=1 in M011'
       ENDIF
      ELSE
       RHOLMG=DBLE(NIT)
      ENDIF
C
      GOTO 99999
C
1000  IER=0
      MT=MT0
      RHOLMG=0.D0
      IF (FD.GE.1D-70) RHOLMG=(DEF/FD)**(1D0/DBLE(ITE))
      IF (MT.GE.9) WRITE (6,'(I15,4D25.16)') 'ite',ITE,DEF,FD,DEF/FD,RHOLMG
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,4D25.16)') 'ite',ITE,DEF,FD,DEF/FD,RHOLMG
       CALL OMSG(72,'M011  ')
       WRITE (CPARAM,'(D25.16)')  RHOLMG
       CALL OMSG(76,'M011  ')
       WRITE (MPROT,10002) 'ite',ITE,DEF,FD,DEF/FD,RHOLMG
      ENDIF
C
99999 MT=MT0 
      IF (BTIME) THEN
       CALL ZTIME(TTMG1)
       TTMG=TTMG+TTMG1-TTMG0
      ENDIF
C
10000 FORMAT (//' CONVERGENCE OF MG-METHOD FAILED'/
     *        ' NORM OF DEFECT AFTER',I6,' ITERATIONS',D12.3/
     *        ' NORM OF INITIAL DEFECT  ',D12.3/
     *        ' NORM OF DEF/INIT. DEF   ',D12.3/
     *        ' KONVERGENZRATE          ',D12.3)
10001 FORMAT (' ITERATION',I6,5X,'  RESIDUUM =',D12.3)
10002 FORMAT (' NUMBER OF ITERATIONS    ',I10/
     *        ' NORM OF DEFECT          ',D12.3/
     *        ' NORM OF INITIAL DEFECT  ',D12.3/
     *        ' NORM OF DEF/INIT. DEF   ',D12.3/
     *        ' KONVERGENZRATE          ',D12.3)
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE   PROLU  (DU1,DV1,DP1,DU2,DV2,DP2,
     *                     KVERT1,KMID1,KADJ1,NEQ1,NEL1,NVT1,
     *                     KVERT2,KMID2,KADJ2,NEQ2,NEL2,NVT2)
************************************************************************
*    Purpose:    Interpolates the coarse grid vector (DU1,DV1,DP1) to
*                the fine grid vector (DU2,DV2,DP2)
*-----------------------------------------------------------------------
*    Input:
*      DU1,DV1,DP1           - coarse grid vector
*      KVERT1,KMID1,..,NVT1  - data of the coarse grid
*      KVERT2,KMID2,..,NVT2  - data of the fine grid
*
*    Output:
*      DU2,DV2,DP2           - interpolated fine grid vector
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      PARAMETER (A1=0.46875D0,A2=-0.09375D0,A3=-0.03125D0,
     *           A4=0.15625D0)
      PARAMETER (A5=0.5625D0,A6=0.1875D0,A7=0.0625D0,A8=0.1875D0)
C
      DIMENSION DU1(*),DV1(*),DP1(*),  DU2(*),DV2(*),DP2(*),
     *          KVERT1(NNVE,*),KMID1(NNVE,*),KADJ1(NNVE,*),
     *          KVERT2(NNVE,*),KMID2(NNVE,*),KADJ2(NNVE,*)
      SAVE
C-----------------------------------------------------------------------
C *** Zero initialization of (DU2,DV2)
      CALL  LCL1 (DU2,NEQ2)
      CALL  LCL1 (DV2,NEQ2)
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      DUH1=DU1(IM1)
      DUH2=DU1(IM2)
      DUH3=DU1(IM3)
      DUH4=DU1(IM4)
C
      DVH1=DV1(IM1)
      DVH2=DV1(IM2)
      DVH3=DV1(IM3)
      DVH4=DV1(IM4)
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C *** Prolongation of pressure
C
      DPH=DP1(IEL1)
      DP2(IELH1)=DPH
      DP2(IELH2)=DPH
      DP2(IELH3)=DPH
      DP2(IELH4)=DPH
C
C *** The edge IM1 and the corresponding fine inner node
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH1)-NVT2
       IB=KMID2(4,IELH2)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4
       DU2(IB)=DU2(IB)+   A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4
       DV2(IA)=DV2(IA)+   A1*DVH1+A2*DVH2+A3*DVH3+A4*DVH4
       DV2(IB)=DV2(IB)+   A1*DVH1+A4*DVH2+A3*DVH3+A2*DVH4
      ELSE
       IA=KMID2(1,IELH1)-NVT2
       IB=KMID2(4,IELH2)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4)
       DV2(IA)=DV2(IA)+2D0*(A1*DVH1+A2*DVH2+A3*DVH3+A4*DVH4)
       DV2(IB)=DV2(IB)+2D0*(A1*DVH1+A4*DVH2+A3*DVH3+A2*DVH4)
      ENDIF
      IC=KMID2(2,IELH1)-NVT2
      DU2(IC)=A5*DUH1+A6*(DUH2+DUH4)+A7*DUH3
      DV2(IC)=A5*DVH1+A6*(DVH2+DVH4)+A7*DVH3
C
C *** The edge IM2 and the corresponding fine inner node
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH2)-NVT2
       IB=KMID2(4,IELH3)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1
       DU2(IB)=DU2(IB)+   A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1
       DV2(IA)=DV2(IA)+   A1*DVH2+A2*DVH3+A3*DVH4+A4*DVH1
       DV2(IB)=DV2(IB)+   A1*DVH2+A4*DVH3+A3*DVH4+A2*DVH1
      ELSE
       IA=KMID2(1,IELH2)-NVT2
       IB=KMID2(4,IELH3)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1)
       DV2(IA)=DV2(IA)+2D0*(A1*DVH2+A2*DVH3+A3*DVH4+A4*DVH1)
       DV2(IB)=DV2(IB)+2D0*(A1*DVH2+A4*DVH3+A3*DVH4+A2*DVH1)
      ENDIF
      IC=KMID2(2,IELH2)-NVT2
      DU2(IC)=A5*DUH2+A6*(DUH3+DUH1)+A7*DUH4
      DV2(IC)=A5*DVH2+A6*(DVH3+DVH1)+A7*DVH4
C
C *** The edge IM3 and the corresponding fine inner node
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH3)-NVT2
       IB=KMID2(4,IELH4)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2
       DU2(IB)=DU2(IB)+   A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2
       DV2(IA)=DV2(IA)+   A1*DVH3+A2*DVH4+A3*DVH1+A4*DVH2
       DV2(IB)=DV2(IB)+   A1*DVH3+A4*DVH4+A3*DVH1+A2*DVH2
      ELSE
       IA=KMID2(1,IELH3)-NVT2
       IB=KMID2(4,IELH4)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2)
       DV2(IA)=DV2(IA)+2D0*(A1*DVH3+A2*DVH4+A3*DVH1+A4*DVH2)
       DV2(IB)=DV2(IB)+2D0*(A1*DVH3+A4*DVH4+A3*DVH1+A2*DVH2)
      ENDIF
      IC=KMID2(2,IELH3)-NVT2
      DU2(IC)=A5*DUH3+A6*(DUH4+DUH2)+A7*DUH1
      DV2(IC)=A5*DVH3+A6*(DVH4+DVH2)+A7*DVH1
C
C *** The edge IM4 and the corresponding fine inner node
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH4)-NVT2
       IB=KMID2(4,IELH1)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3
       DU2(IB)=DU2(IB)+   A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3
       DV2(IA)=DV2(IA)+   A1*DVH4+A2*DVH1+A3*DVH2+A4*DVH3
       DV2(IB)=DV2(IB)+   A1*DVH4+A4*DVH1+A3*DVH2+A2*DVH3
      ELSE
       IA=KMID2(1,IELH4)-NVT2
       IB=KMID2(4,IELH1)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3)
       DV2(IA)=DV2(IA)+2D0*(A1*DVH4+A2*DVH1+A3*DVH2+A4*DVH3)
       DV2(IB)=DV2(IB)+2D0*(A1*DVH4+A4*DVH1+A3*DVH2+A2*DVH3)
      ENDIF
      IC=KMID2(2,IELH4)-NVT2
      DU2(IC)=A5*DUH4+A6*(DUH1+DUH3)+A7*DUH2
      DV2(IC)=A5*DVH4+A6*(DVH1+DVH3)+A7*DVH2
C
C
10    CONTINUE
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE   PROLUM (DU1,DV1,DU2,DV2,
     *                     KVERT1,KMID1,KADJ1,NEQ1,NEL1,NVT1,
     *                     KVERT2,KMID2,KADJ2,NEQ2,NEL2,NVT2)
************************************************************************
*    Purpose:    Interpolates the coarse grid vector (DU1,DV1) to
*                the fine grid vector (DU2,DV2)
*-----------------------------------------------------------------------
*    Input:
*      DU1,DV1               - coarse grid vector
*      KVERT1,KMID1,..,NVT1  - data of the coarse grid
*      KVERT2,KMID2,..,NVT2  - data of the fine grid
*
*    Output:
*      DU2,DV2               - interpolated fine grid vector
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      PARAMETER (A1=0.46875D0,A2=-0.09375D0,A3=-0.03125D0,
     *           A4=0.15625D0)
      PARAMETER (A5=0.5625D0,A6=0.1875D0,A7=0.0625D0,A8=0.1875D0)
C
      DIMENSION DU1(*),DV1(*),DU2(*),DV2(*),
     *          KVERT1(NNVE,*),KMID1(NNVE,*),KADJ1(NNVE,*),
     *          KVERT2(NNVE,*),KMID2(NNVE,*),KADJ2(NNVE,*)
      SAVE
C-----------------------------------------------------------------------
C *** Zero initialization of (DU2,DV2)
      CALL  LCL1 (DU2,NEQ2)
      CALL  LCL1 (DV2,NEQ2)
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      DUH1=DU1(IM1)
      DUH2=DU1(IM2)
      DUH3=DU1(IM3)
      DUH4=DU1(IM4)
C
      DVH1=DV1(IM1)
      DVH2=DV1(IM2)
      DVH3=DV1(IM3)
      DVH4=DV1(IM4)
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C *** The edge IM1 and the corresponding fine inner node
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH1)-NVT2
       IB=KMID2(4,IELH2)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4
       DU2(IB)=DU2(IB)+   A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4
       DV2(IA)=DV2(IA)+   A1*DVH1+A2*DVH2+A3*DVH3+A4*DVH4
       DV2(IB)=DV2(IB)+   A1*DVH1+A4*DVH2+A3*DVH3+A2*DVH4
      ELSE
       IA=KMID2(1,IELH1)-NVT2
       IB=KMID2(4,IELH2)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4)
       DV2(IA)=DV2(IA)+2D0*(A1*DVH1+A2*DVH2+A3*DVH3+A4*DVH4)
       DV2(IB)=DV2(IB)+2D0*(A1*DVH1+A4*DVH2+A3*DVH3+A2*DVH4)
      ENDIF
      IC=KMID2(2,IELH1)-NVT2
      DU2(IC)=A5*DUH1+A6*(DUH2+DUH4)+A7*DUH3
      DV2(IC)=A5*DVH1+A6*(DVH2+DVH4)+A7*DVH3
C
C *** The edge IM2 and the corresponding fine inner node
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH2)-NVT2
       IB=KMID2(4,IELH3)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1
       DU2(IB)=DU2(IB)+   A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1
       DV2(IA)=DV2(IA)+   A1*DVH2+A2*DVH3+A3*DVH4+A4*DVH1
       DV2(IB)=DV2(IB)+   A1*DVH2+A4*DVH3+A3*DVH4+A2*DVH1
      ELSE
       IA=KMID2(1,IELH2)-NVT2
       IB=KMID2(4,IELH3)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1)
       DV2(IA)=DV2(IA)+2D0*(A1*DVH2+A2*DVH3+A3*DVH4+A4*DVH1)
       DV2(IB)=DV2(IB)+2D0*(A1*DVH2+A4*DVH3+A3*DVH4+A2*DVH1)
      ENDIF
      IC=KMID2(2,IELH2)-NVT2
      DU2(IC)=A5*DUH2+A6*(DUH3+DUH1)+A7*DUH4
      DV2(IC)=A5*DVH2+A6*(DVH3+DVH1)+A7*DVH4
C
C *** The edge IM3 and the corresponding fine inner node
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH3)-NVT2
       IB=KMID2(4,IELH4)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2
       DU2(IB)=DU2(IB)+   A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2
       DV2(IA)=DV2(IA)+   A1*DVH3+A2*DVH4+A3*DVH1+A4*DVH2
       DV2(IB)=DV2(IB)+   A1*DVH3+A4*DVH4+A3*DVH1+A2*DVH2
      ELSE
       IA=KMID2(1,IELH3)-NVT2
       IB=KMID2(4,IELH4)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2)
       DV2(IA)=DV2(IA)+2D0*(A1*DVH3+A2*DVH4+A3*DVH1+A4*DVH2)
       DV2(IB)=DV2(IB)+2D0*(A1*DVH3+A4*DVH4+A3*DVH1+A2*DVH2)
      ENDIF
      IC=KMID2(2,IELH3)-NVT2
      DU2(IC)=A5*DUH3+A6*(DUH4+DUH2)+A7*DUH1
      DV2(IC)=A5*DVH3+A6*(DVH4+DVH2)+A7*DVH1
C
C *** The edge IM4 and the corresponding fine inner node
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH4)-NVT2
       IB=KMID2(4,IELH1)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3
       DU2(IB)=DU2(IB)+   A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3
       DV2(IA)=DV2(IA)+   A1*DVH4+A2*DVH1+A3*DVH2+A4*DVH3
       DV2(IB)=DV2(IB)+   A1*DVH4+A4*DVH1+A3*DVH2+A2*DVH3
      ELSE
       IA=KMID2(1,IELH4)-NVT2
       IB=KMID2(4,IELH1)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3)
       DV2(IA)=DV2(IA)+2D0*(A1*DVH4+A2*DVH1+A3*DVH2+A4*DVH3)
       DV2(IB)=DV2(IB)+2D0*(A1*DVH4+A4*DVH1+A3*DVH2+A2*DVH3)
      ENDIF
      IC=KMID2(2,IELH4)-NVT2
      DU2(IC)=A5*DUH4+A6*(DUH1+DUH3)+A7*DUH2
      DV2(IC)=A5*DVH4+A6*(DVH1+DVH3)+A7*DVH2
C
C
10    CONTINUE
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE  RESTRD  (DU1,DV1,DP1,DU2,DV2,DP2,
     *                     KVERT1,KMID1,KADJ1,NEQ1,NEL1,NVT1,
     *                     KVERT2,KMID2,KADJ2,NEQ2,NEL2,NVT2)
************************************************************************
*    Purpose:    Restricts the  fine grid defect vector (DU2,DV2,DP2)
*                to the coarse grid defect vector (DU1,DV1,DP1)
*-----------------------------------------------------------------------
*    Input:
*      DU2,DV2,DP2           - fine grid defect vector
*      KVERT1,KMID1,..,NVT1  - data of the coarse grid
*      KVERT2,KMID2,..,NVT2  - data of the fine grid
*
*    Output:
*      DU1,DV1,DP1           - coarse grid defect vector
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      PARAMETER (A1=0.9375D0,A2=-0.09375D0,A3=-0.03125D0,
     *           A4=0.15625D0)
      PARAMETER (A5=0.5625D0,A6=0.1875D0,A7=0.0625D0,A8=0.1875D0)
C
      DIMENSION DU1(*),DV1(*),DP1(*),  DU2(*),DV2(*),DP2(*),
     *          KVERT1(NNVE,*),KMID1(NNVE,*),KADJ1(NNVE,*),
     *          KVERT2(NNVE,*),KMID2(NNVE,*),KADJ2(NNVE,*)
C
      SAVE
C-----------------------------------------------------------------------
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C *** Restriction of pressure
C
      DP1(IEL1)= DP2(IELH1)+DP2(IELH2)+DP2(IELH3)+DP2(IELH4)
C
C
      I1=KMID2(1,IELH1)-NVT2
      I2=KMID2(4,IELH2)-NVT2
      I3=KMID2(1,IELH2)-NVT2
      I4=KMID2(4,IELH3)-NVT2
      I5=KMID2(1,IELH3)-NVT2
      I6=KMID2(4,IELH4)-NVT2
      I7=KMID2(1,IELH4)-NVT2
      I8=KMID2(4,IELH1)-NVT2
      I9=KMID2(2,IELH1)-NVT2
      I10=KMID2(2,IELH2)-NVT2
      I11=KMID2(2,IELH3)-NVT2
      I12=KMID2(2,IELH4)-NVT2
C
      DUH1= DU2(I1)
      DUH2= DU2(I2)
      DUH3= DU2(I3)
      DUH4= DU2(I4)
      DUH5= DU2(I5)
      DUH6= DU2(I6)
      DUH7= DU2(I7)
      DUH8= DU2(I8)
      DUH9= DU2(I9)
      DUH10=DU2(I10)
      DUH11=DU2(I11)
      DUH12=DU2(I12)
C
      DVH1= DV2(I1)
      DVH2= DV2(I2)
      DVH3= DV2(I3)
      DVH4= DV2(I4)
      DVH5= DV2(I5)
      DVH6= DV2(I6)
      DVH7= DV2(I7)
      DVH8= DV2(I8)
      DVH9= DV2(I9)
      DVH10=DV2(I10)
      DVH11=DV2(I11)
      DVH12=DV2(I12)
C
C *** The edge IM1
C
      IF (KADJ1(1,IEL1).NE.0) THEN
C     case of an inner edge
        IF (KADJ1(1,IEL1).GT.IEL1) THEN
           DU1(IM1)= A1*(DUH1+DUH2)+A2*(DUH4+DUH7)
     *                 +A3*(DUH5+DUH6)+A4*(DUH3+DUH8)
     *                   +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
           DV1(IM1)= A1*(DVH1+DVH2)+A2*(DVH4+DVH7)
     *                 +A3*(DVH5+DVH6)+A4*(DVH3+DVH8)
     *                   +A5*DVH9+A6*(DVH10+DVH12)+A7*DVH11
        ELSE
           DU1(IM1)=DU1(IM1)+A2*(DUH4+DUH7)
     *                        +A3*(DUH5+DUH6)+A4*(DUH3+DUH8)
     *                          +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
           DV1(IM1)=DV1(IM1)+A2*(DVH4+DVH7)
     *                        +A3*(DVH5+DVH6)+A4*(DVH3+DVH8)
     *                          +A5*DVH9+A6*(DVH10+DVH12)+A7*DVH11
        ENDIF
      ELSE
       DU1(IM1)=     A1*(DUH1+DUH2)+2D0*A2*(DUH4+DUH7)
     *          +2D0*A3*(DUH5+DUH6)+2D0*A4*(DUH3+DUH8)
     *          +    A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
       DV1(IM1)=     A1*(DVH1+DVH2)+2D0*A2*(DVH4+DVH7)
     *          +2D0*A3*(DVH5+DVH6)+2D0*A4*(DVH3+DVH8)
     *          +    A5*DVH9+A6*(DVH10+DVH12)+A7*DVH11
      ENDIF
C
C *** The edge IM2
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(2,IEL1).GT.IEL1) THEN
           DU1(IM2)= A1*(DUH3+DUH4)+A2*(DUH6+DUH1)
     *                 +A3*(DUH7+DUH8)+A4*(DUH5+DUH2)
     *                   +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
           DV1(IM2)= A1*(DVH3+DVH4)+A2*(DVH6+DVH1)
     *                 +A3*(DVH7+DVH8)+A4*(DVH5+DVH2)
     *                   +A5*DVH10+A6*(DVH11+DVH9)+A7*DVH12
        ELSE
           DU1(IM2)=DU1(IM2)+A2*(DUH6+DUH1)
     *                        +A3*(DUH7+DUH8)+A4*(DUH5+DUH2)
     *                          +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
           DV1(IM2)=DV1(IM2)+A2*(DVH6+DVH1)
     *                        +A3*(DVH7+DVH8)+A4*(DVH5+DVH2)
     *                          +A5*DVH10+A6*(DVH11+DVH9)+A7*DVH12
        ENDIF
      ELSE
       DU1(IM2)=     A1*(DUH3+DUH4)+2D0*A2*(DUH6+DUH1)
     *          +2D0*A3*(DUH7+DUH8)+2D0*A4*(DUH5+DUH2)
     *          +    A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
       DV1(IM2)=     A1*(DVH3+DVH4)+2D0*A2*(DVH6+DVH1)
     *          +2D0*A3*(DVH7+DVH8)+2D0*A4*(DVH5+DVH2)
     *          +    A5*DVH10+A6*(DVH11+DVH9)+A7*DVH12
      ENDIF
C
C *** The edge IM3
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(3,IEL1).GT.IEL1) THEN
           DU1(IM3)= A1*(DUH5+DUH6)+A2*(DUH8+DUH3)
     *                 +A3*(DUH1+DUH2)+A4*(DUH7+DUH4)
     *                   +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
           DV1(IM3)= A1*(DVH5+DVH6)+A2*(DVH8+DVH3)
     *                 +A3*(DVH1+DVH2)+A4*(DVH7+DVH4)
     *                   +A5*DVH11+A6*(DVH12+DVH10)+A7*DVH9
        ELSE
           DU1(IM3)=DU1(IM3)+A2*(DUH8+DUH3)
     *                        +A3*(DUH1+DUH2)+A4*(DUH7+DUH4)
     *                          +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
           DV1(IM3)=DV1(IM3)+A2*(DVH8+DVH3)
     *                        +A3*(DVH1+DVH2)+A4*(DVH7+DVH4)
     *                          +A5*DVH11+A6*(DVH12+DVH10)+A7*DVH9
        ENDIF
      ELSE
       DU1(IM3)=     A1*(DUH5+DUH6)+2D0*A2*(DUH8+DUH3)
     *          +2D0*A3*(DUH1+DUH2)+2D0*A4*(DUH7+DUH4)
     *          +    A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
       DV1(IM3)=     A1*(DVH5+DVH6)+2D0*A2*(DVH8+DVH3)
     *          +2D0*A3*(DVH1+DVH2)+2D0*A4*(DVH7+DVH4)
     *          +    A5*DVH11+A6*(DVH12+DVH10)+A7*DVH9
      ENDIF
C
C *** The edge IM4
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(4,IEL1).GT.IEL1) THEN
           DU1(IM4)= A1*(DUH7+DUH8)+A2*(DUH2+DUH5)
     *                 +A3*(DUH3+DUH4)+A4*(DUH1+DUH6)
     *                   +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
           DV1(IM4)= A1*(DVH7+DVH8)+A2*(DVH2+DVH5)
     *                 +A3*(DVH3+DVH4)+A4*(DVH1+DVH6)
     *                   +A5*DVH12+A6*(DVH9+DVH11)+A7*DVH10
        ELSE
           DU1(IM4)=DU1(IM4)+A2*(DUH2+DUH5)
     *                        +A3*(DUH3+DUH4)+A4*(DUH1+DUH6)
     *                          +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
           DV1(IM4)=DV1(IM4)+A2*(DVH2+DVH5)
     *                        +A3*(DVH3+DVH4)+A4*(DVH1+DVH6)
     *                          +A5*DVH12+A6*(DVH9+DVH11)+A7*DVH10
        ENDIF
      ELSE
       DU1(IM4)=     A1*(DUH7+DUH8)+2D0*A2*(DUH2+DUH5)
     *          +2D0*A3*(DUH3+DUH4)+2D0*A4*(DUH1+DUH6)
     *          +    A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
       DV1(IM4)=     A1*(DVH7+DVH8)+2D0*A2*(DVH2+DVH5)
     *          +2D0*A3*(DVH3+DVH4)+2D0*A4*(DVH1+DVH6)
     *          +    A5*DVH12+A6*(DVH9+DVH11)+A7*DVH10
      ENDIF
C
C
10    CONTINUE
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE  RESTRM  (DU1,DV1,DU2,DV2,
     *                     KVERT1,KMID1,KADJ1,NEQ1,NEL1,NVT1,
     *                     KVERT2,KMID2,KADJ2,NEQ2,NEL2,NVT2)
************************************************************************
*    Purpose:    Restricts the  fine grid defect vector (DU2,DV2)
*                to the coarse grid defect vector (DU1,DV1)
*-----------------------------------------------------------------------
*    Input:
*      DU2,DV2               - fine grid defect vector
*      KVERT1,KMID1,..,NVT1  - data of the coarse grid
*      KVERT2,KMID2,..,NVT2  - data of the fine grid
*
*    Output:
*      DU1,DV1               - coarse grid defect vector
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      PARAMETER (A1=0.9375D0,A2=-0.09375D0,A3=-0.03125D0,
     *           A4=0.15625D0)
      PARAMETER (A5=0.5625D0,A6=0.1875D0,A7=0.0625D0,A8=0.1875D0)
C
      DIMENSION DU1(*),DV1(*),DU2(*),DV2(*),
     *          KVERT1(NNVE,*),KMID1(NNVE,*),KADJ1(NNVE,*),
     *          KVERT2(NNVE,*),KMID2(NNVE,*),KADJ2(NNVE,*)
C
      SAVE
C-----------------------------------------------------------------------
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C
      I1=KMID2(1,IELH1)-NVT2
      I2=KMID2(4,IELH2)-NVT2
      I3=KMID2(1,IELH2)-NVT2
      I4=KMID2(4,IELH3)-NVT2
      I5=KMID2(1,IELH3)-NVT2
      I6=KMID2(4,IELH4)-NVT2
      I7=KMID2(1,IELH4)-NVT2
      I8=KMID2(4,IELH1)-NVT2
      I9=KMID2(2,IELH1)-NVT2
      I10=KMID2(2,IELH2)-NVT2
      I11=KMID2(2,IELH3)-NVT2
      I12=KMID2(2,IELH4)-NVT2
C
      DUH1= DU2(I1)
      DUH2= DU2(I2)
      DUH3= DU2(I3)
      DUH4= DU2(I4)
      DUH5= DU2(I5)
      DUH6= DU2(I6)
      DUH7= DU2(I7)
      DUH8= DU2(I8)
      DUH9= DU2(I9)
      DUH10=DU2(I10)
      DUH11=DU2(I11)
      DUH12=DU2(I12)
C
      DVH1= DV2(I1)
      DVH2= DV2(I2)
      DVH3= DV2(I3)
      DVH4= DV2(I4)
      DVH5= DV2(I5)
      DVH6= DV2(I6)
      DVH7= DV2(I7)
      DVH8= DV2(I8)
      DVH9= DV2(I9)
      DVH10=DV2(I10)
      DVH11=DV2(I11)
      DVH12=DV2(I12)
C
C *** The edge IM1
C
      IF (KADJ1(1,IEL1).NE.0) THEN
C     case of an inner edge
        IF (KADJ1(1,IEL1).GT.IEL1) THEN
           DU1(IM1)= A1*(DUH1+DUH2)+A2*(DUH4+DUH7)
     *                 +A3*(DUH5+DUH6)+A4*(DUH3+DUH8)
     *                   +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
           DV1(IM1)= A1*(DVH1+DVH2)+A2*(DVH4+DVH7)
     *                 +A3*(DVH5+DVH6)+A4*(DVH3+DVH8)
     *                   +A5*DVH9+A6*(DVH10+DVH12)+A7*DVH11
        ELSE
           DU1(IM1)=DU1(IM1)+A2*(DUH4+DUH7)
     *                        +A3*(DUH5+DUH6)+A4*(DUH3+DUH8)
     *                          +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
           DV1(IM1)=DV1(IM1)+A2*(DVH4+DVH7)
     *                        +A3*(DVH5+DVH6)+A4*(DVH3+DVH8)
     *                          +A5*DVH9+A6*(DVH10+DVH12)+A7*DVH11
        ENDIF
      ELSE
       DU1(IM1)=     A1*(DUH1+DUH2)+2D0*A2*(DUH4+DUH7)
     *          +2D0*A3*(DUH5+DUH6)+2D0*A4*(DUH3+DUH8)
     *          +    A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
       DV1(IM1)=     A1*(DVH1+DVH2)+2D0*A2*(DVH4+DVH7)
     *          +2D0*A3*(DVH5+DVH6)+2D0*A4*(DVH3+DVH8)
     *          +    A5*DVH9+A6*(DVH10+DVH12)+A7*DVH11
      ENDIF
C
C *** The edge IM2
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(2,IEL1).GT.IEL1) THEN
           DU1(IM2)= A1*(DUH3+DUH4)+A2*(DUH6+DUH1)
     *                 +A3*(DUH7+DUH8)+A4*(DUH5+DUH2)
     *                   +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
           DV1(IM2)= A1*(DVH3+DVH4)+A2*(DVH6+DVH1)
     *                 +A3*(DVH7+DVH8)+A4*(DVH5+DVH2)
     *                   +A5*DVH10+A6*(DVH11+DVH9)+A7*DVH12
        ELSE
           DU1(IM2)=DU1(IM2)+A2*(DUH6+DUH1)
     *                        +A3*(DUH7+DUH8)+A4*(DUH5+DUH2)
     *                          +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
           DV1(IM2)=DV1(IM2)+A2*(DVH6+DVH1)
     *                        +A3*(DVH7+DVH8)+A4*(DVH5+DVH2)
     *                          +A5*DVH10+A6*(DVH11+DVH9)+A7*DVH12
        ENDIF
      ELSE
       DU1(IM2)=     A1*(DUH3+DUH4)+2D0*A2*(DUH6+DUH1)
     *          +2D0*A3*(DUH7+DUH8)+2D0*A4*(DUH5+DUH2)
     *          +    A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
       DV1(IM2)=     A1*(DVH3+DVH4)+2D0*A2*(DVH6+DVH1)
     *          +2D0*A3*(DVH7+DVH8)+2D0*A4*(DVH5+DVH2)
     *          +    A5*DVH10+A6*(DVH11+DVH9)+A7*DVH12
      ENDIF
C
C *** The edge IM3
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(3,IEL1).GT.IEL1) THEN
           DU1(IM3)= A1*(DUH5+DUH6)+A2*(DUH8+DUH3)
     *                 +A3*(DUH1+DUH2)+A4*(DUH7+DUH4)
     *                   +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
           DV1(IM3)= A1*(DVH5+DVH6)+A2*(DVH8+DVH3)
     *                 +A3*(DVH1+DVH2)+A4*(DVH7+DVH4)
     *                   +A5*DVH11+A6*(DVH12+DVH10)+A7*DVH9
        ELSE
           DU1(IM3)=DU1(IM3)+A2*(DUH8+DUH3)
     *                        +A3*(DUH1+DUH2)+A4*(DUH7+DUH4)
     *                          +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
           DV1(IM3)=DV1(IM3)+A2*(DVH8+DVH3)
     *                        +A3*(DVH1+DVH2)+A4*(DVH7+DVH4)
     *                          +A5*DVH11+A6*(DVH12+DVH10)+A7*DVH9
        ENDIF
      ELSE
       DU1(IM3)=     A1*(DUH5+DUH6)+2D0*A2*(DUH8+DUH3)
     *          +2D0*A3*(DUH1+DUH2)+2D0*A4*(DUH7+DUH4)
     *          +    A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
       DV1(IM3)=     A1*(DVH5+DVH6)+2D0*A2*(DVH8+DVH3)
     *          +2D0*A3*(DVH1+DVH2)+2D0*A4*(DVH7+DVH4)
     *          +    A5*DVH11+A6*(DVH12+DVH10)+A7*DVH9
      ENDIF
C
C *** The edge IM4
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(4,IEL1).GT.IEL1) THEN
           DU1(IM4)= A1*(DUH7+DUH8)+A2*(DUH2+DUH5)
     *                 +A3*(DUH3+DUH4)+A4*(DUH1+DUH6)
     *                   +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
           DV1(IM4)= A1*(DVH7+DVH8)+A2*(DVH2+DVH5)
     *                 +A3*(DVH3+DVH4)+A4*(DVH1+DVH6)
     *                   +A5*DVH12+A6*(DVH9+DVH11)+A7*DVH10
        ELSE
           DU1(IM4)=DU1(IM4)+A2*(DUH2+DUH5)
     *                        +A3*(DUH3+DUH4)+A4*(DUH1+DUH6)
     *                          +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
           DV1(IM4)=DV1(IM4)+A2*(DVH2+DVH5)
     *                        +A3*(DVH3+DVH4)+A4*(DVH1+DVH6)
     *                          +A5*DVH12+A6*(DVH9+DVH11)+A7*DVH10
        ENDIF
      ELSE
       DU1(IM4)=     A1*(DUH7+DUH8)+2D0*A2*(DUH2+DUH5)
     *          +2D0*A3*(DUH3+DUH4)+2D0*A4*(DUH1+DUH6)
     *          +    A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
       DV1(IM4)=     A1*(DVH7+DVH8)+2D0*A2*(DVH2+DVH5)
     *          +2D0*A3*(DVH3+DVH4)+2D0*A4*(DVH1+DVH6)
     *          +    A5*DVH12+A6*(DVH9+DVH11)+A7*DVH10
      ENDIF
C
C
10    CONTINUE
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE  RESTRU  (DU1,DV1,DU2,DV2,
     *                     KVERT1,KMID1,KADJ1,NEQ1,NEL1,NVT1,
     *                     KVERT2,KMID2,KADJ2,NEQ2,NEL2,NVT2 )
************************************************************************
*    Purpose:  Restricts the  fine grid solution vector (DU2,DV2) to
*              the coarse grid solution vector (DU1,DV1)
*
*    Remark:   RESTRU works for all cases of boundary conditions
*                
*-----------------------------------------------------------------------
*    Input:
*      DU2,DV2               - fine grid solution vector
*      KVERT1,KMID1,..,NVT1  - data of the coarse grid
*      KVERT2,KMID2,..,NVT2  - data of the fine grid
*
*    Output:
*      DU1,DV1               - coarse grid solution vector
*-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      PARAMETER (A1=0.1875D0, A2=0.375D0, A3=-0.0625D0)
      PARAMETER (R1=0.375D0, R2=0.75D0, R3=-0.125D0)
C
      DIMENSION DU1(*),DV1(*),  DU2(*),DV2(*),
     *          KVERT1(NNVE,*),KMID1(NNVE,*),KADJ1(NNVE,*),
     *          KVERT2(NNVE,*),KMID2(NNVE,*),KADJ2(NNVE,*)
      SAVE
C-----------------------------------------------------------------------
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
      I1=KMID2(1,IELH1)-NVT2
      I2=KMID2(4,IELH2)-NVT2
      I3=KMID2(1,IELH2)-NVT2
      I4=KMID2(4,IELH3)-NVT2
      I5=KMID2(1,IELH3)-NVT2
      I6=KMID2(4,IELH4)-NVT2
      I7=KMID2(1,IELH4)-NVT2
      I8=KMID2(4,IELH1)-NVT2
      I9=KMID2(2,IELH1)-NVT2
      I10=KMID2(2,IELH2)-NVT2
      I11=KMID2(2,IELH3)-NVT2
      I12=KMID2(2,IELH4)-NVT2
C
      DUH1= DU2(I1)
      DUH2= DU2(I2)
      DUH3= DU2(I3)
      DUH4= DU2(I4)
      DUH5= DU2(I5)
      DUH6= DU2(I6)
      DUH7= DU2(I7)
      DUH8= DU2(I8)
      DUH9= DU2(I9)
      DUH10=DU2(I10)
      DUH11=DU2(I11)
      DUH12=DU2(I12)
C
      DVH1= DV2(I1)
      DVH2= DV2(I2)
      DVH3= DV2(I3)
      DVH4= DV2(I4)
      DVH5= DV2(I5)
      DVH6= DV2(I6)
      DVH7= DV2(I7)
      DVH8= DV2(I8)
      DVH9= DV2(I9)
      DVH10=DV2(I10)
      DVH11=DV2(I11)
      DVH12=DV2(I12)
C
C *** The edge IM1
C
      IF (KADJ1(1,IEL1).NE.0) THEN
C     case of an inner edge
       IF (KADJ1(1,IEL1).GT.IEL1) THEN
        DU1(IM1)=A1*(DUH1+DUH2)  +A2*DUH9   +A3*(DUH8+DUH3+DUH10+DUH12)
        DV1(IM1)=A1*(DVH1+DVH2)  +A2*DVH9   +A3*(DVH8+DVH3+DVH10+DVH12)
       ELSE
        DU1(IM1)=DU1(IM1) +
     *           A1*(DUH1+DUH2)  +A2*DUH9   +A3*(DUH8+DUH3+DUH10+DUH12)
        DV1(IM1)=DV1(IM1) +
     *           A1*(DVH1+DVH2)  +A2*DVH9   +A3*(DVH8+DVH3+DVH10+DVH12)
       ENDIF
      ELSE
C     case of a boundary edge
        DU1(IM1)=R1*(DUH1+DUH2)  +R2*DUH9   +R3*(DUH8+DUH3+DUH10+DUH12)
        DV1(IM1)=R1*(DVH1+DVH2)  +R2*DVH9   +R3*(DVH8+DVH3+DVH10+DVH12)
      ENDIF
C
C *** The edge IM2
C
      IF (KADJ1(2,IEL1).NE.0) THEN
C     case of an inner edge
       IF (KADJ1(2,IEL1).GT.IEL1) THEN
        DU1(IM2)=A1*(DUH3+DUH4)  +A2*DUH10  +A3*(DUH2+DUH5+DUH9 +DUH11)
        DV1(IM2)=A1*(DVH3+DVH4)  +A2*DVH10  +A3*(DVH2+DVH5+DVH9 +DVH11)
       ELSE
        DU1(IM2)=DU1(IM2) +
     *           A1*(DUH3+DUH4)  +A2*DUH10  +A3*(DUH2+DUH5+DUH9 +DUH11)
        DV1(IM2)=DV1(IM2) +
     *           A1*(DVH3+DVH4)  +A2*DVH10  +A3*(DVH2+DVH5+DVH9 +DVH11)
       ENDIF
      ELSE
C     case of a boundary edge
        DU1(IM2)=R1*(DUH3+DUH4)  +R2*DUH10  +R3*(DUH2+DUH5+DUH9 +DUH11)
        DV1(IM2)=R1*(DVH3+DVH4)  +R2*DVH10  +R3*(DVH2+DVH5+DVH9 +DVH11)
      ENDIF
C
C *** The edge IM3
C
      IF (KADJ1(3,IEL1).NE.0) THEN
C     case of an inner edge
       IF (KADJ1(3,IEL1).GT.IEL1) THEN
        DU1(IM3)=A1*(DUH5+DUH6)  +A2*DUH11  +A3*(DUH4+DUH7+DUH10+DUH12)
        DV1(IM3)=A1*(DVH5+DVH6)  +A2*DVH11  +A3*(DVH4+DVH7+DVH10+DVH12)
       ELSE
        DU1(IM3)=DU1(IM3) +
     *           A1*(DUH5+DUH6)  +A2*DUH11  +A3*(DUH4+DUH7+DUH10+DUH12)
        DV1(IM3)=DV1(IM3) +
     *           A1*(DVH5+DVH6)  +A2*DVH11  +A3*(DVH4+DVH7+DVH10+DVH12)
       ENDIF
      ELSE
C     case of a boundary edge
        DU1(IM3)=R1*(DUH5+DUH6)  +R2*DUH11  +R3*(DUH4+DUH7+DUH10+DUH12)
        DV1(IM3)=R1*(DVH5+DVH6)  +R2*DVH11  +R3*(DVH4+DVH7+DVH10+DVH12)
      ENDIF
C
C *** The edge IM4
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
C     case of an inner edge
       IF (KADJ1(4,IEL1).GT.IEL1) THEN
        DU1(IM4)=A1*(DUH7+DUH8)  +A2*DUH12  +A3*(DUH6+DUH1+DUH9 +DUH11)
        DV1(IM4)=A1*(DVH7+DVH8)  +A2*DVH12  +A3*(DVH6+DVH1+DVH9 +DVH11)
       ELSE
        DU1(IM4)=DU1(IM4) +
     *           A1*(DUH7+DUH8)  +A2*DUH12  +A3*(DUH6+DUH1+DUH9 +DUH11)
        DV1(IM4)=DV1(IM4) +
     *           A1*(DVH7+DVH8)  +A2*DVH12  +A3*(DVH6+DVH1+DVH9 +DVH11)
       ENDIF
      ELSE
C     case of a boundary edge
        DU1(IM4)=R1*(DUH7+DUH8)  +R2*DUH12  +R3*(DUH6+DUH1+DUH9 +DUH11)
        DV1(IM4)=R1*(DVH7+DVH8)  +R2*DVH12  +R3*(DVH6+DVH1+DVH9 +DVH11)
      ENDIF
C
C
10    CONTINUE
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE MR010 (DP1,DP2,KADJ1,KADJ2,NEL1,NEL2)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NNVE=4)
      DIMENSION DP1(*),DP2(*)
      DIMENSION KADJ1(NNVE,*),KADJ2(NNVE,*)
C
      SAVE
C-----------------------------------------------------------------------
C
      DO 10 IEL1=1,NEL1
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C *** Restriction of pressure
C
      DP1(IEL1)= DP2(IELH1)+DP2(IELH2)+DP2(IELH3)+DP2(IELH4)
C
10    CONTINUE
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE MP030(DU1,DU2,KVERT1,KVERT2,KMID1,KMID2,
     *                 KADJ1,KADJ2,NVT1,NVT2,NEL1,NEL2,NMT2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=4)
      PARAMETER (A1=0.5D0,A2=-0.0625D0,A3=0D0,A4=0.0625D0)
      PARAMETER (A5=0.625D0,A6=0.125D0,A7=0.125D0,A8=0.125D0)
      DIMENSION DU1(*),DU2(*),KVERT1(NNVE,*),KVERT2(NNVE,*),
     *          KMID1(NNVE,*),KMID2(NNVE,*),KADJ1(NNVE,*),
     *          KADJ2(NNVE,*)
C
C
C
C *** Zero initialization of (DU2)
      CALL  LCL1 (DU2,NMT2)
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      DUH1=DU1(IM1)
      DUH2=DU1(IM2)
      DUH3=DU1(IM3)
      DUH4=DU1(IM4)
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C *** The edge IM1 and the corresponding fine inner node
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH1)-NVT2
       IB=KMID2(4,IELH2)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4
       DU2(IB)=DU2(IB)+   A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4
      ELSE
       IA=KMID2(1,IELH1)-NVT2
       IB=KMID2(4,IELH2)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4)
      ENDIF
      IC=KMID2(2,IELH1)-NVT2
      DU2(IC)=A5*DUH1+A6*(DUH2+DUH4)+A7*DUH3
C
C *** The edge IM2 and the corresponding fine inner node
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH2)-NVT2
       IB=KMID2(4,IELH3)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1
       DU2(IB)=DU2(IB)+   A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1
      ELSE
       IA=KMID2(1,IELH2)-NVT2
       IB=KMID2(4,IELH3)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1)
      ENDIF
      IC=KMID2(2,IELH2)-NVT2
      DU2(IC)=A5*DUH2+A6*(DUH3+DUH1)+A7*DUH4
C
C *** The edge IM3 and the corresponding fine inner node
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH3)-NVT2
       IB=KMID2(4,IELH4)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2
       DU2(IB)=DU2(IB)+   A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2
      ELSE
       IA=KMID2(1,IELH3)-NVT2
       IB=KMID2(4,IELH4)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2)
      ENDIF
      IC=KMID2(2,IELH3)-NVT2
      DU2(IC)=A5*DUH3+A6*(DUH4+DUH2)+A7*DUH1
C
C *** The edge IM4 and the corresponding fine inner node
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH4)-NVT2
       IB=KMID2(4,IELH1)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3
       DU2(IB)=DU2(IB)+   A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3
      ELSE
       IA=KMID2(1,IELH4)-NVT2
       IB=KMID2(4,IELH1)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3)
      ENDIF
      IC=KMID2(2,IELH4)-NVT2
      DU2(IC)=A5*DUH4+A6*(DUH1+DUH3)+A7*DUH2
C
C
10    CONTINUE
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE MP031(DU1,DU2,KVERT1,KVERT2,KMID1,KMID2,
     *                 KADJ1,KADJ2,NVT1,NVT2,NEL1,NEL2,NMT2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=4)
      PARAMETER (A1=0.46875D0,A2=-0.09375D0,A3=-0.03125D0,A4=0.15625D0)
      PARAMETER (A5=0.5625D0,A6=0.1875D0,A7=0.0625D0,A8=0.1875D0)
      DIMENSION DU1(*),DU2(*),KVERT1(NNVE,*),KVERT2(NNVE,*),
     *          KMID1(NNVE,*),KMID2(NNVE,*),KADJ1(NNVE,*),
     *          KADJ2(NNVE,*)
C
C
C
C *** Zero initialization of (DU2)
      CALL  LCL1 (DU2,NMT2)
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      DUH1=DU1(IM1)
      DUH2=DU1(IM2)
      DUH3=DU1(IM3)
      DUH4=DU1(IM4)
C
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C
C
C *** The edge IM1 and the corresponding fine inner node
C
      IF (KADJ1(1,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH1)-NVT2
       IB=KMID2(4,IELH2)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4
       DU2(IB)=DU2(IB)+   A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4
      ELSE
       IA=KMID2(1,IELH1)-NVT2
       IB=KMID2(4,IELH2)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH1+A2*DUH2+A3*DUH3+A4*DUH4)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH1+A4*DUH2+A3*DUH3+A2*DUH4)
      ENDIF
      IC=KMID2(2,IELH1)-NVT2
      DU2(IC)=A5*DUH1+A6*(DUH2+DUH4)+A7*DUH3
C
C *** The edge IM2 and the corresponding fine inner node
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH2)-NVT2
       IB=KMID2(4,IELH3)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1
       DU2(IB)=DU2(IB)+   A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1
      ELSE
       IA=KMID2(1,IELH2)-NVT2
       IB=KMID2(4,IELH3)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH2+A2*DUH3+A3*DUH4+A4*DUH1)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH2+A4*DUH3+A3*DUH4+A2*DUH1)
      ENDIF
      IC=KMID2(2,IELH2)-NVT2
      DU2(IC)=A5*DUH2+A6*(DUH3+DUH1)+A7*DUH4
C
C *** The edge IM3 and the corresponding fine inner node
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH3)-NVT2
       IB=KMID2(4,IELH4)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2
       DU2(IB)=DU2(IB)+   A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2
      ELSE
       IA=KMID2(1,IELH3)-NVT2
       IB=KMID2(4,IELH4)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH3+A2*DUH4+A3*DUH1+A4*DUH2)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH3+A4*DUH4+A3*DUH1+A2*DUH2)
      ENDIF
      IC=KMID2(2,IELH3)-NVT2
      DU2(IC)=A5*DUH3+A6*(DUH4+DUH2)+A7*DUH1
C
C *** The edge IM4 and the corresponding fine inner node
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
       IA=KMID2(1,IELH4)-NVT2
       IB=KMID2(4,IELH1)-NVT2
       DU2(IA)=DU2(IA)+   A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3
       DU2(IB)=DU2(IB)+   A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3
      ELSE
       IA=KMID2(1,IELH4)-NVT2
       IB=KMID2(4,IELH1)-NVT2
       DU2(IA)=DU2(IA)+2D0*(A1*DUH4+A2*DUH1+A3*DUH2+A4*DUH3)
       DU2(IB)=DU2(IB)+2D0*(A1*DUH4+A4*DUH1+A3*DUH2+A2*DUH3)
      ENDIF
      IC=KMID2(2,IELH4)-NVT2
      DU2(IC)=A5*DUH4+A6*(DUH1+DUH3)+A7*DUH2
C
C
10    CONTINUE
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE MR030(DU2,DU1,KVERT2,KVERT1,KMID2,KMID1,
     *                 KADJ2,KADJ1,NVT2,NVT1,NEL2,NEL1)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
	PARAMETER (NNVE=4)
      PARAMETER (A1=1D0,A2=-0.0625D0,A3=0D0,A4=0.0625D0)
      PARAMETER (A5=0.625D0,A6=0.125D0,A7=0.125D0,A8=0.125D0)
      DIMENSION DU1(1),DU2(1),KVERT2(NNVE,1),KVERT1(NNVE,1),
     *          KMID1(NNVE,1),KMID2(NNVE,1),
     *          KADJ1(NNVE,1),KADJ2(NNVE,1)
C
C
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
      I1=KMID2(1,IELH1)-NVT2
      I2=KMID2(4,IELH2)-NVT2
      I3=KMID2(1,IELH2)-NVT2
      I4=KMID2(4,IELH3)-NVT2
      I5=KMID2(1,IELH3)-NVT2
      I6=KMID2(4,IELH4)-NVT2
      I7=KMID2(1,IELH4)-NVT2
      I8=KMID2(4,IELH1)-NVT2
      I9=KMID2(2,IELH1)-NVT2
      I10=KMID2(2,IELH2)-NVT2
      I11=KMID2(2,IELH3)-NVT2
      I12=KMID2(2,IELH4)-NVT2
C
      DUH1= DU2(I1)
      DUH2= DU2(I2)
      DUH3= DU2(I3)
      DUH4= DU2(I4)
      DUH5= DU2(I5)
      DUH6= DU2(I6)
      DUH7= DU2(I7)
      DUH8= DU2(I8)
      DUH9= DU2(I9)
      DUH10=DU2(I10)
      DUH11=DU2(I11)
      DUH12=DU2(I12)
C
C
C *** The edge IM1
C
      IF (KADJ1(1,IEL1).NE.0) THEN
C     case of an inner edge
        IF (KADJ1(1,IEL1).GT.IEL1) THEN
           DU1(IM1)= A1*(DUH1+DUH2)+A2*(DUH4+DUH7)
     *                 +A3*(DUH5+DUH6)+A4*(DUH3+DUH8)
     *                   +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
        ELSE
           DU1(IM1)=DU1(IM1)+A2*(DUH4+DUH7)
     *                        +A3*(DUH5+DUH6)+A4*(DUH3+DUH8)
     *                          +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
        ENDIF
      ELSE
       DU1(IM1)=     A1*(DUH1+DUH2)+2D0*A2*(DUH4+DUH7)
     *          +2D0*A3*(DUH5+DUH6)+2D0*A4*(DUH3+DUH8)
     *          +    A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
      ENDIF
C
C *** The edge IM2
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(2,IEL1).GT.IEL1) THEN
           DU1(IM2)= A1*(DUH3+DUH4)+A2*(DUH6+DUH1)
     *                 +A3*(DUH7+DUH8)+A4*(DUH5+DUH2)
     *                   +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
        ELSE
           DU1(IM2)=DU1(IM2)+A2*(DUH6+DUH1)
     *                        +A3*(DUH7+DUH8)+A4*(DUH5+DUH2)
     *                          +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
        ENDIF
      ELSE
       DU1(IM2)=     A1*(DUH3+DUH4)+2D0*A2*(DUH6+DUH1)
     *          +2D0*A3*(DUH7+DUH8)+2D0*A4*(DUH5+DUH2)
     *          +    A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
      ENDIF
C
C *** The edge IM3
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(3,IEL1).GT.IEL1) THEN
           DU1(IM3)= A1*(DUH5+DUH6)+A2*(DUH8+DUH3)
     *                 +A3*(DUH1+DUH2)+A4*(DUH7+DUH4)
     *                   +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
        ELSE
           DU1(IM3)=DU1(IM3)+A2*(DUH8+DUH3)
     *                        +A3*(DUH1+DUH2)+A4*(DUH7+DUH4)
     *                          +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
        ENDIF
      ELSE
       DU1(IM3)=     A1*(DUH5+DUH6)+2D0*A2*(DUH8+DUH3)
     *          +2D0*A3*(DUH1+DUH2)+2D0*A4*(DUH7+DUH4)
     *          +    A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
      ENDIF
C
C *** The edge IM4
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(4,IEL1).GT.IEL1) THEN
           DU1(IM4)= A1*(DUH7+DUH8)+A2*(DUH2+DUH5)
     *                 +A3*(DUH3+DUH4)+A4*(DUH1+DUH6)
     *                   +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
        ELSE
           DU1(IM4)=DU1(IM4)+A2*(DUH2+DUH5)
     *                        +A3*(DUH3+DUH4)+A4*(DUH1+DUH6)
     *                          +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
        ENDIF
      ELSE
       DU1(IM4)=     A1*(DUH7+DUH8)+2D0*A2*(DUH2+DUH5)
     *          +2D0*A3*(DUH3+DUH4)+2D0*A4*(DUH1+DUH6)
     *          +    A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
      ENDIF
C
C
10    CONTINUE
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE MP010 (DP1,DP2,KADJ1,KADJ2,NEL1,NEL2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      DIMENSION DP1(*),DP2(*)
      DIMENSION KADJ1(NNVE,*),KADJ2(NNVE,*)
      SAVE
C-----------------------------------------------------------------------
C
      DO 10 IEL1=1,NEL1
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C *** Prolongation of pressure
C
      DPH=DP1(IEL1)
      DP2(IELH1)=DPH
      DP2(IELH2)=DPH
      DP2(IELH3)=DPH
      DP2(IELH4)=DPH
C
10    CONTINUE
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE MR031(DU2,DU1,KVERT2,KVERT1,KMID2,KMID1,
     *                 KADJ2,KADJ1,NVT2,NVT1,NEL2,NEL1)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=4)
      PARAMETER (A1=0.9375D0,A2=-0.09375D0,A3=-0.03125D0,A4=0.15625D0)
      PARAMETER (A5=0.5625D0,A6=0.1875D0,A7=0.0625D0,A8=0.1875D0)
      DIMENSION DU1(1),DU2(1),KVERT2(NNVE,1),KVERT1(NNVE,1),
     *          KMID1(NNVE,1),KMID2(NNVE,1),
     *          KADJ1(NNVE,1),KADJ2(NNVE,1)
C
C
C
      DO 10 IEL1=1,NEL1
C
      IM1=KMID1(1,IEL1)-NVT1
      IM2=KMID1(2,IEL1)-NVT1
      IM3=KMID1(3,IEL1)-NVT1
      IM4=KMID1(4,IEL1)-NVT1
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
C
C
C
      I1=KMID2(1,IELH1)-NVT2
      I2=KMID2(4,IELH2)-NVT2
      I3=KMID2(1,IELH2)-NVT2
      I4=KMID2(4,IELH3)-NVT2
      I5=KMID2(1,IELH3)-NVT2
      I6=KMID2(4,IELH4)-NVT2
      I7=KMID2(1,IELH4)-NVT2
      I8=KMID2(4,IELH1)-NVT2
      I9=KMID2(2,IELH1)-NVT2
      I10=KMID2(2,IELH2)-NVT2
      I11=KMID2(2,IELH3)-NVT2
      I12=KMID2(2,IELH4)-NVT2
C
      DUH1= DU2(I1)
      DUH2= DU2(I2)
      DUH3= DU2(I3)
      DUH4= DU2(I4)
      DUH5= DU2(I5)
      DUH6= DU2(I6)
      DUH7= DU2(I7)
      DUH8= DU2(I8)
      DUH9= DU2(I9)
      DUH10=DU2(I10)
      DUH11=DU2(I11)
      DUH12=DU2(I12)
C
C
C *** The edge IM1
C
      IF (KADJ1(1,IEL1).NE.0) THEN
C     case of an inner edge
        IF (KADJ1(1,IEL1).GT.IEL1) THEN
           DU1(IM1)= A1*(DUH1+DUH2)+A2*(DUH4+DUH7)
     *                 +A3*(DUH5+DUH6)+A4*(DUH3+DUH8)
     *                   +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
        ELSE
           DU1(IM1)=DU1(IM1)+A2*(DUH4+DUH7)
     *                        +A3*(DUH5+DUH6)+A4*(DUH3+DUH8)
     *                          +A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
        ENDIF
      ELSE
       DU1(IM1)=     A1*(DUH1+DUH2)+2D0*A2*(DUH4+DUH7)
     *          +2D0*A3*(DUH5+DUH6)+2D0*A4*(DUH3+DUH8)
     *          +    A5*DUH9+A6*(DUH10+DUH12)+A7*DUH11
      ENDIF
C
C *** The edge IM2
C
      IF (KADJ1(2,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(2,IEL1).GT.IEL1) THEN
           DU1(IM2)= A1*(DUH3+DUH4)+A2*(DUH6+DUH1)
     *                 +A3*(DUH7+DUH8)+A4*(DUH5+DUH2)
     *                   +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
        ELSE
           DU1(IM2)=DU1(IM2)+A2*(DUH6+DUH1)
     *                        +A3*(DUH7+DUH8)+A4*(DUH5+DUH2)
     *                          +A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
        ENDIF
      ELSE
       DU1(IM2)=     A1*(DUH3+DUH4)+2D0*A2*(DUH6+DUH1)
     *          +2D0*A3*(DUH7+DUH8)+2D0*A4*(DUH5+DUH2)
     *          +    A5*DUH10+A6*(DUH11+DUH9)+A7*DUH12
      ENDIF
C
C *** The edge IM3
C
      IF (KADJ1(3,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(3,IEL1).GT.IEL1) THEN
           DU1(IM3)= A1*(DUH5+DUH6)+A2*(DUH8+DUH3)
     *                 +A3*(DUH1+DUH2)+A4*(DUH7+DUH4)
     *                   +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
        ELSE
           DU1(IM3)=DU1(IM3)+A2*(DUH8+DUH3)
     *                        +A3*(DUH1+DUH2)+A4*(DUH7+DUH4)
     *                          +A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
        ENDIF
      ELSE
       DU1(IM3)=     A1*(DUH5+DUH6)+2D0*A2*(DUH8+DUH3)
     *          +2D0*A3*(DUH1+DUH2)+2D0*A4*(DUH7+DUH4)
     *          +    A5*DUH11+A6*(DUH12+DUH10)+A7*DUH9
      ENDIF
C
C *** The edge IM4
C
      IF (KADJ1(4,IEL1).NE.0) THEN 
C     case of an inner edge
        IF (KADJ1(4,IEL1).GT.IEL1) THEN
           DU1(IM4)= A1*(DUH7+DUH8)+A2*(DUH2+DUH5)
     *                 +A3*(DUH3+DUH4)+A4*(DUH1+DUH6)
     *                   +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
        ELSE
           DU1(IM4)=DU1(IM4)+A2*(DUH2+DUH5)
     *                        +A3*(DUH3+DUH4)+A4*(DUH1+DUH6)
     *                          +A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
        ENDIF
      ELSE
       DU1(IM4)=     A1*(DUH7+DUH8)+2D0*A2*(DUH2+DUH5)
     *          +2D0*A3*(DUH3+DUH4)+2D0*A4*(DUH1+DUH6)
     *          +    A5*DUH12+A6*(DUH9+DUH11)+A7*DUH10
      ENDIF
C
C
10    CONTINUE
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE  YAX (DX,DAX,NEQ,A1,A2)  
************************************************************************
*
*   Purpose: - performs the matrix-vector-operation
*
*                   DAX:= A1*(A*DX) + A2*DAX
*
*              of dimension NEQ   (A1,A2 given scalar variables)
*
*            - DX,DAX  have the structure  D=(D1,D2,DP)
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
C *** global constants
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
      DIMENSION DX(*),DAX(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
C-----------------------------------------------------------------------
C
      EXTERNAL MATMUL
C
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C=======================================================================
C     Getting all parameters for MATMUL
C=======================================================================
C *** addresses for the current level ILEV
C
      ISETLV=2
      CALL  SETLEV (ISETLV)
      I2=1+NU
      IP=I2+NU
C
      KVERT=L(LVERT)
      KMID =L(LMID )
      KADJ =L(LADJ )
      KNPR =L(LNPR )
      KMBD =L(KLMBD(ILEV))
      NMBD= KNMBD(ILEV)
C
      KAREA=L(KLAREA(ILEV))
C
C=======================================================================
C
      CALL MATMUL(DAX(1),DAX(I2),DAX(IP),DX(1),DX(I2),DX(IP),A1,A2,
     *            VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *            VWORK(KB1),VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *            NU,NP,KWORK(KMBD),NMBD,INEUM)
C
C
ccc      return
c
ccc       IF (IGRAD.EQ.0) THEN
       IF (ILEV.EQ.NLEV) THEN
        CALL ZNEW (2*NU ,1,LUN,'DUN   ')
        CALL ZNEW (2*NVT,1,LUL,'DUL   ')
        CALL LCP1(DX(1),DWORK(L(LUN)),2*NU)
C
ccc        write(6,*) 'YAX-DX-C',ILEV
        CALL FILT1C(DX(1),DX(I2),
     *              DWORK(L(LUN)),DWORK(L(LUN)+NU),
     *              DWORK(L(LUL)),DWORK(L(LUL)+NVT),
     *              KWORK(KVERT),KWORK(KMID),KWORK(KADJ),
     *              KWORK(KMBD),NMBD,NU,NVT,NEL)
C
        CALL ZDISP(0,LUL,'DUL   ')
        CALL ZDISP(0,LUN,'DUN   ')
       ENDIF
C
ccc       IF (IGRAD.EQ.0) THEN
       IF (ILEV.EQ.NLEV) THEN
        CALL ZNEW (2*NU ,1,LUN,'DUN   ')
        CALL ZNEW (2*NVT,1,LUL,'DUL   ')
        CALL LCP1(DAX(1),DWORK(L(LUN)),2*NU)
C
ccc        write(6,*) 'YAX-DAX-C',ILEV
        CALL FILT1C(DAX(1),DAX(I2),
     *              DWORK(L(LUN)),DWORK(L(LUN)+NU),
     *              DWORK(L(LUL)),DWORK(L(LUL)+NVT),
     *              KWORK(KVERT),KWORK(KMID),KWORK(KADJ),
     *              KWORK(KMBD),NMBD,NU,NVT,NEL)
C
        CALL ZDISP(0,LUL,'DUL   ')
        CALL ZDISP(0,LUN,'DUN   ')
       ENDIF
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE  YDBC (DX,NEQ)  
************************************************************************
*
*   Purpose: - sets Dirichlet boundary components of DX=(D1,D2) to 0
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
C
      KMBD =L(KLMBD(ILEV))
      NMBD= KNMBD(ILEV)
C
      CALL  BDRY0 (DX(1),DX(I2),KWORK(KMBD),NMBD)
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE   YEX (DX,DB,DD,NEQ,RHO)  
************************************************************************
*
*   Purpose: - computes on level ILEV on DX the solution of
*              
*                         A*DX=DB  
*               
*              with a certain accuracy controled by variables on
*              COMMON blocks 
*  
*            - DX,DB,DD have the structure  D=(D1,D2,DP)
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
C
      DIMENSION DX(*),DB(*),DD(*)
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE
C
C=======================================================================
C     Getting all parameters for SMOOTH
C=======================================================================
C
C *** addresses for the current level ILEV
      ISETLV=2
      CALL  SETLEV (ISETLV)
      I2=1+NU
      IP=I2+NU
C
      KVERT=L(LVERT)
      KMID =L(LMID )
      KNPR =L(LNPR )
      KMBD =L(KLMBD(ILEV))
      NMBD= KNMBD(ILEV)
C
      KAREA=L(KLAREA(ILEV))
C
C=======================================================================
C
      IF (ISL.EQ.1) THEN
C
       DO 11 ITE=1,NSL
C
       CALL  VANCAE (DX(1),DX(I2),DX(IP),DD(1),DD(I2),DD(IP),
     *               DB(1),DB(I2),DB(IP),
     *               VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *               VWORK(KB1),VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *               NU,NP,KWORK(KMBD),KWORK(KVERT),
     *               KWORK(KMID),KWORK(KNPR),NMBD,RLXSL,DMAXU,DMAXP)
C
       IF (INEUM.EQ.0) CALL TOL20A(DX(IP),VWORK(KAREA),NP,INEUM)
C
       RES=MAX(DMAXU,DMAXP)
       IF (ITE.EQ.1) RESINI=RES
       RHO=(RES/RESINI)**(1D0/DBLE(ITE))
c       write(6,*) 'VANCAE1',ite,DMAXU,DMAXP
       IF ((RES.LT.EPSSL).AND.(RES.LT.DMPSL*RESINI)) GOTO 99999
C
11     CONTINUE
C
      ENDIF
C
C
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE   YEXA (DX,DB,DD,NEQ,RHO,EPS,ITE)  
************************************************************************
*
*   Purpose: - computes on level ILEV on DX the solution of
*              
*                         A*DX=DB  
*               
*              with a certain accuracy controled by variables on
*              COMMON blocks 
*  
*            - DX,DB,DD have the structure  D=(D1,D2,DP)
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
C
      DIMENSION DX(*),DB(*),DD(*)
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE
C
C=======================================================================
C     Getting all parameters for SMOOTH
C=======================================================================
C
C *** addresses for the current level ILEV
      ISETLV=2
      CALL  SETLEV (ISETLV)
      I2=1+NU
      IP=I2+NU
C
      KVERT=L(LVERT)
      KMID =L(LMID )
      KNPR =L(LNPR )
      KMBD =L(KLMBD(ILEV))
      NMBD= KNMBD(ILEV)
C
      KAREA=L(KLAREA(ILEV))
C
C=======================================================================
C
      IF (ISL.EQ.1) THEN
C
       DO 11 ITE=1,NSL
C
       CALL  VANCAE (DX(1),DX(I2),DX(IP),DD(1),DD(I2),DD(IP),
     *               DB(1),DB(I2),DB(IP),
     *               VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *               VWORK(KB1),VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *               NU,NP,KWORK(KMBD),KWORK(KVERT),
     *               KWORK(KMID),KWORK(KNPR),NMBD,RLXSL,DMAXU,DMAXP)
C
       IF (INEUM.EQ.0) CALL TOL20A(DX(IP),VWORK(KAREA),NP,INEUM)
C
ccc       IF (IGRAD.EQ.0) THEN
        CALL ZNEW (2*NU ,1,LUN,'DUN   ')
        CALL ZNEW (2*NVT,1,LUL,'DUL   ')
        CALL LCP1(DX(1),DWORK(L(LUN)),2*NU)
C
ccc        write(6,*) 'YSM',ILEV,ITE
        CALL FILT1C(DX(1),DX(I2),
     *              DWORK(L(LUN)),DWORK(L(LUN)+NU),
     *              DWORK(L(LUL)),DWORK(L(LUL)+NVT),
     *              KWORK(KVERT),KWORK(KMID),KWORK(KADJ),
     *              KWORK(KMBD),NMBD,NU,NVT,NEL)
C
        CALL ZDISP(0,LUL,'DUL   ')
        CALL ZDISP(0,LUN,'DUN   ')
ccc       ENDIF
       RES=MAX(DMAXU,DMAXP)
       IF (ITE.EQ.1) RESINI=RES
       RHO=(RES/RESINI)**(1D0/DBLE(ITE))
       IF ((RES.LT.EPSMG).AND.(RES.LT.DMPMG*RESINI)) GOTO 99999
C
11     CONTINUE
C
      ENDIF
C
C
C
99999 END
C
C
C
************************************************************************
      SUBROUTINE  YPROL (DUC,DUF)  
************************************************************************
*
*   Purpose: - performs the prolongation   DUF:=p(DUC)
*              with
*                  DUF   - fine correction vector on level ILEV
*                  DUC   - coarse correction vector on level ILEV-1
*  
*            - DUF and DUC have the structure  DU=(DU1,DU2,DP)
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
C
      DIMENSION DUF(*),DUC(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
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
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE 
C
C=======================================================================
C     Getting all parameters for PROLU
C=======================================================================
C
C *** addresses for the fine level ILEV
      ISETLV=2
      CALL  SETLEV (ISETLV)
      I2F=1+NU
      IPF=I2F+NU
C
      KVERTF=L(LVERT)
      KMIDF =L(LMID )
      KADJF =L(LADJ )
C
C *** addresses for the coarse level ILEV-1
      I1=ILEV-1
      NUC=KNU(I1)
      NPC=KNP(I1)
      NVTC=KNVT(I1)
      NMTC=KNMT(I1)
      I2C=1+NUC
      IPC=I2C+NUC
C
      KVERTC=L(KLVERT(I1))
      KMIDC =L(KLMID (I1))
      KADJC =L(KLADJ (I1))
C
C=======================================================================
C
      CALL PROLU (DUC(1),DUC(I2C),DUC(IPC),DUF(1),DUF(I2F),DUF(IPF),
     *            KWORK(KVERTC),KWORK(KMIDC),KWORK(KADJC),NUC,NPC,NVTC,
     *            KWORK(KVERTF),KWORK(KMIDF),KWORK(KADJF),NU,NP,NVT)
C
      IF (ABS(IINT).EQ.2) THEN
       KPLC=L(LD1)
       KPLF=L(LD2)
C
       CALL C2N2DM(DUC(IPC),DWORK(KPLC),KWORK(KMIDC),
     *             KWORK(KADJC),NPC,NMTC,NVTC,0)
C
       IF (IINT.GT.0) THEN
        CALL MP031(DWORK(KPLC),DWORK(KPLF),KWORK(KVERTC),KWORK(KVERTF),
     *             KWORK(KMIDC),KWORK(KMIDF),KWORK(KADJC),KWORK(KADJF),
     *             NVTC,NVT,NPC,NP,NMT)
       ELSE
        CALL MP030(DWORK(KPLC),DWORK(KPLF),KWORK(KVERTC),KWORK(KVERTF),
     *             KWORK(KMIDC),KWORK(KMIDF),KWORK(KADJC),KWORK(KADJF),
     *             NVTC,NVT,NPC,NP,NMT)
       ENDIF
C
       CALL C2N2DM(DUF(IPF),DWORK(KPLF),KWORK(KMIDF),
     *             KWORK(KADJF),NP,NMT,NVT,1)
      ENDIF
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE   YREST (DDF,DDC)  
************************************************************************
*
*   Purpose: - performs the defect restriction   DDC:=r(DDF)
*              with
*                  DDF - fine defect vector on level ILEV+1
*                  DDC - coarse defect vector on level ILEV
*  
*            - DDF and DDC have the structure  DD=(DD1,DD2,DDP)
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
C
      DIMENSION DDF(*),DDC(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
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
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
C=======================================================================
C     Getting all parameters for RESTRD
C=======================================================================
C
C *** addresses for the coarse level ILEV
      ISETLV=2
      CALL  SETLEV (ISETLV)
      I2C=1+NU
      IPC=I2C+NU
C
      KVERTC=L(LVERT)
      KMIDC =L(LMID )
      KADJC =L(LADJ )
C
C *** addresses for the fine level ILEV+1
      I1=ILEV+1
      NUF=KNU(I1)
      NPF=KNP(I1)
      NVTF=KNVT(I1)
      NMTF=KNMT(I1)
      I2F=1+NUF
      IPF=I2F+NUF
C
      KVERTF=L(KLVERT(I1))
      KMIDF =L(KLMID (I1))
      KADJF =L(KLADJ (I1))
C
C=======================================================================
C
      CALL RESTRD(DDC(1),DDC(I2C),DDC(IPC),DDF(1),DDF(I2F),DDF(IPF),
     *            KWORK(KVERTC),KWORK(KMIDC),KWORK(KADJC),NU,NP,NVT,
     *            KWORK(KVERTF),KWORK(KMIDF),KWORK(KADJF),NUF,NPF,NVTF)
C
      IF (ABS(IINT).EQ.2) THEN
       KPLC=L(LD1)
       KPLF=L(LD2)
C
       CALL C2N2DM(DDF(IPF),DWORK(KPLF),KWORK(KMIDF),KWORK(KADJF),
     *             NPF,NMTF,NVTF,0)
C
       IF (IINT.GT.0) THEN
        CALL MR031(DWORK(KPLF),DWORK(KPLC),KWORK(KVERTF),KWORK(KVERTC),
     *             KWORK(KMIDF),KWORK(KMIDC),KWORK(KADJF),KWORK(KADJC),
     *             NVTF,NVT,NPF,NP)
       ELSE
        CALL MR030(DWORK(KPLF),DWORK(KPLC),KWORK(KVERTF),KWORK(KVERTC),
     *             KWORK(KMIDF),KWORK(KMIDC),KWORK(KADJF),KWORK(KADJC),
     *             NVTF,NVT,NPF,NP)
       ENDIF
C
       CALL C2N2DM(DDC(IPC),DWORK(KPLC),KWORK(KMIDC),KWORK(KADJC),
     *             NP,NMT,NVT,1)
      ENDIF
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE  YSM (DX,DB,DD,NEQ,NSMO)  
************************************************************************
*
*   Purpose: - performs NSMO smoothing steps applied to the system
*                          A*DX = DB
*              of dimension NEQ using the auxiliary vector DD
*
*            - DX,DB,DD have the structure  D=(D1,D2,DP)
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
C
      DIMENSION DX(*),DB(*),DD(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE
C
C=======================================================================
C     Getting all parameters for SMOOTH
C=======================================================================
C
C *** addresses for the current level ILEV
      ISETLV=2
      CALL  SETLEV (ISETLV)
      I2=1+NU
      IP=I2+NU
C
      KVERT=L(LVERT)
      KMID =L(LMID )
      KADJ =L(LADJ )
      KNPR =L(LNPR )
      KMBD =L(KLMBD(ILEV))
      NMBD= KNMBD(ILEV)
C
      KAREA=L(KLAREA(ILEV))
C
C=======================================================================
C
      IF (ISM.EQ.1) THEN
       DO 11  ITE=1,NSMO
       CALL VANCAS (DX(1),DX(I2),DX(IP),DD(1),DD(I2),DD(IP),
     *              DB(1),DB(I2),DB(IP),
     *              VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *              VWORK(KB1),VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *              NU,NP,KWORK(KMBD),KWORK(KVERT),KWORK(KMID),
     *              KWORK(KNPR),NMBD)
C
       IF (INEUM.EQ.0) CALL TOL20A(DX(IP),VWORK(KAREA),NP,INEUM)
C
ccc       IF (IGRAD.EQ.0) THEN
        CALL ZNEW (2*NU ,1,LUN,'DUN   ')
        CALL ZNEW (2*NVT,1,LUL,'DUL   ')
        CALL LCP1(DX(1),DWORK(L(LUN)),2*NU)
C
ccc        write(6,*) 'YSM',ILEV,ITE
        CALL FILT1C(DX(1),DX(I2),
     *              DWORK(L(LUN)),DWORK(L(LUN)+NU),
     *              DWORK(L(LUL)),DWORK(L(LUL)+NVT),
     *              KWORK(KVERT),KWORK(KMID),KWORK(KADJ),
     *              KWORK(KMBD),NMBD,NU,NVT,NEL)
C
        CALL ZDISP(0,LUL,'DUL   ')
        CALL ZDISP(0,LUN,'DUN   ')
ccc       ENDIF
C
11     CONTINUE
      ENDIF
C
       return
c
ccc       IF (IGRAD.EQ.0) THEN
       IF (ILEV.EQ.NLEV) THEN
        CALL ZNEW (2*NU ,1,LUN,'DUN   ')
        CALL ZNEW (2*NVT,1,LUL,'DUL   ')
        CALL LCP1(DX(1),DWORK(L(LUN)),2*NU)
C
ccc        write(6,*) 'YSM',ILEV
        CALL FILT1C(DX(1),DX(I2),
     *              DWORK(L(LUN)),DWORK(L(LUN)+NU),
     *              DWORK(L(LUL)),DWORK(L(LUL)+NVT),
     *              KWORK(KVERT),KWORK(KMID),KWORK(KADJ),
     *              KWORK(KMBD),NMBD,NU,NVT,NEL)
C
        CALL ZDISP(0,LUL,'DUL   ')
        CALL ZDISP(0,LUN,'DUN   ')
       ENDIF
C
C
C
      END
C
C
C
************************************************************************
      SUBROUTINE   YSTEP (DX,DD,DB,ALPHA)  
************************************************************************
*
*   Purpose: - performs step size control for prolongation with
*                  DX    - old fine solution vector on level ILEV
*                  DD    - fine correction vector on level ILEV
*                  DB    - fine right hand side vector on level ILEV
*                  ALPHA - relaxation parameter according to some
*                          optimization criterion but within the limits
*                          AMINMG and AMAXMG (COMMON /RPARM/)
*  
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299, NNLEV=9, NNWORK=1)
C
      DIMENSION DX(*),DD(*),DB(*)
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,IMASS,IMASSL,
     *               IUPW,IPRECA,IPRECB,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYC,ILMIN,ILMAX,IINT,ISM,ISL,
     *               NSM,NSL,NSMFAC,IGRAD
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,EPSD,EPSDIV,
     *                EPSUR,EPSPR,DMPD,DMPMG,EPSMG,DMPSL,EPSSL,
     *                RLXSM,RLXSL,AMINMG,AMAXMG
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
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
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SAVE 
C
C=======================================================================
C     Getting all parameters for PROLU
C=======================================================================
C
C
C *** Setting of given ALPHA
      IF (AMINMG.EQ.AMAXMG) THEN
       ALPHA=AMINMG
       RETURN
      ENDIF
C
C
C *** Calculation of optimal ALPHA
      IF (AMINMG.NE.AMAXMG) THEN
C
       ISETLV=2
       CALL SETLEV (ISETLV)
C
       CALL LCP1 (DB,DWORK(L(LD1)),NUP)
       CALL YAX(DX,DWORK(L(LD1)),NUP,-1D0,1D0)
       DO 100 IEQ=NUP,NUP-NEL+1,-1
100    DWORK(L(LD1)+IEQ-1)=-DWORK(L(LD1)+IEQ-1)
       CALL LSP1(DD,DWORK(L(LD1)),NUP,DBY)
C
C
       CALL YAX(DD,DWORK(L(LD1)),NUP,1D0,0D0)
       DO 110 IEQ=NUP,NUP-NEL+1,-1
110    DWORK(L(LD1)+IEQ-1)=-DWORK(L(LD1)+IEQ-1)
       CALL LSP1(DD,DWORK(L(LD1)),NUP,DBX)
C
       ALPHA=DBY/DBX
       IF (ALPHA.LT.AMINMG) ALPHA=AMINMG
       IF (ALPHA.GT.AMAXMG) ALPHA=AMAXMG
ccc       WRITE(*,*)'ALPHA=====',ALPHA,DBY,DBX,ILEV
C
       RETURN
C
      ENDIF
C
C
      END

