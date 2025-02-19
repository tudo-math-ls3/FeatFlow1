***********************************************************************
      SUBROUTINE TDEF(MFILE,MSHOW,BSTOP,BNLEND)
***********************************************************************      
C
C	Solves the linear system for the temperature
C	via multigrid iteration
C
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'bouss.inc'
C-----------------------------------------------------------------------
C
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C     E X T E R N A L S
C-----------------------------------------------------------------------
C *** Coefficient of stiffness matrix
      EXTERNAL COEFFN
C *** definition of finite elements
      EXTERNAL E030,E031,EM30,EM31
C *** Multigrid components
      EXTERNAL  YAXU,YPROLU,YRESTU,YSMU,YEXU,YEXAU,YDBCU,YSTEPU
C
C=======================================================================
C     First generation of the nonlinear block A on level NLMAX
C=======================================================================
C
      CALL ZTIME(TTDEF0)
      BSTOP =.FALSE.
      BNLEND=.FALSE.
C
      ISETLV=2
      ILEV=NLMAX
      CALL  SETLEV (ISETLV)
C
      IF (INLMAX.GT.1) THEN
       IDEFUP=1
      ELSE
       IDEFUP=0
      ENDIF
c
      IF (INLMAX.GT.1) THEN
       CALL LCP1 (DWORK(KFT),DWORK(L(LDT)),NU)
      ENDIF
C
      FAKTOR=ALFA/NY
      THSTEP=THSTEP*FAKTOR
      CALL XMDF2T(KM1,KMASS1,KST1,KA1,KCOLA,KLDA,L(LDT),KT,
     *            NA,NU,KWORK(L(KLMBD(NLEV))),KNMBD(NLEV),THSTEP)
      THSTEP=THSTEP/FAKTOR
c
C
C
      IF ((ISTOK.NE.1).OR.((ISTOK.EQ.1).AND.(IPRECA.EQ.4))) THEN
C
C
        IF (IUPW.EQ.1) THEN
            if (IPRECA.eq.4) then
              IF (IELT.EQ.0) 
     *          CALL TADSTP(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,DWORK(KT),DWORK(KT),
     *          DWORK(L(LDT)),DWORK(L(LDT)),VWORK(KA1),NA,KWORK(KCOLA),
     *          KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *          DWORK(L(LCORVG)),E031,COEFNN,IDEFUP,1,1D0)
              IF (IELT.EQ.1) 
     *          CALL TADSTP(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,DWORK(KT),DWORK(KT),
     *          DWORK(L(LDT)),DWORK(L(LDT)),VWORK(KA1),NA,KWORK(KCOLA),
     *          KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *          DWORK(L(LCORVG)),E030,COEFNN,IDEFUP,1,1D0)
              IF (IELT.EQ.2) 
     *          CALL TADSTN(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,DWORK(KT),DWORK(KT),
     *          DWORK(L(LDT)),DWORK(L(LDT)),VWORK(KA1),NA,KWORK(KCOLA),
     *          KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *          DWORK(L(LCORVG)),EM31,COEFNN,IDEFUP,1,1D0)
              IF (IELT.EQ.3) 
     *          CALL TADSTN(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,DWORK(KT),DWORK(KT),
     *          DWORK(L(LDT)),DWORK(L(LDT)),VWORK(KA1),NA,KWORK(KCOLA),
     *          KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *          DWORK(L(LCORVG)),EM30,COEFNN,IDEFUP,1,1D0)
            endif
            if (ISTOK.NE.1)
     *          CALL TGUPWD(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,
     *              DWORK(KT),DWORK(KT),
     *              DWORK(L(LDT)),DWORK(L(LDT)),
     *              VWORK(KA1),NA,KWORK(KCOLA),
     *              KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *              DWORK(L(LCORVG)),IDEFUP,1)
        ELSE
         IF (IELT.EQ.0) 
     *    CALL TUPWDG(DWORK(KU1),DWORK(KU2),
     *                DWORK(KU1),DWORK(KU2),1.d0,0.d0,
     *                DWORK(KT),DWORK(KT),
     *                DWORK(L(LDT)),DWORK(L(LDT)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                E031,COEFFN,IDEFUP,1,1D0)
         IF (IELT.EQ.1) 
     *    CALL TUPWDG(DWORK(KU1),DWORK(KU2),
     *                DWORK(KU1),DWORK(KU2),1.d0,0.d0,
     *                DWORK(KT),DWORK(KT),
     *                DWORK(L(LDT)),DWORK(L(LDT)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                E030,COEFFN,IDEFUP,1,1D0)
         IF (IELT.EQ.2) 
     *    CALL TUPWNP(DWORK(KU1),DWORK(KU2),
     *                DWORK(KU1),DWORK(KU2),1.d0,0.d0,
     *                DWORK(KT),DWORK(KT),
     *                DWORK(L(LDT)),DWORK(L(LDT)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                EM31,COEFFN,IDEFUP,1,1D0)
         IF (IELT.EQ.3) then
	CALL ZTIME (t3)
          CALL TUPWNP(DWORK(KU1),DWORK(KU2),
     *                DWORK(KU1),DWORK(KU2),1.d0,0.d0,
     *                DWORK(KT),DWORK(KT),
     *                DWORK(L(LDT)),DWORK(L(LDT)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                EM30,COEFFN,IDEFUP,1,1D0)
        ENDIF
        ENDIF
C
       ELSE
        IF (IUPW.EQ.1) THEN
            if (IPRECA.eq.4) then
              IF (IELT.EQ.0) 
     *          CALL TADSTP(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,DWORK(KT),DWORK(KT),
     *          DWORK(L(LDT)),DWORK(L(LDT)),VWORK(KA1),NA,KWORK(KCOLA),
     *          KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *          DWORK(L(LCORVG)),E031,COEFNN,IDEFUP,1,1D0)
              IF (IELT.EQ.1) 
     *          CALL TADSTP(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,DWORK(KT),DWORK(KT),
     *          DWORK(L(LDT)),DWORK(L(LDT)),VWORK(KA1),NA,KWORK(KCOLA),
     *          KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *          DWORK(L(LCORVG)),E030,COEFNN,IDEFUP,1,1D0)
              IF (IELT.EQ.2) 
     *          CALL TADSTN(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,DWORK(KT),DWORK(KT),
     *          DWORK(L(LDT)),DWORK(L(LDT)),VWORK(KA1),NA,KWORK(KCOLA),
     *          KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *          DWORK(L(LCORVG)),EM31,COEFNN,IDEFUP,1,1D0)
              IF (IELT.EQ.3) 
     *          CALL TADSTN(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,DWORK(KT),DWORK(KT),
     *          DWORK(L(LDT)),DWORK(L(LDT)),VWORK(KA1),NA,KWORK(KCOLA),
     *          KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *          DWORK(L(LCORVG)),EM30,COEFNN,IDEFUP,1,1D0)
            endif
            if (ISTOK.NE.1)
     *          CALL TGUPWD(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,
     *              DWORK(KT),DWORK(KT),
     *              DWORK(L(LDT)),DWORK(L(LDT)),
     *              VWORK(KA1),NA,KWORK(KCOLA),
     *              KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *              DWORK(L(LCORVG)),IDEFUP,1)
        ELSE
         IF (IELT.EQ.0) 
     *    CALL TUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KT),DWORK(KT),
     *                DWORK(L(LDT)),DWORK(L(LDT)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                E031,COEFFN,IDEFUP,1,1D0)
         IF (IELT.EQ.1) 
     *    CALL TUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KU1),DWORK(KU2),
     *                DWORK(L(LDT)),DWORK(L(LDT)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                E030,COEFFN,IDEFUP,1,1D0)
         IF (IELT.EQ.2) 
     *    CALL TUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KT),DWORK(KT),
     *                DWORK(L(LDT)),DWORK(L(LDT)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                EM31,COEFFN,IDEFUP,1,1D0)
         IF (IELT.EQ.3) 
     *    CALL TUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DWORK(KT),DWORK(KT),
     *                DWORK(L(LDT)),DWORK(L(LDT)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                EM30,COEFFN,IDEFUP,1,1D0)
        ENDIF
       ENDIF
C
      CALL BDRYAT (VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *            KWORK(L(KLMBD(ILEV))),KNMBD(ILEV),
     *            KWORK(L(LNPR)),DWORK(L(KLDBD(ILEV))))
C
      IF (INLMAX.GT.1) THEN
       CALL BDRY0T (DWORK(L(LDT)),KWORK(L(KLMBD(ILEV))),
     *             KNMBD(ILEV),
     *            KWORK(L(LNPR)),DWORK(L(KLDBD(ILEV))))
      ENDIF
c
C=======================================================================
C	The vector LDT now contains the defect,
C	KA1 the full matrix for the T-equation! 	
C=======================================================================
C
C=======================================================================
C     Calculation of initial defects
C=======================================================================
C
      IF (INLMAX.GT.1) THEN
       CALL  RESDFK(DWORK(L(LDT)),DWORK(L(LDT)),DWORK(KFT),DWORK(KFT),
     *              NU,REST,RDUMMY)
       RESOLD=REST
       RES0=RESOLD
       RES =RESOLD
       EPSRES=DMPUD*RESOLD
      ELSE
       RES   =0D0
       REST  =0D0
       RESOLD=0D0
       EPSRES=0D0
      ENDIF
C
C=======================================================================
C *** Generate the A blocks for all coarse levels
C=======================================================================
C
      IF (NLMAX.GT.NLMIN)  THEN
        DO 22  ILEV=NLMAX-1,NLMIN,-1
        CALL ZTIME(TTT0)
        ISETLV=2
        CALL  SETLEV (ISETLV)
C
        I1=ILEV+1
        KU1F=L(KLUP(I1))
        KU2F=KU1F+KNU(I1)
        KVERTF=L(KLVERT(I1))
        KMIDF =L(KLMID (I1))
        KADJF =L(KLADJ (I1))
        CALL  RESTRU (DWORK(KU1),DWORK(KU2), DWORK(KU1F),DWORK(KU2F),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *                NU,NP,NVT,KWORK(KVERTF),KWORK(KMIDF),KWORK(KADJF),
     *                KNU(I1),KNP(I1),KNVT(I1),2)
c
        KTF=L(KLUP(I1))+2*KNU(I1)+KNP(I1)
c
        CALL  RESTRU (DWORK(KT),DWORK(KT), DWORK(KTF),DWORK(KTF),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *                NU,NP,NVT,KWORK(KVERTF),KWORK(KMIDF),KWORK(KADJF),
     *                KNU(I1),KNP(I1),KNVT(I1),1)
C
C
        FAKTOR=ALFA/NY
        THSTEP=THSTEP*FAKTOR
        CALL XMADF3(KM1,KMASS1,KST1,KA1,KCOLA,KLDA,NA,NU,THSTEP)
        THSTEP=THSTEP/FAKTOR
C
        IF ((ISTOK.NE.1).OR.((ISTOK.EQ.1).AND.(IPRECA.EQ.4))) THEN
         IF (IUPW.EQ.1) THEN
            if (IPRECA.eq.4) then
              IF (IELT.EQ.0) 
     *          CALL TADSTP(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,DWORK(KT),DWORK(KT),
     *          DWORK(L(LDT)),DWORK(L(LDT)),VWORK(KA1),NA,KWORK(KCOLA),
     *          KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *          DWORK(L(LCORVG)),E031,COEFNN,0,1,1D0)
              IF (IELT.EQ.1) 
     *          CALL TADSTP(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,DWORK(KT),DWORK(KT),
     *          DWORK(L(LDT)),DWORK(L(LDT)),VWORK(KA1),NA,KWORK(KCOLA),
     *          KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *          DWORK(L(LCORVG)),E030,COEFNN,0,1,1D0)
              IF (IELT.EQ.2) 
     *          CALL TADSTN(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,DWORK(KT),DWORK(KT),
     *          DWORK(L(LDT)),DWORK(L(LDT)),VWORK(KA1),NA,KWORK(KCOLA),
     *          KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *          DWORK(L(LCORVG)),EM31,COEFNN,0,1,1D0)
              IF (IELT.EQ.3) 
     *          CALL TADSTN(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,DWORK(KT),DWORK(KT),
     *          DWORK(L(LDT)),DWORK(L(LDT)),VWORK(KA1),NA,KWORK(KCOLA),
     *          KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *          DWORK(L(LCORVG)),EM30,COEFNN,0,1,1D0)
            endif
            if (ISTOK.NE.1)
     *          CALL TGUPWD(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,
     *              DWORK(KT),DWORK(KT),
     *              DWORK(L(LDT)),DWORK(L(LDT)),
     *              VWORK(KA1),NA,KWORK(KCOLA),
     *              KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *              DWORK(L(LCORVG)),0,1)
         ELSE
          IF (IELT.EQ.0) 
     *     CALL TUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                 1D0,0D0,DWORK(KT),DWORK(KT),
     *                 DWORK(L(LDT)),DWORK(L(LDT)),
     *                 VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                 KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                 E031,COEFFN,0,1,1D0)
          IF (IELT.EQ.1) 
     *      CALL TUPWDG(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                 1D0,0D0,DWORK(KT),DWORK(KT),
     *                 DWORK(L(LDT)),DWORK(L(LDT)),
     *                 VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                 KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                 E030,COEFFN,0,1,1D0)
          IF (IELT.EQ.2) 
     *     CALL TUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                 1D0,0D0,DWORK(KT),DWORK(KT),
     *                 DWORK(L(LDT)),DWORK(L(LDT)),
     *                 VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                 KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                 EM31,COEFFN,0,1,1D0)
          IF (IELT.EQ.3) 
     *     CALL TUPWNP(DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                 1D0,0D0,DWORK(KT),DWORK(KT),
     *                 DWORK(L(LDT)),DWORK(L(LDT)),
     *                 VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                 KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                 EM30,COEFFN,0,1,1D0)
         ENDIF
        ENDIF
C
        CALL  BDRYAT (VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *               KWORK(L(KLMBD(ILEV))),KNMBD(ILEV),
     *               KWORK(L(LNPR)),DWORK(L(KLDBD(ILEV))))
C
C
22      CONTINUE
      ENDIF
C
C=======================================================================
c	Loesung
C=======================================================================
c
C
       ISETLV=2
       ILEV=NLMAX
       CALL SETLEV(ISETLV)

      CALL  LCP1 (DWORK(KT),DWORK(L(LTOLD)),KNU(nlmax))
      CALL ZTIME(TTDEF1)
      TTILIN=TTILIN+(TTDEF1-TTDEF0)
C
C=======================================================================
C	Multigrid iteration!
C	NO MATRIX RESORTING, at least at the moment
C=======================================================================
c
      DO 12  ILEV=NLMIN,NLMAX
      KOFFX(ILEV)=L(KLUP  (ILEV))-1+2*KNU(ILEV)+KNP(ILEV)
      KOFFB(ILEV)=L(KLF12P(ILEV))-1+2*KNU(ILEV)+KNP(ILEV)
      KOFFD(ILEV)=L(KLAUX (ILEV))-1+2*KNU(ILEV)+KNP(ILEV)
      KNEQ (ILEV)=KNU(ILEV)
      KPRSM(ILEV)=NSMT*NSMUFA**(NLMAX-ILEV)
      KPOSM(ILEV)=NSMT*NSMUFA**(NLMAX-ILEV)
12    CONTINUE
C
      ICYCLE=ICYCU
      IRELMG=1
      EPSUMG=1D99
      IF (ILMAXU.GT.ILMINU) THEN
       IDEFMG=1
      ELSE
       IDEFMG=0
      ENDIF
      CALL  M011 (DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
     *            KNEQ,ILMAXT,ITMG,DMPTMP,EPSUMG,DEFUMG,
     *            YAXU,YPROLU,YRESTU,YSMU,YSMU,YEXU,YEXAU,
     *            YDBCU,YSTEPU,
     *            KIT0,KIT,IRELMG,IDEFMG,RHOMGT,BMGT)
c
      ISETLV=2
      ILEV=NLMAX
      CALL SETLEV(ISETLV)
c	do 27 i=1,20
c27	write (*,*) i,dwork(ku1-1+i),dwork(ku2-1+i)
C
      CALL MGCORT (NU)
      CALL ZTIME(TTDEF2)
      TIMTMG=TIMTMG+(TTDEF2-TTDEF1)
C
C=======================================================================
C *** Calculation of defects
C=======================================================================
C
       CALL LCP1 (DWORK(KFT),DWORK(L(LDT)),NU)
C
      FAKTOR=ALFA/NY
      THSTEP=THSTEP*FAKTOR
      CALL XMDF2T(KM1,KMASS1,KST1,KA1,KCOLA,KLDA,L(LDT),KT,
     *            NA,NU,KWORK(L(KLMBD(NLEV))),KNMBD(NLEV),THSTEP)
      THSTEP=THSTEP/FAKTOR
C
      IF ((ISTOK.NE.1).OR.((ISTOK.EQ.1).AND.(IPRECA.EQ.4))) THEN
C
C
        IF (IUPW.EQ.1) THEN
            if (IPRECA.eq.4) then
              IF (IELT.EQ.0) 
     *          CALL TADSTP(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,DWORK(KT),DWORK(KT),
     *          DWORK(L(LDT)),DWORK(L(LDT)),VWORK(KA1),NA,KWORK(KCOLA),
     *          KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *          DWORK(L(LCORVG)),E031,COEFNN,IDEFUP,1,1D0)
              IF (IELT.EQ.1) 
     *          CALL TADSTP(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,DWORK(KT),DWORK(KT),
     *          DWORK(L(LDT)),DWORK(L(LDT)),VWORK(KA1),NA,KWORK(KCOLA),
     *          KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *          DWORK(L(LCORVG)),E030,COEFNN,IDEFUP,1,1D0)
              IF (IELT.EQ.2) 
     *          CALL TADSTN(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,DWORK(KT),DWORK(KT),
     *          DWORK(L(LDT)),DWORK(L(LDT)),VWORK(KA1),NA,KWORK(KCOLA),
     *          KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *          DWORK(L(LCORVG)),EM31,COEFNN,IDEFUP,1,1D0)
              IF (IELT.EQ.3) 
     *          CALL TADSTN(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,DWORK(KT),DWORK(KT),
     *          DWORK(L(LDT)),DWORK(L(LDT)),VWORK(KA1),NA,KWORK(KCOLA),
     *          KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *          DWORK(L(LCORVG)),EM30,COEFNN,IDEFUP,1,1D0)
            endif
            if (ISTOK.NE.1)
     *          CALL TGUPWD(DWORK(KU1),DWORK(KU2),
     *              DWORK(KU1),DWORK(KU2),1d0,0d0,
     *              DWORK(KT),DWORK(KT),
     *              DWORK(L(LDT)),DWORK(L(LDT)),
     *              VWORK(KA1),NA,KWORK(KCOLA),
     *              KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *              DWORK(L(LCORVG)),IDEFUP,1)
        ELSE
         IF (IELT.EQ.0) 
     *    CALL TUPWDG(DWORK(KU1),DWORK(KU2),
     *                DWORK(KU1),DWORK(KU2),1d0,0d0,
     *                DWORK(KT),DWORK(KT),
     *                DWORK(L(LDT)),DWORK(L(LDT)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                E031,COEFFN,IDEFUP,1,1D0)
         IF (IELT.EQ.1) 
     *    CALL TUPWDG(DWORK(KU1),DWORK(KU2),
     *                DWORK(KU1),DWORK(KU2),1d0,0d0,
     *                DWORK(KT),DWORK(KT),
     *                DWORK(L(LDT)),DWORK(L(LDT)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                E030,COEFFN,IDEFUP,1,1D0)
         IF (IELT.EQ.2) 
     *    CALL TUPWNP(DWORK(KU1),DWORK(KU2),
     *                DWORK(KU1),DWORK(KU2),1d0,0d0,
     *                DWORK(KT),DWORK(KT),
     *                DWORK(L(LDT)),DWORK(L(LDT)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                EM31,COEFFN,IDEFUP,1,1D0)
         IF (IELT.EQ.3) 
     *    CALL TUPWNP(DWORK(KU1),DWORK(KU2),
     *                DWORK(KU1),DWORK(KU2),1d0,0d0,
     *                DWORK(KT),DWORK(KT),
     *                DWORK(L(LDT)),DWORK(L(LDT)),
     *                VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *                KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *                EM30,COEFFN,IDEFUP,1,1D0)
        ENDIF
      ENDIF
C
      CALL BDRYAT (VWORK(KA1),NA,KWORK(KCOLA),KWORK(KLDA),
     *            KWORK(L(KLMBD(ILEV))),KNMBD(ILEV),
     *            KWORK(L(LNPR)),DWORK(L(KLDBD(ILEV))))
C
      CALL BDRY0T (DWORK(L(LDT)),KWORK(L(KLMBD(ILEV))),
     *             KNMBD(ILEV),
     *            KWORK(L(LNPR)),DWORK(L(KLDBD(ILEV))))
C
      CALL  RESDFK(DWORK(L(LDT)),DWORK(L(LDT)),DWORK(KFT),DWORK(KFT),
     *              NU,REST,RDUMMY)
C
C
999      CALL ZTIME(TTDEF3)
      TIMTPO=TIMTPO+(TTDEF3-TTDEF2)
c
C
      WRITE(MTERM,1001)
      WRITE(MFILE,1001)
      WRITE(MTERM,1)
      WRITE(MFILE,1)
      WRITE(MTERM,1003) ITMG,REST,RHOMGT
      WRITE(MFILE,1003) ITMG,REST,RHOMGT
      WRITE(MTERM,1)
      WRITE(MFILE,1)
C
      GOTO 99999
c
   1  FORMAT(80('-'))
   3  FORMAT(80('+'))
1000  FORMAT (6E12.5)
1001  FORMAT('  IT DEF-T',4X,'RHOMGT')
1003  FORMAT(I4,3(D9.2))
10002 FORMAT ('#',I4,'(',I4,')',1X,'TIME=',D10.3,1X,'REL2(P)=',
     *        D10.3,1X,'RELM(P)=',D10.3)
c
99999 END
