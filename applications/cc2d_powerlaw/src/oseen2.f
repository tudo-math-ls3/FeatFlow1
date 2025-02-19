************************************************************************
      SUBROUTINE XOSEE (MFILE,MSHOW,BMG,RHOLMG,INL)  
************************************************************************
C
C-----------------------------------------------------------------------
C PURPOSE:   SOLVER FOR THE GENERALIZED INCOMPRESSIBLE OSEEN EQUATONS
C
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'bouss.inc'
      INCLUDE 'block.inc'
C-----------------------------------------------------------------------
C
      CHARACTER CFILE*15
C
      COMMON /MGFLDS/ KLXS(NNLEV),KLFS(NNLEV),KLAUXS(NNLEV)
C
C *** ARRAYS FOR MULTIGRID MODUSL M011'S
c      DIMENSION  KOFFX(NNLEV),KOFFD(NNLEV),KOFFB(NNLEV),
c     *           KNEQ(NNLEV),KIT(NNLEV),KIT0(NNLEV)
      DIMENSION  KOFFXS(NNLEV),KOFFDS(NNLEV),KOFFBS(NNLEV),KOFFRS(NNLEV)
      DIMENSION  KNEQS(NNLEV),KLM1H(NNLEV),KLRS(NNLEV)
C
C-----------------------------------------------------------------------
C     E X T E R N A L S
C-----------------------------------------------------------------------
C *** COEFFICIENT OF STIFFNESS MATRIX
      EXTERNAL COEFFN
C *** DEFINITION OF FINITE ELEMENTS
      EXTERNAL E030,E031,EM30,EM31
C *** MULTIGRID COMPONENTS
      EXTERNAL  YAXA,YSMA,YEXvk
c      EXTERNAL  YAXS,YPROLS,YRESTS,YSMS,YEXS,YEXAS,YDBCS,YSTEPS
      EXTERNAL  YAX,YPROL2,YREST2,YPROL,YREST,YSM,YEX,YEXA,YDBC,YSTEP
c$$$      EXTERNAL  YAX,YPROL,YREST,YSM,YEX,YEXA,YDBC,YSTEP
c
C=======================================================================
C     INITIALIZATION OF DIRECT SOLVER
C=======================================================================
C
C
       DO 2001 ILEV=NLMIN,NLMAX
       CALL ZTIME(TTT0)
       ISETLV=2
       CALL SETLEV(ISETLV)
C



       CALL ZNEW (NEL+1,3,KLLDPA(ILEV),'KLDPA ')
       CALL ZNEW (NEL  ,3,KLPAT (ILEV),'KPAT  ')
       CALL ZNEW (NEL  ,3,KLPAEL(ILEV),'KPAEL ')
       CALL ZNEW (NEL,  3,LCHECK,'KCHECK')
       IF (IER.NE.0) GOTO 99999
C
       LLDPAT=KLLDPA(ILEV)
       LPATCH=KLPAT (ILEV)
       LPATEL=KLPAEL(ILEV)
C
       LELC1 =KLELC1(ILEV)
       LELC2 =KLELC2(ILEV)
C
       IF (ILEV.NE.NLMIN) THEN
        NBLOCK=NSBLSM
       ELSE
        NBLOCK=NSBLSL
       ENDIF
C
c       nblock=min(nblock,10)
c      goto 998

C=========================patching======================================
c       IF ((ILEV.EQ.NLMIN).AND.(BLOCSL)) THEN
       IF (((ILEV.EQ.NLMIN).AND.(BLOCSL)).OR.(PSBARB.eq.1d0)) THEN
        NPATCH=1
        KWORK(L(LLDPAT))=1
        KWORK(L(LLDPAT)+1)=NEL+1
        DO 2002 I=1,NEL
2002    KWORK(L(LPATCH)-1+I)=I
       else  
C=======================================================================
c
          IF (IBLOCK.EQ.3) THEN
          CALL SGRP3(NPATCH,NBLOCK,KWORK(L(LPATCH)),KWORK(L(LLDPAT)),
     *               KWORK(L(LCHECK)),KWORK(L(LELC1)),KWORK(L(LELC2)),
     *               KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *               DWORK(L(LCORVG)))
          ELSE
c

          CALL SGRP1(NPATCH,NBLOCK,KWORK(L(LPATCH)),KWORK(L(LLDPAT)),
     *               KWORK(L(LCHECK)),KWORK(L(LELC1)),KWORK(L(LELC2)),
     *               KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *               DWORK(L(LCORVG)),DWORK(l(KLASPR(ILEV))))
          ENDIF
        endif
c
       KNPTCH(ILEV)=NPATCH
       CALL ZDISP (0,LCHECK,'KCHECK')
       CALL ZDISP (NPATCH+1,LLDPAT,'LLDPAT')
       IF (IER.NE.0) GOTO 99999
C
c       
       ILDMAX=0
       DO 2010 IPATCH=1,NPATCH
       ILDMAX=MAX(ILDMAX,KWORK(L(LLDPAT)+IPATCH)-
     *                   KWORK(L(LLDPAT)+IPATCH-1))
       DO 2011 ILDPAT=KWORK(L(LLDPAT)+IPATCH-1),
     *                KWORK(L(LLDPAT)+IPATCH)-1
cc       WRITE (6,*) ILEV,IPATCH,ILDPAT,KWORK(L(LPATCH)+ILDPAT-1)
       IEL=KWORK(L(LPATCH)+ILDPAT-1)
       KWORK(L(LPATEL)+IEL-1)=IPATCH
2011   CONTINUE
c
2010   CONTINUE
C


C
       CFILE='#gmv/patch    '
       WRITE(CFILE(11:11),'(I1.1)') ILEV
c

       CALL XGMVPA(70,CFILE,KWORK(L(LVERT)),DWORK(L(LCORVG)),
     *             KWORK(L(LPATEL)))
c

c$$$       CALL XAVSPA(70,CFILE,KWORK(L(LVERT)),DWORK(L(LCORVG)),
c$$$     *             KWORK(L(LPATEL)))
C
       IF (INL.eq.1) THEN
       IF (MSHOW.GE.2) 
     *  WRITE(MTERM,*) 'ILEV,NPATCH,NEL,MAX-PATCH:',ILEV,KNPTCH(ILEV),
     *                  KNEL(ILEV),ILDMAX
       IF (MSHOW.GE.0) 
     *  WRITE(MFILE,*) 'ILEV,NPATCH,NEL,MAX-PATCH:',ILEV,KNPTCH(ILEV),
     *                  KNEL(ILEV),ILDMAX
       ENDIF
cC
cC
2001   CONTINUE
C=======================================================================
c       npatch=1
C=======================================================================

c      ENDIF
C
C=======================================================================
C     INITIALIZATION OF COUPLED MG
C=======================================================================
C
C
      CALL ZTIME(TTT0)
C
C *** INITIALIZATION OF THE OFFSET ARRAYS KOFFX,KOFFB,KOFFD AND KNEQ
c       DO 10 ILEV=NLMIN,NLMAX
c      KOFFX(ILEV)=L(KLUP(ILEV))-1
c      KOFFB(ILEV)=L(KLF12P(ILEV))-1
c      KOFFD(ILEV)=L(KLAUX(ILEV))-1
c      KNEQ (ILEV)=KNUP(ILEV)
c      KPRSM(ILEV)=NPRSMC
c      KPOSM(ILEV)=NPOSMC
c10    CONTINUE
C
c      DO 11  ILEV=NLMIN,NLMAX
c      KPOSM(ILEV)=KPOSM(ILEV)*NCSMC**(NLMAX-ILEV)
c      KPRSM(ILEV)=KPRSM(ILEV)*NCSMC**(NLMAX-ILEV)
c11    CONTINUE
C *** Initialization of the offset arrays KOFFX,KOFFB,KOFFD and KNEQ

      DO 1017 ILV=NLMIN,NLMAX
      KOFFX(ILV)=L(KLUP(ILV))-1
      KOFFB(ILV)=L(KLF12P(ILV))-1
      KOFFD(ILV)=L(KLAUX(ILV))-1       
      KNEQ(ILV)=KNUP(ILV)
      KPOSM(ILV)=NSM*NSMFAC**(NLEV-ILV)
      KPRSM(ILV)=NSM*NSMFAC**(NLEV-ILV)
 1017 CONTINUE
C
      IRELMG=1
      ICYCLE=ICYCle
c      ICYCLE=ICYCC
      ITMG=ILMIN
      IF (ILMAX.GT.ILMIN) THEN
       IDEFMG=1
      ELSE
       IDEFMG=0
      ENDIF

        IF  (MODINT.ge.0) THEN
c           write (*,*) 'klassische Rest/Prola' 
        CALL  M011 (DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
     *               KNEQ,ILMAX,ITMG,DMPMG,EPSMG,
     *               YAX,YPROL,YREST,YSMa,YSMa,YEXvk,YEXA,YDBC,YSTEP,
     *               KIT0,KIT,IRELMG,IDEFMG,RHOLMG,BMG)
        ELSE !   (MODINT.ge.0)
c           write (*,*) 'modifizierte Rest/Prolb'
c
         CALL  M011 (DWORK(1),DWORK(1),DWORK(1),KOFFX,KOFFB,KOFFD,
     *               KNEQ,ILMAX,ITMG,DMPMG,EPSMG,
     *               YAX,YPROL2,YREST2,YSMa,YSMa,YEXvk,YEXA,YDBC,YSTEP,
     *               KIT0,KIT,IRELMG,IDEFMG,RHOLMG,BMG)
         ENDIF!   (MODINT.ge.0)
c
         NMG=NMG+ITMG
c
      TTMGC=TTMGC+TTMG
      TTSC=TTSC+TTS
      TTEC=TTEC+TTE
      TTDC=TTDC+TTD
      TTPC=TTPC+TTP
      TTRC=TTRC+TTR
C
      IF (IER.GT.0) IER=0
      IF (IER.NE.0) GOTO 99999
C
C=======================================================================
C
       DO 3001 ILEV=NLMIN,NLMAX
       CALL ZTIME(TTT0)
       ISETLV=2
       CALL SETLEV(ISETLV)
C
       KNPTCH(ILEV)=0
C      ENDIF

       CALL ZDISP(0,KLLDPA(ILEV),'KLDPA ')
       CALL ZDISP(0,KLPAT (ILEV),'KPAT  ')
       CALL ZDISP(0,KLPAEL(ILEV),'KPAEL ')
       IF (IER.NE.0) GOTO 99999
C
3001   CONTINUE
C
C

C=======================================================================
      CALL  ZTIME(TTT0)
C
      ISETLV=2
      ILEV=NLMAX
      CALL  SETLEV (ISETLV)
C
c      CALL LCP1(DWORK(L(LDEF)),DWORK(KF1),NUP)
c      CALL ZDISP(0,LDEF,'DDEF  ')
c      IF (IER.NE.0) GOTO 99999
C
C       DO 154 II=1,KNU(NLMAX)
C       DO 154 II=1,10
C154    WRITE(6,*) II,DWORK(KU1+II-1),DWORK(KU2+II-1)
C
C       DO 155 II=1,KNP(NLMAX)
C       DO 155 II=1,10
C155    WRITE(6,*) II,DWORK(KP+II-1)
C
       CALL ZTIME(TTT1)
       TTLC=TTLC+TTT1-TTT0
C
C

C
99999 END
C
c========================================================
c     Block Vanca
c========================================================
************************************************************************
      SUBROUTINE  YSMA (DX,DB,DD,NEQ,NSMO)  
************************************************************************
C
C-----------------------------------------------------------------------
*   Purpose: - performs NSMO smoothing steps applied to the system
*                          A*DX = DB
*              of dimension NEQ using the auxiliary vector DD
*            - DX,DB,DD have the structure  D=(D1,D2,DP)
C
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'bouss.inc'
      INCLUDE 'block.inc'
C-----------------------------------------------------------------------
C
      DIMENSION DX(*),DB(*),DD(*)
C
C
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
      LLDPAT=KLLDPA(ILEV)
      LPATCH=KLPAT (ILEV)
      LPATEL=KLPAEL(ILEV) 
      NPATCH=KNPTCH(ILEV)
C
C=======================================================================

ccc$       CALL ZNEW(NUP,-1,LRHS1,'DRHS1 ')
       IF (IER.NE.0) GOTO 99999
C
ccc$       CALL LCP1 (DB,DWORK(L(LRHS1)),NUP)
C
       CALL ZNEW(NUP,-1,LRHS1,'DRHS1 ')
       IF (IER.NE.0) GOTO 99999
C
       CALL LCP1 (DB,DWORK(L(LRHS1)),NUp)
c
c
       DO 1010 ITE=1,NSMO
c
       CALL LCP1(DWORK(L(LRHS1)),DB,NUP)
c
       STOP
       CALL YAX(DX,DB,NUP,-1D0,1D0)

c$$$C
       CALL LCP1(DX,DD,NUP)
       CALL LCL1(DX,NUP)
C


       CALL PATCH (DX(1),DX(I2),DX(IP),DB(1),DB(I2),DB(IP),
     *             KWORK(L(LVERT)),KWORK(L(LMID)),LLDPAT,LPATCH,
     *             LPATEL,NPATCH,ILOES)



c
       CALL  LLC1 (DD,DX,NUP,1D0,RLXSM)


       IF (INEUM.EQ.0) CALL TOL20A(DX(IP),VWORK(KAREA),NP,INEUM)
1010   CONTINUE
c
       CALL LCP1(DWORK(L(LRHS1)),DB,NUP)
C
       CALL ZDISP(0,LRHS1,'DRHS1 ')
       IF (IER.NE.0) GOTO 99999
C

 99    continue
c
       IF (IER.NE.0) GOTO 99999
C
C
C
99999 END
C
C
************************************************************************
      SUBROUTINE PATCH (U1,U2,P,F1,F2,FP,KVERT,KMID,LLDPAT,LPATCH,
     *                  LPATEL,NPATCH,ILOES)
************************************************************************
C
C EINGABE: NPATCH: 	     ZAHL DER PATCHES
C	   KPATCH(NEL):	     ELEMENTNUMMER, GEMAESS PATCHES SORTIERT
C	   KLDPAT(NPATCH+1): STARTADRESSE DES PATCHES AUF KPATCH
C                            (LAST ENTRY:NEL+1)
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'bouss.inc'
      INCLUDE 'block.inc'
      CHARACTER CFILE*6
c
      SAVE
C-----------------------------------------------------------------------
C
c      PARAMETER (NNVE=4)

C
      DIMENSION U1(*),U2(*),P(*)
      DIMENSION F1(*),F2(*),FP(*)
      DIMENSION KMID(NNVE,*),KVERT(NNVE,*)
      EXTERNAL EM31,EM30
C
c      mt=IZBV1
      igold=igrad
c      igrad=0!!! Fuer TEstzwecke
      IAUA=10
      idi1=0
      ifm1=0
c
      KMBD =L(KLMBD(ILEV))  ! this was not expl defined here ??? jaro

      CALL EA00(EM31,EM31,NVE,IELTYP)
C
      DO 111 IPATCH =1,NPATCH
c
C        NELPAT entspricht der Zahl der Elemente im Patch
         NELPAT=(KWORK(L(LLDPAT)-1+IPATCH+1)-
     *            KWORK(L(LLDPAT)-1+IPATCH))
         ILANG=NELPAT*4
c
c===================================================
         IF (NELPAT.EQ.1) THEN
c
            II=KWORK(L(LLDPAT)-1+IPATCH)
            IELE=KWORK(L(LPATCH)+II-1)
c
            IF (ISM.eq.1) ISML=1! klassiche Diagonale
            IF (ISM.eq.2) ISML=2! Volle 8x8 Matrix
            IF (ISM.eq.3) ISML=3! 2 4x4 Matrizen
            IF (ISM.eq.4) THEN
               CALL ARCALC (ARIEL,IELE,KVERT,DWORK(L(LCORVG))) 
               IF (ARIEL.LE.4d0) THEN
                  ISML=1
               ELSE
                  ISML=2
               ENDIF
            ENDIF
c
c            IF ((IGRAD.EQ.1).AND.(ISML.EQ.2)) ISML=2
c           (In der Gradientenform bringt ISML=2 keinen Vorteil)
c
            IF (ISML.eq.1) THEN   
               idi1=idi1+1
            CALL VANCAM (U1,U2,P,F1,F2,FP,
     *           VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *           VWORK(KB1),VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *           KWORK(KMBD),KVERT,KMID,
     *           KWORK(L(LNPR)),NMBD,0,IELE)
c
            ENDIF !isml=1
c
            IF (ISML.eq.2) THEN   
c               IF (iele.eq.1) write (*,*) 'ich mache 8x8', ilev,ite
               ifm1=ifm1+1
            CALL VANCR2 (U1,U2,P,F1,F2,FP,
     *           VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *           VWORK(KB1),VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *           KWORK(KMBD),KVERT,KMID,DWORK(L(LCORVG)),
     *           KWORK(L(LNPR)),NMBD,0,IELE,ite)
c
            ENDIF
c
            IF (ISML.eq.3) THEN   
c               IF (iele.eq.1) write (*,*) 'ich mache 4x4'      IZBV3=0

c         
            CALL VANCAR (U1,U2,P,F1,F2,FP,
     *           VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *           VWORK(KB1),VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *           KWORK(KMBD),KVERT,KMID,
     *           KWORK(L(LNPR)),NMBD,0,IELE)
            ENDIF
c
         ELSE!NELPAT.NE.1
c
            CALL ZTIME (TPAT0)
         CALL ZNEW(ILANG,3,LFHG,'FHG   ')
C         
c------------------------------------------------------
C        		IN DIESEN VEKTOR KOMMEN DIE GLOBALEN 
C			FREIHEITSGRADE DES PATCHES
C        		LAENGE IST 4*ANZAHL ELEMENTE IM PATCH.
C        		ICH HABE FUER DAS PATCH KEIN NMT!
c------------------------------------------------------
C
         CALL ZNEW (NMT,3,LFLAG,'LFLAG ')
C
         LAUF=0
         DO 2221 II=KWORK(L(LLDPAT)-1+IPATCH),
     *                   (KWORK(L(LLDPAT)-1+IPATCH+1)-1)
C
            IELE=KWORK(L(LPATCH)+II-1)
            DO 3331 JJ=1,4
               IH=KMID(JJ,IELE)-NVT
               IF (KWORK(L(LFLAG)-1+IH).EQ.0) THEN
                  KWORK(L(LFLAG)-1+IH)=1
                  KWORK(L(LFHG)+LAUF)=IH
                  LAUF=LAUF+1
               ENDIF
3331        CONTINUE           
2221     CONTINUE
C
         NEQPAT=LAUF
         CALL ZDISP(0,LFLAG,'LFLAG ')
         CALL ZDISP(NEQPAT,LFHG,'FHG   ')
C
c------------------------------------------------------
C        Aufbau der lokalen Matrix A in einem  Vektor der Laenge
C        NEQPAT**2 (bei Grad-Grad) oder 4NEQPAT**2 (bei Def-Def)
c------------------------------------------------------
C
         NAT=NEQPAT*NEQPAT ! NA auf dem Patch
C
         IF (IGRAD.EQ.1) THEN
            NALOC=  NAT
            NVELO=  NEQPAT
            NPRES=  NELPAT
         ELSE
            NALOC=4*NAT
            NVELO=2*NEQPAT
            NPRES=  NELPAT
         ENDIF
c
         CALL ZNEW(NALOC,1,LMATRI,'MATRIX')
C
C
c------------------------------------------------------
         IF ((.NOT.BMATCL).AND.(ILEV.EQ.NLMIN)) GOTO 242
c------------------------------------------------------
c        Routine zur 'COLlection of MATrix elements'
c------------------------------------------------------
         CALL COLMAT(DWORK(L(LMATRI)),LFHG,NAT,NVELO,NEQPAT,NELPAT)
c-------------------------------------------
C
 242    CONTINUE !((.NOT.BMATCL).AND.(ILEV.EQ.NLMIN))

c
C
        CALL ZNEW (NVELO,3,LPIV ,'LPIV  ')
C
c-----------------------------------------------------------------
         IF ((.NOT.BMATCL).AND.(ILEV.EQ.NLMIN)) GOTO 244
         CALL ZNEW(NEQPAt*NEQPAT,1,LINVMA,'LINVMA')
c
       LUMFM1=0
       LUMFI1=0
c     
       NV1=IZBV3*NVELO*NVELO

        CALL ZNEW (NV1,1,LUMFM1,'LUMFM1')! Vektor, den UMFACK fuer die Matrix braucht
        CALL ZNEW (2*NV1,3,LUMFI1,'LUMIN1') ! Hilfsvektor,fuer die Indizes
c

c
        CALL ZTIME (TFACA0)
       IF((INT(OME1).GE.1).and.(ipatch.eq.1).and.(ilev.eq.nlmax))
     *       write (*,*) tfaca0-tpat0, ' Patch bis zu Lu von A'
        CALL FACTOR (NVELO,DWORK(L(LMATRI)),KWORK(L(LPIV)),
     *              KEEP1,CNTL1,ICNTL1,LUMFM1,LUMFI1)
        CALL ZTIME (TFACA1)
       IF((INT(OME1).GE.1).and.(ipatch.eq.1).and.(ilev.eq.nlmax))
     *   print *, TFACA1-TFACA0, ' LU von A' 

c-----------------------------------------------------------------
C   So, jetzt ist die Matrix A1 (LOKAL) AUFGEBAUT,
C   und invertiert. Es geht an die B-Matrizen
c-----------------------------------------------------------------
 244    CONTINUE
C
         CALL ZNEW(NEQPAT*NELPAT,1,LB11,'LB11  ')
         CALL ZNEW(NEQPAT*NELPAT,1,LB22,'LB22  ')
c
         CALL COLMAB(DWORK(L(LB11)),DWORK(L(LB22)),
     *               VWORK(KB1),    VWORK(KB2),            
     *               NELPAT,NEQPAT,LFHG,LPATCH,LLDPAT,IPATCH,KLDB,KCOLB)            
C
         CALL ZNEW (NELPAT*NELPAT,1,LSK,'LSK   ')
         CALL ZNEW (NELPAT,3,LPIV1,'LPIV1 ')
C
Cc----------------------------------
         IF ((.NOT.BMATCL).AND.(ILEV.EQ.NLMIN)) GOTO 246
c----------------------------------
c

        IF (IFACTO.NE.0) THEN
        IF (BSKBLD(ILEV)) THEN!Aufbau der LU-Zerlegung von P, sowie des Pivotvektors
           Call ZTIME (TT31)
           CALL SKBLD(DWORK(L(LB11)),DWORK(L(LB22)),
     *             DWORK(L(LMATRI)),KWORK(L(LPIV )),
     *             DWORK(L(LSK)),NELPAT,NEQPAT,NVELO,
     *             IGRAD,IPATCH,NPATCH)
           Call ZTIME (TT32)
           IF ((INT(OME1).GE.1).and.(ipatch.eq.1).and.(ilev.eq.nlmax))
     *        write (*,*) tt32-tt31 ,' Aufbauen der SK-Matrix:',ilev 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C   Auf LSK steht jetzt die Schurkomplement Matrix
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 

c---------------------------------------------------------
           IF ((.NOT.BMATCL).AND.(ILEV.EQ.NLMIN)) GOTO 248
c---------------------------------------------------------
           LUMFM2=0
           LUMFI2=0
           N1=NELPAT
           NM2=N1*N1*IZBV3
c
           IFC1=IFACTO
           IFACTO=0
           CALL ZTIME (TFACB0)
           CALL FACTOR (NELPAT,DWORK(L(LSK)),KWORK(L(LPIV1)),
     *               KEEP2,CNTL2,ICNTL2, LUMFM2,LUMFI2)
           CALL ZTIME (TFACB1)
           IFACTO=IFC1
           IF ((INT(OME1).GE.1).and.(ipatch.eq.1).and.(ilev.eq.nlmax))
     *        print *, TFACB1-TFACB0, ' LU von P' ,ipatch
c
           IF (BLSKSA) THEN
c             ===================================================
c             Jetzt hamma die Matrix, nu speichern wir sie ab!
c             ===================================================
c
             CFILNM='#mat/ALSK.          '
             WRITE(CFILNM(11:11),'(I1.1)') ILEV
             WRITE(CFILNM(12:12),'(A1)') '.'
c
             IF ((IPATCH.GE.0).AND.(IPATCH.LT.10)) 
     *       WRITE(CFILNM(13:13),'(I1.1)') ipatch
             IF ((IPATCH.GE.10).AND.(IPATCH.LT.100)) 
     *       WRITE(CFILNM(13:14),'(I2.2)') ipatch
             IF ((IPATCH.GE.100).AND.(IPATCH.LT.1000)) 
     *       WRITE(CFILNM(13:15),'(I3.3)') ipatch
             IF ((IPATCH.GE.1000).AND.(IPATCH.LT.10000)) 
     *       WRITE(CFILNM(13:16),'(I4.4)') ipatch
c 
             CALL  OF0 (59,CFILNM,0)!oeffnet file
             CFILE='SKMAT '
             CALL  OWA1 (DWORK(L(LSK)),CFILE,nelpat*nelpat,59,0)
             REWIND(59)
             CLOSE (59)
c
             CFILNM='#mat/PIVO.          '
             WRITE(CFILNM(11:11),'(I1.1)') ILEV
             WRITE(CFILNM(12:12),'(A1)') '.'
             IF ((IPATCH.GE.0).AND.(IPATCH.LT.10)) 
     *        WRITE(CFILNM(13:13),'(I1.1)') ipatch
             IF ((IPATCH.GE.10).AND.(IPATCH.LT.100)) 
     *        WRITE(CFILNM(13:14),'(I1.1)') ipatch
             IF ((IPATCH.GE.100).AND.(IPATCH.LT.1000)) 
     *        WRITE(CFILNM(13:15),'(I1.1)') ipatch
             IF ((IPATCH.GE.1000).AND.(IPATCH.LT.10000)) 
     *        WRITE(CFILNM(13:16),'(I1.1)') ipatch

             CALL  OF0 (59,CFILNM,1)
             CFILE='PIVOT '
c             write (*,*) (kwork(l(lpiv1)-1+i),i=1,5)
c             CALL  OWA1 (DWORK(L(LSK)),CFILE,nelpat*nelpat,59,0)
             CALL  OWA3 (KWORK(L(LPIV1)),CFILE,NELPAT,59,1)
             REWIND(59)
             CLOSE (59)
             IF (IPATCH.EQ.NPATCH) BSKBLD(ilev)=.false.! damit wir das nicht jedes Mal wieder machen
           ENDIF!BLSKSA
c             ===================================================
        ELSE!Einlesen der LU-Zerlegung von P, sowie des Pivotvektors
           Call ZTIME (TT33)
           CFILNM='#mat/ALSK.          '
           WRITE(CFILNM(11:11),'(I1.1)') ILEV
           WRITE(CFILNM(12:12),'(A1)') '.'
           IF ((IPATCH.GE.0).AND.(IPATCH.LT.10)) 
     *        WRITE(CFILNM(13:13),'(I1.1)') ipatch
           IF ((IPATCH.GE.10).AND.(IPATCH.LT.100)) 
     *        WRITE(CFILNM(13:14),'(I1.1)') ipatch
           IF ((IPATCH.GE.100).AND.(IPATCH.LT.1000)) 
     *        WRITE(CFILNM(13:15),'(I1.1)') ipatch
           IF ((IPATCH.GE.1000).AND.(IPATCH.LT.10000)) 
     *        WRITE(CFILNM(13:16),'(I1.1)') ipatch

           CALL  OF0 (59,CFILNM,0)
           CFILE='SKMAT '
           CALL  ORA1 (DWORK(L(LSK)),CFILE,59,0)
           REWIND(59)
           CLOSE (59)
c
           CFILNM='#mat/PIVO.          '
           WRITE(CFILNM(11:11),'(I1.1)') ILEV
           WRITE(CFILNM(12:12),'(A1)') '.'
           IF ((IPATCH.GE.0).AND.(IPATCH.LT.10)) 
     *        WRITE(CFILNM(13:13),'(I1.1)') ipatch
           IF ((IPATCH.GE.10).AND.(IPATCH.LT.100)) 
     *        WRITE(CFILNM(13:14),'(I1.1)') ipatch
           IF ((IPATCH.GE.100).AND.(IPATCH.LT.1000)) 
     *        WRITE(CFILNM(13:15),'(I1.1)') ipatch
           IF ((IPATCH.GE.1000).AND.(IPATCH.LT.10000)) 
     *        WRITE(CFILNM(13:16),'(I1.1)') ipatch

           CALL  OF0 (59,CFILNM,1)
           CFILE='PIVOT '
           CALL  ORA3 (KWORK(L(LPIV1)),CFILE,59,1)
           REWIND(59)
           CLOSE (59)
c             write (*,*) (kwork(l(lpiv1)-1+i),i=1,5)
c

           Call ZTIME (TT34)
           Call ZTIME (TT32)
           IF ((INT(OME1).GE.1).and.(ipatch.eq.1).and.(ilev.eq.nlmax))
     *        write (*,*) tt34-tt33, ' Einlesen  der SK-Matrix:' , ilev
        ENDIF
        ELSE
c==============================================================
c     Das lasse ich jetztz mal ganz weg!!
c==============================================================
        ENDIF
cc
c
 246   CONTINUE!Das ist Abspeichern GG Matrix: Kommt weg!

c
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C   Jetzt ist die Schurkmplement MAtrix auf invertiert. 
c   Es geht ans Loesen
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
 248    CONTINUE
      IF ((.NOT.BMATCL).AND.(ILEV.EQ.NLMIN)) THEN !Falls keine (GG-)Matrix berechnet werden muss
       CALL LCP1(DWORK(L(LBUMU)),DWORK(L(LMATRI)),NU*NU)
       CALL XLCP3(LBUPiU,LPIV,NU)
       CALL LCP1(DWORK(L(LBUMP)),DWORK(L(LSK)),NP*NP)
       CALL XLCP3(LBUPiP,LPIV1,NP)
      ENDIF
c
      IF (BSAVEM) THEN
       IF ((BMATCL).AND.(ILEV.EQ.NLMIN).AND.
     *   ((ISL.EQ.3).OR.(BLOCSL))) THEN
         write (*,*) 'wir machen Den MAtrixtrick'
         CALL LCP1(DWORK(L(LMATRI)),DWORK(L(LBUMU)),NU*NU)
         CALL XLCP3(LPIV,LBUPiU,NU)
         CALL LCP1(DWORK(L(LSK)),DWORK(L(LBUMP)),NP*NP)
         CALL XLCP3(LPIV1,LBUPiP,NP)
         BMATCL=.FALSE.
       ENDIF
      ENDIF
c 250     CONTINUE ! Hier steigen wir ein, wenn Matrix schon berechnet
C
cc        goto 250
        IF (.FALSE.) THEN
c----------------------------------
C        Hier experimentiere ich
c----------------------------------
         CALL NEWDUP (U1,U2,P,F1,F2,FP,LSK,LPIV1,L(LPATCH),L(LLDPAT),
     *             DWORK(L(LINVMA)),DWORK(L(LB11)),DWORK(L(LB22)),
     *             LFHG,NEQPAT,NELPAT,
     *             IPATCH,ILOES)
        ELSE
c----------------------------------
C        Hier ist klassisch
c----------------------------------
        CALL PAUPDN (U1,U2,P,F1,F2,FP,
     *             DWORK(L(LMATRI)),KWORK(L(LPIV)),
     *             DWORK(L(LSK)),KWORK(L(LPIV1)),
     *             VWORK(KB1), VWORK(KB2),
     *             DWORK(L(LB11)),DWORK(L(LB22)),
     *             KWORK(L(LPATCH)),KWORK(L(LLDPAT)),
     *             KWORK(L(LFHG)),NEQPAT,NELPAT,NAT,NVELO,
     *             IPATCH,IGRAD,lsk)
c          write (*,*) nu, 'schur2',neqpat
c          DO 7121 i=1,neqpat
c7121           print * , i,u1(i),U2(i),f1(i),f2(i)
C
      ENDIF
c
C
 250   continue

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
        CALL ZDISP(0,LSK,'LSK   ')
        CALL ZDISP(0,LFHG,'FHG   ')
        CALL ZDISP(0,LMATRI,'MATRIX')
        CALL ZDISP(0,Linvma,'MATRIX')
        CALL ZDISP(0,LPIV,'LPIV  ')
        CALL ZDISP(0,LPIV1,'LPIV1 ')
        CALL ZDISP(0,LB11,'LB11  ')
        CALL ZDISP(0,LB22,'LB22  ')
c
        CALL ZDISP(0,LUMFM1,'LUMFM1')
        CALL ZDISP(0,LUMFI1,'LUMFI1')
c        CALL ZDISP(0,LUMFM2,'LUMFM2')
c        CALL ZDISP(0,LUMFI2,'LUMFI2')
C
C

       ENDIF! (NELPAT.EQ.1)
c
111    CONTINUE
C
       igrad=igold! Fuer Testzwecke
      END
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     ENDE DER HAUPTROUTINE
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C=======================================================================

************************************************************************
      SUBROUTINE  VANMEZ (U1,U2,P,F1,F2,FP,
     *                    A,KACOL,KALD,B1,B2,KBCOL,KBLD,
     *                    KMBD,KVERT,KMID,KNPR,NMBD,IEL,NU)
************************************************************************
C
C-----------------------------------------------------------------------
C Purpose:  Wie VanFUL,
C nur dass auf dem Element 
c lediglich die zwei 4x4-die Matrizen invertiert werden!
C-----------------------------------------------------------------------
c$$$      INCLUDE 'common.inc'
c$$$      INCLUDE 'dwork.inc'
C-----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DOUBLE PRECISION BB1,BB2
      REAL  A,B1,B2
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4)
      PARAMETER (NNARR=299,NNLEV=9,NNAB=21,NNWORK=1)
      DIMENSION U1(*),U2(*),P(*)
      DIMENSION F1(*),F2(*),FP(*)
      DIMENSION A(*),KACOL(*),KALD(*),B1(*),B2(*),KBCOL(*),KBLD(*)
      DIMENSION KMBD(*),KVERT(NNVE,*),KMID(NNVE,*),KNPR(*)
C
C *** Local arrays for informations about one element
      DIMENSION AA1(4,4),AA2(4,4),BB1(4),BB2(4),FF1(4),FF2(4)
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
C
      INCLUDE 'block.inc'
      INCLUDE 'bouss.inc'
      SAVE
C=======================================================================
C     Block Gauss-Seidel on Schur Complement
C=======================================================================
c
      NA=KALD(NU+1)-1 !Einbau der Def-Tensor Matrix 
c
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
ccc      GOTO 111
      IA2=KALD(I+1)-1
      DO 110  IA=IA1+1,IA2
      J=KACOL(IA)
      IF ((J.NE.IU(1)).AND.(J.NE.IU(2)).AND.
     *    (J.NE.IU(3)).AND.(J.NE.IU(4))) THEN!(behaelt 4 Komponenten)
      AOFF1=DBLE(A(IA))
      AOFF2=DBLE(A(3*NA+IA))
      FF1(II)=FF1(II)-AOFF1*U1(J)
      FF2(II)=FF2(II)-AOFF2*U2(J)
      ENDIF
110   CONTINUE
      DO 111  IA=IA1,IA2
      J=KACOL(IA)
c      IF ((J.NE.IU(1)).AND.(J.NE.IU(2)).AND.
c     *    (J.NE.IU(3)).AND.(J.NE.IU(4))) THEN 
      AOFF1=DBLE(A(NA+IA))
      AOFF2=DBLE(A(2*NA+IA))
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
c$$$c       AA(jj,ii)=DBLE(A(IA))
c$$$       AA(II,JJ)=DBLE(A(IA))
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
C *** Calling the subroutine ELUPDN for element update 
      CALL  ELUPDN (U1,U2,P,IEL,IU,BDBC,AA1,AA2,BB1,BB2,FF1,FF2,FFP)
C=======================================================================
C
C
      END
c
c
c
************************************************************************
      SUBROUTINE PAUPDN(U1,U2,P,F1,F2,FP,
     *                 A1,KPIV,ALSK,KPIV1,B1,B2,B11,B22,
     *                 KPATCH,KLDPAT,KFHG,NEQPAT,NELPAT,NAT,NVELO,
     *                 IPATCH,IGRAD,lsk)
************************************************************************
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'
C
      REAL B1,B2
      DOUBLE PRECISION B11,B22
c
      DIMENSION U1(*),U2(*),P(*)
      DIMENSION F1(*),F2(*),FP(*)
      DIMENSION ALSK(NELPAT,NELPAT),A1(NVELO,NVELO)
      DIMENSION KPIV1(NELPAT),KPIV(NVELO)
      DIMENSION B11(NEQPAT,NELPAT), B22(NEQPAT,NELPAT),B1(*),B2(*)
      DIMENSION KPATCH(*),KLDPAT(*),KFHG(*)
      COMMON /ZBV/    IZBV1,IZBV2,IZBV3,IZBV4
c
C	1.) ZUORDNUNG: ALTE <=> NEUE NUMMERN
C       2.) RECHTE SEITE VEKTOR BAUEN
C	3.) B^T A^{-1}G BAUEN
C       4.) DATE UP: P= [B^TA^{-1}B]^{-1} B^T A^{-1}G
C	5.) DATE UP: U= A^{-1}G- A^{-1}BP
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      write (*,*) 'in paup', lumfma, lumfm1,lumfm2,lumfin,lumfi1,lumfi2
c           write (*,*) 'in paup 01', (kwork(l(lumfi1)-1+i),i=2,20)
      LFF =0
      LFF1=0
      LFF2=0 
      LFFP=0
c
c      call transq (ALSK,NELPAT)
c      CALL TRANSQ (A1,nvelo)
      IF (IGRAD.EQ.1) THEN
        CALL ZNEW(NEQPAT,1,LFF1,'FF1   ')
        CALL ZNEW(NEQPAT,1,LFF2,'FF2   ')
        KFF1=L(LFF1)
        KFF2=L(LFF2)
      ELSE
        CALL ZNEW(NVELO ,1,LFF ,'FF    ')
        KFF1=L(LFF)
        KFF2=L(LFF)+NEQPAT
      ENDIF
c
      DO 1 I=1,NEQPAT
           DWORK(KFF1-1+I)=F1(KFHG(I))
           DWORK(KFF2-1+I)=F2(KFHG(I))
 1    CONTINUE
c
      CALL ZNEW(NELPAT,1,LFFP,'FFP   ')
CC
      DO 24 I=1,NELPAT
          DWORK(L(LFFP)-1+I)=
c     *    -FP(KWORK(KPATCH-1+KWORK(KLDPAT-1+IPATCH)+I-1))
     *    -FP(KPATCH(KLDPAT(IPATCH)+I-1))
24    CONTINUE

c----------------------------------------------------------------
C        Hier bringe ich alte P-Werte, die im Patch nicht
C        upgedatet erden, auf die rechte Seite
c----------------------------------------------------------------
         DO 1444 II=1,NEQPAT
          IROW=KFHG(II)
           ILD=KWORK(KLDB-1+IROW)
           ILD2=KWORK(KLDB+IROW)
           IF (ILD.EQ.(ILD2-1)) GOTO 1444
C
          ICOUNT=0
          DO 1555 JJ=1,NELPAT
           JCOL=KPATCH(KLDPAT(IPATCH)-1+JJ)
            DO 1666 LAUF=ILD,(ILD2-1)
             IF (JCOL.EQ.KWORK(KCOLB-1+LAUF)) THEN
               ICOUNT=ICOUNT+1
               ILAUF=LAUF
             ENDIF
1666        CONTINUE 
1555      CONTINUE
           IF (ICOUNT.EQ.1) THEN
             IF (ILAUF.EQ.ILD) THEN
                IRIGHT=ILD2-1
             ELSE
                IRIGHT=ILD
             ENDIF
C
             PDRUCK=P(KWORK(KCOLB-1+IRIGHT))
C
             DWORK(KFF1-1+II)=DWORK(KFF1-1+II)-
     *                         DBLE(B1(IRIGHT))*PDRUCK
             DWORK(KFF2-1+II)=DWORK(KFF2-1+II)-
     *                         DBLE(B2(IRIGHT))*PDRUCK
           ENDIF
1444     CONTINUE
C-----------------------------------------------------------------------
C        Hier werden die alten U-Werte, Ddie im Patch nicht  
C       upgedatet werden, auf die rechte Seite gebracht
C-----------------------------------------------------------------------

      DO 71 II=1,NEQPAT
C
         IROW=KFHG(II)
         ILD =KWORK(KLDA-1+IROW)
         ILD2=KWORK(KLDA+IROW)
C
         DO 72 LAUF=ILD,(ILD2-1)
            ILAUF=0
            NRINA=KWORK(KCOLA-1+LAUF)
            DO 73 JJ=1,NEQPAT                          
              JCOL=KFHG(JJ)
              IF (JCOL.EQ.NRINA) ILAUF=1
73          CONTINUE  
c
            CALL VKADD (AVK1,AVK2)
c
                 AMat1=AVK1*DBLE(VWORK(KA1-1+LAUF))+
     *                 AVK2*DBLE(VWORK(KA1-1+4*NA+LAUF))
                 AMat2=AVK1*DBLE(VWORK(KA1-1+NA+LAUF))
                 AMat3=AVK1*DBLE(VWORK(KA1-1+2*NA+LAUF))
                 AMat4=AVK1*DBLE(VWORK(KA1-1+3*NA+LAUF))+
     *                 AVK2*DBLE(VWORK(KA1-1+4*NA+LAUF))
c
            IF (ILAUF.EQ.0) THEN
              DWORK(KFF1-1+II)=DWORK(KFF1-1+II)-
     *                AMat1*U1(NRINA)-AMat2*U2(NRINA)
              DWORK(KFF2-1+II)=DWORK(KFF2-1+II)-
     *                AMat3*U1(NRINA)-AMat4*U2(NRINA)
c
c$$$              DWORK(KFF1-1+II)=DWORK(KFF1-1+II)-
c$$$     *                DBLE(VWORK(KA1-1+LAUF))*U1(NRINA)
c$$$              DWORK(KFF2-1+II)=DWORK(KFF2-1+II)-
c$$$     *                DBLE(VWORK(KA1-1+LAUF))*U2(NRINA)
            ENDIF                
72       CONTINUE
71    CONTINUE
c
c           write (*,*) 'in paup 06', (kwork(l(lumfi1)-1+i),i=2,20)

c------------------------------------------------------
C
C     Der Aufbau von
C     FP + B1^T A^{-1}F1 + B2^T A^{-1}F2
C     ergibt die rechte Seite fuer die Schurkomplement Gleichung
C
c------------------------------------------------------
      CALL ZTIME (TRSAU0)
       IF (IGRAD.EQ.1) THEN
          CALL RESLU (NEQPAT,A1,DWORK(L(LFF1)),KPIV,LUMFM1,LUMFI1,
     *                KEEP1,CNTL1,ICNTL1)
          CALL RESLU (NEQPAT,A1,DWORK(L(LFF2)),KPIV,LUMFM1,LUMFI1,
     *                KEEP1,CNTL1,ICNTL1)
c	 CALL DGESL (A1,NEQPAT,NEQPAT,KPIV,DWORK(L(LFF1)),0)
c	 CALL DGESL (A1,NEQPAT,NEQPAT,KPIV,DWORK(L(LFF2)),0)
       ELSE

          CALL RESLU (NVELO,A1,DWORK(L(LFF)),KPIV,LUMFM1,LUMFI1,
     *                KEEP1,CNTL1,ICNTL1)
c         CALL DGESL (A1,NVELO,NVELO,KPIV,DWORK(L(LFF)),0)
       ENDIF
      DO 1234 II=1,NELPAT
        DO 2345 JJ=1,NEQPAT
          DWORK(L(LFFP)+II-1)=DWORK(L(LFFP)+II-1)+
     *        B11(JJ,II)*DWORK(KFF1-1+JJ)+
     *        B22(JJ,II)*DWORK(KFF2-1+JJ)
2345    CONTINUE
1234  CONTINUE
C
      CALL ZTIME (TRSAU1)
      TRSAUF=TRSAU1-TRSAU0
      IF ((INT(OME1).GE.1).and.(ipatch.eq.1).and.(ilev.eq.nlmax))
     *      print *, TRSAUF, ' RHS-Aufbau fuer p' 
c------------------------------------------------------
C
C
       IFC1=IFACTO
       IFACTO=1
      LPP=LFFP
      CALL ZTIME (TRESP0)
      CALL RESLU (NELPAT,ALSK,DWORK(L(LPP)),KPIV1,LUMFM2,LUMFI2,
     *            KEEP2,CNTL2,ICNTL2)
      CALL ZTIME (TRESP1)
      TRESP=TRESP1-TRESP0
      IFACTO=IFC1
      IF ((INT(OME1).GE.1).and.(ipatch.eq.1).and.(ilev.eq.nlmax)) THEN
      print *, TRESP, ' RESLU von P ' 
      ENDIF
c     CALL DGESL (ALSK,NELPAT,NELPAT,KPIV1,DWORK(L(LPP)),0)


C
C    Jetzt steht der upgedatete Druck in LPP
C
C
      LUU1=0
      LUU2=0
      LUU=0
      IF (IGRAD.EQ.1) THEN
        CALL ZNEW (NEQPAT,1,LUU1,'LUU1  ')
        CALL ZNEW (NEQPAT,1,LUU2,'LUU2  ')
         KUU1=L(Luu1)
         Kuu2=L(Luu2)
      ELSE
        CALL ZNEW (NVELO,1,LUU,'LUU   ')
         Kuu1=L(Luu)
         Kuu2=L(luu)+NEQPAT
      ENDIF
c
c
c      write (*,*) 'druck', (j,dwork(l(lpp)-1+j),j=1,nelpat)
c      write (*,*) 'B1' 
c      do 99 i=1,neqpat
c         do 99 j=1,nelpat
c99         write (*,*) i,j,b11(i,j)
c
      call ztime (ttsep1)!wwwwww
c
      DO 234 II=1,NEQPAT
        DO 345 JJ=1,NELPAT
          DWORK(KUU1+II-1)=DWORK(KUU1+II-1)+
     *        B11(II,JJ)*DWORK(L(LPP)-1+JJ)
          DWORK(KUU2+II-1)=DWORK(KUU2+II-1)+
     *        B22(II,JJ)*DWORK(L(LPP)-1+JJ)
          if (abs(DWORK(KUU2+II-1)).le.1d-11)
     * DWORK(KUU2+II-1)=0d0
345      CONTINUE
234   CONTINUE

c      write (*,*) 'oi0',(i,dwork(kuu1-1+i),
c     *  dwork(kff1-1+i),'\n',i=1,neqpat)
c
C
      IF (IGRAD.eq.1) THEN
         CALL RESLU (NEQPAT,A1,DWORK(L(LUU1)),KPIV,LUMFM1,LUMFI1,
     *               KEEP1,CNTL1,ICNTL1)
         CALL RESLU (NEQPAT,A1,DWORK(L(LUU2)),KPIV,LUMFM1,LUMFI1,
     *               KEEP1,CNTL1,ICNTL1)
c        CALL DGESL (A1,NEQPAT,NEQPAT,KPIV,DWORK(L(LUU1)),0)
c        CALL DGESL (A1,NEQPAT,NEQPAT,KPIV,DWORK(L(LUU2)),0)
      ELSE
         CALL RESLU (NVELO,A1,DWORK(L(LUU)),KPIV,LUMFM1,LUMFI1,
     *              KEEP1,CNTL1,ICNTL1)
c        CALL DGESL (A1,NVELO,NVELO,KPIV,DWORK(L(LUU)),0)
      ENDIF
C
c      write (*,*) 'oi1',(i,dwork(kuu1-1+i),
c     *  dwork(kff1-1+i),'\n',i=1,neqpat)
c
      DO 11 I=1,NEQPAT
        IF (KWORK(L(LNPR)-1+KFHG(I)+NVT).NE.0) THEN
             DWORK(KUU1-1+I)=0D0
             DWORK(KUU2-1+I)=0D0
        ELSE
             DWORK(KUU1-1+I)=
     *       DWORK(KFF1-1+I)-DWORK(KUU1-1+I)
             DWORK(KUU2-1+I)=
     *       DWORK(KFF2-1+I)-DWORK(KUU2-1+I)
        ENDIF
11    CONTINUE
C
c      write (*,*) 'oi2',(i,dwork(kuu1-1+i),
c     *    dwork(kff1-1+i),'\n',i=1,neqpat)
C
c------------------------------------------------------
C   Jetzt hat man alle Werte, man muss sie nur noch
C   auf den eigentlichen Loesungsvektor kopieren
c------------------------------------------------------
C
      DO 21 I=1,NEQPAT
          IF (KWORK(L(LNPR)-1+KFHG(I)+NVT).NE.0) GOTO 21
             U1(KFHG(I))=DWORK(Kuu1-1+I)
             U2(KFHG(I))=DWORK(Kuu2-1+I)
21    CONTINUE
C
      DO 22 I=1,NELPAT
            P(KPATCH(KLDPAT(IPATCH)+I-1))=DWORK(L(LPP)-1+I)
22    CONTINUE
C
      call ztime (ttsep0)
       IF((INT(OME1).GE.1).and.(ipatch.eq.1).and.(ilev.eq.nlmax)) THEN
      write (*,*) ttsep0-ttsep1, ' der rest der U-Bestimmung'
      print *, '==========================================='
      ENDIF
c
      IF (IGRAD.eq.1) THEN
      CALL ZDISP (0,LFF1,'FF1   ')
      CALL ZDISP (0,LFF2,'FF2   ')
      CALL ZDISP (0,LUU1,'UU1   ')
      CALL ZDISP (0,LUU2,'UU2   ')
      ELSE
      CALL ZDISP (0,LFF ,'FF    ')
      CALL ZDISP (0,LUU ,'UU    ')
      ENDIF
      CALL ZDISP (0,LFFP,'FFP   ')
c
      END
c

c
************************************************************************
      SUBROUTINE NEWDUP(U1,U2,P,F1,F2,FP,LSK,LPIV1,KPATCH,KLDPAT,
     *         AINV,B11,B22,LFHG,NEQPAT,NELPAT,IPATCH,ILOES)
************************************************************************
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
C
      Double Precision B11,B22,AINV
      DIMENSION U1(*),U2(*),P(*)
      DIMENSION F1(*),F2(*),FP(*)
      DIMENSION AINV(*),B11(*),B22(*)
C	1.) ZUORDNUNG: ALTE <=> NEUE NUMMERN
C       2.) RECHTE SEITE VEKTOR BAUEN
C	3.) B^T A^{-1}G BAUEN
C       4.) DATE UP: P= [B^TA^{-1}B]^{-1} B^T A^{-1}G
C	5.) DATE UP: U= A^{-1}G- A^{-1}BP
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      LFF1=0
      LFF2=0 
      LFFP=0
      CALL ZNEW(NEQPAT,1,LFF1,'FF1   ')
      CALL ZNEW(NEQPAT,1,LFF2,'FF2   ')
      CALL ZNEW(NELPAT,1,LFFP,'FFP   ')
C
C 
      DO 1 I=1,NEQPAT
         DWORK(L(LFF1)-1+I)=F1(KWORK(L(LFHG)-1+I))
         DWORK(L(LFF2)-1+I)=F2(KWORK(L(LFHG)-1+I))
1     CONTINUE
C
      DO 24 I=1,NELPAT
      IF (ILOES.EQ.1) THEN
          DWORK(L(LFFP)-1+I)=
     *    -DWORK(KFP-1+KWORK(KPATCH-1+KWORK(KLDPAT-1+IPATCH)+I-1))
      ELSE
          DWORK(L(LFFP)-1+I)=
     *    -FP(KWORK(KPATCH-1+KWORK(KLDPAT-1+IPATCH)+I-1))
      ENDIF
24    CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        HIER BRINGE ICH ALTE P-WERTE, DIE IM PATCH NICHT
C        UPGEDATED WERDEN, AUF DIE RECHTE SEITE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         DO 1444 II=1,NEQPAT
          IROW=KWORK(L(LFHG)+II-1)
           ILD=KWORK(KLDB-1+IROW)
           ILD2=KWORK(KLDB+IROW)
           IF (ILD.EQ.(ILD2-1)) GOTO 1444
C
          ICOUNT=0
          DO 1555 JJ=1,NELPAT
           JCOL=KWORK(KPATCH-1+KWORK(KLDPAT-1+IPATCH)-1+JJ)
            DO 1666 LAUF=ILD,(ILD2-1)
             IF (JCOL.EQ.KWORK(KCOLB-1+LAUF)) THEN
               ICOUNT=ICOUNT+1
               ILAUF=LAUF
             ENDIF
1666        CONTINUE 
1555      CONTINUE
           IF (ICOUNT.EQ.1) THEN
             IF (ILAUF.EQ.ILD) THEN
                IRIGHT=ILD2-1
             ELSE
                IRIGHT=ILD
             ENDIF
C
             IF (ILOES.EQ.1) THEN
                PDRUCK=DWORK(KP-1+KWORK(KCOLB-1+IRIGHT))
             ELSE
                PDRUCK=P(KWORK(KCOLB-1+IRIGHT))
             ENDIF
C
             DWORK(L(LFF1)-1+II)=DWORK(L(LFF1)-1+II)-
     *                         DBLE(VWORK(KB1-1+IRIGHT))*
     *                         PDRUCK
             DWORK(L(LFF2)-1+II)=DWORK(L(LFF2)-1+II)-
     *                         DBLE(VWORK(KB2-1+IRIGHT))*
     *                         PDRUCK
           ENDIF
1444     CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C        HIER BRINGE ICH ALTE U-WERTE, DIE IM PATCH NICHT
C        UPGEDATED WERDEN, AUF DIE RECHTE SEITE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DO 71 II=1,NEQPAT
C
            IROW=KWORK(L(LFHG)+II-1)
            ILD=KWORK(KLDA-1+IROW)
            ILD2=KWORK(KLDA+IROW)
C
            DO 72 LAUF=ILD,(ILD2-1)
                   ILAUF=0
                   NRINA=KWORK(KCOLA-1+LAUF)
                   DO 73 JJ=1,NEQPAT                          
                           JCOL=KWORK(L(LFHG)+JJ-1)
                           IF (JCOL.EQ.NRINA) ILAUF=1
73                 CONTINUE  
                   IF (ILAUF.EQ.0) THEN
                       DWORK(L(LFF1)-1+II)=
     *                 DWORK(L(LFF1)-1+II)-DBLE(VWORK(KA1-1+LAUF))*
     *                 U1(NRINA)
                       DWORK(L(LFF2)-1+II)=
     *                 DWORK(L(LFF2)-1+II)-DBLE(VWORK(KA1-1+LAUF))*
     *                 U2(NRINA)
                   ENDIF                
72          CONTINUE
71    CONTINUE
C------------------------------------------------------
C      Mache A^{-1} g
C------------------------------------------------------
       LSOL1=0
       LSOL2=0
         CALL ZNEW (NEQPAT,1,LSOL1,'LSOL1a')
         CALL ZNEW (NEQPAT,1,LSOL2,'LSOL2b')       
C
         DO 210 IZEI=1,NEQPAT
         DO 210 ISPAL=1,NEQPAT
               dwork(l(lsol1)-1+izei)=
     *            dwork(l(lsol1)-1+izei) +
     *            ainv(NEQPAT*(izei-1)+ispal) *
     *            DWORK(L(LFF1)-1+ispal)
c
               dwork(l(lsol2)-1+izei)=
     *            dwork(l(lsol2)-1+izei) +
     *            ainv(NEQPAT*(izei-1)+ispal) *
     *            DWORK(L(LFF2)-1+ispal)
 210           continue
c
c

         CALL ZCPY (LSOL1,'LSOL1 ',LFF1,'Lff1  ')
         CALL ZCPY (LSOL2,'LSOL2 ',LFF2,'LFF2  ')

c
         CALL ZDISP (0,LSOL1,'LSPO1')
         CALL ZDISP (0,LSOL2,'LSPO2')
C------------------------------------------------------
C      Mache B^{T} (A^{-1} g)
C------------------------------------------------------
C
C	Ergibt die rechte Seite fuer die Schurkomplement Gl.
C
      DO 1234 II=1,NELPAT
        DO 2345 JJ=1,NEQPAT
          DWORK(L(LFFP)+II-1)=DWORK(L(LFFP)+II-1)+
     *        B11(NEQPAT*(II-1)+JJ)*DWORK(L(LFF1)-1+JJ)+
     *        B22(NEQPAT*(II-1)+JJ)*DWORK(L(LFF2)-1+JJ)
2345    CONTINUE
1234  CONTINUE

C
C
      LPP=LFFP
      print *, 'hey, hier kein Loeser!!!'
c      CALL DGESL (DWORK(L(LSK)),NELPAT,NELPAT,KWORK(L(LPIV1)),
c     *              DWORK(L(LPP)),0)

C
C		JETZT STEHT DER UPGEDATETE DRUCK IN LPP
C
C
      LUU1=0
      LUU2=0
      CALL ZNEW (NEQPAT,1,LUU1,'LUU1  ')
      CALL ZNEW (NEQPAT,1,LUU2,'LUU2  ')
      DO 234 II=1,NEQPAT
        DO 345 JJ=1,NELPAT
          DWORK(L(LUU1)+II-1)=DWORK(L(LUU1)+II-1)+
     *        B11(NEQPAT*(JJ-1)+II)*DWORK(L(LPP)-1+JJ)
          DWORK(L(LUU2)+II-1)=DWORK(L(LUU2)+II-1)+
     *        B22(NEQPAT*(JJ-1)+II)*DWORK(L(LPP)-1+JJ)
345      CONTINUE
234   CONTINUE
C

       LSOL1=0
       LSOL2=0
         CALL ZNEW (NEQPAT,1,LSOL1,'LSOL1c')
         CALL ZNEW (NEQPAT,1,LSOL2,'LSOL2d')
c
         DO 310 IZEI=1,NEQPAT
         DO 310 ISPAL=1,NEQPAT
               dwork(l(lsol1)-1+izei)=
     *            dwork(l(lsol1)-1+izei) +
     *            ainv(NEQPAT*(izei-1)+ispal) *
     *            DWORK(L(LUU1)-1+ispal)
c
               dwork(l(lsol2)-1+izei)=
     *            dwork(l(lsol2)-1+izei) +
     *            ainv(NEQPAT*(izei-1)+ispal) *
     *            DWORK(L(LUU2)-1+ispal)
 310           continue
c
      CALL ZDISP (0,LUU1,'UU1   ')
      CALL ZDISP (0,LUU2,'UU2   ')     
         CALL ZCPY (LSOL1,'LSOL1 ',LUU1,'Luu1  ')
         CALL ZCPY (LSOL2,'LSOL2 ',LUU2,'Luu2  ')
c
         CALL ZDISP (0,LSOL1,'LSOL1c')
         CALL ZDISP (0,LSOL2,'LSOL2d')
c
      DO 11 I=1,NEQPAT
        IF (KWORK(L(LNPR)-1+KWORK(L(LFHG)-1+I)+NVT).NE.0) THEN
             DWORK(L(LUU1)-1+I)=0D0
             DWORK(L(LUU2)-1+I)=0D0
        ELSE
             DWORK(L(LUU1)-1+I)=
     *       DWORK(L(LFF1)-1+I)-DWORK(L(LUU1)-1+I)
             DWORK(L(LUU2)-1+I)=
     *       DWORK(L(LFF2)-1+I)-DWORK(L(LUU2)-1+I)
        ENDIF
11    CONTINUE
C

C
C   SO, JETZT HABE ICH ALLE WERTE, JETZT DATE ICH DEN EIGENTLICHEN
C   LOESUNGSVEKTOR UP!
C
      DO 21 I=1,NEQPAT
          IF (KWORK(L(LNPR)-1+KWORK(L(LFHG)-1+I)+NVT).NE.0) GOTO 21
             U1((KWORK(L(LFHG)-1+I)))=DWORK(L(LUU1)-1+I)
             U2((KWORK(L(LFHG)-1+I)))=DWORK(L(LUU2)-1+I)
21    CONTINUE
C
      DO 22 I=1,NELPAT
            P(KWORK(KPATCH-1+KWORK(KLDPAT-1+IPATCH)+I-1))
     *                           =DWORK(L(LPP)-1+I)
22    CONTINUE
C
      CALL ZDISP (0,LFF1,'FF1   ')
      CALL ZDISP (0,LFF2,'FF2   ')
      CALL ZDISP (0,LFFP,'FFP   ')
      CALL ZDISP (0,LUU1,'UU1   ')
      CALL ZDISP (0,LUU2,'UU2   ')
C      CALL ZDISP (0,LPP ,'PP    ')
      END
c
C------------------------------------------------------------------------
      subroutine DUMFFA (N,NNE,AMAT,INDEKs,KEEP,CNTL,ICNTL)
C------------------------------------------------------------------------
c
      PARAMETER (NNARR=299, NNWORK=1, NNLEV=9)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)

C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      DIMENSION AMAT(*),INDEKs(*)
      DIMENSION KEEP(20), ICNTL(20),CNTL(10)
      INCLUDE 'bouss.inc'
      INCLUDE 'block.inc'
      INCLUDE 'dwork.inc'
      SAVE
c  
c     DDIM: Dimension des Problems
c     # Matrixeintraege <> 0
c LVALUE : Laenge des Vektors "Matrix"
C LINDEX : Laenge des Vektors Indexmatrix
C MATRIX :enthaelt anfangs die <>0 Eintraege der Matrix, danach die Faktorisierung
c INDEXMATRIX: Enthaelt Indizes (i,j) der Matrixelemente, erst alle i, dann alle j
C  B ist die rechte Seite
C X : Vektor, auf den die Loesung kommt
c W    :Hilfsvektor
C  
C I.3      LGS: GGU * DELU = MMG
c
c$$$        CALL UMD21I (KEEP, CNTL, ICNTL)
c         ICNTL (3) = 4

         LVALUE=n*N*IZBV4*2
         LINDEX=2*lvalue
c         write (*,*) 'lumfma' ,nne,n, lvalue
         CALL UMD2FA (N, NNE, 0, .FALSE., LVALUE, LINDEX, 
     *                AMAT,INDEKS,
c     $               DWORK(L(LUMFMA)), KWORK(L(LUMFIN)),
     *               KEEP, CNTL, ICNTL, INFO, RINFO)
c         write (*,*)info(1),  'n,nne2', (INDEKS(i),i=2,20)
        IF (INFO (1) .LT. 0) STOP
c        ICNTL (3) = 2

C------------------------------------------------------------------------      
      end
C------------------------------------------------------------------------
      subroutine DUMFSL (N,AMAT,INDEks,B,KEEP,ICNTL,CNTL)
C------------------------------------------------------------------------
c
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299, NNWORK=1, NNLEV=9)
      DOUBLE PRECISION B
      INTEGER INDEKS
      DIMENSION B(*),AMAT(*),INDEks(*)
      DIMENSION KEEP(20),ICNTL(20),CNTL(10)


C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'

      SAVE
c
c
c     DDIM: Dimension des Problems
c     # Matrixeintraege <> 0
c LVALUE : Laenge des Vektors "Matrix"
C LINDEX : Laenge des Vektors Indexmatrix
C MATRIX :enthaelt anfangs die <>0 Eintraege der Matrix, danach die Faktorisierung
c INDEXMATRIX: Enthaelt Indizes (i,j) der Matrixelemente, erst alle i, dann alle j
C  B ist die rechte Seite
C X : Vektor, auf den die Loesung kommt
c W    :Hilfsvektor
C  
C I.3      LGS: GGU * DELU = MMG
c

C Solve Ax = b and print solution.
c       write (*,*) 'in dumfsl 1'!,lumfma,lumfin,cntl,icntl,keep
c       write (*,*) (amat(i),INDEKS(i),i=2,20)
        LXHELP=0
        LWHELP=0
        CALL ZNEW (N,1,LXHELP,'LXHELP')
        CALL ZNEW (4*N,1,LWHELP,'LWHELP')
c
c         print *,lumfma, 'ping', lindex,lvalue
        CALL UMD2SO (N, 0, .FALSE., LVALUE, LINDEX, ! DWORK(L(LUMFMA)),
     $               AMAT,INDEks,! KWORK(L(LUMFIN)),
     $               KEEP, B, DWORK(L(LxHELP)), DWORK(L(LwHELP)),
     $               CNTL, ICNTL, INFO, RINFO)

        DO 10 i=1,n
c          write (*,*) i,DWORK(L(LXHELP)-1+i)
         B(i)=   DWORK(L(LXHELP)-1+i)
 10     continue
         xmaxde=0d0
         do 20 i=1,n
            IF (dwork(l(lwhelp)-1+i).ge.xmaxde)
     *     xmaxde=dwork(l(lwhelp)-1+i)
 20         continue
c            write (*,*) 'xmaxde', xmaxde
        IF (INFO (1) .LT. 0) STOP
C------------------------------------------------------------------------      
        CALL ZDISP (0,LXHELP,'LXHELP')
        CALL ZDISP (0,LWHELP,'LWHELP')
      end

c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE FACTOR (N,AMAT,PIVEC,
     *                   KEEP,CNTL,ICNTL,LUMFMA,LUMFIN)
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    
c     calculates LU-Factorization of Matrix AMAT
c     Returns factorization on AMAT and pivotvector on PIVEC
c     DGTRF1 IST DGETRF SELBSTCOMPILIERT
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'bouss.inc'
      INCLUDE 'block.inc'

c
      INTEGER PIVEC,INFO2,INFO,IIF, LUMFMA, LUMFIN
      DIMENSION AMAT(N,N), PIVEC(N)
      DIMENSION KEEP(*),ICNTL(*),CNTL(*)
C
      SAVE
c
       CALL TRANSQ (AMAT,N)
c
       IF (IFACTO.EQ.0) 
     *    CALL DGETRF (N,N,AMAT,N,PIVEC,INFO2)
c     *    CALL DGEFA (AMAT,N,N,PIVEC,INFO)
       IF (IFACTO.EQ.1) 
     *    CALL DGETRF (N,N,AMAT,N,PIVEC,INFO2)
       IF (IFACTO.eq.2)  THEN 
         NNE=0 !# der MAtrixeintraege <>0
         DO 10 i=1, N
            Do 10 j=1,n

               IF (ABS(AMAT(i,j)).ge.1d-12) THEN
                  NNE=NNe+1
                  DWORK (L(LUMFMA)-1+NNE)=AMAT(i,j)
                  KWORK (L(LUMFIN)-1+NNE)=i
                  KWORK (L(LUMFIN)-1+NNE+N*N)=j
c
cc                  KWORK (L(LUMFIN)-1+NNE)=j
cc                  KWORK (L(LUMFIN)-1+NNE+N*N)=i
                  ENDIF
 10            continue    
         DO 20 K=1,NNE
            KWORK(L(LUMFIN)-1+NNE+k)=KWORK(L(LUMFIN)-1+N*N+k)
 20         CONTINUE!KWORK(L(LUMFIN)-1+N*N+k)=0
        CALL UMD21I (KEEP, CNTL, ICNTL)
c            ICNTL(3)=4

        CALL DUMFFA (N,NNE,DWORK(L(LUMFMA)),KWORK(L(LUMFIN)),
     *  KEEP,CNTL,ICNTL)

       ENDIF!IFACTO=2
c
c
c
      END
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE RESLU (N,AMAT,BVEC,PIVEC,LUMFMA,LUMFIN,keep,cntl,icntl)
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Solves AMAT * X = BVEC , where AMAT is in LU decomposition
c     DGTRS1 IST DGETRS SELBSTCOMPILIERT
c$$$ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'bouss.inc'
      INCLUDE 'block.inc'

      INTEGER PIVEC
      DOUBLE PRECISION BVEC
      DIMENSION AMAT(N,N), PIVEC(N),BVEC(N,*)!,BVEC(N,1)
      DIMENSION KEEP(20),ICNTL(20),CNTL(10)

      SAVE
c
Cc
c      write (*,*) 'vor dumfsl in reslu',lumfma,lumfin
c      write (*,*) (dwork(l(lumfma)-1+i),kwork(l(lumfin)-1+i),i=2,20)
       IF (IFACTO.EQ.0) 
     *   CALL DGETRS ('N',N,n,AMAT,N,PIVEC,BVEC,N,INFO)
c     *   CALL DGESL (AMAT,N,N,PIVEC,BVEC,0)
       IF (IFACTO.EQ.1) 
     *   CALL DGETRS ('N',N,1,AMAT,N,PIVEC,BVEC,N,INFO)
       IF (IFACTO.eq.2) 
     *   CALL DUMFSL (N,DWORK(L(LUMFMA)),KWORK(L(LUMFIN)),BVEC,
     *               KEEP,ICNTL,CNTL)
c      write (*,*) 'nach dumfsl in reslu'
c
      END
      
************************************************************************
      SUBROUTINE BTBLD(LB11,LB22,NELPAT,NEQPAT,LFHG,LPATCH,LLDPAT,
     *                              IPATCH,KLDB,KCOLB,KB1,KB2)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNWORK=1,NNARR=299)
      INCLUDE 'dwork.inc'
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      SAVE
C
C
C	ERSTELLT DIE LOKALEN MATRIZEN B1^{T}, B2^{T}
C	AUS DEN AUSGANGSMATRIZEN.
C	BEACHTE: ES WERDEN DIE TRANSPONIERTEN ERSTELLT,
C	(ENTSPRICHT DER DIVERGENZ!!!)
C
C
         DO 1444 II=1,NEQPAT
          IROW=KWORK(L(LFHG)+II-1)
           ILD=KWORK(KLDB-1+IROW)
           ILD2=KWORK(KLDB+IROW)
C
C 
          ICOUNT=0
          DO 1555 JJ=1,NELPAT
           JCOL=KWORK(L(LPATCH)-1+KWORK(L(LLDPAT)-1+IPATCH)-1+JJ)
C             FINDE ELEMENT B(IROW,JCOL) DER GROSSEN MATRIX
C
            DO 1666 LAUF=ILD,(ILD2-1)
             IF (JCOL.EQ.KWORK(KCOLB-1+LAUF)) THEN
               ICOUNT=ICOUNT+1
               ILAUF=LAUF
               DWORK(L(LB11)-1+NEQPAT*(JJ-1)+II)=DBLE(VWORK(KB1-1+LAUF))
               DWORK(L(LB22)-1+NEQPAT*(JJ-1)+II)=DBLE(VWORK(KB2-1+LAUF))
            ENDIF
1666        CONTINUE
1555       CONTINUE
1444      CONTINUE
c
      END
C
C
CC 
************************************************************************
      SUBROUTINE SKBLD  (B11,B22,A,PIV,ALSK,NELPAT,NEQPAT,
     *                   NVELO,IGRAD,IPATCH,NPATCH)
************************************************************************
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNWORK=1,NNARR=299,NNLEV=9)     
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      INTEGER PIV
      DOUBLE PRECISION B11,B22
      CHARACTER CFILE*6
      DIMENSION  PIV(NVELO)            
      DIMENSION  B11(NEQPAT,NELPAT),B22(NEQPAT,NELPAT)
      DIMENSION  A(NVELO,NVELO), ALSK(NELPAT,NELPAT)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      SAVE
C
c=================================================================
C  		AUFBAU DER SCHURKOMPLEMENT MATRIX
c=================================================================
C

      TRESLU=0
      TBTAB=0
      call ZTIME (TSBLD0)

c
      
c====================================================================
c     sparsameres Abspeichern von B^T
c====================================================================
      LBT1=0
      LBT2=0
      LBT1I=0
      LBT1J=0
      LBT2I=0
      LBT2J=0
      CALL ZNEW (NELPAT*NEQPAT,1,LBT1,'LBT1  ')
      CALL ZNEW (NELPAT*NEQPAT,3,LBT1I,'LBT1I ')
      CALL ZNEW (NELPAT*NEQPAT,3,LBT1J,'LBT1J ')
      NNB1=0
      NNB2=0
      Do 100 IZEI=1,NEQPAT! B^T hat NEQPAT Zeilen
         DO 111 JSPAL=1,NELPAT! und NELPAT Spalten
            IF (ABS(B11(IZEI,JSPAL)).ge.1D-8) THEN
              NNB1=NNB1+1
              DWORK(L(LBT1)-1+NNB1)=B11(IZEI,JSPAL)
              KWORK(L(LBT1I)-1+NNB1)=IZEI
              KWORK(L(LBT1J)-1+NNB1)=JSPAL              
            ENDIF
111        continue
 100        continue
      CALL ZDISP(NNB1,LBT1,'LBT1  ')
      CALL ZDISP(NNB1,LBT1I,'LBT1I ')
      CALL ZDISP(NNB1,LBT1J,'LBT1J ')
c
      CALL ZNEW (NELPAT*NEQPAT,1,LBT2,'LBT2  ')
      CALL ZNEW (NELPAT*NEQPAT,3,LBT2I,'LBT2I ')
      CALL ZNEW (NELPAT*NEQPAT,3,LBT2J,'LBT2J ')

      Do 1000 IZEI=1,NEQPAT! B^T hat NEQPAT Zeilen
         DO 1110 JSPAL=1,NELPAT! und NELPAT Spalten
            IF (ABS(B22(IZEI,JSPAL)).ge.1D-8) THEN
              NNB2=NNB2+1
              DWORK(L(LBT2)-1+NNB2)=B22(IZEI,JSPAL)
              KWORK(L(LBT2I)-1+NNB2)=IZEI
              KWORK(L(LBT2J)-1+NNB2)=JSPAL              
            ENDIF
 1110    continue
 1000 continue 
c           
      CALL ZDISP(NNB2,LBT2,'LBT2  ')
      CALL ZDISP(NNB2,LBT2I,'LBT2I ')
      CALL ZDISP(NNB2,LBT2J,'LBT2J ')
        
c====================================================================
      IF (IGRAD.EQ.1) THEN
         CALL ZNEW (NEQPAT,1,LSPAL1,'LSPL1A')
         CALL ZNEW (NEQPAT,1,LSPAL2,'LSPL2A')
      ELSE
         CALL ZNEW (NVELO,1,LSPAL,'LSPALT')
      ENDIF
c-----------------------------------------------------------------
C	Schleife ueber alle Spalten von B1/B2
c-----------------------------------------------------------------
c
      CALL ZTIME (TVORL1)
      DO 21 ISPALT=1,NELPAT 

         CALL ZTIME (TRELU0)
         IF (IGRAD.EQ.1) THEN
            CALL LCP1(B11(1,ISPALT),DWORK(L(LSPAL1)),NEQPAT)
            CALL LCP1(B22(1,ISPALT),DWORK(L(LSPAL2)),NEQPAT)
c

            CALL RESLU (NEQPAT  ,A,DWORK(L(LSPAL1)),PIV,LUMFM1,LUMFI1,
     *                  KEEP1,CNTL1,ICNTL1)
            CALL RESLU (NEQPAT  ,A,DWORK(L(LSPAL2)),PIV,LUMFM1,LUMFI1,
     *                  KEEP1,CNTL1,ICNTL1)
         ELSE
            CALL LCP1(B11(1,ISPALT),DWORK(L(LSPAL))       ,NEQPAT)
            CALL LCP1(B22(1,ISPALT),DWORK(L(LSPAL)+NEQPAT),NEQPAT)
c
            CALL RESLU (2*NEQPAT,A,DWORK(L(LSPAL)) ,PIV,LUMFM1,LUMFI1,
     *                  KEEP1,CNTL1,ICNTL1)
         ENDIF
c
         CALL ZTIME (TRELU1)
         TRESLU=TRESLU+(TRELU1-TRELU0)
c-----------------------------------------------------------------
C	         In den Vektoren LSPAL* steht jetzt
c                A^{-1}B_{ISPALT}1/2,
C	         also die Loesung von A*X=B, wobei B die
C                ISPALT-te Spalte von B1/2 IST
c-----------------------------------------------------------------
c
c hier wieder dieselbe Frage:
c nota bene: FORTRAN speichert Matrizen SPALTENWEISE,
c     L(LMATRI) hingegen ZEILENWEIESE. Daher hier Alsk(j,i) statt Alsk(i,j)
         CALL ZTIME (TBTMA0)
c
         IF (IGRAD.EQ.1) THEN
         Do 500 ILAUF=1,NNB1
            IZEIL=KWORK(L(LBT1j)-1+ILAUF)
            II   =KWORK(L(LBT1i)-1+ILAUF)
               ALSK(ISPALT,IZEIL)=ALSK(ISPALT,IZEIL)+
     *          DWORK(L(LBT1)-1+ILAUF)*DWORK(l(lspal1)-1+ii)
c              write (*,*) II,IZEIL,B11(II,IZEIL),DWORK(L(LBT1)-1+ILAUF)
 500        CONTINUE
         Do 550 ILAUF=1,NNB2
            IZEIL=KWORK(L(LBT2j)-1+ILAUF)
            II   =KWORK(L(LBT2i)-1+ILAUF)
               ALSK(ISPALT,IZEIL)=ALSK(ISPALT,IZEIL)+
     *          DWORK(L(LBT2)-1+ILAUF)*DWORK(l(lspal2)-1+ii)
 550        CONTINUE
c
         ELSE
c
         Do 600 ILAUF=1,NNB1
            IZEIL=KWORK(L(LBT1j)-1+ILAUF)
            II   =KWORK(L(LBT1i)-1+ILAUF)
               ALSK(ISPALT,IZEIL)=ALSK(ISPALT,IZEIL)+
     *          DWORK(L(LBT1)-1+ILAUF)*DWORK(l(lspal)-1+ii)
 600        CONTINUE
            Do 650 ILAUF=1,NNB2
            IZEIL=KWORK(L(LBT2j)-1+ILAUF)
            II   =KWORK(L(LBT2i)-1+ILAUF)
               ALSK(ISPALT,IZEIL)=ALSK(ISPALT,IZEIL)+
     *         DWORK(L(LBT2)-1+ILAUF)*DWORK(l(lspal)-1+NEQPAT+ii)
 650        CONTINUE   
         ENDIF!IGRAD
c
         CALL ZTIME (TBTMA1)
         TBTAB=TBTAB+(TBTMA1-TBTMA0)
21    CONTINUE
C
      CALL ZTIME (TSBLD1)
c 
c        
c         
      IF ((INT(OME1).GE.1).and.(ipatch.eq.1)) THEN
c      WRITE (*,*) 'ZEITEN:'
c      print *, 'Gesamt: ', (TSBLD1-TSBLD0)
c      print *, 'VORLAUF ', (TVORL1-TSBLD0) , '=', 
c     *           (TVORL1-TSBLD0)*100d0/ (TSBLD1-TSBLD0)
      print *,'    ', TRESLU , '=', (TRESLU*100D0)/(TSBLD1-TSBLD0),'%',
     * 'RESLU'
      print *,'    ', TBTAB  , '=', (TBTAB*100D0)/(TSBLD1-TSBLD0),'%',
     *   'BT-A-B'
      ENDIF
c
      IF (IGRAD.EQ.1) THEN
        CALL ZDISP(0,LSPAL1,'LSPAL1')
        CALL ZDISP(0,LSPAL2,'LSPAL2')
      ELSE
        CALL ZDISP(0,LSPAL,'LSPALT')
      ENDIF
C
            CALL ZDISP(0,LBT1,'LBT1  ')
            CALL ZDISP(0,LBT1I,'LBT1I ')
            CALL ZDISP(0,LBT1J,'LBT1J ')
            CALL ZDISP(0,LBT2,'LBT2  ')
            CALL ZDISP(0,LBT2I,'LBT2I ')
            CALL ZDISP(0,LBT2J,'LBT2J ')
      CALL ZTIME(TSBLD9)
      END
c
************************************************************************
      SUBROUTINE SKBLDN  (B11,B22,A,PIV,ALSK,NELPAT,NEQPAT,
     *                   NVELO,IGRAD)
************************************************************************
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNWORK=1,NNARR=299,NNLEV=9)     
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'
      INTEGER PIV
      DOUBLE PRECISION B11,B22
      DIMENSION  PIV(NVELO)            
      DIMENSION  B11(NEQPAT,NELPAT),B22(NEQPAT,NELPAT)
      DIMENSION  A(NVELO,NVELO), ALSK(NELPAT,NELPAT)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      SAVE
C
c=================================================================
C  		AUFBAU DER SCHURKOMPLEMENT MATRIX
c=================================================================
C
c-----------------------------------------------------------------
C     Bei der Benutzung der BLAS-routinen muss ich nicht
c     A^{-1}b_k rechnen, sondern kann direkt  A^{-1}B machen
c-----------------------------------------------------------------
C

      TRESLU=0
      TBTAB=0
      call ZTIME (TSBLD0)

c
c
c      DO 21 ISPALT=1,NELPAT 

         CALL ZTIME (TRELU0)
         IF (IGRAD.EQ.1) THEN
            CALL ZNEW (NEQPAT*NELPAT,1,LSPAL1,'LSPL1A')
            CALL ZNEW (NEQPAT*NELPAT,1,LSPAL2,'LSPL2A')
c
            CALL LCP1(B11(1,ISPALT),DWORK(L(LSPAL1)+(ispalt-1)*nelpat),
     *               NEQPAT)
            CALL LCP1(B22(1,ISPALT),DWORK(L(LSPAL2)+(ispalt-1)*nelpat),
     *               NEQPAT)
         ELSE
            CALL ZNEW (NVELO*NELPAT,1,LSPAL,'LSPALT')
c
            CALL LCP1(B11(1,ISPALT),DWORK(L(LSPAL)+(ispalt-1)*nelpat),
     *        NEQPAT)
            CALL LCP1(B22(1,ISPALT),DWORK(L(LSPAL)+
     *                             NEQPAT*NELPAT+(ispalt-1)*nelpat),
     *        NEQPAT)
         ENDIF
c
c21    CONTINUE
         IF (IGRAD.EQ.1) THEN
            CALL RESLU (NEQPAT  ,A,DWORK(L(LSPAL1)),PIV,LUMFM1,LUMFI1,
     *                  KEEP1,CNTL1,icntl1)
            CALL RESLU (NEQPAT  ,A,DWORK(L(LSPAL2)),PIV,LUMFM1,LUMFI1,
     *                  KEEP1,CNTL1,icntl1)
         ELSE
c
            CALL RESLU (2*NEQPAT,A,DWORK(L(LSPAL)) ,PIV,LUMFM1,LUMFI1,
     *                  KEEP1,CNTL1,icntl1)
         ENDIF
c
         CALL ZTIME (TRELU1)
         TRESLU=TRESLU+(TRELU1-TRELU0)
c-----------------------------------------------------------------
C	         In den Vektoren LSPAL* steht jetzt
c                A^{-1}B_{ISPALT}1/2,
C	         also die Loesung von A*X=B, wobei B die
C                ISPALT-te Spalte von B1/2 IST
c-----------------------------------------------------------------
c
c hier wieder dieselbe Frage:
c nota bene: FORTRAN speichert Matrizen SPALTENWEISE,
c     L(LMATRI) hingegen ZEILENWEIESE. Daher hier Alsk(j,i) statt Alsk(i,j)
         CALL ZTIME (TBTMA0)

         DO 210 ISPALT=1,NELPAT
         DO 220 IZEIL=1,NELPAT
         DO 230 II=1,NEQPAT
            IF (IGRAD.EQ.1) THEN
               ALSK(ISPALT,IZEIL)=ALSK(ISPALT,IZEIL)+
     *             B11(II,IZEIL)*DWORK(l(lspal1)-1+ii)
     *            +B22(II,IZEIL)*DWORK(l(lspal2)-1+ii)
            ELSE
               ALSK(ISPALT,IZEIL)=ALSK(ISPALT,IZEIL)+
     *             B11(II,IZEIL)*DWORK(l(lspal)-1       +ii)
     *            +B22(II,IZEIL)*DWORK(l(lspal)-1+NEQPAT+ii)
            ENDIF
230      CONTINUE
 220     CONTINUE
210      CONTINUE
c
         CALL ZTIME (TBTMA1)
         TBTAB=TBTAB+(TBTMA1-TBTMA0)
C
      CALL ZTIME (TSBLD1)
c            
      WRITE (*,*) 'ZEITEN:'
      print *, 'Gesamt: ', (TSBLD1-TSBLD0)
      print *, 'RESLU: ', TRESLU , '=', (TRESLU*100D0)/ (TSBLD1-TSBLD0)
     &  ,'%'
      print *, 'BT-A-B: ', TBTAB  , '=', (TBTAB*100D0)/ (TSBLD1-TSBLD0)
     ?  ,'%'

      IF (IGRAD.EQ.1) THEN
        CALL ZDISP(0,LSPAL1,'LSPAL1')
        CALL ZDISP(0,LSPAL2,'LSPAL2')
      ELSE
        CALL ZDISP(0,LSPAL,'LSPALT')
      ENDIF
C
      CALL ZTIME(TSBLD9)
      END
c
************************************************************************
      SUBROUTINE NEWSKB  (B11,B22,A,ALSK,NELPAT,NEQPAT)
************************************************************************
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNWORK=1,NNARR=299)     
      INCLUDE 'dwork.inc'
      DOUBLE PRECISION B11,B22
      DIMENSION  B11(NELPAT*NEQPAT),B22(NELPAT*NEQPAT)
      DIMENSION  A(NEQPAT*NEQPAT),
     *           ALSK(NELPAT*NELPAT)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      SAVE
C
C  		AUFBAU DER SCHURKOMPLEMENT MATRIX
C 

C
         CALL ZNEW (NEQPAT,1,LSPAL1,'LSPL1Y')
         CALL ZNEW (NEQPAT,1,LSPAL2,'LSPL2Y')

      DO 21 ISPALT=1,NELPAT !Ueber alle Spalten von B
c
         CALL ZCLEAR (LSPAL1,'LSPAL1')
         CALL ZCLEAR (LSPAL2,'LSPAL2')
         do 30 izl=1,neqpat   !INV(A) * Spalte ispalt von B1
            do 30 isp=1,neqpat
               dwork(l(lspal1)-1+izl)=
     *            dwork(l(lspal1)-1+izl) +
     *            a(NEQPAT*(izl-1)+isp) *
     *            b11(NEQPAT*(ISPALT-1)+isp)
c
               dwork(l(lspal2)-1+izl)=     
     *            dwork(l(lspal2)-1+izl)+
     *            a(NEQPAT*(izl-1)+isp) *
     *            b22(NEQPAT*(ISPALT-1)+isp)
 30            continue
c
c$$$         do 40 izl=1,neqpat   !INV(A) * Spalte ispalt von B2
c$$$            do 40 isp=1,neqpat
c$$$               dwork(l(lspal2)-1+izl)=     
c$$$     *            dwork(l(lspal2)-1+izl)+
c$$$     *            a(NEQPAT*(izl-1)+isp) *
c$$$     *            b22(NEQPAT*(ISPALT-1)+isp)
c$$$ 40            continue

C
       DO 210 IZEIL=1,NELPAT
           DO 230 II=1,NEQPAT
                ALSK((IZEIL-1)*NELPAT+ISPALT)=
     *          ALSK((IZEIL-1)*NELPAT+ISPALT)+
     *          B11((IZEIL-1)*NEQPAT+II)*
     *          DWORK(L(LSPAL1)-1+II)+
     *          B22((IZEIL-1)*NEQPAT+II)*
     *          DWORK(L(LSPAL2)-1+II)
230        CONTINUE
210     CONTINUE
c
21    CONTINUE
C

      CALL ZDISP(0,LSPAL1,'LSPAL1')
      CALL ZDISP(0,LSPAL2,'LSPAL2')
C
      END
c
************************************************************************
      SUBROUTINE TRANSQ (A,N)
************************************************************************
C
C	BILDET DIE TRANSPONIERTE EINER QUADRATISCHEN
C       MATRIX UND SPEICHERT SIE AUF DENSELBEN VEKTOR
C
       IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
       DOUBLE PRECISION A(N,N)
       DOUBLE PRECISION XHILF
      PARAMETER (NNWORK=1,NNARR=299)
      INCLUDE 'dwork.inc'
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      SAVE
C
      DO 1 I=1,N
         DO 2 J=I+1,N
           XHILF=A(I,J)
           A(I,J)=A(J,I)
           A(J,I)=XHILF
2        CONTINUE
1     CONTINUE
      END

C
C
C
****************************************************************
      SUBROUTINE ARCALC (RATIO,IEL,KVERT,DCORVG)
****************************************************************
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C	BERECHNET DEN ASPECT RATIO DES ELEMENTES IEL
C
      DOUBLE PRECISION NENNER, ZAEHLER
c      PARAMETER (NNBAS=21,NNVE=4)
      DIMENSION DCORVG(2,*),KVERT(4,*)
C
      IECK1=KVERT(1,IEL)    
      IECK2=KVERT(2,IEL)    
      IECK3=KVERT(3,IEL)    
      IECK4=KVERT(4,IEL)    
C
      XE1=DCORVG(1,IECK1)
      YE1=DCORVG(2,IECK1)
      XE2=DCORVG(1,IECK2)
      YE2=DCORVG(2,IECK2)
      XE3=DCORVG(1,IECK3)
      YE3=DCORVG(2,IECK3)
      XE4=DCORVG(1,IECK4)
      YE4=DCORVG(2,IECK4)
C
      XM1=(XE2+XE1)/2D0
      YM1=(YE2+YE1)/2D0
      XM2=(XE3+XE2)/2D0
      YM2=(YE3+YE2)/2D0
      XM3=(XE4+XE3)/2D0
      YM3=(YE4+YE3)/2D0
      XM4=(XE1+XE4)/2D0
      YM4=(YE1+YE4)/2D0
C
      ZAEHLER=(XM3-XM1)**2+(YM3-YM1)**2
      NENNER =(XM4-XM2)**2+(YM4-YM2)**2
C
      RATIO=SQRT(ZAEHLER/NENNER)
C
      END
C
****************************************************************
      SUBROUTINE ARCLC2 (ARMAX,HMIN,HMAX,NELEM,KVERT,DCORVG,ASPRAT)
****************************************************************
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C	BERECHNET DEN GROESSTEN ASPECT RATIO ALLER ELEMENTE
C       SOWIE HMIN und HMAX
C       Schreibt die AR in den Vektor ASPRAT
C
      DOUBLE PRECISION NENNER, ZAEHLER,ASPRAT
c      PARAMETER (NNBAS=21,NNVE=4)
      DIMENSION DCORVG(2,*),KVERT(4,*),ASPRAT(*)
C
      SAVE

      ARMAX=0D0
      HMAX=0D0
      HMIN=1D99
      DO 1 IEL=1,NELEM
C
      IECK1=KVERT(1,IEL)    
      IECK2=KVERT(2,IEL)    
      IECK3=KVERT(3,IEL)    
      IECK4=KVERT(4,IEL)    
C
      XE1=DCORVG(1,IECK1)
      YE1=DCORVG(2,IECK1)
      XE2=DCORVG(1,IECK2)
      YE2=DCORVG(2,IECK2)
      XE3=DCORVG(1,IECK3)
      YE3=DCORVG(2,IECK3)
      XE4=DCORVG(1,IECK4)
      YE4=DCORVG(2,IECK4)
C
      XM1=(XE2+XE1)/2D0
      YM1=(YE2+YE1)/2D0
      XM2=(XE3+XE2)/2D0
      YM2=(YE3+YE2)/2D0
      XM3=(XE4+XE3)/2D0
      YM3=(YE4+YE3)/2D0
      XM4=(XE1+XE4)/2D0
      YM4=(YE1+YE4)/2D0
C
      ZAEHLER=(XM3-XM1)**2+(YM3-YM1)**2
      NENNER =(XM4-XM2)**2+(YM4-YM2)**2
C
      RATIO=SQRT(ZAEHLER/NENNER)
C
c
      IF (RATIO.LT.1D0) RATIO=1D0/RATIO
      ASPRAT(IEL)=(RATIO)
      ARMAX = MAX(ARMAX,RATIO)
      HMIN=MIN(HMIN,ZAEHLER,NENNER)
      HMAX=MAX(HMAX,ZAEHLER,NENNER)
C
 1    CONTINUE
      END
c
****************************************************************
      SUBROUTINE GETIEL (IEL,INR1,INR2,RATIO,KMID,NVT)
****************************************************************
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KMID(4,*)
C
      IF (RATIO.LT.1D0) THEN
          INR1=KMID(1,IEL)-NVT
          INR2=KMID(3,IEL)-NVT
      ELSE
          INR1=KMID(2,IEL)-NVT
          INR2=KMID(4,IEL)-NVT
      ENDIF
C
      END
****************************************************************
      SUBROUTINE GETNEL (IELNR,INR,IELNEW,INRNEW,NVT,
     *                   KADJ,KMID,BRAND)
****************************************************************
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KMID(4,*),KADJ(4,*)
C
C	IN INR STEHT DIE (GLOBALE) NUMMER DES KNOTENS,
C   	DESSEN NACHBARELEMENT BESTIMMT WERDEN SOLL.
C	IN IELNR STEH DIE (ALTE) ELEMENTNUMMER
C	DIE ROUTINE LIEFERT IN IELNEW DIE NUMMER DES
C	NACHBARELEMENTES UND IN INRNEW DIE NUMMER DES KNOTENS,
C	IN DESSEN RICHTUNG GGF. WEITERGEMACHT WIRD.
C  
      DO 1 II=1,4
        IF ((KMID(II,IELNR)-NVT).EQ.INR) LOKNR=II
1     CONTINUE
C
      LOK2=MOD(LOKNR+1,4)+1
      IELNEW=KADJ(LOK2,IELNR)
      INRNEW=KMID(LOK2,IELNR)-NVT
      IF (IELNEW.EQ.0) BRAND=.TRUE.
99    END

C
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      SUBROUTINE LOEDOU (IDIM,NEQPAT,LFHG)
C
C     LOESCHT DOUBLETTEN IM VEKTOR LFHG
C     BESTIMMT NEQPAT, DIE ANZAHL DER UNBEKANNTEN IM PATCH
C
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299,NNBAS=21,NNVE=4,NNWORK=1)
C
C *** STANDARD DIMENSIONING FOR WORKSPACE CONCEPT
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
C *** STANDARD COMMON BLOCKS
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
C
C
      NEQPAT=0
C
      DO 123 ILAUF=IDIM,1,-1
       IF (KWORK(L(LFHG)+ILAUF).EQ.KWORK(L(LFHG)+ILAUF-1)) THEN
         NEQPAT=NEQPAT+1
C
         DO 124 JLAUF=ILAUF,IDIM-1
           KWORK(L(LFHG)+JLAUF)=KWORK(L(LFHG)+JLAUF+1)
124      CONTINUE
C
       ENDIF
123   CONTINUE
C
      NEQPAT=IDIM-NEQPAT
CCCCCC KAPPEN DER LETZEN, 'DOPPELTEN' EINTRAEGE

      CALL ZDISP(NEQPAT,LFHG,'FHG   ')
C
       END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE MATVEC (LLMAT,LVECKI,LSOL,LZEIL,LSPALT)
C
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      SAVE
C
C
      CALL ZLEN(LVECKI,LENG1)
      IF (LENG1.NE.LSPALT)
     *   WRITE (*,*) 'ERROR IN MATVEC! INCOMPATIBLE LENGHT LVEC'
C
      LSOL=0
      CALL ZNEW(LZEIL,1,LSOL,'LSOL  ')
C
      DO 1234 I=1,LZEIL
        DO 2345 J=1,LSPALT
          DWORK(L(LSOL)+I-1)=DWORK(L(LSOL)+I-1)+
     *        DWORK(L(LLMAT)-1+LSPALT*(I-1)+J)*DWORK(L(LVECKI)-1+J)
2345    CONTINUE
1234  CONTINUE
      END
C
C
C
************************************************************************
      SUBROUTINE SGRP1(NPATCH,NBLOCK,KPATCH,KLDPAT,KCHECK,KELC1,KELC2,
     *                 KVERT,KMID,KADJ,DCORVG,ASPRAT)
************************************************************************
C
C-----------------------------------------------------------------------
C PURPOSE:   FASST GEWISSE ELEMENTE ZU GRUPPEN (PATCHES) ZUSAMMEN
C
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'
C-----------------------------------------------------------------------
C
c      PARAMETER (NNVE=4)
      DIMENSION KPATCH(*),KLDPAT(*),KCHECK(*),KELC1(*),KELC2(*),
     *          KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*),DCORVG(2,*),
     *          ASPRAT(*)
C
C
C
      NPATCH=0
      LAUF=0
      LAUF2=1
      KLDPAT(1)=1
C
C
      DO 1 IEL=1,NEL
C
      IF (KCHECK(IEL).NE.0) GOTO 1
C



      CALL ARCALC(RATIO,IEL,KVERT,DCORVG)
c$$$      RATIO=ASPRAT(IEL)
      RATIOH=RATIO
      IF (RATIO.LT.1D0) RATIO=1D0/RATIO
      IF (RATIO.LT.PSBARB) GOTO 1
C
      NPATCH=NPATCH+1
C
      LAUF=LAUF+1
      LAUFH=LAUF-1
      KCHECK(IEL)=1
      KPATCH(LAUF)=IEL
C
      CALL GETIEL (IEL,INR1,INR2,RATIOH,KMID,NVT)
C
C
      BRAND=.FALSE.
      INR=INR1
      IELNR=IEL
C
2     CALL GETNEL (IELNR,INR,IELNEW,INRNEW,NVT,KADJ,KMID,BRAND)
      IF ((KCHECK(IELNEW).EQ.1).OR.BRAND) GOTO 3
C
      CALL ARCALC (RATIO,IELNEW,KVERT,DCORVG) 
c$$$      RATIO=ASPRAT(IELNEW)
      IF (RATIO.LT.1D0) RATIO=1D0/RATIO
      IF ((LAUF-LAUFH).GE.NBLOCK) GOTO 3
C
      KCHECK(IELNEW)=1
      LAUF=LAUF+1
      KPATCH(LAUF)=IELNEW
      IELNR=IELNEW
      INR=INRNEW
      IF (RATIO.LT.PSBARB) GOTO 3
      GOTO 2
C
3     CONTINUE
C
      BRAND=.FALSE.
      INR=INR2
      IELNR=IEL
C
20    CALL GETNEL (IELNR,INR,IELNEW,INRNEW,NVT,KADJ,KMID,BRAND)
      IF ((KCHECK(IELNEW).EQ.1).OR.BRAND) GOTO 30
C
      CALL ARCALC (RATIO,IELNEW,KVERT,DCORVG) 
c$$$      RATIO=ASPRAT(IELNEW)

      IF (RATIO.LT.1D0) RATIO=1D0/RATIO
      IF ((LAUF-LAUFH).GE.NBLOCK) GOTO 30
C
      KCHECK(IELNEW)=1
      LAUF=LAUF+1
      KPATCH(LAUF)=IELNEW
      IELNR=IELNEW
      INR=INRNEW
      IF (RATIO.LT.PSBARB) GOTO 30
      GOTO 20
C
30    CONTINUE
C
      LAUF2=LAUF2+1         
      KLDPAT(LAUF2)=LAUF+1
1     CONTINUE
C
      DO 200 IEL=1,NEL
      IF (KCHECK(IEL).NE.0) GOTO 200
C
      NPATCH=NPATCH+1

      KCHECK(IEL)=1
      LAUF=LAUF+1
      LAUF2=LAUF2+1                      
      KPATCH(LAUF)=IEL
      KLDPAT(LAUF2)=LAUF+1
200   CONTINUE
C
C
C
99999 END
C



c====================================================================
c  Subroutinen fuer die Blockung nach GG-Elementen
c  i.e. alle Tochterelemente eines GG-Elementes bilden ein Patch
c====================================================================
************************************************************************
      SUBROUTINE SGRP3(NPATCH,NBLOCK,KPATCH,KLDPAT,KCHECK,KELC1,KELC2,
     *                 KVERT,KMID,KADJ,DCORVG)
************************************************************************
C
C-----------------------------------------------------------------------
C PURPOSE:   
C
C-----------------------------------------------------------------------
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'block.inc'
C-----------------------------------------------------------------------24
C
      INTEGER NELC
c      PARAMETER (NNVE=4)
      DIMENSION KPATCH(*),KLDPAT(*),KCHECK(*),KELC1(*),KELC2(*),
     *          KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*),DCORVG(2,*)
C
C
C
      IF (NEL.GE.NELMAC) THEN
        NPATCH=NELMAC
        NELC=NEL/NELMAC
        IPATCH=1
        KLDPAT(IPATCH)=1
C
        DO 100 IEL=1,NEL
        KPATCH(IEL)=KELC2(IEL)
        IF (MOD(IEL,NELC).EQ.0) THEN
         IPATCH=IPATCH+1
         KLDPAT(IPATCH)=(IPATCH-1)*NELC+1
        ENDIF
100   CONTINUE
      ELSE !Falls ich spaeter mit dem Blocken anfange
         NPATCH=NEL
         NELC=NELC
         KLDPAT(1)=1
         do 200 iel=1,nel
            KPATCH(IEL)=iel
            KLDPAT(IEL+1)=iel+1
 200     continue
      ENDIF
C
C
C
99999 END
************************************************************************
      SUBROUTINE MAC2EL (KELC1,KELC2,KELC1H,KELC2H,KADJ,NEL,NELMAC,
     *                   INDLEV)
************************************************************************
C
C-----------------------------------------------------------------------
C PURPOSE:   
C
C-----------------------------------------------------------------------
C
      PARAMETER (NNVE=4)
      DIMENSION KELC1(*),KELC2(*),KELC1H(*),KELC2H(*),KADJ(NNVE,*)
C
C
C
      IF (INDLEV.EQ.1) THEN
C      
       DO 100 IEL=1,NEL
       KELC1(IEL)=IEL
       KELC2(IEL)=IEL
100    CONTINUE
C
      ELSE
C
       DO 200 IELC=1,NELMAC*4**(INDLEV-2)
       IEL1=IELC
       IEL2=KADJ(2,IEL1)
       IEL3=KADJ(2,IEL2)
       IEL4=KADJ(2,IEL3)
       KELC1(IEL1)=KELC1H(IELC)
       KELC1(IEL2)=KELC1H(IELC)
       KELC1(IEL3)=KELC1H(IELC)
       KELC1(IEL4)=KELC1H(IELC)
200    CONTINUE
C

       INDC2=0
       DO 210 IELC=1,NELMAC
       DO 220 IEL =1,NEL
       IF (KELC1(IEL).EQ.IELC) THEN
        INDC2=INDC2+1
        KELC2(INDC2)=IEL
       ENDIF
220    CONTINUE
210    CONTINUE
C
      ENDIF
C
C
C
C
99999 END
c--------------------------------------------------
      SUBROUTINE COLMAT(A1,LFHG,NAT,NVELO,NEQPAT,NELPAT)
c--------------------------------------------------
c
c     Sammelt aus der 'grossen' MAtrix A1 die
c     Eintraege, die zu der Patchmatrix gehoeren
c
c--------------------------------------------------
c
c      DOUBLE PRECISION A1
c
      INCLUDE 'common.inc'
      INCLUDE 'dwork.inc'
      INCLUDE 'bouss.inc'
      INCLUDE 'block.inc'
      DIMENSION A1(NVELO,NVELO)
      SAVE

c-----------------------------------------------------------------
C      Finde Element A(irow,icol) in der grossen Matrix
c      und setze es auf A1(ii,jj). Wird nicht mit der
c      Gradientenform gerechnet, kommen auch noch die 
c      entsprechenden Eintraege in die drei anderen
c      Bloecke der MAtrix
c-----------------------------------------------------------------
      CALL VKADD (AVK1,AVK2)

      DO 100 II=1,NEQPAT
        IROW=KWORK(L(LFHG)+II-1)
        ILD =KWORK(KLDA-1+IROW)
        ILD2=KWORK(KLDA  +IROW)
C
        DO 200 JJ=1,NEQPAT
C
          JCOL=KWORK(L(LFHG)+JJ-1)
C                                        
          DO 300 LAUF=ILD,(ILD2-1)
            IF (JCOL.EQ.KWORK(KCOLA-1+LAUF)) THEN
              A1(jj,ii)           =AVK1*dble(vwork(KA1     -1+lauf))+
     *                             AVK2*dble(vwork(KA1+4*NA-1+lauf))
c nota bene: FORTRAN speichert Matrizen SPALTENWEISE,
c     L(LMATRI) hingegen ZEILENWEIESE. Daher hier A1(j,i) statt A1(i,j)
             IF (IGRAD.NE.1) THEN
              A1(jj+NEQPAT,ii+NEQPAT)=AVK1*dble(vwork(KA1+3*NA-1+lauf))+
     *                                AVK2*dble(vwork(KA1+4*NA-1+lauf))
              A1(jj+NEQPAT,ii)       =AVK1*dble(vwork(KA1+  NA-1+lauf))
              A1(jj,ii+NEQPAT)       =AVK1*dble(vwork(KA1+2*NA-1+lauf))
             ENDIF
              GOTO 200
            ENDIF
 300      CONTINUE           
 200    CONTINUE
 100  CONTINUE
c
c-----------------------------------------------------------------
c     Behandlung des Dirichletrandes!
c-----------------------------------------------------------------
        DO 400 I=1,NEQPAT
          IROW=KWORK(L(LFHG)+I-1)
          IF (KWORK(L(LNPR)-1+IROW+NVT).NE.0) THEN
            DO 500 J=1,NEQPAT
              A1(j,i)              =0d0
c              A1(j,i)              =0d0
              IF (IGRAD.NE.1) THEN
              A1(j+NEQPAT,i)       =0d0
              A1(j,i+NEQPAT)       =0d0
c              A1(j,i+NEQPAT)       =0d0
c              A1(j+NEQPAT,i)       =0d0
              A1(j+NEQPAT,i+NEQPAT)=0d0
              ENDIF
 500        CONTINUE
c
            A1(i,i)              =1D99
            IF (IGRAD.NE.1) 
     *      A1(i+NEQPAT,i+NEQPAT)=1d99
          ENDIF
 400    CONTINUE
c
c
c$$$           write (*,*) nu, 'Matraze2'
c$$$          DO 7712 i=1,nvelo
c$$$             do 7712 j=1,nvelo
c$$$ 7712           print * , i,j,a1(j,i)    


99999    END
c
c
************************************************************************
      SUBROUTINE COLMAB(BB1,BB2,B1,B2,NELPAT,NEQPAT,
     *                  LFHG,LPATCH,LLDPAT,IPATCH,KLDB,KCOLB)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      DOUBLE PRECISION BB1,BB2
      REAL    B1,B2         
      PARAMETER (NNWORK=1,NNARR=299)
      DIMENSION BB1(NEQPAT,NELPAT), BB2(NEQPAT,NELPAT)
      DIMENSION B1(*),B2(*)
      INCLUDE 'dwork.inc'
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      SAVE
c-----------------------------------------------------------------
C      Finde Element B1(irow,icol) [sowie B2(irow,icol)]
C      in der grossen Matrix und setze es auf BB1(ii,jj).
c      [respektive BB2(irow,icol)] (Merke: Aufbau von B, *nicht* B^T
c-----------------------------------------------------------------
C
         DO 1444 II=1,NEQPAT
          IROW=KWORK(L(LFHG)+II-1)
           ILD=KWORK(KLDB-1+IROW)
           ILD2=KWORK(KLDB+IROW)
C
C 
          ICOUNT=0
          DO 1555 JJ=1,NELPAT
           JCOL=KWORK(L(LPATCH)-1+KWORK(L(LLDPAT)-1+IPATCH)-1+JJ)
C             Finde Element B(IROW,JCOL) der grossen Matrix
C
            DO 1666 LAUF=ILD,(ILD2-1)
             IF (JCOL.EQ.KWORK(KCOLB-1+LAUF)) THEN
               ICOUNT=ICOUNT+1
               ILAUF=LAUF
               BB1(II,JJ) =DBLE(B1(LAUF))
               BB2(II,JJ) =DBLE(B2(LAUF))
             ENDIF
1666        CONTINUE
1555       CONTINUE
1444      CONTINUE
c
      END

************************************************************************
      SUBROUTINE  YAXA (DX,DAX,NEQ,A1,A2)  
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
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
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
      INCLUDE 'bouss.inc'
      INCLUDE 'jump.inc'
C-----------------------------------------------------------------------
C
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      EXTERNAL MATML
      SAVE
C=======================================================================
C     Getting all parameters for MATML
C=======================================================================
C *** addresses for the current level ILEV
C
      ISETLV=2
      CALL  SETLEV (ISETLV)
      I2=1+NU
      IP=I2+NU
C
      KMBD =L(KLMBD(ILEV))
      NMBD= KNMBD(ILEV)
      KNY=L(KLNY(NLEV))
C
C=======================================================================
c
       IF(IJUMP .GE. 1)THEN 
         CALL JUMPL (DWORK(KU1),DWORK(KU2),DWORK(KU1),DWORK(KU2),
     *                1D0,0D0,DX(1),DX(I2),
     *                DAX(1),DAX(I2),
     *                VWORK(KA1),NA,KWORK(KCOLA),
     *                KWORK(KLDA),KWORK(L(LVERT)),KWORK(L(LMID)),
     *                KWORK(L(LADJ)),KWORK(L(LMEL)),
     *                DWORK(L(LCORVG)),E031,COEFFN,IDEFUP,1D0,
     *                DWORK(KNY))
c
      ENDIF
c
      IF (IPRECO.NE.1) then!Voller Deftensor als VK
        CALL MATML(DAX(1),DAX(I2),DAX(IP),DX(1),DX(I2),DX(IP),A1,A2,
     *            VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *            VWORK(KB1),VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *            NU,NP,KWORK(KMBD),NMBD,INEUM)
C
      ELSE
        CALL MATML (DAX(1),DAX(I2),DAX(IP),DX(1),DX(I2),DX(IP),A1,A2,
     *            VWORK(KA1),KWORK(KCOLA),KWORK(KLDA),
     *            VWORK(KB1),VWORK(KB2),KWORK(KCOLB),KWORK(KLDB),
     *            NU,NP,KWORK(KMBD),NMBD,INEUM)
      ENDIF
C
      END
C
C
C
      SUBROUTINE VKADD(AVK1,AVK2)
c
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNLEV=9)
      DOUBLE PRECISION AVK1,AVK2
      INCLUDE 'bouss.inc'
c
c      IF (IPRECO.EQ.0) then!Voller Deftensor als VK
        AVK1=1d0
         AVK2=0d0
c      ELSE !Nur Grad-Grad als Vorkonditionierer
c         AVK1=0d0
c         AVK2=1D0
c      endif
      END
