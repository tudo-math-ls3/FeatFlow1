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
* XABM7                                                                *
*                                                                      *
* Purpose  Allocation of the blocks in vector DA , KOFF and COECON     *
*          Call of ABM7                                                *
*                                                                      *
* Subroutines/functions called  ABM7, ZNEW, ZDISP, ZLEN, ZTYPE, ZCLEAR *
*                                                                      *
* Version from  12/08/94                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LA       I*4    Vector of numbers of the arrays                      *
*                 New arrays are allocated on DWORK for LA(IBLOC)=0    *
* LCOL     I*4    Numbers of the pointer arrays KCOL and KLD           *
* LLD      I*4    calculated by AP7                                    *
* ICLEAR   I*4    ICLEAR=1  all matrices are cleared                   *
* ARR      C*6    Names of blocks (for error messages only)            *
* BSNGL    L*4    NBLOC logical values                                 *
*                 .TRUE. means that corresponding array is converted   *
*                 to single precision after completion                 *
* For the description of the remaining parameters see AB07             *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* LA       I*4    Vector of numbers of arrays on DA                    *
* IER      I*4    Error indicator                                      *
*                 -114 Type of at least one array is not double prec.  *
*                 -115 Length of at least one array is < NA            *
*                                                                      *
************************************************************************
*                                                                      *
* ABM7                                                                 *
*                                                                      *
* Purpose  Calculation of NBLOC matrices corresponding to              *
*          bilinear forms  a(v,u)                                      *
*          Nonparametric quadrilateral elements                        *
*          Same trial and test functions                               *
*          Storage technique 7/8                                       *
*          Each matrix makes use of the same pointer vectors           *
*          Double precision version                                    *
*                                                                      *
* Subroutines/functions called NDFL, NDFGL, CB2Q                       *
*                                                                      *
* Version from  12/08/94                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* KCOL     I*4    Pointer vectors for each matrix on DA corresponding  *
* KLD      I*4    to storage technique 7/8 (calculated by AP7)         *
* NA       I*4    Number of elements per matrix                        *
* NEQ      I*4    Number of equations per matrix                       *
* NBLOC    I*4    Number of matrices stored on DA                      *
* KOFF     I*4    Matrix IBLOC starts at position KOFF(IBLOC)+1 on DA  *
* KVERT    I*4    Arrays describing the triangulation                  *
* KMID     I*4                                                         *
* DCORVG   R*8                                                         *
* ELE      SUBR   EXTERNAL SUBROUTINE - values of basis functions      *
* COEFF    R*8    EXTERNAL FUNCTION - coefficients of the bilinear form*
* BCON     LOG    BCON(IBLOC)=.TRUE. means constant coefficients in    *
*                 matrix IBLOC                                         *
* KAB      I*4    Pairs of multiindices occuring in the bilinear forms *
*                 specified separately for each matrix                 *
* KABN     I*4    Number of additive terms in each bilinear form       *
* ICUB     I*4    Number of cubature formula in CB2Q                   *
* ISYMM    I*4    0  storage technique 7                               *
*                 1  storage technique 8                               *
*                 2  storage technique 7 but for symmetric matrices    *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DA       R*8    Calculated matrices                                  *
* IER      I*4    Error indicator                                      *
*                 -116 Wrong value in array KAB                        *
*                                                                      *
************************************************************************
C
      SUBROUTINE XABM7(LA,LCOL,LLD,NA,NEQ,NBLOC,ICLEAR,ELE,
     *                 COEFF,BCON,KAB,KABN,ICUB,ISYMM,ARR,BSNGL)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*12
C
      PARAMETER (NNARR=299,NNAB=21,NNDER=6)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION LA(*),KAB(2,NNAB,*),KABN(*),BCON(*),ARR(*),BSNGL(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL COEFF,ELE
      SAVE /TRIAA/,/OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='XABM7'
      IF (ICHECK.GE.997) CALL OTRC('XABM7 ','12/08/94')
      IER=0
C
      DO 1 IBLOC=1,NBLOC
      IF (LA(IBLOC).EQ.0) THEN
       CALL ZNEW(NA,1,LA(IBLOC),ARR(IBLOC))
       IF (IER.NE.0) GOTO 99999
      ELSE
C ***  Check input parameter
       CALL ZTYPE(LA(IBLOC),ITYPE)
       IF (ITYPE.NE.1) THEN
        WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
        CALL WERR(-114,'XABM7 ')
        GOTO 99999
       ENDIF
       CALL ZLEN(LA(IBLOC),ILEN)
       IF (ILEN.LT.NA) THEN
        WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
        CALL WERR(-115,'XABM7 ')
        GOTO 99999
       ENDIF
       IF (ICLEAR.EQ.1) CALL ZCLEAR(LA(IBLOC),ARR(IBLOC))
      ENDIF
1     CONTINUE
C
      CALL ZNEW(NBLOC,-3,LOFF,'KOFF  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NBLOC*NNDER**2,1,LOECON,'COECON')
      IF (IER.NE.0) GOTO 99999
C
      DO 2 IBLOC=1,NBLOC
      KWORK(L(LOFF)+IBLOC-1)=L(LA(IBLOC))-1
2     CONTINUE
C
      CALL ABM7(DWORK(1),KWORK(L(LCOL)),KWORK(L(LLD)),NA,NEQ,
     *          NBLOC,KWORK(L(LOFF)),KWORK(L(LVERT)),KWORK(L(LMID)),
     *          DWORK(L(LCORVG)),ELE,COEFF,
     *          BCON,DWORK(L(LOECON)),KAB,KABN,ICUB,ISYMM)
      IF (IER.NE.0) GOTO 99999
C
      CALL ZDISP(0,LOECON,'COECON')
      CALL ZDISP(0,LOFF,'KOFF  ')
C
      DO 3 IBLOC=NBLOC,1,-1
      IF (BSNGL(IBLOC)) THEN
       CALL ZCTYPE(2,LA(IBLOC),ARR(IBLOC))
       IF (IER.NE.0) GOTO 99999
      ENDIF
3     CONTINUE
C
99999 END
C
C
C
      SUBROUTINE ABM7(DA,KCOL,KLD,NA,NEQ,NBLOC,KOFF,KVERT,KMID,
     *                DCORVG,ELE,COEFF,BCON,COECON,
     *                KAB,KABN,ICUB,ISYMM)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21,NNCOF=6)
      DIMENSION KCOL(*),KVERT(NNVE,*),KMID(NNVE,*),KOFF(*)
      DIMENSION KLD(*),KDFG(NNBAS),KDFL(NNBAS),DA(*),DCORVG(2,*)
      DIMENSION BCON(*),KAB(2,NNAB,*),KABN(*),COECON(NNDER,NNDER,*)
      DIMENSION KENTRY(NNBAS,NNBAS)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      COMMON /COFBAS/ COB(NNBAS,NNCOF)
      EXTERNAL ELE,COEFF
      SAVE
C
      SUB='ABM7'
      IF (ICHECK.GE.997) CALL OTRC('ABM7  ','12/08/94')
C
C *** Preparation - evaluation of parameters
      IER=0
      BSYMM=ISYMM.GE.1
C *** Which derivatives of basis functions are needed?
      DO 1 I = 1,NNDER
1     BDER(I)=.FALSE.
      DO 2 IBLOC=1,NBLOC
      DO 3 I=1,KABN(IBLOC)
      DO 3 J=1,2
      I1=KAB(J,I,IBLOC)
      IF (I1.LE.0.OR.I1.GT.NNDER) THEN
       WRITE (CPARAM,'(I15)') IBLOC
       CALL WERR(-116,'AB07  ')
       GOTO 99999
      ENDIF
3     BDER(I1)=.TRUE.
2     CONTINUE
C *** Dummy call of ELE sets number of element
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      BFIRST=.TRUE.
C *** Dummy call of COEFF for nonlinear problems
C *** COEFF must set BDER(IDER)=.TRUE. if derivative IDER is needed
      AUX=COEFF(0D0,0D0,-1,-1,0,BFIRST)
      BCON0=.TRUE.
C
      DO 4 IBLOC=1,NBLOC
      IF (BCON(IBLOC)) THEN
       DO 5 I=1,KABN(IBLOC)
       IA=KAB(1,I,IBLOC)
       IB=KAB(2,I,IBLOC)
5      COECON(IA,IB,IBLOC)=COEFF(0D0,0D0,IA,IB,IBLOC,BFIRST)
      ELSE
       BCON0=.FALSE.
      ENDIF
4     CONTINUE
C********************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
C********************************************************************
C
C *** Loop over all elements
      JDOFE1=1
      DO 100 IEL=1,NEL
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Determine entry positions in matrix
      DO 110 JDOFE=1,IDFL
      ILD=KLD(KDFG(JDOFE))
      KENTRY(JDOFE,JDOFE)=ILD
      IF (BSYMM) JDOFE1=JDOFE
      JCOL0=ILD
      DO 111 IDOFE=JDOFE1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 111
      IDFG=KDFG(IDOFE)
      DO 112 JCOL=JCOL0,NA
      IF (KCOL(JCOL).EQ.IDFG) GOTO 113
112   CONTINUE
113   JCOL0=JCOL+1
      KENTRY(JDOFE,IDOFE)=JCOL
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
C *** Dummy call - ELE may save arithmetic operations
      CALL ELE(0D0,0D0,-2)
      IF (IER.LT.0) GOTO 99999
C
C *** Loop over all cubature points
      DO 200 ICUBP = 1,NCUBP
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
C *** Cubature points on the actual element + weights
      XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *  +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
      YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1+0.5D0*
     *         (DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2
      OM=DOMEGA(ICUBP)*DETJ
C
      CALL ELE(XX,YY,-3)
      IF (IER.LT.0) GOTO 99999
C
C *** Summing up over all pairs of multiindices
      BFIRST=.TRUE.
      DO 300 IBLOC=1,NBLOC
C
      DO 301 IALBET=1,KABN(IBLOC)
      IA=KAB(1,IALBET,IBLOC)
      IB=KAB(2,IALBET,IBLOC)
      IF (.NOT.BCON(IBLOC)) THEN
       AUX=COEFF(XX,YY,IA,IB,IBLOC,BFIRST)*OM
      ELSE
       AUX=COECON(IA,IB,IBLOC)*OM
      ENDIF
C
      DO 310 JDOFE=1,IDFL
      DB=DBAS(KDFL(JDOFE),IA)
      IF (BSYMM) JDOFE1=JDOFE
      DO 310 IDOFE=JDOFE1,IDFL
      JCOLB=KENTRY(JDOFE,IDOFE)+KOFF(IBLOC)
      DA(JCOLB)=DA(JCOLB)+DB*DBAS(KDFL(IDOFE),IB)*AUX
310   CONTINUE
301   CONTINUE
C
      BFIRST=.FALSE.
300   CONTINUE
C
200   CONTINUE
100   CONTINUE
C
      IF (ISYMM.GT.1) THEN
       DO 400 IBLOC=1,NBLOC
C
       DO 401 IROW=1,NEQ
       DO 402 ICOL=KLD(IROW)+1,KLD(IROW+1)-1
       J1=KCOL(ICOL)
       IF (J1.GE.IROW) GOTO 401
       DO 403 ICOL1=KLD(J1+1)-1,KLD(J1),-1
       IF (KCOL(ICOL1).EQ.IROW) GOTO 404
403    CONTINUE
       GOTO 402
404    ICOLB=ICOL+KOFF(IBLOC)
       ICOL1B=ICOL1+KOFF(IBLOC)
C ***  Elements ICOL and ICOL1 in block IBLOC
       DA(ICOLB)=DA(ICOL1B)
402    CONTINUE
401    CONTINUE
C
400    CONTINUE
      ENDIF
C
99999 END
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
* XMABM7                                                               *
*                                                                      *
* Purpose  Calculate matrices corresponding to bilinear forms          *
*          (multigrid version)                                         *
*          Successive call of XABM7                                    *
*                                                                      *
* Subroutines/functions called  XABM7                                  *
*                                                                      *
* Version from  12/08/94                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEV     I*4    Number of levels                                     *
*                 for which matrices are generated                     *
* KLA      I*4    NBLOC*NLEV numbers of arrays                         *
* KLCOL    I*4    NLEV numbers of pointer vectors                      *
* KLLD     I*4                                                         *
* NBLOC    I*4                                                         *
* ICLEAR   I*4                                                         *
* ELE      SUBR                                                        *
* COEFF    SUBR                                                        *
* BCON     L*4                                                         *
* KAB      I*4    See parameters of XAB07 and AB07                     *
* KABN     I*4                                                         *
* ICUB     I*4                                                         *
* ISYMM    I*4                                                         *
* CARR     C*6    Names of matrices for all levels                     *
* BSNGL    L*4    NBLOC*NLEV logical values                            *
*                 .TRUE. meand that corresponding array is converted   *
*                 to single precision after completion                 *
* Meshes on COMMON /MGTRD/ and /MGTRA/                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KLA      I*4    NBLOC*NLEV numbers of arrays on DA                   *
*                                                                      *
************************************************************************
C
      SUBROUTINE XMABM7(KLA,KLCOL,KLLD,KNA,KNEQ,NBLOC,ICLEAR,ELE,
     *                  COEFF,BCON,KAB,KABN,ICUB,ISYMM,CARR,BSNGL)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CARR*12
C
      PARAMETER (NNARR=299,NNAB=21,NNLEV=9)
      DIMENSION KLA(NBLOC,*),KLCOL(*),KLLD(*)
      DIMENSION KNA(*),KNEQ(*),KAB(2,NNAB,*),KABN(*),BCON(*)
      DIMENSION CARR(NBLOC,*),BSNGL(NBLOC,*)
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
      EXTERNAL ELE,COEFF
      SAVE /TRIAA/,/TRIAD/,/MGTRD/,/MGTRA/,/MGPAR/,/ERRCTL/,/CHAR/
C
      SUB='XMABM7'
      IF (ICHECK.GE.997) CALL OTRC('XMABM7','12/08/94')
      IER=0
C
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
      CALL XABM7(KLA(1,ILEV),KLCOL(ILEV),KLLD(ILEV),KNA(ILEV),
     *           KNEQ(ILEV),NBLOC,ICLEAR,ELE,COEFF,BCON,
     *           KAB,KABN,ICUB,ISYMM,CARR(1,ILEV),BSNGL(1,ILEV))
      IF (IER.NE.0) GOTO 99999
C
10    CONTINUE
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
* XMAB09                                                               *
*                                                                      *
* Purpose  Calculate matrices corresponding to bilinear forms          *
*          (multigrid version)                                         *
*          Successive call of XAB09                                    *
*                                                                      *
* Subroutines/functions called  XAB09                                  *
*                                                                      *
* Version from  07/21/92                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEV     I*4    Number of levels                                     *
*                 for which matrices are generated                     *
* KLA      I*4    NBLOC*NLEV numbers of arrays                         *
* KLCOL    I*4    NLEV numbers of pointer vectors                      *
* KLLD     I*4                                                         *
* NBLOC    I*4                                                         *
* ICLEAR   I*4                                                         *
* ELE1     SUBR                                                        *
* ELE2     SUBR                                                        *
* ELE3     SUBR                                                        *
* COEFF    SUBR                                                        *
* BCON     L*4                                                         *
* KAB      I*4    See parameters of XAB07 and AB07                     *
* KABN     I*4                                                         *
* ICUB     I*4                                                         *
* CARR     C*6    Names of matrices for all levels                     *
* BSNGL    L*4    NBLOC*NLEV logical values                            *
*                 .TRUE. meand that corresponding array is converted   *
*                 to single precision after completion                 *
* Meshes on COMMON /MGTRD/ and /MGTRA/                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KLA      I*4    NBLOC*NLEV numbers of arrays on DA                   *
*                                                                      *
************************************************************************
C
      SUBROUTINE XMAB09(KLA,KLCOL,KLLD,KNA,NBLOC,ICLEAR,ELE1,ELE2,ELE3,
     *                  COEFF,BCON,KAB,KABN,ICUB,CARR,BSNGL)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CARR*12
C
      PARAMETER (NNARR=299,NNAB=21,NNLEV=9)
      DIMENSION KLA(NBLOC,*),KLCOL(*),KLLD(*)
      DIMENSION KNA(*),KAB(2,NNAB,*),KABN(*),BCON(*)
      DIMENSION CARR(NBLOC,*),BSNGL(NBLOC,*)
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
      EXTERNAL ELE1,ELE2,ELE3,COEFF
      SAVE /TRIAA/,/TRIAD/,/MGTRD/,/MGTRA/,/MGPAR/,/ERRCTL/,/CHAR/
C
      SUB='XMAB09'
      IF (ICHECK.GE.997) CALL OTRC('XMAB09','07/21/92')
      IER=0
C
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
      CALL XAB09(KLA(1,ILEV),KLCOL(ILEV),KLLD(ILEV),KNA(ILEV),
     *           NBLOC,ICLEAR,ELE1,ELE2,ELE3,COEFF,BCON,
     *           KAB,KABN,ICUB,CARR(1,ILEV),BSNGL(1,ILEV))
      IF (IER.NE.0) GOTO 99999
C
10    CONTINUE
C     
99999 END
c
c
c
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
* XMAP9                                                                *
*                                                                      *
* Purpose  Calculate pointer vectors                                   *
*          (multigrid version)                                         *
*          Successive call of XAP9                                     *
*                                                                      *
* Subroutines/functions called  XAP9                                   *
*                                                                      *
* Version from  07/21/92                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEV     I*4    Number of levels                                     *
* ELE      SUBR                                                        *
* ISYMM    I*4                                                         *
* Meshes on COMMON /MGTRD/ and /MGTRA/                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KNA      I*4    Number of entries for each level                     *
* KNEQ     I*4    Number of unknowns for each level                    *
* IER      I*4    Error indicator                                      *
*                 Set by ZNEW                                          *
*                                                                      *
************************************************************************
      SUBROUTINE XMAP9(KLCOL,KLLD,KNA,KNEQ,ELE1,ELE2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNARR=299,NNLEV=9)
      DIMENSION KLCOL(*),KLLD(*),KNA(*),KNEQ(*)
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
      SAVE /TRIAA/,/TRIAD/,/MGTRD/,/MGTRA/,/MGPAR/,/ERRCTL/,/CHAR/
C
      SUB='XMAP9 '
      IF (ICHECK.GE.997) CALL OTRC('XMAP9 ','07/21/92')
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
      CALL XAP9(KLCOL(ILEV),KLLD(ILEV),KNA(ILEV),KNEQ(ILEV),ELE1,ELE2)
      IF (IER.NE.0) GOTO 99999
C     
10    CONTINUE
C     
99999 END
c
c
c
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
* XMSB2                                                                *
*                                                                      *
* Purpose  Generate sequence of meshes for multigrid applications      *
*          by successive calls of XSB0X  (two-level ordering)          *
*                                                                      *
* Subroutines/functions called  XSB0X, WERR, ZCPY                      *
*                                                                      *
* Version from  05/11/91                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEV     I*4    desired number of l                                  *
* IMID     I*4                                                         *
* IADJ     I*4                                                         *
* IVEL     I*4                                                         *
* IDISP    I*4                                                         *
* IBDP     I*4     See parameter list of XSB0X                         *
* S2DI     SUBR                                                        *
* S2DB     SUBR                                                        *
* PARX     FNCT                                                        *
* PARY     FNCT                                                        *
* TMAX     FNCT                                                        *
* Coarse grid on /TRIAD and /TRIAA/                                    *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DCORVG   R*8   Cartesian coordinates of vertices                     *                                                                      *                                                                      *
************************************************************************
C
      SUBROUTINE XMSB2(IMID,IADJ,IVEL,IDISP,IBDP,
     *                 S2DI,S2DB,PARX,PARY,TMAX)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
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
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      EXTERNAL S2DI,S2DB,PARX,PARY,TMAX
      SAVE 
C
      SUB='XMSB2 '
      IF (ICHECK.GE.997) CALL OTRC('XMSB2 ','11/05/91')
C
      NLEV=NLMAX
C
      IF (NLEV.GT.NNLEV) THEN
       CALL WERR(-180,'XMSB2 ')
       GOTO 99999
      ENDIF
C
C
      IADJ0=1
      IBDP0=2
C
C
      DO 10 ILEV=1,NLEV
      IF (ILEV.EQ.NLEV) THEN
       IADJ0=IADJ
       IBDP0=IBDP
      ENDIF
      IF (ILEV.EQ.1) THEN
       CALL XSB0X(0,IMID,IADJ0,IVEL,IDISP,IBDP0,
     *            S2DI,S2DB,PARX,PARY,TMAX)
      ELSE
       CALL XSB0X(1,IMID,IADJ0,IVEL,IDISP,IBDP0,
     *            S2DI,S2DB,PARX,PARY,TMAX)
      ENDIF
C
************************************************************************
C
C      GREPS=0.20D0
C      IF (ILEV.EQ.NLEV) CALL GRDIST(DWORK(L(LCORVG)),KWORK(L(LNPR)),
C     *                              NVT,GREPS)
C
************************************************************************
C
      IF (ILEV.GE.2) 
     * CALL CHCOOR(DWORK(L(LCORVG)),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             KWORK(L(LADJ)),KWORK(L(LNPR)),NEL,NVT)
************************************************************************
C
      IF (ILEV.GE.NLMIN) THEN
C
C *** Save dimensions for all levels between NLMIN and NLMAX
C
      KNEL(ILEV) =NEL
      KNVT(ILEV) =NVT
      KNMT(ILEV) =NMT
      KNVEL(ILEV)=NVEL
      KNVBD(ILEV)=NVBD     
C
C
C *** Save arrays for all levels between NLMIN and NLMAX
C
C
      IF (ILEV.LT.NLEV) THEN
C
C       CALL ZCPY(LCORVG,'DCORVG',KLCVG(ILEV),'DCVG0 ')
C       IF (IER.NE.0) GOTO 99999
       IF (LCORMG.NE.0) THEN
        CALL ZCPY(LCORMG,'DCORMG',KLCMG(ILEV),'DCMG0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       CALL ZCPY(LVERT,'KVERT ',KLVERT(ILEV),'KVERT0')
       IF (IER.NE.0) GOTO 99999
       IF (LMID.NE.0) THEN
        CALL ZCPY(LMID,'KMID  ',KLMID(ILEV),'KMID0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (IADJ.EQ.1) THEN
        CALL ZCPY(LADJ,'KADJ  ',KLADJ(ILEV),'KADJ0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (LVEL.NE.0) THEN
        CALL ZCPY(LVEL,'KVEL  ',KLVEL(ILEV),'KVEL0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (LMEL.NE.0) THEN
        CALL ZCPY(LMEL,'KMEL  ',KLMEL(ILEV),'KMEL0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       CALL ZCPY(LNPR,'KNPR  ',KLNPR(ILEV),'KNPR0 ')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LMM,'KMM   ',KLMM(ILEV),'KMM0  ')
       IF (IER.NE.0) GOTO 99999
       IF (IBDP.GE.0) THEN
        CALL ZCPY(LVBD,'KVBD  ',KLVBD(ILEV),'KVBD0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (IBDP.GE.1) THEN
       CALL ZCPY(LEBD,'KEBD  ',KLEBD(ILEV),'KEBD0 ')
       IF (IER.NE.0) GOTO 99999
       ENDIF
       CALL ZCPY(LBCT,'KLBCT ',KLBCT(ILEV),'KBCT0 ')
       IF (IER.NE.0) GOTO 99999
       IF (IBDP.GE.2) THEN
        CALL ZCPY(LVBDP,'DVBDP ',KLVBDP(ILEV),'DVBDP0')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (IBDP.GE.2.AND.LMBDP.NE.0) THEN
        CALL ZCPY(LMBDP,'DMBDP ',KLMBDP(ILEV),'DMBDP0')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
      ELSE
C
       KLCVG(ILEV) =LCORVG
       KLCMG(ILEV) =LCORMG
       KLVERT(ILEV)=LVERT
       KLMID(ILEV) =LMID
       KLADJ(ILEV) =LADJ
       KLVEL(ILEV) =LVEL
       KLMEL(ILEV) =LMEL
       KLNPR(ILEV) =LNPR
       KLMM(ILEV)  =LMM
       KLVBD(ILEV) =LVBD
       KLEBD(ILEV) =LEBD
       KLBCT(ILEV) =LBCT
       KLVBDP(ILEV)=LVBDP
       KLMBDP(ILEV)=LMBDP
C
      ENDIF
C
      DO 20 ILEVH=1,NLEV-1
      KLCVG(ILEVH)=LCORVG
20    CONTINUE
C
      ENDIF
C
10    CONTINUE
C
C
99999 END
C
C
C
************************************************************************
C
      SUBROUTINE XMORS2(IRMESH)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299,NNLEV=9,NNWORK=1)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CTRIF(NNLEV)*15
C
      DATA CTRIF/ '#tries/#tria1  ','#tries/#tria2  ','#tries/#tria3  ',
     *            '#tries/#tria4  ','#tries/#tria5  ','#tries/#tria6  ',
     *            '#tries/#tria7  ','#tries/#tria8  ','#tries/#tria9  '/
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
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
      SAVE 
C
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      SUB='XMORS2'
      IF (ICHECK.GE.997) CALL OTRC('XMORS2','06/05/95')
C
C
C
      NLEV=NLMAX
C
      IF (IRMESH.GT.1) THEN
       IFMT=1
      ELSE
       IFMT=0
      ENDIF
C 
      DO 10 ILEV=NLMIN,NLMAX
C
C *** Read refined mesh in file CTRIF(ILEV) with format IFMT
      CALL XORS(65,CTRIF(ILEV),IFMT)
      REWIND(65)
      CLOSE(65)
C
C *** Save dimensions for all levels between NLMIN and NLMAX
      KNEL(ILEV) =NEL
      KNVT(ILEV) =NVT
      KNMT(ILEV) =NMT
      KNVEL(ILEV)=NVEL
      KNVBD(ILEV)=NVBD     
C
C
C *** Save arrays for all levels between NLMIN and NLMAX
      IF (ILEV.LT.NLMAX) THEN
C
       IF (LCORMG.NE.0) THEN
        CALL ZCPY(LCORMG,'DCORMG',KLCMG(ILEV),'DCMG0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       CALL ZCPY(LVERT,'KVERT ',KLVERT(ILEV),'KVERT0')
       IF (IER.NE.0) GOTO 99999
C
       IF (LMID.NE.0) THEN
        CALL ZCPY(LMID,'KMID  ',KLMID(ILEV),'KMID0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LADJ.NE.0) THEN
        CALL ZCPY(LADJ,'KADJ  ',KLADJ(ILEV),'KADJ0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LVEL.NE.0) THEN
        CALL ZCPY(LVEL,'KVEL  ',KLVEL(ILEV),'KVEL0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LMEL.NE.0) THEN
        CALL ZCPY(LMEL,'KMEL  ',KLMEL(ILEV),'KMEL0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       CALL ZCPY(LNPR,'KNPR  ',KLNPR(ILEV),'KNPR0 ')
       IF (IER.NE.0) GOTO 99999
C
       CALL ZCPY(LMM,'KMM   ',KLMM(ILEV),'KMM0  ')
       IF (IER.NE.0) GOTO 99999
C
       IF (LVBD.NE.0) THEN
        CALL ZCPY(LVBD,'KVBD  ',KLVBD(ILEV),'KVBD0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LEBD.NE.0) THEN
        CALL ZCPY(LEBD,'KEBD  ',KLEBD(ILEV),'KEBD0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       CALL ZCPY(LBCT,'KLBCT ',KLBCT(ILEV),'KBCT0 ')
       IF (IER.NE.0) GOTO 99999
C
       IF (LVBDP.NE.0) THEN
        CALL ZCPY(LVBDP,'DVBDP ',KLVBDP(ILEV),'DVBDP0')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (LMBDP.NE.0) THEN
        CALL ZCPY(LMBDP,'DMBDP ',KLMBDP(ILEV),'DMBDP0')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
C
      ELSE
C
C
       KLCVG(ILEV) =LCORVG
       KLCMG(ILEV) =LCORMG
       KLVERT(ILEV)=LVERT
       KLMID(ILEV) =LMID
       KLADJ(ILEV) =LADJ
       KLVEL(ILEV) =LVEL
       KLMEL(ILEV) =LMEL
       KLNPR(ILEV) =LNPR
       KLMM(ILEV)  =LMM
       KLVBD(ILEV) =LVBD
       KLEBD(ILEV) =LEBD
       KLBCT(ILEV) =LBCT
       KLVBDP(ILEV)=LVBDP
       KLMBDP(ILEV)=LMBDP
C
      ENDIF
C
10    CONTINUE
C
      DO 20 ILEV=1,NLMAX-1
      KLCVG(ILEV)=LCORVG
20    CONTINUE
C
C
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
* ORALCn                                                               *
*                                                                      *
* Purpose  Get arrays from file, previously stored by XOWA or OWAn     *
*                                                                      *
* Subroutines/functions called  ORALC0                                 *
*                                                                      *
* Version from  10/11/94                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----    ----                                                        *
* MFILE   I*4    I/O unit - previouly opened (by OF0)                  *
* IFMT    I*4    0  Format free input                                  *
*                1  Formatted input                                    *
* A1      R*8    Scaling factor: DX=DX+A1*DX(readin)                   *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
* DX      R*8    Double precision array - version ORA1                 *
* VX      R*4    Single precision array - version ORA2                 *
* ARR     C*6    Name of array read from data file                     *
* IER     I*4    Error indicator                                       *
*                -110  Error while reading from unit  MFILE            *
*                -113  Wrong value of  ITYPE  read from  MFILE         *
*                                                                      *
************************************************************************
C
      SUBROUTINE ORALC1(DX,A1,ARR,MFILE,IFMT)
C
C *** Double precision version
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6
      DIMENSION DX(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='ORALC1'
      IF (ICHECK.GE.997) CALL OTRC('ORALC1','10/11/94')
C
      IER=0
C
      IF (IFMT.NE.0) GOTO 99997
C
      READ (MFILE,ERR=99997,END=99997) ARR,ITYPE,ILEN,ILEN8,JRECL
C
      IF (ITYPE.NE.1) GOTO 99997
C
      CALL ORALCD(DX,A1,ILEN8,MFILE,JRECL)
      IF (IER.NE.0) GOTO 99999
C
      LNR=0
      WRITE (CPARAM,'(A6,2I15)') ARR,LNR,MFILE
      CALL OMSG(7,'ORALC1')
      GOTO 99999
C
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'XORALC')
C
99999 END
C
C
C
      SUBROUTINE ORALC2(VX,A1,ARR,MFILE,IFMT)
C
C *** Single precision version
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6
      DIMENSION VX(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='ORALC2'
      IF (ICHECK.GE.997) CALL OTRC('ORALC2','10/11/94')
C
      IER=0
C
      IF (IFMT.NE.0) GOTO 99997
C
      READ (MFILE,ERR=99997,END=99997) ARR,ITYPE,ILEN,ILEN8,JRECL
C
      IF (ITYPE.NE.2) GOTO 99997
C
      CALL ORALCV(VX,A1,ILEN8,MFILE,JRECL)
      IF (IER.NE.0) GOTO 99999
C
      LNR=0
      WRITE (CPARAM,'(A6,2I15)') ARR,LNR,MFILE
      CALL OMSG(7,'ORALC2')
      GOTO 99999
C
99997 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'XORALC')
C
99999 END
C
C
C
      SUBROUTINE ORALCD(DX,A1,N,MFILE,IRECL8)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DX(*),DXH(512)
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /CHAR/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('ORALCD','10/11/94')
C
      IREC=N/IRECL8
      DO 1 JREC=1,IREC
      J1=(JREC-1)*IRECL8
      READ (MFILE,ERR=99998,END=99998) (DXH(I),I=1,IRECL8)
      DO 1 I=1,IRECL8
      DX(J1+I)=DX(J1+I)+DXH(I)
1     CONTINUE
      IF (MOD(N,IRECL8).EQ.0) GOTO 99999
      J1=IREC*IRECL8
      READ (MFILE,ERR=99998,END=99998) (DXH(I),I=1,MOD(N,IRECL8))
      DO 2 I=1,MOD(N,IRECL8)
      DX(J1+I)=DX(J1+I)+DXH(I)
2     CONTINUE
      GOTO 99999
C
99998 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'ORALCD')
C
99999 END
C
C
C
      SUBROUTINE ORALCV(VX,A1,N,MFILE,IRECL8)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VX(*),VXH(512)
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /CHAR/,/ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('ORALCV','10/11/94')
C
      IREC=N/IRECL8
      DO 1 JREC=1,IREC
      J1=(JREC-1)*IRECL8
      READ (MFILE,ERR=99998,END=99998) (VXH(I),I=1,IRECL8)
      DO 1 I=1,IRECL8
      VX(J1+I)=VX(J1+I)+VXH(I)
1     CONTINUE
      IF (MOD(N,IRECL8).EQ.0) GOTO 99999
      J1=IREC*IRECL8
      READ (MFILE,ERR=99998,END=99998) (VXH(I),I=1,MOD(N,IRECL8))
      DO 2 I=1,MOD(N,IRECL8)
      VX(J1+I)=VX(J1+I)+VXH(I)
2     CONTINUE
      GOTO 99999
C
99998 WRITE (CPARAM,'(I15)') MFILE
      CALL WERR(-110,'ORALCV')
C
99999 END
C
C
C
      SUBROUTINE OWA2V(VX,ARR,NX,MFILE,IFMT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6,CFORM*15
      DIMENSION VX(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='OWA2V'
      IF (ICHECK.GE.997) CALL OTRC('OWA2V ','10/11/94')
C
      IER=0
      BFMT=IFMT.EQ.1
C
      IF (BFMT) THEN
       CFORM=FMT(1)
       WRITE (MFILE,'(2A15,2I15)') ARR,CFORM,1,NX
      ELSE
       WRITE (MFILE) ARR,2,NX,NX,IRECL8
      ENDIF
C
      IF (BFMT) THEN
       WRITE (MFILE,CFORM) (VX(I),I=1,NX)
      ELSE
       CALL OWA0V(VX,NX,MFILE,IRECL8)
      ENDIF
C
      LNR=0
      WRITE (CPARAM,'(A6,2I15)') ARR,LNR,MFILE
      CALL OMSG(8,'XOWA2V')
C
99999 END
C
C
C
      SUBROUTINE OWA0V(VX,N,MFILE,IRECL8)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('OWA0V ','10/11/94')
C
      IREC=N/IRECL8
      DO 1 JREC=1,IREC
      J1=(JREC-1)*IRECL8
1     WRITE (MFILE) (VX(J1+I),I=1,IRECL8)
      IF (MOD(N,IRECL8).EQ.0) GOTO 99999
      J1=IREC*IRECL8
      WRITE (MFILE) (VX(J1+I),I=1,MOD(N,IRECL8))
C
99999 END
C
C
C
      SUBROUTINE XVBM0(LB,NEQ,NBLOC,ICLEAR,ELE,COEFF,BCON,KB,KBN,ICUB,
     *                 ARR,BSNGL)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6
C
      PARAMETER (NNARR=299,NNAB=21,NNDER=6)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION LB(*),KB(NNAB,*),KBN(*),BCON(*),ARR(*),BSNGL(*)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL COEFF,ELE
      SAVE
C
      SUB='XVBM0'
      IF (ICHECK.GE.997) CALL OTRC('XVBM0 ','12/04/90')
      IER=0
C
      DO 1 IBLOC=1,NBLOC
      IF (LB(IBLOC).EQ.0) THEN
       CALL ZNEW(NEQ,1,LB(IBLOC),ARR(IBLOC))
       IF (IER.NE.0) GOTO 99999
      ELSE
C ***  Check input parameter
       CALL ZTYPE(LB(IBLOC),ITYPE)
       IF (ITYPE.NE.1) THEN
        WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
        CALL WERR(-114,'XVBM0 ')
        GOTO 99999
       ENDIF
       CALL ZLEN(LB(IBLOC),ILEN)
       IF (ILEN.LT.NEQ) THEN
        WRITE (CPARAM,'(A6,I15)') ARR(IBLOC),IBLOC
        CALL WERR(-115,'XVBM0 ')
        GOTO 99999
       ENDIF
       IF (ICLEAR.EQ.1) CALL ZCLEAR(LB(IBLOC),ARR(IBLOC))
      ENDIF
1     CONTINUE
C
      CALL ZNEW(NBLOC,-3,LOFF,'KOFF  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NBLOC*NNDER,1,LOECON,'COECON')
      IF (IER.NE.0) GOTO 99999
C
      DO 2 IBLOC=1,NBLOC
      KWORK(L(LOFF)+IBLOC-1)=L(LB(IBLOC))-1
2     CONTINUE
C
      CALL VBM0(DWORK(1),NBLOC,KWORK(L(LOFF)),KWORK(L(LVERT)),
     *          KWORK(L(LMID)),DWORK(L(LCORVG)),
     *          ELE,COEFF,BCON,DWORK(L(LOECON)),KB,KBN,ICUB)
      IF (IER.NE.0) GOTO 99999
C
      CALL ZDISP(0,LOECON,'COECON')
      CALL ZDISP(0,LOFF,'KOFF  ')
C
      DO 3 IBLOC=NBLOC,1,-1
      IF (BSNGL(IBLOC)) THEN
       CALL ZCTYPE(2,LB(IBLOC),ARR(IBLOC))
       IF (IER.NE.0) GOTO 99999
      ENDIF
3     CONTINUE
C
99999 END
C
C
C
      SUBROUTINE VBM0(DB,NBLOC,KOFF,KVERT,KMID,DCORVG,
     *                ELE,COEFF,BCON,COECON,KB,KBN,ICUB)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NNAB=21,NNCOF=6)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KOFF(*)
      DIMENSION KDFG(NNBAS),KDFL(NNBAS),DB(*),DCORVG(2,*)
      DIMENSION BCON(*),KB(NNAB,*),KBN(*),COECON(NNDER,*)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG,KDFL,IDFL
      COMMON /COFBAS/ COB(NNBAS,NNCOF)
      EXTERNAL ELE,COEFF
      SAVE
C
      SUB='VB0M'
      IF (ICHECK.GE.997) CALL OTRC('VB0M  ','12/04/90')
C
C *** Preparation - evaluation of parameters
      IER=0
C
C *** Which derivatives of basis functions are needed?
      DO 1 IDER=1,NNDER
1     BDER(IDER)=.FALSE.
      DO 2 IBLOC=1,NBLOC
      DO 3 IBN=1,KBN(IBLOC)
      IB=KB(IBN,IBLOC)
      IF (IB.LE.0.OR.IB.GT.NNDER) THEN
       WRITE (CPARAM,'(I15)') IBLOC
       CALL WERR(-117,'VB0M  ')
       GOTO 99999
      ENDIF
3     BDER(IB)=.TRUE.
2     CONTINUE
C
C *** Dummy call of ELE sets number of element
      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C
      BFIRST=.TRUE.
C *** Dummy call of COEFF for nonlinear problems
C *** COEFF must set BDER(IDER)=.TRUE. if derivative IDER is needed
      AUX=COEFF(0D0,0D0,-1,0,BFIRST)
      BCON0=.TRUE.
C
      DO 4 IBLOC=1,NBLOC
      IF (BCON(IBLOC)) THEN
       DO 5 IBN=1,KBN(IBLOC)
       IB=KB(IBN,IBLOC)
5      COECON(IB,IBLOC)=COEFF(0D0,0D0,IB,IBLOC,BFIRST)
      ELSE
       BCON0=.FALSE.
      ENDIF
4     CONTINUE
************************************************************************
C *** Calculation of the linear form
************************************************************************
C
C *** Loop over all elements
      DO 100 IEL=1,NEL
      CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.LT.0) GOTO 99999
C
C *** Evaluation of coordinates of the vertices
      DO 110 IVE=1,NVE
      JP=KVERT(IVE,IEL)
      KVE(IVE)=JP
      DX(IVE)=DCORVG(1,JP)
      DY(IVE)=DCORVG(2,JP)
110   CONTINUE
C
      DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
      DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
      DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
      DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
C
C *** Dummy call - ELE may save arithmetic operations
      CALL ELE(0D0,0D0,-2)
      IF (IER.LT.0) GOTO 99999
C
C *** Loop over all cubature points
      DO 200 ICUBP=1,NCUBP
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
C *** Cubature points on the actual element + weights
      XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *  +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
      YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1+0.5D0*
     *         (DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2
      OM=DOMEGA(ICUBP)*DETJ
C
      CALL ELE(XX,YY,-3)
      IF (IER.LT.0) GOTO 99999
C
C *** Summing up over all multiindices
      BFIRST=.TRUE.
      DO 300 IBLOC=1,NBLOC
C
      DO 301 IBN=1,KBN(IBLOC)
      IB=KB(IBN,IBLOC)
      IF (.NOT.BCON(IBLOC)) THEN
       AUX=COEFF(XX,YY,IB,IBLOC,BFIRST)*OM
      ELSE
       AUX=COECON(IB,IBLOC)*OM
      ENDIF
C
      DO 310 JDOFE=1,IDFL
      IGLOB=KDFG(JDOFE)+KOFF(IBLOC)
      DB(IGLOB)=DB(IGLOB)+DBAS(KDFL(JDOFE),IB)*AUX
310   CONTINUE
301   CONTINUE
C
      BFIRST=.FALSE.
300   CONTINUE
C
200   CONTINUE
100   CONTINUE
C
99999 END
