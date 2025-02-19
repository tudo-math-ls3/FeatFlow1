************************************************************************
      SUBROUTINE INIT1 (MDATA,CDATA,MFILE,MSHOW,IPROJ,IWORKG,IWMAXG)
************************************************************************
*
*   Purpose: - generates geometry for all levels
*            - allocates all arrays
*            - reads a start vector if ISTART=1 or 2
*            - generates linear matrices for all levels
*            - stores pointers for all arrays on COMMON blocks
*            - sets boundary parameters for all levels
*            - sets Dirichlet bc's for the finest level
*            - generates rhs for the finest level
*            - opens a file with unit number MFILE for user output
*            - etc
*
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299,NNAB=21,NNLEV=9,NNWORK=1)
      PARAMETER (NBLOCA=1,NBLOCB=2,NBLOCF=3)
      PARAMETER (NBLA1=NBLOCA*NNLEV,NBLB1=NBLOCB*NNLEV,
     *           NBLF1=NBLOCF*NNLEV)
C
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CDATA*60,CFILE*60
C
C *** Names of matrices and vectors (for messages only)
      CHARACTER CARRST*12,CARRM*12,CARRDB*12
      CHARACTER CFILST*12,CFILM*12
      CHARACTER ARRDA*6,ARRDST*6,ARRDB*6,ARRDF*6,ARRDFP*6,ARRDM*6
      DIMENSION CARRST(NNLEV),CARRM(NNLEV),CARRDB(NNLEV)
      DIMENSION CFILST(NNLEV),CFILM(NNLEV)
      DIMENSION ARRDA(NBLOCA),ARRDST(NBLOCA),ARRDM(NBLOCA)
      DIMENSION ARRDB(NBLOCB),ARRDF(NBLOCF)
C
C *** Names of matrices and vectors (for output and tracing)
      DATA CARRST/'#data/VST.1 ','#data/VST.2 ','#data/VST.3 ',
     *            '#data/VST.4 ','#data/VST.5 ','#data/VST.6 ',
     *            '#data/VST.7 ','#data/VST.8 ','#data/VST.9 '/
      DATA CARRM /'#data/VM.1  ','#data/VM.2  ','#data/VM.3  ',
     *            '#data/VM.4  ','#data/VM.5  ','#data/VM.6  ',
     *            '#data/VM.7  ','#data/VM.8  ','#data/VM.9  '/
      DATA CARRDB/'#data/VB.1  ','#data/VB.2  ','#data/VB.3  ',
     *            '#data/VB.4  ','#data/VB.5  ','#data/VB.6  ',
     *            '#data/VB.7  ','#data/VB.8  ','#data/VB.9  '/
      DATA CFILST/'#ns/ST1     ','#ns/ST2     ','#ns/ST3     ',
     *            '#ns/ST4     ','#ns/ST5     ','#ns/ST6     ',
     *            '#ns/ST7     ','#ns/ST8     ','#ns/ST9     '/
      DATA CFILM /'#ns/MA1     ','#ns/MA2     ','#ns/MA3     ',
     *            '#ns/MA4     ','#ns/MA5     ','#ns/MA6     ',
     *            '#ns/MA7     ','#ns/MA8     ','#ns/MA9     '/
      DATA ARRDA/'VA    '/,ARRDST/'DST   '/
      DATA ARRDF/'DF1   ','DF2   ','DTF   '/
      DATA ARRDM/'DM    '/
      DATA ARRDB/'DB1   ','DB2   '/,ARRDFP/'DFP   '/
C
C *** Structure of bilinear and linear forms 
      DIMENSION KABSTN(NBLOCA),KABST(2,NNAB,NBLOCA)
      DIMENSION KABBN(NBLOCB),KABB(2,NNAB,NBLOCB)
      DIMENSION KFN(NBLOCF),KF(NNAB,NBLOCF)
      DATA KABSTN/3/, KABBN/1,1/, KFN/1,1,1/
C
      DIMENSION BCONST(NBLOCA),BCONB(NBLOCB),BCONF(NBLOCF)
      DATA BCONST/.TRUE./, BCONB/.TRUE.,.TRUE./
      DATA BCONF/.FALSE.,.FALSE.,.FALSE./
C
      DIMENSION BSNGLA(NBLOCA,NNLEV),BSNGLB(NBLOCB,NNLEV),BSNGLF(NBLOCF)
      DATA BSNGLA/NBLA1*.FALSE./,BSNGLB/NBLB1*.FALSE./
      DATA BSNGLF/.FALSE.,.FALSE.,.FALSE./
C
      DIMENSION KLST(NBLOCA,NNLEV),KLMASS(NBLOCA,NNLEV)
      DIMENSION KLB(NBLOCB,NNLEV),KLF(NBLOCF,NNLEV),LF(NBLOCF)
      DIMENSION KLCOLA(NNLEV),KLLDA(NNLEV),KNA(NNLEV),KNEQ(NNLEV)
      DIMENSION KLCOLB(NNLEV),KLLDB(NNLEV),KNB(NNLEV),KNEQB(NNLEV)
      DIMENSION KLMH(NNLEV)
      DATA KLB/NBLB1*0/,KLF/NBLF1*0/,LF(1)/0/,LF(2)/0/,LF(3)/0/
      DATA KLMH/NNLEV*0/
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
      COMMON /NSTIME/ TTGRID,TTPOST,TTADF,TTUPW,TTBDR,TTLC,TTILU,
     *                TTMGU,TTSU,TTEU,TTDU,TTPU,TTRU,
     *                TTMGP,TTSP,TTEP,TTDP,TTPP,TTRP
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
      COMMON /MGBDRY/ KLMBD(NNLEV),KLDBD(NNLEV),KNMBD(NNLEV),
     *                KLNPRO(NNLEV),INEUM
      COMMON /MGFLD/  KLA(NNLEV),KLST,KLMASS,KLM(NNLEV),KLCOLA,KLLDA,
     *                KLB1(NNLEV),KLB2(NNLEV),KLCOLB,KLLDB,
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGDIM/  KNA,KNB,KNU(NNLEV),KNP(NNLEV),
     *                KNUP(NNLEV)
      COMMON /MGPROJ/ KLC(NNLEV),KLCOLC(NNLEV),KLLDC(NNLEV),KNC(NNLEV)
      COMMON /MGILUU/ ISORTU,KLAILU(NNLEV),KLTRA1(NNLEV),KLTRA2(NNLEV)
      COMMON /MGILUP/ ISORTP,KLCILU(NNLEV),KLTRC1(NNLEV),KLTRC2(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KMASS1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
C
C *** user COMMON block
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
      CHARACTER CPARM1*60,CMESH1*60,CFILE1*60,CSTART*60,CSOL*60
      COMMON /FILES/ IMESH1,MMESH1,CPARM1,CMESH1,MFILE1,CFILE1,
     *               ISTART,MSTART,CSTART,ISOL,MSOL,CSOL
C=======================================================================
      INCLUDE 'bouss.inc'
C=======================================================================
C
      SAVE 
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C     E X T E R N A L S
C-----------------------------------------------------------------------
C *** Parametrization of the domain
      EXTERNAL PARX,PARY,TMAX
C *** Control of refinement - here: regular refinement
      EXTERNAL S2DI0,S2DB0
C *** Coefficient of stiffness matrix, right hand side, exact solution
      EXTERNAL COEFST,COEFFB,RHS,UE
C *** definition of finite elements
      EXTERNAL E030,E031,EM30,EM31,E010
C=======================================================================
C     Initialization
C=======================================================================
      SUB='INIT1 '
C
      CALL ZTIME(TTT0)
C
C *** Structure of bilinear and linear forms
C     B1 block
      KABB(1,1,1)=2
      KABB(2,1,1)=1
C     B2 block
      KABB(1,1,2)=3
      KABB(2,1,2)=1
C
      KF(1,1)=1
      KF(1,2)=1
      KF(1,3)=1
C
      ISYMMA=0
C
C=======================================================================
C     Data input
C=======================================================================
C
C      CALL  OF0  (MDATA,CDATA,1)
      CALL  OPNDAT (MDATA,MSHOW,CDATA)
      CALL  GDAT (MDATA,MSHOW,IRMESH,IPROJ)
      CLOSE (MDATA)
      CFILE=CFILE1
      MFILE=MFILE1
C
C=======================================================================
C     Grid generation
C=======================================================================
C
      IF (IMESH1.EQ.1) THEN
       CALL RDPARM (CPARM1,MDATA)
       CLOSE(MDATA)
      ENDIF
C
      CALL  XORSC (MMESH1,CMESH1)
      IF (IER.NE.0) GOTO 99999
      CLOSE(MMESH1)
C
C
      IMID=1
      IADJ=1
      IVEL=0
      IDISP=1
      IBDP=2
C
      IF (IRMESH.EQ.0) THEN
       CALL XMSB2(IMID,IADJ,IVEL,IDISP,IBDP,S2DI0,S2DB0,PARX,PARY,TMAX)
       IF (IER.NE.0) GOTO 99999
      ELSE
       CALL XMORS2(IRMESH)
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
C
      DO 11  II=NLMIN,NLMAX
      ILEV=II
      NVT=KNVT(II)
      NMT=KNMT(II)
      NEL=KNEL(II)
      NVBD=KNVBD(II)
      KNU(II)=NMT
      KNP(II)=NEL
      NUP=NMT+NMT+NEL
      NUPT=NUP+NMT
      KNUP(II)=NUP
      KNUPT(II)=NUPT
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'ILEV,NVT,NMT,NEL,NVBD: ', II,NVT,NMT,NEL,NVBD
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'ILEV,NVT,NMT,NEL,NVBD: ', II,NVT,NMT,NEL,NVBD
C
      IF (IBDR.GE.2) THEN
       CALL ZNEW(NVT+NMT,-3,KLNPRO(II),'KNPRO ')
       IF (IER.NE.0) GOTO 99999
       CALL XLCP3(KLNPR(II),KLNPRO(II),NVT+NMT)
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      CALL ZNEW(NVBD,3,KLMBD(II),'KMBD  ')
      CALL ZNEW(NVBD,1,KLDBD(II),'DDBD  ')
      IF (IER.NE.0) GOTO 99999
      CALL BDRNEU(KWORK(L(KLMBD(II))),KWORK(L(KLVBD(II))),
     *            KWORK(L(KLEBD(II))),KWORK(L(KLVERT(II))),
     *            KWORK(L(KLMID(II))),KWORK(L(KLNPR(II))),
     *            DWORK(L(KLDBD(II))),DWORK(L(KLMBDP(II))),NMBD,INEUM)
      IF (IER.NE.0) GOTO 99999
      KNMBD(II)=NMBD      
11    CONTINUE
C
C
      CALL ZTIME(TTT1)
      TTGRID=TTT1-TTT0
C
      IF (MSHOW.GE.2) WRITE(MTERM,*)
      IF (MSHOW.GE.0) WRITE(MFILE,*)
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'time for grid generation : ',
     *                                TTGRID
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'time for grid generation : ', 
     *                                TTGRID
      IF (MSHOW.GE.2) WRITE(MTERM,*)
      IF (MSHOW.GE.0) WRITE(MFILE,*)
      IWORKG=IWORK
      IWMAXG=IWMAX
C
C=======================================================================
C     Generation of: - pointer structures
C                    - STOKES,B1,B2 blocks 
C=======================================================================
C
C *** Generation of Laplace/mass block
C
      CALL ZTIME(TTT0L)
      CALL ZTIME(TTT0)
      CALL XMAP7(KLCOLA,KLLDA,KNA,KNEQ,E031,ISYMMA)
      IF (IER.NE.0) GOTO 99999
C
C=======================================================================
C *** Generation of lumped + real mass matrix
C=======================================================================
C
      KABST (1,1,1)=1
      KABST (2,1,1)=1
      KABSTN(1)    =1
C
      ICLRA=1
      IF (IMASSL.EQ.0)
     *  CALL XMAB07(KLMASS,KLCOLA,KLLDA,KNA,KNEQ,NBLOCA,ICLRA,E031,
     *              COEFST,BCONST,KABST,KABSTN,ABS(ICUBML),
     *              ISYMMA,CARRM,BSNGLA)
      IF (IMASSL.EQ.1)
     *  CALL XMAB07(KLMASS,KLCOLA,KLLDA,KNA,KNEQ,NBLOCA,ICLRA,E030,
     *              COEFST,BCONST,KABST,KABSTN,ABS(ICUBML),
     *              ISYMMA,CARRM,BSNGLA)
      IF (IMASSL.EQ.2)
     *  CALL XMABM7(KLMASS,KLCOLA,KLLDA,KNA,KNEQ,NBLOCA,ICLRA,EM31,
     *              COEFST,BCONST,KABST,KABSTN,ABS(ICUBML),
     *              ISYMMA,CARRM,BSNGLA)
      IF (IMASSL.EQ.3) 
     *  CALL XMABM7(KLMASS,KLCOLA,KLLDA,KNA,KNEQ,NBLOCA,ICLRA,EM30,
     *              COEFST,BCONST,KABST,KABSTN,ABS(ICUBML),
     *              ISYMMA,CARRM,BSNGLA)
      IF (IER.NE.0) GOTO 99999
C
C=======================================================================
C
      DO 20  II=NLMIN,NLMAX
      CALL  ZNEW (KNU(II),1,KLM(II),'DMASS ')
      IF (IER.NE.0) GOTO 99999
C
      DO 21 INU=1,KNU(II)
C
      IF (ICUBML.GT.0) THEN
       DMH=0D0
       DO 22 ILD=KWORK(L(KLLDA(II))+INU-1),KWORK(L(KLLDA(II))+INU)-1
22     DMH=DMH+DWORK(L(KLMASS(1,II))+ILD-1)
      ELSE
       ILD=KWORK(L(KLLDA(II))+INU-1)
       DMH=DWORK(L(KLMASS(1,II))+ILD-1)
      ENDIF
C
      DWORK(L(KLM(II))+INU-1)=DMH
21    CONTINUE
C
      IF ( (IMASS.EQ.0).OR.(IPRECA.EQ.4).OR.
     *    ((IMASS.EQ.1).AND.((IMASSL.NE.IELT).OR.
     *                       (ABS(ICUBML).NE.ICUBM)))) THEN
       CALL ZDISP (0,KLMASS(1,II),'MASMAT')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
20    CONTINUE
C
C=======================================================================
C
      IF ( (IPRECA.NE.4).AND.
     *    ((IMASS.EQ.1 ).AND.((IMASSL.NE.IELT).OR.
     *                        (ABS(ICUBML).NE.ICUBM)))) THEN
       ICLRA=1
       IF (IELT.EQ.0)
     *  CALL XMAB07(KLMASS,KLCOLA,KLLDA,KNA,KNEQ,NBLOCA,ICLRA,E031,
     *              COEFST,BCONST,KABST,KABSTN,ICUBM,ISYMMA,CARRM,
     *              BSNGLA)
       IF (IELT.EQ.1)
     *  CALL XMAB07(KLMASS,KLCOLA,KLLDA,KNA,KNEQ,NBLOCA,ICLRA,E030,
     *              COEFST,BCONST,KABST,KABSTN,ICUBM,ISYMMA,CARRM,
     *              BSNGLA)
       IF (IELT.EQ.2)
     *  CALL XMABM7(KLMASS,KLCOLA,KLLDA,KNA,KNEQ,NBLOCA,ICLRA,EM31,
     *              COEFST,BCONST,KABST,KABSTN,ICUBM,ISYMMA,CARRM,
     *              BSNGLA)
       IF (IELT.EQ.3) 
     *  CALL XMABM7(KLMASS,KLCOLA,KLLDA,KNA,KNEQ,NBLOCA,ICLRA,EM30,
     *              COEFST,BCONST,KABST,KABSTN,ICUBM,ISYMMA,CARRM,
     *              BSNGLA)
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
C=======================================================================
C *** Generation of block ST
C=======================================================================
C
      IF (IPRECA.NE.4) THEN
       KABST (1,1,1)=2
       KABST (2,1,1)=2
       KABST (1,2,1)=3
       KABST (2,2,1)=3
       KABSTN(1)    =2
C
       ICLRA=1
       IF (IELT.EQ.0) 
     *  CALL XMAB07(KLST,KLCOLA,KLLDA,KNA,KNEQ,NBLOCA,ICLRA,E031,COEFST,
     *              BCONST,KABST,KABSTN,ICUBA,ISYMMA,CARRST,BSNGLA)
       IF (IELT.EQ.1) 
     *  CALL XMAB07(KLST,KLCOLA,KLLDA,KNA,KNEQ,NBLOCA,ICLRA,E030,COEFST,
     *              BCONST,KABST,KABSTN,ICUBA,ISYMMA,CARRST,BSNGLA)
       IF (IELT.EQ.2)
     *  CALL XMABM7(KLST,KLCOLA,KLLDA,KNA,KNEQ,NBLOCA,ICLRA,EM31,COEFST,
     *              BCONST,KABST,KABSTN,ICUBA,ISYMMA,CARRST,BSNGLA)
       IF (IELT.EQ.3)
     *  CALL XMABM7(KLST,KLCOLA,KLLDA,KNA,KNEQ,NBLOCA,ICLRA,EM30,COEFST,
     *              BCONST,KABST,KABSTN,ICUBA,ISYMMA,CARRST,BSNGLA)
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
C=======================================================================
C *** Generation of blocks B1,B2
C=======================================================================
C
      IF (IPRECB.LT.2) THEN
       IF ((IELT.EQ.0).OR.(IELT.EQ.2)) 
     *      CALL XMAP9(KLCOLB,KLLDB,KNB,KNEQB,E031,E010)
       IF ((IELT.EQ.1).OR.(IELT.EQ.3)) 
     *      CALL XMAP9(KLCOLB,KLLDB,KNB,KNEQB,E030,E010)
       IF (IER.NE.0) GOTO 99999
C
       ICLRB=1
       IF ((IELT.EQ.0).OR.(IELT.EQ.2)) 
     *  CALL XMAB09(KLB,KLCOLB,KLLDB,KNB,NBLOCB,ICLRB,E031,E010,E010,
     *              COEFFB,BCONB,KABB,KABBN,ICUBB,CARRDB,BSNGLB)
       IF ((IELT.EQ.1).OR.(IELT.EQ.3)) 
     *  CALL XMAB09(KLB,KLCOLB,KLLDB,KNB,NBLOCB,ICLRB,E030,E010,E010,
     *              COEFFB,BCONB,KABB,KABBN,ICUBB,CARRDB,BSNGLB)
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
C
      DO 13 ILEV=NLMIN,NLMAX
C
      IF (IPRECB.GE.2) THEN
       KNEQB(ILEV)=KNU(ILEV)
       KNB  (ILEV)=2*KNU(ILEV)
       CALL ZNEW(KNB(ILEV)  ,1,KLB(1,ILEV) ,ARRDB(1))
       CALL ZNEW(KNB(ILEV)  ,1,KLB(2,ILEV) ,ARRDB(2))
       CALL ZNEW(KNB(ILEV)  ,3,KLCOLB(ILEV),'KCOLB ')
       CALL ZNEW(KNU(ILEV)+1,3,KLLDB (ILEV),'LLDB  ')
       IF (IER.NE.0) GOTO 99999
C
       ISETLV=2
       CALL SETLEV(ISETLV)
       CALL BBUILD(KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *             DWORK(L(LCORVG)),
     *             DWORK(L(KLB(1,ILEV))),DWORK(L(KLB(2,ILEV))),
     *             KWORK(KCOLB),KWORK(KLDB),KNB(ILEV),NEL,NVT,NMT)
C
       CALL ZDISP (KNB(ILEV),KLB(1,ILEV) ,ARRDB(1))
       CALL ZDISP (KNB(ILEV),KLB(2,ILEV) ,ARRDB(2))
       CALL ZDISP (KNB(ILEV),KLCOLB(ILEV),'KCOLB ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'ILEV,NU,NA,NB:',ILEV,KNU(ILEV),
     *                                KNA(ILEV),KNB(ILEV)
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'ILEV,NU,NA,NB:',ILEV,KNU(ILEV),
     *                                KNA(ILEV),KNB(ILEV)
C
      NEQ =KNEQ(ILEV)
      NEQB=KNEQB(ILEV)
      IF (NEQB.NE.NEQ.OR.NEQ.NE.KNU(ILEV).OR.KNEL(ILEV).NE.KNP(ILEV)) 
     *    THEN
        WRITE(MTERM,*) 'ERROR in INIT1: NEQ.NE.NEQB'
        STOP
      ENDIF
C
      KLB1(ILEV)=KLB(1,ILEV)
      KLB2(ILEV)=KLB(2,ILEV)
13    CONTINUE
C
C
      IF (MSHOW.GE.2) WRITE(MTERM,*)
      IF (MSHOW.GE.0) WRITE(MFILE,*)
C 
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
C=======================================================================
C     C matrix generation
C=======================================================================
C
      DO 15 ILEV=NLMIN,NLMAX
       CALL ZTIME(TTT0)
       ISETLV=2
       CALL SETLEV(ISETLV)
C
       KNC(ILEV)=5*NP
       CALL ZNEW(KNC(ILEV)  ,1,KLC(ILEV)   ,'CC    ')
       CALL ZNEW(KNC(ILEV)  ,3,KLCOLC(ILEV),'KCOLC ')
       CALL ZNEW(NP+1       ,3,KLLDC(ILEV) ,'KLDC  ')
       IF (IER.NE.0) GOTO 99999
C
       CALL PROJST(KWORK(L(KLCOLC(ILEV))),KWORK(L(KLLDC(ILEV))),
     *             KWORK(L(LADJ)),NEL,KNC(ILEV))
C
       CALL ZDISP (KNC(ILEV),KLC(ILEV)   ,'CC    ')
       CALL ZDISP (KNC(ILEV),KLCOLC(ILEV),'KCOLC  ')
       IF (IER.NE.0) GOTO 99999
C
       CALL PROJMA(DWORK(L(KLC(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *             KWORK(L(KLLDC(ILEV))),KWORK(L(LNPR)),
     *             KWORK(L(LMID)),KWORK(L(LADJ)),
     *             DWORK(KM1),DWORK(KB1),DWORK(KB2),
     *             KWORK(KCOLB),KWORK(KLDB),NEL,NVT,NMT)
C
       IF (MSHOW.GE.2) 
     *  WRITE(MTERM,*) 'ILEV,NP,NC:',ILEV,KNP(ILEV),KNC(ILEV)
       IF (MSHOW.GE.0) 
     *  WRITE(MFILE,*) 'ILEV,NP,NC:',ILEV,KNP(ILEV),KNC(ILEV)
C
       CALL ZTIME(TTT1)
       TTADF=TTADF+TTT1-TTT0
C
C=======================================================================
C     ILU(C) matrix generation + grid sorting for P
C=======================================================================
C
       CALL ZTIME(TTT0)
       IF ((ISMP.EQ.4).OR.(ISLP.EQ.3).OR.(ISLP.EQ.4).OR.
     *     (ISORTP.GT.0)) THEN
C
        IF (ISORTP.GT.0) THEN
         KLTRC1(ILEV)=0
         KLTRC2(ILEV)=0
         CALL ZNEW(KNP(ILEV), 3,KLTRC1(ILEV),'KKTRC1')
         CALL ZNEW(KNP(ILEV), 3,KLTRC2(ILEV),'KKTRC2')
         IF (IER.NE.0) GOTO 99999
C
C=======================================================================
C     coordinate sorting for ISORTP = 1 or 2
C=======================================================================
C
         IF ((ISORTP.EQ.1).OR.(ISORTP.EQ.2)) THEN
          LCOREC=0
          CALL ZNEW(2*KNP(ILEV),-1,LCOREC,'DCOREC')
          IF (IER.NE.0) GOTO 99999
C
          CALL TRCORE(KWORK(L(LVERT)),DWORK(L(LCORVG)),DWORK(L(LCOREC)),
     *                KNEL(ILEV))
C
          INUM=ISORTP
          CALL TRSRT(KWORK(L(KLTRC1(ILEV))),KWORK(L(KLTRC2(ILEV))),
     *               INUM,KNEL(ILEV),DWORK(L(LCOREC)))
          IF (IER.NE.0) GOTO 99999
C
          CALL ZDISP (0,LCOREC,'DCOREC')
          IF (IER.NE.0) GOTO 99999
         ENDIF
C
C=======================================================================
C     coordinate sorting for ISORTP = 3
C=======================================================================
C
         IF (ISORTP.EQ.3) THEN
          CALL ZNEW(KNC(ILEV),-3,LCOLH,'KCOLH ')
          IF (IER.NE.0) GOTO 99999
C
          CALL ZCPY(KLCOLC(ILEV),'KCOLC ',LCOLH,'KCOLH ')
          IF (IER.NE.0) GOTO 99999
C
          NDEGP=5
          CALL ZNEW(NDEGP,-3,LDEGP,'KDEGP ')
          IF (IER.NE.0) GOTO 99999
C
          CALL CUTCE0(KWORK(L(KLLDC(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *                KWORK(L(LCOLH)),KWORK(L(LDEGP)),KNEL(ILEV),NDEGP)
C
          CALL CUTCE1(KWORK(L(KLLDC(ILEV))),KWORK(L(LCOLH)),KNEL(ILEV),
     *                KWORK(L(KLTRC1(ILEV))),KWORK(L(KLTRC2(ILEV))))
C
          CALL ZDISP(0,LDEGP,'KDEGP ')
          CALL ZDISP(0,LCOLH,'KCOLH ')
          IF (IER.NE.0) GOTO 99999
         ENDIF
C
C=======================================================================
C     matrix sorting for ISORTP > 0
C=======================================================================
C
         CALL ZNEW(KNC(ILEV)   ,-1,LCH  ,'DCH   ')
         CALL ZNEW(KNC(ILEV)   ,-3,LCOLH,'KCOLH ')
         CALL ZNEW(KNEL(ILEV)+1,-3,LLDH ,'KLDH  ')
         IF (IER.NE.0) GOTO 99999
C
         CALL ZCPY(KLC(ILEV)   ,'CC    ',LCH,  'DCH   ')
         CALL ZCPY(KLCOLC(ILEV),'KCOLC ',LCOLH,'KCOLH ')
         CALL ZCPY(KLLDC(ILEV) ,'KLDC  ',LLDH, 'KLDH  ')
         IF (IER.NE.0) GOTO 99999
C
         CALL MTSRTD(DWORK(L(KLC(ILEV)))   ,DWORK(L(LCH)),
     *               KWORK(L(KLCOLC(ILEV))),KWORK(L(LCOLH)),
     *               KWORK(L(KLLDC(ILEV))) ,KWORK(L(LLDH)),
     *               KWORK(L(KLTRC1(ILEV))),KWORK(L(KLTRC2(ILEV))),
     *               KNEL(ILEV))
C
         CALL ZDISP(0,LLDH ,'KLDH  ')
         CALL ZDISP(0,LCOLH,'KCOLH ')
         CALL ZDISP(0,LCH  ,'DCH   ')
         IF (IER.NE.0) GOTO 99999
C
        ENDIF
C
C=======================================================================
C     ILU(C) matrix
C=======================================================================
C
        IF ((ISMP.EQ.4).OR.(ISLP.EQ.3).OR.(ISLP.EQ.4)) THEN
         CALL ZNEW(KNC(ILEV),-1,KLCILU(ILEV),'DDCILU')
         IF (IER.NE.0) GOTO 99999
         CALL ZCPY(KLC(ILEV),'CC    ',KLCILU(ILEV),'DDCILU')
         IF (IER.NE.0) GOTO 99999
C
         TOLILU=1D-12
         ALPILU=0D0
         INDILU=1
         CALL IFD17(DWORK(L(KLCILU(ILEV))),KWORK(L(KLCOLC(ILEV))),
     *              KWORK(L(KLLDC(ILEV))),NEL,INDILU,ALPILU,TOLILU)
         IF (IER.NE.0) GOTO 99999
C
         CALL ZCTYPE(2,KLCILU(ILEV),'DDCILU')
         IF (IER.NE.0) GOTO 99999
        ENDIF
C
       ENDIF
C
C=======================================================================
C     grid sorting for U
C=======================================================================
C
       IF (ISORTU.GT.0) THEN
        KLTRA1(ILEV)=0
        KLTRA2(ILEV)=0
        CALL ZNEW(KNU(ILEV), 3,KLTRA1(ILEV),'KKTRA1')
        CALL ZNEW(KNU(ILEV), 3,KLTRA2(ILEV),'KKTRA2')
        IF (IER.NE.0) GOTO 99999
C
C=======================================================================
C     coordinate sorting ISORTU = 1 or 2
C=======================================================================
C
        IF ((ISORTU.EQ.1).OR.(ISORTU.EQ.2)) THEN
         LCORMC=0
         CALL ZNEW(2*KNU(ILEV),-1,LCORMC,'DCORMC')
         IF (IER.NE.0) GOTO 99999
C
         CALL TRCORM(KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *               DWORK(L(LCORMC)),KNEL(ILEV),KNVT(ILEV))
C
         INUM=ISORTU
         CALL TRSRT(KWORK(L(KLTRA1(ILEV))),KWORK(L(KLTRA2(ILEV))),
     *              INUM,KNU(ILEV),DWORK(L(LCORMC)))
         IF (IER.NE.0) GOTO 99999
C
         CALL ZDISP (0,LCORMC,'DCORMC')
         IF (IER.NE.0) GOTO 99999
        ENDIF
C
C=======================================================================
C     coordinate sorting ISORTU = 3
C=======================================================================
C
        IF (ISORTU.EQ.3) THEN
         CALL ZNEW(KNA(ILEV),-3,LCOLH,'KCOLH ')
         IF (IER.NE.0) GOTO 99999
C
         CALL ZCPY(KLCOLA(ILEV),'KCOLA ',LCOLH,'KCOLH ')
         IF (IER.NE.0) GOTO 99999
C
         NDEGU=7
         CALL ZNEW(NDEGU,-3,LDEGU,'KDEGU ')
         IF (IER.NE.0) GOTO 99999
C
         CALL CUTCE0(KWORK(L(KLLDA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *               KWORK(L(LCOLH)),KWORK(L(LDEGU)),KNU(ILEV),NDEGU)
C
         CALL CUTCE1(KWORK(L(KLLDA(ILEV))),KWORK(L(LCOLH)),KNU(ILEV),
     *               KWORK(L(KLTRA1(ILEV))),KWORK(L(KLTRA2(ILEV))))
C
         CALL ZDISP(0,LDEGU,'KDEGU ')
         CALL ZDISP(0,LCOLH,'KCOLH ')
         IF (IER.NE.0) GOTO 99999
        ENDIF
C
       ENDIF
C
       CALL ZTIME(TTT1)
       TTILU=TTILU+TTT1-TTT0
C
C=======================================================================
C     matrix restructuring
C=======================================================================
C
       CALL ZTIME(TTT0)
C
       IF (IPRECA.EQ.0) THEN
        CALL ZCTYPE(2,KLST(1,ILEV),CARRST)
        IF (IER.NE.0) GOTO 99999
C
        IF (IMASS.EQ.1) THEN
         CALL ZCTYPE(2,KLMASS(1,ILEV),CARRM)
         IF (IER.NE.0) GOTO 99999
        ENDIF
       ENDIF
C
C
       IF (IPRECA.EQ.2) THEN
        CALL ZCTYPE(2,KLST(1,ILEV),CARRST)
        IF (IER.NE.0) GOTO 99999
        CALL  OF0 (59,CFILST(ILEV),0)
        CFILE='STMAT '
        CALL  OWA2 (VWORK(L(KLST(1,ILEV))),CFILE,NA,59,0)
        REWIND(59)
        CLOSE (59)
        IF (IER.NE.0) GOTO 99999
        CALL ZDISP (0,KLST(1,ILEV),'STMAT ')
        IF (IER.NE.0) GOTO 99999
C
        IF (IMASS.EQ.1) THEN
         CALL ZCTYPE(2,KLMASS(1,ILEV),CARRM)
         IF (IER.NE.0) GOTO 99999
         CALL  OF0 (59,CFILM(ILEV),0)
         CFILE='MASMAT'
         CALL  OWA2V (VWORK(L(KLMASS(1,ILEV))),CFILE,NA,59,0)
         REWIND(59)
         CLOSE (59)
         IF (IER.NE.0) GOTO 99999
         CALL ZDISP (0,KLMASS(1,ILEV),'MASMAT')
         IF (IER.NE.0) GOTO 99999
        ENDIF
       ENDIF
C
C
       IF (IPRECA.EQ.3) THEN
        CALL  OF0 (59,CFILST(ILEV),0)
        CFILE='STMAT '
        CALL  OWA1 (DWORK(L(KLST(1,ILEV))),CFILE,NA,59,0)
        REWIND(59)
        CLOSE (59)
        IF (IER.NE.0) GOTO 99999
        CALL ZDISP (0,KLST(1,ILEV),'STMAT ')
        IF (IER.NE.0) GOTO 99999
C
        IF (IMASS.EQ.1) THEN
         CALL  OF0 (59,CFILM(ILEV),0)
         CFILE='MASMAT'
         CALL  OWA2V (VWORK(L(KLMASS(1,ILEV))),CFILE,NA,59,0)
         REWIND(59)
         CLOSE (59)
         IF (IER.NE.0) GOTO 99999
         CALL ZDISP (0,KLMASS(1,ILEV),'MASMAT')
         IF (IER.NE.0) GOTO 99999
        ENDIF
       ENDIF
C
C
       IF (IPRECB.EQ.0) THEN
        CALL ZCTYPE(2,KLB(1,ILEV),CARRDB)
        IF (IER.NE.0) GOTO 99999
        CALL ZCTYPE(2,KLB(2,ILEV),CARRDB)
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (IPRECB.EQ.2) THEN
        CALL ZDISP (0,KLB(1,ILEV) ,ARRDB(1))
        CALL ZDISP (0,KLB(2,ILEV) ,ARRDB(2))
        CALL ZDISP (0,KLCOLB(ILEV),'KCOLB ')
        CALL ZDISP (0,KLLDB (ILEV),'KLDB  ')
        IF (IER.NE.0) GOTO 99999
        KLB1(ILEV)=0
        KLB2(ILEV)=0
       ENDIF
C
       IF (IPRECB.EQ.3) THEN
        CALL ZCTYPE(2,KLB(1,ILEV),CARRDB)
        IF (IER.NE.0) GOTO 99999
        CALL ZCTYPE(2,KLB(2,ILEV),CARRDB)
        IF (IER.NE.0) GOTO 99999
        IPRECB=0
       ENDIF
C
       IF (IPRECB.EQ.4) THEN
        IPRECB=1
       ENDIF
C
       CALL ZTIME(TTT1)
       TTADF=TTADF+TTT1-TTT0
C
15    CONTINUE
C
      IF (MSHOW.GE.2) WRITE(MTERM,*)
      IF (MSHOW.GE.0) WRITE(MFILE,*)
C
C=======================================================================
C    Allocation of:  - solution vector with boundary conditions on NLMAX
C                    - RHS on NLMAX and auxiliary vectors
C                    - Auxiliary vectors 
C=======================================================================
C
      DO 30  ILEV=NLMIN,NLMAX
C
      CALL ZTIME(TTT0)
      ISETLV=1
      CALL  SETLEV (ISETLV)
C
C=======================================================================
C *** Allocation of solution and defect vectors and right hand side 
C=======================================================================
C
      NEQ=KNEQ(ILEV)
      NUP=2*NEQ+NEL
      NUPT=NUP+NEQ
C
      CALL ZNEW(NUPT,1,LUP,'DU12PT')
      IF (IER.NE.0) GOTO 99999
      KLUP(ILEV)=LUP
C
      IF (ILEV.EQ.NLMAX)  THEN
        IF ((OMGMIN.GT.0D0).OR.(OMGMAX.GT.0D0)) THEN
         CALL ZNEW(NEQ,1,LU1OLD,'DU1OLD')
         IF (IER.NE.0) GOTO 99999
         CALL ZNEW(NEQ,1,LU2OLD,'DU2OLD')
         IF (IER.NE.0) GOTO 99999
        ENDIF
C
        CALL ZNEW(NEL,1,LPOLD,'DPOLD ')
        IF (IER.NE.0) GOTO 99999
C
        CALL ZNEW(NEQ,1,LTOLD,'DTOLD ')
        IF (IER.NE.0) GOTO 99999
C
        CALL ZNEW(NEQ,1,LD1,'DD1   ')
        IF (IER.NE.0) GOTO 99999
        CALL ZNEW(NEQ,1,LD2,'DD2   ')
        IF (IER.NE.0) GOTO 99999
        CALL ZNEW(NEQ,1,LDT,'DDT   ')
        IF (IER.NE.0) GOTO 99999
        CALL ZNEW(NEL,1,LDP,'DDP   ')
        IF (IER.NE.0) GOTO 99999
      ENDIF
C
      CALL ZNEW(NUPT,1,LF12P,'DF12P ')
      IF (IER.NE.0) GOTO 99999
      KLF12P(ILEV)=LF12P
C
C=======================================================================
C *** Allocation of iteration matrix A 
C=======================================================================
C
      NA=KNA(ILEV)
      CALL  ZNEW (NA,2,LA1,'VA    ')
      IF (IER.NE.0) GOTO 99999
      KLA(ILEV)=LA1
C
C=======================================================================
C *** calculation of a vector with the areas of all finite elements
C=======================================================================
C
      CALL ZNEW(NEL+1,2,LAREA,'VAREA ')
      IF (IER.NE.0) GOTO 99999
      KLAREA(ILEV)=LAREA
      CALL SETARE(VWORK(L(LAREA)),NEL,KWORK(L(LVERT)),DWORK(L(LCORVG))) 
C
C=======================================================================
C *** calculation of rhs and start vectors on finest level only
C=======================================================================
C
      IF (ILEV.EQ.NLMAX)  THEN
C
      LF(1)=0
      LF(2)=0
      LF(3)=0
      ICLRF=1
      IF (IELT.EQ.0) 
     * CALL  XVB0 (LF,NEQ,NBLOCF,ICLRF,E031,
     *             RHS,BCONF,KF,KFN,ICUBF,ARRDF,BSNGLF)
      IF (IELT.EQ.1) 
     * CALL  XVB0 (LF,NEQ,NBLOCF,ICLRF,E030,
     *             RHS,BCONF,KF,KFN,ICUBF,ARRDF,BSNGLF)
      IF (IELT.EQ.2) 
     * CALL  XVBM0(LF,NEQ,NBLOCF,ICLRF,EM31,
     *             RHS,BCONF,KF,KFN,ICUBF,ARRDF,BSNGLF)
      IF (IELT.EQ.3) 
     * CALL  XVBM0(LF,NEQ,NBLOCF,ICLRF,EM30,
     *             RHS,BCONF,KF,KFN,ICUBF,ARRDF,BSNGLF)
      IF (IER.NE.0) GOTO 99999
      LF1=LF(1)
      LF2=LF(2)
      LF3=LF(3)
C
C=======================================================================
C *** start vector as prolongation from level NLMAX
C=======================================================================
C
      IF (ABS(ISTART).EQ.1) THEN
       IF (ISTART.EQ.1) THEN
        IFMTS=0
       ELSE
        IFMTS=1
       ENDIF
       CFILE='DU12P '
       CALL  OF0 (MSTART,CSTART,IFMTS)
       CALL  ORA11 (DWORK(L(LUP)),CFILE,MSTART,IFMTS)
       CLOSE(MSTART)
      ENDIF
C
      IF (ABS(ISTART).EQ.3) THEN
       IF (ISTART.EQ.3) THEN
        IFMTS=0
       ELSE
        IFMTS=1
       ENDIF
       CFILE='DU12P '
       CALL  OF0 (MSTART,CSTART,IFMTS)
       CALL  ORA11 (DWORK(L(LUP)),CFILE,MSTART,IFMTS)
       CLOSE(MSTART)
       DO 8000 IU=1,NU
8000   DWORK(L(LUP)+NUP-1+IU)=0D0
      ENDIF
C

C=======================================================================
C *** start vector as prolongation from level NLMAX-1
C=======================================================================
C
      IF (ABS(ISTART).EQ.2) THEN
        I1=NLMAX-1
        KU1C=L(KLUP(I1))
        IF (ISTART.EQ.2) THEN
         IFMTS=0
        ELSE
         IFMTS=1
        ENDIF
        CFILE='DU12P '
        CALL  OF0 (MSTART,CSTART,IFMTS)
        CALL  ORA11 (DWORK(KU1C),CFILE,MSTART,IFMTS)
        CLOSE(MSTART)
        NUC=KNU(I1)
        NPC=KNP(I1)
        NVTC=KNVT(I1)
        KU2C=KU1C+NUC
        KPC=KU2C+NUC
        KVERTC=L(KLVERT(I1))
        KMIDC=L(KLMID(I1))
        KADJC=L(KLADJ(I1))
        KU1=L(LUP)
        KU2=KU1+NEQ
        KP=KU2+NEQ
C
        CALL  PROLU (DWORK(KU1C),DWORK(KU2C),DWORK(KPC),
     *               DWORK(KU1),DWORK(KU2),DWORK(KP),
     *               KWORK(KVERTC),KWORK(KMIDC),KWORK(KADJC),
     *               NUC, NPC, NVTC,
     *               KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *               NEQ, NEL, NVT)
      ENDIF
C
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
C=======================================================================
C *** Dirichlet boundary updates of solution and rhs vector
C=======================================================================
C
      CALL ZTIME(TTT0)
      KU1=L(LUP)
      KU2=KU1+NEQ
      CALL BDRSET (DWORK(KU1),DWORK(KU2),DWORK(L(LF1)),
     *             DWORK(L(LF2)),KWORK(L(KLMBD(NLEV))),
     *             DWORK(L(KLDBD(NLEV))),KWORK(L(LNPR)),
     *             KNMBD(NLEV),NVT,PARX,PARY,UE,2)
C
      CALL  LCP1 (DWORK(L(LF1)), DWORK(L(LF12P)),     NEQ)
      CALL  LCP1 (DWORK(L(LF2)), DWORK(L(LF12P)+NEQ), NEQ)
C
      KT=KU1+2*NEQ+NEL
      CALL BDSETT (DWORK(KT),DWORK(L(LF3)),
     *             KWORK(L(KLMBD(NLEV))),DWORK(L(KLDBD(NLEV))),
     *             KWORK(L(LVERT)),KWORK(L(LMID)),
     *             KWORK(L(LNPR)),DWORK(L(LCORVG)),
     *             KNMBD(NLEV),PARX,PARY,UE,0D0)
C
      CALL  LCP1 (DWORK(L(LF3)), DWORK(L(LF12P)+2*NEQ+NEL), NEQ)
      CALL ZTIME(TTT1)
      TTBDR=TTBDR+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      CALL  ZDISP (0, LF1, ARRDF(1) )
      CALL  ZDISP (0, LF2, ARRDF(2) )
      CALL  ZDISP (0, LF3, ARRDF(3) )
      IF (IER.NE.0) GOTO 99999
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
      ENDIF
C
C=======================================================================
C *** Auxiliary vector on all levels
C=======================================================================
C
      CALL ZTIME(TTT0)
      CALL ZNEW(NUPT,1,LAUX,'DAUX  ')
      IF (IER.NE.0) GOTO 99999
      KLAUX(ILEV)=LAUX
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
30    CONTINUE
C
C
      CALL ZTIME(TTT1)
      TTLIN=TTT1-TTT0L
C
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'time for initialization of linear operators : ', 
     *                TTLIN
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'time for initialization of linear operators : ', 
     *                TTLIN
      IF (MSHOW.GE.2) WRITE(MTERM,*)
      IF (MSHOW.GE.0) WRITE(MFILE,*)
C
C
C
99999 END

************************************************************************
* Open file CDATA, read the MSHOW variable from the data file CDATA.
* Set the file pointer to the beginning of the file to allow complete
* parsing of the file.
************************************************************************
      SUBROUTINE OPNDAT (MDATA, MSHOW, CDATA)
      IMPLICIT NONE
      INTEGER MDATA, MSHOW
      CHARACTER*(*) CDATA
      
      INTEGER i

      CALL OF0 (MDATA,CDATA,1)
      DO 10 i=1,20
        READ (MDATA,*)
10    CONTINUE
      READ (MDATA,*) MSHOW
      REWIND (MDATA)
      
      END

************************************************************************
      SUBROUTINE GDAT (MDATA,MSHOW,IRMESH,IPROJ)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNLEV=9)
C
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------
C *** Standard COMMON blocks
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
C
C *** COMMON blocks for multigrid data management
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTIME/ TTMG,TTS,TTE,TTD,TTP,TTR,IMTIME
      COMMON /MGILUU/ ISORTU,KLAILU(NNLEV),KLTRA1(NNLEV),KLTRA2(NNLEV)
      COMMON /MGILUP/ ISORTP,KLCILU(NNLEV),KLTRC1(NNLEV),KLTRC2(NNLEV)
C
C *** user COMMON blocks
      INTEGER  VIPARM 
      DIMENSION VIPARM(100)                     
      EQUIVALENCE (IAUSAV,VIPARM)   
      COMMON /IPARM/ IAUSAV,IELT,ISTOK,IRHS,IBDR,IERANA,
     *               IMASS,IMASSL,IUPW,IPRECA,IPRECB,
     *               ICUBML,ICUBM,ICUBA,ICUBN,ICUBB,ICUBF,
     *               INLMIN,INLMAX,ICYCU,ILMINU,ILMAXU,IINTU,
     *               ISMU,ISLU,NSMU,NSLU,NSMUFA,ICYCP,ILMINP,ILMAXP,
     *               IINTP,ISMP,ISLP,NSMP,NSLP,NSMPFA
C
      DOUBLE PRECISION VRPARM,NY
      DIMENSION VRPARM(100)
      EQUIVALENCE (NY,VRPARM)                          
      COMMON /RPARM/  NY,RE,UPSAM,OMGMIN,OMGMAX,OMGINI,
     *                EPSUR,EPSUD,DMPUD,DMPUMG,DMPUSL,RLXSMU,RLXSLU,
     *                AMINU,AMAXU,EPSP,DMPPMG,DMPPSL,RLXSMP,RLXSLP,
     *                AMINP,AMAXP
C
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
      COMMON /NSADAT/ TIMEMX,DTMIN,DTMAX,DTFACT,TIMEIN,EPSADI,EPSADL,
     *                EPSADU,IEPSAD,IADIN,IREPIT,IADTIM,PRDIF1,PRDIF2
      COMMON /NSSAV/  INSAV,INSAVN
      COMMON /NSSAVF/ DTFILM,DTFILO,DTAVS,DTAVSO,DTGMV,DTGMVO,
     *                IFUSAV,IFPSAV,IFXSAV,IAVS,IGMV,IFINIT
      COMMON /NSTIME/ TTGRID,TTPOST,TTADF,TTUPW,TTBDR,TTLC,TTILU,
     *                TTMGU,TTSU,TTEU,TTDU,TTPU,TTRU,
     *                TTMGP,TTSP,TTEP,TTDP,TTPP,TTRP
      COMMON /NSEXL/  ITEXL,LTML,TIML11,TIML12,TIML31,TIML32
C
      CHARACTER CPARM1*60,CMESH1*60,CFILE1*60,CSTART*60,CSOL*60
      COMMON /FILES/ IMESH1,MMESH1,CPARM1,CMESH1,MFILE1,CFILE1,
     *               ISTART,MSTART,CSTART,ISOL,MSOL,CSOL
C
      INCLUDE 'bouss.inc'
      SAVE 
C
C-----------------------------------------------------------------------
C *** Input file
C-----------------------------------------------------------------------
C
      READ(MDATA,*)
      READ(MDATA,*)
      READ(MDATA,*) 
      IF (MSHOW.GE.2) WRITE(MTERM,1)
   1  FORMAT(80('-'))
      IF (MSHOW.GE.2) WRITE(MTERM,*) '        INPUT DATA'
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      READ(MDATA,*)
      READ(MDATA,*)
C
      READ(MDATA,*) IMESH1
      IF (IMESH1.NE.1) IMESH1=0
      READ(MDATA,*) IRMESH
      IRMESH=ABS(IRMESH)
      READ(MDATA,*) CPARM1
      READ(MDATA,*) CMESH1
      MMESH1=61
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'Parametrization file = ',CPARM1
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'Coarse grid file     = ',CMESH1
C
      READ(MDATA,*) CFILE1
      MFILE1=62
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'MFILE,CFILE: ', MFILE1,CFILE1
      MFILE=MFILE1
C
      READ(MDATA,*) ISTART
      IF (ABS(ISTART).GT.3) ISTART=0
      READ(MDATA,*) CSTART
      MSTART=63
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*)'ISTART,MSTART,CSTART: ', ISTART,MSTART,CSTART
C
      READ(MDATA,*) ISOL
      IF (ABS(ISOL).GT.1) ISOL=1
      READ(MDATA,*) CSOL
      MSOL=64
      IF (MSHOW.GE.2) WRITE(MTERM,*)'ISOL,MSOL,CSOL: ', ISOL,MSOL,CSOL
C
C-----------------------------------------------------------------------
C *** Open file for user output
C-----------------------------------------------------------------------
C
      CALL  OF0 (MFILE1,CFILE1,1)
C
      IF (MSHOW.GE.0) WRITE(MFILE,1)
      IF (MSHOW.GE.0) WRITE(MFILE,*) '        INPUT DATA'
      IF (MSHOW.GE.0) WRITE(MFILE,1)
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'Parametrization file = ',CPARM1
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'Coarse grid file     = ',CMESH1
C
C-----------------------------------------------------------------------
C *** Values for /OUTPUT/
C-----------------------------------------------------------------------
C
      READ(MDATA,*) 
      READ(MDATA,*)
      READ(MDATA,*)
C
      READ(MDATA,*) M
      M=ABS(M)
C
      READ(MDATA,*) MT
      MT=ABS(MT)
C
      READ(MDATA,*) ICHECK
      ICHECK=ABS(ICHECK)
C
      READ(MDATA,*) MSHOW
      MSHOW=ABS(MSHOW)
C
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'message level for file output:  M = ',M
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'message level for terminal output:  MT = ',MT
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'level for tracing:  ICHECK = ',ICHECK
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'output level:  MSHOW = ',MSHOW
      IF (MSHOW.GE.2) WRITE(MTERM,1)
C
C-----------------------------------------------------------------------
C *** Values for /IPARM/,etc.
C-----------------------------------------------------------------------
C
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'Integer parameters of /IPARM/ :'
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Integer parameters of /IPARM/,etc. :'
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      WRITE(MFILE,1)
      READ(MDATA,*) 
      READ(MDATA,*)
      READ(MDATA,*)
C
      READ(MDATA,*)  NLMIN
      NLMIN=ABS(NLMIN)
      IF (NLMIN.EQ.0) NLMIN=1
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'minimum mg-level:  NLMIN = ', NLMIN
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'minimum mg-level:  NLMIN = ', NLMIN
C
      READ(MDATA,*)  NLMAX
      NLMAX=ABS(NLMAX)
      IF (NLMAX.LT.NLMIN) NLMAX=NLMIN
      NLEV=NLMAX
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'maximum mg-level:  NLMAX = ', NLMAX
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'maximum mg-level:  NLMAX = ', NLMAX
C
      READ(MDATA,*) IELT
      IF ((IELT.LT.0).OR.(IELT.GT.3)) IELT=3
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'element type   = ',IELT
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'element type   = ',IELT
C
      READ(MDATA,*)  ISTOK
      IF (ISTOK.NE.1) ISTOK=0
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Stokes calculation:  ISTOK = ', ISTOK
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Stokes calculation:  ISTOK = ', ISTOK
C
      READ(MDATA,*) IRHS
      IF ((IRHS.LT.0).OR.(IRHS.GT.1)) IRHS=1
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'RHS    generation   = ',IRHS
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'RHS    generation   = ',IRHS
C
      READ(MDATA,*) IBDR
      IF ((IBDR.LT.0).OR.(IBDR.GT.2)) IBDR=2
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'Boundary generation = ',IBDR
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'Boundary generation = ',IBDR
C
      READ(MDATA,*) IERANA
      IERANA=ABS(IERANA)
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'Error evaluation    = ',IERANA
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'Error evaluation    = ',IERANA
C
      READ(MDATA,*) IMASS
      IF (IMASS.NE.1) IMASS=0
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'mass evaluation     = ',IMASS
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'mass evaluation     = ',IMASS
C
      READ(MDATA,*) IMASSL
      IF ((IMASSL.LT.0).OR.(IMASSL.GT.3)) IMASSL=3
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'lumped mass eval.   = ',IMASSL
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'lumped mass eval.   = ',IMASSL
C
      READ(MDATA,*) IUPW
      IF ((IUPW.GT.1).OR.(IUPW.LT.0)) IUPW=1
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'convective part     = ',IUPW
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'convective part     = ',IUPW
C
      READ(MDATA,*) IPRECA
      IF ((IPRECA.LT.0).OR.(IPRECA.GT.4)) IPRECA=1
      IF ((IPRECA.EQ.4).AND.(IUPW.EQ.1))  IPRECA=1
      IF ((IPRECA.EQ.4).AND.(ISTOK.EQ.1)) IPRECA=1
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'Accuracy for ST     = ',IPRECA
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'Accuracy for ST     = ',IPRECA
C
      READ(MDATA,*) IPRECB
      IF ((IPRECB.LT.0).OR.(IPRECB.GT.4)) IPRECB=2
      IF ((IBDR.EQ.2).AND.(IPRECB.EQ.0)) IPRECB=1
      IF ((IBDR.EQ.2).AND.(IPRECB.GE.2)) IPRECB=4
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'Accuracy for B      = ',IPRECB
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'Accuracy for B      = ',IPRECB
C
      READ(MDATA,*) ICUBML
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'ICUB lumped mass matrix = ',ICUBML
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'ICUB lumped mass matrix = ',ICUBML
C
      READ(MDATA,*) ICUBM
      ICUBM=ABS(ICUBM)
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'ICUB mass matrix        = ',ICUBM
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'ICUB mass matrix        = ',ICUBM
C
      READ(MDATA,*) ICUBA
      ICUBA=ABS(ICUBA)
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'ICUB diff. matrix       = ',ICUBA
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'ICUB diff. matrix       = ',ICUBA
C
      READ(MDATA,*) ICUBN
      ICUBN=ABS(ICUBN)
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'ICUB conv. matrix       = ',ICUBN
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'ICUB conv. matrix       = ',ICUBN
C
      READ(MDATA,*) ICUBB
      ICUBB=ABS(ICUBB)
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'ICUB matrices B1,B2     = ',ICUBB
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'ICUB matrices B1,B2     = ',ICUBB
C
      READ(MDATA,*) ICUBF
      ICUBF=ABS(ICUBF)
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'ICUB right hand side    = ',ICUBF
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'ICUB right hand side    = ',ICUBF
C
      READ(MDATA,*) INLMIN
      IF (INLMIN.LT.-1) INLMIN=-1
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*)'minimum of nonlinear iterations: INLMIN = ',
     *  INLMIN
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*)'minimum of nonlinear iterations: INLMIN = ',
     *  INLMIN
C
      READ(MDATA,*) INLMAX
      IF (INLMAX.LT.-1) INLMAX=-1
      IF (INLMAX.LT.INLMIN) INLMAX=INLMIN
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*)'maximum of nonlinear iterations: INLMAX = ',
     *  INLMAX
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*)'maximum of nonlinear iterations: INLMAX = ',
     *  INLMAX
C
      ITEXL=0
      IF ((INLMIN.EQ.INLMAX).AND.(INLMIN.EQ. 1)) THEN
       ITEXL=1
       TIML31=1D0
       TIML32=1D0
       TIML11=1D0
       TIML12=1D0
      ENDIF
      IF ((INLMIN.EQ.INLMAX).AND.(INLMIN.EQ.-1)) THEN
       INLMIN=1
       INLMAX=1
      ENDIF
C
C
      READ(MDATA,*)  ISORTU
      ISORTU=ABS(ISORTU)
      IF (ISORTU.GT.3) ISORTU=3
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'matrix resorting:  ISORTU = ', ISORTU
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'matrix resorting:  ISORTU = ', ISORTU
C
      READ(MDATA,*)  ICYCU
      ICYCU=ABS(ICYCU)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'type of mg-cycle:  ICYCU = ', ICYCU
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'type of mg-cycle:  ICYCU = ', ICYCU
C
      READ(MDATA,*) ILMINU
      ILMINU=ABS(ILMINU)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'minimum of linear mg steps :  ILMINU = ', ILMINU
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'minimum of linear mg steps :  ILMINU = ', ILMINU
C
      READ(MDATA,*) ILMAXU
      ILMAXU=ABS(ILMAXU)
      IF (ILMAXU.LT.ILMINU) ILMAXU=ILMINU
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'maximum of linear mg steps :  ILMAXU = ', ILMAXU
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'maximum of linear mg steps :  ILMAXU = ', ILMAXU
C
      READ(MDATA,*)  IINTU
      IF (IINTU.NE.2) IINTU=1
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'type of interpolation:  IINTU = ', IINTU
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'type of interpolation:  IINTU = ', IINTU
C
      READ(MDATA,*) ISMU
      IF ((ISMU.LT.1).OR.(ISMU.GT.4)) ISMU=2
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'type of smoother :  ISMU = ',ISMU
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'type of smoother :  ISMU = ',ISMU
C
      READ(MDATA,*) ISLU
      IF ((ISLU.LT.1).OR.(ISLU.GT.4)) ISLU=1
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'type of solver :  ISLU = ',ISLU
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'type of solver :  ISLU = ',ISLU
C
      READ(MDATA,*) NSMU
      NSMU=ABS(NSMU)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'number of smoothing steps :  NSMU = ', NSMU
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'number of smoothing steps :  NSMU = ', NSMU
C
      READ(MDATA,*) NSLU
      NSLU=ABS(NSLU)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'maximum of solver-iterations :  NSLU = ',NSLU
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'maximum of solver-iterations :  NSLU = ',NSLU
C
      READ(MDATA,*) NSMUFA
      IF (NSMUFA.LT.1) NSMUFA=1
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'factor sm. steps on coarser lev.:NSMUFA=',NSMUFA
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'factor sm. steps on coarser lev.:NSMUFA=',NSMUFA
C
C
      READ(MDATA,*)  ISORTP
      ISORTP=ABS(ISORTP)
      IF (ISORTP.GT.3) ISORTP=3
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'matrix resorting:  ISORTP = ', ISORTP
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'matrix resorting:  ISORTP = ', ISORTP
C
      READ(MDATA,*)  ICYCP
      ICYCP=ABS(ICYCP)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'type of mg-cycle:  ICYCP = ', ICYCP
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'type of mg-cycle:  ICYCP = ', ICYCP
C
      READ(MDATA,*) ILMINP
      ILMINP=ABS(ILMINP)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'minimum of linear mg steps :  ILMINP = ', ILMINP
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'minimum of linear mg steps :  ILMINP = ', ILMINP
C
      READ(MDATA,*) ILMAXP
      ILMAXP=ABS(ILMAXP)
      IF (ILMAXP.LT.ILMINP) ILMAXP=ILMINP
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'maximum of linear mg steps :  ILMAXP = ', ILMAXP
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'maximum of linear mg steps :  ILMAXP = ', ILMAXP
C
      READ(MDATA,*)  IINTP
      IF ((IINTP.LT.1).OR.(IINTP.GT.4)) IINTP=4
      IF (MSHOW.GE.2)  
     * WRITE(MTERM,*) 'type of interpolation:  IINTP = ', IINTP
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'type of interpolation:  IINTP = ', IINTP
C
      READ(MDATA,*) ISMP
      IF ((ISMP.LT.1).OR.(ISMP.GT.4)) ISMP=2
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'type of smoother :  ISMP = ',ISMP
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'type of smoother :  ISMP = ',ISMP
C
      READ(MDATA,*) ISLP
      IF ((ISLP.LT.1).OR.(ISLP.GT.4)) ISLP=1
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'type of solver :  ISLP = ',ISLP
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'type of solver :  ISLP = ',ISLP
C
      READ(MDATA,*) NSMP
      NSMP=ABS(NSMP)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'number of smoothing steps :  NSMP = ', NSMP
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'number of smoothing steps :  NSMP = ', NSMP
C
      READ(MDATA,*) NSLP
      NSLP=ABS(NSLP)
      IF (MSHOW.GE.2)  
     * WRITE(MTERM,*) 'maximum of solver-iterations :  NSLP = ',NSLP
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'maximum of solver-iterations :  NSLP = ',NSLP
C
      READ(MDATA,*) NSMPFA
      IF (NSMPFA.LT.1) NSMPFA=1
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'factor sm. steps on coarser lev.:NSMPFA=',NSMPFA
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'factor sm. steps on coarser lev.:NSMPFA=',NSMPFA
C
      READ(MDATA,*) ILMINT
      ILMINT=ABS(ILMINT)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'minimum of linear mg steps :  ILMINT = ', ILMINT
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'minimum of linear mg steps :  ILMINT = ', ILMINT
C
      READ(MDATA,*) ILMAXT
      ILMAXT=ABS(ILMAXT)
      IF (ILMAXT.LT.ILMINT) ILMAXT=ILMINT
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'maximum of linear mg steps :  ILMAXT = ', ILMAXT
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'maximum of linear mg steps :  ILMAXT = ', ILMAXT
C
      READ(MDATA,*)  IINTT
      IF ((IINTT.LT.1).OR.(IINTT.GT.4)) IINTT=4
      IF (MSHOW.GE.2)  
     * WRITE(MTERM,*) 'type of interpolation:  IINTT = ', IINTT
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'type of interpolation:  IINTT = ', IINTT
C
      READ(MDATA,*) ISMT
      IF ((ISMT.LT.1).OR.(ISMT.GT.4)) ISMT=2
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'type of smoother :  ISMT = ',ISMT
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'type of smoother :  ISMT = ',ISMT
C
      READ(MDATA,*) ISLT
      IF ((ISLT.LT.1).OR.(ISLT.GT.4)) ISLT=1
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'type of solver :  ISLT = ',ISLT
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'type of solver :  ISLT = ',ISLT
C
      READ(MDATA,*) NSMT
      NSMT=ABS(NSMT)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'number of smoothing steps :  NSMT = ', NSMT
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'number of smoothing steps :  NSMT = ', NSMT
C
      READ(MDATA,*) NSLT
      NSLT=ABS(NSLT)
      IF (MSHOW.GE.2)  
     * WRITE(MTERM,*) 'maximum of solver-iterations :  NSLT = ',NSLT
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'maximum of solver-iterations :  NSLT = ',NSLT
C
C-----------------------------------------------------------------------
C *** Values for /RPARM/,etc.
C-----------------------------------------------------------------------
C
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'Real parameters of /RPARM/,etc. :'
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'Real parameters of /RPARM/,etc. :'
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.0) WRITE(MFILE,1)
      READ(MDATA,*) 
      READ(MDATA,*)
      READ(MDATA,*)
C
      READ(MDATA,*) RE
      RE=ABS(RE)
      IF (RE.LT.1D-8) RE=1D-8
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'Viscosity parameter:  1/NU = ', RE
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'Viscosity parameter:  1/NU = ', RE
C
      NY=1.D0/RE
C
      READ(MDATA,*) ALFA
      IF (MSHOW.GE.2) WRITE(MTERM,*) 
     *               'Temperature parameter:  Alfa = ', ALFA
      IF (MSHOW.GE.0) WRITE(MFILE,*) 
     *               'Temperature parameter:  Alfa = ', ALFA
C
      READ(MDATA,*) VOLEXK
      IF (MSHOW.GE.2) WRITE(MTERM,*) 
     *               ' parameter:  Volexk = ', VOLEXK
      IF (MSHOW.GE.0) WRITE(MFILE,*) 
     *               ' parameter:  Volexk = ', VOLEXK
C
      READ(MDATA,*) G(1)
      READ(MDATA,*) G(2)
c
      READ(MDATA,*) UPSAM
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*)'parameter for Samarskij-upwind:  UPSAM = ', UPSAM
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*)'parameter for Samarskij-upwind:  UPSAM = ', UPSAM
C
      READ(MDATA,*) OMGMIN
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'lower limit for optimal OMEGA: OMGMIN = ', OMGMIN
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'lower limit for optimal OMEGA: OMGMIN = ', OMGMIN
C
      READ(MDATA,*) OMGMAX
      IF (OMGMAX.LT.OMGMIN) OMGMAX=OMGMIN
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'upper limit for optimal OMEGA: OMGMAX = ', OMGMAX
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'upper limit for optimal OMEGA: OMGMAX = ', OMGMAX
C
      READ(MDATA,*) OMGINI
      IF (OMGINI.LT.OMGMIN) OMGINI=OMGMIN
      IF (OMGINI.GT.OMGMAX) OMGINI=OMGMAX
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'start value for optimal OMEGA: OMGINI = ', OMGINI
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'start value for optimal OMEGA: OMGINI = ', OMGINI
C
      READ(MDATA,*) OMEGAT
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'fixed OMEGA for Temp: OMEGAT = ', OMEGAT
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'fixed OMEGA for Temp: OMEGAT = ', OMEGAT
C
      READ(MDATA,*) EPSUR
      EPSUR=ABS(EPSUR)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'limit for U-changes :          EPSUR  = ', EPSUR
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'limit for U-changes :          EPSUR  = ', EPSUR
C
      READ(MDATA,*) EPSUD
      EPSUD=ABS(EPSUD)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'limit for U-defects   :        EPSUD  = ', EPSUD
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'limit for U-defects   :        EPSUD  = ', EPSUD
C
      READ(MDATA,*) DMPUD
      DMPUD=ABS(DMPUD)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'damping of U-defects :         DMPUD  = ', DMPUD
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'damping of U-defects :         DMPUD  = ', DMPUD
C
      READ(MDATA,*) DMPUMG
      DMPUMG=ABS(DMPUMG)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'damping of U-MG residuals :    DMPUMG = ', DMPUMG
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'damping of U-MG residuals :    DMPUMG = ', DMPUMG
C
      READ(MDATA,*) DMPUSL
      DMPUSL=ABS(DMPUSL)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'damping of solver residuals :  DMPUSL = ', DMPUSL
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'damping of solver residuals :  DMPUSL = ', DMPUSL
C
      READ(MDATA,*) RLXSMU
      RLXSMU=ABS(RLXSMU)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'relaxation for the U-smoother: RLXSMU = ', RLXSMU
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'relaxation for the U-smoother: RLXSMU = ', RLXSMU
C
      READ(MDATA,*) RLXSLU
      RLXSLU=ABS(RLXSLU)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'relaxation for the U-solver :  RLXSLU = ', RLXSLU
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'relaxation for the U-solver :  RLXSLU = ', RLXSLU
C
      READ(MDATA,*) AMINU
      IF (MSHOW.GE.2)  
     * WRITE(MTERM,*)'lower limit for opt. U-ALPHA:   AMINU  = ', AMINU
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*)'lower limit for opt. U-ALPHA:   AMINU  = ', AMINU
C
      READ(MDATA,*) AMAXU
      IF (AMAXU.LT.AMINU) AMAXU=AMINU
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*)'upper limit for opt. U-ALPHA:   AMAXU  = ', AMAXU
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*)'upper limit for opt. U-ALPHA:   AMAXU  = ', AMAXU
C
C
      READ(MDATA,*) EPSP
      EPSP=ABS(EPSP)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'error reduction for p :        EPSP   = ', EPSP 
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'error reduction for p :        EPSP   = ', EPSP 
C
      READ(MDATA,*) DMPPMG
      DMPPMG=ABS(DMPPMG)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'damping of P-MG residuals :    DMPPMG = ', DMPPMG
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'damping of P-MG residuals :    DMPPMG = ', DMPPMG
C
      READ(MDATA,*) DMPPSL
      DMPPSL=ABS(DMPPSL)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'damping of P-solver residuals: DMPPSL = ', DMPPSL
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'damping of P_solver residuals: DMPPSL = ', DMPPSL
C
      READ(MDATA,*) RLXSMP
      RLXSMP=ABS(RLXSMP)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'relaxation for the P-smoother: RLXSMP = ', RLXSMP
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'relaxation for the P-smoother: RLXSMP = ', RLXSMP
C
      READ(MDATA,*) RLXSLP
      RLXSLP=ABS(RLXSLP)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'relaxation for the P-solver :  RLXSLP = ', RLXSLP
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'relaxation for the P-solver :  RLXSLP = ', RLXSLP
C
      READ(MDATA,*) AMINP
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*)'lower limit for opt. P-ALPHA:   AMINP  = ', AMINP
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*)'lower limit for opt. P-ALPHA:   AMINP  = ', AMINP
C
      READ(MDATA,*) AMAXP
      IF (AMAXP.LT.AMINP) AMAXP=AMINP
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*)'upper limit for opt. P-ALPHA:   AMAXP  = ', AMAXP
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*)'upper limit for opt. P-ALPHA:   AMAXP  = ', AMAXP
C
      READ(MDATA,*) DMPTMP
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*)'damping of T-residuals:   DMPTMP  = ', DMPTMP
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*)'damping of T-residuals:   DMPTMP  = ', DMPTMP
C
C-----------------------------------------------------------------------
C *** Values for /NS.../
C-----------------------------------------------------------------------
C
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'Parameters of /NS.../ :'
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'Parameters of /NS.../ :'
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.0) WRITE(MFILE,1)
      READ(MDATA,*) 
      READ(MDATA,*)
      READ(MDATA,*)
C
      READ(MDATA,*) IPROJ
      IF (IPROJ.GT.1) IPROJ=1
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'projection type          : IPROJ  = ', IPROJ
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'projection type          : IPROJ  = ', IPROJ
C
      READ(MDATA,*) NITNS
      NITNS=ABS(NITNS)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Number of time steps     : NITNS  = ', NITNS
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Number of time steps     : NITNS  = ', NITNS
C
      READ(MDATA,*) EPSNS
      EPSNS=ABS(EPSNS)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'limit for time derivative: EPSNS  = ', EPSNS
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'limit for time derivative: EPSNS  = ', EPSNS
C
      READ(MDATA,*) TIMENS
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Total time               : TIMENS = ', TIMENS
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Total time               : TIMENS = ', TIMENS
C
      READ(MDATA,*) THETA
      IF (THETA.LT.0D0) THETA=0D0
      IF (THETA.GT.1D0) THETA=1D0
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Theta                    : THETA  = ', THETA
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Theta                    : THETA  = ', THETA
C
      READ(MDATA,*) TSTEP
      EPSNS=ABS(EPSNS)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Time step                : TSTEP  = ', TSTEP
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Time step                : TSTEP  = ', TSTEP
C
      READ(MDATA,*) IFRSTP
      IF (IFRSTP.NE.1) IFRSTP=0
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Fractional step          : IFRSTP = ', IFRSTP
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Fractional step          : IFRSTP = ', IFRSTP
C
      READ(MDATA,*) INSAV
      IF (INSAV.LT.0) INSAV=0
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Stepsize for nonsteady savings: INSAV = ', INSAV
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Stepsize for nonsteady savings: INSAV = ', INSAV
C
      READ(MDATA,*) INSAVN
      IF (INSAVN.LT. 0) INSAVN=0
      IF (INSAVN.GT.10) INSAVN=10
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Number of files               : INSAVN = ',
     *  INSAVN
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Number of files               : INSAVN = ',
     *  INSAVN
C
      READ(MDATA,*) DTFILM
      DTFILM=ABS(DTFILM)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Time step for Film            : DTFILM = ',
     *  DTFILM
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Time step for Film            : DTFILM = ',
     *  DTFILM
      DTFILO=TIMENS
C
      READ(MDATA,*) DTAVS
      DTAVS=ABS(DTAVS)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Time step for AVS             : DTAVS = ',
     *  DTAVS
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Time step for AVS             : DTAVS = ',
     *  DTAVS
      DTAVSO=TIMENS
C
      READ(MDATA,*) DTGMV
      DTGMV=ABS(DTGMV)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Time step for GMV             : DTGMV = ',
     *  DTGMV
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Time step for GMV             : DTGMV = ',
     *  DTGMV
      DTGMVO=TIMENS
C
      READ(MDATA,*) IFUSAV
      IFUSAV=MIN(IFUSAV,NLMAX)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Level for velocity            : IFUSAV = ', 
     * IFUSAV
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Level for velocity            : IFUSAV = ',
     *  IFUSAV
C
      READ(MDATA,*) IFPSAV
      IFPSAV=MIN(IFPSAV,NLMAX)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Level for pressure            : IFPSAV = ',
     *  IFPSAV
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Level for pressure            : IFPSAV = ',
     *  IFPSAV
C
      READ(MDATA,*) IFXSAV
      IFXSAV=MIN(IFXSAV,NLMAX)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Level for streamlines         : IFXSAV = ',
     *  IFXSAV
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Level for streamlines         : IFXSAV = ',
     *  IFXSAV
C
      READ(MDATA,*) IAVS
      IAVS=MIN(IAVS,NLMAX)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Level for AVS                 : IAVS = ',
     *  IAVS
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Level for AVS                 : IAVS = ',
     *  IAVS
C
      READ(MDATA,*) IGMV
      IGMV=MIN(IGMV,NLMAX)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Level for GMV                 : IGMV = ',
     *  IGMV
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Level for GMV                 : IGMV = ',
     *  IGMV
C
      READ(MDATA,*) IFINIT
      IFINIT=ABS(IFINIT)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Start file                    : IFINIT = ',
     *  IFINIT
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Start file                    : IFINIT = ',
     *  IFINIT
C
      READ(MDATA,*) IADTIM
      IF (ABS(IADTIM).GT.3) IADTIM=0
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Type of adaptivity            : IADTIM = ',
     *  IADTIM
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Type of adaptivity            : IADTIM = ',
     *  IADTIM
C
      READ(MDATA,*) TIMEMX
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Max. Time                     : TIMEMX = ',
     *  TIMEMX
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Max. Time                     : TIMEMX = ',
     *  TIMEMX
C
      READ(MDATA,*) DTMIN
      DTMIN=ABS(DTMIN)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Min. Timestep                 : DTMIN  = ',
     *  DTMIN
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Min. Timestep                 : DTMIN  = ',
     *  DTMIN
C
      READ(MDATA,*) DTMAX
      DTMAX=ABS(DTMAX)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Max. Timestep                 : DTMAX  = ',
     *  DTMAX
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Max. Timestep                 : DTMAX  = ',
     *  DTMAX
C
      READ(MDATA,*) DTFACT
      DTFACT=ABS(DTFACT)
      IF (DTFACT.LT.1D0) DTFACT=1D0
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Max. Timestep change          : DTFACT = ',
     *  DTFACT
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Max. Timestep change          : DTFACT = ',
     *  DTFACT
C
      READ(MDATA,*) TIMEIN
      TIMEIN=ABS(TIMEIN)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Time for start procedure      : TIMEIN = ',
     *  TIMEIN
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Time for start procedure      : TIMEIN = ',
     *  TIMEIN
C
      READ(MDATA,*) EPSADI
      EPSADI=ABS(EPSADI)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'EPS for start procedure       : EPSADI = ',
     *  EPSADI
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'EPS for start procedure       : EPSADI = ',
     *  EPSADI
C
      READ(MDATA,*) EPSADL
      EPSADL=ABS(EPSADL)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'EPS for acceptance            : EPSADL = ',
     *  EPSADL
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'EPS for acceptance            : EPSADL = ',
     *  EPSADL
C
      READ(MDATA,*) EPSADU
      EPSADU=ABS(EPSADU)
      IF (EPSADU.GT.1D0) EPSADU=1D0
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'EPS for not acceptance        : EPSADU = ',
     *  EPSADU
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'EPS for not acceptance        : EPSADU = ',
     *  EPSADU
C
      READ(MDATA,*) IEPSAD
      IF ((IEPSAD.LT.0).OR.(IEPSAD.GT.8)) IEPSAD=1
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Acceptance criterion          : IEPSAD = ',
     *  IEPSAD
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Acceptance criterion          : IEPSAD = ',
     *  IEPSAD
C
      READ(MDATA,*) IADIN
      IF ((IADIN.LT.0).OR.(IADIN.GT.2)) IADIN=0
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Start procedure               : IADIN  = ',
     *  IADIN
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Start procedure               : IADIN  = ',
     *  IADIN
C
      READ(MDATA,*) IREPIT
      IF (IREPIT.LT.1) IREPIT=1
      IF (IREPIT.GT.9) IREPIT=9
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Max.numbers of repetitions    : IREPIT = ',
     *  IREPIT
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Max.numbers of repetitions    : IREPIT = ',
     *  IREPIT
C
      READ(MDATA,*) PRDIF1
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Parameter for reactive prec.  : PRDIF1 = ',
     *  PRDIF1
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Parameter for reactive prec.  : PRDIF1 = ',
     *  PRDIF1
C
      READ(MDATA,*) PRDIF2
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Parameter for diffusive prec. : PRDIF2 = ',
     *  PRDIF2
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Parameter for diffusive prec. : PRDIF2 = ',
     *  PRDIF2
C
C
C
      IAUSAV=0
C
      IMTIME=1
      TTMG =0D0
      TTS  =0D0
      TTE  =0D0
      TTD  =0D0
      TTP  =0D0
      TTR  =0D0
C
      TTGRID=0D0
      TTPOST=0D0
      TTADF =0D0
      TTUPW =0D0
      TTBDR =0D0
      TTLC  =0D0
      TTILU =0D0
C
      TTMGU =0D0
      TTSU  =0D0
      TTEU  =0D0
      TTDU  =0D0
      TTPU  =0D0
      TTRU  =0D0
      TTMGP =0D0
      TTSP  =0D0
      TTEP  =0D0
      TTDP  =0D0
      TTPP  =0D0
      TTRP  =0D0
C
C
      IF (IFRSTP.EQ.1) THEN
       THETA =1D0-SQRT(0.5D0)
       THETAP=1D0-2D0*THETA
       FALPHA=THETAP/(1D0-THETA)
       FBETA =THETA /(1D0-THETA)
       THSTEP=3D0*TSTEP*FALPHA*THETA
      ELSE
       THSTEP=TSTEP*THETA
      ENDIF
C
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.0) WRITE(MFILE,1)
C
C=======================================================================
C     Point-value files
C=======================================================================
C
      OPEN (UNIT=40,FILE='#points/tf_0')
      OPEN (UNIT=41,FILE='#points/u1_0')
      OPEN (UNIT=42,FILE='#points/u2_0')
      OPEN (UNIT=43,FILE='#points/u3_0')
      OPEN (UNIT=44,FILE='#points/u4_0')
      OPEN (UNIT=45,FILE='#points/p1_0')
      OPEN (UNIT=46,FILE='#points/p2_0')
      OPEN (UNIT=47,FILE='#points/p3_0')
      OPEN (UNIT=48,FILE='#points/p4_0')
      OPEN (UNIT=49,FILE='#points/p5_0')
      OPEN (UNIT=50,FILE='#points/p6_0')
      OPEN (UNIT=51,FILE='#points/p7_0')
      OPEN (UNIT=52,FILE='#points/p8_0')
      OPEN (UNIT=53,FILE='#points/f1_0')
C
C
C
      END
