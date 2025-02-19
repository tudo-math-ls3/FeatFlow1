************************************************************************
      SUBROUTINE INIT1 (MDATA,CDATA,MFILE,MSHOW,IWORKG,IWMAXG)
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
      PARAMETER (NBLOCA=1,NBLOCB=2,NBLOCF=2)
      PARAMETER (NBLA1=NBLOCA*NNLEV,NBLB1=NBLOCB*NNLEV,NBLF1=
     *           NBLOCF*NNLEV)
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
      DATA ARRDA/'VA    '/,ARRDST/'DST   '/,ARRDF/'DF1   ','DF2   '/
      DATA ARRDM/'DM    '/
      DATA ARRDB/'DB1   ','DB2   '/,ARRDFP/'DFP   '/
C
C *** Structure of bilinear and linear forms 
      DIMENSION KABSTN(NBLOCA),KABST(2,NNAB,NBLOCA)
      DIMENSION KABBN(NBLOCB),KABB(2,NNAB,NBLOCB)
      DIMENSION KFN(NBLOCF),KF(NNAB,NBLOCF)
      DATA KABSTN/3/, KABBN/1,1/, KFN/1,1/
C
      DIMENSION BCONST(NBLOCA),BCONB(NBLOCB),BCONF(NBLOCF)
      DATA BCONST/.TRUE./, BCONB/.TRUE.,.TRUE./, BCONF/.FALSE.,.FALSE./
C
      DIMENSION BSNGLA(NBLOCA,NNLEV),BSNGLB(NBLOCB,NNLEV),BSNGLF(NBLOCF)
      DATA BSNGLA/NBLA1*.FALSE./,BSNGLB/NBLB1*.FALSE./
      DATA BSNGLF/.FALSE.,.FALSE./
C
      DIMENSION KLB(NBLOCB,NNLEV),KLF(NBLOCF,NNLEV),LF(NBLOCF)
      DIMENSION KNA(NNLEV),KNEQ(NNLEV),KNB(NNLEV),KNEQB(NNLEV)
      DIMENSION KLMH(NNLEV)
      DATA KLB/NBLB1*0/,KLF/NBLF1*0/,LF(1)/0/,LF(2)/0/
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
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
      COMMON /NSTIME/ TTGRID,TTPOST,TTADF,TTUPW,TTBDR,TTLC
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
      COMMON /MGFLD/  KLA(NNLEV),KLST(NNLEV),KLM(NNLEV),KLCOLA(NNLEV),
     *                KLLDA(NNLEV),KLB1(NNLEV),KLB2(NNLEV),
     *                KLCOLB(NNLEV),KLLDB(NNLEV),
     *                KLUP(NNLEV),KLF12P(NNLEV),KLAUX(NNLEV),
     *                KLAREA(NNLEV),LU1OLD,LU2OLD,LPOLD,LD1,LD2,LDP
      COMMON /MGDIM/  KNA,KNB,KNU(NNLEV),KNP(NNLEV),KNUP(NNLEV)
C
      COMMON /LEVDIM/ NA,NB,NU,NP,NUP
      COMMON /ADRFLD/ KA1,KST1,KM1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                KU1,KU2,KP,KF1,KF2,KFP,KAUX1,KAUX2,KAUXP
c
      INCLUDE 'jump.inc'
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
      CHARACTER CPARM1*60,CMESH1*60,CFILE1*60,CSTART*60,CSOL*60
      COMMON /FILES/ IMESH1,MMESH1,CPARM1,CMESH1,MFILE1,CFILE1,
     *               ISTART,MSTART,CSTART,ISOL,MSOL,CSOL
C

      INCLUDE 'block.inc'
      INCLUDE 'bouss.inc'
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
C
      ISYMMA=0
C
C=======================================================================
C     Data input
C=======================================================================
C
C      CALL  OF0  (MDATA,CDATA,1)
      CALL OPNDAT (MDATA,MSHOW,CDATA)
      CALL  GDAT (MDATA,MSHOW,IRMESH)
      CLOSE (MDATA)
      CFILE=CFILE1
      MFILE=MFILE1
C
c     Daten fuer das Blocken etc
C
      CDATA='#data/patchpara.dat'
      CALL  OF0  (MDATA,CDATA,1)
      CALL  GDAT2 (MDATA,MSHOW,IRMESH)
c
c      write (*,*) 'iinv', iinver
c
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
      IMID=2
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
c      DO 11  II=NLMIN,NLMAX
      DO 11  II=1,NLMAX
      ILEV=II
      NVT=KNVT(II)
      NMT=KNMT(II)
      NEL=KNEL(II)
      NVBD=KNVBD(II)
      KNU(II)=NMT
      KNP(II)=NEL
      NUP=NMT+NMT+NEL
      KNUP(II)=NUP
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'ILEV,NVT,NMT,NEL,NVBD: ', II,NVT,NMT,NEL,NVBD
      IF (MSHOW.GE.1) 
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
C-----------------------------------------------------------------
c   Fuer die Blockung nach GG-Elementen
C-----------------------------------------------------------------
c
      CALL ZNEW(NEL,3,KLELC1(II),'KELC1 ')
      CALL ZNEW(NEL,3,KLELC2(II),'KELC2 ')
C
c      IF (II.LE.1) THEN
      IF (II.LE.NSBLSL) THEN
       INDLEV=1
       NELMAC=NEL
      ELSE
       INDLEV=II-NSBLSL+1
      ENDIF
c$$$      IF (II.LT.NLMIN) goto 4
c$$$      IF (II.LE.nlmin) THEN
c$$$       INDLEV=1
c$$$       NELMAC=NEL
c$$$      ELSE
c$$$       INDLEV=II
c$$$      ENDIF
C
      CALL MAC2EL(KWORK(L(KLELC1(II)))  ,KWORK(L(KLELC2(II))),
     *            KWORK(L(KLELC1(II-1))),KWORK(L(KLELC2(II-1))),
     *            KWORK(L(KLADJ(II))),NEL,NELMAC,INDLEV)
C
c
 4    continue
c$$$      DO 4444 IEL=1,NEL
c$$$      WRITE(6,*) II,NELMAC,IEL,KWORK(L(KLELC1(II))+IEL-1),
c$$$     *                         KWORK(L(KLELC2(II))+IEL-1)
c$$$4444  CONTINUE
C-----------------------------------------------------------------
c
11    CONTINUE
C
C
      CALL ZTIME(TTT1)
      TTGRID=TTT1-TTT0
C
      IF (MSHOW.GE.2) WRITE(MTERM,*)
      IF (MSHOW.GE.1) WRITE(MFILE,*)
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'time for grid generation : ',
     *                                TTGRID
      IF (MSHOW.GE.1) WRITE(MFILE,*) 'time for grid generation : ', 
     *                                TTGRID
      IF (MSHOW.GE.2) WRITE(MTERM,*)
      IF (MSHOW.GE.1) WRITE(MFILE,*)
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
      IF(IJUMP .EQ. 0) CALL XMAP7(KLCOLA,KLLDA,KNA,KNEQ,E031,ISYMMA)
      IF(IJUMP .GE. 1)CALL XMAPP7(KLCOLA,KLLDA,KNA,KNEQ,KMEL,E031,ISYMM)
      IF (IER.NE.0) GOTO 99999
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
C *** Generation of mass matrix
C=======================================================================
C
      IF (ISTAT.NE.0) THEN
       KABST (1,1,1)=1
       KABST (2,1,1)=1
       KABSTN(1)    =1
C
       ICLRA=1
       IF (IELT.EQ.0)
     *  CALL XMAB07(KLMH,KLCOLA,KLLDA,KNA,KNEQ,NBLOCA,ICLRA,E031,
     *              COEFST,BCONST,KABST,KABSTN,ABS(ICUBM),ISYMMA,CARRM,
     *              BSNGLA)
       IF (IELT.EQ.1)
     *  CALL XMAB07(KLMH,KLCOLA,KLLDA,KNA,KNEQ,NBLOCA,ICLRA,E030,
     *              COEFST,BCONST,KABST,KABSTN,ABS(ICUBM),ISYMMA,CARRM,
     *              BSNGLA)
       IF (IELT.EQ.2)
     *  CALL XMABM7(KLMH,KLCOLA,KLLDA,KNA,KNEQ,NBLOCA,ICLRA,EM31,
     *              COEFST,BCONST,KABST,KABSTN,ABS(ICUBM),ISYMMA,CARRM,
     *              BSNGLA)
       IF (IELT.EQ.3) 
     *  CALL XMABM7(KLMH,KLCOLA,KLLDA,KNA,KNEQ,NBLOCA,ICLRA,EM30,
     *              COEFST,BCONST,KABST,KABSTN,ABS(ICUBM),ISYMMA,CARRM,
     *              BSNGLA)
       IF (IER.NE.0) GOTO 99999
C
C=======================================================================
C
       DO 20  II=NLMIN,NLMAX
C
       IF (IMASS.EQ.0) THEN
        CALL  ZNEW (KNU(II),1,KLM(II),'VMASS ')
        IF (IER.NE.0) GOTO 99999
C
        DO 21 INU=1,KNU(II)
C
        IF (IMASSL.EQ.0) THEN
         DMH=0D0
         DO 22 ILD=KWORK(L(KLLDA(II))+INU-1),KWORK(L(KLLDA(II))+INU)-1
22       DMH=DMH+DWORK(L(KLMH(II))+ILD-1)
        ELSE
         ILD=KWORK(L(KLLDA(II))+INU-1)
         DMH=DWORK(L(KLMH(II))+ILD-1)
        ENDIF
C
        DWORK(L(KLM(II))+INU-1)=DMH
c        WRITE(6,*) II,INU,DWORK(L(KLM(II))+INU-1),KLM(II),L(KLM(II))
21      CONTINUE
C
        CALL  ZDISP (0, KLMH(II),'VMASSH')
        IF (IER.NE.0) GOTO 99999
       ELSE
        KLM(II) =KLMH(II)
        KLMH(II)=0
       ENDIF
C
20     CONTINUE
C
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
      IF (IPRECB.GE.2) THEN
       KNEQB(ILEV)=KNU(ILEV)
       KNB  (ILEV)=4*KNEL(ILEV)
       CALL ZNEW(2*KNB(ILEV)  ,1,KLB(1,ILEV) ,ARRDB(1))
       CALL ZNEW(2*KNB(ILEV)  ,1,KLB(2,ILEV) ,ARRDB(2))
       CALL ZNEW(KNB(ILEV)  ,3,KLCOLB(ILEV),'KCOLB ')
       CALL ZNEW(KNU(ILEV)+1,3,KLLDB (ILEV),'LLDB  ')
       IF (IER.NE.0) GOTO 99999
C
c
       ISETLV=2
       CALL SETLEV(ISETLV)
C
       CALL BBUILD(KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *             DWORK(L(LCORVG)),
     *             DWORK(L(KLB(1,ILEV))),DWORK(L(KLB(2,ILEV))),
     *             KWORK(KCOLB),KWORK(KLDB),KNB(ILEV),NEL,NVT,NMT)
c
       CALL BBUILD(KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *             DWORK(L(LCORVG)),
     *             DWORK(L(KLB(1,ILEV))+KNB  (ILEV)),
     *             DWORK(L(KLB(2,ILEV))+KNB  (ILEV)),
     *             KWORK(KCOLB),KWORK(KLDB),KNB(ILEV),NEL,NVT,NMT)
C

c       CALL ZDISP (2*KNB(ILEV),KLB(1,ILEV) ,ARRDB(1))
c       CALL ZDISP (2*KNB(ILEV),KLB(2,ILEV) ,ARRDB(2))
c       CALL ZDISP (KNB(ILEV),KLCOLB(ILEV),'KCOLB ')
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
C
13    CONTINUE
C

C
      IF (MSHOW.GE.2) WRITE(MTERM,*)
      IF (MSHOW.GE.1) WRITE(MFILE,*)
C 
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
C=======================================================================
C     matrix restructuring
C=======================================================================
C
      DO 15 ILEV=NLMIN,NLMAX
       CALL ZTIME(TTT0)
       ISETLV=2
       CALL SETLEV(ISETLV)
C
       IF (IPRECA.EQ.0) THEN
        CALL ZCTYPE(2,KLST(ILEV),CARRST)
        IF (IER.NE.0) GOTO 99999
C
        IF (ISTAT.EQ.1) THEN
c        write (*,*) 'typenumwandlung: HIER'
         CALL ZCTYPE(2,KLM(ILEV),CARRM)
         IF (IER.NE.0) GOTO 99999
        ENDIF
       ENDIF
C
C


       IF (IPRECA.EQ.2) THEN
        CALL ZCTYPE(2,KLST(ILEV),CARRST)
        IF (IER.NE.0) GOTO 99999
        CALL  OF0 (59,CFILST(ILEV),0)
        CFILE='STMAT '
        CALL  OWA2 (VWORK(L(KLST(ILEV))),CFILE,NA,59,0)
        REWIND(59)
        CLOSE (59)
        IF (IER.NE.0) GOTO 99999
        CALL ZDISP (0,KLST(ILEV),'STMAT ')
        IF (IER.NE.0) GOTO 99999
C
        IF (ISTAT.EQ.1) THEN
         CALL ZCTYPE(2,KLM(ILEV),CARRM)
         IF (IER.NE.0) GOTO 99999
         CALL  OF0 (59,CFILM(ILEV),0)
         CFILE='MASMAT'
         CALL  OWA2V (VWORK(L(KLM(ILEV))),CFILE,NA,59,0)
         REWIND(59)
         CLOSE (59)
         IF (IER.NE.0) GOTO 99999
         CALL ZDISP (0,KLM(ILEV),'MASMAT')
         IF (IER.NE.0) GOTO 99999
        ENDIF
       ENDIF
C
C
       IF (IPRECA.EQ.3) THEN
        CALL  OF0 (59,CFILST(ILEV),0)
        CFILE='STMAT '
        CALL  OWA1 (DWORK(L(KLST(ILEV))),CFILE,NA,59,0)
        REWIND(59)
        CLOSE (59)
        IF (IER.NE.0) GOTO 99999
        CALL ZDISP (0,KLST(ILEV),'STMAT ')
        IF (IER.NE.0) GOTO 99999
C
        IF (ISTAT.EQ.1) THEN
         CALL  OF0 (59,CFILM(ILEV),0)
         CFILE='MASMAT'
         CALL  OWA1 (DWORK(L(KLM(ILEV))),CFILE,NA,59,0)
         REWIND(59)
         CLOSE (59)
         IF (IER.NE.0) GOTO 99999
         CALL ZDISP (0,KLM(ILEV),'MASMAT')
         IF (IER.NE.0) GOTO 99999
        ENDIF
       ENDIF
C
c333       IF (IPRECA.EQ.4) THEN
c        write (*,*) 'typenumwandlung: HIER'
c333         CALL ZCTYPE(2,KLM(ILEV),CARRM)
c333     IF (IER.NE.0) GOTO 99999
c333       ENDIF
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
C        KLB1(ILEV)=0
C        KLB2(ILEV)=0
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
C                    - UOLD-vector on level NLMAX
C=======================================================================
C
      DO 25  ILEV=NLMIN,NLMAX
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
C
      CALL ZNEW(NUP,1,LUP,'DU12P ')
      IF (IER.NE.0) GOTO 99999
      KLUP(ILEV)=LUP
C
      IF (ILEV.EQ.NLMAX)  THEN
c        IF ((OMGMIN.GT.0D0).OR.(OMGMAX.GT.0D0)) THEN
         CALL ZNEW(NEQ,1,LU1OLD,'DU1OLD')
         IF (IER.NE.0) GOTO 99999
         CALL ZNEW(NEQ,1,LU2OLD,'DU2OLD')
         IF (IER.NE.0) GOTO 99999
         CALL ZNEW(NEL,1,LPOLD,'DPOLD ')
         IF (IER.NE.0) GOTO 99999
c        ENDIF
C
        CALL ZNEW(NEQ,1,LD1,'DD1   ')
        IF (IER.NE.0) GOTO 99999
        CALL ZNEW(NEQ,1,LD2,'DD2   ')
        IF (IER.NE.0) GOTO 99999
        CALL ZNEW(NEL,1,LDP,'DDP   ')
        IF (IER.NE.0) GOTO 99999
      ENDIF
C
      CALL ZNEW(NUP,1,LF12P,'DF12P ')
      IF (IER.NE.0) GOTO 99999
      KLF12P(ILEV)=LF12P
C
C allocation of a vector for local viscosities !CHANGEDEV !
C
      CALL ZNEW(NEL,1,LNY,'VLNY  ')
      IF (IER.NE.0) GOTO 99999
      KLNY(ILEV)=LNY
c
C=======================================================================
C *** Allocation of iteration matrix A 
C=======================================================================
C
      NA=KNA(ILEV)
      CALL  ZNEW (5*NA,2,LA1,'VA    ') !CHANGEDEV !
      IF (IER.NE.0) GOTO 99999
      KLA(ILEV)=LA1
C=======================================================================
C *** Allocation of  matrix AV for the preconditioner 
C=======================================================================
C
      NAV=KNA(ILEV)
      CALL  ZNEW (NA,2,LAV1,'VAV   ') !CHANGEDEV !
      IF (IER.NE.0) GOTO 99999
      KLAV(ILEV)=LAV1
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
C=======================================================================
C *** start vector as prolongation from level NLMAX
C=======================================================================
C
c$$$      IF (ISTART.EQ.0) CALL UINIT(KWORK(L(LVERT)),KWORK(L(LMID)),
c$$$     *    DWORK(L(LCORVG)),DWORK(L(LUP)),DWORK(L(LUP)+NEQ),NEL,NVT)  
C
      IF (ABS(ISTART).EQ.1) THEN
       IF (ISTART.EQ.1) THEN
        IFMTS=0
       ELSE
        IFMTS=1
       ENDIF
       CFILE='DU12P '
       CALL  OF0 (MSTART,CSTART,IFMTS)
       CALL  ORA1 (DWORK(L(LUP)),CFILE,MSTART,IFMTS)
       CLOSE(MSTART)
      ENDIF
C
C=======================================================================
C *** start vector as prolongation from level NLMAX
C=======================================================================

      IF (ABS(ISTART).EQ.1) THEN
       IF (ISTART.EQ.1) THEN
        IFMTS=0
       ELSE
        IFMTS=1
       ENDIF
       CFILE='DU12P '
       CALL  OF0 (MSTART,CSTART,IFMTS)
       CALL  ORA1 (DWORK(L(LUP)),CFILE,MSTART,IFMTS)
       CLOSE(MSTART)
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
        CALL  ORA1 (DWORK(KU1C),CFILE,MSTART,IFMTS)
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
C	Generation of Vector KMEL, which contains the numbers of
C       the elements adjacent to each midpoint. 
C=======================================================================
C
      DO 17 ILV=NLMIN,NLMAX
        CALL ZNEW (2*KNMT(ILV),3,KLMEL(ILV),'KMEL  ')
C
        CALL MKKMEL(KWORK(L(KLMEL(ILV))),KWORK(L(KLMID(ILV))),
     *              KNEL(ILV),KNMT(ILV),KNVT(ILV))
17    CONTINUE
C=======================================================================
C *** start vector is the inflow profile in the whole OMEGA
C=======================================================================
C
      IF (ABS(ISTART).EQ.3) THEN
      NU  =KNU(NLMAX)
      NVT =KNVT(NLMAX)
      KU1 =L(KLUP(NLMAX))
      KU2 =KU1+NU
c
c$$$      CALL INLOES (DWORK(KU1),DWORK(KU2),
c$$$     *            NU,NVT,KNVBD(NLMAX),KWORK(L(KLMEL(NLEV))),
c$$$     *            KWORK(L(KLVERT(NLEV))), DWORK(L(KLMBDP (NLEV))),
c$$$     *            KWORK(L(KLMBD (NLEV))),
c$$$     *            kwork(l(lnpr)),DWORK(L(LCORVG)),1)
c$$$

c$$$      write (*,*) (i,dwork(ku1-1+i),dwork(ku2-1+i),'\n',i=1,50)
c$$$      print *
       ENDIF
C=======================================================================
C *** Dirichlet boundary updates of solution and rhs vector
C=======================================================================
C
      CALL ZTIME(TTT0)
      IF ((IBDR.EQ.1).AND.(ISTAT.EQ.0)) THEN
       CALL PDSET (KWORK(L(KLVBD(NLEV))),KWORK(L(KLEBD(NLEV))),
     *             KWORK(L(KLVERT(NLEV))),KWORK(L(KLMID(NLEV))),
     *             KWORK(L(KLNPR(NLEV))),DWORK(L(LCORVG)),
     *             DWORK(L(KLMBDP(NLEV))),DWORK(L(LF1)),DWORK(L(LF2)),
     *             1D0)
      ENDIF
C
      KU1=L(LUP)
      KU2=KU1+NEQ
      CALL BDRSET (DWORK(KU1),DWORK(KU2),DWORK(L(LF1)),
     *             DWORK(L(LF2)),KWORK(L(KLMBD(NLEV))),
     *             DWORK(L(KLDBD(NLEV))),KWORK(L(LNPR)),
     *             KNMBD(NLEV),NVT,PARX,PARY,UE)
C
      CALL  LCP1 (DWORK(L(LF1)), DWORK(L(LF12P)),     NEQ)
      CALL  LCP1 (DWORK(L(LF2)), DWORK(L(LF12P)+NEQ), NEQ)
      CALL ZTIME(TTT1)
      TTBDR=TTBDR+TTT1-TTT0
C
      CALL ZTIME(TTT0)
      CALL  ZDISP (0, LF1, ARRDF(1) )
      CALL  ZDISP (0, LF2, ARRDF(2) )
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
      CALL ZNEW(NUP,1,LAUX,'DAUX  ')
      IF (IER.NE.0) GOTO 99999
      KLAUX(ILEV)=LAUX
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0
C
25    CONTINUE
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
      SUBROUTINE GDAT (MDATA,MSHOW,IRMESH)
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
C *** COMMON blocks for multigrid routine M010
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTIME/ TTMG,TTS,TTE,TTD,TTP,TTR,IMTIME
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
C *** COMMON blocks for time discretization
      COMMON /NSPAR/  TSTEP,THETA,THSTEP,TIMENS,EPSNS,NITNS,ITNS,ISTAT
      COMMON /NSFRAC/ THETAP,FALPHA,FBETA,IFRSTP
      COMMON /NSADAT/ TIMEMX,DTMIN,DTMAX,DTFACT,TIMEIN,EPSADI,EPSADL,
     *                EPSADU,IEPSAD,IADIN,IREPIT,IADTIM
      COMMON /NSSAV/  INSAV,INSAVN
      COMMON /NSSAVF/ DTFILM,DTFILO,DTAVS,DTAVSO,DTGMV,DTGMVO,
     *                IFUSAV,IFPSAV,IFXSAV,IAVS,IGMV,IFINIT
      COMMON /NSTIME/ TTGRID,TTPOST,TTADF,TTUPW,TTBDR,TTLC
      COMMON /NSEXL/  ITEXL,LTML,TIML11,TIML12,TIML31,TIML32
C
      CHARACTER CPARM1*60,CMESH1*60,CFILE1*60,CSTART*60,CSOL*60
      COMMON /FILES/ IMESH1,MMESH1,CPARM1,CMESH1,MFILE1,CFILE1,
     *               ISTART,MSTART,CSTART,ISOL,MSOL,CSOL
C
      INCLUDE 'jump.inc'
C
      INCLUDE 'bouss.inc'
      INCLUDE 'block.inc'
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
      READ(MDATA,*) IELTNEWTON
      IF ((IELTNEWTON.LT.0).OR.(IELTNEWTON.GT.3)) IELTNEWTON=3
      IF (MSHOW.GE.2)
     *  WRITE(MTERM,*) 'element type for newton term    =', IELTNEWTON
      IF (MSHOW.GE.0) 
     *  WRITE(MFILE,*) 'element type for newton term   = ', IELTNEWTON
C
      READ(MDATA,*) IELTJUMP
      IF ((IELTJUMP.LT.0).OR.(IELTJUMP.GT.3)) IELTJUMP=3
      IF (MSHOW.GE.2)
     *  WRITE(MTERM,*) 'element type for newton term    =', IELTJUMP
      IF (MSHOW.GE.0) 
     *  WRITE(MFILE,*) 'element type for newton term   = ', IELTJUMP
C
      READ(MDATA,*)  ISTOK
      IF (ISTOK.NE.1) ISTOK=0
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'Stokes calculation:  ISTOK = ', ISTOK
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'Stokes calculation:  ISTOK = ', ISTOK
C
      READ(MDATA,*) IRHS
      IF ((IRHS.LT.0).OR.(IRHS.GT.2)) IRHS=2
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
      READ(MDATA,*) ICHDEF
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'Error evaluation    = ',ICHDEF
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'Error evaluation    = ',ICHDEF
C
      READ(MDATA,*) IMASS
      IF (IMASS.NE.1) IMASS=0
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'mass evaluation     = ',IMASS
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'mass evaluation     = ',IMASS
C
      READ(MDATA,*) IMASSL
      IF (IMASSL.NE.1) IMASSL=0
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'lumped mass eval.   = ',IMASSL
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'lumped mass eval.   = ',IMASSL
C
      READ(MDATA,*) IUPW
      IF ((IUPW.GT.1).OR.(IUPW.LT.0)) IUPW=1
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'convective part     = ',IUPW
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'convective part     = ',IUPW
C
      READ(MDATA,*) IPRECA !CHANGEDEV !
      IF ((IPRECA.LT.0).OR.(IPRECA.GT.4)) IPRECA=4
c      IF ((IPRECA.EQ.4).AND.(IUPW.EQ.1))  IPRECA=1
c      IF ((IPRECA.EQ.4).AND.(ISTOK.EQ.1)) IPRECA=1
      IF (IPRECA.EQ.3)                    IPRECA=4
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'Accuracy for ST     = ',IPRECA
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'Accuracy for ST     = ',IPRECA
C
      READ(MDATA,*) IPRECB
      IF ((IPRECB.LT.0).OR.(IPRECB.GT.4)) IPRECB=2
      IF (IPRECB.EQ.1) IPRECB=0
      IF (IPRECB.EQ.2) IPRECB=3
      IF (IPRECB.EQ.4) IPRECB=3
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'Accuracy for B      = ',IPRECB
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'Accuracy for B      = ',IPRECB
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
      READ(MDATA,*) ICUBNEWTON
      ICUBF=ABS(ICUBF)
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'ICUB FOR NEWTON TERM=',ICUBNEWTON
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'ICUB FOR NEWTON TERM=',ICUBNEWTON
C
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
      READ(MDATA,*)  ICYC
      ICYCLE=ABS(ICYC)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'type of mg-cycle:  ICYCLE = ', ICYCLE
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'type of mg-cycle:  ICYCLE = ', ICYCLE
C
      READ(MDATA,*) ILMIN
      ILMIN=ABS(ILMIN)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'minimum of linear mg steps :  ILMIN = ', ILMIN
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'minimum of linear mg steps :  ILMIN = ', ILMIN
C
      READ(MDATA,*) ILMAX
      ILMAX=ABS(ILMAX)
      IF (ILMAX.LT.ILMIN) ILMAX=ILMIN
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'maximum of linear mg steps :  ILMAX = ', ILMAX
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'maximum of linear mg steps :  ILMAX = ', ILMAX
C
      READ(MDATA,*)  IINT
      IF ((ABS(IINT).NE.1).AND.(ABS(IINT).NE.2)) IINT=2
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'type of interpolation:  IINT = ', IINT
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'type of interpolation:  IINT = ', IINT
C
      READ(MDATA,*) ISM
      IF ((ISM.LT.1).AND.(ISM.GT.4)) ISM=1
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'type of smoother :  ISM = ',ISM
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'type of smoother :  ISM = ',ISM
C
      READ(MDATA,*) ISL
      IF ((ISL.LT.1).AND.(ISL.GT.3)) ISL=1
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'type of solver :  ISL = ',ISL
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'type of solver :  ISL = ',ISL
C
      READ(MDATA,*) NSM
      NSM=ABS(NSM)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'number of smoothing steps :  NSM = ', NSM
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'number of smoothing steps :  NSM = ', NSM
C
      READ(MDATA,*) NSL
      NSL=ABS(NSL)
      IF (ISL.EQ.3) then
         NSLOLD=NSL
         NSL=1
         ENDIF
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'number of solver steps :  NSL = ', NSL
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'number of solver steps :  NSL = ', NSL
C
      READ(MDATA,*) NSMFAC
      IF (NSMFAC.LT.1) NSMFAC=1
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'factor sm. steps on coarser lev.:NSMFAC=',NSMFAC
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'factor sm. steps on coarser lev.:NSMFAC=',NSMFAC
C
      DO 11  II=1,NNLEV
      KPOSM(II)=NSM
11    KPRSM(II)=NSM
C
      DO 12  II=1,NLEV
      KPOSM(II)=KPOSM(II)*NSMFAC**(NLEV-II)
      KPRSM(II)=KPRSM(II)*NSMFAC**(NLEV-II)
12    CONTINUE
C
      DO 13  II=1,NNLEV
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'KPRSM,KPOSM ON LEVEL: ',II,KPRSM(II),KPOSM(II)
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'KPRSM,KPOSM ON LEVEL: ',II,KPRSM(II),KPOSM(II)
13    CONTINUE
C
C
      READ(MDATA,*) IGRAD
      IF (MSHOW.GE.2) THEN
         IF (IGRAD.eq.0) THEN
          WRITE(MTERM,*) 'Use just Def-Tensor for OUTER it.: IGRAD = ', 
     *                   IGRAD
         ELSE
          WRITE(MTERM,*) 'Use the gradient for OUTER it.: IGRAD = ', 
     *                   IGRAD
         ENDIF
      ENDIF
      IF (MSHOW.GE.0) THEN
         IF (IGRAD.eq.0) THEN
          WRITE(MFILE,*) 'Use just Def-Tensor for OUTER it.:  IGRAD = ', 
     *                   IGRAD
         ELSE
          WRITE(MFILE,*) 'Use the gradient for OUTER it.:: IGRAD = ', 
     *                   IGRAD
         ENDIF
       ENDIF
C
      READ(MDATA,*) IFIXMIN
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*)'minimum number of fixed point iteration 
     *                                          : IFIXMIN = ',IFIXMIN
      IF (MSHOW.GE.0) 
     * WRITE(MTERM,*)'minimum number of fixed point iteration 
     *                                          : IFIXMIN = ',IFIXMIN
      READ(MDATA,*) IFIXMAX
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*)'maximum number of fixed point iteration
     *                                          : IFIXMAX = ',IFIXMAX
      IF (MSHOW.GE.0) 
     * WRITE(MTERM,*)'maximum number of iteration of fixed point
     *                                          : IFIXMAX = ',IFIXMAX
C
      READ(MDATA,*) INEWTON
      IF (MSHOW.GE.2) THEN
         IF (INEWTON.eq.0) THEN
          WRITE(MTERM,*) 'NEWTON METHOD ARE NOT USED.: INEWTON = ', 
     *                   INEWTON
         ELSE
          WRITE(MTERM,*) 'NEWTON METHOD ARE NOT USED.: INEWTON = ', 
     *                   INEWTON
         ENDIF
      ENDIF
C
      READ(MDATA,*) ISPGRAD
      IF (MSHOW.GE.2) THEN
         IF (ISPGRAD.eq.0) THEN
          WRITE(MTERM,*) 'THE GRADIENT TENSOR ARE USED.: ISPGRAD = ', 
     *                                           ISPGRAD
         ELSE
          WRITE(MTERM,*) 'THE SYMETRIC PART OF GRADIENT ARE USED.:, 
     *                                           ISPGRAD = ', ISPGRAD 
         ENDIF
      ENDIF
C
      READ(MDATA,*) IJUMP
      IF (MSHOW.GE.2) THEN
         IF (ISPGRAD.eq.0) THEN
          WRITE(MTERM,*) 'THE JUMP                  .:IJUMP  = ',IJUMP
         ELSE
          WRITE(MTERM,*) 'THE JUMP                 .:,IJUMP = ', IJUMP 
         ENDIF
      ENDIF
C
      READ(MDATA,*) IPRECO
      IF (IGRAD.EQ.1) IPRECO=1
c
      IF (MSHOW.GE.2) THEN
         IF (IPRECO.eq.0) THEN
          WRITE(MTERM,*) 'Use just Def-Tensor for PRECOND.:  IPRECO = ', 
     *                   IPRECO
         ELSEIF (IPRECO.EQ.1) THEN
          WRITE(MTERM,*) 'Use the gradient for PRECOND.: IPRECO = ', 
     *                   IPRECO
         ELSE
          WRITE(MTERM,*)  
     *     'Mixed version preconditioner:  '
          ENDIF
      ENDIF
      IF (MSHOW.GE.0) THEN
         IF (IPRECO.eq.0) THEN
          WRITE(MFILE,*) 'Use just Def-Tensor for PRECOND.:  IPRECO = ', 
     *                   IPRECO
         ELSEIF (IPRECO.EQ.1) THEN
          WRITE(MFILE,*) 'Use the gradient for PRECOND.: IPRECO = ', 
     *                   IPRECO
         ELSE
          WRITE(MFILE,*) 
     *     'Mixed version preconditioner: '
         ENDIF
       ENDIF
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
      READ(MDATA,*) PLA
      IF (MSHOW.GE.2) 
     *     WRITE(MTERM,*) 'Power-law parameter:  ALPHA = ', PLA
      IF (MSHOW.GE.0) 
     *     WRITE(MFILE,*) 'Power-law parameter:  ALPHA = ', PLA
C
      READ(MDATA,*) PLE
      IF (MSHOW.GE.2) 
     *     WRITE(MTERM,*) 'Power-law parameter:  EPSILON = ', PLE
      IF (MSHOW.GE.0) 
     *     WRITE(MFILE,*) 'Power-law parameter:  EPSILON = ', PLE
C
      READ(MDATA,*) DCU1
      IF (MSHOW.GE.2) 
     *     WRITE(MTERM,*) 'EXTERIOR CYLIDER :  DCU1 = ', DCU1
      IF (MSHOW.GE.0) 
     *     WRITE(MFILE,*) 'EXTERIOR CYLIDER :  DCU1 = ', DCU1
C
      READ(MDATA,*) DCU2 
      IF (MSHOW.GE.2) 
     *     WRITE(MTERM,*) 'INTERIOR CYLIDER :  DCU1 = ', DCU2 
      IF (MSHOW.GE.0) 
     *     WRITE(MFILE,*) 'INTERIOR CYLIDER :  DCU1 = ', DCU2
C
      READ(MDATA,*) DMEAN
      IF (MSHOW.GE.2) 
     *     WRITE(MTERM,*) ' THE MEAN VALUE OF THE PRESSURE = ',DMEAN
      IF (MSHOW.GE.0) 
     *     WRITE(MFILE,*) ' THE MEAN VALUE OF THE PRESSURE = ',DMEAN 
C
      READ(MDATA,*) IMODEL
      IF (MSHOW.GE.2) 
     *     WRITE(MTERM,*) 'DIFFERENT TYPE OF MODEL:  IMODEL = ',IMODEL
      IF (MSHOW.GE.0) 
     *     WRITE(MFILE,*) 'DIFFERENT TYPE OF MODEL:  IMODEL = ',IMODEL 
C
      READ(MDATA,*) INEWPRE
      IF (MSHOW.GE.2) 
     *     WRITE(MTERM,*) 'NEWTON FOR DEPPENDENT PRESSURE  = ',INEWPRE
      IF (MSHOW.GE.0) 
     *     WRITE(MFILE,*) 'NEWTON FOR DEPPENDENT PRESSURE  = ',INEWPRE 
C
      READ(MDATA,*) DNP
      IF (MSHOW.GE.2) 
     *     WRITE(MTERM,*) 'RELAXATION PARA. FOR NEWT. PRESSURE  = ',DNP
      IF (MSHOW.GE.0) 
     *     WRITE(MFILE,*) 'RELAXATION PARA. FOR NEWT. PRESSURE  = ',DNP 
C
C
      READ(MDATA,*) ALPHA
      IF (MSHOW.GE.2) 
     *     WRITE(MTERM,*) 'COEFF FOR PRESSURE :  ALPHA = ',ALPHA
      IF (MSHOW.GE.0) 
     *     WRITE(MFILE,*) 'COEFF FOR PRESSURE :  ALPHA = ',ALPHA 
C
      READ(MDATA,*)DBETA
      IF (MSHOW.GE.2) 
     *     WRITE(MTERM,*) 'COEFF FOR PRESSURE :  BETA = ',DBETA
      IF (MSHOW.GE.0) 
     *     WRITE(MFILE,*) 'COEFF FOR PRESSURE :  BETA = ',DBETA 
C
      READ(MDATA,*)DEPSI1
      IF (MSHOW.GE.2) 
     *     WRITE(MTERM,*) 'FOR PRESSURE :  DEPSI1 = ',DEPSI1
      IF (MSHOW.GE.0) 
     *     WRITE(MFILE,*) 'FOR PRESSURE :  DEPSI1 = ',DEPSI1 
C
      READ(MDATA,*)DEPSI2
      IF (MSHOW.GE.2) 
     *     WRITE(MTERM,*) 'FOR PRESSURE :  DEPSI2 = ',DEPSI2
      IF (MSHOW.GE.0) 
     *     WRITE(MFILE,*) 'FOR PRESSURE :  DEPSI2 = ',DEPSI2 
C
C
      READ(MDATA,*) DNPWSTART
      IF (MSHOW.GE.2) 
     *  WRITE(MTERM,*) 'Newton method for power_law term = ', DNPWSTART
      IF (MSHOW.GE.0) 
     *   WRITE(MFILE,*) 'Newton method for Power-law  = ', DNPWSTART
C
      READ(MDATA,*) DNNS
      IF (MSHOW.GE.2) 
     *  WRITE(MTERM,*) 'Newton method for convective term: DNNS=', DNNS
      IF (MSHOW.GE.0) 
     *  WRITE(MFILE,*) 'Newton method for convective term: DNNS=', DNNS
C
      READ(MDATA,*) DTOLERNEW
      IF (MSHOW.GE.2) 
     *  WRITE(MTERM,*) 'THE TOLERNECE FOR NEWTON :DTOLERNEW=', DTOLERNEW
      IF (MSHOW.GE.0) 
     *  WRITE(MFILE,*) 'THE TOLERNECE FOR NEWTON :DTOLERNEW=', DTOLERNEW 
C
      READ(MDATA,*) DJUMP
      IF (MSHOW.GE.2) 
     *  WRITE(MTERM,*) 'Relaxation parameter for the jump: DJUMP=',DJUMP
      IF (MSHOW.GE.0) 
     *  WRITE(MFILE,*) 'Relaxation parameter for the jump: DJUMP=',DJUMP

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
      READ(MDATA,*) EPSD
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'limit for U-defects   :        EPSD   = ', EPSD 
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'limit for U-defects   :        EPSD   = ', EPSD 
      READ(MDATA,*) EPSDIV
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'limit for DIV-defects :        EPSDIV = ', EPSDIV
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'limit for DIV-defects :        EPSDIV = ', EPSDIV
      READ(MDATA,*) EPSUR
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'limit for U-changes :          EPSUR  = ', EPSUR
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'limit for U-changes :          EPSUR  = ', EPSUR
      READ(MDATA,*) EPSPR
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'limit for P-changes :          EPSPR  = ', EPSPR
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'limit for P-changes :          EPSPR  = ', EPSPR
      READ(MDATA,*) DMPD
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'defect improvement  :          DMPD   = ', DMPD
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'defect improvement  :          DMPD   = ', DMPD
      READ(MDATA,*) DMPMG
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'damping of MG residuals     :  DMPMG  = ', DMPMG
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'damping of MG residuals     :  DMPMG  = ', DMPMG
      READ(MDATA,*) EPSMG
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'limit for MG residuals      :  EPSMG  = ', EPSMG
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'limit for MG residuals      :  EPSMG  = ', EPSMG
      READ(MDATA,*) DMPSL
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*)'damping of residuals for solving:  DMPSL = ',DMPSL
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*)'damping of residuals for solving:  DMPSL = ',DMPSL
      READ(MDATA,*) EPSSL
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*)'limit of changes for solving:    EPSSL = ',EPSSL
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*)'limit of changes for solving:    EPSSL = ',EPSSL
C
C
      READ(MDATA,*) RLXSM
      RLXSM=ABS(RLXSM)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'relaxation for the U-smoother: RLXSM = ', RLXSM
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'relaxation for the U-smoother: RLXSM = ', RLXSM
C
      READ(MDATA,*) RLXSL
      RLXSL=ABS(RLXSL)
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'relaxation for the U-solver :  RLXSL = ', RLXSL
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'relaxation for the U-solver :  RLXSL = ', RLXSL
C
      READ(MDATA,*) AMINMG
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*)'lower limit optimal MG-ALPHA: AMINMG = ', AMINMG
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*)'lower limit optimal MG-ALPHA: AMINMG = ', AMINMG
      READ(MDATA,*) AMAXMG
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*)'upper limit optimal MG-ALPHA: AMAXMG = ', AMAXMG
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*)'upper limit optimal MG-ALPHA: AMAXMG = ', AMAXMG
C
C
C-----------------------------------------------------------------------
C *** Values for /NS.../
C-----------------------------------------------------------------------
C
      MSOLD=MSHOW
      MSHOW=0
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'Parameters of /NS.../ :'
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'Parameters of /NS.../ :'
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.0) WRITE(MFILE,1)
      READ(MDATA,*) 
      READ(MDATA,*)
      READ(MDATA,*)
C
      READ(MDATA,*) ISTAT
      IF ((ISTAT.LT.0).OR.(ISTAT.GT.1)) ISTAT=0
      WRITE(MTERM,*) 'Time dependency          : ISTAT  = ', ISTAT
      WRITE(MFILE,*) 'Time dependency          : ISTAT  = ', ISTAT
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
      MSHOW=MSOLD
      READ(MDATA,*) 
c
      READ(MDATA,*) IZBV1
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'ZBV- Variable 1    : IZBV1 = ',
     *  IZBV1
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*)  'ZBV- Variable 1    : IZBV1 = ',
     *  IZBV1
C
c
      READ(MDATA,*) IZBV2
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'ZBV- Variable 2    : IZBV2 = ',
     *  IZBV2
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*)  'ZBV- Variable 2    : IZBV2 = ',
     *  IZBV2
      IZBV2=0! Dummy!!!!
c
      READ(MDATA,*) IZBV3
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'ZBV- Variable 3    : IZBV3 = ',
     *  IZBV3
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*)  'ZBV- Variable 3    : IZBV3 = ',
     *  IZBV3
C
      READ(MDATA,*) IZBV4
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'ZBV- Variable 4    : IZBV3 = ',
     *  IZBV4
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*)  'ZBV- Variable 4    : IZBV3 = ',
     *  IZBV4
C
      IF (ISTAT.EQ.0) THEN
       IAUSAV=INSAV
      ELSE
       IAUSAV=0
      ENDIF
C
      IMTIME=2
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
C
C
      IF (ISTAT.EQ.0) THEN
       IADTIM=0
       TSTEP =1D0
       THETA =1D0
       IFRSTP=0
      ENDIF
C
      IF ((ISTAT.NE.0).AND.(IFRSTP.EQ.1)) THEN
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

C=======================================================================
      SUBROUTINE GDAT2 (MDATA,MSHOW,IRMESH)
C=======================================================================
C  Liest die Blockingdaten ein
C=======================================================================
c      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
c     PARAMETER (NNLEV=9)
      INCLUDE 'common.inc'
      INCLUDE 'block.inc'
      SAVE
c
      READ (MDATA,*) 
      READ (MDATA,*)
      READ (MDATA,*)
      READ (MDATA,*)
      READ (MDATA,*)
      READ (MDATA,*)
c
      READ (MDATA,*) BANISO
      IF (MSHOW.GE.2) THEN
         IF (BANISO) WRITE (*,*) 'Es erfolgt anisotrope Verfeinerung'
      ENDIF
c
      READ (MDATA,*) IREFLV
      IF (MSHOW.GE.2) THEN
         IF (BANISO) WRITE (*,*) 'Ab Level ', IREFLV
      ENDIF
c
      READ (MDATA,*) OME1
      IF ((MSHOW.GE.2) .AND. (BANISO))
     *   WRITE (*,*) 'OME1    =  ', OME1
c
      READ (MDATA,*) OME2
      IF ((MSHOW.GE.2) .AND. (BANISO))
     *   WRITE (*,*) 'OME2    =  ', OME2
c
      READ (MDATA,*)
      READ (MDATA,*)
      READ (MDATA,*)
c
      READ (MDATA,*) IBLOCK
      IF (IBLOCK.gt.0) THEN
         IF (MSHOW.GE.2) 
     *   WRITE (*,*) 'IBLOCK    =  ', IBLOCK
        IF (ISL.EQ.3)  NSL=NSLOLD
         ISL=10!Ueberschreibt den Wert aus cc2d.dat
         write (*,*) 'nslold = ', nslold
      ENDIF 
c
      READ (MDATA,*) BLOCSL
      IF ((BLOCSL).AND.(IBLOCK.gt.0)) then
       WRITE (*,*) 'Auf dem Grobgitter wird exakt geloest'
       ISL=10 ! Ueberschreibt den Wert aus cc2d.dat
       NSL=1 !Bei ex. Invertierung brauche ich nur 1 Schritt
       ENDIF
c
      READ (MDATA,*) BSAVEM
      IF (BSAVEM) THEN
       IF (((BLOCSL).AND.(IBLOCK.gt.0)).OR.(ISL.EQ.3)) then
        WRITE (*,*) 'Die GG-Matrix wird Abgespeichert'
       ENDIF
      ENDIF
C
      READ (MDATA,*) BLSKSA
      IF (BLSKSA) THEN
        WRITE (*,*) 'Die Schurkomplementmatrix wird auf
     * Platte geschrieben'
      ENDIF
c
      READ (MDATA,*) IFACTO
      IF (IFACTO.LE.1) IFACTO=1
      IF (IFACTO.GE.2) IFACTO=2
      IF (MSHOW.GE.2) 
     *   WRITE (*,*) 'IFACTO    =  ', IFACTO
c
      READ (MDATA,*) IINVER
      IF (MSHOW.GE.2)      
     *   WRITE (*,*) 'IINVER    =  ', IINVER
c
      READ (MDATA,*)
      READ (MDATA,*)
      READ (MDATA,*)
c
      READ (MDATA,*) PSBARM
      IF (MSHOW.GE.2) 
     *   WRITE (*,*) 'PSBARM    =  ', PSBARM
c
      READ (MDATA,*) PSBARP
      IF (MSHOW.GE.2) 
     *   WRITE (*,*) 'PSBARP    =  ', PSBARP
C
      READ (MDATA,*) PSBARR
      IF (MSHOW.GE.2) 
     *   WRITE (*,*) 'PSBARR    =  ', PSBARR
c
      READ (MDATA,*) PSBARB
      IF (MSHOW.GE.2) 
     *   WRITE (*,*) 'PSBARB    =  ', PSBARB
C
      READ (MDATA,*) PSBARS
      IF (MSHOW.GE.2) 
     *   WRITE (*,*) 'PSBARS    =  ', PSBARS
c
      READ (MDATA,*) PSBVOL
      IF (MSHOW.GE.2) 
     *   WRITE (*,*) 'PSBVOL    =  ', PSBVOL
C
      READ (MDATA,*)
      READ (MDATA,*)
      READ (MDATA,*)
c
C
      READ (MDATA,*) NSBLSM
      IF (MSHOW.GE.2) 
     *   WRITE (*,*) 'NSBLSM    =  ', NSBLSM
C
      READ (MDATA,*) NSBLSL
      IF (MSHOW.GE.2) 
     *   WRITE (*,*) 'NSBLSL    =  ', NSBLSL
c
C=======================================================================
C  Jetzt die Daten fuer die modifizierte Rest/Prol
C=======================================================================
c
      READ (MDATA,*) 
      READ (MDATA,*)
      READ (MDATA,*)
c
      READ (MDATA,*) MODINT
      IF (MSHOW.GE.2) 
     *   WRITE (*,*) 'MODINT    =  ', MODINT
C
      READ (MDATA,*) IINTU
      IF (MSHOW.GE.2) 
     *   WRITE (*,*) 'IINTU    =  ', IINTU
C
      READ (MDATA,*) IINTP
      IF (MSHOW.GE.2) 
     *   WRITE (*,*) 'IINTP    =  ', IINTP
c
      END
C     #######################################################################
C  
C     Initialization of the velocity field
C
      SUBROUTINE UINIT(KVERT,KMID,DCORVG,U1,U2,NEL,NVT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
      DIMENSION KVERT(4,*),KMID(4,*),DCORVG(2,*),U1(*),U2(*)
C
      NVE=4
      NMT=NVT
C
      DO 10 IEL=1,NEL
      DO 20 IVE=1,NVE
C
      IMID=KMID(IVE,IEL)
      IF (IMID.LE.NMT) GOTO 20
C
      NMT =NMT+1
      IVEH=IVE+1
      IF (IVEH.EQ.NVE+1) IVEH=1
C
      IVT1=KVERT(IVE,IEL)
      PX1 =DCORVG(1,IVT1)
      PY1 =DCORVG(2,IVT1)
C
      IVT2=KVERT(IVEH,IEL)
      PX2 =DCORVG(1,IVT2)
      PY2 =DCORVG(2,IVT2)
C
      PX=0.5D0*(PX1+PX2)
      PY=0.5D0*(PY1+PY2)
      U1(IMID-NVT)=FDATIN(0,1,PX,PY,0D0,0D0)
      U2(IMID-NVT)=FDATIN(0,2,PX,PY,0D0,0D0)
C
20    CONTINUE
10    CONTINUE
C
      END

C
C     #######################################################################
C   

