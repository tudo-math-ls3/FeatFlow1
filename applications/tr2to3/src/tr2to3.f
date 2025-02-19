      PROGRAM TR2TO3
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      INCLUDE 'tr2to3.inc'
      PARAMETER (NNARR=299,NNLEV=9,NNPZ=100)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CFILEI*60,CFILEO*60,CDATA*60,CPARM*60
      CHARACTER CTRIA(NNLEV)*15,CTRIB(NNLEV)*15,CTRIF(NNLEV)*15
C
      DIMENSION PZARR(NNPZ)
C
      DATA CTRIA/ '#avs/tria1.inp ','#avs/tria2.inp ','#avs/tria3.inp ',
     *            '#avs/tria4.inp ','#avs/tria5.inp ','#avs/tria6.inp ',
     *            '#avs/tria7.inp ','#avs/tria8.inp ','#avs/tria9.inp '/
      DATA CTRIB/ '#gmv/tria1.gmv ','#gmv/tria2.gmv ','#gmv/tria3.gmv ',
     *            '#gmv/tria4.gmv ','#gmv/tria5.gmv ','#gmv/tria6.gmv ',
     *            '#gmv/tria7.gmv ','#gmv/tria8.gmv ','#gmv/tria9.gmv '/
      DATA CTRIF/ '#tries/#tria1  ','#tries/#tria2  ','#tries/#tria3  ',
     *            '#tries/#tria4  ','#tries/#tria5  ','#tries/#tria6  ',
     *            '#tries/#tria7  ','#tries/#tria8  ','#tries/#tria9  '/
C
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
C
C *** EQUIVALENCE statement needed for DWORK concept
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
C *** Parametrization of the domain
      EXTERNAL PARX,PARY,TMAX
C *** Control of refinement - here, regular refinement is used
      EXTERNAL S2DI0,S2DB0
C
C
C
C=======================================================================
C *** Initialization
C=======================================================================
      CALL ZINIT(NNWORK,'feat.msg',
     *           '#data/TRI2TO3.ERR','#data/TRI2TO3.PRT',
     *           '#data/TRI2TO3.SYS','#data/TRI2TO3.TRC') 
C
      MDATA=80
      CDATA='#data/tr2to3.dat'
      CALL  OF0 (MDATA,CDATA,1)
      CALL  GDAT(MDATA,IMESH,CPARM,IGMV,IAVS,CFILEI,CFILEO,NPZ,
     *           PZMIN,PZMAX,PZARR)
      CLOSE (MDATA)
C
C=======================================================================
C *** Data preparation for OMEGA
C=======================================================================
      IF (IMESH.EQ.1) THEN
       CALL RDPARM (CPARM,MDATA)
       CLOSE(MDATA)
      ENDIF
C
C=======================================================================
C *** Read in coarse grid from file CFILEI
C=======================================================================
      MUNITT=50
      CALL XORSC(MUNITT,CFILEI)
      IF (IER.NE.0) GOTO 99998
      CLOSE(50)
C
      NFINE=0
      IMID=0
      IADJ=1
      IVEL=0
      IDISP=1
      IBDP=2
      CALL XSB0X(NFINE,IMID,IADJ,IVEL,IDISP,IBDP, 
     *           S2DI0,S2DB0,PARX,PARY,TMAX)
      IF (IER.NE.0) GOTO 99998
C
C=======================================================================
C *** Save old arrays and constants
C=======================================================================
C
      NELOLD=NEL
      NVTOLD=NVT
C
      CALL ZNEW(2*NVT,1,LCORH,'DCORH')
      IF (IER.NE.0) GOTO 99998
      CALL ZCPY(LCORVG,'DCORVG',LCORH,'DCORH ')
      IF (IER.NE.0) GOTO 99998
      CALL ZDISP(0,LCORVG,'DCORVG')      
      IF (IER.NE.0) GOTO 99998
C
      CALL ZNEW(NVT,3,LNPRH,'KNPRH ')
      IF (IER.NE.0) GOTO 99998
      CALL ZCPY(LNPR,'KNPR  ',LNPRH,'KNPRH ')
      IF (IER.NE.0) GOTO 99998
      CALL ZDISP(0,LNPR,'KNPR  ')      
      IF (IER.NE.0) GOTO 99998
C
      CALL ZNEW(4*NEL,3,LVERTH,'KVERTH')
      IF (IER.NE.0) GOTO 99998
      CALL ZCPY(LVERT,'KVERT ',LVERTH,'KVERTH')
      IF (IER.NE.0) GOTO 99998
      CALL ZDISP(0,LVERT,'KVERT ')      
      IF (IER.NE.0) GOTO 99998
C
C=======================================================================
C *** Generate new arrays and constants
C=======================================================================
C
      NEL=   NPZ *NELOLD
      NVT=(NPZ+1)*NVTOLD
C
      CALL ZNEW(3*NVT,1,LCORVG,'DCORVG')
      IF (IER.NE.0) GOTO 99998
C
      CALL ZNEW(NVT,3,LNPR,'KNPR  ')
      IF (IER.NE.0) GOTO 99998
C
      CALL ZNEW(8*NEL,3,LVERT,'KVERT ')
      IF (IER.NE.0) GOTO 99998
C
C=======================================================================
C
      CALL GENC3(DWORK(L(LCORVG)),DWORK(L(LCORH)),KWORK(L(LVERT)),
     *           KWORK(L(LVERTH)),KWORK(L(LNPR)),KWORK(L(LNPRH)),
     *           NEL,NELOLD,NVT,NVTOLD,NPZ,PZMIN,PZMAX,PZARR)

C
      WRITE(MTERM,*) 'Coarse elements 2d ',NELOLD
      WRITE(MTERM,*) 'Coarse elements 3d ',NEL
      WRITE(MTERM,*)
C
C=======================================================================
C *** Write in coarse grid from file CFILEO, opened as unit 51
C=======================================================================
      MUNITT=51
      CALL XOWSC3(MUNITT,CFILEO)
      IF (IER.NE.0) GOTO 99998
C
C=======================================================================
C ***  graphical output of new coarse mesh
C=======================================================================
C
C *** Write in coarse grid as GMV file
      IF (IGMV.GT.0) THEN
       CALL GMVTR(52,'#gmv/c3d.gmv',NEL,NVT,KWORK(L(LVERT)),
     *            DWORK(L(LCORVG)))
       IF (IER.NE.0) GOTO 99998
      ENDIF
C
C *** Write in coarse grid as AVS file
      IF (IAVS.GT.0) THEN
       CALL AVSTR(53,'#avs/c3d.inp',NEL,NVT,KWORK(L(LVERT)),
     *            DWORK(L(LCORVG)))
       IF (IER.NE.0) GOTO 99998
      ENDIF
C
C *****************************************************************
C
      WRITE(MTERM,*) 'NWORK ',NWORK
      WRITE(MTERM,*) 'IWMAX ',IWMAX
C
      GOTO 99999
C
C
C
99998 WRITE(MTERM,*) 'IER', IER
      WRITE(MTERM,*) 'IN SUBROUTINE ',SUB
C
99999 END
C
C
C-----------------------------------------------------------------------
C
C
      SUBROUTINE GDAT(MDATA,IMESH,CPARM,IGMV,IAVS,CFILEI,CFILEO,NPZ,
     *                PZMIN,PZMAX,PZARR)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER CFILEI*60,CFILEO*60,CPARM*60
C
      DIMENSION PZARR(*)
C
C *** Standard COMMON blocks
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      SAVE 
C
C *****************************************************************
C *** Input file
C *****************************************************************
C
      READ(MDATA,*)
      READ(MDATA,*)
      READ(MDATA,*)
      READ(MDATA,*) M
      READ(MDATA,*) MT
      READ(MDATA,*) ICHECK
C
      READ(MDATA,*) IMESH
      IF (IMESH.NE.1) IMESH=0
C      
      READ(MDATA,*) CPARM
C
      READ(MDATA,*) IGMV
      IGMV=ABS(IGMV)
C
      READ(MDATA,*) IAVS
      IAVS=ABS(IAVS)
C
      READ(MDATA,*) CFILEI
C
      READ(MDATA,*) CFILEO
C
      READ(MDATA,*) NPZ
      NPZ=ABS(NPZ)
C
      READ(MDATA,*) PZMIN
C
      READ(MDATA,*) PZMAX
      IF (PZMAX.LE.PZMIN) THEN
       WRITE(MTERM,*) 'WRONG VALUE FOR PZMAX'
       STOP
      ENDIF
C
      DO 10 IPZ=1,NPZ-1
10    READ(MDATA,*) PZARR(IPZ)
C
C
      PZH=PZMIN
      DO 20 IPZ=1,NPZ-1
      PZ=PZARR(IPZ)
      IF ((PZ.LE.PZMIN).OR.(PZ.GE.PZMAX).OR.(PZ.LE.PZH)) THEN
       WRITE(MTERM,*) 'WRONG VALUE FOR PZARR',IPZ
       STOP
      ENDIF
      PZH=PZ
20    CONTINUE
C
C
C
      END
