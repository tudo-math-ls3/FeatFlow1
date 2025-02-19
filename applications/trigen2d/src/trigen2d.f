      PROGRAM TRIGEN2D
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      INCLUDE 'trigen2d.inc'
      PARAMETER (NNARR=299,NNLEV=9)
      PARAMETER (NNELM=1000)
C
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CFILEI*60,CFILEO*60,CDATA*60,CPARM*60
      CHARACTER CTRIA(NNLEV)*15,CTRIB(NNLEV)*15,CTRIF(NNLEV)*15
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
      COMMON /ELEMOD/ ITYPEL,NELMOD,DELMA,DELMB,KELMOD(NNELM)
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
     *           '#data/TRIGEN2D.ERR','#data/TRIGEN2D.PRT',
     *           '#data/TRIGEN2D.SYS','#data/TRIGEN2D.TRC') 
C
      MDATA=80
      CDATA='#data/trigen2d.dat'
      CALL  OF0 (MDATA,CDATA,1)
      CALL  GDAT(MDATA,IMESH,CPARM,IBDCHK,IGMV,IAVS,NLEV,IFMT,
     *           CFILEI,CFILEO)
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
C *** Prepare for mesh modifications
C=======================================================================
C
      IF ((ITYPEL.NE.0).AND.(NELMOD.NE.0)) THEN
       CALL ZNEW(NEL,3,LNELM,'KNELM ')
       IF (IER.NE.0) GOTO 99998
C
       DO 10 IELMOD=1,NELMOD
       IEL=KELMOD(IELMOD)
       IF (ITYPEL.EQ.1) ITYP= 1
       IF (ITYPEL.EQ.2) ITYP=-1
       IF (ITYPEL.EQ.3) ITYP=-1
       IF (ITYPEL.EQ.4) ITYP= 1
       IF (ITYPEL.EQ.5) ITYP= 10
       IF (ITYPEL.EQ.6) ITYP= 20
10     KWORK(L(LNELM)+IEL-1)=ITYP
C
      ENDIF
C
C=======================================================================
C *** Perform checks and find elements to modify
C=======================================================================
C
      IF ((ITYPEL.EQ.1).AND.(NELMOD.NE.0)) THEN
C
110    CALL CHEC1A(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),ICH1A)
       IF (IER.NE.0) GOTO 99998
       IF (ICH1A.GT.0) GOTO 110
C
       CALL CHEC1B(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),ICH1B)
       IF (IER.NE.0) GOTO 99998
       NFINEL=ICH1B
C
       IF (M.GE.1) THEN      
        WRITE(MPROT,*)         'TYP1 NEW ELEMENTS: ',NFINEL,ICH1A,ICH1B
        IF (M.GE.2) WRITE(6,*) 'TYP1 NEW ELEMENTS: ',NFINEL,ICH1A,ICH1B
        DO 120 IEL=1,NEL
        WRITE(MPROT,*)         IEL,KWORK(L(LNELM)+IEL-1)
        IF (M.GE.3) WRITE(6,*) IEL,KWORK(L(LNELM)+IEL-1)
120     CONTINUE
       ENDIF
C
      ENDIF
C
C
      IF ((ITYPEL.EQ.2).AND.(NELMOD.NE.0)) THEN
C
       CALL CHEC2A(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),ICH2A)
       IF (IER.NE.0) GOTO 99998
C
       CALL CHEC2B(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),ICH2B)
       IF (IER.NE.0) GOTO 99998
       NFINEL=ICH2B
C
       IF (M.GE.1) THEN      
        WRITE(MPROT,*)         'TYP2 NEW ELEMENTS: ',NFINEL,ICH2A,ICH2B
        IF (M.GE.2) WRITE(6,*) 'TYP2 NEW ELEMENTS: ',NFINEL,ICH2A,ICH2B
        DO 220 IEL=1,NEL
        WRITE(MPROT,*)         IEL,KWORK(L(LNELM)+IEL-1)
        IF (M.GE.3) WRITE(6,*) IEL,KWORK(L(LNELM)+IEL-1)
220     CONTINUE
       ENDIF
C
      ENDIF
C
C
      IF ((ITYPEL.EQ.3).AND.(NELMOD.NE.0)) THEN
C
       CALL CHEC3A(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),ICH3A)
       IF (IER.NE.0) GOTO 99998
C
       CALL CHEC3B(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),ICH3B)
       IF (IER.NE.0) GOTO 99998
       NFINEL=ICH3B
C
       IF (M.GE.1) THEN      
        WRITE(MPROT,*)         'TYP3 NEW ELEMENTS: ',NFINEL,ICH3A,ICH3B
        IF (M.GE.2) WRITE(6,*) 'TYP3 NEW ELEMENTS: ',NFINEL,ICH3A,ICH3B
        DO 330 IEL=1,NEL
        WRITE(MPROT,*)         IEL,KWORK(L(LNELM)+IEL-1)
        IF (M.GE.3) WRITE(6,*) IEL,KWORK(L(LNELM)+IEL-1)
330     CONTINUE
       ENDIF
C
      ENDIF
C
C
      IF ((ITYPEL.EQ.4).AND.(NELMOD.NE.0)) THEN
C
       CALL CHEC4A(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),ICH4A)
       IF (IER.NE.0) GOTO 99998
C
       CALL CHEC4B(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),ICH4B)
       IF (IER.NE.0) GOTO 99998
       NFINEL=ICH4B
C
       IF (M.GE.1) THEN      
        WRITE(MPROT,*)         'TYP4 NEW ELEMENTS: ',NFINEL,ICH4A,ICH4B
        IF (M.GE.2) WRITE(6,*) 'TYP4 NEW ELEMENTS: ',NFINEL,ICH4A,ICH4B
        DO 440 IEL=1,NEL
        WRITE(MPROT,*)         IEL,KWORK(L(LNELM)+IEL-1)
        IF (M.GE.3) WRITE(6,*) IEL,KWORK(L(LNELM)+IEL-1)
440     CONTINUE
       ENDIF
C
      ENDIF
C
C
      IF ((ITYPEL.EQ.5).AND.(NELMOD.NE.0)) THEN
C
       NFINEL=4*NELMOD
C
       IF (M.GE.1) THEN      
        WRITE(MPROT,*)         'TYP5 NEW ELEMENTS: ',NFINEL
        IF (M.GE.2) WRITE(6,*) 'TYP5 NEW ELEMENTS: ',NFINEL
        DO 550 IEL=1,NEL
        WRITE(MPROT,*)         IEL,KWORK(L(LNELM)+IEL-1)
        IF (M.GE.3) WRITE(6,*) IEL,KWORK(L(LNELM)+IEL-1)
550     CONTINUE
       ENDIF
C
      ENDIF
C
C
      IF ((ITYPEL.EQ.6).AND.(NELMOD.NE.0)) THEN
C
       CALL CHEC6A(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),ICH6A)
       IF (IER.NE.0) GOTO 99998
       NFINEL=2*ICH6A
C
       IF (M.GE.1) THEN      
        WRITE(MPROT,*)         'TYP6 NEW ELEMENTS: ',NFINEL,ICH6A
        IF (M.GE.2) WRITE(6,*) 'TYP6 NEW ELEMENTS: ',NFINEL,ICH6A
        DO 660 IEL=1,NEL
        WRITE(MPROT,*)         IEL,KWORK(L(LNELM)+IEL-1)
        IF (M.GE.3) WRITE(6,*) IEL,KWORK(L(LNELM)+IEL-1)
660     CONTINUE
       ENDIF
C
      ENDIF
C
C=======================================================================
C *** Generate new arrays
C=======================================================================
C
      CALL ZNEW(2*4*(NFINEL+NEL),1,LCORH,'DCORVG')
      IF (IER.NE.0) GOTO 99998
      CALL ZCPY(LCORVG,'DCORVG',LCORH,'DCORVG')
      IF (IER.NE.0) GOTO 99998
      CALL ZDISP(0,LCORVG,'DCORVG')      
      IF (IER.NE.0) GOTO 99998
      LCORVG=LCORH
C
      CALL ZNEW(4*(NFINEL+NEL),1,LVBDPH,'DVBDP ')
      IF (IER.NE.0) GOTO 99998
      CALL ZCPY(LVBDP,'DVBDP ',LVBDPH,'DVBDP ')
      IF (IER.NE.0) GOTO 99998
      CALL ZDISP(0,LVBDP,'DVBDP ')      
      IF (IER.NE.0) GOTO 99998
      LVBDP=LVBDPH
C
      CALL ZNEW(4*(NFINEL+NEL),3,LVBDH,'KVBD  ')
      IF (IER.NE.0) GOTO 99998
      CALL ZCPY(LVBD,'KVBD  ',LVBDH,'KVBD  ')
      IF (IER.NE.0) GOTO 99998
      CALL ZDISP(0,LVBD,'KVBD  ')      
      IF (IER.NE.0) GOTO 99998
      LVBD=LVBDH
C
      CALL ZNEW(4*(NFINEL+NEL),3,LNPRH,'KNPR  ')
      IF (IER.NE.0) GOTO 99998
      CALL ZCPY(LNPR,'KNPR  ',LNPRH,'KNPR  ')
      IF (IER.NE.0) GOTO 99998
      CALL ZDISP(0,LNPR,'KNPR  ')      
      IF (IER.NE.0) GOTO 99998
      LNPR=LNPRH
C
      CALL ZNEW(4*(NFINEL+NEL),3,LVERTH,'KVERT ')
      IF (IER.NE.0) GOTO 99998
      CALL ZCPY(LVERT,'KVERT ',LVERTH,'KVERT ')
      IF (IER.NE.0) GOTO 99998
      CALL ZDISP(0,LVERT,'KVERT ')      
      IF (IER.NE.0) GOTO 99998
      LVERT=LVERTH
C
      CALL ZNEW(4*(NFINEL+NEL),3,LADJH,'KADJ  ')
      IF (IER.NE.0) GOTO 99998
      CALL ZCPY(LADJ,'KADJ  ',LADJH,'KADJ  ')
      IF (IER.NE.0) GOTO 99998
      CALL ZDISP(0,LADJ,'KADJ  ')      
      IF (IER.NE.0) GOTO 99998
      LADJ=LADJH
C
      CALL ZNEW(NEL,3,LINDEL,'KINDEL')
      IF (IER.NE.0) GOTO 99998
C
      NELOLD=NEL
      NVTOLD=NVT
C
C=======================================================================
C *** perform element modifications
C=======================================================================
C
      IF ((ITYPEL.EQ.1).AND.(NELMOD.NE.0)) THEN
C
       CALL CREF1A(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),KWORK(L(LINDEL)),
     *             NELOLD,NVTOLD)
       IF (IER.NE.0) GOTO 99998
C
       CALL CREF1B(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),KWORK(L(LINDEL)),
     *             NELOLD,NVTOLD)
       IF (IER.NE.0) GOTO 99998
C
       CALL CREF1C(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),KWORK(L(LINDEL)),
     *             NELOLD,NVTOLD)
       IF (IER.NE.0) GOTO 99998
C
       CALL CREF1D(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),KWORK(L(LINDEL)),
     *             NELOLD,NVTOLD)
       IF (IER.NE.0) GOTO 99998
C
      ENDIF
C
C
C
      IF ((ITYPEL.EQ.2).AND.(NELMOD.NE.0)) THEN
C
       CALL CREF2A(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),KWORK(L(LINDEL)),
     *             NELOLD,NVTOLD)
       IF (IER.NE.0) GOTO 99998
C
       CALL CREF2B(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),KWORK(L(LINDEL)),
     *             NELOLD,NVTOLD)
       IF (IER.NE.0) GOTO 99998
C
       CALL CREF2D(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),KWORK(L(LINDEL)),
     *             NELOLD,NVTOLD)
       IF (IER.NE.0) GOTO 99998
C
       CALL CREF2E(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),KWORK(L(LINDEL)),
     *             NELOLD,NVTOLD)
       IF (IER.NE.0) GOTO 99998
C
      ENDIF
C
C
C
      IF ((ITYPEL.EQ.3).AND.(NELMOD.NE.0)) THEN
C
       CALL CREF3A(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),KWORK(L(LINDEL)),
     *             NELOLD,NVTOLD)
       IF (IER.NE.0) GOTO 99998
C
       CALL CREF3B(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),KWORK(L(LINDEL)),
     *             NELOLD,NVTOLD)
       IF (IER.NE.0) GOTO 99998
C
       CALL CREF3D(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),KWORK(L(LINDEL)),
     *             NELOLD,NVTOLD)
       IF (IER.NE.0) GOTO 99998
C
       CALL CREF3E(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),KWORK(L(LINDEL)),
     *             NELOLD,NVTOLD)
       IF (IER.NE.0) GOTO 99998
C
      ENDIF
C
C
C
      IF ((ITYPEL.EQ.4).AND.(NELMOD.NE.0)) THEN
C
       CALL CREF4A(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),KWORK(L(LINDEL)),
     *             NELOLD,NVTOLD)
       IF (IER.NE.0) GOTO 99998
C
       CALL CREF4B(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),KWORK(L(LINDEL)),
     *             NELOLD,NVTOLD)
       IF (IER.NE.0) GOTO 99998
C
      ENDIF
C
C
C
      IF ((ITYPEL.EQ.5).AND.(NELMOD.NE.0)) THEN
C
       CALL CREF5A(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),KWORK(L(LINDEL)),
     *             NELOLD,NVTOLD)
       IF (IER.NE.0) GOTO 99998
C
      ENDIF
C
C
C
      IF ((ITYPEL.EQ.6).AND.(NELMOD.NE.0)) THEN
C
       CALL CREF6A(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),KWORK(L(LNELM)),KWORK(L(LINDEL)),
     *             NELOLD,NVTOLD)
       IF (IER.NE.0) GOTO 99998
C
      ENDIF
C
C=======================================================================
C
      CALL ZDISP(0,LINDEL,'KINDEL')
      IF (IER.NE.0) GOTO 99998
C
      CALL ZDISP(2*NVT,LCORVG,'DCORVG')      
      IF (IER.NE.0) GOTO 99998
      CALL ZDISP(NVT,LNPR,'KNPR  ')      
      IF (IER.NE.0) GOTO 99998
C
C=======================================================================
C *** output of new coarse mesh
C=======================================================================
C
C *** Write in coarse grid as GMV file
      IF (IGMV.GT.0) THEN
       CALL GMVTR(52,'#gmv/coarse.gmv',NEL,NVT,
     *            KWORK(L(LVERT)),DWORK(L(LCORVG)))
       IF (IER.NE.0) GOTO 99998
      ENDIF
C
C *** Write in coarse grid as AVS file
      IF (IAVS.GT.0) THEN
       CALL AVSTR(53,'#avs/coarse.inp',NEL,NVT,
     *            KWORK(L(LVERT)),DWORK(L(LCORVG)))
       IF (IER.NE.0) GOTO 99998
      ENDIF
C
C *** Write in coarse grid from file CFILEO, opened as unit 51
      MUNITT=51
      CALL XOWSC1(MUNITT,CFILEO)
      IF (IER.NE.0) GOTO 99998
C
C=======================================================================
C *** perform boundary check for further refinements
C=======================================================================
C
      IF (IBDCHK.GT.0) THEN
C
       CALL SBC2P(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),NVBD)
       
       CALL ZDISP(0,LVBD,'KVBD  ')      
       IF (IER.NE.0) GOTO 99998
       CALL ZDISP(0,LEBD,'KEBD  ')      
       IF (IER.NE.0) GOTO 99998
       CALL ZDISP(0,LVBDP,'DVBDP ')      
       IF (IER.NE.0) GOTO 99998
C
       CALL XS2A
C
       CALL ZNEW(NVBD,-3,LVBD,'KVBD  ')
       IF (IER.NE.0) GOTO 99998
       CALL ZNEW(NVBD,-1,LVBDP,'DVBDP ')
       IF (IER.NE.0) GOTO 99998
       LMBDP=0
C
       CALL SBD01(DWORK(L(LCORVG)),KWORK(L(LNPR)),
     *            KWORK(L(LVBD)),DWORK(L(LVBDP)),DWORK(L(LMBDP)))
C
       NFINE=0
       IMID=0
       IADJ=1
       IVEL=0
       IDISP=1
       IBDP=2
       CALL XSB0X(NFINE,IMID,IADJ,IVEL,IDISP,IBDP, 
     *            S2DI0,S2DB0,PARX,PARY,TMAX)
       IF (IER.NE.0) GOTO 99998
C
       IBDCHH=IBDCHK
       CALL BDRCHK(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *             KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *             KWORK(L(LMM)),IBDCHH,PARX,PARY,TMAX)
       IF (IER.NE.0) GOTO 99998
C
      ENDIF
C
C=======================================================================
C *** output of modifications into new coarse mesh file
C=======================================================================
C
      IF ((ITYPEL.NE.0).AND.(NELMOD.NE.0)) THEN
       WRITE (MUNITT,*) '***'
       WRITE (MUNITT,*) '***'
       WRITE (MUNITT,*) '***'
       WRITE (MUNITT,*) '*** MESH GENERATED BY TRIADC'
       WRITE (MUNITT,*) '*** IMESH :',IMESH
       WRITE (MUNITT,*) '*** CFILEI:',CFILEI
       WRITE (MUNITT,*) '*** CFILEO:',CFILEO
       WRITE (MUNITT,*) '*** ITYPEL:',ITYPEL
       WRITE (MUNITT,*) '*** NELMOD:',NELMOD
       WRITE (MUNITT,*) '*** DELMA:',DELMA
       WRITE (MUNITT,*) '*** DELMB:',DELMB
C
       WRITE (MUNITT,*) '*** KELMOD'
       DO 900 IELMOD=1,NELMOD
900    WRITE (MUNITT,*) KELMOD(IELMOD)
C
       WRITE (MUNITT,*) '*** NFINEL:',NFINEL
       WRITE (MUNITT,*) '*** KNELM'
       DO 910 IEL=1,NELOLD
       INELM=KWORK(L(LNELM)+IEL-1)
       IF (INELM.NE.0) WRITE (MUNITT,*) IEL,INELM
910    CONTINUE
      ENDIF
C
C=======================================================================
C
      IF ((ITYPEL.NE.0).AND.(NELMOD.NE.0)) THEN
       CALL ZDISP(0,LNELM,'KNELM ')
       IF (IER.NE.0) GOTO 99998
      ENDIF
C
C=======================================================================
C *** further refinements
C=======================================================================
C
      IF (NLEV.GT.0) THEN
C
       WRITE(MPROT,10000)
       IF (M.GE.1) WRITE(6,10000)
C
       CALL SBC2P(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),NVBD)
       
       CALL ZDISP(0,LVBD,'KVBD  ')      
       IF (IER.NE.0) GOTO 99998
       CALL ZDISP(0,LEBD,'KEBD  ')      
       IF (IER.NE.0) GOTO 99998
       CALL ZDISP(0,LVBDP,'DVBDP ')      
       IF (IER.NE.0) GOTO 99998
C
       CALL XS2A
C
       CALL ZNEW(NVBD,-3,LVBD,'KVBD  ')
       IF (IER.NE.0) GOTO 99998
       CALL ZNEW(NVBD,-1,LVBDP,'DVBDP ')
       IF (IER.NE.0) GOTO 99998
       LMBDP=0
C
       CALL SBD01(DWORK(L(LCORVG)),KWORK(L(LNPR)),
     *            KWORK(L(LVBD)),DWORK(L(LVBDP)),DWORK(L(LMBDP)))
C
       IMID =1
       IADJ =1
       IVEL =0
       IDISP=1
       IBDP =2
C
       DO 1000 ILEV=1,NLEV
C
       IF (ILEV.EQ.1) THEN
        CALL XSB0X(0,IMID,IADJ,IVEL,IDISP,IBDP, 
     *             S2DI0,S2DB0,PARX,PARY,TMAX)
        IF (IER.NE.0) GOTO 99998
       ELSE
        CALL XSB0X(1,IMID,IADJ,IVEL,IDISP,IBDP, 
     *             S2DI0,S2DB0,PARX,PARY,TMAX)
        IF (IER.NE.0) GOTO 99998
       ENDIF
C
       WRITE(MPROT,10000) ILEV,NEL,NVT,NVBD
       IF (M.GE.1)
     *     WRITE(6,10000) ILEV,NEL,NVT,NVBD
C
C
      IF (ILEV.GE.2) 
     * CALL CHCOOR(DWORK(L(LCORVG)),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             KWORK(L(LADJ)),KWORK(L(LNPR)),NEL,NVT)
C
C=======================================================================
C *** boundary check of refined mesh
C=======================================================================
C
       IF (IBDCHK.GT.0) THEN
        IBDCHH=IBDCHK
        CALL BDRCHK(DWORK(L(LCORVG)),DWORK(L(LVBDP)),KWORK(L(LVBD)),
     *              KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *              KWORK(L(LMM)),IBDCHH,PARX,PARY,TMAX)
        IF (IER.NE.0) GOTO 99998
C
       ENDIF
C
C=======================================================================
C *** output of refined mesh
C=======================================================================
C ***  Write refined mesh in file CTRIF(ILEV)
       IF ((IFMT.NE.0).AND.(ILEV.LE.ABS(IFMT))) THEN
        IF (IFMT.GT.0) THEN
         CALL XOWS(54,CTRIF(ILEV),1)
        ELSE
         CALL XOWS(54,CTRIF(ILEV),0)
        ENDIF
        REWIND(54)
        CLOSE(54)
       ENDIF
C
C ***  Write refined mesh as GMV file
       IF ((IGMV.GT.0).AND.(ILEV.LE.IGMV)) THEN
        CALL GMVTR(55,CTRIB(ILEV),NEL,NVT,
     *             KWORK(L(LVERT)),DWORK(L(LCORVG)))
        IF (IER.NE.0) GOTO 99998
        REWIND(55)
        CLOSE(55)
       ENDIF
C
C ***  Write refined mesh as AVS file
       IF ((IAVS.GT.0).AND.(ILEV.LE.IAVS)) THEN
        CALL AVSTR(56,CTRIA(ILEV),NEL,NVT,
     *             KWORK(L(LVERT)),DWORK(L(LCORVG)))
        IF (IER.NE.0) GOTO 99998
        REWIND(56)
        CLOSE(56)
       ENDIF
C
1000   CONTINUE
C
      ENDIF
C
C *****************************************************************
C
      WRITE(MTERM,*) 'NWORK ',NWORK
      WRITE(MTERM,*) 'IWMAX ',IWMAX
C
      GOTO 99999
C
10000 FORMAT('LEVEL ',I2,' NEL ',I8,' NVT ',I8,' NVBD ',I6)
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
      SUBROUTINE GDAT(MDATA,IMESH,CPARM,IBDCHK,IGMV,IAVS,NLEV,IFMT,
     *                CFILEI,CFILEO)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER CFILEI*60,CFILEO*60,CPARM*60
      PARAMETER (NNELM=1000)
C
C *** Standard COMMON blocks
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /ELEMOD/ ITYPEL,NELMOD,DELMA,DELMB,KELMOD(NNELM)
      SAVE 
C
C *****************************************************************
C *** Input file
C *****************************************************************
C
      DO 10 IM=1,NNELM
10    KELMOD(IM)=0
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
      READ(MDATA,*) IBDCHK
      IBDCHK=ABS(IBDCHK)
C
      READ(MDATA,*) IGMV
      IGMV=ABS(IGMV)
C
      READ(MDATA,*) IAVS
      IAVS=ABS(IAVS)
C
      READ(MDATA,*) NLEV
      NLEV=ABS(NLEV)
C
      READ(MDATA,*) IFMT
C      
      READ(MDATA,*) CFILEI
C
      READ(MDATA,*) CFILEO
C
      READ(MDATA,*) ITYPEL
      ITYPEL=ABS(ITYPEL)
C
      READ(MDATA,*) NELMOD
      NELMOD=ABS(NELMOD)
C
      READ(MDATA,*) DELMA
      DELMA=ABS(DELMA)
C
      READ(MDATA,*) DELMB
      DELMB=ABS(DELMB)
C
C
      IF (ITYPEL.GT.6) THEN
       WRITE(6,*) 'WRONG TYPE !!!'
       ITYPEL=0
       NELMOD=0
      ENDIF
C
      IF ((ITYPEL.EQ.1).OR.(ITYPEL.EQ.4)
     *                 .OR.(ITYPEL.EQ.5).OR.(ITYPEL.EQ.6)) THEN
       IF (DELMA.LT.1D-5   ) DELMA=1D-5
       IF (DELMA.GT.0.490D0) DELMA=0.49D0
       IF (DELMB.LT.1D-3   ) DELMB=1D-3
       IF (DELMB.GT.0.495D0) DELMB=0.495D0
       DELMB=1D0-DELMA
      ENDIF
C
      IF ((ITYPEL.EQ.2).OR.(ITYPEL.EQ.3)) THEN
       IF (DELMA.LT.1D-5   ) DELMA=1D-5
       IF (DELMA.GT.0.990D0) DELMA=0.99D0
       IF (DELMB.LT.1D-3   ) DELMB=1D-3
       IF (DELMB.GT.0.995D0) DELMB=0.995D0
      ENDIF
C
      DO 100 IELMOD=1,NELMOD
      READ(MDATA,*) KELMOD(IELMOD)
      KELMOD(IELMOD)=ABS(KELMOD(IELMOD))
100   CONTINUE
C
C
C
      END
