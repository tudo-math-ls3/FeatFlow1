      PROGRAM TRIGEN3D
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      INCLUDE 'trigen3d.inc'
      PARAMETER (NNARR=299,NNLEV=9)
C
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CFILEI*60,CFILEO*60,CDATA*60
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
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NET,NAT,NVE,NEE,NAE,NVEL,NEEL,NVED,
     *                NVAR,NEAR,NBCT,NVBD,NEBD,NABD
      COMMON /TRIAA/  LCORVG,LCORMG,LCORAG,LVERT,LEDGE,LAREA,LADJ,
     *                LVEL,LEEL,LAEL,LVED,LAED,LVAR,LEAR,LEVE,LAVE,
     *                LNPR,LBCT,LVBD,LEBD,LABD
C
C *** EQUIVALENCE statement needed for DWORK concept
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
C *** Parametrization of the domain
      EXTERNAL PARX,PARY,PARZ
C *** Control of refinement - here: regular refinement
      EXTERNAL SEDB,SADB
C
C
C
C=======================================================================
C *** Initialization
C=======================================================================
      CALL ZINIT(NNWORK,'feat.msg',
     *           '#data/TRIGEN3D.ERR','#data/TRIGEN3D.PRT',
     *           '#data/TRIGEN3D.SYS','#data/TRIGEN3D.TRC') 
C
      MDATA=80
      CDATA='#data/trigen3d.dat'
      CALL  OF0 (MDATA,CDATA,1)
      CALL  GDAT(MDATA,IGMV,IAVS,NLEV,IFMT,CFILEI,CFILEO)
      CLOSE (MDATA)
C
C=======================================================================
C *** Read in coarse grid from file CFILEI
C=======================================================================
      MUNITT=50
      CALL XORSC(MUNITT,CFILEI)
      IF (IER.NE.0) GOTO 99998
      CLOSE(50)
C
C=======================================================================
C *** Refinements
C=======================================================================
C
      ISE=0
      ISA=2
      ISVEL=1
      ISEEL=0
      ISAEL=0
      ISVED=0
      ISAED=0
      ISVAR=0
      ISEAR=0
      ISEVE=0
      ISAVE=0
      ISVBD=0
      ISEBD=0
      ISABD=0
      IDISP=1
C
      IF (NLEV.GT.0) THEN
C
       DO 1000 ILEV=1,NLEV
C
       IF (ILEV.EQ.1) THEN
        CALL XSB0X(0,1,MAX(1,ISE),ISA,ISVEL,ISEEL,ISAEL,
     *             ISVED,ISAED,ISVAR,ISEAR,ISEVE,ISAVE,
     *             ISVBD,ISEBD,ISABD,IDISP,PARX,PARY,PARZ,
     *             SEDB,SADB)
       ELSE
        CALL XSB0X(1,0,MAX(1,ISE),ISA,ISVEL,ISEEL,ISAEL,
     *             ISVED,ISAED,ISVAR,ISEAR,ISEVE,ISAVE,
     *             ISVBD,ISEBD,ISABD,IDISP,PARX,PARY,PARZ,
     *             SEDB,SADB)      
       ENDIF   
C
************************************************************************
C
       CALL TRPARV(DWORK(L(LCORVG)),KWORK(L(LNPR)),KWORK(L(LVEL)),
     *             NVT,NVEL)
C
************************************************************************
C
       IF (ILEV.GE.2) THEN
        CALL CHCOOR(DWORK(L(LCORVG)),KWORK(L(LVERT)),KWORK(L(LAREA)),
     *              KWORK(L(LADJ)),KWORK(L(LNPR)),NEL,NVT)
        CALL TRPARV(DWORK(L(LCORVG)),KWORK(L(LNPR)),KWORK(L(LVEL)),
     *              NVT,NVEL)
       ENDIF
C
************************************************************************
C
       WRITE(MPROT,10000) ILEV,NEL,NVT,NAT,NET
       IF (M.GE.1)
     *     WRITE(6,10000) ILEV,NEL,NVT,NAT,NET
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
        CALL GMVTR(55,CTRIB(ILEV),NEL,NVT,KWORK(L(LVERT)),
     *       DWORK(L(LCORVG)))
        IF (IER.NE.0) GOTO 99998
        REWIND(55)
        CLOSE(55)
       ENDIF
C
C ***  Write refined mesh as AVS file
       IF ((IAVS.GT.0).AND.(ILEV.LE.IAVS)) THEN
        CALL AVSTR(56,CTRIA(ILEV),NEL,NVT,KWORK(L(LVERT)),
     *       DWORK(L(LCORVG)))
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
10000 FORMAT('LEVEL ',I2,' NEL ',I8,' NVT ',I8,' NAT ',I8,' NET ',I8)
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
      SUBROUTINE GDAT(MDATA,IGMV,IAVS,NLEV,IFMT,CFILEI,CFILEO)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER CFILEI*60,CFILEO*60
      PARAMETER (NNELM=1000)
C
C *** Standard COMMON blocks
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
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
C
C
      END
