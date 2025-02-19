************************************************************************
      PROGRAM  CC3DI 
************************************************************************
*
*   Purpose: - estimator for storage amount for cc3d
*
************************************************************************
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER CDATA*60
      CHARACTER CPARM1*60,CMESH1*60,CFILE1*60,CSTART*60,CSOL*60
C
C
C=======================================================================
C     MTERM,MKEYB machine dependent
C=======================================================================
C
      MTERM=6
      MKEYB=5
C
      MDATA=79
      CDATA='#data/cc3d.dat'
C
      WRITE(MTERM,*) ' Input: Number of coarse mesh elements? '
      READ (MKEYB,*)   NELC
C
      OPEN (UNIT=MDATA,FILE=CDATA)
C
C=======================================================================
C     Read #data/pp3d.dat
C=======================================================================
C
      READ(MDATA,*)
      READ(MDATA,*)
      READ(MDATA,*) 
      READ(MDATA,*)
      READ(MDATA,*)
C
      READ(MDATA,*) IMESH1
      READ(MDATA,*) IRMESH
      READ(MDATA,*) CPARM1
      READ(MDATA,*) CMESH1
      READ(MDATA,*) CFILE1
      READ(MDATA,*) ISTART
      READ(MDATA,*) CSTART
      READ(MDATA,*) ISOL
      READ(MDATA,*) CSOL
C
      READ(MDATA,*) 
      READ(MDATA,*)
      READ(MDATA,*)
C
      READ(MDATA,*) M
      READ(MDATA,*) MT
      READ(MDATA,*) ICHECK
      READ(MDATA,*) MSHOW
C
      READ(MDATA,*) 
      READ(MDATA,*)
      READ(MDATA,*)
C
      READ(MDATA,*)  NLMIN
      NLMIN=ABS(NLMIN)
      IF (NLMIN.EQ.0) NLMIN=1
C
      READ(MDATA,*)  NLMAX
      NLMAX=ABS(NLMAX)
      IF (NLMAX.LT.NLMIN) NLMAX=NLMIN
      NLEV=NLMAX-NLMIN+1
C
      READ(MDATA,*) IELT
      READ(MDATA,*) ISTOK
C
      READ(MDATA,*) IRHS
      IF ((IRHS.LT.0).OR.(IRHS.GT.2)) IRHS=2
      READ(MDATA,*) IBDR
      IF ((IBDR.LT.0).OR.(IBDR.GT.2)) IBDR=2
C
      READ(MDATA,*) IERANA
C
      READ(MDATA,*) IMASS
      IF (IMASS.NE.1) IMASS=0
      READ(MDATA,*) IMASSL
      IF (IMASSL.NE.1) IMASSL=0
C
      READ(MDATA,*) IUPW
      IF ((IUPW.GT.1).OR.(IUPW.LT.0)) IUPW=1
C
      READ(MDATA,*) IPRECA
      IF ((IPRECA.LT.0).OR.(IPRECA.GT.4)) IPRECA=1
      IF ((IPRECA.EQ.4).AND.(IUPW.EQ.1))  IPRECA=1
      IF ((IPRECA.EQ.4).AND.(ISTOK.EQ.1)) IPRECA=1
      IF (IPRECA.EQ.3)                    IPRECA=2
      READ(MDATA,*) IPRECB
      IF ((IPRECB.LT.0).OR.(IPRECB.GT.4)) IPRECB=2
      IF (IPRECB.EQ.1) IPRECB=0
      IF (IPRECB.EQ.2) IPRECB=3
      IF (IPRECB.EQ.4) IPRECB=3
C
      READ(MDATA,*) ICUBM
      READ(MDATA,*) ICUBA
      READ(MDATA,*) ICUBN
      READ(MDATA,*) ICUBB
      READ(MDATA,*) ICUBF
C
      READ(MDATA,*) INLMIN
      IF (INLMIN.LT.-1) INLMIN=-1
      READ(MDATA,*) INLMAX
      IF (INLMAX.LT.-1) INLMAX=-1
      IF (INLMAX.LT.INLMIN) INLMAX=INLMIN
C
      ITEXL=0
      IF ((INLMIN.EQ.INLMAX).AND.(INLMIN.EQ. 1)) THEN
       ITEXL=1
      ENDIF
C
      READ(MDATA,*) ICYCLE
      READ(MDATA,*) ILMIN
      READ(MDATA,*) ILMAX
      READ(MDATA,*) IINT
      READ(MDATA,*) ISM
      READ(MDATA,*) ISL
      READ(MDATA,*) NSM
      READ(MDATA,*) NSL
      READ(MDATA,*) NSMFAC
C
      READ(MDATA,*) 
      READ(MDATA,*)
      READ(MDATA,*)
C
      READ(MDATA,*) RE
      READ(MDATA,*) UPSAM
C
      READ(MDATA,*) OMGMIN
      READ(MDATA,*) OMGMAX
      IF (OMGMAX.LT.OMGMIN) OMGMAX=OMGMIN
C
      READ(MDATA,*) OMGINI
      READ(MDATA,*) EPSD
      READ(MDATA,*) EPSDIV
      READ(MDATA,*) EPSUR
      READ(MDATA,*) EPSPR
      READ(MDATA,*) DMPD
      READ(MDATA,*) DMPMG
      READ(MDATA,*) EPSMG
      READ(MDATA,*) DMPSL
      READ(MDATA,*) EPSSL
      READ(MDATA,*) RLXSM
      READ(MDATA,*) RLXSL
      READ(MDATA,*) AMINMG
      READ(MDATA,*) AMAXMG
C
      READ(MDATA,*) 
      READ(MDATA,*)
      READ(MDATA,*)
C
      READ(MDATA,*) ISTAT
      READ(MDATA,*) NITNS
      READ(MDATA,*) EPSNS
      READ(MDATA,*) TIMENS
      READ(MDATA,*) THETA
      READ(MDATA,*) TSTEP
      READ(MDATA,*) IFRSTP
      READ(MDATA,*) INSAV
      READ(MDATA,*) INSAVN
      READ(MDATA,*) DTFILM
      READ(MDATA,*) DTAVS
      READ(MDATA,*) DTGMV
      READ(MDATA,*) IFUSAV
      READ(MDATA,*) IFPSAV
      READ(MDATA,*) IFXSAV
      READ(MDATA,*) IAVS
      READ(MDATA,*) IGMV
      READ(MDATA,*) IFINIT
C
      READ(MDATA,*) IADTIM
      IF (ABS(IADTIM).GT.3) IADTIM=0
C
      READ(MDATA,*) TIMEMX
      READ(MDATA,*) DTMIN
      READ(MDATA,*) DTMAX
      READ(MDATA,*) DTFACT
      READ(MDATA,*) TIMEIN
      READ(MDATA,*) EPSADI
      READ(MDATA,*) EPSADL
      READ(MDATA,*) EPSADU
      READ(MDATA,*) IEPSAD
      READ(MDATA,*) IADIN
      READ(MDATA,*) IREPIT
C
C
      CLOSE(MDATA)
C
C=======================================================================
C     Statistics
C=======================================================================
C      
      WRITE (MTERM,1001) 'Number of levels           : ',NLEV
      NFAC=INT(8D0**(NLEV-1))
      WRITE (MTERM,1001)
C
      WRITE (MTERM,1001) 'Number of coarse mesh cells: ',NELC
      NFACEL=NELC*NFAC
      WRITE (MTERM,1001) 'Number of   fine mesh cells: ',NFACEL
      WRITE (MTERM,1001)
C
      NWORKF=102
C
      IOMEGA=0
      IF ((OMGMIN.GT.0D0).OR.(OMGMAX.GT.0D0)) IOMEGA=1
C
      IF (ISTAT.EQ.0) THEN
       IMASS=0
       IRHS=0
       IADTIM=0
       ITEXL=0
       IF (IBDR.GE.2) IBDR=1
      ENDIF
C
      WRITE (MTERM,1001) 'ISTAT : ',ISTAT 
      WRITE (MTERM,1001) 'IMASS : ',IMASS 
      WRITE (MTERM,1001) 'IUPW  : ',IUPW 
      WRITE (MTERM,1001) 'IPRECA: ',IPRECA
      WRITE (MTERM,1001) 'IPRECB: ',IPRECB
      WRITE (MTERM,1001) 'IOMEGA: ',IOMEGA
      WRITE (MTERM,1001) 'IBDR  : ',IBDR
      WRITE (MTERM,1001) 'IRHS  : ',IRHS
      WRITE (MTERM,1001) 'IADTIM: ',IADTIM
      WRITE (MTERM,1001) 'ITEXL : ',ITEXL
C
C
      IF (IUPW.EQ.0) NWORKF=NWORKF+1
C
      IF (ISTAT.EQ.1) THEN
       IF ((IMASS.EQ.1).AND.(IPRECA.LT.2)) THEN
        IF (IPRECA.EQ.1) NWORKF=NWORKF+38
        IF (IPRECA.EQ.0) NWORKF=NWORKF+19
       ENDIF
       IF ((IMASS.EQ.0).AND.(IPRECA.LT.2)) THEN
        IF (IPRECA.EQ.1) NWORKF=NWORKF+3
        IF (IPRECA.EQ.0) NWORKF=NWORKF+2
       ENDIF
      ENDIF
C
      IF (IPRECA.LT.2) THEN
       IF (IPRECA.EQ.1) NWORKF=NWORKF+38
       IF (IPRECA.EQ.0) NWORKF=NWORKF+19
      ENDIF
C
      IF (IPRECB.NE.2) THEN
       IF ((IPRECB.EQ.1).OR.(IPRECB.EQ.4)) NWORKF=NWORKF+25
       IF ((IPRECB.EQ.0).OR.(IPRECB.EQ.3)) NWORKF=NWORKF+15
      ENDIF
C
      IF (IOMEGA.EQ.1) NWORKF=NWORKF+10
C
      IF (IBDR  .EQ.2) NWORKF=NWORKF+2
C
      IF (IRHS  .GE.1) NWORKF=NWORKF+10
C
      IF (ISTAT.EQ.1) NWORKF=NWORKF+10
C
      IF (IADTIM.NE.0) NWORKF=NWORKF+5
C
      IF (ITEXL .NE.0) THEN
       NWORKF=NWORKF+10
       IF (ABS(IADTIM).EQ.2) NWORKF=NWORKF+10
      ENDIF
C
C
C
      WRITE (MTERM,1001) 
      WRITE (MTERM,1002) 'FACTOR: ',NWORKF,DBLE(NWORKF)/128D0
      WRITE (MTERM,1001) 'NNWORK: ',NWORKF*NFACEL
      WRITE (MTERM,1001) 
C
      WRITE (MTERM,1003) 'MBYTE : ',DBLE(8*NWORKF*NFACEL)/1D6
C
C
C
1001  FORMAT(A40,3X,I10)
1002  FORMAT(A40,3X,I10,3X,1F12.5)
1003  FORMAT(A40,3X,1F12.5)
C
C
C
      END
