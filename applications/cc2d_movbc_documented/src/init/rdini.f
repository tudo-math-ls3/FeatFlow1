************************************************************************
* RDOUT
*
* Read the parameters from the .DAT file concerning the OUTPUT
* parameters. This is necessary on program start to prepare all the
* output. It will open the file CFNAME, read the parameters and
* close the file. 
*
* In:
*   MDATA  - Unit number to use for reading process
*   CFNAME - name of the file
*   MFILE  - A unit filename for output to a file
*
* Out:
*   MSHOW  - Message level for initialisation routines (from file).
*            Controls if the initialisation routines print out the
*            readed parameters to screen
*
* Modifies the parameters in the /OUTPUT/ COMMON according to
* the .DAT file. Opens the file unit MFILE for user output to a file.
* The filename that should be used for file output is taken from the
* DAT file. The unit number is stored in the COMMON block as well.
************************************************************************

      SUBROUTINE RDOUT (MDATA,CFNAME,MFILE,MSHOW)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cfiles.inc'
      
C parameters
      
      CHARACTER CFNAME*(*)
      INTEGER MSHOW,MDATA,MFILE
      
C local variables

      INTEGER I,IFMTS

C Open DAT file
      
      IFMTS = 1
      CALL  OF0 (MDATA,CFNAME,IFMTS)
      
C Ignore all lines up to the name of the protocol file variables

      DO I=1,9
        READ (MDATA,*)
      END DO
      
C Read the name of protocol file

      READ (MDATA,*) CFILE1

C Ignore all lines up to the /OUTPUT/ variables

      DO I=11,17
        READ (MDATA,*)
      END DO

C This data must be read in directly without GETINT.
      
      READ(MDATA,*) M
      M=ABS(M)

      READ(MDATA,*) MT
      MT=ABS(MT)

      READ(MDATA,*) ICHECK
      ICHECK=ABS(ICHECK)

      READ(MDATA,*) MSHOW
      MSHOW=ABS(MSHOW)
      
      IF (MSHOW.GE.2) THEN
        WRITE(MTERM,*) 'Message level for file output:      M      = ',M
        WRITE(MTERM,*) 'Message level for terminal output:  MT     = ',
     *                  MT
        WRITE(MTERM,*) 'Level for tracing:                  ICHECK = ',
     *                  ICHECK
        WRITE(MTERM,*) 'Output level in initialisation:     MSHOW  = ',
     *                  MSHOW
      END IF
      
      CLOSE (MDATA)

C Open file for user output

      MFILE1 = MFILE
      CALL  OF0 (MFILE,CFILE1,IFMTS)
      
      END


************************************************************************
* RDDAT
*
* Read the parameters in the .DAT file and store the information in
* the variables of the Common blocks /IPARAM/ and /RPARAM/.
*
* This routine is used in the routine INIT1 to parse the
* initialization file. It will open the file CFNAME, read the parameters 
* and close the file. 
*
* The routein RDOUT has to be called previously to this!
*
* In:
*  MDATA  - Handle/IO number to the .DAT file to use. 
*  CFNAME - name of the file
*  MSHOW  - Controls the output of the subroutine GDAT. The higher the
*           number, the more output. 0 gives basic output, 2 gives
*           most information, also about the content of the parameters.
* 
* Out:
*  IRMESH - Parameter IMESH in the .DAT file. This information is not
*           stored in the Common blocks for now.
************************************************************************

************************************************************************
      SUBROUTINE RDDAT (MDATA,CFNAME,MSHOW)
************************************************************************
      
      IMPLICIT NONE

C *** Output/error variables
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'

C *** multigrid data
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgtria.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      
C *** variables for Navier Stokes equations
      INCLUDE 'cns.inc'
      INCLUDE 'cnsparfrac.inc'
      
C *** input/output files
      INCLUDE 'cfiles.inc'      

C *** user COMMON blocks
      INCLUDE 'cinidat.inc'
      INCLUDE 'cinidat2.inc'

C parameter

      INTEGER MDATA, MSHOW
      INTEGER MFILE
      CHARACTER*(*) CFNAME

C local variables

      INTEGER II,IFMTS

C-----------------------------------------------------------------------
C *** Input file
C-----------------------------------------------------------------------

C Open the file

      IFMTS = 1
      CALL  OF0 (MDATA,CFNAME,IFMTS)

C Read it

   1  FORMAT(80('-'))
   
C Some lines to ignore
   
      READ(MDATA,*)
      READ(MDATA,*)
      READ(MDATA,*) 
      READ(MDATA,*)
      READ(MDATA,*)
C
      READ(MDATA,*) IMESH1
      IF (IMESH1.NE.1) IMESH1=0
      READ(MDATA,*) IRMESH
      IRMESH=ABS(IRMESH)
      READ(MDATA,*) CPARM1
      READ(MDATA,*) CMESH1

      IF (MSHOW.GE.2) WRITE(MTERM,*) 'Parametrization file = ',CPARM1
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'Coarse grid file     = ',CMESH1

C Filename for file output was already opened earlier in RDOUT, so here 
C ignore it

      READ(MDATA,*) 

      IF (MSHOW.GE.2) WRITE(MTERM,*) 'MFILE,CFILE: ', MFILE1,CFILE1
      MFILE=MFILE1

      READ(MDATA,*) ISTART
      IF (.NOT.(((ABS(ISTART).GE.0).AND.(ABS(ISTART).LE.NNLEV)).OR.
     *     ((ABS(ISTART).GT.100).AND.(ABS(ISTART).LE.(100+NNLEV)))) )
     *    ISTART=0
      READ(MDATA,*) CSTART
      MSTART=63
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*)'ISTART,MSTART,CSTART: ', ISTART,MSTART,CSTART

      READ(MDATA,*) ISOL
      IF (ABS(ISOL).GT.1) ISOL=1
      READ(MDATA,*) CSOL
      MSOL=64
      IF (MSHOW.GE.2) WRITE(MTERM,*)'ISOL,MSOL,CSOL: ', ISOL,MSOL,CSOL

      IF (MSHOW.GE.0) WRITE(MFILE,1)
      IF (MSHOW.GE.0) WRITE(MFILE,*) '        INPUT DATA'
      IF (MSHOW.GE.0) WRITE(MFILE,1)
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'Parametrization file = ',CPARM1
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'Coarse grid file     = ',CMESH1
C
C-----------------------------------------------------------------------
C *** Values for /OUTPUT/
C-----------------------------------------------------------------------

C Can be ignored. Had been read earlier in RDOUT

      READ(MDATA,*) 
      READ(MDATA,*)
      READ(MDATA,*)
C M
      READ(MDATA,*) 
C MT
      READ(MDATA,*) 
C ICHECK
      READ(MDATA,*) 
C MSHOW
      READ(MDATA,*) 
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
      READ(MDATA,*) IPRECA
      IF ((IPRECA.LT.0).OR.(IPRECA.GT.4)) IPRECA=1
      IF ((IPRECA.EQ.4).AND.(IUPW.EQ.1))  IPRECA=1
      IF ((IPRECA.EQ.4).AND.(ISTOK.EQ.1)) IPRECA=1
      IF (IPRECA.EQ.3)                    IPRECA=2
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
      IF ((ABS(IINT).LT.1).OR.(ABS(IINT).GT.4)) IINT=2
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'type of interpolation:  IINT = ', IINT
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'type of interpolation:  IINT = ', IINT

C *MK*: extended prolongation/restriction/matrix generation
      READ(MDATA,*)  IAVPR
      IF ((IAVPR.LT.0)) IAVPR=3
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'type of averaging in pr/rest: IAVPR  = ', IAVPR
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'type of averaging in pr/rest: IAVPR  = ', IAVPR

      READ(MDATA,*)  IAPRM
      IF ((IAPRM.LT.0)) IAPRM=3
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'type of averaging:             IAPRM  = ', IAPRM
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'type of averaging:             IAPRM  = ', IAPRM
C [*MK*] 

      READ(MDATA,*) ISM
      IF (ISM.NE.1) ISM=1
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'type of smoother :  ISM = ',ISM
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'type of smoother :  ISM = ',ISM
C
      READ(MDATA,*) ISL
      IF (ISL.NE.1) ISL=1
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

C *MK*: extended prolongation/restriction/matrix generation

      READ(MDATA,*) DPREP
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'treshold param. pr/rest:       DPREP  = ', DPREP
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'treshold param. pr/rest:       DPREP  = ', DPREP 

      READ(MDATA,*) DMTEP
      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'treshold param. mat. gen.:     DMTEP  = ', DMTEP 
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'treshold param. mat. gen.:     DMTEP  = ', DMTEP

C [*MK*]

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
     
      IF (ISTAT.EQ.0) THEN
       IAUSAV=INSAV
      ELSE
       IAUSAV=0
      ENDIF

C And finally close the DAT file we read from

      CLOSE (MDATA)
      
      END
      
************************************************************************
* Initialise Navier Stokes time stepping scheme
*
* This routine initialises the variables that are used by the
* (instationary) Navier-Stokes solver for time stepping.
* It will modify the necessary common blocks directly according to
* the parameters in the DAT file.
*
* The DAT file must have been read in previously!
************************************************************************
      
      SUBROUTINE INNSTT
      
      IMPLICIT NONE

      INCLUDE 'cout.inc'
      INCLUDE 'cfiles.inc'

      INCLUDE 'cns.inc'
      INCLUDE 'cnsparfrac.inc'
      
      INCLUDE 'cinidat.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'

C     ITEXL - ITeration EXtrapolation Linear.
C             Status flag during the calculation.
C     This flag indicates during the calculation whether the
C     linear extrapolation scheme is active or not. For a documentation
C     about how this works, see the main algorithm.
C
C     Linear extrapolation replaces nonlinear iteration.
C     It's active as soon as INLMIN=INLMAX=1, as then there is
C     no nonlinear iteration. We indicate this to the main
C     algorithm by setting ITEXL to 1.
C
C     ITEXL describes a status machine that corresponds to the
C     content of LTML as well as TIML11,TIML12,TIML31,TIML32,
C     and defines the behaviour of the extrapolation handling.
C     The values are defined as following:
C      = -1 : Predictor-step with no substeps, i.e. time-step
C             with step-size 3xTSTEP. Use TIM1x for adaptive time-
C             stepping.
C             No information about previous time-step available,
C             i.e. LTML does not define u^(n-1)
C      =  1 : Anywhere during the simulation, predictor-step
C             with step-size 3xTSTEP. Use TIM1x for adaptive time-
C             stepping.
C             LTML identifies the solution u^(n-1)
C      = -3 : Standard macro time step consisting of 3 substeps, 
C             i.e. time-step with step-size TSTEP. 
C             Use TIM3x for adaptive time-stepping.
C             No information about previous time-step available,
C             i.e. LTML does not define u^(n-1).
C      =  3 : Anywhere during the simulation, small/standard time-step
C             with step-size TSTEP. 
C             Use TIM3x for adaptive time-stepping.
C             LTML identifies the solution u^(n-1).
C
C     When linear extrapolation should be used, ITEXL must be = 1
C     on entry of NONSTL. It will change during NONSTL in that case,
C     but will have the same value as it had before when NONSTL ends.

      ITEXL=0
      
      IF ((INLMIN.EQ.INLMAX).AND.(INLMIN.EQ. 1)) THEN
      
        ITEXL=1
        
C       Initialize adaptive time-step control weights for adaptive
C       time stepping with linear extrapolation (cf. p. 204/205, 
C       (187/188) Turek's book).
C       TIML3x describes the two weights of the velocity vectors
C       that are used for each of the three small substeps with
C       stepsize TSTEP of a macrostep.
C       TIML1x on the other side describes the weights that
C       are used in front of the velocity vectors for the
C       single predictor step with step size 3xTSTEP
C       that is done previously to the substeps.

        TIML31=1D0
        TIML32=1D0
        TIML11=1D0
        TIML12=1D0
        
C       The linear extrapolation approximates the velocity u^(l+1) with
C       u~ computed by (cf. p. 204f (187f), Turek's book)
C
C       u~ = [ dT2 / (dT2-dT1) ] * u(t-dT1)  +  [ -dT1 / (dT2-dT1) ] * u(t-dT2)
C       
C       The TIMLxy-parameter are used as:
C
C         TIMLx1 = dT1       
C         TIMLx2 = (dT2-dT1) 
C
C       In the begining of the iteration, we have step size TSTEP (at least 
C       with 1-step scheme). Nevertheless, we initialize TIMLx1=TIMLx1=1 !
C       Why? This is just a dummy! The time stepping routines will later
C       directly before the calculation initialize TIMLx1=TSTEP. On the
C       other hand, the solution corresponding to TIMLx2=1 will be 0.
C       Therefore initializing TIMLx1=TIMLx2=1 will not harm here!
C       Indeed, it's the correct initialization, because then the fractions
C       for u~ (see above) will calculate correctly in the first two time
C       steps!
        
      ENDIF

      IF ((INLMIN.EQ.INLMAX).AND.(INLMIN.EQ.-1)) THEN
        INLMIN=1
        INLMAX=1
      ENDIF

      IMTIME=2
      TTMG =0D0
      TTS  =0D0
      TTE  =0D0
      TTD  =0D0
      TTP  =0D0
      TTR  =0D0

      TTGRID=0D0
      TTPOST=0D0
      TTADF =0D0
      TTUPW =0D0
      TTBDR =0D0
      TTLC  =0D0

C     If we perform a stationary simutaion, switch the adaptive
C     time-stepping (Theta-Scheme) off.

      IF (ISTAT.EQ.0) THEN
        IADTIM=0
        TSTEP =1D0
        THETA =1D0
        IFRSTP=0
      ENDIF

C     In instationary simulation we either use a standard 1-step scheme
C     or Fractional-Step. Depending on what is configured by the
C     parameter IFRSTP, we have to initialize the Theta scheme.
C
C     The time-stepping algorithm always collects 3 steps with stepsize
C     TSTEP (configured in the DAT-file) to one macro-step of stepsize
C     3xTSTEP (note that every Theta-scheme can be decomposed into
C     three substeps) - so our logical stepsize we are working everywhere
C     through the algorithm is 3xTSTEP and only TSTEP if we are 
C     calculating the substeps.

      IF ((ISTAT.NE.0).AND.(IFRSTP.EQ.1)) THEN
      
C       The FS Theta-Scheme uses by theory 4 parameters: 
C
C         Theta   = 1 - sqrt(2) / 2
C         Theta'  = 1 - 2 * Theta
C         alpha   = ( 1 - 2 * Theta ) / ( 1 - Theta )
C         beta    = 1 - alpha
C
C       The parameter THETA in the DAT-file is ignored and replaced
C       by a hard-coded setting.
      
        THETA =1D0-SQRT(0.5D0)
        THETAP=1D0-2D0*THETA
        FALPHA=THETAP/(1D0-THETA)
        FBETA =THETA /(1D0-THETA)

C       Each macro-step u^n -> u^(n+1) is then split up in three
C       substeps: 
C         u^n         -> u^(n+1/3),
C         u^(n+1/3)   -> u^(n+2/3),
C         u^(n+2/3)   -> u^(n+1)    = result for the next time-step
C
C       A general form of the Theta-scheme is given by:
C
C         [ I + Theta~ N(u^new) ] u^new   +   ... grad(p^new)   =   RHS
C
C       in 1-step scheme, Theta~ measures the step-length of the
C       substeps (see below). In FS-Theta-Scheme, Theta~ is defined by
C
C         Theta~   =   alpha * Theta * K   =   beta * Theta' * K
C
C       with K = step-length of the macro-step. Remember, our
C       macro-stepsize is 3xTSTEP.

        THSTEP=3D0*TSTEP*FALPHA*THETA
        
      ELSE
      
C       The standard 1-step scheme uses Theta=0.5 (Crank Nicolson)
C       or Theta=1 (Backward Euler). This is configured in the 
C       DAT-file. 
C       The other case (stationary simulation) also sets Theta
C       this way, but Theta is overwritten in the simulation
C       subroutine again.
      
        THSTEP=TSTEP*THETA
        
      ENDIF
      
C     The THSTEP-parameter is set here only for sure. It is set in the
C     MG-step procedure again according to the current time stepping
C     mode.
      
      END
      
      
************************************************************************
* Initialise Point value files
*
* This routine opens output channels 40..53 for output of points
* that are calculated during the solution process.
************************************************************************
      
      SUBROUTINE INPWFS
      
      IMPLICIT NONE

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

      END
