************************************************************************
      PROGRAM  CC2D_MOVBC  
************************************************************************
*
*   Explanation: works similar as standard CC2D 
*
*   Differences: bndry.f contains position (and movement) of "fictitious
*                boundary components" 
*
*                indat2d.f contains position (and movement) of 
*                "fictitious boundary components" and additionally the 
*                Dirichlet values on the surface
*
*   Perform the given example and modify w.r.t. your application !!!
*
*
*   UNOFFICIAL TEST VERSION ONLY !!!!
*
*   The file FICTBDRY.F defines the shape, position,... of the
*   fictitious boundary object(s). This file can be configured
*   by the user.
*
************************************************************************

      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cfiles.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cns.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      
      INCLUDE 'cinidat.inc'
      INCLUDE 'cinigeometry.inc'

      CHARACTER CDATA*60

C externals

C *** Exact solution (for  error analysis)
      DOUBLE PRECISION UE
      EXTERNAL UE
      
C local variables

      INTEGER MDATA,MFILE,MSHOW,IWORKG,IWMAXG,IWORKI,IWMAXI
      INTEGER IFILEN,ITFILM
      INTEGER IFPOST

      DOUBLE PRECISION TTTSUM,TTT0,TTT1,TTTH,TTLIN
      DOUBLE PRECISION DFWVI,DAWVI

C=======================================================================
C     Initialization
C=======================================================================

C At the very beginning initialise the pseudodynamic memory management.
C Nothing must be done before that!!!

      CALL ZINIT(NNWORK,'feat.msg','#data/CC2D.ERR','#data/CC2D.PRT',
     *                             '#data/CC2D.SYS','#data/CC2D.TRC') 

C Now read parameters and initialise the output structures.

C Use unit number 62 for file output and open the file as prescribed in
C the DAT file.

      MDATA=79
      CDATA='#data/cc2d.dat'
      
      CALL RDOUT (MDATA,CDATA,62,MSHOW)
      MFILE = MFILE1
      
C Read parameters of the DAT file

      CALL RDDAT (MDATA,CDATA,MSHOW)
      
C Set up geometry data

      CALL RDGEOM ()

      WRITE (MTERM,'(A)') '----------------------------------------'//
     *                    '---------------------------------------'

1091  CALL ZTIME(TTTSUM)

      CALL ZTIME(TTT0)

C Initialisation of the solver structures
      
      CALL  INIT1 (MFILE,IWORKG,IWMAXG)
      IWORKI=IWORK
      IWMAXI=IWMAX

      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0

C=======================================================================
C *** Call the (non)steady loop to solve
C=======================================================================

C Open output channels

      CALL INPWFS ()
      
C Initialise NS time stepping

      CALL INNSTT ()
      
C initialise output variables for nonsteady loop:

      ITFILM=IFINIT-1
      IFILEN=0

C Initialise counter variables in /NSCOUN/; these are incremented 
C during the loop

      NSUBST = 0
      NNONL  = 0
      NMGU   = 0

C Start the nonsteady loop

      CALL NONSTL (MFILE,MSHOW,ITFILM,IFILEN)
      IF (IER.NE.0) GOTO 99998

C=======================================================================
C *** Postprocessing
C=======================================================================

      CALL ZTIME(TTT0)
C
      IFPOST=0
      IFILEN=0
      CALL FPOST(IFPOST,IFILEN,ITFILM,UE,MSHOW,DFWVI,DAWVI,0,'')
      
C *MK*: DFWVI/DAWVI added to return drag/lift
C
      CALL ZTIME(TTT1)
      TTPOST=TTPOST+TTT1-TTT0
C
C=======================================================================
C     Statistics
C=======================================================================
C
      IF (MSHOW.GE.2) WRITE (MTERM,*)
      IF (MSHOW.GE.0) WRITE (MFILE,*)
      IF (MSHOW.GE.2) WRITE (MTERM,*) 'STATISTICS :'
      IF (MSHOW.GE.0) WRITE (MFILE,*) 'STATISTICS :'
C
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'NWORK :      ',NWORK
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'NWORK :      ',NWORK
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'IWORKG:      ',IWORKG
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'IWORKG:      ',IWORKG
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'IWMAXG:      ',IWMAXG
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'IWMAXG:      ',IWMAXG
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'IWORKI:      ',IWORKI
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'IWORKI:      ',IWORKI
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'IWMAXI:      ',IWMAXI
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'IWMAXI:      ',IWMAXI
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'IWORK :      ',IWORK
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'IWORK :      ',IWORK
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'IWMAX :      ',IWMAX
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'IWMAX :      ',IWMAX
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.0) WRITE(MFILE,1)
C
      CALL ZTIME(TTT1)
      TTTSUM=TTT1-TTTSUM
      TTTH  =TTGRID+TTPOST+TTADF+TTUPW+TTBDR+TTLC+TTMG
      TTLIN =TTADF+TTUPW+TTBDR+TTLC
C
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'total time : ', TTTSUM
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'total time : ', TTTSUM
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'appr. time : ', TTTH
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'appr. time : ', TTTH
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'grid  time : ', TTGRID
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'grid  time : ', TTGRID
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'post  time : ', TTPOST
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'post  time : ', TTPOST
C
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'lin.  time : ', TTLIN
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'lin.  time : ', TTLIN
      IF (MSHOW.GE.2) WRITE(MTERM,*) '-> mavec time : ', TTADF
      IF (MSHOW.GE.0) WRITE(MFILE,*) '-> mavec time : ', TTADF
      IF (MSHOW.GE.2) WRITE(MTERM,*) '-> konv. time : ', TTUPW
      IF (MSHOW.GE.0) WRITE(MFILE,*) '-> konv. time : ', TTUPW
      IF (MSHOW.GE.2) WRITE(MTERM,*) '-> bdry  time : ', TTBDR
      IF (MSHOW.GE.0) WRITE(MFILE,*) '-> bdry  time : ', TTBDR
      IF (MSHOW.GE.2) WRITE(MTERM,*) '-> LC    time : ', TTLC
      IF (MSHOW.GE.0) WRITE(MFILE,*) '-> LC    time : ', TTLC
C
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'mg    time : ', TTMG
      IF (MSHOW.GE.1) WRITE(MFILE,*) 'mg    time : ', TTMG
C
      IF (MSHOW.GE.2) WRITE(MTERM,*) '#substeps  : ', NSUBST
      IF (MSHOW.GE.0) WRITE(MFILE,*) '#substeps  : ', NSUBST
      IF (MSHOW.GE.2) WRITE(MTERM,*) '#nonlinear : ', NNONL
      IF (MSHOW.GE.1) WRITE(MFILE,*) '#nonlinear : ', NNONL
      IF (MSHOW.GE.2) WRITE(MTERM,*) '#mg        : ', NMGU
      IF (MSHOW.GE.1) WRITE(MFILE,*) '#mg        : ', NMGU
      IF (MSHOW.GE.2) WRITE(MTERM,1)
      IF (MSHOW.GE.1) WRITE(MFILE,1)
C
      IF (MSHOW.GE.2) THEN
      WRITE(MTERM,*)
      WRITE(MTERM,*) ' MULTIGRID COMPONENTS [in percent]:'
      WRITE(MTERM,*) ' smoothing     :', 1.D2*TTS/TTMG
      WRITE(MTERM,*) ' solver        :', 1.D2*TTE/TTMG
      WRITE(MTERM,*) ' defect calc.  :', 1.D2*TTD/TTMG
      WRITE(MTERM,*) ' prolongation  :', 1.D2*TTP/TTMG
      WRITE(MTERM,*) ' restriction   :', 1.D2*TTR/TTMG
      WRITE(MTERM,1)
      ENDIF
C
      IF (MSHOW.GE.0) THEN
      WRITE(MFILE,*)
      WRITE(MFILE,*) ' MULTIGRID COMPONENTS [in percent]:'
      WRITE(MFILE,*) ' smoothing     :', 1.D2*TTS/TTMG
      WRITE(MFILE,*) ' solver        :', 1.D2*TTE/TTMG
      WRITE(MFILE,*) ' defect calc.  :', 1.D2*TTD/TTMG
      WRITE(MFILE,*) ' prolongation  :', 1.D2*TTP/TTMG
      WRITE(MFILE,*) ' restriction   :', 1.D2*TTR/TTMG
      WRITE(MFILE,1)
      ENDIF
C
      GOTO 99999
C=======================================================================
C     Error case
C=======================================================================
99998 WRITE(MTERM,*) 'IER', IER
      WRITE(MTERM,*) 'IN SUBROUTINE ',SUB
C
99999 CLOSE(MFILE)

C Clean up geometries

      CALL DONGEO ()

C Clean up all handles that are allocated in INIT1.
C If this produces an error, check the file done.f that handles
C the deallocation of all arrays!

      CALL DONE1 (0)

   1  FORMAT(80('-'))

99997 END
