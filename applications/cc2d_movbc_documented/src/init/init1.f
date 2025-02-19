************************************************************************
      SUBROUTINE INIT1 (MSHOW,IWORKG,IWMAXG)
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
*            - etc
*
* The routine is mainly a collection of calls to ALCxxx and GENxxx-
* (and perhaps DONxxx-) routines. These subroutines can be found
* in ALCGEN.F and DONE.F.
*
* Before this routine is called, RDOUT and RDDAT must have been
* executed!
*
************************************************************************
      IMPLICIT NONE
      
C-----------------------------------------------------------------------
C     C O M M O N S 
C-----------------------------------------------------------------------

C *** Standard COMMON blocks
      INCLUDE 'cmem.inc'
      INCLUDE 'cmem2.inc'

C *** Output/error variables
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'

C *** multigrid variables      
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

C *** element definitions
      INCLUDE 'cbasicelem.inc'      
      
C *** user COMMON blocks
      INCLUDE 'cinidat.inc'
      
C *** Names of matrices and vectors (for output and tracing)
      INCLUDE 'fnstokes.inc'
      INCLUDE 'fnmassmat.inc'
      
C *** Common blocks of triangulation

      INCLUDE 'ctria.inc'

C define element E031 as external; this is once used for building
C the matrix structure for all elements of the Ex3x-set

      EXTERNAL E031

C Fictitious boundary management
      
      INTEGER NFBDYC
      EXTERNAL NFBDYC

C=======================================================================
C     local variables and constants
C=======================================================================

C parameters

      INTEGER MSHOW,MFILE
      INTEGER IWORKG,IWMAXG

C further local variables, parameters

      CHARACTER CFILE*60

      DOUBLE PRECISION TTT0,TTT1,TTT0L,TTLIN
      
      INTEGER ISYMMA
      INTEGER ISETLV,NEQ,NEQB
      
      INTEGER LF12P, LF1, LF2, IFMTS
      INTEGER I1, II, KU1C

C=======================================================================
C     Initialization
C=======================================================================
      SUB='INIT1 '

      CALL ZTIME(TTT0)

C=======================================================================
C     Data input
C=======================================================================

      CFILE=CFILE1
      MFILE=MFILE1

C=======================================================================
C     Grid generation
C=======================================================================

C parametrisation

      CALL GENPAR (.TRUE.,IMESH1,CPARM1)
      
C triangulation. II specifies the number of pre-refinements of
C the coarse grid if the mesh is not read in from a file. If NLMAX is
C set > 9, we shift NLMIN and NLMAX by some levels so that NLMAX is 9.

      II = 0
      IF (NLMAX.GT.9) THEN
        II = NLMAX-9
        NLMIN = MAX(1,NLMIN-II)
        NLMAX = 9
        WRITE (MTERM,'(A)') 'Warning: NLMAX too large.'
        WRITE (MTERM,'(A)') 'Performing level shift for calculation '//
     *                      'on finer levels.'
        WRITE (MTERM,'(A,I3,A,I3)') 'Calculation will be on level ',
     *                      NLMIN+II,'..',NLMAX+II
      END IF
      
      CALL GENORS (.TRUE.,IRMESH,II,CMESH1)

C print mesh statistics

      DO II=NLMIN,NLMAX
        IF (MSHOW.GE.2) 
     *    WRITE(MTERM,*) 'ILEV,NVT,NMT,NEL,NVBD: ', II,
     *          KNVT(II),KNMT(II),KNEL(II),KNVBD(II)
        IF (MSHOW.GE.1) 
     *    WRITE(MFILE,*) 'ILEV,NVT,NMT,NEL,NVBD: ', II,
     *          KNVT(II),KNMT(II),KNEL(II),KNVBD(II)
      END DO

C=======================================================================
C     Solver preparation
C=======================================================================
   
C Allocate memory for boundary information backup

      DO II=NLMIN,NLMAX
        IF ((IBDR.GE.2).OR.(NFBDYC().GT.0)) THEN
C make a backup of KNPR in case of time dependent boundary conditions
C or moving boundaries
          CALL ALCKNP (II)
          IF (IER.NE.0) GOTO 99999
        ENDIF

C allocate memory for boundary-midpoint-information

        CALL ALCMBD(II)
        IF (IER.NE.0) GOTO 99999
      END DO

C implement boundary information; update NMBD

      DO II=NLMIN,NLMAX
        CALL IMPBDG (II)
      END DO

      CALL ZTIME(TTT1)
      TTGRID=TTT1-TTT0

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
      
C allocate pointer structures for matrix using one of
C E030/E031/EM30/EM31-element
C --> all have the same structure

C don't exploit any symmetry

      ISYMMA=0
      CALL XMAP7(KLCOLA,KLLDA,KNA,KNEQA,E031,ISYMMA)
      IF (IER.NE.0) GOTO 99999
C
C=======================================================================
C *** Generation of block ST
C --> Laplace matric int(grad(phi_j),grad(phi_i))
C     as alternative to building the full nonlinear A-matrix
C     in every iteration step.
C     --> by adding a nonlinearity to a in every step.
C
C This is only used if IPRECA!=4, because in the case =4 the
C real matrix A is build in every iteration step.
C=======================================================================
C
      IF (IPRECA.NE.4) THEN
       CALL GENSTM 
      ENDIF
C
C=======================================================================
C *** Generation of mass matrices in instationary case
C=======================================================================
C
      IF (ISTAT.NE.0) THEN
        CALL GENMM (IMASS.EQ.0,IMASSL.NE.0)
      END IF
C
C=======================================================================
C *** Generation of blocks B1,B2
C=======================================================================

C generate B-blocks by usual quadrature or by exact evaluation

      CALL GENMTB(IPRECB.GE.2)

C show statistics and quick-check generated information

      DO ILEV=NLMIN,NLMAX

        IF (MSHOW.GE.2) WRITE(MTERM,*) 'ILEV,NU,NA,NB:',ILEV,KNU(ILEV),
     *                                  KNA(ILEV),KNB(ILEV)
        IF (MSHOW.GE.0) WRITE(MFILE,*) 'ILEV,NU,NA,NB:',ILEV,KNU(ILEV),
     *                                  KNA(ILEV),KNB(ILEV)

        NEQ =KNEQA(ILEV)
        NEQB=KNEQB(ILEV)
        IF (NEQB.NE.NEQ.OR.NEQ.NE.KNU(ILEV).OR.KNEL(ILEV).NE.KNP(ILEV)) 
     *    THEN
          WRITE(MTERM,*) 'ERROR in INIT1: NEQ.NE.NEQB'
          STOP
        ENDIF
        
      END DO

      IF (MSHOW.GE.2) WRITE(MTERM,*)
      IF (MSHOW.GE.1) WRITE(MFILE,*)
 
      CALL ZTIME(TTT1)
      TTADF=TTADF+TTT1-TTT0
C
C=======================================================================
C     matrix restructuring
C=======================================================================
C
      DO ILEV=NLMIN,NLMAX
        CALL ZTIME(TTT0)
       
C directly change the data type of Stokes-/mass matrices as given in
C the parameters in the .DAT file:
C * code removed because all matrices are now double! *
C        CALL CHMATP (ILEV,2,(IPRECA.EQ.0).OR.(IPRECA.GE.2),
C     *                      (IPRECB.EQ.0).OR.(IPRECB.EQ.3),
C     *                      ((IPRECA.EQ.0).OR.(IPRECA.GE.2)).AND.
C     *                       (ISTAT.EQ.1))
       
       IF ((IPRECA.EQ.2).OR.(IPRECA.EQ.3)) THEN

C Write the matrix to the hard disc

        CALL  OF0 (59,CFILST(ILEV),0)
        CFILE='STMAT '
C double precision
        CALL OWA1 (DWORK(L(KLST(ILEV))),CFILE,KNA(ILEV),59,0)
        REWIND(59)
        CLOSE (59)
        IF (IER.NE.0) GOTO 99999
        
C Release the Stokes matrix; it's read in later

        CALL DISSTM(ILEV)
        IF (IER.NE.0) GOTO 99999

        IF (ISTAT.EQ.1) THEN
        
C Write the mass matrix to hard disc

         CALL  OF0 (59,CFILM(ILEV),0)
         CFILE='MASMAT'
C double precision
         CALL  OWA1 (DWORK(L(KLM(ILEV))),CFILE,KNA(ILEV),59,0)
         REWIND(59)
         CLOSE (59)
         IF (IER.NE.0) GOTO 99999
         
C Release the mass matrix; it's read in later

         CALL DISMM (ILEV)
         IF (IER.NE.0) GOTO 99999
        ENDIF
       ENDIF

       IF (IPRECB.EQ.2) THEN
C release all handles
        CALL DISMTB(ILEV)
       ENDIF

       CALL ZTIME(TTT1)
       TTADF=TTADF+TTT1-TTT0

      END DO

      IF (MSHOW.GE.2) WRITE(MTERM,*)
      IF (MSHOW.GE.0) WRITE(MFILE,*)
C
C=======================================================================
C    Allocation of:  - solution vector with boundary conditions on NLMAX
C                    - RHS on NLMAX and auxiliary vectors
C                    - UOLD-vector on level NLMAX
C=======================================================================
C
      DO ILEV=NLMIN,NLMAX

        CALL ZTIME(TTT0)
        ISETLV=1
        CALL  SETLEV (ISETLV)

C=======================================================================
C *** Allocation of solution and defect vectors and right hand side 
C=======================================================================

C Full right hand side vector is allocated on every level
C (for U,V and P in one vector)
C Solution vectors for velocities are only allocated on maximum level.
C Solution vectors for pressures are allocated on every level.
C Backups of solution vectors are only allocated on finest level
C if necessary for the calculation.

        CALL ALCSOL (ILEV,.TRUE.,.TRUE.,(ILEV.EQ.NLMAX),
     *               (ILEV.EQ.NLMAX).AND.
     *               ((OMGMIN.GT.0D0).OR.(OMGMAX.GT.0D0)))

C=======================================================================
C *** Allocation of iteration matrix A 
C This is build with the Stokes matris ST+correction term in every step
C=======================================================================

        CALL ALCMTA (ILEV)
        IF (IER.NE.0) GOTO 99999

C=======================================================================
C *** calculation of a vector with the areas of all finite elements
C=======================================================================

        CALL ALCARE (ILEV)
        CALL SETARE(DWORK(L(KLAREA(ILEV))),KNEL(ILEV),
     *       KWORK(L(KLVERT(ILEV))),DWORK(L(KLCVG(ILEV)))) 

C=======================================================================

        CALL ZTIME(TTT1)
        TTLC=TTLC+TTT1-TTT0
        TTT0 = TTT1

C=======================================================================
C *** calculation of rhs and start vectors on finest level only
C=======================================================================

C RHS/start vector handling has only to be performed on 
C the finest level:

        IF (ILEV.EQ.NLMAX)  THEN

C At first generate the right hand side on the finest level
C for U and V-velocity; the return values LF1 and LF2 are new
C handles to these vectors.

          CALL GENRHS(ILEV,LF1,LF2)

C Copy these vectors to our global right hand side on the
C finest level. Afterwards the handles can be released,
C we continue working with our global representation.

          LF12P = KLF12P(ILEV)
          CALL  LCP1 (DWORK(L(LF1)), DWORK(L(LF12P)),KNEQA(ILEV))
          CALL  LCP1 (DWORK(L(LF2)), DWORK(L(LF12P)+KNEQA(ILEV)), 
     *                KNEQA(ILEV))
          CALL  ZDISP (0, LF1, 'LF1TMP' )
          CALL  ZDISP (0, LF2, 'LF2TMP' )
          IF (IER.NE.0) GOTO 99999

C=======================================================================
C *** start vector as prolongation from level NLMAX-(ISTART-1) or
C     ISTART-100, read from disc
C=======================================================================
        
          IF (ABS(ISTART).GT.0) THEN
            
C calculate the level of the start vector
            
            IF (ISTART.LE.100) THEN
              I1=NLMAX-(ISTART-1)
            ELSE
              I1=ISTART-100
            END IF
          
            IF (I1.GT.NLMAX) THEN
            
C Oops, we can't go backwards in levels; restriction only works for a 
C defect vector, whereas we can use prolongation to go forward in 
C levels with solution vectors.

C Theoretically a restriction of the start solution can be performed
C by RESTRU. Unfortunately we don't know how long this vector is
C and so it's not implemented up until now...

              WRITE(MTERM,*) 'Warning: level of predefined start'//
     *                       ' vector too high!'
              WRITE(MTERM,*) 'Using zero vector as start vector.'
              WRITE(MFILE,*) 'Warning: level of predefined start'//
     *                       ' vector too high!'
              WRITE(MFILE,*) 'Using zero vector as start vector.'
            ELSE
            
              KU1C=L(KLUP(I1))
              IF (ISTART.GT.0) THEN
                IFMTS=0
              ELSE
                IFMTS=1
              END IF
              
C read the vector from the file

              CFILE='DU12P '
              CALL  OF0 (MSTART,CSTART,IFMTS)
              CALL  ORA1 (DWORK(KU1C),CFILE,MSTART,IFMTS)
              CLOSE(MSTART)
              
C prolongate it to the current level

              CALL MLTPRL (I1,NLMAX)
            END IF
          END IF

          CALL ZTIME(TTT1)
          TTLC=TTLC+TTT1-TTT0
C
C=======================================================================
C *** Dirichlet boundary updates of solution and rhs vector
C=======================================================================

          CALL ZTIME(TTT0)
          
          CALL IMPSLR (ILEV,(IBDR.EQ.1).AND.(ISTAT.EQ.0),.TRUE.,.TRUE.)

          CALL ZTIME(TTT1)
          TTBDR=TTBDR+TTT1-TTT0

        ENDIF
C
C=======================================================================
C *** Auxiliary vector on all levels
C=======================================================================
C
        CALL ZTIME(TTT0)
        
        CALL ALCAUX(ILEV)
        IF (IER.NE.0) GOTO 99999
        
        CALL ZTIME(TTT1)
        TTLC=TTLC+TTT1-TTT0

      END DO

      CALL ZTIME(TTT1)
      TTLIN=TTT1-TTT0L

      IF (MSHOW.GE.2) 
     * WRITE(MTERM,*) 'time for initialization of linear operators : ', 
     *                TTLIN
      IF (MSHOW.GE.0) 
     * WRITE(MFILE,*) 'time for initialization of linear operators : ', 
     *                TTLIN
      IF (MSHOW.GE.2) WRITE(MTERM,*)
      IF (MSHOW.GE.0) WRITE(MFILE,*)
      
99999 END
