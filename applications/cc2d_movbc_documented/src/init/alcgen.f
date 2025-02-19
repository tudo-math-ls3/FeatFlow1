************************************************************************
* The following subroutines do the real work in the initialization.
* They are typically called by INIT1, but can also be called
* independently if necessary.
*
* All routines work by using the parameters that are provided in the
* COMMON blocks, so before calling them, the necessary information
* (e.g. NVT, NMT,...) must be provided there! This is normally
* done in the INIT1-subroutine.
*
* General rule for ALCGEN.F and DONE.F:
* - These routines handle the information in the COMMON blocks,
*   especially for multigrid.
* - ALCxxx-routines only allocate the memory, perhaps depending on
*   their parameters
* - GENxxx-routines build information structures in the COMMON blocks.
*   Depending on the type of information, memory is either allocated
*   automatically or has to be allocated previously by a call to
*   an ALCxxx-routine.
*   Some GENxxx-routines (or the subroutines they call) need to work
*   with grid structures in /TRIAx/. If that's the case, these
*   routines automatically switch the current level if necessary.
* - DONxxx-routines release the memory, ALCxxx/GENxxx-routines
*   allocate.
* - IMPxxx-routines modify existing structures by implementing 
*   necessary information into these.
*
************************************************************************

C Switch level, but only if necessary. -> Local helper routine

      SUBROUTINE SWTLEV (ILV)
      
      IMPLICIT NONE
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      
      INTEGER ILV,ISETLV
      
      IF (ILEV.NE.ILV) THEN
        ILEV=ILV
        ISETLV=1
        CALL SETLEV (ISETLV)
      END IF
      
      END      

C ----------------------------------------------------------------------
C The following routine allocates the vector KLNPRO(ILV)
C and stores a copy of the vector KLNPR(ILV) there. 
C This can be used in case of time dependent boundary conditions
C to make a backup of KLNPR(ILV) so it can be restored every time.
C
C The routine does not check the parameter IBDRY of the .DAT file,
C so the caller must take care about if this routine should be called
C or not.
C ----------------------------------------------------------------------

      SUBROUTINE ALCKNP (ILV)
      
      IMPLICIT NONE

      INCLUDE 'cerr.inc'

      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgtria.inc'
      
      INCLUDE 'cnsparfrac.inc'
      
      INTEGER ILV
      INTEGER NVT, NMT

C Allocate the memory and make a backup of KNPR.
C Make copies of the MG-information to prevent a possible destruction by ZNEW!

      NVT = KNVT(ILV)
      NMT = KNMT(ILV)
      
      CALL ZNEW(NVT+NMT,-3,KLNPRO(ILV),'KNPRO ')
      IF (IER.EQ.0) THEN
        CALL XLCP3(KLNPR(ILV),KLNPRO(ILV),NVT+NMT)
      END IF

      END 


C ----------------------------------------------------------------------
C Allocate vector with area of elements on level ILV.
C ----------------------------------------------------------------------

      SUBROUTINE ALCARE (ILV)
      
      IMPLICIT NONE

      INCLUDE 'cerr.inc'

      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgtria.inc'
      
      INCLUDE 'cnsparfrac.inc'

      INTEGER ILV
      INTEGER LAREA,NEL

C Allocate the memory and make a backup of KNPR.
C Make copies of the MG-information to prevent a possible destruction by ZNEW!

      NEL = KNEL(ILV)

      CALL ZNEW(NEL+1,1,LAREA,'DAREA ')
      IF (IER.EQ.0) THEN
        KLAREA(ILV)=LAREA
      END IF

      END 

C ----------------------------------------------------------------------
C Allocate matrix A on level ILV.
C ----------------------------------------------------------------------

      SUBROUTINE ALCMTA (ILV)
      
      IMPLICIT NONE

      INCLUDE 'cerr.inc'

      INCLUDE 'cmem.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      
      INCLUDE 'cnsparfrac.inc'

      INTEGER ILV
      INTEGER NA,LA1

C Make copies of the MG-information to prevent a possible destruction by ZNEW!

      NA=KNA(ILV)
      CALL  ZNEW (NA,1,LA1,'DA    ')
      IF (IER.EQ.0) THEN
        KLA(ILV)=LA1
      END IF

      END 


C ----------------------------------------------------------------------
C Allocate auxiliary vectors on level ILV.
C ----------------------------------------------------------------------

      SUBROUTINE ALCAUX (ILV)
      
      IMPLICIT NONE

      INCLUDE 'cerr.inc'

      INCLUDE 'cmem.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      
      INCLUDE 'cnsparfrac.inc'

      INTEGER ILV
      INTEGER LAUX,NUP

C Allocate the memory and make a backup of KNPR.
C Make copies of the MG-information to prevent a possible destruction by ZNEW!

      NUP = KNUP(ILV)
      
      CALL ZNEW(NUP,1,LAUX,'DAUX  ')
      IF (IER.EQ.0) THEN
        KLAUX (ILV)=LAUX
      END IF

      END 


C ----------------------------------------------------------------------
C Allocate structures describing the midpoints on the
C boundary on level ILV:
C ----------------------------------------------------------------------

      SUBROUTINE ALCMBD (ILV)
      
      IMPLICIT NONE

      INCLUDE 'cerr.inc'

      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgtria.inc'
      
      INTEGER ILV
      INTEGER NMT

C Allocate the memory and make a backup of KNPR.
C Make copies of the MG-information to prevent a possible destruction by ZNEW!

      NMT=KNMT(ILV)
      
      CALL ZNEW(NMT,3,KLMBD(ILV),'KMBD  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NMT,1,KLDBD(ILV),'DDBD  ')
      IF (IER.NE.0) GOTO 99999

99999 END 

C ----------------------------------------------------------------------
C Allocate memory for solution vectors and right hand sides.
C
C In:
C  ILV    - the level
C  BRHS   - allocate memory for right hand sides level ILV
C  BSOL   - allocate memory for solution vectors level ILV
C  BSLD   - allocate memory for defect correction vectors
C           (should be done only once on finest level)
C  BSLOLD - allocate memory for backup of solution vectors
C           (should be done only once on finest level)
C
C The handles of the solution vector is normally only relevant on 
C the finest level, so the caller should switch to the finest level 
C before calling this subroutine! Handles of solution vectors are
C stored in the LD1, LD2,...-variables of the COMMON block /MGFLD/.
C ----------------------------------------------------------------------

      SUBROUTINE ALCSOL (ILV, BRHS, BSOL, BSLD, BSLOLD)
      
      IMPLICIT NONE

      INCLUDE 'cerr.inc'

      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgtria.inc'
      
      INTEGER ILV
      LOGICAL BRHS, BSLD, BSOL, BSLOLD
      INTEGER NEQ, NUP, NEL, LF12P, LUP

C Allocate the memory and make a backup of KNPR.
C Make copies of the MG-information to prevent a possible destruction by ZNEW!

      NEQ=KNEQA(ILV)
      NEL=KNEL(ILV)
      NUP=2*NEQ+NEL

C backup of defect correction vectors on finest level

      IF (BSLOLD) THEN
        CALL ZNEW(2*NEQ+NEL,1,LU1OLD,'DU1OLD')
        IF (IER.NE.0) GOTO 99999
      END IF

C defect correction vectors on finest level

      IF (BSLD) THEN
        CALL ZNEW(2*NEQ+NEL,1,LD1,'DD12P ')
        IF (IER.NE.0) GOTO 99999
      END IF

C solution vectors... on each level for use in nonlinear iteration

      IF (BSOL) THEN
        CALL ZNEW(NUP,1,LUP,'DU12P ')
        IF (IER.NE.0) GOTO 99999
        KLUP(ILV)=LUP
      END IF 
      
C right hand side vectors
      
      IF (BRHS) THEN
        CALL ZNEW(NUP,1,LF12P,'DF12P ')
        IF (IER.NE.0) GOTO 99999
        KLF12P(ILV)=LF12P
      END IF

99999 END 

C ----------------------------------------------------------------------
C Generate the Stokes-/Laplace matrices ST on all levels.
C
C This matrix can be used as a replacement for building the whole
C system matrix A in every iteration step. For this purpose
C in every step a nonlinear term is added to ST to produce an
C approximation for the A-block.
C
C The caller has to decide whether the matrix ST is used or not;
C the parameters in the .DAT file (whether ST should be built or not)
C are not considered. This routine only cares about the parameters
C in the .DAT file that configure the matrix, e.g. used element,...
C
C The resulting matrices will have double precision.
C ----------------------------------------------------------------------

      SUBROUTINE GENSTM 
      
      IMPLICIT NONE

      INCLUDE 'cerr.inc'

      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgtria.inc'
      
      INCLUDE 'cbasicelem.inc'
      
      INCLUDE 'cinidat.inc'

C *** Names of matrices and vectors (for output and tracing)

      INCLUDE 'arnstokes.inc'

C externals

C *** Coefficient of stiffness matrix, right hand side, exact solution
      DOUBLE PRECISION COEFST,COEFFB,RHS,UE
      EXTERNAL COEFST,COEFFB,RHS,UE
C *** definition of finite elements
      EXTERNAL E030,E031,EM30,EM31,E010
      
      INTEGER ICLRA, ISYMM

C configuration of the matrix blocks

      INTEGER NBLOCA
      PARAMETER (NBLOCA=1)
      INTEGER NBLA1      
      PARAMETER (NBLA1=NBLOCA*NNLEV)

C local variables: structure of bilinear form of Laplacian matrix:
      
      INTEGER KABSTN(NBLOCA),KABST(2,NNAB,NBLOCA)
      
      LOGICAL BCONST(NBLOCA)      
      DATA BCONST/.TRUE./
      SAVE BCONST
      
C structure constants for matrices      
      
      LOGICAL BSNGLA(NBLOCA,NNLEV)
      DATA BSNGLA/NBLA1*.FALSE./      
      SAVE BSNGLA 
      
C initialise generation-descriptors

       KABST (1,1,1)=2
       KABST (2,1,1)=2
       KABST (1,2,1)=3
       KABST (2,2,1)=3
       KABSTN(1)    =2

C Don't exploit any symmetry
       
       ISYMM = 0

       ICLRA=1
       IF (IELT.EQ.0) 
     *  CALL XMAB07(KLST,KLCOLA,KLLDA,KNA,KNEQA,NBLOCA,ICLRA,E031,
     *              COEFST,BCONST,KABST,KABSTN,ICUBA,ISYMM,CARRST,
     *              BSNGLA)
       IF (IELT.EQ.1) 
     *  CALL XMAB07(KLST,KLCOLA,KLLDA,KNA,KNEQA,NBLOCA,ICLRA,E030,
     *              COEFST,BCONST,KABST,KABSTN,ICUBA,ISYMM,CARRST,
     *              BSNGLA)
       IF (IELT.EQ.2)
     *  CALL XMABM7(KLST,KLCOLA,KLLDA,KNA,KNEQA,NBLOCA,ICLRA,EM31,
     *              COEFST,BCONST,KABST,KABSTN,ICUBA,ISYMM,CARRST,
     *              BSNGLA)
       IF (IELT.EQ.3)
     *  CALL XMABM7(KLST,KLCOLA,KLLDA,KNA,KNEQA,NBLOCA,ICLRA,EM30,
     *              COEFST,BCONST,KABST,KABSTN,ICUBA,ISYMM,CARRST,
     *              BSNGLA)

99999 END 

C ----------------------------------------------------------------------
C Allocate and generate mass-matrices on all levels
C (i.e. level NLMIN..NLMAX).
C
C This will generate the mass matrices. The necessary memory is
C automatically allocated.
C 
C BLUMP  - whether to lump the mass matrices
C BDILMP - whether to perform diagonal lumping
C          (-> copy the diagonal of the real matrix instead of summing
C              all entries to form the diagonal)
C
C The resulting matrices will have double precision.
C ----------------------------------------------------------------------

      SUBROUTINE GENMM (BLUMP, BDILMP)
      
      IMPLICIT NONE

      INCLUDE 'cerr.inc'

      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgtria.inc'
      
      INCLUDE 'cbasicelem.inc'
      
      INCLUDE 'cinidat.inc'

C *** Names of matrices and vectors (for output and tracing)
      
      INCLUDE 'arnmassmat.inc'
      
C externals

C *** Coefficient of stiffness matrix, right hand side, exact solution
      DOUBLE PRECISION COEFST
      EXTERNAL COEFST
      
C *** definition of finite elements

      EXTERNAL E030,E031,EM30,EM31,E010
      
C configuration of the matrix blocks

      INTEGER NBLOCA
      PARAMETER (NBLOCA=1)
      INTEGER NBLA1      
      PARAMETER (NBLA1=NBLOCA*NNLEV)

C local variables: structure of bilinear form of Laplacian matrix:
      
      INTEGER KABSTN(NBLOCA),KABST(2,NNAB,NBLOCA)
      
      LOGICAL BCONST(NBLOCA)      
      DATA BCONST/.TRUE./
      SAVE BCONST
      
      LOGICAL BSNGLA(NBLOCA,NNLEV)
      DATA BSNGLA/NBLA1*.FALSE./      
      SAVE BSNGLA 
      
C parameters

      LOGICAL BLUMP, BDILMP
      
C local variables: backup of mass matrix handles

      INTEGER KLMH

C some more local variables...
      
      INTEGER ICLRA, ISYMM, II, INU, NU, ILD
      DOUBLE PRECISION DMH
      
C initialise generation-descriptors for generation of real mass matrices

      KABST (1,1,1)=1
      KABST (2,1,1)=1
      KABSTN(1)    =1

C Don't exploit any symmetry
       
      ISYMM = 0
      ICLRA=1

C Call the XMAB07 routine that builds all the matrices...

      IF (IELT.EQ.0)
     *  CALL XMAB07(KLM,KLCOLA,KLLDA,KNA,KNEQA,NBLOCA,ICLRA,E031,
     *              COEFST,BCONST,KABST,KABSTN,ABS(ICUBM),ISYMM,CARRM,
     *              BSNGLA)
      IF (IELT.EQ.1)
     *  CALL XMAB07(KLM,KLCOLA,KLLDA,KNA,KNEQA,NBLOCA,ICLRA,E030,
     *              COEFST,BCONST,KABST,KABSTN,ABS(ICUBM),ISYMM,CARRM,
     *              BSNGLA)
      IF (IELT.EQ.2)
     *  CALL XMABM7(KLM,KLCOLA,KLLDA,KNA,KNEQA,NBLOCA,ICLRA,EM31,
     *              COEFST,BCONST,KABST,KABSTN,ABS(ICUBM),ISYMM,CARRM,
     *              BSNGLA)
      IF (IELT.EQ.3) 
     *  CALL XMABM7(KLM,KLCOLA,KLLDA,KNA,KNEQA,NBLOCA,ICLRA,EM30,
     *              COEFST,BCONST,KABST,KABSTN,ABS(ICUBM),ISYMM,CARRM,
     *              BSNGLA)
      
C generate lumped mass matrices?

      IF (BLUMP) THEN

C perform the lumping level for level:
      
        DO II=NLMIN,NLMAX

C make a backup of the matrix handle, reinit the KLM array

          KLMH    = KLM(II)
          KLM(II) = 0
          
C allocate memory -> a vector for diagonal entries of the lumped matrices
          
          NU=KNU(II)
          CALL  ZNEW (NU,1,KLM(II),'VMASS ')
          IF (IER.NE.0) THEN
            KLM(II) = KLMH
            GOTO 99999
          END IF
          
C Build the elements of the matrix (--> diagonal entries), either
C by taking the diagonal of the mass matrix or by summing the elements
C in a line

          DO INU=1,KNU(II)

            IF (BDILMP) THEN
              DMH=0D0
              DO ILD=KWORK(L(KLLDA(II))+INU-1),KWORK(L(KLLDA(II))+INU)-1
                DMH=DMH+DWORK(L(KLMH)+ILD-1)
              END DO
            ELSE
              ILD=KWORK(L(KLLDA(II))+INU-1)
              DMH=DWORK(L(KLMH)+ILD-1)
            ENDIF

C store the diagonal entry

            DWORK(L(KLM(II))+INU-1)=DMH
            
          END DO

C lumped matrix has been formed, we can release the origional one

          CALL  ZDISP (0, KLMH,'VMASSH')
          IF (IER.NE.0) GOTO 99999

        END DO

      END IF

99999 END 

C ----------------------------------------------------------------------
C Allocate and generate/calculate velocity right-hand-sides 
C on level ILV.
C
C The memory for the RHS vector(s) is automatically allocated.
C The handles for the U and V velocity components of the vector
C are returned in LF1 and LF2. The caller is responsible for
C freeing the memory allocated by this.
C KNEQA(ILV) has to be initialized! The current level will be switched
C to ILV if necessary!
C ----------------------------------------------------------------------

      SUBROUTINE GENRHS (ILV,LF1,LF2)
      
      IMPLICIT NONE

      INCLUDE 'cerr.inc'

      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgtria.inc'
      
      INCLUDE 'cbasicelem.inc'
      
      INCLUDE 'cinidat.inc'

C externals

C *** Coefficient of stiffness matrix, right hand side, exact solution
      
      DOUBLE PRECISION RHS
      EXTERNAL RHS
      
C *** definition of finite elements
      
      EXTERNAL E030,E031,EM30,EM31,E010
      
C parameters
       
      INTEGER ILV,LF1,LF2

C constants for proper operation

      INTEGER NBLOCF
      PARAMETER (NBLOCF=2)

C Names of vectors

      CHARACTER ARRDF*6
      DIMENSION ARRDF(NBLOCF)
      DATA ARRDF/'DF1   ','DF2   '/
      SAVE ARRDF

C constants for forming the linear form

      INTEGER KFN(NBLOCF),KF(NNAB,NBLOCF)
      
      LOGICAL BSNGLF(NBLOCF)      
      DATA BSNGLF /.FALSE.,.FALSE./
      SAVE BSNGLF

      LOGICAL BCONF(NBLOCF)
      DATA BCONF /.FALSE.,.FALSE./
      SAVE BCONF
      
C local variables

      INTEGER LF(NBLOCF),ICLRF,NEQ,I,J

C init the temp. array
      LF(1)=0
      LF(2)=0
      ICLRF=1
      
C get the number of equations in each vector from the Common block

      NEQ = KNEQA(ILV)
      
C initialize bilinear form
      
      KFN(1)=1
      KFN(2)=1
      
      DO J=1,NBLOCF
        DO I = 1,NNAB
          KF(I,J) = 0
        END DO
      END DO
      KF(1,1)=1
      KF(1,2)=1

C switch to level ILV if necessary; to ensure proper operation of
C XVB0 -> that uses variables in Common blocks /TRIAx/ !

      CALL SWTLEV (ILV)

C generate right hand side(s); store their handles in LF(..).
      
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
      
C return the handles
      
      LF1=LF(1)
      LF2=LF(2)
      
99999 END 

C ----------------------------------------------------------------------
C Allocate and generate the B-matrices of the system matrix
C on all levels.
C
C BEXACT=false - the matrix is constructed by using quadrature rules.
C BEXACT=true  - use exact evaluation for constructing the matrix
C                (only possible when using constant pressure/linear
C                 velocities as there is an exact quadrature rule!)
C Is BEXACT=true, the grid arrays KNU,KNB,KNEL,KNVT,KNMT,KVERT,...
C have to be initialised.
C
C The resulting matrices will have double precision.
C ----------------------------------------------------------------------

      SUBROUTINE GENMTB (BEXACT)
      
      IMPLICIT NONE

      INCLUDE 'cerr.inc'

      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgtria.inc'
      
      INCLUDE 'cbasicelem.inc'
      
      INCLUDE 'cinidat.inc'
      
C names of matrices/vectors

      INCLUDE 'arnbmat.inc'

C constants

      INTEGER NBLOCB,NBLB1
      PARAMETER (NBLOCB=2,NBLB1=NBLOCB*NNLEV)

C parameters

      LOGICAL BEXACT

C externals

C *** Coefficient of stiffness matrix, right hand side, exact solution
      DOUBLE PRECISION COEFST,COEFFB,RHS,UE
      EXTERNAL COEFST,COEFFB,RHS,UE
      
C *** definition of finite elements

      EXTERNAL E030,E031,EM30,EM31,E010
      
C local variables to define the bilinear form

      INTEGER KLB(NBLOCB,NNLEV)
      DATA KLB /NBLB1*0/
      SAVE KLB

      LOGICAL BCONB(NBLOCB)
      DATA BCONB /.TRUE.,.TRUE./
      SAVE BCONB
      
      INTEGER KABBN(NBLOCB),KABB(2,NNAB,NBLOCB)
      DATA KABBN/1,1/
      SAVE KABBN

      LOGICAL BSNGLB(NBLOCB,NNLEV)
      DATA BSNGLB/NBLB1*.FALSE./
      SAVE BSNGLB

C Constants for output

      CHARACTER ARRDB*6      
      DIMENSION ARRDB(NBLOCB)           
      DATA ARRDB/'DB1   ','DB2   '/
      SAVE ARRDB
      
C local variables

      INTEGER ICLRB,I,J,K
      
      IF (.NOT.BEXACT) THEN

C Use standard integration technique to build the B-matrices.
C Initialize bilinear form...

        DO K=1,NBLOCB
          DO J=1,NNAB
            DO I=1,2
              KABB(I,J,K)=0
            END DO
          END DO
        END DO

C ( A         B1 )
C (      A    B2 )
C ( B1^T B2^T 0  )


C B1-block
      
        KABB(1,1,1)=2
        KABB(2,1,1)=1
      
C B2 block

        KABB(1,1,2)=3
        KABB(2,1,2)=1

C As the B-matrices of the system matrix is not square, we can't
C use matrix storage format 7; so we use the more general 
C matrix format 9 where 7 is a special case of.

C Call the appropriate generation routines to generate the matrix;
C first generate the structure...

        IF ((IELT.EQ.0).OR.(IELT.EQ.2)) 
     *     CALL XMAP9(KLCOLB,KLLDB,KNB,KNEQB,E031,E010)
        IF ((IELT.EQ.1).OR.(IELT.EQ.3)) 
     *     CALL XMAP9(KLCOLB,KLLDB,KNB,KNEQB,E030,E010)
        IF (IER.NE.0) GOTO 99999

C then allocate and generate the matrix itself

        ICLRB=1
      
        IF ((IELT.EQ.0).OR.(IELT.EQ.2)) 
     *   CALL XMAB09(KLB,KLCOLB,KLLDB,KNB,NBLOCB,ICLRB,E031,E010,E010,
     *               COEFFB,BCONB,KABB,KABBN,ICUBB,CARRDB,BSNGLB)
        IF ((IELT.EQ.1).OR.(IELT.EQ.3)) 
     *   CALL XMAB09(KLB,KLCOLB,KLLDB,KNB,NBLOCB,ICLRB,E030,E010,E010,
     *               COEFFB,BCONB,KABB,KABBN,ICUBB,CARRDB,BSNGLB)
        IF (IER.NE.0) GOTO 99999
        
C write the handles of the matrices to the Common blocks for later use
        
        DO I=NLMIN,NLMAX
          
          KLB1(I)=KLB(1,I)
          KLB2(I)=KLB(2,I)

        END DO

      ELSE
      
C Use "quick and (not so) dirty" exact quadrature rules for
C constant pressure / linear velocities to calculate the B-matrices
C directly.
      
C loop through all levels

        DO I=NLMIN,NLMAX
      
          KNEQB(I)=KNU(I)
          KNB  (I)=2*KNU(I)
        
C allocate vectors for the matrix structure
        
          CALL ZNEW(KNB(I)  ,1,KLB(1,I) ,ARRDB(1))
          IF (IER.NE.0) GOTO 99999
          CALL ZNEW(KNB(I)  ,1,KLB(2,I) ,ARRDB(2))
          IF (IER.NE.0) GOTO 99999
          CALL ZNEW(KNB(I)  ,3,KLCOLB(I),'KCOLB ')
          IF (IER.NE.0) GOTO 99999
          CALL ZNEW(KNU(I)+1,3,KLLDB (I),'LLDB  ')
          IF (IER.NE.0) GOTO 99999

C build the matrix directly
          
          CALL BBUILD(KWORK(L(KLVERT(I))),KWORK(L(KLMID(I))),
     *             KWORK(L(KLADJ(I))),DWORK(L(KLCVG(I))),
     *             DWORK(L(KLB(1,I))),DWORK(L(KLB(2,I))),
     *             KWORK(L(KLCOLB(I))),KWORK(L(KLLDB(I))),
     *             KNB(I),KNEL(I),KNVT(I),KNMT(I))

C Free unused memory; there are only KNB(I) elements in KLB, KLCOLB,...

          CALL ZDISP (KNB(I),KLB(1,I) ,ARRDB(1))
          CALL ZDISP (KNB(I),KLB(2,I) ,ARRDB(2))
          CALL ZDISP (KNB(I),KLCOLB(I),'KCOLB ')
          IF (IER.NE.0) GOTO 99999
        
C write the handles of the matrices to the Common blocks for later use
          
          KLB1(I)=KLB(1,I)
          KLB2(I)=KLB(2,I)

        END DO

      END IF
      
99999 END 

C ----------------------------------------------------------------------
C Change matrix precision on level ILV.
C
C This will change the data type of the stokes-/mass-/B-matrix
C on level ILV depending on the parameters:
C
C IDTP  - the desired data type: 0=double prec, 1=single prec.,
C                                2=integer
C BST   - change data type of Stokes matrix to single precision
C BBMAT - change data type of B-matrices to single precision
C BMASS - change data type of mass matrix to single precision
C ----------------------------------------------------------------------

      SUBROUTINE CHMATP (ILV,IDTP,BST,BBMAT,BMASS)
      
      IMPLICIT NONE

      INCLUDE 'cerr.inc'

      INCLUDE 'cmem.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      
      INCLUDE 'cnsparfrac.inc'

      INTEGER ILV,IDTP
      LOGICAL BST,BBMAT,BMASS
      
C *** Names of matrices and vectors (for output and tracing)

      INCLUDE 'arnstokes.inc'
      INCLUDE 'arnmassmat.inc'
      INCLUDE 'arnbmat.inc'

      IF (BST) THEN
        CALL ZCTYPE(IDTP,KLST(ILV),CARRST)
        IF (IER.NE.0) GOTO 99999
      END IF

      IF (BMASS) THEN
        CALL ZCTYPE(IDTP,KLM(ILV),CARRM)
        IF (IER.NE.0) GOTO 99999
      END IF

      IF (BBMAT) THEN
        CALL ZCTYPE(IDTP,KLB1(ILV),CARRDB)
        IF (IER.NE.0) GOTO 99999
        CALL ZCTYPE(IDTP,KLB2(ILV),CARRDB)
        IF (IER.NE.0) GOTO 99999
      END IF

99999 END 

C ----------------------------------------------------------------------
C Update solution vector and/or RHS vector on level ILV.
C
C Implement pressure drop and dirichlet boundary conditions.
C
C BPDBCR  - whether to implement pressure drop into RHS-vectors
C BDBCR   - whether to implement dirichlet boundary conditions
C           into RHS vectors
C BDBCS   - whether to implement dirichlet boundary conditions
C           into solution vectors
C ----------------------------------------------------------------------

      SUBROUTINE IMPSLR (ILV, BPDBCR, BDBCR, BDBCS)
      
      IMPLICIT NONE

      INCLUDE 'cerr.inc'

      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgtria.inc'
      
C parameters

      INTEGER ILV
      LOGICAL BPDBCR, BDBCR, BDBCS
      
C externals

      DOUBLE PRECISION PARX,PARY,TMAX
      EXTERNAL PARX,PARY,TMAX
      EXTERNAL UE
      
C local variables
      
      INTEGER KU1,KU2,LF12P
      
      LF12P=KLF12P(ILV)

      IF (BPDBCR) THEN
C implement pressure drop values into right hand side
        CALL PDSET (KWORK(L(KLVBD(ILV))),KWORK(L(KLEBD(ILV))),
     *             KWORK(L(KLVERT(ILV))),KWORK(L(KLMID(ILV))),
     *             KWORK(L(KLNPR(ILV))),DWORK(L(KLCVG(ILV))),
     *             DWORK(L(KLMBDP(ILV))),
     *             DWORK(L(LF12P)),DWORK(L(LF12P)+KNEQA(ILV)),
     *             1D0,KNVBD(ILV))
      END IF

C implement Dirichlet values into solution vectors and/or vectors of
C right hand side:

      KU1=L(KLUP(ILV))
      KU2=KU1+KNEQA(ILV)
      
      IF (BDBCS) THEN
        CALL BDRSET (DWORK(KU1),DWORK(KU2),
     *             KWORK(L(KLNPR(ILV))),KNVT(ILV),PARX,PARY,UE,
     *             DWORK(L(KLCVG(ILV))),KWORK(L(KLADJ(ILV))),KNEL(ILV),
     *             KWORK(L(KLMID(ILV))),KWORK(L(KLVERT(ILV))))
      END IF
      
      IF (BDBCR) THEN
        CALL BDRSET (DWORK(L(LF12P)),
     *             DWORK(L(LF12P)+KNEQA(ILV)),
     *             KWORK(L(KLNPR(ILV))),KNVT(ILV),PARX,PARY,UE,
     *             DWORK(L(KLCVG(ILV))),KWORK(L(KLADJ(ILV))),KNEL(ILV),
     *             KWORK(L(KLMID(ILV))),KWORK(L(KLVERT(ILV))))
      END IF

99999 END 

C ----------------------------------------------------------------------
C Multiple prolongation
C
C This routine will take the solution vektor KLUP(ISTRT) on level
C ISTRT and prolongate it through the levels up to level IEND.
C The resulting vector is written to KLUP(IEND)
C (or DWORK(L(KLUP(IEND))) respectively).
C ----------------------------------------------------------------------

      SUBROUTINE MLTPRL (ISTART, IEND)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgtria.inc'

C parameters 

      INTEGER ISTART,IEND

C local variables

      INTEGER ILV
      INTEGER KU1C,KU2C,KPC
      INTEGER KVERTC,KMIDC ,KADJC,KU1,KU2,KP 

C perform the loop for the levels up to IEND
      
      DO ILV=ISTART+1,IEND
      
        KVERTC=L(KLVERT(ILV-1))
        KMIDC =L(KLMID(ILV-1))
        KADJC =L(KLADJ(ILV-1))
        
        KU1C=L(KLUP(ILV-1))
        KU2C=KU1C+KNEQA(ILV-1)
        KPC =KU2C+KNEQA(ILV-1)

        KU1   =L(KLUP(ILV))
        KU2   =KU1+KNEQA(ILV)
        KP    =KU2+KNEQA(ILV)

        CALL  PROLU (DWORK(KU1C),DWORK(KU2C),DWORK(KPC),
     *               DWORK(KU1),DWORK(KU2),DWORK(KP),
     *               KWORK(KVERTC),KWORK(KMIDC),KWORK(KADJC),
     *               KNEQA(ILV-1), KNEL(ILV-1), KNVT(ILV-1),
     *               KWORK(L(KLVERT(ILV))),KWORK(L(KLMID(ILV))),
     *               KWORK(L(KLADJ(ILV))),
     *               KNEQA(ILV), KNEL(ILV), KNVT(ILV),
     *               DWORK(L(KLCVG(ILV-1))),DWORK(L(KLAREA(ILV-1))))
     
      END DO
      
      END
      
C ----------------------------------------------------------------------
C Generate/read parametrisation and triangulation
C
C BPAR   - Generate information about parametrisation
C IPAR   - 0=FEAT parametrisation, 1=OMEGA parametrisation from file
C CPARM  - Filename of the file containing the parametrisation
C          (if IPAR>0)
C ----------------------------------------------------------------------

      SUBROUTINE GENPAR (BPAR, IPAR, CPARM)
      
      IMPLICIT NONE

      INCLUDE 'cerr.inc'

      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgtria.inc'

      INCLUDE 'cinidat.inc'
      
C parameters

      LOGICAL BPAR
      INTEGER IPAR
      CHARACTER*(*) CPARM
      
C externals
      
      INTEGER MFILE
      
      MFILE = 65
      
      IF (BPAR.AND.(IPAR.EQ.1)) THEN
C Read the domain parametrisation from file
       CALL RDPARM (CPARM,MFILE)
       CLOSE(MFILE)
      ENDIF

99999 END

C ----------------------------------------------------------------------
C Generate/read triangulation
C
C BTRI   - Generate information about triangulation
C IRGMSH - 0=create mesh information for all levels
C            This will also calculate information about midpoints and
C            elements meeting on each vertex.
C          1=read mesh information for all levels from file(s)
C            in unformatted form
C          2=read mesh information for all levels from file(s)
C            in formatted form
C IPREF  - If IRGMSH=0, this specifies a number of pre-refinements
C          of the coarsest level to obtain level 1.
C CMESH  - Filename of the file containing information about triang.
C          (if IRGMSH>0)
C ----------------------------------------------------------------------

      SUBROUTINE GENORS (BTRI, IRGMSH, IPREF, CMESH)
      
      IMPLICIT NONE

      INCLUDE 'cerr.inc'

      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgtria.inc'

      INCLUDE 'cinidat.inc'
      
C parameters

      LOGICAL BTRI
      INTEGER IRGMSH,IPREF
      CHARACTER*(*) CMESH
      
C externals
      
C *** Parametrization of the domain
      DOUBLE PRECISION PARX,PARY,TMAX
      EXTERNAL PARX,PARY,TMAX
C *** Control of refinement - here: regular refinement
      EXTERNAL S2DI0,S2DB0

C local variables
      
      INTEGER IMID,IADJ,IVEL,IDISP,IBDP
      INTEGER NVT,NMT,NEL,NVBD,NUP
      INTEGER MFILE,II
      
      MFILE = 65
      
C Always read coarse grid data from unit MMESH1, filename CMESH1

      IF (BTRI) THEN

        CALL  XORSC (MFILE,CMESH)
        IF (IER.NE.0) GOTO 99999
        CLOSE(MFILE)

C Create subdivisions of mesh for level 1..NLEV. Either generate them
C from the coarse grid or read them from hard disc in formatted
C or unfoirmatted form.

        IF (IRGMSH.EQ.0) THEN
        
C Generate subdivisions directly.
C IMID=2 to calculate midpoints,
C IVEL=1 to calculate the elements meeting on a vertex.
C Important for force calculations and grid deformation!

          IMID=2
          IADJ=1
          IVEL=1
          IDISP=1
          IBDP=2
          
C IPREF specifies the global number of pre-refinements of the coarse-grid
C to obtain level 1:

          CALL XMSB2(IMID,IADJ,IVEL,IDISP,IBDP,S2DI0,S2DB0,
     *               PARX,PARY,TMAX,IPREF)
          IF (IER.NE.0) GOTO 99999
          
        ELSE
          CALL XMORS2(IRGMSH)
          IF (IER.NE.0) GOTO 99999
        ENDIF

        DO II=NLMIN,NLMAX
C Update additional data in the multigrid COMMON blocks describing
C the triangulation        
          NVT=KNVT(II)
          NMT=KNMT(II)
          NEL=KNEL(II)
          NVBD=KNVBD(II)
          
          KNU(II)=NMT
          KNP(II)=NEL
          NUP=NMT+NMT+NEL
          KNUP(II)=NUP
        END DO
      END IF
      
99999 END

C ----------------------------------------------------------------------
C Implement boundary conditions into grid on level ILV
C
C This is a wrapper routine for implementing boundary conditions into
C the grid information. This includes e.g. the points where Neumann
C boundary condition occur, moving boundaries and so on.
C Basically this routine is only a wrapper routine for BDRNEU,
C but it updates the grid information in the Common blocks, too,
C depending on the results of BDRNEU.
C
C For proper operation this will switch the current level to ILV
C if necessary.
C ----------------------------------------------------------------------

      SUBROUTINE IMPBDG (ILV)
      
      IMPLICIT NONE

      INCLUDE 'cerr.inc'

      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgtria.inc'

      INCLUDE 'cinidat.inc'
      
C parameters

      INTEGER ILV

C local variables

      INTEGER NMBD

      CALL SWTLEV (ILV)

C Implement boundary information; calculate NMBD.

      CALL BDRNEU(KWORK(L(KLMBD(ILV))),KWORK(L(KLVBD(ILV))),
     *              KWORK(L(KLEBD(ILV))),KWORK(L(KLVERT(ILV))),
     *              KWORK(L(KLMID(ILV))),KWORK(L(KLNPR(ILV))),
     *              DWORK(L(KLDBD(ILV))),DWORK(L(KLMBDP(ILV))),
     *              DWORK(L(KLCVG(ILV))),KWORK(L(KLADJ(ILV))),
     *              KNVBD(ILV),KNVT(ILV),KNEL(ILV),
     *              NMBD,INEUM)
      IF (IER.NE.0) GOTO 99999
     
C update the multigrid Common block information
      
      KNMBD(ILV)=NMBD      

99999 END


C ----------------------------------------------------------------------
C Implement matrix restriction
C
C This routine capsules matrix restriction for the system matrix.
C The matriux entries belonging to cells with too high aspect ratio 
C on level ILV are rebuilt by restriction of the matrix of level
C ILV+1. The routine has no function when called on level NLMAX.
C ----------------------------------------------------------------------

      SUBROUTINE IMPMRS (ILV)
      
      IMPLICIT NONE

      INCLUDE 'cerr.inc'

      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmgtria.inc'

      INCLUDE 'cinidat.inc'
      INCLUDE 'cinidat2.inc'
      
C parameters

      INTEGER ILV
      
C local variables

      INTEGER KA1,KA2,KCOLA1,KCOLA2,KLDA1,KLDA2,LLV2,LLV1,LLA2,LLA1
      INTEGER LLM2,LLM1,NVT2,NVT1,NEL2,NEL1,NMT2,NMT1,LLAR1,LLAR2
      INTEGER LCORVG
      
      INTEGER IADM1
      
C Call matrix restriction routine to rebuild the matrix
C entrys belonging to anisotropic cells. This procedure is applied to all
C matrices except for the finest level

      IF (ILV.NE.NLMAX) THEN

C Initialize some helper variables
        KA1   =L(KLA(ILV))
        KA2   =L(KLA(ILV+1))
        KCOLA1=L(KLCOLA(ILV))
        KCOLA2=L(KLCOLA(ILV+1))
        KLDA1 =L(KLLDA(ILV))
        KLDA2 =L(KLLDA(ILV+1))
        LLV2=L(KLVERT(ILV+1))
        LLV1=L(KLVERT(ILV))
        LLA2=L(KLADJ(ILV+1))
        LLA1=L(KLADJ(ILV))
        LLM2=L(KLMID(ILV+1))
        LLM1=L(KLMID(ILV))
        NVT2=KNVT(ILV+1)
        NVT1=KNVT(ILV)
        NEL2=KNEL(ILV+1)
        NEL1=KNEL(ILV)
        NMT2=KNMT(ILV+1)
        NMT1=KNMT(ILV)
        LLAR2=L(KLAREA(ILV+1))
        LLAR1=L(KLAREA(ILV))
        LCORVG=KLCVG(ILV)

C Evaluate the configuration bitfield how to perform the 
C matrix modification

        IADM1=0
        IF(IAND(IAPRM,4).EQ.4) IADM1=1
        IF(IAND(IAPRM,36).EQ.36) IADM1=2

        IF ((IADM1.NE.0).AND.(DMTEP.GE.0D0))
     *    CALL  MAREST(KWORK(LLV1),KWORK(LLV2),KWORK(LLM1),KWORK(LLM2),
     *               KWORK(LLA1),KWORK(LLA2),
     *               KWORK(KLDA1),KWORK(KLDA2),
     *               KWORK(KCOLA1),KWORK(KCOLA2),
     *               DWORK(KA1),DWORK(KA2),
     *               DWORK(L(LCORVG)),DWORK(LLAR1),
     *               NEL1,NEL2,NVT1,NVT2,
     *               DMTEP,IADM1)

      END IF
      
99999 END
     
      