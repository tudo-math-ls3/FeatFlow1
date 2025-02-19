************************************************************************
* RHS and matrix construction routine 1 
*
* Recombinate the linear part of the RHS for the Theta-scheme.
* for one macro-step. Combinate Stokes- and Mass-matrix to build the
* linear part of the nonlinear system matrix.
*
* This routine is called in MGSTP to assemble mass-/system
* matrices for the nonlinear iteration.
*
* In:
*   KST1   - starting address in DWORK of array [1..?] of double.
*            This "Stokes"-matrix represents the linear / Laplacian
*            part of the nonlinear matrix
*              N(u) = -nu * Laplace (.)  +  u * grad(.)
*            This precalculated matrix is used as a bases for N that
*            way that other terms are added to it to build
*            the nonlinear matrix.
*   KM1    - starting address in DWORK of array [1..?] of double.
*            This is the precalculated mass-matrix of the system
*   KA1    - starting address in DWORK of array [1..NA] of double.
*            Will be filled with data about the linear part of the
*            nonlinear system matrix for the Theta-scheme:
*                     [ M*I + THSTEP * (-nu*Laplace) ] 
*   KCOLA,
*   KLDA   - Structure of the system matrix
*   KF1,
*   KF2    - starting address in DWORK of array [1..NU] of double.
*            This is the the X- and Y-component of the right-hand-side
*            vector of the equation
*   KU1,
*   KU2    - starting address in DWORK of array [1..NU] of double.
*            This is the the X- and Y-component of the current solution
*            vector u^n of the equation
*   THSTEP - double. Theta-scheme identifier.
*            =0: Forward Euler, =1: Backward Euler, =0.5: Crank Nicolson
*            For Fractional-Step theta scheme this is set to different
*            values, depending on the current step in the scheme.
*
* In (from COMMON-blocks):
*   IMASS  - =1: if we have a real mass matrix in the same structure
*                as the system matrix
*            =0: if we have a lumped mass matrix with entries only on
*                the diagonal
*   IPRECA - Determines whether the linear part of the system matrix/RHS
*            is constructed with precalculated data.
*            If IPRECA=4, the Stokes-part -nu*Laplace(u) is not added
*             to KA1. If there is a lumped mass matrix, that part of KA1
*             is build and a part of the RHS-vector is constructed:
*               RHS = M*u + KF,
*               KA1 = M*I
*             If the mass matrix is not lumped, KA1 is filled with 0 and
*             the RHS-vector KFx is not touched.
*            If IPRECA=2,3, KM1/KST1 is read from disc before
*             reconstruction.   
*
* Out:
*   DWORK(KA1..) - array for the nonlinear system matrix, 
*                  filled with data about the linear part
*   DWORK(KF1..),
*   DWORK(KF2..)
*          - By linear combination with the system matrix, these vectors
*            are updated now to represent the "linear part" of the right 
*            hand side of the Theta-scheme:
*                  KF = [M*u-THSTEP*{-nu*Laplace(u)}] + KF
************************************************************************
      SUBROUTINE XMADF1(KM1,KST1,KA1,KCOLA,KLDA,KF1,KF2,KU1,KU2,
     *                  NA,NU,THSTEP)

      IMPLICIT NONE

C main COMMON blocks

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      
      INCLUDE 'cinidat.inc'

C parameters

      INTEGER KM1,KST1,KA1,KCOLA,KLDA
      INTEGER KF1,KF2,KU1,KU2,NA,NU
      DOUBLE PRECISION THSTEP

C constants

      CHARACTER CFILST*12,CFILM*12,CFILE*12
      DIMENSION CFILST(NNLEV),CFILM(NNLEV)
      DATA CFILST/'#ns/ST1     ','#ns/ST2     ','#ns/ST3     ',
     *            '#ns/ST4     ','#ns/ST5     ','#ns/ST6     ',
     *            '#ns/ST7     ','#ns/ST8     ','#ns/ST9     '/
      DATA CFILM /'#ns/MA1     ','#ns/MA2     ','#ns/MA3     ',
     *            '#ns/MA4     ','#ns/MA5     ','#ns/MA6     ',
     *            '#ns/MA7     ','#ns/MA8     ','#ns/MA9     '/
      SAVE CFILST,CFILM,CFILE

C local variables

      INTEGER INU, INUA, INA

C     Standard matrix assembling branch:

      IF ((IPRECA.EQ.0).OR.(IPRECA.EQ.1)) THEN
      
C         We have to build the linear part of the nonlinear system 
C         matrix with the help of the Stokes- and mass-matrices:
C
C       KA1 = [  M*I  +  THSTEP * N(u)  ]
C           = [  M*I  +  THSTEP * (-nu * Laplace(.))  +  u * grad(.)  ]
C                ^^^              ^^^^^^^^^^^^^^^^^     ^^^^^^^^^^^
C               ->KM1                ->KST1             -> ignored here
C
C       Check if we have a real or a lumped mass matrix.

        IF (IMASS.EQ.1) THEN
        
C         We have a real mass matrix in the structure of the 
C         system matrix.
C         Now check the Theta-parameter. If it's 0, we can skip the 
C         linear combination with the Stokes-matrix. Otherwise
C         simply add the matrices:
C
C                        KM1 + THSTEP*KST1
        
          IF (THSTEP.NE.0D0) THEN
            CALL LLC1(DWORK(KST1),DWORK(KA1),NA,THSTEP,0D0)
            CALL LLC1(DWORK(KM1),DWORK(KA1),NA,1D0,1D0)
          ELSE
            CALL LCP1(DWORK(KM1),DWORK(KA1),NA)
          ENDIF
          
        ELSE
        
C         Using lumped mass matrices we have only to tackle
C         the diagonal when adding it to the system matrix.
C         Again add the matrices:
C
C                        KM1 + THSTEP*KST1
        
          IF (THSTEP.NE.0D0) THEN
            CALL LLC1(DWORK(KST1),DWORK(KA1),NA,THSTEP,0D0)
            DO INU=1,NU
              INUA=KWORK(KLDA+INU-1)-1
              DWORK(KA1+INUA)=DWORK(KA1+INUA)+DWORK(KM1+INU-1)
            END DO
          ELSE
            CALL LCL1(DWORK(KA1),NA)
            DO INU=1,NU
              INUA=KWORK(KLDA+INU-1)-1
              DWORK(KA1+INUA)=DWORK(KM1+INU-1)
            END DO
          ENDIF
        ENDIF

C       Multiply the current solution- and RHS-vector to build
C       the "linear" part of the RHS-vector of the Theta-scheme:
C
C         RHS := [ M*I - THSTEP*(-nu*Laplace(u^n)) ]  +  RHS

        CALL  LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *               DWORK(KU1),DWORK(KF1),1D0,1D0)
        CALL  LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *               DWORK(KU2),DWORK(KF2),1D0,1D0)

      ENDIF

C=======================================================================

C     Check if the matrices are so big that they have to be read in
C     from a file. The duty of this branch is the same as above.

      IF ((IPRECA.EQ.2).OR.(IPRECA.EQ.3)) THEN

        IF (THSTEP.NE.0D0) THEN
          CALL  OF0 (59,CFILST(ILEV),0)
          CFILE='STMAT '
          CALL  ORA1 (DWORK(KA1),CFILE,59,0)
          REWIND(59)
          CLOSE (59)
          IF (IER.NE.0) RETURN
          CALL LLC1(DWORK(KA1),DWORK(KA1),NA,THSTEP,0D0)
        ELSE
          CALL LCL1(DWORK(KA1),NA)
        ENDIF

        IF (IMASS.EQ.1) THEN
          CALL  OF0 (59,CFILM(ILEV),0)
          CFILE='MASMAT'
          CALL  ORALC1 (DWORK(KA1),1D0,CFILE,59,0)
          REWIND(59)
          CLOSE (59)
          IF (IER.NE.0) RETURN
        ELSE
          IF (THSTEP.NE.0D0) THEN
            DO INU=1,NU
              INUA=KWORK(KLDA+INU-1)-1
              DWORK(KA1+INUA)=DWORK(KA1+INUA)+DWORK(KM1+INU-1)
            END DO
          ELSE
            DO INU=1,NU
              INUA=KWORK(KLDA+INU-1)-1
              DWORK(KA1+INUA)=DWORK(KM1+INU-1)
            END DO
          ENDIF
        ENDIF

        CALL  LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *               DWORK(KU1),DWORK(KF1),1D0,1D0)
        CALL  LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *               DWORK(KU2),DWORK(KF2),1D0,1D0)

      ENDIF

C=======================================================================

C     In this branch we don't rely on any precalculated information
C     in KST1. If we have a lumped mass matrix, we can build a part of
C     the RHS-vector and of KA1 - otherwise simply clear KA1 as is
C     will be calculated later by the caller.

      IF (IPRECA.EQ.4) THEN

        IF (IMASS.EQ.0) THEN
          DO INU=1,NU
            DWORK(KF1+INU-1)= DWORK(KF1+INU-1)
     *                       +DWORK(KM1+INU-1)*DWORK(KU1+INU-1)
            DWORK(KF2+INU-1)= DWORK(KF2+INU-1)
     *                       +DWORK(KM1+INU-1)*DWORK(KU2+INU-1)
          END DO

          CALL LCL1(DWORK(KA1),NA)

          DO INU=1,NU
            INUA=KWORK(KLDA+INU-1)-1
            DWORK(KA1+INUA)=DWORK(KM1+INU-1)
          END DO
        ELSE
          CALL LCL1(DWORK(KA1),NA)
        ENDIF

      ENDIF

99998 END

************************************************************************
* Defect-RHS and matrix construction routine 2
*
* Build the linear part of the defect vector for the nonlinear 
* iteration. Combinate Stokes- and Mass-matrix to build the
* linear part of the nonlinear system matrix.
*
* This routine is called in NSDEF to assemble mass-/system
* matrices for the nonlinear iteration.
*
* In:
*   KST1   - starting address in DWORK of array [1..?] of double.
*            This "Stokes"-matrix represents the linear / Laplacian
*            part of the nonlinear matrix
*              N(u) = -nu * Laplace (.)  +  u * grad(.)
*            This precalculated matrix is used as a bases for N that
*            way that other terms are added to it to build
*            the nonlinear matrix.
*   KM1    - starting address in DWORK of array [1..?] of double.
*            This is the precalculated mass-matrix of the system
*   KA1    - starting address in DWORK of array [1..NA] of double.
*            Will be filled with data about the linear part of the
*            nonlinear system matrix for the Theta-scheme:
*                     [ M*I + THSTEP * (-nu*Laplace) ] 
*   KCOLA,
*   KLDA   - Structure of the system matrix
*   KB1,
*   KB2    - starting address in DWORK of array [1..*] of double.
*            This are the two "pressure" parts of the system matrix.
*                [ KST1        KB1 ] ( u1 )   ( d1 )
*                [       KST1  KB2 ] ( u2 ) = ( d2 )
*                [ KB1^T KB2^T  0  ] ( p  )   ( dp )
*   KU1,
*   KU2    
*   KP     - starting address in DWORK of array [1..NU] of double.
*            This is the the X-,Y- and pressure-component of the current 
*            solution vector u^n of the equation
*   KD1,
*   KD2,
*   KDP    - starting address in DWORK of array [1..*] of double.
*            KD1,KD2,KDP must contain the current RHS-vector of the linear 
*            system for velocity and pressure. The resulting defect
*            vector is written info KD1,KD2,KDP, thus overwriting 
*            KD1,KD2,KDP. 
*   THSTEP - double. Theta-scheme identifier.
*            =0: Forward Euler, =1: Backward Euler, =0.5: Crank Nicolson
*            For Fractional-Step theta scheme this is set to different
*            values, depending on the current step in the scheme.
*            For stationary simulations this parameter must be set
*            to 1 to include the full Laplacian matrix into the
*            system matrix:
*                    [alpha*M + THETA*KST1 ] 
*   ISTAT  - Whether the simulation is stationary or not.
*            As this routine generally builds
*                    [alpha*M + Theta_1*nu*k*L ] u
*            (the linear part of S(u)), this parameter decides whether
*            ALPHA=0 or ALPHA=1, so whether the mass matrix is added 
*            or not.   
*   INEUM  - Whether there are Neumann boundary components in our
*            problem or not.
*
* In (from COMMON-blocks):
*   IMASS  - =1: if we have a real mass matrix in the same structure
*                as the system matrix
*            =0: if we have a lumped mass matrix with entries only on
*                the diagonal
*   IPRECA - Determines whether the linear part of the system matrix/RHS
*            is constructed with precalculated data.
*            If IPRECA=4, the Stokes-part -nu*Laplace(u) is not added
*             to KA1. If there is a lumped mass matrix, that part of KA1
*             is build and a part of the RHS-vector is constructed:
*               RHS = M*u + KF,
*               KA1 = M*I
*             If the mass matrix is not lumped, KA1 is filled with 0 and
*             the RHS-vector KFx is not touched.
*            If IPRECA=2,3, KM1/KST1 is read from disc before
*             reconstruction.   
*
* Out:
*   DWORK(KA1..) - array for the nonlinear system matrix, 
*                  filled with data about the linear part
*   DWORK(KD1..),
*   DWORK(KD2..),
*   DWORK(KDP..)
*          - These arrays will receive the defect vector, calculated
*            with the RHS, the solution and the linear part of the
*            nonlinear system matrix:
*                  KD = KD - [alpha*M*u-THSTEP*{-nu*Laplace(u)}]*u
************************************************************************

      SUBROUTINE XMADF2(KM1,KST1,KA1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *                  KD1,KD2,KDP,KU1,KU2,KP,NA,NU,NP,KMBD,NMBD,INEUM,
     *                  THSTEP,ISTAT)
     
************************************************************************

      IMPLICIT NONE

C main COMMON blocks

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      
      INCLUDE 'cinidat.inc'

C parameters

      INTEGER KM1,KST1,KA1,KCOLA,KLDA,KB1,KB2,KCOLB,KLDB,
     *        KD1,KD2,KDP,KU1,KU2,KP,NA,NU,NP,KMBD(*),NMBD,INEUM,ISTAT
      INTEGER KF1,KF2
      DOUBLE PRECISION THSTEP

C constants

      CHARACTER CFILST*12,CFILM*12,CFILE*12
      DIMENSION CFILST(NNLEV),CFILM(NNLEV)
      DATA CFILST/'#ns/ST1     ','#ns/ST2     ','#ns/ST3     ',
     *            '#ns/ST4     ','#ns/ST5     ','#ns/ST6     ',
     *            '#ns/ST7     ','#ns/ST8     ','#ns/ST9     '/
      DATA CFILM /'#ns/MA1     ','#ns/MA2     ','#ns/MA3     ',
     *            '#ns/MA4     ','#ns/MA5     ','#ns/MA6     ',
     *            '#ns/MA7     ','#ns/MA8     ','#ns/MA9     '/
      SAVE CFILST,CFILM,CFILE

C local variables

      INTEGER INU, INUA, INA

C     Standard matrix assembling branch:

      IF ((IPRECA.EQ.0).OR.(IPRECA.EQ.1)) THEN

C       There's a slight difference to XMADF1 here - we have to
C       perform sifferent tasks whether we have a stationary
C       or an instationary simulation. In case that the simulation
C       is instationary, we go the same way as XMADF1:

        IF (ISTAT.EQ.1) THEN

C         We have to build the linear part of the nonlinear system 
C         matrix with the help of the Stokes- and mass-matrices:
C
C         KA1 = [  M*I  +  THSTEP * N(u)  ]
C             = [  M*I  +  THSTEP * (-nu * Laplace(.))  +  u * grad(.)  ]
C                  ^^^              ^^^^^^^^^^^^^^^^^     ^^^^^^^^^^^
C                 ->KM1                ->KST1             -> ignored here
C
C         Check if we have a real or a lumped mass matrix.

          IF (IMASS.EQ.1) THEN

C           We have a real mass matrix in the structure of the 
C           system matrix.
C           Now check the Theta-parameter. If it's 0, we can skip the 
C           linear combination with the Stokes-matrix. Otherwise
C           simply add the matrices:
C
C                          KM1 + THSTEP*KST1

            IF (THSTEP.NE.0D0) THEN
              CALL LLC1(DWORK(KST1),DWORK(KA1),NA,THSTEP,0D0)
              CALL LLC1(DWORK(KM1 ),DWORK(KA1),NA,1D0,1D0)
            ELSE
              CALL LCP1(DWORK(KM1),DWORK(KA1),NA)
            ENDIF
            
          ELSE
          
C           Using lumped mass matrices we have only to tackle
C           the diagonal when adding it to the system matrix.
C           Again add the matrices:
C
C                          KM1 + THSTEP*KST1

            IF (THSTEP.NE.0D0) THEN
              CALL LLC1(DWORK(KST1),DWORK(KA1),NA,THSTEP,0D0)
              DO INU=1,NU
                INUA=KWORK(KLDA+INU-1)-1
                DWORK(KA1+INUA)=DWORK(KA1+INUA)+DWORK(KM1+INU-1)
              END DO
            ELSE
              CALL LCL1(DWORK(KA1),NA)
              DO INU=1,NU
                INUA=KWORK(KLDA+INU-1)-1
                DWORK(KA1+INUA)=DWORK(KM1+INU-1)
              END DO
            ENDIF
          ENDIF
        ELSE
        
C         The case of an stationary simulation is much easier. 
C         Remember, we have to solve a system of the form
C
C             S(u_h)u_h  +  k B p_h  =  g_h,   B^T u_h = 0
C
C         where
C
C             S(u) = [alpha*M + Theta_1*nu*k*L + Theta_2*k*K(u)] u
C         
C         is the matrix whose linear parts we are building here.
C         In the stationary approach, there is alpha=0, so adding
C         the mass matrix can be dropped completely!
C
C         Therefore we only have to build
C
C         KA1 = [  THSTEP * (-nu * Laplace(.))  +  u * grad(.)  ]
C                           ^^^^^^^^^^^^^^^^^     ^^^^^^^^^^^
C                              ->KST1             -> ignored here
C
C         what can simply be done by copying KST1 to KA1 with the 
C         correct scaling factor.
        
          CALL LLC1(DWORK(KST1),DWORK(KA1),NA,THSTEP,0D0)
          
        ENDIF

C       Multiply the solution vector with the built matrix, and 
C       subtract it from the RHS-vector (the current setting of KDx)
C       to build the defect vector.  
C
C         DEF := DEF - [ ALPHA*M*I - THSTEP*(-nu*Laplace(u^n)) ] u

        IF (INLMAX.GT.1) THEN
          CALL  LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DWORK(KU1),DWORK(KD1),-1D0,1D0)
          CALL  LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DWORK(KU2),DWORK(KD2),-1D0,1D0)
        ENDIF

      ENDIF

C=======================================================================

C     Check if the matrices are so big that they have to be read in
C     from a file. The duty of this branch is the same as above.

      IF ((IPRECA.EQ.2).OR.(IPRECA.EQ.3)) THEN

        IF (ISTAT.EQ.1) THEN
          IF (THSTEP.NE.0D0) THEN
            CALL  OF0 (59,CFILST(ILEV),0)
            CFILE='STMAT '
            CALL  ORA1 (DWORK(KA1),CFILE,59,0)
            REWIND(59)
            CLOSE (59)
            IF (IER.NE.0) RETURN
            CALL LLC1(DWORK(KA1),DWORK(KA1),NA,THSTEP,0D0)
          ELSE
            CALL LCL1(DWORK(KA1),NA)
          ENDIF

          IF (IMASS.EQ.1) THEN
            CALL  OF0 (59,CFILM(ILEV),0)
            CFILE='MASMAT'
            CALL  ORALC1 (DWORK(KA1),1D0,CFILE,59,0)
            REWIND(59)
            CLOSE (59)
            IF (IER.NE.0) RETURN
          ELSE
            IF (THSTEP.NE.0D0) THEN
            DO INU=1,NU
              INUA=KWORK(KLDA+INU-1)-1
              DWORK(KA1+INUA)=DWORK(KA1+INUA)+DWORK(KM1+INU-1)
            END DO
            ELSE
            DO INU=1,NU
              INUA=KWORK(KLDA+INU-1)-1
              DWORK(KA1+INUA)=DWORK(KM1+INU-1)
            END DO
            ENDIF
          ENDIF
        ELSE
          CALL  OF0 (59,CFILST(ILEV),0)
          CFILE='STMAT '
          CALL  ORA1 (DWORK(KA1),CFILE,59,0)
          REWIND(59)
          CLOSE (59)
          IF (IER.NE.0) RETURN
          CALL LLC1(DWORK(KA1),DWORK(KA1),NA,THSTEP,0D0)
        ENDIF

        IF (INLMAX.GT.1) THEN
          CALL  LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DWORK(KU1),DWORK(KD1),-1D0,1D0)
          CALL  LAX17 (DWORK(KA1),KWORK(KCOLA),KWORK(KLDA),NU,
     *                 DWORK(KU2),DWORK(KD2),-1D0,1D0)
        ENDIF

      ENDIF

C=======================================================================

C     In this branch we don't rely on any precalculated information
C     in KST1. If we have a lumped mass matrix, we can build a part of
C     the RHS-vector and of KA1 - otherwise simply clear KA1 as is
C     will be calculated later by the caller.

      IF (IPRECA.EQ.4) THEN

        IF (ISTAT.EQ.1) THEN
          IF (IMASS.EQ.0) THEN
            IF (INLMAX.GT.1) THEN
              DO INU=1,NU
                DWORK(KD1+INU-1)= DWORK(KD1+INU-1)
     *                           -DWORK(KM1+INU-1)*DWORK(KU1+INU-1)
                DWORK(KD2+INU-1)= DWORK(KD2+INU-1)
     *                           -DWORK(KM1+INU-1)*DWORK(KU2+INU-1)
              END DO
            ENDIF

            CALL LCL1(DWORK(KA1),NA)

            DO INU=1,NU
              INUA=KWORK(KLDA+INU-1)-1
              DWORK(KA1+INUA)=DWORK(KM1+INU-1)
            END DO
          ELSE
            CALL LCL1(DWORK(KA1),NA)
          ENDIF
        ELSE
          CALL LCL1(DWORK(KA1),NA)
        ENDIF

      ENDIF

C     That's it for the velocity part of the system matrix
C
C=======================================================================
C
C     Now turn over to handle the pressure part. Up to now we
C     constructed the defect vector KD1/KD2 only using the velocity
C     part - but of course that's not enough. Also the pressure-part
C     has to be included into the defect calculation.
C
C     So now update the defect vector by
C
C         KD1 = KD1 - KB1*KP
C
C     to get the correct defect for the velocity. Remember, the
C     whole thing we are calculating here is:
C
C        (def1)   (d1)   [ KST1        KB1 ] ( u1 ) 
C        (def2) = (d2) - [       KST1  KB2 ] ( u2 ) 
C        (defp)   (dp)   [ KB1^T KB2^T  0  ] ( p  ) 
C
C     Standard handling: IPRECB<>2

      IF (IPRECB.NE.2) THEN
      
C       If there are Neumann boundary components, we first have
C       to modify the current defect vector: All entries corresponding
C       to Neumann nodes are multiplied by 1/2.
C
C       This compensates the BUG in the VANCA smoother, which accidently 
C       introduced the factor 2.0 in the smoothing/solving procedure!
      
        IF (INEUM.EQ.1) CALL  BDRDEF(DWORK(KD1),KMBD,NMBD,0.5D0)
        
C       Subtract KB1*p from KD1 to get the correct defect for U
        
        CALL  LAX19 (DWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,
     *               DWORK(KP),DWORK(KD1),-1D0,1D0)
     
C       Multiply the Neumann nodes by 2 again to compensate 
C       the above multiplication by 1/2.
     
        IF (INEUM.EQ.1) CALL  BDRDEF(DWORK(KD1),KMBD,NMBD,2.0D0)
        
C       The same treatment now for the V-variable of the defect:
        
        IF (INEUM.EQ.1) CALL  BDRDEF(DWORK(KD2),KMBD,NMBD,0.5D0)
        CALL  LAX19 (DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,
     *               DWORK(KP),DWORK(KD2),-1D0,1D0)
        IF (INEUM.EQ.1) CALL  BDRDEF(DWORK(KD2),KMBD,NMBD,2.0D0)
        
      ENDIF

C     Special treatment: IPRECB=2=elementwise application.
C     Might not be implemented completely...

      IF (IPRECB.EQ.2) THEN
        CALL BMUL1 (KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *              DWORK(L(LCORVG)),DWORK(KP),DWORK(KD1),DWORK(KD2),
     *              NEL,NVT,NMT,-1D0,1D0)
      ENDIF

C=======================================================================

C     Implement Dirichlet-values into the defect.
C     On all Dirichlet-nodes, the defect must be exactly 0.

      IF (INLMAX.GT.1) CALL BDRY0 (DWORK(KD1),DWORK(KD2),KMBD,NMBD)

C=======================================================================

C     Finally calculate the defect vector of the pressure by subtracting
C     the velocity components (multiplied with KB) from KDP.
C
C       KDP = KDP  -  KB1^T KU1  -  KB2^T KU2

      IF (IPRECB.NE.2) THEN
       CALL  LTX19 (DWORK(KB1),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU1),DWORK(KDP),-1D0,1D0)
       CALL  LTX19 (DWORK(KB2),KWORK(KCOLB),KWORK(KLDB),NU,NP,
     *              DWORK(KU2),DWORK(KDP),-1D0,1D0)
      ENDIF

C     Special treatment: IPRECB=2=elementwise application.
C     Might not be implemented completely...

      IF (IPRECB.EQ.2) THEN
       WRITE(6,*) 'BTMUL ???'
       CALL BTMUL1(KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LADJ)),
     *             DWORK(L(LCORVG)),DWORK(KU1),DWORK(KU2),DWORK(KDP),
     *             NEL,NVT,NMT,-1D0)
      ENDIF

99998 END

************************************************************************
* Matrix construction routine 3
*
* Combinate Stokes- and Mass-matrix to build the
* linear part of the nonlinear system matrix.
*
* This routine is called in NGSTP to assemble mass-/system
* matrices for the nonlinear iteration on all lower levels than NLMAX.
*
* In:
*   KST1   - starting address in DWORK of array [1..?] of double.
*            This "Stokes"-matrix represents the linear / Laplacian
*            part of the nonlinear matrix
*              N(u) = -nu * Laplace (.)  +  u * grad(.)
*            This precalculated matrix is used as a bases for N that
*            way that other terms are added to it to build
*            the nonlinear matrix.
*   KM1    - starting address in DWORK of array [1..?] of double.
*            This is the precalculated mass-matrix of the system
*   KA1    - starting address in DWORK of array [1..NA] of double.
*            Will be filled with data about the linear part of the
*            nonlinear system matrix for the Theta-scheme:
*                     [ M*I + THSTEP * (-nu*Laplace) ] 
*   KCOLA,
*   KLDA   - Structure of the system matrix
*   THSTEP - double. Theta-scheme identifier.
*            =0: Forward Euler, =1: Backward Euler, =0.5: Crank Nicolson
*            For Fractional-Step theta scheme this is set to different
*            values, depending on the current step in the scheme.
*            For stationary simulations this parameter must be set
*            to 1 to include the full Laplacian matrix into the
*            system matrix:
*                    [alpha*M + THETA*KST1 ] 
*   ISTAT  - Whether the simulation is stationary or not.
*            As this routine generally builds
*                    [alpha*M + Theta_1*nu*k*L ] u
*            (the linear part of S(u)), this parameter decides whether
*            ALPHA=0 or ALPHA=1, so whether the mass matrix is added 
*            or not.   
*
* In (from COMMON-blocks):
*   IMASS  - =1: if we have a real mass matrix in the same structure
*                as the system matrix
*            =0: if we have a lumped mass matrix with entries only on
*                the diagonal
*   IPRECA - Determines whether the linear part of the system matrix/RHS
*            is constructed with precalculated data.
*            If IPRECA=4, the Stokes-part -nu*Laplace(u) is not added
*             to KA1. If there is a lumped mass matrix, that part of KA1
*             is build and a part of the RHS-vector is constructed:
*               RHS = M*u + KF,
*               KA1 = M*I
*             If the mass matrix is not lumped, KA1 is filled with 0 and
*             the RHS-vector KFx is not touched.
*            If IPRECA=2,3, KM1/KST1 is read from disc before
*             reconstruction.   
*
* Out:
*   DWORK(KA1..) - array for the nonlinear system matrix, 
*                  filled with data about the linear part
************************************************************************

      SUBROUTINE XMADF3(KM1,KST1,KA1,KCOLA,KLDA,NA,NU,THSTEP,ISTAT)
      
************************************************************************
      
      IMPLICIT NONE

C main COMMON blocks

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      
      INCLUDE 'cinidat.inc'

C parameters

      INTEGER KM1,KST1,KA1,KCOLA,KLDA,NA,NU,ISTAT
      DOUBLE PRECISION THSTEP

C constants

      CHARACTER CFILST*12,CFILM*12,CFILE*12
      DIMENSION CFILST(NNLEV),CFILM(NNLEV)
      DATA CFILST/'#ns/ST1     ','#ns/ST2     ','#ns/ST3     ',
     *            '#ns/ST4     ','#ns/ST5     ','#ns/ST6     ',
     *            '#ns/ST7     ','#ns/ST8     ','#ns/ST9     '/
      DATA CFILM /'#ns/MA1     ','#ns/MA2     ','#ns/MA3     ',
     *            '#ns/MA4     ','#ns/MA5     ','#ns/MA6     ',
     *            '#ns/MA7     ','#ns/MA8     ','#ns/MA9     '/
      SAVE CFILST, CFILM, CFILE

C local variables

      INTEGER INU, INA, INUA

C     Standard matrix assembling branch:

      IF ((IPRECA.EQ.0).OR.(IPRECA.EQ.1)) THEN

C       There's a slight difference to XMADF1 here - we have to
C       perform sifferent tasks whether we have a stationary
C       or an instationary simulation. In case that the simulation
C       is instationary, we go the same way as XMADF1:

        IF (ISTAT.EQ.1) THEN

C         We have to build the linear part of the nonlinear system 
C         matrix with the help of the Stokes- and mass-matrices:
C
C         KA1 = [  M*I  +  THSTEP * N(u)  ]
C             = [  M*I  +  THSTEP * (-nu * Laplace(.))  +  u * grad(.)  ]
C                  ^^^              ^^^^^^^^^^^^^^^^^     ^^^^^^^^^^^
C                 ->KM1                ->KST1             -> ignored here
C
C         Check if we have a real or a lumped mass matrix.

          IF (IMASS.EQ.1) THEN

C           We have a real mass matrix in the structure of the 
C           system matrix.
C           Now check the Theta-parameter. If it's 0, we can skip the 
C           linear combination with the Stokes-matrix. Otherwise
C           simply add the matrices:
C
C                          KM1 + THSTEP*KST1

            IF (THSTEP.NE.0D0) THEN
              CALL LLC1(DWORK(KST1),DWORK(KA1),NA,THSTEP,0D0)
              CALL LLC1(DWORK(KM1),DWORK(KA1),NA,1D0,1D0)
            ELSE
              CALL LCP1(DWORK(KM1),DWORK(KA1),NA)
            ENDIF
          ELSE

C           Using lumped mass matrices we have only to tackle
C           the diagonal when adding it to the system matrix.
C           Again add the matrices:
C
C                          KM1 + THSTEP*KST1

            IF (THSTEP.NE.0D0) THEN
              CALL LLC1(DWORK(KST1),DWORK(KA1),NA,THSTEP,0D0)
              DO INU=1,NU
                INUA=KWORK(KLDA+INU-1)-1
                DWORK(KA1+INUA)=DWORK(KA1+INUA)+DWORK(KM1+INU-1)
              END DO
            ELSE
              CALL LCL1(DWORK(KA1),NA)
              DO INU=1,NU
                INUA=KWORK(KLDA+INU-1)-1
                DWORK(KA1+INUA)=DWORK(KM1+INU-1)
              END DO
            ENDIF
          ENDIF
        ELSE
        
C         The case of a stationary simulation is much easier. 
C         Remember, we have to solve a system of the form
C
C             S(u_h)u_h  +  k B p_h  =  g_h,   B^T u_h = 0
C
C         where
C
C             S(u) = [alpha*M + Theta_1*nu*k*L + Theta_2*k*K(u)] u
C         
C         is the matrix whose linear parts we are building here.
C         In the stationary approach, there is alpha=0, so adding
C         the mass matrix can be dropped completely!
C
C         Therefore we only have to build
C
C         KA1 = [  THSTEP * (-nu * Laplace(.))  +  u * grad(.)  ]
C                           ^^^^^^^^^^^^^^^^^     ^^^^^^^^^^^
C                              ->KST1             -> ignored here
C
C         what can simply be done by copying KST1 to KA1 with the 
C         correct scaling factor.

          CALL LLC1(DWORK(KST1),DWORK(KA1),NA,THSTEP,0D0)
          
        ENDIF

      ENDIF

C=======================================================================

C     Check if the matrices are so big that they have to be read in
C     from a file. The duty of this branch is the same as above.

      IF ((IPRECA.EQ.2).OR.(IPRECA.EQ.3)) THEN

        IF (THSTEP.NE.0D0) THEN
          CALL  OF0 (59,CFILST(ILEV),0)
          CFILE='STMAT '
          CALL  ORA1 (DWORK(KA1),CFILE,59,0)
          REWIND(59)
          CLOSE (59)
          IF (IER.NE.0) RETURN
          CALL LLC1(DWORK(KA1),DWORK(KA1),NA,THSTEP,0D0)
        ELSE
          CALL LCL1(DWORK(KA1),NA)
        ENDIF

        IF (IMASS.EQ.1) THEN
          CALL  OF0 (59,CFILM(ILEV),0)
          CFILE='MASMAT'
          CALL  ORALC1 (DWORK(KA1),1D0,CFILE,59,0)
          REWIND(59)
          CLOSE (59)
          IF (IER.NE.0) RETURN
        ELSE
          IF (THSTEP.NE.0D0) THEN
            DO INU=1,NU
              INUA=KWORK(KLDA+INU-1)-1
              DWORK(KA1+INUA)=DWORK(KA1+INUA)+DWORK(KM1+INU-1)
            END DO
          ELSE
            DO INU=1,NU
              INUA=KWORK(KLDA+INU-1)-1
              DWORK(KA1+INUA)=DWORK(KM1+INU-1)
            END DO
          ENDIF
        ENDIF

      ENDIF
C=======================================================================

C     In this branch we don't rely on any precalculated information
C     in KST1. If we have a lumped mass matrix, we can build a part of
C     the RHS-vector and of KA1 - otherwise simply clear KA1 as is
C     will be calculated later by the caller.

      IF (IPRECA.EQ.4) THEN

        IF (IMASS.EQ.0) THEN
          CALL LCL1(DWORK(KA1),NA)
          DO INU=1,NU
            INUA=KWORK(KLDA+INU-1)-1
            DWORK(KA1+INUA)=DWORK(KM1+INU-1)
          END DO
        ELSE
          CALL LCL1(DWORK(KA1),NA)
        ENDIF

      ENDIF

C=======================================================================

99998 END

************************************************************************
* Computes the norms RESU, RESDIV
*
* This routine computes the "U-residual" RESU and the "divergence
* residual" RESDIV.
*
* In:
*   U1,
*   U2,
*   P      - Velocity/pressure vector of the solution of 
*            the linearised system    A * (U1,U2,P) = (F1,F2,FP)
*   D1,
*   D2,
*   DP     - Defect vector   (D1,D2,DP) = (F1,F2,FP) - A*(U1,U2,P)
*   F1,
*   F2,
*   FP     - RHS vector of the system   A * (U1,U2,P) = (F1,F2,FP)
*   NU     - length of the velocity vectors
*   NP     - length of the pressure vectors
* 
* Out:
*   RESU   - Norm of the velocity residual:
*                               || (D1,D2) ||_E
*               RESU   = -----------------------------
*                        max ( ||F1||_E , ||F2||_E )
*
*   RESDIV - Norm of the pressure residual
*                           || P ||_E
*               RESDIV = ----------------
*                        || (U1,U2) ||_E
************************************************************************

      SUBROUTINE RESDFK(U1,U2,P,D1,D2,DP,F1,F2,FP,NU,NP,RESU,RESDIV)

      IMPLICIT NONE
      
      INCLUDE 'ctria.inc'

C parameters      

      INTEGER NU,NP
      DOUBLE PRECISION U1(*),U2(*),P(*),D1(*),D2(*),DP(*),F1(*),F2(*),
     *                 FP(*)
      DOUBLE PRECISION RESU, RESDIV

C local variables
      
      DOUBLE PRECISION RESF, RESF1, RESF2, RESU1, RESU2
      DOUBLE PRECISION DNORMU, DNRMU1, DNRMU2

C-----------------------------------------------------------------------
C     Compute the relative l2-norms  RESU,RESDIV
C-----------------------------------------------------------------------

C     RESF := max ( ||F1||_E , ||F2||_E )

      CALL LL21 (F1,NU,RESF1)
      CALL LL21 (F2,NU,RESF2)
      RESF=MAX(RESF1,RESF2)
      IF (ABS(RESF).LT.1D-8) RESF=1D0

C                   || (D1,D2) ||_E
C     RESU = -----------------------------
C            max ( ||F1||_E , ||F2||_E )

      CALL LL21 (D1,NU,RESU1)
      CALL LL21 (D2,NU,RESU2)
      RESU=SQRT(RESU1*RESU1+RESU2*RESU2)/RESF

C     DNORMU = || (U1,U2) ||_l2 

      CALL LL21 (U1,NU,DNRMU1)
      CALL LL21 (U2,NU,DNRMU2)
      DNORMU=SQRT(DNRMU1*DNRMU1+DNRMU2*DNRMU2)
      IF (ABS(DNORMU).LT.1D-8) DNORMU=1D0

C                 || P ||_E
C     RESDIV = ----------------
C              || (U1,U2) ||_E

      CALL LL21 (DP,NP,RESDIV)
      RESDIV=RESDIV/DNORMU
C
      END

