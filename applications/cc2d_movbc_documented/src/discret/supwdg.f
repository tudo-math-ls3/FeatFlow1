************************************************************************
* Implementation of the streamline diffusion technique
************************************************************************

      SUBROUTINE SUPWDG(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                  KCOLA,KLDA,KVERT,KMID,DCORVG,ELE,COEFFN,
     *                  IDEF,DCMASS)

************************************************************************
*     Purpose: -  Adds the SUPG-part on matrix block A after
*                 it was initialized by the linear part or
*              -  builds up the complete nonlinear system matrix
*              -  The input vector Ui is the old velocity field
*              -  The input vectors UjLi are the transport directions
*
*     PARAMETRIC VERSION
*
* In:
*   A, NA,
*   KCOLA,
*   LLDA   - array [1..*] of double/integer
*            Structure arrays of system matrix, 
*            maybe initialised with linear parts.  
*   KVERT,
*   KMID,
*   DCORVG - usual geometry information
*
*   U1L1,
*   U1L2   - array [1..NU] of double
*            main velocity field used for assembling of 
*            the nonlinearity. Can be undefined if A1L=0.
*   U2L1,
*   U2L2   - array [1..NU] of double
*            secondary velocity field, used for the assembling
*            of the nonlinearity. Can be undefined if A2L=0.
*   A1L    - double; weighting factor for U1L1/U1L2
*   A2L    - double; weighting factor for U2L1/U2L2
*
*   IDEF   - Controls the behaviour of this routine.
*            =0: modify system matrix, add nonlinearity.
*                Defect vectors D1,D2 and velocity vectors U1,U2
*                can be undefined.
*            =1: modify both, system matrix and defect vector
*            =2: modify defect vectors, include nonlinearity.
*                A can be undefined (but not KCOLA,KLDA!)
*
*   U1,
*   U2     - array [1..NU] of double
*            Solution vector for modifying the defect vector.
*            Only if IDEF=1,2, otherwise it can be undefined.
*   D1,
*   D2     - array [1..NU] of double
*            Defect vector, modified by U1,U2.
*            Only if IDEF=1,2, otherwise it can be undefined.
*   
*   DCMASS - This defines how to modify the system matrix - if the
*            matrix is modified at all (see IDEF!).
*             =0: subtract the mass matrix from the system matrix;
*                 The nonlinearity is not build.
*            <>0: Add the nonlinear term including the stabilization
*                 to the system matrix.
*                 If the full matrix with not-lumped mass matrix 
*                 is to be build (IPRECA=4 and IMASS=1),
*                 add DCMASS*(Mass matrix) to the system matrix.
*
*   ELE    - used element for the discretisation
*   COEFFN - Coefficient function defining the coefficients of
*            the bilinear form of the convectuve block 
*            (for future enhancements; currently not used)
*
* In (from COMMON-blocks in former times):
*
*   THSTEP - Current Step-size of the Theta-scheme
*
* In (from COMMON-blocks):
*
*   UPSAM  - control parameter
*   RE     - 1/nu = viscosity
*   ISTOK  - =1 if the Stokes equation is to be discretized
*            rather than Navier-Stokes.
*   IPRECA - if =4, the complete system matrix (including the linear
*            parts) is build
*   IMASS  - only if IPRECA=4; must be 1 if the real mass matrix
*            should be used
*   Remark: IPRECA=4 and IMASS=0 does not work! In this case, no mass
*   matrix is added to the system matrix!
*
* Out:
*   A      - system matrix;
*            the nonlinearity is added to that matrix
*   D1,
*   D2     - Modified defect vector; only if IDEF=1,2.
*
* Remarks:
*  
* 1.) In a typical call of the upwinding, the caller can use:
*     A1L = 1, U1L1/U1L2 = velocity field
*     A2L = 0, U2L1/U2L2 = undefined
*   So the upwinding scheme only uses one velocity field.
*   Such a call e.g. adds the integral
*                ( U1Lx*grad(.) , v )_Omega
*   to the system matrix.
*
*  2.) In case that there are two velocity fields representing
*   the solution (may happen in a nonstationary simulation where
*   U1L1/U1L2 represents the solution in the current and U2L1/U2L2
*   that of the previous time step), A1L/A2L defines how these both
*   velocity vectors should be weighted to compute the actual
*   velocity field for the assembling:
*                U_act = A1L*U1Lx + A2L*U2Lx
*   This is e.g. used for the linear extrapolation technique to
*   reconstruct a velocity from two previous time steps...
*
*  3.) In the nonlinear iteration, as a right hand side there arises
*   a defect vector D, which linear part can easily being assembled.
*   However, there is a nonlinearity to be included into that vector,
*   too. By setting IDEF=1,2, this routine incorporates the nonlinearity
*   into that vector, using the formula
*
*             D = D - THSTEP * UUx * grad (Ux)
*
*   If IPRECA=4, IMASS=1, DCMASS<>0, D is updated according to
*
*             D = D - DCMASS*M*UUx - THSTEP * UUx * grad (Ux)
*   
*   If IPRECA=4, IMASS=1, DCMASS=0, D is updated according to
*
*             D = D + M*UUx - THSTEP * UUx * grad (Ux)
*
*  4.) Currently, NEQ=NMT is assumed, thus this routine only works
*      with E030/E031. This can easily be corrected if necessary...
************************************************************************

      IMPLICIT NONE
      
C include necessary COMMON blocks

      INCLUDE 'cout.inc' 
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      INCLUDE 'ccub.inc'
      
      INCLUDE 'casmbly.inc'
      
      INCLUDE 'cinidat.inc'
      
      INCLUDE 'cnsparfrac.inc'
      
C parameters
      
      DOUBLE PRECISION A(*)
      DOUBLE PRECISION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*)
      DOUBLE PRECISION D1(*),D2(*),DCORVG(2,*)
      DOUBLE PRECISION DENTRY(NNBAS,NNBAS)
      DOUBLE PRECISION A1L,A2L,DCMASS
      INTEGER NA
      INTEGER KCOLA(*),KLDA(*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*)
      INTEGER KENTRY(NNBAS,NNBAS)
      INTEGER IDEF
      
      DOUBLE PRECISION COEFFN
      EXTERNAL ELE, COEFFN

C externals

      INTEGER NDFL
      EXTERNAL NDFL

C local variables

      INTEGER ICUB,I,IELTYP,IEQ,JDOFE,ILD,JCOL,JCOL0,IDOFE,IDFG
      INTEGER IVE,JP,JDFL,JDFG,IDOFEH,JDOFEH,IA
      DOUBLE PRECISION DNY,CT0,DUMAX,DU1,DU2,DUNORM,DELTA
      DOUBLE PRECISION DJ1,DJ2,DJ3,DJ4
      DOUBLE PRECISION XI1,XI2
      DOUBLE PRECISION OM,HBAS,HBASJ1,HBASJ2,HBASJ3,HSUMJ,AH
      DOUBLE PRECISION HBASI1,HBASI2,HBASI3,HSUMI,DENTH

C     In case that IPRECA=4, create the Laplace-part of the system
C     and include it into the system matrix each time this routine
C     is called. Otherwise don't create the Laplace part, just
C     add the nonlinearity.
C     This handling is simply realized by setting the factor NU in
C     front of the Laplace part to 0.

      IF (IPRECA.EQ.4) THEN
        DNY=NY
      ELSE
        DNY=0D0
      ENDIF

C     A similar handling holds for the case that the (full!) mass
C     matrix is to be included into the system matrix while the
C     matrix is to be rebuild completely. 
C     DCMASS represents the factor in front of the mass matrix.
C     Note that we cannot use DCMASS directly! In the Theta scheme,
C     we typically calculate 
C         [ DCMASS*M + THSTEP*(Laplace) + ... ]
C     so we have to weight everything except for the mass matrix!
C     We make a little trick here to realize that. We always weight
C     everything by THSTEP including the mass matrix - but divide
C     DCMASS by THSTEP before to compensate that!
C         THSTEP * [ CT0*M + (Laplace) + ... ]
C     with CT0 = DCMASS/THSTEP.
C
C     If only the nonlinear part is to be build (which is the normal
C     case, as rebuilding the comlete matrix including the full mass
C     matrix is normally not necessary), simply set CT0=0 which 
C     prevents calculating the mass matrix.

      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
        CT0=DCMASS/THSTEP
      ELSE
        CT0=0D0
      ENDIF

C     Initialize BDER for our element. We want the element to calculate
C     function values as well as first X/Y derivatives:

      DO I = 1,NNDER
        BDER(I)=.FALSE.
      END DO

      DO I=1,3
        BDER(I)=.TRUE.
      END DO

C     Ask the element about its type:

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      
C     Get the local number of degrees of freedom:
      
      IDFL=NDFL(IELTYP)
      
C     Initialize the cubature formula identifyer in the COMMON block
C     with our cubature formula we have to use here:

      ICUB=ICUBN
      CALL CB2Q(ICUB)
      
      IF (IER.NE.0) GOTO 99999

************************************************************************
C *** Calculation of the matrix - storage technique 7 or 8
************************************************************************

C     Initialize the values in the cubature points. Allowes the element
C     to precalculate some information for faster access later

      ICUBP=ICUB
      CALL ELE(0D0,0D0,-2)

C     Calculate the maximum norm of the actual velocity field
C     U = A1*U1 + A2*U2 into DUMAX. 
C     Round up the norm to 1D-8 if it's too small...

      DUMAX=0D0
      IF (A2L.EQ.0D0) THEN
        DO IEQ=1,NMT
          DU1=U1L1(IEQ)
          DU2=U1L2(IEQ)
          DUNORM=SQRT(DU1**2+DU2**2)
          DUMAX=MAX(DUMAX,DUNORM)
        END DO
      ELSE       
        DO IEQ=1,NMT
          DU1=A1L*U1L1(IEQ)+A2L*U2L1(IEQ)
          DU2=A1L*U1L2(IEQ)+A2L*U2L2(IEQ)
          DUNORM=SQRT(DU1**2+DU2**2)
          DUMAX=MAX(DUMAX,DUNORM)
        END DO
      ENDIF       

      IF (DUMAX.LT.1D-8) DUMAX=1D-8

C *** Loop over all elements

      DO IEL=1,NEL
      
C       Calculate the local and global degrees of freedom on the
C       current element IEL:
      
        CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
      	IF (IER.LT.0) GOTO 99999

C       Determine local DELTA for streamline-diffusion
C       (cf. p. 131/132 (121) in Turek's CFD book).
C
C       For Stokes flow, we have the equation
C
C                -nu*div(u) + grad(p) = f
C       
C       not containing any convective part. Therefore we don't need
C       any stabilization technique. So in this case, switch of the
C       stabilization by setting DELTA to 0:

        IF (ISTOK.NE.1) THEN
          CALL DELTSD(U1L1,U1L2,U2L1,U2L2,A1L,A2L,IEL,DUMAX,DELTA,
     *                KVERT,KMID,DCORVG)
        ELSE
          DELTA=0D0
        END IF

C       Determine entry positions in matrix.
C
C       Here we build a small matrix DENTRY/KENTRY which
C       corresponds to the current element. We will assemble the
C       contributions of our current element into this matrix and
C       will add the result into the main matrix later.
C
C       To successfully integrate the contributions into the main
C       matrix, we compute in advance the positions in the main
C       matrix where we have to add the contribution to.
C       KENTRY(X,Y) corresponds to the index in the array A
C       of the entry in line Y, row KCOL(line start + X).

        DO JDOFE=1,IDFL
          ILD=KLDA(KDFG(JDOFE))
          KENTRY(JDOFE,JDOFE)=ILD
          DENTRY(JDOFE,JDOFE)=0D0
          JCOL0=ILD
          DO IDOFE=1,IDFL
            IF (IDOFE.NE.JDOFE) THEN
              IDFG=KDFG(IDOFE)
              DO JCOL=JCOL0,NA
                IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
              END DO
113           JCOL0=JCOL+1
              KENTRY(JDOFE,IDOFE)=JCOL
              DENTRY(JDOFE,IDOFE)=0D0
            END IF
          END DO
        END DO

C       Calculate auxiliary Jacobian factors for the transformation 
C       onto the reference element. See QTRAFO.F for details...

        DO IVE = 1, NVE
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
        END DO

        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))

C       Loop over the cubature points in our element to calculate
C       its contribution to each of the degrees of freedom on our
C       element

        DO ICUBP = 1, NCUBP
        
C         Get the coordinates of the cubature point on the reference
C         element
        
          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)

C         Calculate the Jacobian of the bilinear mapping onto the
C         reference element and the weight OM of the cubature formula:

          DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
          DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
          DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
          DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
          DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
          OM=DOMEGA(ICUBP)*DETJ

C         Call the element to calculate the values in the current
C         cubature point on the reference element:

          CALL ELE(XI1,XI2,-3)
          IF (IER.LT.0) GOTO 99999

C         Now we have to assemble the "local" matrix DENTRY/KENTRY.
C         This assembling decomposes now into different parts,
C         depending on what has to me assembled.
C
C         We want to set up the nonlinear part of the matrix
C
C           n~_h (u_h, u_h, v_h) 
C
C         = n_h (u_h, u_h, v_h) + sum_T ( delta_T ( u_h*grad u_h, u_h*grad v_h)_T )
C           ^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C          standard nonlin. part                  stabilization
C
C         More precisely, as we want to assemble the matrix which is 
C         later multiplied with coefficient vectors, we have to insert
C         basis functions in the above terms instead of u_h and v_h.
C         Assuming the representation u_h=sum_j(u_j*Phi_j) and 
C         v_h=sum_i(u_i,Phi_i), the above term is evaluated in the
C         DOF's as:
C         
C           n_h (u_h, Phi_j, Phi_i) 
C         + sum_T ( delta_T ( u_h*grad Phi_j, u_h*grad Phi_i )_T )
C
C         In nonstationary simulations, the system matrix typically
C         contains a mass matrix to respect the time derivative.
C         The matrix has the form
C
C         [  dcmass*M*I  +  THWEIG * (-nu * Laplace(.))  ] + THWEIG * u grad(.)
C
C         so if DCMASS<>0, incorporate the (real, not lumped!) mass matrix
C         into the local matrix.

          IF (DCMASS.NE.0D0) THEN

C           Calculate the actual velocity in the current cubature point
C           into (DU1,DU2). If we only have a primary velocity field
C           (A2L=0), we can calculate that only by summing up the
C           velocities in U1Lx, otherwise we have to sum up
C           A1*U1Lx + A2*U2Lx.

            DU1=0D0
            DU2=0D0
            IF (A2L.EQ.0D0) THEN
              DO JDFL=1,IDFL
                HBAS=DBAS(KDFL(JDFL),1)
                IF (ABS(HBAS).GE.1D-8) THEN
                  JDFG=KDFG(JDFL)
                  DU1=DU1+U1L1(JDFG)*HBAS
                  DU2=DU2+U1L2(JDFG)*HBAS
                ENDIF
              END DO
            ELSE
              DO JDFL=1,IDFL
                HBAS=DBAS(KDFL(JDFL),1)
                IF (ABS(HBAS).GE.1D-8) THEN
                  JDFG=KDFG(JDFL)
                  DU1=DU1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
                  DU2=DU2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
                ENDIF
              END DO
            ENDIF

C           We take a more detailed look onto the last scalar product
C           of n~_h (u_h, u_h, v_h) what we want to calculate here.
C
C           The vector u_h=(DU1,DU2) contains both velocity components,
C           for the X as well as for the Y velocity. On the other hand
C           the system matrix we want to build here will be designed for 
C           one velocity component only! Therefore, Phi_i and Phi_j
C           are scalar functions, so grad(Phi_i), grad(Phi_j) are vectors
C           with two components. Therefore, the last scalar product is more 
C           in detail:
C
C               ( u_h*grad Phi_j, u_h*grad Phi_i )_T
C
C           =   ( < (DU1) , (grad(Phi_j)_1) > , < (DU1) , (grad(Phi_i)_1) > )_T
C                   (DU2) , (grad(Phi_j)_2)       (DU2) , (grad(Phi_i)_2)  
C
C           =   < (DU1) , (grad(Phi_j)_1) >  *  < (DU1) , (grad(Phi_j)_1) >
C                 (DU2) , (grad(Phi_j)_2)         (DU2) , (grad(Phi_j)_2)
C
C           =   HSUMJ * HSUMI
C
C           i.e. a product of two scalar values!
C
C           Summing up over all pairs of multiindices.
C
C           Outer loop over the DOF's j=1..IDFL on our current element, 
C           which corresponds to the basis functions Phi_j:

            DO JDOFE=1,IDFL
            
C             Fetch the contributions of the basis functions Phi_j for
C             function value and first derivatives for the current
C             DOF into HBASxy:
            
              JDOFEH=KDFL(JDOFE)
              HBASJ1=DBAS(JDOFEH,1)
              HBASJ2=DBAS(JDOFEH,2)
              HBASJ3=DBAS(JDOFEH,3)

C             Calculate 
C
C                 U * grad(Phi_j)  =  < grad(Phi_j), U >
C     
C               = ( grad(Phi_j)_1 , (DU1) )
C                 ( grad(Phi_j)_2   (DU2) )
              
              HSUMJ = HBASJ2*DU1 + HBASJ3*DU2

C             Inner loop over the DOF's i=1..IDFL, which corresponds to
C             the basis function Phi_i:

              DO IDOFE=1,IDFL
              
                IF (IDOFE.EQ.JDOFE) THEN
                
C                 Short version of the evaluation of the matrix
C                 contribution - see below for a more detailed
C                 description what is added together here!
                
                  AH = HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *               + DNY*(HBASJ2**2+HBASJ3**2)
     *               + CT0*HBASJ1**2
     
                ELSE
                
C                 Fetch the contributions of the basis function Phi_i for
C                 function value and first derivatives for the current
C                 DOF into HBASIy:
                
                  IDOFEH=KDFL(IDOFE)
                  HBASI1=DBAS(IDOFEH,1)
                  HBASI2=DBAS(IDOFEH,2)
                  HBASI3=DBAS(IDOFEH,3)

C                 Calculate 
C
C                     U * grad(Phi_i)  =  < grad(Phi_i), U >
C     
C                   = ( grad(Phi_i)_1 , (DU1) )
C                     ( grad(Phi_i)_2   (DU2) )

                  HSUMI=HBASI2*DU1+HBASI3*DU2
     
C                 Finally calculate the controbution to the system
C                 matrix. Depending on the configuration of DNY,
C                 IPRECA, DCMASS,... this decomposes into three
C                 different parts:
C
C                 AH = n~_h(u_h,phi_j,phi_i)        | nonlinear part
C                    + DNY*(grad(phi_j,grad(phi_i)) | -nu*Laplace(u)
C                    + CT0*(phi_j*phi_i)            | Mass matrix
C
C                 The last two parts are probably not added to the
C                 matrix by setting DNY or CT0 to 0, respectively.
C
C                 For saving some numerical operations, we write:
C     
C                     HSUMI * (Delta * HSUMJ + HBASJ1)
C
C                 =   Delta * HSUMI * HSUMJ  
C                   + HSUMI * HBASJ1
C     
C                 =   Delta * ( U*grad(Phi_j), U*grad(Phi_i) )
C                   + (grad(Phi_i)*U,Phi_j)
C
C               <->   sum_T ( delta_T ( u_h*grad Phi_j, u_h*grad Phi_i) )_T
C                   + n_h (u_h, Phi_j, Phi_i)
                  
                  AH = HSUMI*(DELTA*HSUMJ+HBASJ1)
     *               + DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *               + CT0*HBASJ1*HBASI1
     
                ENDIF ! (IDOFE.EQ.JDOFE)

C               Weighten the calculated value AH by the cubature
C               weight OM and add it to the local matrix. After the
C               loop over all DOF's is finished, each entry contains
C               the calculated integral.

                DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
                
              END DO ! IDOFE
              
            END DO ! JDOFE

          ELSE

C           Coefficient in front of the mass matrix is 0.
C           Subtract the mass matrix from the system matrix.
C
C           Outer loop over the DOF's j=1..IDFL corresponding
C           to Phi_j:

            DO JDOFE=1,IDFL
            
C             Fetch the function value in the current DOF:
            
              HBASJ1=DBAS(KDFL(JDOFE),1)

C             Inner loop over the DOF's i=1..IDFL corresponding
C             to Phi_i:

              DO IDOFE=1,IDFL
                
C               Fetch the function value of the other DOF of the
C               current element
                
                HBASI1=DBAS(KDFL(IDOFE),1)
                
C               Calculate the contribution for the entry. The factor
C               THSTEP is compensated later when the local matrix
C               is included into the global matrix and/or in the
C               modification of the RHS vector.
                
                AH=-1D0/THSTEP*HBASJ1*HBASI1
                
C               ...and incorporate it into the local matrix.
                
                DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
                
              END DO ! IDOFE
              
            END DO ! JDOFE

          END IF ! (DCMASS.NE.0D0)

        END DO ! ICUBP
      
C       Now we have set up a "local" system matrix. We can either
C       include it into the real matrix or we can use it to simply
C       modify the RHS vector to create a defect vector (throwing
C       away the information about the matrix afterwards, which would
C       result in a matrix free modification of the RHS vector).
      
        DO JDOFE=1,IDFL
        
          DO IDOFE=1,IDFL
          
C           Get the entry from the local matrix and weight it according
C           to the current THETA in the Theta scheme, given by the
C           parameter. 
C           (Remark: For stationary simulations, THSTEP is typically
C            1D0 which includes the local matrix into the global one
C            directly).
          
            DENTH=THSTEP*DENTRY(JDOFE,IDOFE)

C           For IDEF=0/1, incorporate our "local" system matrix into 
C           the global matrix. The position of each entry DENTRY(X,Y) 
C           in the global matrix array A was saved in element KENTRY(X,Y)
C           before.

            IF (IDEF.LT.2) THEN
              IA   =KENTRY(JDOFE,IDOFE)
              A(IA)=A(IA)+DENTH
            ENDIF

C           For IDEF=1,2, build the defect vector
C               D = RHS - A*U
C           This is done matrix free, only with the help of the local 
C           matrix.
C           In this case, D=(D1,D2) is expected to be the RHS on
C           entry and will be updated to be the defect vector when
C           this routine is left.

            IF (IDEF.GT.0) THEN 
              IDFG=KDFG(IDOFE)
              JDFG=KDFG(JDOFE)
              D1(JDFG)= D1(JDFG)-DENTH*U1(IDFG)
              D2(JDFG)= D2(JDFG)-DENTH*U2(IDFG)
            ENDIF 

          END DO ! IDOFE
          
        END DO ! JDOFE

C       Matrix/defect vector updated for that element; proceed to the
C       next one...

      END DO ! IEL

99999 END

************************************************************************

      SUBROUTINE SUPWNP(U1L1,U1L2,U2L1,U2L2,A1L,A2L,U1,U2,D1,D2,A,NA,
     *                  KCOLA,KLDA,KVERT,KMID,DCORVG,ELE,COEFFN,
     *                  IDEF,DCMASS)

************************************************************************
*     Purpose: -  Adds the SUPG-part on matrix block A after
*                 it was initialized by the linear part or
*              -  builds up the complete nonlinear system matrix
*              -  The input vector Ui is the old velocity field
*              -  The input vectors UjLi are the transport directions
*
*     NONPARAMETRIC VERSION
*
* In:
*   A, NA,
*   KCOLA,
*   LLDA   - array [1..*] of double/integer
*            Structure arrays of system matrix, 
*            maybe initialised with linear parts.  
*   KVERT,
*   KMID,
*   DCORVG - usual geometry information
*
*   U1L1,
*   U1L2   - array [1..NU] of double
*            main velocity field used for assembling of 
*            the nonlinearity. Can be undefined if A1L=0.
*   U2L1,
*   U2L2   - array [1..NU] of double
*            secondary velocity field, used for the assembling
*            of the nonlinearity. Can be undefined if A2L=0.
*   A1L    - double; weighting factor for U1L1/U1L2
*   A2L    - double; weighting factor for U2L1/U2L2
*
*   IDEF   - Controls the behaviour of this routine.
*            =0: modify system matrix, add nonlinearity.
*                Defect vectors D1,D2 and velocity vectors U1,U2
*                can be undefined.
*            =1: modify both, system matrix and defect vector
*            =2: modify defect vectors, include nonlinearity.
*                A can be undefined (but not KCOLA,KLDA!)
*
*   U1,
*   U2     - array [1..NU] of double
*            Solution vector for modifying the defect vector.
*            Only if IDEF=1,2, otherwise it can be undefined.
*   D1,
*   D2     - array [1..NU] of double
*            Defect vector, modified by U1,U2.
*            Only if IDEF=1,2, otherwise it can be undefined.
*   
*   DCMASS - This defines how to modify the system matrix - if the
*            matrix is modified at all (see IDEF!).
*             =0: subtract the mass matrix from the system matrix;
*                 The nonlinearity is not build.
*            <>0: Add the nonlinear term including the stabilization
*                 to the system matrix.
*                 If the full matrix with not-lumped mass matrix 
*                 is to be build (IPRECA=4 and IMASS=1),
*                 add DCMASS*(Mass matrix) to the system matrix.
*
*   ELE    - used element for the discretisation
*   COEFFN - Coefficient function defining the coefficients of
*            the bilinear form of the convectuve block 
*            (for future enhancements; currently not used)
*
* In (from COMMON-blocks in former times):
*
*   THSTEP - Current Step-size of the Theta-scheme
*
* In (from COMMON-blocks):
*
*   UPSAM  - control parameter
*   RE     - 1/nu = viscosity
*   ISTOK  - =1 if the Stokes equation is to be discretized
*            rather than Navier-Stokes.
*   IPRECA - if =4, the complete system matrix (including the linear
*            parts) is build
*   IMASS  - only if IPRECA=4; must be 1 if the real mass matrix
*            should be used
*   Remark: IPRECA=4 and IMASS=0 does not work! In this case, no mass
*   matrix is added to the system matrix!
*
* Out:
*   A      - system matrix;
*            the nonlinearity is added to that matrix
*   D1,
*   D2     - Modified defect vector; only if IDEF=1,2.
*
* Remarks:
*  
* 1.) In a typical call of the upwinding, the caller can use:
*     A1L = 1, U1L1/U1L2 = velocity field
*     A2L = 0, U2L1/U2L2 = undefined
*   So the upwinding scheme only uses one velocity field.
*   Such a call e.g. adds the integral
*                ( U1Lx*grad(.) , v )_Omega
*   to the system matrix.
*
*  2.) In case that there are two velocity fields representing
*   the solution (may happen in a nonstationary simulation where
*   U1L1/U1L2 represents the solution in the current and U2L1/U2L2
*   that of the previous time step), A1L/A2L defines how these both
*   velocity vectors should be weighted to compute the actual
*   velocity field for the assembling:
*                U_act = A1L*U1Lx + A2L*U2Lx
*   This is e.g. used for the linear extrapolation technique to
*   reconstruct a velocity from two previous time steps...
*
*  3.) In the nonlinear iteration, as a right hand side there arises
*   a defect vector D, which linear part can easily being assembled.
*   However, there is a nonlinearity to be included into that vector,
*   too. By setting IDEF=1,2, this routine incorporates the nonlinearity
*   into that vector, using the formula
*
*             D = D - THSTEP * UUx * grad (Ux)
*
*   If IPRECA=4, IMASS=1, DCMASS<>0, D is updated according to
*
*             D = D - DCMASS*M*UUx - THSTEP * UUx * grad (Ux)
*   
*   If IPRECA=4, IMASS=1, DCMASS=0, D is updated according to
*
*             D = D + M*UUx - THSTEP * UUx * grad (Ux)
*
*  4.) Currently, NEQ=NMT is assumed, thus this routine only works
*      with EM30/EM31. This can easily be corrected if necessary...
************************************************************************

      IMPLICIT NONE
      
C include necessary COMMON blocks

      INCLUDE 'cout.inc' 
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'
      INCLUDE 'ccub.inc'
      
      INCLUDE 'casmbly.inc'
      
      INCLUDE 'cinidat.inc'

      INCLUDE 'cnsparfrac.inc'
      
C parameters
      
      DOUBLE PRECISION A(*)
      DOUBLE PRECISION U1L1(*),U1L2(*),U2L1(*),U2L2(*),U1(*),U2(*)
      DOUBLE PRECISION D1(*),D2(*),DCORVG(2,*)
      DOUBLE PRECISION DENTRY(NNBAS,NNBAS)
      DOUBLE PRECISION A1L,A2L,DCMASS
      INTEGER NA
      INTEGER KCOLA(*),KLDA(*)
      INTEGER KVERT(NNVE,*),KMID(NNVE,*)
      INTEGER KENTRY(NNBAS,NNBAS)
      INTEGER IDEF
      
      DOUBLE PRECISION COEFFN
      EXTERNAL ELE, COEFFN

C externals

      INTEGER NDFL
      EXTERNAL NDFL

C local variables

      INTEGER ICUB,I,IELTYP,IEQ,JDOFE,ILD,JCOL,JCOL0,IDOFE,IDFG
      INTEGER IVE,JP,JDFL,JDFG,IDOFEH,JDOFEH,IA
      DOUBLE PRECISION DNY,CT0,DUMAX,DU1,DU2,DUNORM,DELTA
      DOUBLE PRECISION DJ1,DJ2,DJ3,DJ4
      DOUBLE PRECISION XI1,XI2,XX,YY
      DOUBLE PRECISION OM,HBAS,HBASJ1,HBASJ2,HBASJ3,HSUMJ,AH
      DOUBLE PRECISION HBASI1,HBASI2,HBASI3,HSUMI,DENTH

C     In case that IPRECA=4, create the Laplace-part of the system
C     and include it into the system matrix each time this routine
C     is called. Otherwise don't create the Laplace part, just
C     add the nonlinearity.
C     This handling is simply realized by setting the factor NU in
C     front of the Laplace part to 0.

      IF (IPRECA.EQ.4) THEN
        DNY=NY
      ELSE
        DNY=0D0
      ENDIF

C     A similar handling holds for the case that the (full!) mass
C     matrix is to be included into the system matrix while the
C     matrix is to be rebuild completely. 
C     DCMASS represents the factor in front of the mass matrix.
C     Note that we cannot use DCMASS directly! In the Theta scheme,
C     we typically calculate 
C         [ DCMASS*M + THSTEP*(Laplace) + ... ]
C     so we have to weight everything except for the mass matrix!
C     We make a little trick here to realize that. We always weight
C     everything by THSTEP including the mass matrix - but divide
C     DCMASS by THSTEP before to compensate that!
C         THSTEP * [ CT0*M + (Laplace) + ... ]
C     with CT0 = DCMASS/THSTEP.
C
C     If only the nonlinear part is to be build (which is the normal
C     case, as rebuilding the comlete matrix including the full mass
C     matrix is normally not necessary), simply set CT0=0 which 
C     prevents calculating the mass matrix.

      IF ((IPRECA.EQ.4).AND.(IMASS.EQ.1)) THEN
        CT0=DCMASS/THSTEP
      ELSE
        CT0=0D0
      ENDIF

C     Initialize BDER for our element. We want the element to calculate
C     function values as well as first X/Y derivatives:

      DO I = 1,NNDER
        BDER(I)=.FALSE.
      END DO

      DO I=1,3
        BDER(I)=.TRUE.
      END DO

C     Ask the element about its type:

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      
C     Get the local number of degrees of freedom:
      
      IDFL=NDFL(IELTYP)
      
C     Initialize the cubature formula identifyer in the COMMON block
C     with our cubature formula we have to use here:

      ICUB=ICUBN
      CALL CB2Q(ICUB)
      
      IF (IER.NE.0) GOTO 99999

************************************************************************
*** Calculation of the matrix - storage technique 7 or 8
************************************************************************

C     In contrast to the parametric case, we don't initialize
C     the element here - it wuld be useless...
C
C     Calculate the maximum norm of the actual velocity field
C     U = A1*U1 + A2*U2 into DUMAX. 
C     Round up the norm to 1D-8 if it's too small...

      DUMAX=0D0
      IF (A2L.EQ.0D0) THEN
        DO IEQ=1,NMT
          DU1=U1L1(IEQ)
          DU2=U1L2(IEQ)
          DUNORM=SQRT(DU1**2+DU2**2)
          DUMAX=MAX(DUMAX,DUNORM)
        END DO
      ELSE       
        DO IEQ=1,NMT
          DU1=A1L*U1L1(IEQ)+A2L*U2L1(IEQ)
          DU2=A1L*U1L2(IEQ)+A2L*U2L2(IEQ)
          DUNORM=SQRT(DU1**2+DU2**2)
          DUMAX=MAX(DUMAX,DUNORM)
        END DO
      ENDIF       

      IF (DUMAX.LT.1D-8) DUMAX=1D-8

C *** Loop over all elements

      DO IEL=1,NEL
      
C       Calculate the local and global degrees of freedom on the
C       current element IEL:
      
        CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999

C       Determine local DELTA for streamline-diffusion
C       (cf. p. 131/132 (121) in Turek's CFD book).
C
C       For Stokes flow, we have the equation
C
C                -nu*div(u) + grad(p) = f
C       
C       not containing any convective part. Therefore we don't need
C       any stabilization technique. So in this case, switch of the
C       stabilization by setting DELTA to 0:

        IF (ISTOK.NE.1) THEN
          CALL DELTSD(U1L1,U1L2,U2L1,U2L2,A1L,A2L,IEL,DUMAX,DELTA,
     *                KVERT,KMID,DCORVG)
        ELSE
          DELTA=0D0
        END IF

C       Determine entry positions in matrix.
C
C       Here we build a small matrix DENTRY/KENTRY which
C       corresponds to the current element. We will assemble the
C       contributions of our current element into this matrix and
C       will add the result into the main matrix later.
C
C       To successfully integrate the contributions into the main
C       matrix, we compute in advance the positions in the main
C       matrix where we have to add the contribution to.
C       KENTRY(X,Y) corresponds to the index in the array A
C       of the entry in line Y, row KCOL(line start + X).

        DO JDOFE=1,IDFL
          ILD=KLDA(KDFG(JDOFE))
          KENTRY(JDOFE,JDOFE)=ILD
          DENTRY(JDOFE,JDOFE)=0D0
          JCOL0=ILD
          DO IDOFE=1,IDFL
            IF (IDOFE.NE.JDOFE) THEN
              IDFG=KDFG(IDOFE)
              DO JCOL=JCOL0,NA
                IF (KCOLA(JCOL).EQ.IDFG) GOTO 113
              END DO
113           JCOL0=JCOL+1
              KENTRY(JDOFE,IDOFE)=JCOL
              DENTRY(JDOFE,IDOFE)=0D0
            END IF
          END DO
        END DO

C       Calculate auxiliary Jacobian factors for the transformation 
C       onto the reference element. See QTRAFO.F for details...

        DO IVE = 1, NVE
          JP=KVERT(IVE,IEL)
          KVE(IVE)=JP
          DX(IVE)=DCORVG(1,JP)
          DY(IVE)=DCORVG(2,JP)
        END DO

        DJ1=0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2=0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3=0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4=0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))

C       Now that we know the coordinates of the corners of the
C       element, make the dummy call to ELE to precalculate
C       information if possible (may save arithmetic operations).
C       This is different to the parametric case, where it was
C       possible to do it only once outside of the loop over
C       the elements...

        CALL ELE(0D0,0D0,-2)
        IF (IER.LT.0) GOTO 99999

C       Loop over the cubature points in our element to calculate
C       its contribution to each of the degrees of freedom on our
C       element

        DO ICUBP = 1, NCUBP
        
C         Get the coordinates of the cubature point on the reference
C         element
        
          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)

C         Calculate the Jacobian of the bilinear mapping onto the
C         reference element and the weight OM of the cubature formula:

          DJAC(1,1)=0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
          DJAC(1,2)=0.5D0*DJ1+0.5D0*DJ2*XI1
          DJAC(2,1)=0.5D0*DJ4-0.5D0*DJ3*XI2
          DJAC(2,2)=0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
          DETJ=DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
          OM=DOMEGA(ICUBP)*DETJ

C         On difference to the parametric case, calculate the
C         coordinates of the cubature point on the actual element.
C         We'll need them...

          XX=0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *      +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
          YY=0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1
     *      +0.5D0*(DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2

C         Call the element to calculate the values in the current
C         cubature point on the reference element:

          CALL ELE(XX,YY,-3)
          IF (IER.LT.0) GOTO 99999

C         Now we have to assemble the "local" matrix DENTRY/KENTRY.
C         This assembling decomposes now into different parts,
C         depending on what has to me assembled.
C
C         We want to set up the nonlinear part of the matrix
C
C           n~_h (u_h, u_h, v_h) 
C
C         = n_h (u_h, u_h, v_h) + sum_T ( delta_T ( u_h*grad u_h, u_h*grad v_h)_T )
C           ^^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
C          standard nonlin. part                  stabilization
C
C         More precisely, as we want to assemble the matrix which is 
C         later multiplied with coefficient vectors, we have to insert
C         basis functions in the above terms instead of u_h and v_h.
C         Assuming the representation u_h=sum_j(u_j*Phi_j) and 
C         v_h=sum_i(u_i,Phi_i), the above term is evaluated in the
C         DOF's as:
C         
C           n_h (u_h, Phi_j, Phi_i) 
C         + sum_T ( delta_T ( u_h*grad Phi_j, u_h*grad Phi_i )_T )
C
C         In nonstationary simulations, the system matrix typically
C         contains a mass matrix to respect the time derivative.
C         The matrix has the form
C
C         [  dcmass*M*I  +  THWEIG * (-nu * Laplace(.))  ] + THWEIG * u grad(.)
C
C         so if DCMASS<>0, incorporate the (real, not lumped!) mass matrix
C         into the local matrix.

          IF (DCMASS.NE.0D0) THEN

C           Calculate the actual velocity in the current cubature point
C           into (DU1,DU2). If we only have a primary velocity field
C           (A2L=0), we can calculate that only by summing up the
C           velocities in U1Lx, otherwise we have to sum up
C           A1*U1Lx + A2*U2Lx.

            DU1=0D0
            DU2=0D0
            IF (A2L.EQ.0D0) THEN
              DO JDFL=1,IDFL
                HBAS=DBAS(KDFL(JDFL),1)
                IF (ABS(HBAS).GE.1D-8) THEN
                  JDFG=KDFG(JDFL)
                  DU1=DU1+U1L1(JDFG)*HBAS
                  DU2=DU2+U1L2(JDFG)*HBAS
                ENDIF
              END DO
            ELSE
              DO JDFL=1,IDFL
                HBAS=DBAS(KDFL(JDFL),1)
                IF (ABS(HBAS).GE.1D-8) THEN
                  JDFG=KDFG(JDFL)
                  DU1=DU1+(A1L*U1L1(JDFG)+A2L*U2L1(JDFG))*HBAS
                  DU2=DU2+(A1L*U1L2(JDFG)+A2L*U2L2(JDFG))*HBAS
                ENDIF
              END DO
            ENDIF

C           We take a more detailed look onto the last scalar product
C           of n~_h (u_h, u_h, v_h) what we want to calculate here.
C
C           The vector u_h=(DU1,DU2) contains both velocity components,
C           for the X as well as for the Y velocity. On the other hand
C           the system matrix we want to build here will be designed for 
C           one velocity component only! Therefore, Phi_i and Phi_j
C           are scalar functions, so grad(Phi_i), grad(Phi_j) are vectors
C           with two components. Therefore, the last scalar product is more 
C           in detail:
C
C               ( u_h*grad Phi_j, u_h*grad Phi_i )_T
C
C           =   ( < (DU1) , (grad(Phi_j)_1) > , < (DU1) , (grad(Phi_i)_1) > )_T
C                   (DU2) , (grad(Phi_j)_2)       (DU2) , (grad(Phi_i)_2)  
C
C           =   < (DU1) , (grad(Phi_j)_1) >  *  < (DU1) , (grad(Phi_j)_1) >
C                 (DU2) , (grad(Phi_j)_2)         (DU2) , (grad(Phi_j)_2)
C
C           =   HSUMJ * HSUMI
C
C           i.e. a product of two scalar values!
C
C           Summing up over all pairs of multiindices.
C
C           Outer loop over the DOF's j=1..IDFL on our current element, 
C           which corresponds to the basis functions Phi_j:

            DO JDOFE=1,IDFL
            
C             Fetch the contributions of the basis functions Phi_j for
C             function value and first derivatives for the current
C             DOF into HBASxy:
            
              JDOFEH=KDFL(JDOFE)
              HBASJ1=DBAS(JDOFEH,1)
              HBASJ2=DBAS(JDOFEH,2)
              HBASJ3=DBAS(JDOFEH,3)

C             Calculate 
C
C                 U * grad(Phi_j)  =  < grad(Phi_j), U >
C     
C               = ( grad(Phi_j)_1 , (DU1) )
C                 ( grad(Phi_j)_2   (DU2) )
              
              HSUMJ = HBASJ2*DU1 + HBASJ3*DU2

C             Inner loop over the DOF's i=1..IDFL, which corresponds to
C             the basis function Phi_i:

              DO IDOFE=1,IDFL
              
                IF (IDOFE.EQ.JDOFE) THEN
                
C                 Short version of the evaluation of the matrix
C                 contribution - see below for a more detailed
C                 description what is added together here!
                
                  AH = HSUMJ*(DELTA*HSUMJ+HBASJ1)
     *               + DNY*(HBASJ2**2+HBASJ3**2)
     *               + CT0*HBASJ1**2
     
                ELSE
                
C                 Fetch the contributions of the basis function Phi_i for
C                 function value and first derivatives for the current
C                 DOF into HBASIy:
                
                  IDOFEH=KDFL(IDOFE)
                  HBASI1=DBAS(IDOFEH,1)
                  HBASI2=DBAS(IDOFEH,2)
                  HBASI3=DBAS(IDOFEH,3)

C                 Calculate 
C
C                     U * grad(Phi_i)  =  < grad(Phi_i), U >
C     
C                   = ( grad(Phi_i)_1 , (DU1) )
C                     ( grad(Phi_i)_2   (DU2) )

                  HSUMI=HBASI2*DU1+HBASI3*DU2
     
C                 Finally calculate the controbution to the system
C                 matrix. Depending on the configuration of DNY,
C                 IPRECA, DCMASS,... this decomposes into three
C                 different parts:
C
C                 AH = n~_h(u_h,phi_j,phi_i)        | nonlinear part
C                    + DNY*(grad(phi_j,grad(phi_i)) | -nu*Laplace(u)
C                    + CT0*(phi_j*phi_i)            | Mass matrix
C
C                 The last two parts are probably not added to the
C                 matrix by setting DNY or CT0 to 0, respectively.
C
C                 For saving some numerical operations, we write:
C     
C                     HSUMI * (Delta * HSUMJ + HBASJ1)
C
C                 =   Delta * HSUMI * HSUMJ  
C                   + HSUMI * HBASJ1
C     
C                 =   Delta * ( U*grad(Phi_j), U*grad(Phi_i) )
C                   + (grad(Phi_i)*U,Phi_j)
C
C               <->   sum_T ( delta_T ( u_h*grad Phi_j, u_h*grad Phi_i) )_T
C                   + n_h (u_h, Phi_j, Phi_i)
                  
                  AH = HSUMI*(DELTA*HSUMJ+HBASJ1)
     *               + DNY*(HBASJ2*HBASI2+HBASJ3*HBASI3)
     *               + CT0*HBASJ1*HBASI1
     
                ENDIF ! (IDOFE.EQ.JDOFE)

C               Weighten the calculated value AH by the cubature
C               weight OM and add it to the local matrix. After the
C               loop over all DOF's is finished, each entry contains
C               the calculated integral.

                DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
                
              END DO ! IDOFE
              
            END DO ! JDOFE

          ELSE

C           Coefficient in front of the mass matrix is 0.
C           Subtract the mass matrix from the system matrix.
C
C           Outer loop over the DOF's j=1..IDFL corresponding
C           to Phi_j:

            DO JDOFE=1,IDFL
            
C             Fetch the function value in the current DOF:
            
              HBASJ1=DBAS(KDFL(JDOFE),1)

C             Inner loop over the DOF's i=1..IDFL corresponding
C             to Phi_i:

              DO IDOFE=1,IDFL
                
C               Fetch the function value of the other DOF of the
C               current element
                
                HBASI1=DBAS(KDFL(IDOFE),1)
                
C               Calculate the contribution for the entry. The factor
C               THSTEP is compensated later when the local matrix
C               is included into the global matrix and/or in the
C               modification of the RHS vector.
                
                AH=-1D0/THSTEP*HBASJ1*HBASI1
                
C               ...and incorporate it into the local matrix.
                
                DENTRY(JDOFE,IDOFE)=DENTRY(JDOFE,IDOFE)+OM*AH
                
              END DO ! IDOFE
              
            END DO ! JDOFE

          END IF ! (DCMASS.NE.0D0)

        END DO ! ICUBP
      
C       Now we have set up a "local" system matrix. We can either
C       include it into the real matrix or we can use it to simply
C       modify the RHS vector to create a defect vector (throwing
C       away the information about the matrix afterwards, which would
C       result in a matrix free modification of the RHS vector).
      
        DO JDOFE=1,IDFL
        
          DO IDOFE=1,IDFL
          
C           Get the entry from the local matrix and weight it according
C           to the current THETA in the Theta scheme, given by the
C           parameter. 
C           (Remark: For stationary simulations, THSTEP is typically
C            1D0 which includes the local matrix into the global one
C            directly).
          
            DENTH=THSTEP*DENTRY(JDOFE,IDOFE)

C           For IDEF=0/1, incorporate our "local" system matrix into 
C           the global matrix. The position of each entry DENTRY(X,Y) 
C           in the global matrix array A was saved in element KENTRY(X,Y)
C           before.

            IF (IDEF.LT.2) THEN
              IA   =KENTRY(JDOFE,IDOFE)
              A(IA)=A(IA)+DENTH
            ENDIF

C           For IDEF=1,2, build the defect vector
C               D = RHS - A*U
C           This is done matrix free, only with the help of the local 
C           matrix.
C           In this case, D=(D1,D2) is expected to be the RHS on
C           entry and will be updated to be the defect vector when
C           this routine is left.

            IF (IDEF.GT.0) THEN 
              IDFG=KDFG(IDOFE)
              JDFG=KDFG(JDOFE)
              D1(JDFG)= D1(JDFG)-DENTH*U1(IDFG)
              D2(JDFG)= D2(JDFG)-DENTH*U2(IDFG)
            ENDIF 

          END DO ! IDOFE
          
        END DO ! JDOFE

C       Matrix/defect vector updated for that element; proceed to the
C       next one...

      END DO ! IEL

99999 END

************************************************************************
* Calculation of a local DELTA
*
* This routine calculates a local DELTA=DELTA_T for a finite element
* T=IEL. This can be used by the streamline diffusion stabilization
* technique as a multiplier of the (local) bilinear form.
*
* In:
*   IEL    - Element where the Delta should be calculated
*
*   KVERT,
*   KMID,
*   DCORVG - Usual geometry information
*
*   U1L1,
*   U1L2   - array [1..NU] of double
*            Main velocity field. 
*   U2L1,
*   U2L2   - array [1..NU] of double
*            Secondary velocity field. 
*   A1L    - double; weighting factor for U1L1/U1L2
*   A2L    - double; weighting factor for U2L1/U2L2
*   UPSAM  - double; user defined parameter for configuring the 
*            streamline diffusion.
*            < 0: Simple calculation of Delta, using
*                 DELTA = |UPSAM| * h_T
*            > 0: usually UPSAM = 0.1 .. 2;
*                 Samarskji-like calculation of DELTA using:
*                 DELTA = UPSAM * h_t/||u||_T * 2*Re_T/(1+Re_T)
*   DUMAX  - Maximum norm of velocity in the domain:
*            DUMAXR = ||u||_Omega
*
* Out:
*   DELTA  - local Delta
*
* Used COMMON block variables:
*   IDFL   - Number of degrees of freedom on element IEL
*   KDFL   - array [1..IDFL] of integer
*            Array with global degrees of freedom, corresponding to
*            local degrees of freedom 1..IDFL on element IEL.
*   NY     - Coefficient NU in front of the Laplacian term 
*            of the Navier-Stokes equation
*               NU * Laplace(u) + u*grad(u) + ...
* 
* Remarks:
*
*   The effective velocity that is used for calculating the DELTA
*   is combined by a weighted mean of the two velocity fields U1,U2
*   by:
*                     Ux = A1*U1Lx + A2*U2Lx
*   The coefficients A1,A2 allow the caller to take influence on which
*   velocity field to weight more.
*
************************************************************************

      SUBROUTINE  DELTSD  (U1L1,U1L2,U2L1,U2L2,A1L,A2L,IEL,DUMAX,DELTA,
     *                     KVERT,KMID,DCORVG)

      IMPLICIT NONE
      
C include necessary COMMON blocks

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc' 
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cbasicelem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      
      INCLUDE 'casmbly.inc'
      INCLUDE 'cnsparfrac.inc'
      
      INCLUDE 'cinidat.inc'

C parameters
      
      DOUBLE PRECISION U1L1(*),U1L2(*),U2L1(*),U2L2(*)
      DOUBLE PRECISION A1L,A2L,DUMAX,DELTA
      INTEGER IEL
      INTEGER KVERT(4,*),KMID(4,*)
      DOUBLE PRECISION DCORVG(2,*)

C local variables
      
      DOUBLE PRECISION HLOCAL,DU1,DU2,UNORM,RELOC
      INTEGER IDOF

C     Loop through the local degrees of freedom on element IEL.
C     Sum up the velocities on these DOF's. This will result
C     in the vector (DU1,DU2) representing the (mean) X/Y-velocity
C     through element IEL.
C
C     For Ex30/Ex31 element, U1/U2 represent the mean velocity
C     along an egde/on the midpoint of each edge, so U1/U2 is
C     clearly an approximation to the velocity in element T.

      DU1=0D0
      DU2=0D0
      DO IDOF=1,IDFL
        DU1=DU1+(A1L*U1L1(KDFG(IDOF))+A2L*U2L1(KDFG(IDOF)))
        DU2=DU2+(A1L*U1L2(KDFG(IDOF))+A2L*U2L2(KDFG(IDOF)))
      END DO

C     Calculate the norm of that local velocity:

      UNORM = SQRT(DU1**2+DU2**2) / DBLE(IDFL)
      
C     Now we have:   UNORM = ||u||_T
C     and:           u_T = a1*u1_T + a2*u2_T
C
C     If the norm of the velocity is small, we choose DELTA = 0,
C     which results in central difference in the streamline diffusion
C     matrix assembling:

      IF (UNORM.LE.1D-8) THEN
        DELTA=0D0
      ELSE

C       u_T defines the "slope" of the velocity through
C       the element T. At next, calculate the local mesh width
C       HLOCAL = h = h_T on our element T=IEL:

         CALL HWAHL (HLOCAL,UNORM, DU1, DU2, IEL, KVERT,KMID,DCORVG)

C       Calculate DELTA... (cf. p. 131/132 (121) in Turek's CFD book)

        IF (UPSAM.LT.0D0) THEN

C         For UPSAM<0, we use simple calculation of Delta:        

          DELTA=ABS(UPSAM)*HLOCAL

        ELSE
        
C         For Delta >= 0, we use standard Samarskji-like calculation
C         of Delta. At first calculate the local Reynolds number
C         RELOC = Re_T = ||u||_T * h_T / NU
          
          RELOC=UNORM*HLOCAL/NY
          
C         and then the DELTA = UPSAM * h_t/||u|| * 2*Re_T/(1+Re_T)
          
          DELTA=UPSAM*HLOCAL/DUMAX*2D0*(RELOC/(1D0+RELOC))

        ENDIF ! (UPSAM.LT.0D0)
        
      ENDIF ! (UNORM.LE.1D-8)

      END

************************************************************************
* Determine the local mesh width for an element JEL of a triangulation.
* 
* In:
*   JEL    - Element where the local h should be calculated
*   UNORM  - norm ||u||_T = mean velocity through element T=JEL
*   XBETA1,
*   XBETA2 - mean velocity u_T = (xbeta1,xbeta2) through element T=JEL
*
*   KVERT,
*   KMID,
*   DCORVG - Usual geometry information
*
* Out:
*   HLOCAL - local mesh width
************************************************************************

      SUBROUTINE HWAHL (HLOCAL, UNORM,  XBETA1, 
     *                  XBETA2, JEL,KVERT,KMID,DCORVG)

      IMPLICIT NONE
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'cbasicelem.inc'
      INCLUDE 'celem.inc'

C parameters

      INTEGER JEL, KVERT(NNVE,*),KMID(NNVE,*)
      DOUBLE PRECISION DCORVG(2,*), HLOCAL, UNORM, XBETA1, XBETA2
      
C local variables
      
      DOUBLE PRECISION LAMBDA
      INTEGER NECK1,NECK2,NECK3,NECK4
      DOUBLE PRECISION X1,Y1,X2,Y2,X3,Y3,X4,Y4
      DOUBLE PRECISION ALPMAX, ALPHA


C     Fetch the numbers of the four corners of element JEL

      NECK1=KVERT(1,JEL)
      NECK2=KVERT(2,JEL)
      NECK3=KVERT(3,JEL)
      NECK4=KVERT(4,JEL)

C     Fetch the coordinates of these corners

      X1=DCORVG(1,NECK1)
      Y1=DCORVG(2,NECK1)
      X2=DCORVG(1,NECK2)
      Y2=DCORVG(2,NECK2)
      X3=DCORVG(1,NECK3)
      Y3=DCORVG(2,NECK3)
      X4=DCORVG(1,NECK4)
      Y4=DCORVG(2,NECK4)

C     Scale: (deactivated)
C
C      skal=max(xbeta1,xbeta2)
C      xbeta1=xbeta1
C      xbeta2=xbeta2

      ALPMAX=0D0

C     Loop through the four corners of element JEL and check
C     of a line with slope BETA=(xbeta1,xbeta2) starting in this
C     corner really intersects with one of the edges of the element.
C     Remark that we only have to check the two opposite edges
C     to the current corner!
C
C     -----------------------------------------------------------------
C     Check the first corner:

      CALL SCHNITT(X1,Y1,ALPHA,XBETA1,XBETA2,
     *            X3,Y3,LAMBDA,X2,Y2)
      ALPMAX=MAX(ALPHA,ALPMAX)

      CALL SCHNITT(X1,Y1,ALPHA,XBETA1,XBETA2,
     *            X3,Y3,LAMBDA,X4,Y4)
      ALPMAX=MAX(ALPHA,ALPMAX)
 
C     -----------------------------------------------------------------
C     The second one...
      
      CALL SCHNITT(X2,Y2,ALPHA,XBETA1,XBETA2,
     *            X4,Y4,LAMBDA,X1,Y1)
      ALPMAX=MAX(ALPHA,ALPMAX)

      CALL SCHNITT(X2,Y2,ALPHA,XBETA1,XBETA2,
     *            X4,Y4,LAMBDA,X3,Y3)
      ALPMAX=MAX(ALPHA,ALPMAX)

C     -----------------------------------------------------------------
C     The third one...

      CALL SCHNITT(X3,Y3,ALPHA,XBETA1,XBETA2,
     *            X1,Y1,LAMBDA,X2,Y2)
      ALPMAX=MAX(ALPHA,ALPMAX)

      CALL SCHNITT(X3,Y3,ALPHA,XBETA1,XBETA2,
     *            X1,Y1,LAMBDA,X4,Y4)
      ALPMAX=MAX(ALPHA,ALPMAX)

C     -----------------------------------------------------------------
C     And the fourth=last one...
      
      CALL SCHNITT(X4,Y4,ALPHA,XBETA1,XBETA2,
     *            X2,Y2,LAMBDA,X1,Y1)
      ALPMAX=MAX(ALPHA,ALPMAX)

      CALL SCHNITT(X4,Y4,ALPHA,XBETA1,XBETA2,
     *            X2,Y2,LAMBDA,X3,Y3)
      ALPMAX=MAX(ALPHA,ALPMAX)

C     -----------------------------------------------------------------
C     finally determine the local h=h_T

      HLOCAL=ALPMAX*4D0*unorm

      END


*******************************************************************
* Intersect two lines in R^2
*******************************************************************

      SUBROUTINE SCHNITT (XO,YO,ALPHA,BETA1,BETA2,
     *                    XA,YA,LAMBDA,XB,YB)
C
C   Schnitt von zwei Geraden im R^2
C
      IMPLICIT NONE
      
C parameters

      DOUBLE PRECISION XO, YO, BETA1, BETA2, ALPHA, XA, YA, XB, YB
      DOUBLE PRECISION LAMBDA
      
C local variables

      DOUBLE PRECISION SKAL
      LOGICAL BFLAG

      skal=beta2*(xb-xa)-beta1*(yb-ya)
C      
      if (skal.eq.0D0) then
C  
C        beta und der Richtungsvektor sind parallel
C
         alpha=0D0
c         write(*,*) 'eins'
         bflag=.false.
      else  
           lambda=(beta1*(ya-yo)-beta2*(xa-xo))/skal
           bflag=.true.      
      endif
C
C     Ist der Schnittpunkt innerhalb des Elements?
C
      if (bflag) then
      if ((lambda.ge.-1e-1).and.(lambda.le.1.11e0)) then
         if (beta1.ne.0D0) then
            alpha=((xa-xo)+lambda*(xb-xa))/beta1
         else
            if (beta2.ne.0D0) then
               alpha=((ya-yo)+lambda*(yb-ya))/beta2
            else
               alpha=0D0
           endif
         endif
      else
c         write(*,*) 'drei'
         alpha=0D0
      endif
      endif
      END

