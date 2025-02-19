************************************************************************
* NFBDYC - Number of fictitious boundary components
*
* User provided function. Has to return the number of fictitious
* boundary components in the given geometry.
*
* Remark: If 0 is returned, fictitious boundary handling is disabled !!!
************************************************************************

      INTEGER FUNCTION NFBDYC ()
      
      IMPLICIT NONE
      
      INCLUDE 'cinigeometry.inc'
      
C Normally only one f.b.c.
C If we have a composed geometry, we have more perhaps. This can be seen
C in the IGEOCN variable.
      
C The standard implementation sets NFBDYC=0, so deactivating fictitious
C boundary handling:

      NFBDYC = 0
      
      END
      
************************************************************************
* ISFBDY - Test for being in a fictitious boundary object
*
* This is a user provided routine. It has to check whether the given
* point (x,y) is in an object surrounded by a fictitious boundary.
* The return value must be = 0 if the point is a normal point and
* != 0 if the point is in an object with a fictitious boundary.
*
* in:
*  X    - X-coordinate of the point to test
*  Y    - Y-coordinate of the point to test
*  IFCB - =0: test, if (X,Y) is in any fictitious boundary component;
*             the return value is either 0 or the number of any
*             fictitious boundary component containing this point
*         >0: test, if (X,Y) is contained especially in boundary
*             component IFBC. The return value is either 0 or
*             +-IFBC, depending on whether IFBC contains (X,Y) or not.
*
* result:
*   0 - If (x,y) is not in a fictitious boundary object
*  >0 - If (x,y) is in a rigid fictitious boundary object.
*  <0 - If (x,y) if in a virtual fictitious boundary object
*
*  The absolute value defines the number of a/the fictitious boundary
*  component containing (X,Y). If IFCB=0, this is the number
*  of any of the fictitious boundary components containing (X,Y).
*  If IFCB>0, the return value is 0 or +-IFCB, depending on if
*  (X,Y) is contained especially in this component or not.
*
*  A positive number indicates a rigid fictitious boundary object (wall)
*  which must be treated as dirichlet-boundary in the framework.
*  A negative number indicates a "virtual" fictitious boundary
*  object that is treated as Neumann boundary. In such points the
*  framework does not perform any changes to the matrix/RHS/...
*  to implement any type of boundary conditions. This can be used
*  e.g. for defining a region where to refine the grid without
*  changing the type of geometry in that region.
*
* This routine has to consider the number of fixed boundary components,
* too. As all fixed boundary conponents take the numbers 1..NBCT,
* all moving boundary components take numbers NBCT+1..NBCT+NFBDYC.
* So the lowest number != 0 returned by this function has to be
* +-(NBCT+1) !
************************************************************************

      INTEGER FUNCTION ISFBDY (X,Y,IFBC)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      INCLUDE 'cinigeometry.inc'
      
C COMMON-blocks with given parameters of the fictitious boundary
C component(s)

      INCLUDE 'cinidat.inc'
      
C parameters

      DOUBLE PRECISION X,Y
      INTEGER IFBC

C externals

      INTEGER NFBDYC
      EXTERNAL NFBDYC
      
      INTEGER GINCMP,GINWRP
      EXTERNAL GINCMP,GINWRP

C local variables

      DOUBLE PRECISION XM1,YM1,DIST

C quick parameter check

      IF ((IFBC.LT.0).OR.(IFBC.GT.NFBDYC())) THEN
        WRITE (*,'(A,I,A)') 'ERROR in call to ISFBDY()! IFBC=',IFBC,
     *                      ' Program halted!'
        STOP
      END IF

C standard return value

      ISFBDY = 0

C In this implementation we only have one boundary component, so we
C don't have to take care of IFBC for now...

      IF (ICIRTP.EQ.0) THEN

C ----------------------------------------------------------------------
C fictitious boundary component: circle at DCXPOS/DCYPOS, radius DCRAD
C ----------------------------------------------------------------------

        XM1=DCXPOS
        YM1=DCYPOS

        DIST=DSQRT((X-XM1)**2+(Y-YM1)**2)

        IF (DIST.LE.DCRAD) THEN
          ISFBDY = NBCT+1
        END IF

      ELSE IF (ICIRTP.EQ.1) THEN

C ----------------------------------------------------------------------
C fictitious boundary component: square with midpoint at DCXPOS/DCYPOS, 
C width/height=2xDCRAD
C ----------------------------------------------------------------------

        XM1=DCXPOS
        YM1=DCYPOS

        DIST=MAX(ABS(X-XM1),ABS(Y-YM1))

        IF (DIST.LE.DCRAD) THEN
          ISFBDY = NBCT+1
        END IF
        
      ELSE IF (ICIRTP.GE.1) THEN

C User-defined complex geometry.
C nothing implemented here.

      END IF

C Every other value for ICIRTP disables the fictitious boundary 
C management, i.e. there is no fict. boundary component at all.

      END
      
************************************************************************
* FBDVOL - Get the (analytical) volume of the fictitious
*          boundary component
*
* This routine computes the reference volume of the fictitious boundary
* component IFBC. This is used as the reference volume in comparisons.
*
* In:
*  IFBC   - Number of fictitious boundary component to calculate
*           the volume for.
*           0 = calculate summed volume of all fictitious boundary
*               components in the domain
*           NBCT+1..NBCT+NFBDYC = calculate volume only for
*               fictitious boundary component IFBC
*
* Return: calculated volume or -1D0, if the calculation is not supported
************************************************************************

      DOUBLE PRECISION FUNCTION FBDVOL (IFBC)

      IMPLICIT NONE
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      INCLUDE 'cinigeometry.inc'
      
C COMMON-blocks with given parameters of the fictitious boundary
C component(s)

      INCLUDE 'cinidat.inc'
      
C parameters

      INTEGER IFBC

      FBDVOL=-1D0

      IF ((IFBC.EQ.0).OR.
     *    ((FBDVOL.GT.0).AND.(IFBC.EQ.1))) THEN

        IF (ICIRTP.EQ.0) THEN
C circle
          FBDVOL = DCRAD**2*PI
          
        END IF
        
      END IF
      
      END
      
      
************************************************************************
* FBDGMV - Write analytical fictitious boundaries to GMV file
*
* This routine is called by the framework in the postprocessing
* part. It should write out an analytical description of all
* fictitious boundary parts to the file unit FUNIT in syntax
* of GMV (as far as this is possible).
*
* The routine is allowed to use material number MMAT..MMAT+NFBDYC()
* for the fictitious bondary components.
*
* The caller of this routine is responsible for the administration
* of the sections in the GMV file. This routine only has to write out
* the content for the POLYGON's section to approximate the boundary
* components. The header/footer of this section is written to the
* file by the caller.
*
* In:
*  MUNIT  - File unit where to write the GMV output to.
*  MMAT   - The first material number to use.
************************************************************************

      SUBROUTINE FBDGMV (MUNIT,MMAT)

      IMPLICIT NONE
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      INCLUDE 'cmem.inc'

C COMMON-blocks with given parameters of the fictitious boundary
C component(s)

      INCLUDE 'cinidat.inc'
      INCLUDE 'cinigeometry.inc'
      
C parameters

      INTEGER MUNIT,MMAT
      
C externals

      INTEGER NFBDYC
      EXTERNAL NFBDYC
      
      INTEGER GCRELL
      EXTERNAL GCRELL
      
C local variables

      INTEGER I,NL
      DOUBLE PRECISION DEG
      
      IF (NFBDYC().EQ.0) RETURN

C One circle:

      IF (ICIRTP.EQ.0) THEN
        
C As we only have one boundary component, this routine is only 
C called once...
      
      
C Emulate the circle by a number of lines...

C Number of circle segments is dependent of radius

        NL = MIN(10000,MAX(1,NINT(DCRAD*2*PI*1000)))

C Material number MMAT=moving boundary, NL+1 nodes

        WRITE(MUNIT,*) MMAT,NL+1

C X-coordinates

        DO I=0,NL
          DEG = DBLE(I)*2D0*PI/DBLE(NL)
          WRITE(MUNIT,*) DCXPOS+DCRAD*SIN(DEG)
        END DO

C Y-coordinates

        DO I=0,NL
          DEG = DBLE(I)*2D0*PI/DBLE(NL)
          WRITE(MUNIT,*) DCYPOS+DCRAD*COS(DEG)
        END DO

C Z-coordinates

        DO I=0,NL
          WRITE(MUNIT,*) 0D0
        END DO
      
      ELSE IF (ICIRTP.GE.1) THEN

C User-defined geometry

      END IF
      
      END

************************************************************************
* Prescribe Dirichlet boundary conditions on fictitious boundary
*
* This routine is called by the framework when implementing Dirichlet 
* boundary conditions in a point of the fictitious boundary domain.
* It can be used to prescribe e.g. tangential velocity initiated
* by rotation.
* ITYP describes the information that has to be returned in the
* return value.
*
* In:
*  ITYP   - 1=X-velocity, 2=Y-velocity
*  X,Y    - the point
* 
* Out:
*  Return value = the desired information
************************************************************************

      DOUBLE PRECISION FUNCTION FBINDT (ITYP,X,Y)
      
      IMPLICIT NONE
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      
      INCLUDE 'cinigeometry.inc'
      
C parameters

      INTEGER ITYP
      DOUBLE PRECISION X,Y
      
C local variables

      DOUBLE PRECISION XM1,YM1,XR1,YR1,XN1,YN1,SPD
      
C Standard return value: no velocity.
      
      FBINDT = 0D0
      
      GOTO 99999
      
C In case of a circle we might prescribe the rotation

      IF (ICIRTP.EQ.0) THEN

C fictitious boundary component: circle at DCXPOS/DCYPOS, radius DCRAD

C Midpoint:

        XM1 = DCXPOS
        YM1 = DCYPOS
        
C difference vector

        XR1 = X-XM1
        YR1 = Y-YM1
        
C create the normal vector to that; this is at the same time the tangential
C vector on the circle

        XN1 = YR1
        YN1 = -XR1
        
C normalise it, such that points on the boundary have velocity=STD

        SPD = 0.05D0

        XN1 = SPD * XN1/DCRAD**2
        YN1 = SPD * YN1/DCRAD**2
        
        IF (ITYP.EQ.1) FBINDT=XN1
        IF (ITYP.EQ.2) FBINDT=YN1

      END IF
      
99999 END
      
      
************************************************************************
* FBDNML - Get normal vector of fictitious boundary
*
* This routine is called by the framework to compute analytical normal
* vectors of the fictitious boundary. To a given point (X,Y) -
* which may be inside the fictitious boundary or not, depending
* on the approximation - this routine has to calculate a normalised
* vector (XN,YN) perpendicular to the fictitious boundary domain IFBC.
*
* In:
*  X,Y    - Point where to calculate the normal vector
*  IFBC   - Number of fictitious boundary component to calculate
*           the volume for.
*           0 = automatically determine a fictitious boundary
*               component which boundary is taken for the calculation
*               of the normal vector (i.e. by analytically analysing
*               the distance)
*           NBCT+1..NBCT+NFBDYC = Number of fictitious boundary
*               component which is used for the calculation of the
*               normal vector.
*
* Out:
*  XN,YN  - normal vector of the fictitious boundary component IFBC
*           in the point (X,Y); must be normalised to length=1D0 !
* 
* If calculation of (XN,YN) is possible, this routine should set
* (XN,YN)=(0D0,0D0) to indicate this to the framework.
************************************************************************

      SUBROUTINE FBDNML (X,Y,IFBC,XN,YN)

      IMPLICIT NONE
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'

C COMMON-blocks with given parameters of the fictitious boundary
C component(s)

      INCLUDE 'cinidat.inc'
      INCLUDE 'cinigeometry.inc'
      
C parameters

      DOUBLE PRECISION X,Y,XN,YN,A
      INTEGER IFBC

      XN = 0D0
      YN = 0D0

C For the moment we only have one boundary component, so IFBC is unused.

C Circle case:

      IF (ICIRTP.EQ.0) THEN
C Calculate the normal by the taking the difference vector
C to the midpoint

        XN = X-DCXPOS
        YN = Y-DCYPOS
      END IF

C normalise the vector

      IF ((XN.NE.0D0).OR.(YN.NE.0D0)) THEN
        A = 1D0/DSQRT(XN*XN+YN*YN)
        XN=A*XN
        YN=A*YN
      END IF
      
      END
      
************************************************************************
* FBDDST - Get distance of point to fictitious boundary interface  
*
* This routine is called by the framework to compute the minimum 
* distance of a point (X,Y) to the interface of a fictitious
* boundary component. The parameter IFBC determines the number
* of the boundary component to determine the distance to.
* When called with IFBC=0, the routine should calculate the minimum
* distance to any of the fictitious boundary interfaces.
* When called with IFBC>0 (i.e. IFBC of NBCT+1..NBCT+NFBDYC), 
* the routine should calculate the distance to that specific
* fictitious boundary component.
*
* The distance is returned in the parameter DIST. For points inside 
* of a fictitious boundary component, this value is <0 whereas 
* points outside the fictitious boundary are described by a value >0.
*
* If the distance was successfully calculated, the routine must
* return a value > 0. If it's not possible to determine the 
* distance (either if it's not implemented or not possible 
* due to complexity) the routine must return a value <= 0 to indicate
* this. The framework will then take care of the fact that determining
* the distance to this fictitious boundary component is not
* possible.
*
* In:
*  X,Y    - Point to determine the distance to the fictitious
*           boundary interface
*  IFBC   - Number of fictitious boundary component to calculate
*           the distance to.
*           0 = Compute the minimum distance to all fictitious
*               boundaty interfaces
*           NBCT+1..NBCT+NFBDYC = Number of fictitious boundary
*               component where to calculate the distance to.
* Out:
*  DIST   - Distance of (X,Y) to the interface of a/the fictitious
*           boundary
* 
* Return:
*  >  0   - If the distance was successfully computed
*  <= 0   - If the computation is not possible
* 
************************************************************************

      INTEGER FUNCTION FBDDST (X,Y,IFBC, DIST)

      IMPLICIT NONE

      INCLUDE 'cinigeometry.inc'
      
      DOUBLE PRECISION X,Y,DIST
      INTEGER IFBC
      
      FBDDST = 0

C We are only able to handle the circle case here:

      IF (ICIRTP.EQ.0) THEN
        DIST = DSQRT((X-DCXPOS)**2+(Y-DCYPOS)**2) - DCRAD
        FBDDST = 1
      END IF

      END

**********************************************************************
* Cubature point projection method, analytical version
*
* This subroutine must project cubature points on the reference 
* interval [-1,1] to real cubature point coordinates on a given
* quadrilateral.
*
* The cubature points on the reference interval are given
* in the COMMON-block variable DXI(1..NCUBP,1), with /CUB/.NCUBP the
* number of cubature points. The points (X1,Y1) and (X2,Y2) are
* reconstructed intersection points of the fictitious boundary
* interface with the edges of the quadrilateral IEL. DCOORD
* denote the corner vertices of this element in counterclockwise
* order. If IEL=0 there is no element belonging to the coordinates
* DCOORD, but still DCOORD describes the four corner points of
* a quadrilateral.
*
* The routine must translates the parameter values in DXI(.,.) into
* real coordinates on the interface of the fictitious boundary
* interface, which intersects with the edges of the element in
* the coordinates (X1,Y1) and (X2,Y2). The coordinates of these 
* points have to be stored in 
*              (DCUBP(1,1..NCUBP),DCUBP(2,1..NCUBP)).
*
* DJAC will receive the Jacobian determinant of the mapping between
* the reference interval [-1,1] to the interface segment which
* intersects with the edges in the coordates (X1,Y1) / (X2,Y2).
* This coordinate pair is ordered in that way that the vector
* (X2-X1,Y2-Y1) denotes the tangential and (-(Y2-Y1),X2-X1)
* the normal vector of the fictitious boundary component pointing
* "to the outside".
*
* If the mapping between the interval [-1,1] to the
* interface segment does not have a constant Jacobian determinant,
* the value of DETJ has to be set to 1D0. Furthermore the
* values of the norm of the Jacobian determinant have to be
* included into the weights of the cubature formula, which are
* stored in the DOMEGA-variables of the COMMON block /CUB/ !
*
* Depending on the type of the fictitious boundary component this
* routine can perform an exact reconstruction of the quadrature
* points on the fictitious boundary interface. In the case that
* these points are too complicated to be reconstructed analytically,
* this routine can use the subroutine CBLPRQ, which creates a linear
* reconstruction of the interface to distribute the cubature points
* appropriately.
*
* In:
*  DCOORD - array [1..2,1..4] of double
*           Coordinates of the four corners of the real quadrilateral.
*           DCOORD(1,.) saves the X-, DCOORD(2,.) the Y-coordinates.
*  IEL    - The number of current element or 0, if there is no current
*           element belonging to the coordinates in DCOORD.
*  (X1,Y1),
*  (X2,Y2) - Start- and endpoint of the line segment intersecting
*            the edges of the element. These coordinates are given
*            in "real" coordinates on the "real" element, not on the 
*            reference element.
*  IFBC   - Number of fictitious boundary component currently being
*           treated
*           0 = no explicit handling of a special component.
*               Automatically determine the best approximation of
*               any fictitious boundary component that intersects
*               with the element
*           NBCT+1..NBCT+NFBDYC = Number of fictitious boundary
*               component whose interface is being meant. 
* 
* Out:
*  DJAC    - Norm of the Jacobi determinant of the mapping between 
*            [-1,1] and the interface segment in real coordinates, 
*            or 1D0 if the Jacobian determinant determinant is not
*            constant.
*  DCUBP   - array [1..2,1..NCUBP] of double
*            Coordinates of cubature points in real coordinates
*            on the quadrilateral given by DCOORD.
*
* If the norm of the Jacobian determinant os not constant, the
* integration weights in /CUB/.DOMEGA have to be modified properly.
**********************************************************************
            
      SUBROUTINE FBCBPR (DCOORD,IEL,IFBC,X1,Y1,X2,Y2,DCUBP,DJAC)
      
      IMPLICIT NONE
      
      INCLUDE 'ccub.inc'
      INCLUDE 'cinigeometry.inc'
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

C parameters
      
      INTEGER IEL,IFBC
      DOUBLE PRECISION X1,Y1,X2,Y2,DJAC,P
      DOUBLE PRECISION DCOORD(2,4),DCUBP(2,NNCUBP)
      
C local variables
      
      DOUBLE PRECISION A1,A2,L1,L2,V1,V2,W1,W2,CURANG
      INTEGER I
      
C For now we only have the case of exactly one fictitious bondary
C component, so we don't care about IFBC.
      
C In case of a circle we know the interface:

      IF (ICIRTP.EQ.0) THEN

C Ok, quadrilateral IEL intersects with our circle in two points.
C We must know the angle of these two points. Therefore we transform
C them back onto a "reference circle" and calculate the angles
C with trigonometric functions:

        L1 = 1D0/DSQRT((X1-DCXPOS)**2+(Y1-DCYPOS)**2)
        V1 = L1*(X1-DCXPOS)
        W1 = L1*(Y1-DCYPOS)

        L2 = 1D0/DSQRT((X2-DCXPOS)**2+(Y2-DCYPOS)**2)
        V2 = L2*(X2-DCXPOS)
        W2 = L2*(Y2-DCYPOS)
        
C Use the ARCCOS in the upper or lower half circle
        
        A1 = DACOS(V1)
        IF (W1.LT.0D0) A1 = 2*PI - A1 

        A2 = DACOS(V2)
        IF (W2.LT.0D0) A2 = 2*PI - A2 
        
C Make sure angle A2 is "before" A1 in the orientation
        
        IF (A2.GT.A1) A2 = A2 - 2*PI;
        
C Ok, we parametrize for the arc length. [-1,1] has to be mapped
C to the arc between the angles A1 and A2. Loop through the
C cubature points and map them to DCUBP. 

        DO I=1,NCUBP
        
C Calculate the angle of the cubature point
        
          P = 0.5D0+(0.5D0*DXI(I,1))
          CURANG = (1D0-P)*A1 + P*A2
          
C Use that to calculate the cubature point itself

          DCUBP (1,I) = DCXPOS + DCRAD*COS(CURANG)
          DCUBP (2,I) = DCYPOS + DCRAD*SIN(CURANG)

        END DO
        
C Fortran uses (as usual) the arc length parametrization 0..2*PI
C for trigonometrical functions. Therefore the difference between
C A1 and A2 is the arc length of our circle segment - up to
C a constant factor, depending on the radius. This is then the
C Jacobian "determinant" of the mapping. We have to divide it
C by 2 because our reference interval is [-1,1], not [0,1].

        DJAC = 0.5D0*DCRAD*ABS(A2-A1)
        
      ELSE
      
C General handling: linear reconstruction of the interface and
C distribution of the cubature points

        CALL CBLPRQ (DCOORD,IEL,X1,Y1,X2,Y2,DCUBP,DJAC)
        
      END IF
      
      END
      
***********************************************************************
* Initialize geometry configuration 
*
* This initializes the COMMON-block variables for the fictitious
* boundary geometry - usually a circle.
***********************************************************************

      SUBROUTINE RDGEOM ()
      
      IMPLICIT NONE
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      
      INCLUDE 'cinigeometry.inc'
      
C     Set up a circle

      ICIRTP = 0
      DCXPOS = 0.2
      DCYPOS = 0.2
      DCRADX = 0.05
      DCRADY = 0.05
      DCROT  = 0
      
C initialise sin/cos(rot)

      DCRSIN = SIN(DCROT*PI/180D0)
      DCRCOS = COS(DCROT*PI/180D0)

      END
      
***********************************************************************
* Clean up geometry configuration 
*
* Cleans up all information that was initialized in RDGEOM.
* This can e.g. be used to dispose arrays. In the standard
* implementation there's nothing to do because the circle is a
* simple direct object that does not need and additional information
* on the heap...
***********************************************************************

      SUBROUTINE DONGEO ()

      IMPLICIT NONE

      END
      