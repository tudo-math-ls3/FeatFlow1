************************************************************************
* Specify input data for inflow/outflow/...
*
* User specified routine. Prescribed data for files coeff.f and bndry.f.
* With this routine the user can specify input data for the flow
* simulation like inflow, outflow, ...
* The framework repeatly calls this routine with different
* ITYP and IBLOC parameters. Depending on these parameters this
* routine has to return the desired infomation about the flow
* the framework needs.
*
* In:
*  ITYP   - Type of information the framework needs.
*           = 1: velocity dirichlet value
*           = 2: velocity x-derivative
*           = 3: velocity y-derivative
*           = 4: exact pressure
*           = 5: Right hand side for momentum equation
*           = 6: Right hand side for continuity equation
*           = 7: Mean pressure value
*  IBLOC  - Current matrix block; corresponds to the solution component
*           that is currently being assembled by the framework. This
*           must be used by this routine to determine the "direction"
*           of the desired information
*           = 1: matrix block for x-velocity
*           = 2: matrix block for y-velocity
*  X,Y    - Coordinates of the point where the framework needs
*           the information
*  TIMENS - Current time in instationary Navier-Stokes calculation
*  RE     - Reynolds number; from the DAT file
*
*
* Out:
*  Return value = desired information
************************************************************************

      DOUBLE PRECISION FUNCTION FDATIN(ITYP,IBLOC,X,Y,TIMENS,RE)

      IMPLICIT NONE
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgpar.inc'
      
C parameters

      DOUBLE PRECISION X,Y,TIMENS,RE
      INTEGER ITYP, IBLOC
      INTEGER INODE
      
C       external routines
      
      INTEGER ISFBDY,ISFBDM,NFBDYC
      EXTERNAL ISFBDY,ISFBDM,NFBDYC
      
C       local variables

      INTEGER INPR
      DOUBLE PRECISION DPAR
      
C     The following function implements a parabolic inflow profile with
C     maximum velocity = UMAX. LEN denotes the length of the profile
C     (e.g. length of an edge, depending on the configuration).
C     The parameter value t in [0,LEN] denotes the position 
C     where the profile should be evaluated.

      DOUBLE PRECISION UMAX,T,LEN,PPROF
      PPROF (T,LEN,UMAX) = 4*UMAX*T*(LEN-T)/(LEN*LEN)
      
C -------------------------------------------------
C set the standard result of this function to zero:
C -------------------------------------------------

      FDATIN=0D0

C --------------------------------------------------------------------
C Now analyse the input parameters and return the corresponding result
C for the current point in the geometry. This has to be specified
C by the user depending in the current geometry and simulation!
C
C In short words: Find out, which DIRICHLET information FEATFLOW 
C wants to know and return it!
C
C The type of information is specified by ITYP, the current point
C by X/Y, the current boundary component on IBLOC, the current time
C of the Navier Stokes simulation in TIMENS,...
C --------------------------------------------------------------------

C=======================================================================
C *** Case 1: Velocity boundary values and/or exact solution
C=======================================================================
C
      IF (ITYP.EQ.1) THEN

C       benchmark-like inflow

        IF (IBLOC.EQ.1) THEN

C         parabolic inflow profile from the left - horizontal channel
C         (bench1) fluid speed = 0.3

C          IF (X.EQ.0.0D0)     FDATIN= 4D0*0.3D0/0.1681D0*Y*(0.41D0-Y)

          IF (X.EQ.0.0D0)     FDATIN = PPROF(Y,0.41D0,0.3D0) 

        END IF
          
        IF (IBLOC.EQ.2) THEN
        END IF
      
C       Set the velocity on all boundary nodes of all moving boundary
C       objects to zero:

        IF (NFBDYC().GT.0) THEN
          IF (ISFBDY (X,Y,0).GT.0) THEN
            FDATIN=0D0
          END IF
        END IF

      END IF

C=======================================================================
C *** Case 2: Velocity x-derivative of exact solution
C=======================================================================

      IF (ITYP.EQ.2) THEN

       IF (IBLOC.EQ.1) THEN
        IF (X.EQ.0.0D0) FDATIN=0D0
       ENDIF

       IF (IBLOC.EQ.2) THEN
        IF (X.EQ.0.0D0) FDATIN=0D0
       ENDIF

      ENDIF

C=======================================================================
C *** Case 3: Velocity y-derivative of exact solution
C=======================================================================

      IF (ITYP.EQ.3) THEN

        IF (IBLOC.EQ.1) THEN
          FDATIN=0D0
        ENDIF

        IF (IBLOC.EQ.2) THEN
          FDATIN=0D0
        ENDIF

      ENDIF
C
C
C=======================================================================
C *** Case 4: Exact pressure solution
C=======================================================================

      IF (ITYP.EQ.4) THEN

       FDATIN=0D0

      ENDIF


C=======================================================================
C *** Case 5: Right hand side for momentum equation
C=======================================================================

      IF (ITYP.EQ.5) THEN

        IF (IBLOC.EQ.1) THEN
          FDATIN=0D0
        ENDIF

        IF (IBLOC.EQ.2) THEN
          FDATIN=0D0
        ENDIF

      ENDIF

C=======================================================================
C *** Case 6: Right hand side for continuity equation
C=======================================================================
C
      IF (ITYP.EQ.6) THEN

        FDATIN=0D0

      ENDIF

C=======================================================================
C *** Case 7: Mean pressure values
C=======================================================================

      IF (ITYP.EQ.7) THEN
        DPAR=X
        INPR=IBLOC

        IF ((DPAR.GT.1D0).AND.(DPAR.LT.2D0).AND.(INPR.EQ.1)) THEN
          FDATIN=0D0
        ENDIF

        IF ((DPAR.GT.3D0).AND.(DPAR.LT.4D0).AND.(INPR.EQ.1)) THEN
          FDATIN=0D0
        ENDIF

      ENDIF

99999 END

************************************************************************
* Specify Neumann boundary parts
*
* User specified routine. Prescribed data for files coeff.f and bndry.f.
* With this routine the user can specify boundary segments as
* Neumann boundary.
* The framework repeatly calls this routine with different
* parameters. Depending on these parameters this
* routine has to return the desired infomation about the flow
* the framework needs.
*
* In:
*  INPART - Type of information the framework needs.
*           = 0: Number of Neumann boundary segments in the simulation
*           > 0: Information about boundary segment INPART
*  TIMENS - Current time in instationary Navier-Stokes calculation
*
* Out:
*  If INPART = 0:
*    INPART - Number of boundary segments
*
*  If INPART > 0:
*    INPRN  - Number of the boundary component, DPAR1/DPAR2 refer to.
*    DPAR1,DPAR2
*           - Parameter values of the boundary segments that should be
*             treated as Neumann boudary
************************************************************************
      SUBROUTINE NEUDAT(INPART,INPRN,DPARN1,DPARN2,TIMENS)

      IMPLICIT NONE
      
C     parameters

      INTEGER INPART, INPRN
      DOUBLE PRECISION DPARN1, DPARN2, TIMENS
      
C --------------------------------------------------------------------
C Now analyse the input parameters and return the corresponding result
C for the current point in the geometry. This has to be specified
C by the user depending in the current geometry and simulation!
C
C In short words: Find out, which NEUMANN information FEATFLOW 
C wants to know and return it!
C --------------------------------------------------------------------

C=======================================================================
C *** Case 0: Specify number of Neumann-boundary parts
C=======================================================================

      IF (INPART.EQ.0) THEN

C       Benchmark-like configuration; one Neumann part on the right side

        INPART=1      
          
C=======================================================================
C *** Case <>0: Specify Neumann-boundary parts
C=======================================================================

      ELSE IF (INPART.GT.0) THEN

C       Neumann boundary on outflow
C
C       INPRN specifies the number of the boundary component which is
C       specified by the setting of DPRN1/DPRN2.
C       DPRN1/DPRN2 specify the minimum/maximum parameter value
C       on the boundary that determine the Neumann boundary component.
C
C       Benchmark-like configuration: Outflow on the right between
C       boundary parameter value 1.0 and 2.0.

        IF (INPART.EQ.1) THEN
          INPRN =1
          DPARN1=1D0
          DPARN2=2D0
        END IF
          
      ENDIF

99999 END

************************************************************************
* Data for Point-output (for fpost and bdpres and bdforc)
*
* This routine allowes to specify vertices whose velocity/pressure
* values are measured over time and written into external files for
* postprocessing by the user.
* Furthermore this routine allowes to specify coefficients in the
* integral for the drag- and lift-calculation, to customize that to the
* corresponding geometry.
*
* In:
*   TIMENS - Current simulation time
*   DNU    - Current viscosity coefficient NU=1/RE
*
* Out:
*   Modifies the COMMON block variables in /NSPTS/.
************************************************************************

      SUBROUTINE PTSDAT(TIMENS,DNU)
      
      IMPLICIT NONE
      
      INCLUDE 'cnspts.inc'
      INCLUDE 'cinigeometry.inc'

C externals

      INTEGER NFBDYC
      EXTERNAL NFBDYC

C parameters

      DOUBLE PRECISION TIMENS, DNU

C local variables

      DOUBLE PRECISION RHO, DIST, UMEAN

C --------------------------------------------------------------------
C Now analyse the input parameters and return the corresponding result
C for the current point in the geometry. This has to be specified
C by the user depending in the current geometry and simulation!
C
C In short words: Find out, which NS-information FEATFLOW 
C wants to know and put this to the variables in the /NSPTS/
C COMMON block!
C --------------------------------------------------------------------

C=======================================================================
C *** Points for velocity, pressure and flux-tracing
C=======================================================================

C     In KPU(1) and KPU(2) it's possible to specify two points where the
C     velocities are measured during nonstationary simulations.
C     By setting KPU(1) / KPU(2) to the number of a vertex of the coarse
C     grid, the values of the velocities in these vertices are tracked
C     over time. 
C     According to INPWFS, the X/Y-velocities in KPU(1) is written into
C     the files "#points/u1_0" / "#points/u2_0".
C     The X/Y-velocities in the vertex KPU(2) is written into the files 
C     "#points/u3_0" / "#points/u4_0".

      KPU(1)=77
      KPU(2)=298

C     Like KPU, KPP(1..4) allowes to track the value of the pressure in
C     four different vertices over time. The pressure values are written
C     to external files. When KPP(.) is set to a vertex number on the
C     coarse grid, the pressure values are written out as follows:
C     Pressure in KPP(1) is written to "#points/p1_0".
C     Pressure in KPP(2) is written to "#points/p2_0".
C     Pressure in KPP(3) is written to "#points/p3_0".
C     Pressure in KPP(4) is written to "#points/p4_0".

      KPP(1)=77
      KPP(2)=298
      KPP(3)=72
      KPP(4)=74

C     KPX allowes to track the values of the streamfunction in two
C     vertieces. By setting KPX(1) and KPX(2) to two numbers of vertices
C     on the coarse grid, the difference KPX(1)-KPX(2) of the
C     streamfunction values in these points is written into the
C     external file "#points/f1_0".

      KPX(1)=7
      KPX(2)=298

C=======================================================================
C *** Parameters for 2 integral pressures in bdpres.f
C=======================================================================

C     By setting KPI(1) and/or KPI(2) to a value <> 0, it's possible
C     to calculate integrals on two boundary components.
C     Set rule is as follows:
C       KPI(1)   = number of boundary component where a boundary integral
C                  should be computed.
C       DPI(1,1) = Minimum parameter value on the boundary
C       DPI(2,1) = Maximum parameter value on the boundary
C
C     -> Boundary integram computed on boundary component KPI(1) in
C        parameter value interval [DPI(1,1),DPI(2,1)]
C
C     The similar holds for the second set KPI(2)/DPI(1,2)/DPI(2,2) which
C     allow to specify a second boundary component where to calculate a
C     boundary integral.
C
C     The integral values of the first boundary component specified
C     by KPI(1)/DPI(1,1)/DPI(2,1) are written into the file
C     "#points/p5_0".

      KPI(1)  =2
      DPI(1,1)=0D0
      DPI(2,1)=1D0

C     The integral values of the first boundary component specified
C     by KPI(2)/DPI(1,2)/DPI(2,2) are written into the file
C     "#points/p6_0".

      KPI(2)  =2
      DPI(1,2)=0D0
      DPI(2,2)=0.5D0

C=======================================================================
C *** Parameters for lift (DFW) and drag (DAW) in bdforc.f (INPR=2)
C ***
C *** dfw=2 int_s [dpf(1) dut/dn n_y - p n_x] ds / dpf(2)
C *** daw=2 int_s [dpf(1) dut/dn n_x + p n_y] ds / dpf(2)
C ***
C=======================================================================

C     These parameters are very much dependent on the current
C     configuration. We implement the standard cases here.
C
C     Benchmark-like configuration
C
C     RHO    = density of the fluid; normally normed to 1D0
C     UMEAN  = mean velocity of the parabolic profile
C     DIST   = length of the obstacle that is facing the flow;
C              depending on the direction of the flow!

C      RHO  =1.0D0
C      DIST =0.1D0
C      UMEAN=0.2D0

      RHO   = 1.0D0
      DIST  = DCRADY*2D0
      
C     Inflow velocity = 0.3
      
      UMEAN = 0.3D0 * 2D0/3D0

C     Write the coefficients into DPF(1)/DPF2. 

      DPF(1) = RHO*DNU
      DPF(2) = RHO*DIST*UMEAN**2

C     The values of Drag and Lift that are calculated using these
C     constants are written into the file "#points/p7_0".

99999 END
