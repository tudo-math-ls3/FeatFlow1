************************************************************************
* This file contains the parametrisation routines for the
* computational domain:
*   PARX - Obtain the X-coordinate of a point from its parameter value
*   PARX - Obtain the Y-coordinate of a point from its parameter value
*   TMAX - Obtain the maximum parameter value on a boundary component
*
* The routines here are basically wrapper routines choosing the
* "real" parametrisation routines. The governing variable is the
* variable IMESH of the DAT file:
* If IMESH=1, these routines will use the standard OMEGA implementation 
*             routines for defining the computational domain by an 
*             external file.
* If IMESH=0, these routines will use the user defined FEAT
*             implementation FPARX, FPARY, FTMAX that is provided
*             in the file PARQ2D.F.
* Compatibility remark: If importing an "old" PARQ2D.F into the
* new implementation, the "old" PARX/PARY/TMAX-routines have simply
* to be renamed to FPARX/FPARY/FTMAX in order to work.
************************************************************************

      DOUBLE PRECISION FUNCTION PARX(T,IBCT)

      IMPLICIT NONE
      
      INCLUDE 'cfiles.inc'
      
C parameters

      DOUBLE PRECISION T
      INTEGER IBCT
      
C externals

      DOUBLE PRECISION FPARX,OPARX
      EXTERNAL FPARX,OPARX
      
C call the "correct" PARX-implementation

      IF (IMESH1.EQ.0) THEN
        PARX = FPARX(T,IBCT)
      ELSE IF (IMESH1.EQ.1) THEN
        PARX = OPARX (T,IBCT)
      END IF
      
      END


************************************************************************
      DOUBLE PRECISION FUNCTION PARY(T,IBCT)

      IMPLICIT NONE
      
      INCLUDE 'cfiles.inc'

C parameters

      DOUBLE PRECISION T
      INTEGER IBCT
      
C externals

      DOUBLE PRECISION FPARY,OPARY
      EXTERNAL FPARY,OPARY

C call the "correct" PARY-implementation

      IF (IMESH1.EQ.0) THEN
        PARY = FPARY(T,IBCT)
      ELSE IF (IMESH1.EQ.1) THEN
        PARY = OPARY (T,IBCT)
      END IF

      END


************************************************************************
      DOUBLE PRECISION FUNCTION TMAX(IBCT)

      IMPLICIT NONE
      
      INCLUDE 'cfiles.inc'

C parameters

      INTEGER IBCT
      
C externals

      DOUBLE PRECISION FTMAX,OTMAX
      EXTERNAL FTMAX,OTMAX

C call the "correct" TMAX-implementation

      IF (IMESH1.EQ.0) THEN
        TMAX = FTMAX(IBCT)
      ELSE IF (IMESH1.EQ.1) THEN
        TMAX = OTMAX (IBCT)
      END IF

      END
