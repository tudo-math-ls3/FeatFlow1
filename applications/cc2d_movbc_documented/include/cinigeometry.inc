************************************************************************
* This file contains the variable representing the user provided
* information in the DAT file of the application.
* More specifically the variables here define the behavior of
* the (moving-boundary-) geometry.
************************************************************************

      INTEGER ICIRTP, IFBMNT
      COMMON /GEOINT/ ICIRTP,IFBMNT
      
      DOUBLE PRECISION DCXPOS,DCYPOS,DCRPOS,DCRADX,DCRADY,DCRAD,
     *                 DCROT, DCRSIN, DCRCOS
      COMMON /GEODBL/ DCXPOS,DCYPOS,DCRPOS,DCRADX,DCRADY,DCRAD,
     *                DCROT, DCRSIN, DCRCOS
      
      SAVE /GEOINT/, /GEODBL/
      
C ICIRTP - Type of obstacle (circle, square, ellipse, ...)
C DCRAD  - Radius of a circle / half length of an edge on a square
C DCRADX - X-Radius of an ellipse
C DCRADY - Y-Radius of an ellipse
C DCROT  - Roration of the obstacle, 0..360°.
C
C For saving computational time the sin/cos of the rotation is
C precalculated:
C DCRSIN - SIN(DCROT); must always be updated if DCROT is changed
C DCRCOS - COS(DCROT); must always be updated if DCROT is changed
