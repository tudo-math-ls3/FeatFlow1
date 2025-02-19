***********************************************************************
* This file introduces a wrapper element EXXX ehich acts as a
* replacement for E030, E031, EM30, EM31. It calls the correct
* element routine by choice of the parameter IELT in the
* COMMON block.
***********************************************************************

***********************************************************************
* Nonconforming FE wrapper routine.
*
* Calls nonconforming element routine depending on parameter
* in the DAT file.
***********************************************************************

      SUBROUTINE EXXXNC(XI1,XI2,IPAR)

      IMPLICIT NONE
      
      DOUBLE PRECISION XI1,XI2
      INTEGER IPAR
      
      INCLUDE 'cinidat.inc'
      
      IF (IELT.EQ.0) CALL E031(XI1,XI2,IPAR)
      IF (IELT.EQ.1) CALL E030(XI1,XI2,IPAR)
      IF (IELT.EQ.2) CALL EM31(XI1,XI2,IPAR)
      IF (IELT.EQ.3) CALL EM30(XI1,XI2,IPAR)

      END
      