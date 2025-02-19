************************************************************************
* This file contains the COMMON block /XYPAR/ in which the parameters
* for the Y-routines have to be stored before calling them.
* As this COMMON block is an extension to the multigrid common blocks,
* the file CBASICMG.INC has to be included before!
************************************************************************

      DOUBLE PRECISION DXYPAR
      INTEGER KXYPAR

      COMMON /XYPAR/  DXYPAR(NNARR),KXYPAR(NNARR)
      
      SAVE /XYPAR/
