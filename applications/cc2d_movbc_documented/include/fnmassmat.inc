***********************************************************************
* This file contains the filenames of the stokes-matrices.
*
* As this file defines the filenames for all levels, the file
* CBASICMG.INC has to be included before!
***********************************************************************

      CHARACTER CFILM*12
      DIMENSION CFILM(NNLEV)
      
      DATA CFILM /'#ns/MA1     ','#ns/MA2     ','#ns/MA3     ',
     *            '#ns/MA4     ','#ns/MA5     ','#ns/MA6     ',
     *            '#ns/MA7     ','#ns/MA8     ','#ns/MA9     '/
      
      SAVE CFILM     