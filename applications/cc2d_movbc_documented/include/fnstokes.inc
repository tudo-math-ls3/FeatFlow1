***********************************************************************
* This file contains the filenames of the stokes-matrices.
*
* As this file defines the filenames for all levels, the file
* CBASICMG.INC has to be included before!
***********************************************************************

      CHARACTER CFILST*12
      DIMENSION CFILST(NNLEV)
      
      DATA CFILST/'#ns/ST1     ','#ns/ST2     ','#ns/ST3     ',
     *            '#ns/ST4     ','#ns/ST5     ','#ns/ST6     ',
     *            '#ns/ST7     ','#ns/ST8     ','#ns/ST9     '/
      
      SAVE CFILST      