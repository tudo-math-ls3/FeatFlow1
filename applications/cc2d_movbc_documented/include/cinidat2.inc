************************************************************************
* This file contains additional variables the user provided by
* the DAT file. It's an extension to the file CINIDAT.INC.
************************************************************************

C Treshold parameters for prolongation/Restriction

      DOUBLE PRECISION DPREP, DMTEP
      INTEGER IAVPR, IAPRM

      COMMON /RPARM2/ DPREP, DMTEP
      COMMON /IPARM2/ IAVPR, IAPRM

      SAVE /RPARM2/,/IPARM2/
