************************************************************************
* This file contains COMMON blocks that are used during the assembly
* of matrices. It's an extension to CBASICELEM.INC, so this has to be
* included prior to this file!
************************************************************************

      INTEGER KDFG, KDFL, IDFL
	
      COMMON /COAUX1/ KDFG(NNBAS),KDFL(NNBAS),IDFL

      SAVE /COAUX1/
	  