************************************************************************
* This file contains the COMMON blocks that are used by elements to
* return calculated values.
*
* It's an extension to files CBASICELEM.INC and CBASICTRIA.INC, which
* have to be included prior to this file!
************************************************************************

      DOUBLE PRECISION DX, DY, DBAS, DJAC, DETJ
      DOUBLE PRECISION COB
      INTEGER KVE, IEL
      LOGICAL BDER

      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     *                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /COFBAS/ COB(NNBAS,NNCOF)

      SAVE /ELEM/, /COFBAS/