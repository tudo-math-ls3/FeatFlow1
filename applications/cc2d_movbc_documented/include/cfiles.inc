************************************************************************
* This file contain COMMON blocks describing the input/output
* files used by the solver - start vectors, solution vectors,...
************************************************************************

      CHARACTER CPARM1*60,CMESH1*60,CFILE1*60,CSTART*60,CSOL*60
      INTEGER IMESH1,MFILE1,ISTART,MSTART,ISOL,MSOL
      
      COMMON /FILES/ IMESH1,CPARM1,CMESH1,MFILE1,CFILE1,
     *               ISTART,MSTART,CSTART,ISOL,MSOL,CSOL

      SAVE /FILES/
