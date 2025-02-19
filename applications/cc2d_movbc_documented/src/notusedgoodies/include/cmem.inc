************************************************************************
* This file contains the public COMMON block and the definitions
* responsible for pseudodynamic memory management. 
*
* The parameter NNWORK is to be modified by the user. It describes
* the total amount of memory (in doubles) that is used by the
* main problem.
*
* NNARR describes the total amount of handles available to the program.
************************************************************************

      INTEGER NNWORK, NNARR

C NNWORK = number of DOUBLE-values in the DWORK-block
C
C 64 Bit:
C --------
C 500 MB:
C      PARAMETER (NNWORK=35000000)
C 1 GB:
C      PARAMETER (NNWORK=70000000)
C 2 GB:
C      PARAMETER (NNWORK=140000000)
C 3 GB:
C      PARAMETER (NNWORK=210000000)
C 4 GB:
C      PARAMETER (NNWORK=280000000)
C fast 6 GB:
C      PARAMETER (NNWORK=400000000)

C 32 Bit:
C --------
C Dummy
      PARAMETER (NNWORK=1)
C 12.8 MB
C      PARAMETER (NNWORK=2000000)
C 64 MB
C      PARAMETER (NNWORK=10000000)
C 128 MB
C      PARAMETER (NNWORK=20000000)
C 256 MB
C      PARAMETER (NNWORK=40000000)
C 500 MB:
C      PARAMETER (NNWORK=70000000)
C 1 GB:
C      PARAMETER (NNWORK=140000000)
C 1.5 GB:
C      PARAMETER (NNWORK=210000000)
C 2 GB:
C      PARAMETER (NNWORK=280000000)
C nearly 3 GB:
C      PARAMETER (NNWORK=400000000)

C NNARR = max. number of handles that can be used:

      PARAMETER (NNARR=299)

C Now the declaration of the main global structures follow

      INTEGER NWORK, IWORK, IWMAX, L
      
      DOUBLE PRECISION DWORK
      REAL VWORK(1)
      INTEGER KWORK(1)
      
      COMMON  NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)

C Equivalence statement; essential for pseudodyn. memory management

      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
