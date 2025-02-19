C ********************************************************************
C This file introduces subroutines for easy reading of parameters from
C a DAT file.
C
C Every parameter read from the file is written to screen and/or
C to the output file.
C ********************************************************************
      
C Read INTEGER-value from file IUNIT to INT. Write it to the
C screen/file after the output string NAME.

      SUBROUTINE GETINT(IUNIT,NAME,INT)
      IMPLICIT NONE
      INCLUDE 'cout.inc'
      INCLUDE 'cfiles.inc'

      CHARACTER NAME*(*)
      INTEGER IUNIT
      INTEGER INT
      
      READ(IUNIT,*) INT
      IF (MT.GE.0) WRITE(MFILE1,'(A,I10)') NAME,INT
      IF (MT.GE.2) WRITE(MTERM,'(A,I10)') NAME,INT
      
      END

C Read Integer-value from file IUNIT. Interpret it as
C boolean: 0=false, <> 0 = true.
C (So boolean values are treated as integers in the file!)
C
C Write the value to the screen/file after the output string NAME.

      SUBROUTINE GETBOL(IUNIT,NAME,BOL)
      IMPLICIT NONE
      INCLUDE 'cout.inc'
      INCLUDE 'cfiles.inc'

      CHARACTER NAME*(*)
      INTEGER IUNIT
      LOGICAL BOL
      
      INTEGER INT
      
      READ(IUNIT,*) INT
      BOL = INT.NE.0
      IF (MT.GE.0) WRITE(MFILE1,'(A,I10)') NAME,INT
      IF (MT.GE.2) WRITE(MTERM,'(A,I10)') NAME,INT
      
      END

C Read DOUBLE-value from file IUNIT to DBL. Write it to the
C screen/file after the output string NAME.

      SUBROUTINE GETDBL(IUNIT,NAME,DBL)
      IMPLICIT NONE
      INCLUDE 'cout.inc'
      INCLUDE 'cfiles.inc'

      CHARACTER NAME*(*)
      INTEGER IUNIT
      DOUBLE PRECISION DBL
      
      READ(IUNIT,*) DBL
      IF (MT.GE.0) WRITE(MFILE1,'(A,D25.10)') NAME,DBL
      IF (MT.GE.2) WRITE(MTERM,'(A,D25.10)') NAME,DBL
      
      END

C Read STRING-value from file IUNIT to STR. Write it to the
C screen/file after the output string NAME.

      SUBROUTINE GETSTR(IUNIT,NAME,STR)
      IMPLICIT NONE
      INCLUDE 'cout.inc'
      INCLUDE 'cfiles.inc'

      CHARACTER NAME*(*),STR*(*)
      INTEGER IUNIT
      
      READ(IUNIT,*) STR
      IF (MT.GE.0) WRITE(MFILE1,'(A,A)') NAME,STR
      IF (MT.GE.2) WRITE(MTERM,'(A,A)') NAME,STR
      
      END
