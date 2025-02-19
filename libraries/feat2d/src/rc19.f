************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* XRC19                                                                *
*                                                                      *
* Purpose  Call of RC19                                                *
*          Compress DWORK if desired                                   *
*          Double precision version                                    *
*                                                                      *
* Subroutines/functions called   RC19, ZDISP, ZNEW                     *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LA       I*4    Numbers of the corresponding arrays                  *
* LCOL     I*4    set by ZNEW                                          *
* LLD      I*4                                                         *
* IDISP    I*4    1  means call of ZDISP after deletion of zero entries*
* ARR1     C*6    Names of blocks (for error messages only)            *
* ARR2     C*6    Name of array corresponding to LCOL                  *
* For the description of the remaining parameters see RC19             *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                 -114 Type of at least one array is not double prec.  *
*                 -115 Length of at least one array is < NA            *
*                                                                      *
************************************************************************
*                                                                      *
* RC19                                                                 *
*                                                                      *
* Purpose  Removing entries of small modulus from an array             *
*          Storage technique 9                                         *
*          Double precision version                                    *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DA       R*8    Field of matrices                                    *
* KCOL     I*4    Pointer vectors                                      *
* KLD      I*4                                                         *
* NA       I*4    length of DA                                         *
* NEQ      I*4    number of equations                                  *
* NBLOC    I*4    Number of matrices stored on DA                      *
* KOFF     I*4    Matrix IBLOC starts at position KOFF(IBLOC)+1 on DA  *
* TOL      R*8    entries of  modulus <= TOL  are deleted              *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DA       R*8    Compressed matrices                                  *
* KCOL     I*4    described in the same storage technique              *
* KLD      I*4                                                         *
* NA       I*4                                                         *
*                                                                      *
************************************************************************
C
      SUBROUTINE XRC19 (LA,LCOL,LLD,NA,NEQ,NBLOC,TOL,IDISP,ARR1,ARR2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR1*6,ARR2*6
      PARAMETER (NNARR=299)
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION ARR1(*),LA(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='XRC19'
      IF (ICHECK.GE.997) CALL OTRC('XRC19 ','01/02/89')
      IER=0
C
      DO 1 IBLOC=1,NBLOC
C *** Check input parameter
      CALL ZTYPE(LA(IBLOC),ITYPE)
      IF (ITYPE.NE.1) THEN
       WRITE (CPARAM,'(A6,I15)') ARR1(IBLOC),IBLOC
       CALL WERR(-114,'XRC19 ')
       GOTO 99999
      ENDIF
      CALL ZLEN(LA(IBLOC),ILEN)
      IF (ILEN.LT.NA) THEN
       WRITE (CPARAM,'(A6,I15)') ARR1(IBLOC),IBLOC
       CALL WERR(-115,'XRC19 ')
       GOTO 99999
      ENDIF
1     CONTINUE
C
      CALL ZNEW(NBLOC,-3,LOFF,'KOFF  ')
      IF (IER.NE.0) GOTO 99999
C
      DO 2 IBLOC=1,NBLOC
      KWORK(L(LOFF)+IBLOC-1)=L(LA(IBLOC))-1
2     CONTINUE
C
      CALL RC19(DWORK(1),KWORK(L(LCOL)),KWORK(L(LLD)),NA,NEQ,NBLOC,
     *          KWORK(L(LOFF)),TOL)
C
      CALL ZDISP(0,LOFF,'KOFF  ')
      IF (IDISP.EQ.1) THEN
       DO 3 IBLOC=1,NBLOC
3      CALL ZDISP(NA,LA(IBLOC),ARR1(IBLOC))
       IF (IER.NE.0) GOTO 99999
       CALL ZDISP(NA,LCOL,ARR2)
      ENDIF
C
99999 END
C
C
C
      SUBROUTINE RC19(DA,KCOL,KLD,NA,NEQ,NBLOC,KOFF,TOL)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL (B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DA(*),KLD(*),KCOL(*),KOFF(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      IF (ICHECK.GE.997) CALL OTRC('RC19  ','01/02/89')
C
      BMSG2=M.GE.2.OR.MT.GE.2
C
      IP=1
      NAOLD=NA
C
      DO 10 IROW=1,NEQ
      BZERO=.TRUE.
      ILD=KLD(IROW)
      KLD(IROW)=IP
      DO 20 J=ILD,KLD(IROW+1)-1
      DO 21 IBLOC=1,NBLOC
      IF (ABS(DA(J+KOFF(IBLOC))).GT.TOL) GOTO 22
21    CONTINUE
      GOTO 20
22    BZERO=.FALSE.
      DO 23 IBLOC=1,NBLOC
23    DA(IP+KOFF(IBLOC))=DA(J+KOFF(IBLOC))
      KCOL(IP)=KCOL(J)
      IP=IP+1
20    CONTINUE
      IF (BZERO) THEN
C ***  Row IROW contains only zero entries ***
       IF (BMSG2) THEN
        WRITE (CPARAM,'(I15)') IROW
        CALL OMSG(54,'RC19  ')
       ENDIF
       DO 24 IBLOC=1,NBLOC
24     DA(IP+KOFF(IBLOC))=0D0
       KCOL(IP)=KCOL(ILD)
       IP=IP+1
      ENDIF
10    CONTINUE
C
      KLD(NEQ+1)=IP
      NA=IP-1
      WRITE (CPARAM,'(2I15)') NAOLD,NA
      CALL OMSG(21,'RC19  ')
C
      END
