************************************************************************
      DOUBLE PRECISION FUNCTION PARX(T,IBCT)
************************************************************************
*   parxypre.f:
*     - contains the parametrization functions:
*           PARX, PARY, TMAX
*       for the preprocessing variant ('rdparm')
*     - the COMMON-block /TDATA/ is used frequently and has to be
*       initialized by the modul 'initxx.f'
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=299, NNWORK=1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
C
      DIMENSION KWORK(1)
      EQUIVALENCE (DWORK(1),KWORK(1))
C
      COMMON /TDATA/ NBCT_T,NNCOMP,NNNPAR,
     *               LNCOMP,LICPTR,LITYP,LNSPLN,LNPAR,
     *               LIPPTR,LXPAR,LYPAR
C
      DATA INIT  /0/
      SAVE
C=======================================================================
C *** IBCT dependent quantities
      NCOMP=KWORK(L(LNCOMP)+IBCT-1)
      ICPTR=KWORK(L(LICPTR)+IBCT-1)
C
C *** determine boundary component ICOMP of part IBCT and
C     the ICOMP-dependent quantities
C
      ICOMP=T+1.d0
      tdiff=T-ICOMP+1.d0
      if(tdiff.lt.0d0 .or. tdiff.ge.1d0) then
        write(*,*)'error in PARX (userpre): conflict with ICOMP'
        stop
      endif
      iicomp=ICPTR+ICOMP-1
      ITYP=KWORK(L(LITYP)+iicomp)
      NSPLIN=KWORK(L(LNSPLN)+iicomp)
      NPAR=KWORK(L(LNPAR)+iicomp)
      IPPTR=KWORK(L(LIPPTR)+iicomp)
C
C *** pointer for the LXPAR/LYPAR-arrays
      KXPAR=L(LXPAR)+IPPTR
      KYPAR=L(LYPAR)+IPPTR      
C=======================================================================
C   ITYP=1:  line      
C=======================================================================
      if(ITYP.eq.1) then
C      
        if(NSPLIN.ne.1 .or. NPAR.ne.2) then
          write(*,*)'error in PARX (userpre): conflict at ITYP=1'
          stop
        endif
        x1=DWORK(KXPAR)
        xdiff=DWORK(KXPAR+1)
        PARX=x1 + xdiff*tdiff        
C=======================================================================
C   ITYP=2:  part of a circle      
C=======================================================================        
      else if(ITYP.eq.2) then
C      
        if(NSPLIN.ne.1 .or. NPAR.ne.3) then
          write(*,*)'error in PARX (userpre): conflict at ITYP=2'
          stop
        endif
        xm=DWORK(KXPAR)
        r=DWORK(KXPAR+1)
        phi1=DWORK(KXPAR+2)
        phi2=DWORK(KYPAR+2)
        phi=phi1 + tdiff*(phi2-phi1)
        PARX=xm + r*COS(phi)
C=======================================================================
C   ITYP=3:  spline
C=======================================================================        
      else if(ITYP.eq.3) then
C      
C *** determine "number of the subspline -1"=ISPLIN and tt=t_tilde
        ttdif=tdiff*NSPLIN
        ISPLIN=ttdif
        tt=ttdif-ISPLIN
        if(tt.lt.0d0 .or. tt.ge.1d0 .or. ISPLIN.gt.NSPLIN-1) then
          write(*,*)'error in PARX (userpre): conflict with ISPLIN'
          stop
        endif
C      
C *** determine cubic spline-function values at 'tt'
        phi1=(2d0*tt*tt - 3d0*tt)*tt + 1d0
        phi2=(-2d0*tt + 3d0)*tt*tt
        phi3=(tt*tt - 2d0*tt + 1d0)*tt
        phi4=(tt - 1d0)*tt*tt
C      
C *** determine parameter constants
        kk=KXPAR + 4*ISPLIN
        c1=DWORK(kk)
        c2=DWORK(kk+1)
        c3=DWORK(kk+2)
        c4=DWORK(kk+3)
C
        PARX=c1*phi1+c2*phi2+c3*phi3+c4*phi4
C=======================================================================      
      else
        write(*,*)'error in PARX (userpre): wrong ITYP-value'
        stop
      endif
C
      END
************************************************************************
      DOUBLE PRECISION FUNCTION PARY(T,IBCT)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNARR=299, NNWORK=1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
C
      DIMENSION KWORK(1)
      EQUIVALENCE (DWORK(1),KWORK(1))
C
      COMMON /TDATA/ NBCT_T,NNCOMP,NNNPAR,
     *               LNCOMP,LICPTR,LITYP,LNSPLN,LNPAR,
     *               LIPPTR,LXPAR,LYPAR
C
      save
C=======================================================================
C *** IBCT dependent quantities
      NCOMP=KWORK(L(LNCOMP)+IBCT-1)
      ICPTR=KWORK(L(LICPTR)+IBCT-1)
C
C *** determine boundary component ICOMP of part IBCT and
C     the ICOMP-dependent quantities
C
      ICOMP=T+1.d0
      tdiff=T-ICOMP+1.d0
      if(tdiff.lt.0d0 .or. tdiff.ge.1d0) then
        write(*,*)'error in PARY (userpre): conflict with ICOMP'
        stop
      endif
      iicomp=ICPTR+ICOMP-1
      ITYP=KWORK(L(LITYP)+iicomp)
      NSPLIN=KWORK(L(LNSPLN)+iicomp)
      NPAR=KWORK(L(LNPAR)+iicomp)
      IPPTR=KWORK(L(LIPPTR)+iicomp)
C
C *** pointer for the LXPAR/LYPAR-arrays
      KXPAR=L(LXPAR)+IPPTR
      KYPAR=L(LYPAR)+IPPTR      
C=======================================================================
C   ITYP=1:  line      
C=======================================================================
      if(ITYP.eq.1) then
C      
        if(NSPLIN.ne.1 .or. NPAR.ne.2) then
          write(*,*)'error in PARY (userpre): conflict at ITYP=1'
          stop
        endif
        y1=DWORK(KYPAR)
        ydiff=DWORK(KYPAR+1)
        PARY=y1 + ydiff*tdiff        
C=======================================================================
C   ITYP=2:  part of a circle      
C=======================================================================        
      else if(ITYP.eq.2) then
C      
        if(NSPLIN.ne.1 .or. NPAR.ne.3) then
          write(*,*)'error in PARY (userpre): conflict at ITYP=2'
          stop
        endif
        ym=DWORK(KYPAR)
        r=DWORK(KXPAR+1)
        phi1=DWORK(KXPAR+2)
        phi2=DWORK(KYPAR+2)
        phi=phi1 + tdiff*(phi2-phi1)
        PARY=ym + r*SIN(phi)
C=======================================================================
C   ITYP=3:  spline
C=======================================================================        
      else if(ITYP.eq.3) then
C      
C *** determine "number of the subspline -1"=ISPLIN and tt=t_tilde
        ttdif=tdiff*NSPLIN
        ISPLIN=ttdif
        tt=ttdif-ISPLIN
        if(tt.lt.0d0 .or. tt.ge.1d0 .or. ISPLIN.gt.NSPLIN-1) then
          write(*,*)'error in PARY (userpre): conflict with ISPLIN'
          stop
        endif
C      
C *** determine cubic spline-function values at 'tt'
        phi1=(2d0*tt*tt - 3d0*tt)*tt + 1d0
        phi2=(-2d0*tt + 3d0)*tt*tt
        phi3=(tt*tt - 2d0*tt + 1d0)*tt
        phi4=(tt - 1d0)*tt*tt
C      
C *** determine parameter constants
        kk=KYPAR + 4*ISPLIN
        c1=DWORK(kk)
        c2=DWORK(kk+1)
        c3=DWORK(kk+2)
        c4=DWORK(kk+3)
C
        PARY=c1*phi1+c2*phi2+c3*phi3+c4*phi4
C=======================================================================      
      else
        write(*,*)'error in PARY (userpre): wrong ITYP-value'
        stop
      endif
C
      END
************************************************************************
      DOUBLE PRECISION FUNCTION TMAX(IBCT)
************************************************************************
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNARR=299, NNWORK=1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
C
      DIMENSION KWORK(1)
      EQUIVALENCE (DWORK(1),KWORK(1))
C
      COMMON /TDATA/ NBCT_T,NNCOMP,NNNPAR,
     *               LNCOMP,LICPTR,LITYP,LNSPLN,LNPAR,
     *               LIPPTR,LXPAR,LYPAR
C
      save
C=======================================================================
      NCOMP=KWORK(L(LNCOMP)+IBCT-1)
C
      TMAX=NCOMP
      END
