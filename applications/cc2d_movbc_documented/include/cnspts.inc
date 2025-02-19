************************************************************************
* This file contains additional COMMON blocks for the Navier Stokes
* solver. They are mainly used for implementing boundary conditions.
************************************************************************

      INTEGER KPU, KPP, KPX, KPI
      DOUBLE PRECISION DPI, DPF

      COMMON /NSPTS/  KPU(2),KPP(4),KPX(4),KPI(2),DPI(2,2),DPF(2)
	
      SAVE /NSPTS/
	