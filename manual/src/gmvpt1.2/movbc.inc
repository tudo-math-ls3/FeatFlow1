int Test_MovBC(Daten_Typ K[DIM])
{
  const Daten_Typ DVel=M_PI_4;
  const Daten_Typ DRad=0.226;
  const Daten_Typ DEps=0.05;
  Daten_Typ Phi1,Phi2,Phi3,Phi4,Phi5,Phi6,Phi7;
  Daten_Typ Dist1,Dist2,DRad2;

  if(fabs(K[3])<=0.226) {
    Phi1=Global_Time*DVel;
    Phi2=Phi1+M_PI_2;
    Phi3=Phi2+M_PI_2;
    Phi4=Phi3+M_PI_2;
    DRad2=sqrt(K[0]*K[0]+K[1]*K[1]);
    if(DRad2<=DRad) {
      if(DRad2<=1E-6) return 1; /**/
      Phi5=acos(K[0]/DRad2);
      if(K[1]<0.0) Phi5=2.0*M_PI-Phi5;
      while(Phi5<Phi1) Phi5+=2.0*M_PI;
      if((Phi5>=Phi1)&&(Phi5<=Phi2)) {
        Phi6=Phi1;
        Phi7=Phi2;
      } else if((Phi5>=Phi2)&&(Phi5<=Phi3)){
        Phi6=Phi2;
        Phi7=Phi3;
      } else if((Phi5>=Phi3)&&(Phi5<=Phi4)){
        Phi6=Phi3;
        Phi7=Phi4;
      } else if(Phi5>=Phi4){
        Phi6=Phi4;
        Phi7=Phi1;
      } else {
        fprintf(stderr,"Schrott!!!!!! %.8f, %.8f\n",Phi5,Phi1);
	return 0;/**/
      }
      Dist1=DRad2*sin(Phi5-Phi6);
      Dist2=DRad2*sin(Phi7-Phi5);
      if((fabs(Dist1)<=DEps)||(fabs(Dist2)<=DEps)) return 1; 
    }
  }
  return 0;
}
