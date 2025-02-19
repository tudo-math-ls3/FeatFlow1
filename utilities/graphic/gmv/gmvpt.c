#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

/*
 Schalter:
   PART3D: Aktivierung der 3D Unterstuetzung
   BOUSS,POWERLAW : Liest zusaetzliche Eingabefelder
   MOVBC: Aktiviere "Moving Boundaries" Unterstützung
   TIN: Schreibt 'Tracer Identification Numbers' heraus
   AGE: Schreibt das Alter der einzelnen Partikel heraus
   LDOUBLE: Verwendet intern Long Doubles als Standardtyp
   RTD: Schreibt Informationen über gelöschte Tracer in
        eine Logdatei mit dem Namen "${präfix}_rtd.csv".
        RTD impliziert TIN und AGE
 */

#ifdef RTD
#define TIN
#define AGE
#endif

#ifndef LDOUBLE
typedef double Daten_Typ;
#else
typedef long double Daten_Typ;
#endif
typedef long double LDaten_Typ;
typedef unsigned int Index_Typ;
typedef unsigned int Zaehler_Typ;
typedef Zaehler_Typ Anz_Typ;

#define BSize 1024 
char NOut[BSize],NIn[BSize];
#ifdef RTD
char Nrtd[BSize];
#ifdef LDOUBLE
const char *rtd_data="%.10Lg";
#else
const char *rtd_data="%.8g";
#endif
const char *rtd_idx="%u";
FILE *rtd_log;
#endif

Daten_Typ Global_Time=0.0;

#ifdef MOVBC
Daten_Typ Input_Time=0.0;
#endif

#ifdef PART3D
#define DIM 3 
#define NpC 8 
#define CpM 8 
const char CellStr[]="hex 8"; 
#else
#define DIM 2 
#define NpC 4 
#define CpM 4 
const char CellStr[]="quad 4";
#endif

#ifdef PART3D
#ifdef BOUSS
  #ifdef POWERLAW
    #define Num_Field 3
  #else
    #define Num_Field 2
  #endif
#else
  #ifdef POWERLAW
    #define Num_Field 2
  #else
    #define Num_Field 1
  #endif
#endif
#else
#ifdef BOUSS
  #ifdef POWERLAW
    #define Num_Field 4
  #else
    #define Num_Field 3
  #endif
#else
  #ifdef POWERLAW
    #define Num_Field 3
  #else
    #define Num_Field 2
  #endif
#endif
#endif  
const char *Field_Names[Num_Field]=
{
  "pressure 1"
#ifndef PART3D
  ,"streamfunction 1"
#endif
#ifdef BOUSS
  ,"temperature 1"
#endif
#ifdef POWERLAW
  ,"ny 1"
#endif
};

const char *Tracer_Names[Num_Field]=
{
  "Druck"
#ifndef PART3D  
  ,"Stromfkt"
#endif
#ifdef BOUSS
  ,"Temp"
#endif
#ifdef POWERLAW
  ,"ny"
#endif
};

typedef struct
{
  Anz_Typ Anz_Index;
  Index_Typ *I;
} Index_List;

typedef struct
{
  Anz_Typ Anz_Nodes;
  Daten_Typ *K[DIM];
} Nodes_Typ;

/* Moegliche Zustaende fuer Reg_Flag[i] */
#define REG_UNKNOWN 0 
#define REG_REGULAR 1 
#define REG_IRREGULAR 2 

typedef struct
{
  Anz_Typ Anz_Cells;
  Index_Typ *Vert;
  unsigned char *Reg_Flag;
  /* Formfaktoren fuer verbesserte Interpolation */
  /* Daten_Typ *Form; */
} Cells_Typ;

typedef struct 
{
  Nodes_Typ Nodes;
  Cells_Typ Cells;	
} Gitter_Typ;

typedef struct
{
  Gitter_Typ Gitter;
  Daten_Typ *V[DIM];
  Daten_Typ *Field[Num_Field];
} Input_Typ;

typedef struct
{
  Daten_Typ L[NpC];
} B_Coord_Typ;

typedef struct
{
  Anz_Typ Anz_Tracer;
  Daten_Typ *K[DIM];
  Daten_Typ *Farbe,*Expire;
#ifdef AGE
  Daten_Typ *Alter;
#endif
#ifdef TIN 
  Index_Typ *tin;
#endif
  Index_Typ *InCell;
  B_Coord_Typ *B_Coords;
} Tracer_Typ;

typedef struct
{
  Tracer_Typ Tracers;
  Daten_Typ *VecNorm;
  Daten_Typ *Field[Num_Field];
} Output_Typ;

typedef struct
{
  Daten_Typ P[2][DIM],Farbe,Haltbarkeit;
  Anz_Typ N[DIM];
  Anz_Typ Modul,Offset,Start,Stop;
} Block_Typ;

#ifdef TIN 
Index_Typ New_tin=1;
#endif

Output_Typ O;
Input_Typ I;
/* Die Einzelnen Ebenen des Mehrgitters */
Anz_Typ Anz_Level=0;
Gitter_Typ *MG;
/* Liste der "entarteten" Zellen */
Index_List *CL;

const Daten_Typ EPS=0.1E-4;
const Daten_Typ Small_EPS=0.5E-7;
Daten_Typ EpsScale=1.0;

/* Statistische Datenerhebung zu Debuggin-Zwecken */
Anz_Typ Tracer_Verworfen=0,Tracer_Gerettet=0;
Anz_Typ Tracer_Verworfen_Inc=0,Tracer_Gerettet_Inc=0;
/**/
 
void Lese_Daten(FILE *f,Anz_Typ Anz,Daten_Typ *Daten)
{
  Zaehler_Typ i;
  if(Daten==NULL)
    {
      for(i=0;i<Anz;i++) fscanf(f,"%*f");
    }
  else
    {
      for(i=0;i<Anz;i++,Daten++)
#ifndef LDOUBLE 
	fscanf(f,"%lf",Daten);
#else
	fscanf(f,"%Lf",Daten);
#endif
    }
}
void Schreibe_Daten(FILE *f,Anz_Typ Anz,Daten_Typ *Daten)
{
  Zaehler_Typ i;
  if(Daten==NULL)
    {
      for(i=0;i<Anz;i++) fprintf(f,"%.8f\n",0.0);
    }
  else
    {
      for(i=0;i<Anz;i++,Daten++)
#ifndef LDOUBLE
	fprintf(f,"%.8g\n",*Daten);
#else
	fprintf(f,"%.10Lg\n",*Daten);
#endif
    }
}
void Lese_Index(FILE *f,Anz_Typ Anz,Index_Typ *Index)
{
  Zaehler_Typ i;
  if(Index==NULL)
    {
      for(i=0;i<Anz;i++) fscanf(f,"%*u");
    }
  else
    {
      for(i=0;i<Anz;i++,Index++) fscanf(f,"%u",Index);
    }
}
void Schreibe_Index(FILE *f,Anz_Typ Anz,Index_Typ *Index)
{
  Zaehler_Typ i;
  if(Index==NULL)
    {
      for(i=0;i<Anz;i++) fprintf(f,"%u\n",0);
    }
  else
    {
      for(i=0;i<Anz;i++,Index++) fprintf(f,"%u\n",*Index);
    }
}
void Lese_Anz(FILE *f,Anz_Typ *Anz)
{
  if(Anz==NULL) fscanf(f,"%*u");
  else fscanf(f,"%u",Anz);
}
void Schreibe_Anz(FILE *f,Anz_Typ Anz)
{
  fprintf(f,"%u\n",Anz);
}

void Scan_Name(FILE *f,const char *S)
{
  fscanf(f,"\n");
  fscanf(f,S);
}

void Schreibe_Vertreter(FILE *f)
{
  fprintf(f,"gmvinput ascii\n");
  fprintf(f,"nodes fromfile \"coarse.gmv\"\n");
  fprintf(f,"cells fromfile \"coarse.gmv\"\n");
}

char * strip_string(char * str) {
  size_t i;
  if(!str) return str;
  // Entferne Leerzeichen am Anfang
  while((str[0]==' ')||(str[0]=='\t')) str ++;
  // Entferne das Zeilenende
  i=strlen(str);
  str[--i]=(char)(0);
  // Entferne Leerzeichen am Ende
  while((str[i-1]==' ')||(str[i-1]=='\t')) str[--i]=(char)(0);
  // Entferne Anführungszeichen  falls vorhanden
  if((str[0]=='\'')||(str[0]=='\"')) {
    str++;i-=2;
    str[i]=(char)(0);
  }
  return str;
}

void Lese_Knoten(FILE *f,Nodes_Typ *Nodes) {
   Anz_Typ A;
   Zaehler_Typ i;
   char Buffer[1024];
   char * Name;
   FILE *New_f;
   
   Scan_Name(f,"nodes");
   fgets(Buffer,1024,f);
   Name=strstr(Buffer,"fromfile");
   if(Name) {
     if(!Nodes) return;
     Name+=strlen("fromfile");
     Name=strip_string(Name);
     fprintf(stderr,"Verzweige nach %s !\n",Name);
     New_f=fopen(Name,"r");
     if(!New_f) {
       fprintf(stderr,"Verzweigen nach %s mißlungen!\n",Name);
       exit(EXIT_FAILURE);
     }
     Lese_Knoten(New_f,Nodes);
     fclose(New_f);
     return;
   }
   sscanf(Buffer,"%u",&A);
   if(Nodes!=NULL)
   {
     printf("Ich alloziere nun Speicher für %u Knoten...\n",A);
     Nodes->Anz_Nodes=A;
     for(i=0;i<DIM;i++) 
       Nodes->K[i]=
	 (Daten_Typ *)calloc(A+1,sizeof(Daten_Typ));
     printf("Ich lese nun %u Knoten ein...\n",A);
     for(i=0;i<3;i++)
       {
	 if(i<DIM) Lese_Daten(f,A,Nodes->K[i]+1);
	 else Lese_Daten(f,A,NULL);
       }
   }
   else
   {
     printf("Ich überlese nun %u Knoten...\n",A);
     Lese_Daten(f,3*A,NULL);
   }

}

void Lese_Zellen(FILE *f,Cells_Typ *Cells) {
   Anz_Typ A;
   Zaehler_Typ i;
   char Buffer[1024];
   char * Name;
   FILE *New_f;

   Scan_Name(f,"cells");
   fgets(Buffer,1024,f);
   Name=strstr(Buffer,"fromfile");
   if(Name) {
     if(!Cells) return;
     Name+=strlen("fromfile");
     Name=strip_string(Name);
     fprintf(stderr,"Verzweige nach %s !\n",Name);
     New_f=fopen(Name,"r");
     if(!New_f) {
       fprintf(stderr,"Verzweigen nach %s mißlungen!\n",Name);
       exit(EXIT_FAILURE);
     }
     Lese_Zellen(New_f,Cells);
     fclose(New_f);
     return;
   }
   sscanf(Buffer,"%u",&A);
   if(Cells!=NULL)
   {
     printf("Ich alloziere nun Speicher für %u Zellen...\n",A);
     Cells->Anz_Cells=A;
     Cells->Vert=(Index_Typ *)calloc((A+1)*NpC,sizeof(Index_Typ));
     Cells->Reg_Flag=(unsigned char *)calloc(A+1,sizeof(unsigned char));
     printf("Ich lese nun %u Zellen ein...\n",A);
     for(i=1;i<=A;i++)
     {
        Scan_Name(f,CellStr);
        Lese_Index(f,NpC,Cells->Vert+i*NpC);
     }
   }
   else
   {
     printf("Ich überlese nun %u Zellen...\n",A);
     for(i=1;i<=A;i++)
     {
       Scan_Name(f,CellStr);
       Lese_Index(f,NpC,NULL);
     }
   }
}
 
void Lese_Gitter(FILE *f,Gitter_Typ *Gitter)
{
  fscanf(f,"gmvinput ascii");
  if(Gitter) {
    Lese_Knoten(f,&(Gitter->Nodes));
    Lese_Zellen(f,&(Gitter->Cells));
  } else {
    Lese_Knoten(f,NULL);
    Lese_Zellen(f,NULL);
  }
}

void Update_Input(const char *Name,Input_Typ *Input)
{
  FILE *f;
  Anz_Typ A,i;
  f=fopen(Name,"r");
  Lese_Gitter(f,NULL);
  A=Input->Gitter.Nodes.Anz_Nodes;
  printf("Ich lese nun die Datensätze ein...\n");
  Scan_Name(f,"velocity 1");
  for(i=0;i<3;i++)
    {
      if(i<DIM) Lese_Daten(f,A,Input->V[i]+1);
      else Lese_Daten(f,A,NULL);
    }
  printf("Geschwindigkeitsfeld gelesen!\n");
  Scan_Name(f,"variable");
  for(i=0;i<Num_Field;i++)
    {
      Scan_Name(f,Field_Names[i]);
      Lese_Daten(f,A,Input->Field[i]+1);
      printf("Zusatzdatenfeld %d gelesen!\n",i+1);
    }
  fclose(f);
}

void Startup_Input(const char *Name,Input_Typ *Input)
{
  Anz_Typ A,i;
  Input->Gitter=MG[Anz_Level-1];
  A=Input->Gitter.Nodes.Anz_Nodes;
  printf("Ich alloziere Speicher für diverse Datenfelder\n");
  for(i=0;i<DIM;i++) 
    Input->V[i]=(Daten_Typ *)calloc(A+1,sizeof(Daten_Typ));
  for(i=0;i<Num_Field;i++) 
    Input->Field[i]=(Daten_Typ *)calloc(A+1,sizeof(Daten_Typ));
  if (Name!=NULL) Update_Input(Name,Input);
}
 
void Cleanup_Nodes(Nodes_Typ *N)
{
  Zaehler_Typ i;
  N->Anz_Nodes=0;
  for(i=0;i<DIM;i++) 
    if(N->K[i]!=NULL)
    {
      free(N->K[i]);
      N->K[i]=NULL;
    }
}

void Cleanup_Cells(Cells_Typ *C)
{
  C->Anz_Cells=0;
  if(C->Vert!=NULL)
  {
    free(C->Vert);
    C->Vert=NULL;
  }
  if(C->Reg_Flag!=NULL)
  {
    free(C->Reg_Flag);
    C->Reg_Flag=NULL;
  }
}

void Cleanup_Gitter(Gitter_Typ *G)
{
  Cleanup_Nodes(&(G->Nodes));
  Cleanup_Cells(&(G->Cells));
}

void Cleanup_Input(Input_Typ *Input)
{
  Anz_Typ i;
   if(Input!=NULL)
   {
     printf("Ich gebe nun den allozierten Speicher wieder frei...\n");
     for(i=0;i<DIM;i++)
       free(Input->V[i]);
     for(i=0;i<Num_Field;i++)
       free(Input->Field[i]);
     /* Gitter wird schon bei Cleanup_Suchen freigegeben! */
   }
}

Daten_Typ Dist_Max(Daten_Typ *V1, Daten_Typ *V2)
{
  Zaehler_Typ i;
  Daten_Typ d=0.0,temp;
  for(i=DIM;i>0;i--)
    {
      temp=fabs((*V1)-(*V2));
      V1++;V2++;
      if(d<temp) d=temp;
    }
  return d;
}

Daten_Typ Norm_Max(Daten_Typ *V)
{
  Zaehler_Typ i;
  Daten_Typ d,temp;
  for(i=DIM,d=0.0;i>0;i--)
    {
      temp=fabs(*V);
      V++;
      if(d<temp) d=temp;
    }
  return d;
}

Daten_Typ Dist_Euklid(Daten_Typ *V1, Daten_Typ *V2)
{
  Zaehler_Typ i;
  Daten_Typ d,temp;
  for(i=DIM,d=0.0;i>0;i--)
    {
      temp=(*V1)-(*V2);
      V1++;V2++;
      d+=temp*temp;
    }
  return sqrt(d);
}

Daten_Typ Norm_Euklid(Daten_Typ *V)
{
  Zaehler_Typ i;
  Daten_Typ d,temp;
  for(i=DIM,d=0.0;i>0;i--)
    {
      temp=*V;
      V++;
      d+=temp*temp;
    }
  return sqrt(d);
}

void Vec_Scale (Daten_Typ Faktor, Daten_Typ *V)
{
  Zaehler_Typ i;
  for(i=DIM;i>0;i--) *(V++)*=Faktor;
}

Daten_Typ Vec_Dot(Daten_Typ *V1,Daten_Typ *V2)
{
  Zaehler_Typ i;
  Daten_Typ Erg;
  for(i=DIM,Erg=0.0;i>0;i--) 
    {
      Erg+=(*V1)*(*V2);
      V1++;V2++;
    }
  return Erg;
}

#ifdef PART3D
Index_Typ transform3d(Daten_Typ X[8],Daten_Typ Y[8],Daten_Typ Z[8],
		      Daten_Typ P[3],Daten_Typ UVW[3]);

int Test_Box(const Daten_Typ X[NpC], const Daten_Typ Y[NpC],const Daten_Typ Z[NpC],
	     const Daten_Typ K[3])
{
  const Daten_Typ *H[3];
  Daten_Typ Hmin[3],Hmax[3],Dummy;
  Zaehler_Typ i,j;

  H[0]=X;H[1]=Y;H[2]=Z;
  for(j=0;j<DIM;j++)
    {
      Hmin[j]=Hmax[j]=(H[j])[0];
      for(i=1;i<NpC;i++)
	{
	  Dummy=(H[j])[i];
	  if(Dummy<Hmin[j]) Hmin[j]=Dummy;
	  if(Dummy>Hmax[j]) Hmax[j]=Dummy;
	}
      if(!(((Hmin[j]-EPS*EpsScale)<=K[j])&&
	   ((Hmax[j]+EPS*EpsScale)>=K[j]))) 
	return 0;
    }
  return 1;
}

int Test_Zelle(Index_Typ C,Gitter_Typ *G,Daten_Typ K[DIM])
{
  const Daten_Typ EX[NpC]={-1, 1, 1,-1,-1, 1, 1,-1};
  const Daten_Typ EY[NpC]={-1,-1, 1, 1,-1,-1, 1, 1};
  const Daten_Typ EZ[NpC]={-1,-1,-1,-1, 1, 1, 1, 1};
  Daten_Typ X[NpC],Y[NpC],Z[NpC],UVW[DIM];
  Zaehler_Typ i;
  for(i=0;i<NpC;i++)
    {
      X[i]=G->Nodes.K[0][G->Cells.Vert[i+NpC*C]];
      Y[i]=G->Nodes.K[1][G->Cells.Vert[i+NpC*C]];
      Z[i]=G->Nodes.K[2][G->Cells.Vert[i+NpC*C]];      
    }
  if(!Test_Box(X,Y,Z,K)) return 0;
  transform3d(X,Y,Z,K,UVW);
  return Test_Box(EX,EY,EZ,UVW);
}

unsigned char Is_Regular(Index_Typ C,Gitter_Typ *G)
{
  const Index_Typ T[8][3]=
  {
    {1,3,4},
    {0,2,5},
    {1,6,3},
    {0,2,7},
    {0,5,7},
    {1,6,4},
    {2,5,7},
    {3,6,4},
  };
  Daten_Typ Erg;
  Zaehler_Typ i,j,k;
  if(G->Cells.Reg_Flag[C]!=REG_UNKNOWN) return G->Cells.Reg_Flag[C];
  for(i=0;i<8;i++)
  {
    for(k=0;k<3;k++)
    {
      for(j=0,Erg=0.0;j<3;j++)
      {
        
        Erg+=(G->Nodes.K[j][G->Cells.Vert[NpC*C+T[i][k]]]
              -G->Nodes.K[j][G->Cells.Vert[NpC*C+i]])
             *(G->Nodes.K[j][G->Cells.Vert[NpC*C+T[i][(k+1)%3]]]
             -G->Nodes.K[j][G->Cells.Vert[NpC*C]+i]);
      }
      if(fabs(Erg)>=EPS)
      {
        G->Cells.Reg_Flag[C]=REG_IRREGULAR;
        return REG_IRREGULAR;
      }
    }
  }
  G->Cells.Reg_Flag[C]=REG_REGULAR;
  return REG_REGULAR;
}
#else

int Test_Box(const Daten_Typ X[NpC], const Daten_Typ Y[NpC],
	     const Daten_Typ K[DIM])
{
  const Daten_Typ *H[DIM];
  Daten_Typ Hmin[DIM],Hmax[DIM],Dummy;
  Zaehler_Typ i,j;

  H[0]=X;H[1]=Y;
  for(j=0;j<DIM;j++)
    {
      Hmin[j]=Hmax[j]=(H[j])[0];
      for(i=1;i<NpC;i++)
	{
	  Dummy=(H[j])[i];
	  if(Dummy<Hmin[j]) Hmin[j]=Dummy;
	  if(Dummy>Hmax[j]) Hmax[j]=Dummy;
	}
      if(!(((Hmin[j]-EPS*EpsScale)<=K[j])&&
	   ((Hmax[j]+EPS*EpsScale)>=K[j]))) 
	return 0;
    }
  return 1;
}

int Test_Linie(Gitter_Typ *Gitter,Index_Typ P1,Index_Typ P2,Daten_Typ K[2])
{
  Daten_Typ Normal[DIM],Vect[DIM],Erg;
  Zaehler_Typ i;
  for(i=0;i<DIM;i++)
    Vect[i]=K[i]-Gitter->Nodes.K[i][P1];
  /* Nach innen zeigender Normalenvektor */
  Normal[0]=-(Gitter->Nodes.K[1][P2]-Gitter->Nodes.K[1][P1]);
  Normal[1]=Gitter->Nodes.K[0][P2]-Gitter->Nodes.K[0][P1];
  /**/
  Erg=1.0/Norm_Euklid(Normal);
  Vec_Scale(Erg,Normal);
  /**/
  Erg=Vec_Dot(Normal,Vect);
  if(Erg>=-EPS*EpsScale) return 1;
  return 0;
}

int Test_Zelle(Index_Typ C,Gitter_Typ *Gitter,Daten_Typ K[DIM])
{
  Zaehler_Typ i;
  Index_Typ *D;
  Daten_Typ X[NpC],Y[NpC];
  for(i=0;i<NpC;i++)
    {
      X[i]=Gitter->Nodes.K[0][Gitter->Cells.Vert[i+NpC*C]];
      Y[i]=Gitter->Nodes.K[1][Gitter->Cells.Vert[i+NpC*C]];  
    }
  if(!Test_Box(X,Y,K)) return 0;
  D=Gitter->Cells.Vert+NpC*C;
  for(i=0;i<NpC;i++) 
  {
    if(!Test_Linie(Gitter,*(D+i),*(D+(i+1)%NpC),K)) return 0;
  }
  return 1;
}

unsigned char Is_Regular(Index_Typ C,Gitter_Typ *G)
{
  Daten_Typ V[NpC][2],Erg;
  Zaehler_Typ i,j;
  if(G->Cells.Reg_Flag[C]!=REG_UNKNOWN) return G->Cells.Reg_Flag[C];
  for(i=0;i<NpC;i++)
  {
    for(j=0;j<DIM;j++)
      {
	V[i][j]=G->Nodes.K[j][G->Cells.Vert[NpC*C+((i+1)%NpC)]]-
	  G->Nodes.K[j][G->Cells.Vert[NpC*C+i]];
      }
  }
  for(i=0;i<NpC;i++)
  {
    Erg=V[i][0]*V[(i+1)%NpC][0]+V[i][1]*V[(i+1)%NpC][1];
    if(fabs(Erg)>=EPS)
    {
      G->Cells.Reg_Flag[C]=REG_IRREGULAR;
      return REG_IRREGULAR;
    }
  }
  G->Cells.Reg_Flag[C]=REG_REGULAR;
  return REG_REGULAR;
}
#endif

#ifdef MOVBC
#include "movbc.inc"
#endif

int Datei_Existiert(const char *Name)
{
  FILE *f;
  f=fopen(Name,"r");
  if(f==NULL) return 0;
  fclose(f);
  return 1;
}

void File_Error(const char *Name)
{
  fprintf(stderr,"Datei %s nicht gefunden!\n",Name);
  exit(EXIT_FAILURE);
}

void Reset_EpsScale(void)
{
  EpsScale=1.0;
}
#ifdef PART3D
#define EpsScaleFac 4.0 
#else
#define EpsScaleFac 2.0
#endif

void Raise_EpsScale(void)
{
  EpsScale*=EpsScaleFac;
}
void Lower_EpsScale(void)
{
  EpsScale/=EpsScaleFac;
  if(EpsScale<1.0) Reset_EpsScale();
}
void Compute_EpsScale(Zaehler_Typ Level)
{
  Zaehler_Typ i;
  Reset_EpsScale();
  for(i=2;i<=Level;i++) Raise_EpsScale();
}

/* Suche Zelle auf Mehrgitterweise */
Index_Typ Suche_Zelle(Daten_Typ K[DIM],Index_Typ Old_Cell)
{
  Index_Typ i,Akt_Level,Akt_Cell,n;

  Akt_Level=Anz_Level;
  Reset_EpsScale();
  Akt_Cell=Old_Cell;
  if(Akt_Cell)
  {
    /* Vergröbere zuerst das Gitter nach bedarf max. bis Level 1*/
    while((Akt_Level>1)&&!Test_Zelle(Akt_Cell,MG+Akt_Level-1,K))
    {
      Akt_Level--;
      if(Akt_Cell>MG[Akt_Level-1].Cells.Anz_Cells)
        Akt_Cell=(Akt_Cell-MG[Akt_Level-1].Cells.Anz_Cells-1)/(CpM-1)+1;
    }
    /* Vollstaendige Suche, falls erfolglos */
    if(!Test_Zelle(Akt_Cell,MG+Akt_Level-1,K)) Akt_Cell=0;
  }
  else Akt_Level=1;
  /* Mache solange vollständige Suchen in den entarteten Zellen, 
     bis ein Einstiegsmakro gefunden wurde oder feinster Level erreicht ist.
     Worst Case Fall tritt ein, wenn ein Tracer sich aus dem Gebiet
     herausbewegt!
     */
 loop1:
  while((Akt_Level<=Anz_Level)&&(!Akt_Cell))
  {
    Compute_EpsScale(Akt_Level);
    for(i=0;i<CL[Akt_Level-1].Anz_Index;i++)
    {
      if(Test_Zelle(CL[Akt_Level-1].I[i],MG+Akt_Level-1,K))
      {
        Akt_Cell=CL[Akt_Level-1].I[i];
        break;
      }
    }
    if(!Akt_Cell) Akt_Level++;
  }
  /* Tracer befindet sich in keiner Zelle */
  if(Akt_Level>Anz_Level)
    {
      Akt_Cell=0;
      goto mark1;
    }
  /* Verfeinere das Gitter Schrittweise bis zum feinsten Level */
  while(Akt_Level<Anz_Level)
  {
    n=0;
    Akt_Level++;
    Raise_EpsScale();
    if(Test_Zelle(Akt_Cell,MG+Akt_Level-1,K)) n=1;
    else
      {
	for(i=1;i<CpM;i++)
	  {
	    if(Test_Zelle((Akt_Cell-1)*(CpM-1)+MG[Akt_Level-2].Cells.Anz_Cells+i,
			  MG+Akt_Level-1,K))
	      {
		Akt_Cell=(Akt_Cell-1)*(CpM-1)+MG[Akt_Level-2].Cells.Anz_Cells+i;
		n=1;
		break;
	      }
	  }
      }
    /* Wenn es nicht in einer Zelle der aktuellen Hirarchie ist, dann kann
       immer noch in einer entarteten Zelle des naechstfeineren Gitters sein.
    */
    if(!n)
      {
	Akt_Cell=0;
	Akt_Level++;
	goto loop1;
      }
  }
  if(Akt_Level!=Anz_Level) fprintf(stderr,"Schwerer Fehler beim Suchen!\n");
  /* Dies hier ist sowohl alt, als auch neu.
     Um einen seltsamen Fehler mit verschwindenden Partikeln abzufangen,
     mache ich eine vollstaendige Suche fuer Partikel, die angeblich in
     keiner Zelle liegen.
  */
  if (Akt_Cell==0)
    {
    mark1:
      Tracer_Verworfen_Inc++;
      Reset_EpsScale();
      for(i=1;i<=MG[Anz_Level-1].Cells.Anz_Cells;i++)
	{
	  if(Test_Zelle(i,MG+Anz_Level-1,K))
	    {
	      Akt_Cell=i;
	      Tracer_Gerettet_Inc++;
	      break;
	    }
	}
    }
  /* Ende des neuen Teiles. */
  return Akt_Cell;
}

void Startup_Suchen(void)
{
  Zaehler_Typ i,j,k,l;
  Index_Typ OldCell;
  Anz_Typ *A;
  Anz_Typ n=0;
  Daten_Typ K[DIM];
  FILE *f;
  Index_List L;
  L.I=NULL;L.Anz_Index=0;
  /* Zuerst Bestimme ich die Anzahl der Verfeinerungstufen des MG */
  i=0;
  do {sprintf(NIn,"#gmv/tria%d.gmv",++i);} while(Datei_Existiert(NIn));
  Anz_Level=--i;
  MG=(Gitter_Typ *)calloc(Anz_Level,sizeof(Gitter_Typ));
  A=(Anz_Typ *)calloc(Anz_Level,sizeof(Anz_Typ));
  for(i=1;i<=Anz_Level;i++)
  {
    sprintf(NIn,"#gmv/tria%d.gmv",i);
    f=fopen(NIn,"r");
    Lese_Gitter(f,MG+i-1);
    fclose(f);
    A[i-1]=MG[i-1].Nodes.Anz_Nodes;
    if(i<Anz_Level) Cleanup_Nodes(&(MG[i-1].Nodes));
  }
  for(i=0;i<Anz_Level-1;i++) 
    {
      MG[i].Nodes=MG[Anz_Level-1].Nodes;
      MG[i].Nodes.Anz_Nodes=A[i];
    }
  free(A);
  /* Debug */
  printf("Bestimme entartete Zellen...\n");
  /**/
  /* Bestimme entartete Zellen,d.h. alle Zellen, die nicht vollständig
     in einem Makro der naechstgröberen Verfeinerungstufe liegen */
  /* Debug */
  Reset_EpsScale();
  printf("Level 1\n");
  /**/
  CL=(Index_List *)calloc(Anz_Level,sizeof(Index_List));
  /* Im gröbsten Level sind alle Zellen entartet */
  CL[0].I=(Index_Typ *)calloc(MG[0].Cells.Anz_Cells,sizeof(Index_Typ));
  CL[0].Anz_Index=MG[0].Cells.Anz_Cells;
  for(i=0;i<CL[0].Anz_Index;i++) CL[0].I[i]=i+1;
  /* Bearbeite Level 2..Anz_Level */
  for(i=2;i<=Anz_Level;i++)
  {
    /* Debug */
    printf("Level %d\n",i);
    /**/
    /* Lösche Knotenliste */
    L.Anz_Index=0;
    L.I=(Index_Typ *)realloc(L.I,0);
    /* Teste alle neuen Knoten im aktuellen Level daraufhin, ob sie nicht
       innerhalb des nächstgröberen Gitters liegen mittels
       vollständiger Suche.
    */
    OldCell=0;
    for(j=MG[i-2].Nodes.Anz_Nodes+1;j<=MG[i-1].Nodes.Anz_Nodes;j++)
    {
      for(k=0;k<DIM;k++) K[k]=MG[i-1].Nodes.K[k][j];
      /*
      n=0;
      for(k=1;k<=MG[i-2].Cells.Anz_Cells;k++)
        if(Test_Zelle(k,MG+i-2,K))
        {
          n=1;
          break;
        }
	*/
      /* Experimentell */ 
      l=Anz_Level;
      Anz_Level=i-1;
      OldCell=n=Suche_Zelle(K,OldCell);
      Anz_Level=l;
      /**/
      /* Wenn nicht gefunden, dann entarteter Knoten */
      if(!n)
      {
        L.I=(Index_Typ *)realloc(L.I,(++L.Anz_Index)*sizeof(Index_Typ));
        L.I[L.Anz_Index-1]=j;
      }
    }
    /* Debug */
    printf("%d entartete Knoten\n",L.Anz_Index);
    /**/
    /* Nun haben wir eine vollständige Liste der entarteten Knoten für
       diesen Verfeinerungslevel. Nun suchen wir die Zellen, die 
       mindestens einen dieser Knoten enthalten.
    */
    /* Für alle Zellen des aktuellen Levels */
    for(j=1;j<=MG[i-1].Cells.Anz_Cells;j++)
    {
      n=0;
      /* Für jeden Knotenindex dieser Zelle */
      for(k=0;k<NpC;k++)
      /* Suche, ob in Liste */
      {
        for(l=0;l<L.Anz_Index;l++)
        {
          if (*(MG[i-1].Cells.Vert+j*NpC+k)==L.I[l])
          {
            n=1;
            break;
          }
        }
        if(n==1) break;
      }
      if(n==1)
      {
	CL[i-1].I=(Index_Typ *) realloc(CL[i-1].I,(++(CL[i-1].Anz_Index))
					*sizeof(Index_Typ));
	CL[i-1].I[CL[i-1].Anz_Index-1]=j;
      }
    }
    /* Debug */
    printf("%d entartete Zellen\n",CL[i-1].Anz_Index);
    /**/
    Raise_EpsScale();
  }
  free(L.I);
  /* Debug */
  printf("Startup_Suchen beendet!\n");
  Reset_EpsScale();
  /**/
}

void Cleanup_Suchen(void)
{
  Zaehler_Typ i;
  for(i=0;i<Anz_Level-1;i++) Cleanup_Cells(&(MG[i].Cells));
  Cleanup_Gitter(&(MG[Anz_Level-1]));
  Anz_Level=0;
  free(MG);
  for(i=0;i<Anz_Level;i++) 
  {
    CL[i].Anz_Index=0;
    if(CL[i].I!=NULL)
    {
      free(CL[i].I);
      CL[i].I=NULL;
    }
  }
}

Index_Typ Append_Tracers(Output_Typ *O,Anz_Typ Anz)
{
   Anz_Typ Neu,i,j;
   Neu=O->Tracers.Anz_Tracer+Anz;
   for(i=0;i<DIM;i++)
     O->Tracers.K[i]=(Daten_Typ *)realloc(O->Tracers.K[i],Neu*sizeof(Daten_Typ));
   O->Tracers.Farbe=(Daten_Typ *)realloc(O->Tracers.Farbe,Neu*sizeof(Daten_Typ));
   O->Tracers.Expire=(Daten_Typ *)realloc(O->Tracers.Expire,Neu*sizeof(Daten_Typ));
#ifdef AGE
   O->Tracers.Alter=(Daten_Typ *)realloc(O->Tracers.Alter,Neu*sizeof(Daten_Typ));
#endif
#ifdef TIN
   O->Tracers.tin=(Index_Typ *)realloc(O->Tracers.tin,Neu*sizeof(Index_Typ));
#endif
   O->Tracers.InCell=(Index_Typ *)realloc(O->Tracers.InCell,Neu*sizeof(Index_Typ));
   O->Tracers.B_Coords=(B_Coord_Typ *)realloc(O->Tracers.B_Coords,Neu*sizeof(B_Coord_Typ));
   O->VecNorm=(Daten_Typ *)realloc(O->VecNorm,Neu*sizeof(Daten_Typ));
   for(i=0;i<Num_Field;i++)
     {
       O->Field[i]=(Daten_Typ *)realloc(O->Field[i],Neu*sizeof(Daten_Typ));
       for(j=0;j<Anz;j++) O->Field[i][O->Tracers.Anz_Tracer+j]=0.0;
     }
   Neu=O->Tracers.Anz_Tracer;
   O->Tracers.Anz_Tracer+=Anz;
   return Neu;
}

void Clear_Tracer(Output_Typ *O,Index_Typ Pos)
{
  Index_Typ i,j;
  O->Tracers.Anz_Tracer--;
#ifdef RTD
  fprintf(rtd_log,rtd_data,Global_Time);
  fprintf(rtd_log,",");
  fprintf(rtd_log,rtd_idx,O->Tracers.tin[Pos]);
  for(j=0;j<DIM;j++) {
    fprintf(rtd_log,",");
    fprintf(rtd_log,rtd_data,O->Tracers.K[j][Pos]);
  }
  fprintf(rtd_log,",");
  fprintf(rtd_log,rtd_data,O->Tracers.Farbe[Pos]);
  fprintf(rtd_log,",");
  fprintf(rtd_log,rtd_data,O->Tracers.Alter[Pos]);
  fprintf(rtd_log,"\n");
  fflush(rtd_log);
#endif
  for(i=Pos;i<O->Tracers.Anz_Tracer;i++)
  {
    for(j=0;j<DIM;j++)
      O->Tracers.K[j][i]=O->Tracers.K[j][i+1];
    O->Tracers.Farbe[i]=O->Tracers.Farbe[i+1];
    O->Tracers.Expire[i]=O->Tracers.Expire[i+1];
#ifdef AGE
    O->Tracers.Alter[i]=O->Tracers.Alter[i+1];
#endif
#ifdef TIN
    O->Tracers.tin[i]=O->Tracers.tin[i+1];
#endif
    O->Tracers.InCell[i]=O->Tracers.InCell[i+1];
    O->Tracers.B_Coords[i]=O->Tracers.B_Coords[i+1];
    O->VecNorm[i]=O->VecNorm[i+1];
    for(j=0;j<Num_Field;j++)
      O->Field[j][i]=O->Field[j][i+1];
  }
}

/* Maximale Anzahl von Newtonschritten */
#define MAXITER 1024
/* Genauigkeit der Interpolation */
#define GdI 1E-7

Daten_Typ Bounce(const Daten_Typ I)
{
  if(I<=0.0) return 0.0;
  if(I>=1.0) return 1.0;
  return I;
}

void Norm_Barc(Daten_Typ L[NpC])
{
  Zaehler_Typ i;
  Daten_Typ Erg;
  for(Erg=0.0,i=0;i<NpC;i++) Erg+=L[i];
  if(fabs(Erg)<=GdI) for(i=0;i<NpC;i++) L[i]=1.0/NpC;
  else for(i=0;i<NpC;i++) L[i]/=Erg;
}

#ifdef PART3D

/* Berechne die Determinante eine Spaltenweise angegebenen 3x3 Matrix */
LDaten_Typ Mat_Det(LDaten_Typ S1[3],LDaten_Typ S2[3],LDaten_Typ S3[3])
{
  return S1[0]*S2[1]*S3[2]+S2[0]*S3[1]*S1[2]+S3[0]*S1[1]*S2[2]-
    (S1[2]*S2[1]*S3[0]+S2[2]*S3[1]*S1[0]+S3[2]*S1[1]*S2[0]);
}

/* Lösen eines 3x3 Gleichungssystemes, der Rückgabewert ist im Fehlerfall
   ungleich Null; M wird Spaltenweise angegeben
*/
int Solve_Linear(LDaten_Typ M[9],Daten_Typ X[3],LDaten_Typ B[3])
{
  LDaten_Typ Det;
  Det=Mat_Det(M,M+3,M+6);
  if(fabs(Det)==0.0) return -1;
  X[0]=Mat_Det(B,M+3,M+6)/Det;
  X[1]=Mat_Det(M,B,M+6)/Det;
  X[2]=Mat_Det(M,M+3,B)/Det;
  return 0;
}

/* Berechnet die Koordinaten auf dem Einheitselement */
Zaehler_Typ Solve_Newton(LDaten_Typ F[24],Daten_Typ P[3],Daten_Typ UVW[3])
{
  Zaehler_Typ i,j;
  LDaten_Typ JG[9],G[3];
  Daten_Typ Old_UVW[3];
  /* Startlösung setzen */
  for(i=0;i<3;i++) UVW[i]=0.0;
  i=0;
  do {
    i++;
    /* Abspeichern der alten Lösung */
    for(j=0;j<3;j++) Old_UVW[j]=UVW[j];
    /* Aufbau der Jacobimatrix der Abbildung */
    JG[0]=F[3]+F[12]*UVW[1]+F[15]*UVW[2]+F[21]*UVW[1]*UVW[2];
    JG[3]=F[6]+F[12]*UVW[0]+F[18]*UVW[2]+F[21]*UVW[0]*UVW[2];
    JG[6]=F[9]+F[15]*UVW[0]+F[18]*UVW[1]+F[21]*UVW[0]*UVW[1];
    JG[1]=F[4]+F[13]*UVW[1]+F[16]*UVW[2]+F[22]*UVW[1]*UVW[2];
    JG[4]=F[7]+F[13]*UVW[0]+F[19]*UVW[2]+F[22]*UVW[0]*UVW[2];
    JG[7]=F[10]+F[16]*UVW[0]+F[19]*UVW[1]+F[22]*UVW[0]*UVW[1];
    JG[2]=F[5]+F[14]*UVW[1]+F[17]*UVW[2]+F[23]*UVW[1]*UVW[2];
    JG[5]=F[8]+F[14]*UVW[0]+F[20]*UVW[2]+F[23]*UVW[0]*UVW[2];
    JG[8]=F[11]+F[17]*UVW[0]+F[20]*UVW[1]+F[23]*UVW[0]*UVW[1];
    G[0]=F[0]+JG[0]*UVW[0]+F[6]*UVW[1]+F[9]*UVW[2]+
      F[18]*UVW[1]*UVW[2]-P[0];
    G[1]=F[1]+F[4]*UVW[0]+JG[4]*UVW[1]+F[10]*UVW[2]+
      F[16]*UVW[0]*UVW[2]-P[1];
    G[2]=F[2]+F[5]*UVW[0]+F[8]*UVW[1]+JG[8]*UVW[2]+
      F[14]*UVW[0]*UVW[1]-P[2];
    /* Berechne den Defekt */
    if(Solve_Linear(JG,UVW,G))
      {
	fprintf(stderr,"Singularitaet gefunden => Abbruch!\n");
	return -1;
      }
    /* Erzeuge die neue Lösung */
    for(j=0;j<3;j++) UVW[j]=Old_UVW[j]-UVW[j];
  } while((i<MAXITER)&&(Dist_Max(Old_UVW,UVW)>=GdI));
  if (Dist_Max(Old_UVW,UVW)>=GdI)
    {
      fprintf(stderr,"Maximale Iterationszahl ueberschritten!\n");
      return 0;
    }
  return i;
}

Index_Typ transform3d(Daten_Typ X[8],Daten_Typ Y[8],Daten_Typ Z[8],
		      Daten_Typ P[3],Daten_Typ UVW[3])
{
  LDaten_Typ F[24];
  const LDaten_Typ Q8=0.125;
  /* Ich berechne nun gewisse Formfaktoren des Elementes */
  F[ 0]=( X[0]+X[1]+X[2]+X[3]+X[4]+X[5]+X[6]+X[7])*Q8;
  F[ 1]=( Y[0]+Y[1]+Y[2]+Y[3]+Y[4]+Y[5]+Y[6]+Y[7])*Q8;
  F[ 2]=( Z[0]+Z[1]+Z[2]+Z[3]+Z[4]+Z[5]+Z[6]+Z[7])*Q8;

  F[ 3]=(-X[0]+X[1]+X[2]-X[3]-X[4]+X[5]+X[6]-X[7])*Q8;
  F[ 4]=(-Y[0]+Y[1]+Y[2]-Y[3]-Y[4]+Y[5]+Y[6]-Y[7])*Q8;
  F[ 5]=(-Z[0]+Z[1]+Z[2]-Z[3]-Z[4]+Z[5]+Z[6]-Z[7])*Q8;

  F[ 6]=(-X[0]-X[1]+X[2]+X[3]-X[4]-X[5]+X[6]+X[7])*Q8;
  F[ 7]=(-Y[0]-Y[1]+Y[2]+Y[3]-Y[4]-Y[5]+Y[6]+Y[7])*Q8;
  F[ 8]=(-Z[0]-Z[1]+Z[2]+Z[3]-Z[4]-Z[5]+Z[6]+Z[7])*Q8;  

  F[ 9]=(-X[0]-X[1]-X[2]-X[3]+X[4]+X[5]+X[6]+X[7])*Q8;
  F[10]=(-Y[0]-Y[1]-Y[2]-Y[3]+Y[4]+Y[5]+Y[6]+Y[7])*Q8;
  F[11]=(-Z[0]-Z[1]-Z[2]-Z[3]+Z[4]+Z[5]+Z[6]+Z[7])*Q8;

  F[12]=( X[0]-X[1]+X[2]-X[3]+X[4]-X[5]+X[6]-X[7])*Q8;
  F[13]=( Y[0]-Y[1]+Y[2]-Y[3]+Y[4]-Y[5]+Y[6]-Y[7])*Q8;
  F[14]=( Z[0]-Z[1]+Z[2]-Z[3]+Z[4]-Z[5]+Z[6]-Z[7])*Q8;

  F[15]=( X[0]-X[1]-X[2]+X[3]-X[4]+X[5]+X[6]-X[7])*Q8;
  F[16]=( Y[0]-Y[1]-Y[2]+Y[3]-Y[4]+Y[5]+Y[6]-Y[7])*Q8;
  F[17]=( Z[0]-Z[1]-Z[2]+Z[3]-Z[4]+Z[5]+Z[6]-Z[7])*Q8;

  F[18]=( X[0]+X[1]-X[2]-X[3]-X[4]-X[5]+X[6]+X[7])*Q8;
  F[19]=( Y[0]+Y[1]-Y[2]-Y[3]-Y[4]-Y[5]+Y[6]+Y[7])*Q8;
  F[20]=( Z[0]+Z[1]-Z[2]-Z[3]-Z[4]-Z[5]+Z[6]+Z[7])*Q8;

  F[21]=(-X[0]+X[1]-X[2]+X[3]+X[4]-X[5]+X[6]-X[7])*Q8;
  F[22]=(-Y[0]+Y[1]-Y[2]+Y[3]+Y[4]-Y[5]+Y[6]-Y[7])*Q8;
  F[23]=(-Z[0]+Z[1]-Z[2]+Z[3]+Z[4]-Z[5]+Z[6]-Z[7])*Q8;
/* Nun löse ich mit dem Newton-Verfahren nach den Koordinaten
     auf dem Einheitselement auf
  */
  return Solve_Newton(F,P,UVW); 
}

/* Berechnet die Gewichtungen L der Ecken des Einheitselementes */
void evaluate3d(Daten_Typ UVW[3],Daten_Typ L[8])
{
  Daten_Typ A,B,C,T;
  /* A,B,C sind die Einzelkomponentengewichte im ersten Punkt */
  A=(1.0-UVW[0])/2.0;
  B=(1.0-UVW[1])/2.0;
  C=(1.0-UVW[2])/2.0;
  L[0]=L[3]=L[4]=L[7]=A;
  L[1]=L[2]=L[5]=L[6]=1.0-A;
  L[0]*=B;L[1]*=B;L[4]*=B;L[5]*=B;T=1.0-B;
  L[2]*=T;L[3]*=T;L[6]*=T;L[7]*=T;
  L[0]*=C;L[1]*=C;L[2]*=C;L[3]*=C;T=1.0-C;
  L[4]*=T;L[5]*=T;L[6]*=T;L[7]*=T;
}

void Calc_Barc(Gitter_Typ *G,Tracer_Typ *T,Index_Typ P)
{
  Index_Typ C;
  Zaehler_Typ i;
  Daten_Typ Alpha,Beta,Gamma,Alph,Bet,Gam;
  Daten_Typ X[NpC],Y[NpC],Z[NpC],XYZ[DIM],UVW[DIM];
  C=T->InCell[P];
  for(i=0;i<NpC;i++)
    {
      X[i]=G->Nodes.K[0][G->Cells.Vert[i+NpC*C]];
      Y[i]=G->Nodes.K[1][G->Cells.Vert[i+NpC*C]];
      Z[i]=G->Nodes.K[2][G->Cells.Vert[i+NpC*C]];
    }
  for(i=0;i<DIM;i++) XYZ[i]=T->K[i][P];
  if (Is_Regular(C,G)==REG_REGULAR){
    XYZ[0]-=X[0];
    XYZ[1]-=Y[0];
    XYZ[2]-=Z[0];
    UVW[0]=X[1]-X[0];
    UVW[1]=Y[1]-Y[0];
    UVW[2]=Z[1]-Z[0];
    Alpha=Vec_Dot(XYZ,UVW)/Vec_Dot(UVW,UVW);
    UVW[0]=X[3]-X[0];
    UVW[1]=Y[3]-Y[0];
    UVW[2]=Z[3]-Z[0];
    Beta=Vec_Dot(XYZ,UVW)/Vec_Dot(UVW,UVW);
    UVW[0]=X[4]-X[0];
    UVW[1]=Y[4]-Y[0];
    UVW[2]=Z[4]-Z[0];
    Gamma=Vec_Dot(XYZ,UVW)/Vec_Dot(UVW,UVW);
    Alph=1.0-Alpha;
    Bet=1.0-Beta;
    Gam=1.0-Gamma;
    T->B_Coords[P].L[0]=Alph*Bet*Gam;
    T->B_Coords[P].L[1]=Alpha*Bet*Gam;
    T->B_Coords[P].L[2]=Alpha*Beta*Gam;
    T->B_Coords[P].L[3]=Alph*Beta*Gam;
    T->B_Coords[P].L[4]=Alph*Bet*Gamma;
    T->B_Coords[P].L[5]=Alpha*Bet*Gamma;
    T->B_Coords[P].L[6]=Alpha*Beta*Gamma;
    T->B_Coords[P].L[7]=Alph*Beta*Gamma;
  } else {
    if(transform3d(X,Y,Z,XYZ,UVW)<0)
      for(i=0;i<DIM;i++) UVW[i]=0.0;
    evaluate3d(UVW,T->B_Coords[P].L);
  }
  Norm_Barc(T->B_Coords[P].L);
}
#else

Zaehler_Typ newton_solver(LDaten_Typ A[2], LDaten_Typ B[2], LDaten_Typ C[2], 
		  LDaten_Typ D[2], Daten_Typ XY[2], Daten_Typ UV[2])
{
  Zaehler_Typ i,j;
  LDaten_Typ F[2],JF[2][2],det;
  Daten_Typ Old_UV[2];
  /* Initialisierung */
  UV[0]=UV[1]=0.0;
  i=0;
  do {
    i++;
    for(j=0;j<2;j++) Old_UV[j]=UV[j];
    for(j=0;j<2;j++)
      {
	F[j]=A[j]*UV[0]+B[j]*UV[1]+C[j]*UV[0]*UV[1]+D[j]-XY[j];
	JF[j][0]=C[j]*UV[1]+A[j];
	JF[j][1]=C[j]*UV[0]+B[j];
      }
    det=JF[0][0]*JF[1][1]-JF[1][0]*JF[0][1];
    if (fabs(det)==0.0)
      {
	fprintf(stderr,"Singularitaet gefunden => Abbruch!\n");
	return -1;
      }
    det=1.0/det;
    for(j=0;j<2;j++) 
      UV[j]-=det*(F[j]*JF[1-j][1-j]-F[1-j]*JF[j][1-j]);
      
  } while((i< MAXITER) && (Dist_Max(Old_UV,UV)>=GdI));
  if (Dist_Max(Old_UV,UV)>=GdI)
    {
      fprintf(stderr,"Maximale Iterationszahl ueberschritten!\n");
      return 0;
    }
  return i;
}

/* K[0][] enthaelt die X- und K[1][] die Y-Koordinaten der Eckpunkte */
Zaehler_Typ transform2d(Daten_Typ X[4],Daten_Typ Y[4],
			Daten_Typ XY[2],Daten_Typ UV[2])
{
  LDaten_Typ A[2],B[2],C[2],D[2],DJ[4];
  Zaehler_Typ i;
  DJ[0]=-X[0]-X[1]+X[2]+X[3];
  DJ[1]= X[0]-X[1]+X[2]-X[3];
  DJ[2]=-Y[0]+Y[1]-Y[2]+Y[3];
  DJ[3]=-Y[0]+Y[1]+Y[2]-Y[3];
  for(i=0;i<4;i++) DJ[i]*=0.5;
  A[0]=X[1]-X[0]+DJ[1];
  A[1]=DJ[3];
  B[0]=DJ[0];
  B[1]=Y[2]-Y[0]-DJ[3];
  C[0]=DJ[1];
  C[1]=-DJ[2];
  D[0]=X[0]+X[1]+DJ[0];
  D[1]=Y[0]+Y[2]+DJ[2];
  for(i=0;i<2;i++)
    {
      A[i]*=0.5;
      B[i]*=0.5;
      C[i]*=0.5;
      D[i]*=0.5;
    }
  i=newton_solver(A,B,C,D,XY,UV);
  return i;
}

void evaluate2d(Daten_Typ UV[2],Daten_Typ L[4])
{
  Daten_Typ Alpha,Beta;
  Alpha=(1.0-UV[0])/2.0;
  Beta=(1.0-UV[1])/2.0;
  L[0]=Alpha*Beta;
  L[1]=(1.0-Alpha)*Beta;
  L[2]=(1.0-Alpha)*(1.0-Beta);
  L[3]=Alpha*(1-Beta);
}

void Calc_Barc(Gitter_Typ *G,Tracer_Typ *T,Index_Typ P)
{
  Index_Typ C;
  Zaehler_Typ i;
  Daten_Typ Alpha,Beta,Alph,Bet,X[4],Y[4],XY[2],UV[2];
  C=T->InCell[P];
  for(i=0;i<4;i++)
  {
    X[i]=G->Nodes.K[0][G->Cells.Vert[i+4*C]];
    Y[i]=G->Nodes.K[1][G->Cells.Vert[i+4*C]];
  }
  for(i=0;i<2;i++) XY[i]=T->K[i][P];
  if (Is_Regular(C,G)==REG_REGULAR)
  {
    XY[0]-=X[0];
    XY[1]-=Y[0];
    UV[0]=X[1]-X[0];
    UV[1]=Y[1]-Y[0];
    Alpha=Vec_Dot(XY,UV)/Vec_Dot(UV,UV);
    UV[0]=X[3]-X[0];
    UV[1]=Y[3]-Y[0];
    Beta=Vec_Dot(XY,UV)/Vec_Dot(UV,UV);
    Alph=1.0-Alpha;
    Bet=1.0-Beta;
    T->B_Coords[P].L[0]=Alph*Bet;
    T->B_Coords[P].L[1]=Alpha*Bet;
    T->B_Coords[P].L[2]=Alpha*Beta;
    T->B_Coords[P].L[3]=Alph*Beta;  
  }
  else
  {
    if(transform2d(X,Y,XY,UV)<0)
      for(i=0;i<2;i++) UV[i]=0.0;
    evaluate2d(UV,T->B_Coords[P].L);
  }
  Norm_Barc(T->B_Coords[P].L);
}
#endif

void Gen_Tracers(Output_Typ *O,Block_Typ *B)
{
  Index_Typ Offset,OldCell;
  Zaehler_Typ i,j,l,Lev[DIM+1];
  Daten_Typ H[DIM],P[DIM];
  for(i=0,j=1;i<DIM;i++)
    j*=B->N[i];
  l=j;
  Offset=Append_Tracers(O,l);
  for(i=0;i<DIM;i++)
    H[i]=(B->P[1][i]-B->P[0][i])/B->N[i];
 
  Lev[0]=1;
  for(i=1;i<=DIM;i++) Lev[i]=Lev[i-1]*B->N[i-1];
  OldCell=0;
  for(i=0;i<l;i++)
    {
      for(j=0;j<DIM;j++)
	O->Tracers.K[j][Offset]=B->P[0][j]+(0.5+((i%Lev[j+1])/Lev[j]))*H[j];
      O->Tracers.Farbe[Offset]=B->Farbe;
      O->Tracers.Expire[Offset]=B->Haltbarkeit;
#ifdef AGE
      O->Tracers.Alter[Offset]=0.0;
#endif
      for(j=0;j<DIM;j++)P[j]=O->Tracers.K[j][Offset];
#ifdef MOVBC
  /* Teste zuerst, ob das Partikel im bewegten Rand ist */
      if(Test_MovBC(P)) OldCell=O->Tracers.InCell[Offset]=0;
      else 
#endif
	OldCell=O->Tracers.InCell[Offset]=Suche_Zelle(P,OldCell);
	if(!O->Tracers.InCell[Offset]) Clear_Tracer(O,Offset);
	else
	  {
	    Calc_Barc(MG+Anz_Level-1,&(O->Tracers),Offset);
#ifdef TIN
	    O->Tracers.tin[Offset]=New_tin++;
#endif
	    Offset++;
	  }
    }
}

void Bewege_Tracer(Output_Typ *O,Input_Typ *I,Daten_Typ t)
{
  Zaehler_Typ i,j,k,pos,stop;
  Daten_Typ Dummy,P[DIM];
  stop=O->Tracers.Anz_Tracer;
  for(i=0,pos=0;i<stop;i++)
  {
    for(k=0;k<DIM;k++)
      {
	for(j=0,Dummy=0.0;j<NpC;j++)
	  Dummy+=O->Tracers.B_Coords[pos].L[j]*
	    I->V[k][I->Gitter.Cells.Vert[j+NpC*O->Tracers.InCell[pos]]];
	P[k]=O->Tracers.K[k][pos]+=t*Dummy;
      }
#ifdef MOVBC
  /* Teste zuerst, ob das Partikel im bewegten Rand ist */
    if(Test_MovBC(P)) O->Tracers.InCell[pos]=0;
    else {
#endif
      O->Tracers.InCell[pos]=Suche_Zelle(P,O->Tracers.InCell[pos]);
#ifdef MOVBC
    }
#endif
    O->Tracers.Expire[pos]-=t;
#ifdef AGE
    O->Tracers.Alter[pos]+=t;
#endif
    if((!O->Tracers.InCell[pos])||(O->Tracers.Expire[pos]<-EPS)) Clear_Tracer(O,pos);
    else {
      Calc_Barc(&(I->Gitter),&(O->Tracers),pos);
      pos++;
    }
  }
}

void Calculate_Values(Output_Typ *O,Input_Typ *I)
{
  Zaehler_Typ i,j,l;
  Daten_Typ Dummy,Vel[DIM];
  /* Zuerst berechne ich die Norm der Geschwindigkeiten */
  for (i=0;i<O->Tracers.Anz_Tracer;i++) {
    for(l=0;l<DIM;l++) {
      Vel[l]=0.0;
      for(j=0;j<NpC;j++) {
	Vel[l]+=O->Tracers.B_Coords[i].L[j]*
	  I->V[l][I->Gitter.Cells.Vert[j+NpC*O->Tracers.InCell[i]]];
      }
    }
    O->VecNorm[i]=Norm_Euklid(Vel);
  }
  /* Nun fuelle ich die anderen Datenfelder */
  for(l=0;l<Num_Field;l++)
    {
      for(i=0;i<O->Tracers.Anz_Tracer;i++)
	{
	  for(j=0,Dummy=0.0;j<NpC;j++)
	    Dummy+=O->Tracers.B_Coords[i].L[j]*
	      I->Field[l][I->Gitter.Cells.Vert[j+NpC*O->Tracers.InCell[i]]];
	  O->Field[l][i]=Dummy;
	}
    }
}

const Daten_Typ DELTA=0.5E-2; 
void Schreibe_Output(const char *Name,Output_Typ *O,Gitter_Typ *F)
{
  FILE *f;
  Anz_Typ A,k;

#ifdef MOVBC 
  Anz_Typ i,j,l;
  Index_Typ n;
#ifdef PART3D
#define ANZ_POLY 6 
 const Index_Typ Index[ANZ_POLY][4]=
 {
   {0,1,5,4},
   {3,0,4,7},
   {2,3,7,6},
   {1,2,6,5},
   {4,7,6,5},
   {0,1,2,3}
 };
#else
#define ANZ_POLY 1 
 const Index_Typ Index[ANZ_POLY][4]=
 {
   {0,1,2,3}
 };
#endif
 Daten_Typ K[DIM],P[4];
#endif
  f=fopen(Name,"w");
  Schreibe_Vertreter(f);
#ifdef MOVBC
  /* Test (Polygone) */
  fprintf(f,"polygons\n");
  A=F->Cells.Anz_Cells;
  for(i=1;i<=A;i++) {
    n=1;
    for(k=0;k<NpC;k++) {
      for(j=0;j<DIM;j++) K[j]=F->Nodes.K[j][F->Cells.Vert[i*NpC+k]];
      if(!Test_MovBC(K)) { n=0; break; }
    }
    if(n) {
      for(j=0;j<ANZ_POLY;j++) {
	fprintf(f,"2 4\n");
	for(k=0;k<DIM;k++) {
	  for(l=0;l<4;l++) P[l]=F->Nodes.K[k][F->Cells.Vert[i*NpC+Index[j][l]]];
	  Schreibe_Daten(f,4,P);
	}
	if(DIM<3) Schreibe_Daten(f,4,NULL);
      }
    }
  }
  fprintf(f,"endpoly\n");
  /**/
#endif
  fprintf(f,"tracers ");
  A=O->Tracers.Anz_Tracer;
  Schreibe_Anz(f,A);
  /* Schreibe die Tracercoordinaten */
  for(k=0;k<DIM;k++) Schreibe_Daten(f,A,O->Tracers.K[k]);
  if(DIM<3) Schreibe_Daten(f,A,NULL);
  /* Schreibe die Norm der Geschwindigkeiten */
  fprintf(f,"VecNorm\n");
  Schreibe_Daten(f,A,O->VecNorm);
  /* Schreibe die variablen Datenfelder */
  for(k=0;k<Num_Field;k++)
  {
    fprintf(f,"%s\n",Tracer_Names[k]);
    Schreibe_Daten(f,A,O->Field[k]);
  }
  /* Die Felder "Farbe" und Alter fallen aus der Rolle, da sie nicht aus den
     Simulationsdateien interpoliert werden und deshalb gesondert
     behandelt werden muessen.
  */
  /* Schreibe Zusatzdatenfelder */
  fprintf(f,"Farbe\n");
  Schreibe_Daten(f,A,O->Tracers.Farbe);
#ifdef AGE
  fprintf(f,"Alter\n");
  Schreibe_Daten(f,A,O->Tracers.Alter);
#endif
  fprintf(f,"endtrace\n");
#ifdef TIN
  fprintf(f,"traceids\n");
  Schreibe_Index(f,A,O->Tracers.tin);
#endif
  fprintf(f,"probtime ");
  Schreibe_Daten(f,1,&Global_Time);
  fprintf(f,"endgmv\n");
  fclose(f);
}

void Cleanup_Output(Output_Typ *O)
{
  Zaehler_Typ i;
  for(i=0;i<DIM;i++) free(O->Tracers.K[i]);
  free(O->Tracers.Farbe);
  free(O->Tracers.Expire);
#ifdef AGE
  free(O->Tracers.Alter);
#endif
#ifdef TIN
  free(O->Tracers.tin);
#endif
  free(O->Tracers.InCell);
  free(O->Tracers.B_Coords);
  free(O->VecNorm);
  O->Tracers.Anz_Tracer=0;
  for(i=0;i<Num_Field;i++) free(O->Field[i]);
}

Zaehler_Typ Modul_Input,Modul_Output;
Anz_Typ Anz_Zeitschritte,Anz_Block;
Daten_Typ Zeitschritt;
Block_Typ *Blocks;

void Lese_Parameter(const char *Name)
{
  FILE *f;
  Zaehler_Typ i,j;
  f=fopen(Name,"r");
  if(f==NULL) File_Error(Name);
  Lese_Anz(f,&Anz_Zeitschritte);Lese_Daten(f,1,&Zeitschritt);
  Lese_Anz(f,&Modul_Input);Lese_Anz(f,&Modul_Output);
  Lese_Anz(f,&Anz_Block);
  Blocks=(Block_Typ *)calloc(Anz_Block,sizeof(Block_Typ));
  for(i=0;i<Anz_Block;i++)
  {
    for(j=0;j<2;j++)
    {
      Lese_Daten(f,DIM,Blocks[i].P[j]);
    }
    Lese_Index(f,DIM,Blocks[i].N);
    Lese_Anz(f,&(Blocks[i].Modul));
    Lese_Anz(f,&(Blocks[i].Offset));
    Lese_Anz(f,&(Blocks[i].Start));
    Lese_Anz(f,&(Blocks[i].Stop));
    Lese_Daten(f,1,&(Blocks[i].Farbe));
    Lese_Daten(f,1,&(Blocks[i].Haltbarkeit));
  }
}

#define CONFIGFILE "#data/gmvpt.dat" 
#define COARSE "#gmv/coarse.gmv" 
#define GROBGITTER "#gmv/tria1.gmv" 

int main(int argc,char *argv[])
{
  Zaehler_Typ In_Pos,Out_Pos,T_Pos,i;
  char * Praefix;
  if (argc>1) Praefix=argv[1];
  else Praefix="t";
#ifdef PART3D
  printf("GMVPT-3D");
#else
  printf("GMVPT-2D");
#endif
  printf(" Version 1.30");
#ifdef MOVBC
  printf(" mit bewegtem Rand,");
#endif
#ifdef AGE
  printf(" mit Herausschreiben des Alters der Partikel,");
#endif
#ifdef TIN
  printf(" mit Tracer Indentification Numbers,");
#endif
#ifdef RTD
  printf(" mit RTD-Logdatei,");
#endif
#ifdef LDOUBLE
  printf(" mit erhöhter Genauigkeit,");
#endif
  printf(" angepasst fuer ");
#ifdef BOUSS
  printf("BOUSS");
#else
#ifdef PART3D
  printf("PP3D,CC3D,CP3D");
#else
  printf("PP2D,CC2D,CP2D");
#endif
#endif
#ifdef POWERLAW
  printf(" (POWERLAW)");
#endif
  puts("...");
  if (!Datei_Existiert(CONFIGFILE)) File_Error(CONFIGFILE);
  if (!Datei_Existiert(GROBGITTER)) File_Error(GROBGITTER);
  if (!Datei_Existiert(COARSE)) File_Error(COARSE);
  Lese_Parameter(CONFIGFILE);
  In_Pos=0;
  Out_Pos=0;
  for(i=0;i<DIM;i++) O.Tracers.K[i]=NULL;
  for(i=0;i<Num_Field;i++) O.Field[i]=NULL;
  O.Tracers.Anz_Tracer=0;
  O.Tracers.Farbe=NULL;
  O.Tracers.Expire=NULL;
#ifdef AGE
  O.Tracers.Alter=NULL;
#endif
#ifdef TIN
  O.Tracers.tin=NULL;
#endif
  O.Tracers.InCell=NULL;
  O.Tracers.B_Coords=NULL;
  Startup_Suchen();
  Startup_Input(NULL,&I);
  /* Debugging */
  Tracer_Verworfen_Inc=Tracer_Gerettet_Inc=0;
  /**/
#ifdef RTD
  sprintf(Nrtd,"%s_rtd.csv",Praefix);
  if(!(rtd_log=fopen(Nrtd,"w"))) {
    fprintf(stderr,"Datei %s konnte nicht zum Schreiben geöffnet werden!\n",Nrtd);
    exit(EXIT_FAILURE);
  }
#endif
  for(T_Pos=0;T_Pos<=Anz_Zeitschritte;T_Pos++) {
    if(((T_Pos+Modul_Input-1)%Modul_Input==0)&&(T_Pos>0)){
      do {
	In_Pos++;
	sprintf(NIn,"#gmv/u.%d.gmv",In_Pos);
      } while(!Datei_Existiert(NIn));
      Update_Input(NIn,&I);
#ifdef MOVBC
      Input_Time=Global_Time;
#endif
    }
    if(T_Pos>0) {
      printf("Bewege Tracer um einen Zeitschritt...\n");
      Bewege_Tracer(&O,&I,Zeitschritt);
    }
    for(i=0;i<Anz_Block;i++){
      if((T_Pos>=Blocks[i].Start)&&(T_Pos<=Blocks[i].Stop)&&
	 ((T_Pos+(Blocks[i].Modul-Blocks[i].Offset))%Blocks[i].Modul==0)){
	printf("Erzeuge Block %d ...\n",i+1); 
	Gen_Tracers(&O,Blocks+i);
      }
    }
    Append_Tracers(&O,0);
    if(T_Pos%Modul_Output==0){
      printf("Ausgabe der Tracer!\n");
      sprintf(NOut,"#gmv/%s.%d.gmv",Praefix,Out_Pos++);
      Calculate_Values(&O,&I);
      Schreibe_Output(NOut,&O,&(I.Gitter));
    }
    /* Debugging */
    printf("Das Programm versuchte %u Tracer zu verwerfen.\n",
	   Tracer_Verworfen_Inc);
    printf("Es konnten aber %u Tracer gerettet werden.\n",
	   Tracer_Gerettet_Inc);
    Tracer_Verworfen+=Tracer_Verworfen_Inc;
    Tracer_Gerettet+=Tracer_Gerettet_Inc;
    Tracer_Verworfen_Inc=Tracer_Gerettet_Inc=0;
    /**/
    Global_Time+=Zeitschritt;
  }
    /* Debugging */
    printf("Das Programm versuchte insgesammt %u Tracer zu verwerfen.\n",
	   Tracer_Verworfen);
    printf("Es konnten aber insgesammt %u Tracer gerettet werden.\n",
	   Tracer_Gerettet);
    /**/
#ifdef RTD
    fclose(rtd_log);
#endif  
  Cleanup_Input(&I);
  Cleanup_Output(&O);
  Cleanup_Suchen();
  return EXIT_SUCCESS;
}
