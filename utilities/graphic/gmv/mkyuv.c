#include <stdlib.h>
#include <stdio.h>

#define MYBUF 4000

int xres,yres,start,stop,add,i,adaptive;
char S[MYBUF];
int main(int argc,char **argv)
{
   if((argc<7)||(argc>8))
     {
        printf("Aufruf: %s <xres> <yres> <attribut> <prefix> \
<start> <stop> [<add>]\nWenn <add> nicht explizit angegeben \
wird gilt <add>:=<start>\n<add>:=0 bei adaptiver Zeitschrittweite\n",argv[0]);
	return EXIT_SUCCESS;
     }
   xres=atoi(argv[1]);
   yres=atoi(argv[2]);
   start=atoi(argv[5]);
   stop=atoi(argv[6]);
   adaptive=0;
   if(argc==7)add=start;
   else
     {
       add=atoi(argv[7]);
       if(add==0)
	 {
	   adaptive=1;
	   add=1;
	 }
     }
   for(i=1;start<=stop;start+=add)
     {
        sprintf(S,"test -r u.%d.gmv",start);
	if(system(S))
	  {
	     if(!adaptive)printf("Datei u.%d.gmv nicht auffindbar!\n",start);
	     continue;
	  } 
	sprintf(S,"gmv -m -a %s -w 0 0 %d %d -i u.%d.gmv -s\n",argv[3],xres,yres,start);
	printf(S);
	system(S);
	sprintf(S,"sgitopnm AzsnapgmvAz | ppmtoyuvsplit %s%d\n",argv[4],i);
	printf(S);
	system(S);
	printf("Frame %d wurde erstellt!\n",i);
	i++;
     }
   remove("AzsnapgmvAz");
   return EXIT_SUCCESS;
}



