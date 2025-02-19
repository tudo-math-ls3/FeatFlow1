#include <stdlib.h>
#include <stdio.h>

#define DEFAULT_BITRATE 5000000.0
#define MYBUF 5000

char s[MYBUF];
int i,number,xres,yres;
long int max_bits;
int main(int argc,char **argv)
{
   if((argc<5)||(argc>6))
     {
	printf("Aufruf: %s <xres> <yres> <prefix> <frames> [<max_MB>]\n",argv[0]);
	return EXIT_SUCCESS;
     }
   xres=atoi(argv[1]);
   yres=atoi(argv[2]);
   number=atoi(argv[4]);
   if(argc==5) max_bits=(DEFAULT_BITRATE*number)/25;
   else max_bits=atof(argv[5])*1024*1024*8;
   puts("Starte MPEG Encoder...");
   sprintf(s,"mpeg -PF -p 3 -a 1 -b %d -h %d -v %d -x %ld -s %s.mpeg %s >/dev/null",number,xres,yres,max_bits,argv[3],argv[3]);
   puts(s);
   system(s);
   return EXIT_SUCCESS;
}
