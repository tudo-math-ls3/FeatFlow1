#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <avs/flow.h>
#include <avs/avs.h>
#include <avs/field.h>

#define cb "\
  panel \"$Module\"  -p ui -xy 684,36 -wh 273,397\n\
  manip \"$Module:File name\" -w browser -p \"$Module\" -xy 13,9 -wh 248,194\n\
  manip \"$Module:Animation On\/Off\" -w toggle -p \"$Module\" -xy 27,220 -wh 118,22\n\
  manip \"$Module:Single Step On\/Off\" -w toggle -p \"$Module\" -xy 27,200 -wh 118,22\n\
  manip \"$Module:Nr. frames\" -w idial -p \"$Module\" -xy 25,254 -wh 90,130\n\
  manip \"$Module:Current frame\" -w idial -p \"$Module\" -xy 154,254 -wh 90,130\n\
"

#define FORMATS ".anim"
/* ..inp.1 */

static int  anim=1,cur_anim=1,singlestep=1,sync=0,maxstep,step;
static char file[100];
static unsigned nr_digits=0,first_step=1;
static char line[200],string[20],name[200],*end_step;
static char format[40];

int  read_seq_desc()
{

  int  param;

  AVSset_module_name("ReadAnimMax",MODULE_DATA);

  AVScreate_output_port("File","string");
  param=AVSadd_parameter("File name","string",0,0,FORMATS);
  AVSconnect_widget(param,"browser");

  param=AVSadd_parameter("Animation On/Off","boolean",1,0,1);
  param=AVSadd_parameter("Single Step On/Off","boolean",1,0,1);
  param=AVSadd_parameter("Nr. files","integer",256,1,256);
  AVSconnect_widget(param,"idial");
  param=AVSadd_parameter("Current file","integer",1,1,256);
  AVSconnect_widget(param,"idial");
}



char  *get_seq(file)
char  *file;
{
   char  *ptr;
   int   i=0;


   ptr=file+(strlen(file)-1);

   while((*ptr !='.'))
       if((*ptr=='/')|| (ptr==file))
	   return NULL;
       else ptr--;
   end_step=ptr+1;
   ptr--;

   for(i=0; (*ptr !='.') && (isdigit(*ptr)) ; i++,--ptr)
	 if(ptr==file) return NULL;
   if(*ptr=='.'){ /*OK*/
     first_step=atoi(++ptr);
     nr_digits=i;
     return ptr;
   }else
      return NULL ;

}

void get_num_files(ptr,steps)
char *ptr;
int  *steps;
{
  FILE *file;
  int count=first_step;

  *steps=1;
  while(count<maxstep) {
    sprintf(ptr,"%d.",count);
    strcat(line,string);
    
    /*    fprintf(stderr,"ReadAnim:name:%s \n",line);*/
    file=fopen(line,"r");
    if(file) { 
      /*      fprintf(stderr,"ReadAnim: open file:%s \n",line);*/
      fclose(file);
      (*steps)++;
    }
    count++;
  }
}

void get_dir(file)
char *file;
{
  char *ptr;
  ptr=file+(strlen(file)-1);
  
  while((*ptr !='/'))
    if((ptr==file))
      break;
    else 
      ptr--;
  ptr++;
  *ptr=0;
}


int main(argc,argv)
int argc;
char *argv[];

{  
    int  read_seq_desc();
    char *get_seq();
    int  c_step,curstep=1;
    char *file_name;
    int  nrframes,i,reset;
    char *ob,*eb ;
    char  *ptr,*ptr1;
    int nrsteps=0;
    FILE *FilePtr;


    AVScorout_init(argc,argv,read_seq_desc);
    AVScorout_set_sync(0);
/*    AVScommand("kernel",cb,&ob,&eb);*/
    cur_anim=1;
    reset=0;

    while(1){
       if(cur_anim) AVScorout_wait();

       AVScorout_input(&file_name,&anim,&singlestep,&nrframes,&c_step);

       nrsteps=nrframes;
       cur_anim=anim ;
       if(!singlestep)
	 cur_anim=0;


       if(AVSparameter_changed("Nr. files"))
         AVSmodify_parameter("Current file",AVS_VALUE|AVS_MINVAL|AVS_MAXVAL,
			      first_step,first_step,(maxstep==1)? maxstep+1:maxstep);

       if(AVSparameter_changed("Current file")) {
	 if(curstep>c_step) {
	   reset=1;
	 }
	 curstep=step=c_step;
       }

       if(AVSparameter_changed("File name")){
	 FilePtr=fopen(file_name,"r");
	 strcpy(line,file_name);
	 get_dir(line);
	 fscanf(FilePtr,"%s",file_name);
	 fscanf(FilePtr,"%d",&maxstep);
	 fclose(FilePtr);
	 printf(stderr,"ReadAnim: filename %s maxstep %d\n",file_name,maxstep);
         strcat(line,file_name);

         if((ptr=get_seq(line))==NULL){
            AVSmessage("ReadAnim :",AVS_Warning,AVSmodule,"ReadAnim",
	              "OK","Can't detect the seq\n");
            AVScorout_wait() ;
	 } else {
	   strcpy(string,end_step);
	   *ptr='\0';
	   strcpy(name,line);
	   
	   get_num_files(ptr,&nrsteps);
	   fprintf(stderr,"ReadAnim: nrsteps %d\n",nrsteps);
	   AVSmodify_parameter("Nr. files",AVS_MINVAL | AVS_MAXVAL | AVS_VALUE,
			       nrsteps,1,(nrsteps==1)? nrsteps+1:nrsteps);

           AVSmodify_parameter("Current file",AVS_MINVAL | AVS_MAXVAL | AVS_VALUE,
			       first_step,first_step,(maxstep==1)? maxstep+1:maxstep);
	   curstep=step=first_step;
	   reset=1;
	 }
       }

       if(reset==1) {
	 reset=0;
	 AVScorout_mark_changed();
	 AVSmark_output_unchanged("File");
	 AVScorout_output(line);
	 AVScorout_wait();
       }

       if(cur_anim || !file_name ) continue ;


       if(!singlestep) {
	 while(step<=maxstep) {
	   sprintf(ptr,"%d.",step);
	   strcat(line,string);
	  
	   FilePtr=fopen(line,"r");
	   if(FilePtr) { 
	     /*fprintf(stderr,"ReadAnim: cant open file:%s \n",line);*/
	     fclose(FilePtr);
	     step++;
	     break;
	   } else
	     step++;
	 }

/*	 AVScorout_exec();*/
	 AVScorout_output(line);
	 
	 AVSmodify_parameter("Current file",AVS_VALUE|AVS_MINVAL|AVS_MAXVAL,
			       step-1,first_step,(maxstep==1)? maxstep+1: maxstep);
	 curstep++;
	 /*fprintf(stderr,"Single Step: nr_steps:%d file:%s step:%d\n",nrsteps,line,curstep);*/
       } else {
	 while(step<=maxstep) {
	   sprintf(ptr,"%d.",step);
	   strcat(line,string);
	   
	   FilePtr=fopen(line,"r");
	   if(FilePtr) { 
	     /*fprintf(stderr,"ReadAnim: open file:%s \n",line);*/
	     fclose(FilePtr);
	     step++;
	     
	     AVScorout_exec();
	     AVScorout_mark_changed();
	     AVScorout_output(line);
	     AVScorout_exec();
	     
	     AVSmodify_parameter("Current file",AVS_VALUE|AVS_MINVAL|AVS_MAXVAL,
				 step-1,first_step,(maxstep==1)? maxstep+1: maxstep);
	     
	     /*fprintf(stderr,"Anim nrsteps:%d file:%s step:%d\n",nrsteps,line,curstep);*/
	     
	     AVScorout_wait();
	     AVScorout_input(&file_name,&anim,&singlestep,&nrframes,&c_step);
	     if(anim)
	       break;
	   } else 
	     step++;	   
	 }
       }
       AVSmodify_parameter("Single Step On/Off",AVS_VALUE,1,0,1);
       AVSmodify_parameter("Animation On/Off",AVS_VALUE,1,0,1);
       cur_anim=1;
      
  }
}

