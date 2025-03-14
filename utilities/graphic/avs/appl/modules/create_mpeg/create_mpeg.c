/****************************************************************************
                  INTERNATIONAL AVS CENTER
	(This disclaimer must remain at the top of all files)

WARRANTY DISCLAIMER

This module and the files associated with it are distributed free of charge.
It is placed in the public domain and permission is granted for anyone to use,
duplicate, modify, and redistribute it unless otherwise noted.  Some modules
may be copyrighted.  You agree to abide by the conditions also included in
the AVS Licensing Agreement, version 1.0, located in the main module
directory located at the International AVS Center ftp site and to include
the AVS Licensing Agreement when you distribute any files downloaded from 
that site.

The International AVS Center, MCNC, the AVS Consortium and the individual
submitting the module and files associated with said module provide absolutely
NO WARRANTY OF ANY KIND with respect to this software.  The entire risk as to
the quality and performance of this software is with the user.  IN NO EVENT
WILL The International AVS Center, MCNC, the AVS Consortium and the individual
submitting the module and files associated with said module BE LIABLE TO
ANYONE FOR ANY DAMAGES ARISING FROM THE USE OF THIS SOFTWARE, INCLUDING,
WITHOUT LIMITATION, DAMAGES RESULTING FROM LOST DATA OR LOST PROFITS, OR ANY
SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES.

This AVS module and associated files are public domain software unless
otherwise noted.  Permission is hereby granted to do whatever you like with
it, subject to the conditions that may exist in copyrighted materials. Should
you wish to make a contribution toward the improvement, modification, or
general performance of this module, please send us your comments:  why you
liked or disliked it, how you use it, and most important, how it helps your
work. We will receive your comments at avs@ncsc.org.

Please send AVS module bug reports to avs@ncsc.org.

******************************************************************************/


/*

 AUTHOR: Alex Knowles

 For further information please contact:

  The Computer Graphics Unit
  Manchester Computing Centre
  The University of Manchester
  Manchester M13 9PL
  United Kingdom

  tel: +44 61 275 6095
  fax: +44 61 275 6040

  email: cgu-info@mcc.ac.uk
*/

/*__________________________________________________________________________
Alex Knowles
Manchester Computing Centre, Computer Graphics Unit	   Tel. 061 275 6095
E-Mail	ark@hpa.cgu.mcc.ac.uk
        alex@ed.ac.uk
--------------------------------------------------------------------------*/

/* a module to create mpeg movies of images passed to it. They must be
   of consistant size. Small images look good. Large Images are slow and
   jerky on replay.
*/
/* version 1.1
   April 1994 small bug fixes to fix memory allocation
   May 1994 added destroy routine toget rid of temp files
   May 1994 make it work with an SGI
*/

/* use AVS's memory stuff */
#define MEM_DEFS_ENABLE 1

#include <stdio.h>
#include <string.h>
#include <avs/avs.h>
#include <avs/port.h>
#include <avs/field.h>
 
/* IAC CODE CHANGE
#ifdef  _INCLUDE_POSIX_SOURCE
#include <unistd.h>
#include <symlink.h>
#else
#define  _INCLUDE_POSIX_SOURCE
#include <unistd.h>
#include <symlink.h>
#undef _INCLUDE_POSIX_SOURCE
#endif
   IAC CODE CHANGE */

#define TRUE 1
#define FALSE 0
#define MIN(a,b)	((a)>(b) ? (b) : (a))	/* minimum of a and b        */
#define MAX(a,b)	((a)>(b) ? (a) : (b))	/* maximum of a and b        */
#define CORRECT(x)	MIN( MAX((x),0), 255)	/* limits x to range 0->255  */
#define PIXELS(x)	(width*height)		/* num of pixels in an image */
#define BYTES(x)	(PIXELS(x)*4)		/* num of bytes in an image  */
#define EVEN(x)		(((x)%2==1) ? (x)-1:(x))/* round down to an even num */
/* checks to see if `a` points to something and isn`t the empty string */
#define SOMETHING(a)	(a && strcmp(a,""))

/* button names os that they are easy to change */
#define FILENAME "filename"
#define FRAME "Frame"
#define SIZE "Size"

typedef unsigned char byte;

/* path max is the max length of a path from the file browser */
/* the number 1024 was from the browser.h include file this is the*/
/* maximun filename length from the file browser */
#ifndef PATH_MAX
#define PATH_MAX 1024
#endif

char	mpegfile[PATH_MAX];
int	mpegcount=0;
int	width, height;

/* function prototypes */
void make_mpeg();
void do_mpeg();
void clear_mpegtemp();
void repeat_mpeg();

/* *****************************************/
/*  Module Description                     */
/* *****************************************/
int Create_MPEG_desc()
{
  int in_port, out_port, param, iresult;
  extern int Create_MPEG_compute();
  extern int Create_MPEG_destroy();

  AVSset_module_name("Create MPEG", MODULE_RENDER);
  
  /* Input Port Specifications               */
  in_port = AVScreate_input_port("input", 
				 "field 2D 4-vector uniform byte", REQUIRED);
  
  /* Parameter Specifications                */
  param = AVSadd_parameter(FILENAME, "string", "", "", "*");
  AVSconnect_widget(param, "browser");
  AVSadd_parameter_prop(param,"height","integer",10);

  param = AVSadd_parameter("Pause", "boolean", 1, 0, 1);
  AVSconnect_widget(param, "toggle");

  param = AVSadd_parameter("Forget it", "oneshot", 0, 0, 1);
  AVSconnect_widget(param, "oneshot");

  param = AVSadd_parameter("Make it", "oneshot", 0, 0, 1);
  AVSconnect_widget(param, "oneshot");

  param = AVSadd_parameter("Repeat", "oneshot", 0, 0, 1);
  AVSconnect_widget(param, "oneshot");

  param = AVSadd_parameter(SIZE,"string","Width : 000 Height : 000", "", "");
  AVSconnect_widget(param, "text");
  AVSadd_parameter_prop(param,"width","integer",4);

  param = AVSadd_parameter(FRAME, "string", "Frame : 000", "", "");
  AVSconnect_widget(param, "text");
  AVSadd_parameter_prop(param,"width","integer",4);
  
  AVSset_compute_proc(Create_MPEG_compute);
  AVSset_destroy_proc(Create_MPEG_destroy);
  return(1);
}

int
Create_MPEG_destroy ()
{
  if( mpegcount>0 ) 	
    clear_mpegtemp();	/* erase all the temp files */
}
 
/* *****************************************/
/* Module Compute Routine                  */
/* *****************************************/
int
Create_MPEG_compute(input, filename, pause, forget, make, repeat, size, frame)
     AVSfield_char *input;	/* the input image to be added to the movie */
     char	*filename;	/* filename of mpeg movie.. with directory  */
     int 	pause;		/* pausing means add no more images	    */
     int 	forget;		/* cancel this movie and erase temp files   */
     int 	make;		/* shall i make it   and erase temp files   */
     int 	repeat;		/* repeat the last frame ???		    */
     char	*size;		/* string with the size of the image in it  */
     char	*frame;		/* string with the frame number in it	    */
{
  char line[128]; 		/* the frame count is built up in this	    */

  /* make it so that you can't turn off pause unless a file is chosen */
  if( AVSparameter_changed("Pause") && !SOMETHING(filename) ) {
    AVSwarning("You can't start a movie without specifying a file.\n");
    AVSmodify_parameter("Pause",AVS_VALUE,1,0,1);
    return(1);
  }
  
  /* make it so that you cannot change the filename halfway through a movie */
  if( AVSparameter_changed(FILENAME) ){
    if( mpegcount==0 ){
      strcpy(mpegfile,filename);
    } else {
      AVSwarning("You can't change the filename while\n"
		 "you're making a movie!!!\n"
		 "Either finish this one or start a new one\n");
      /* now make filename what it should be      */
      AVSmodify_parameter(FILENAME,AVS_VALUE,mpegfile,"","");
    }
  }
  
  if( AVSinput_changed("input",0) && !pause && SOMETHING(filename) )
    /* if new input & we're not paused & filename is something */
    do_mpeg(input);			/* add this frame to the sequence   */
  else
    if(mpegcount>0)
      if( repeat && !pause )
	repeat_mpeg();			/* make links to the last frame     */
      else{
	AVSmodify_parameter("Pause",AVS_VALUE,1,0,1);	/* turn pause on    */
	if( make ){
	  make_mpeg();				/* make the mpeg movie	    */
	  clear_mpegtemp();			/* erase all the temp files */
	} else
	  if( forget ){
	    clear_mpegtemp();	/* erase all the temp files */
	    putchar('\n');
	  }
      }

  sprintf(line, "Width : %03d   Height : %03d",
	  EVEN( MAXX( input ) ),  EVEN( MAXY( input ) ) );
  AVSmodify_parameter(SIZE,AVS_VALUE,line,"","");

  sprintf( line, "Frame : %03d",mpegcount);
  AVSmodify_parameter(FRAME,AVS_VALUE,line,"","");

  return(1);
}

/* ***********************************************************************/
/* Initialization for modules contained in this file.                    */
/* ***********************************************************************/
static int ((*mod_list[])()) = {  Create_MPEG_desc  };
#define NMODS (sizeof(mod_list) / sizeof(char *))
     
AVSinit_modules()
{
  AVSinit_from_module_list(mod_list, NMODS);
}
 
/* ***********************************************************************/
/* subroutines by ARK... subroutines 'r' us!!! 				 */
/* ***********************************************************************/
  
/* gives Y luminance value in range 0->255	*/
#define Yval(ptr)	(0.2989*((float)*(ptr+1))+	/* red   */	\
			 0.5866*((float)*(ptr+2))+	/* green */	\
			 0.1144*((float)*(ptr+3)) )	/* blue  */

/* gives U -- blue chroma value in range -127->128	*/
#define Uval(ptr)       (0.493*(0.8856*((float)*(ptr+3))-   /* blue  */	\
				0.5866*((float)*(ptr+2))-   /* green */	\
				0.2989*((float)*(ptr+1)) ) )/* red   */

/* gives V -- red chroma value in range -127->128	*/
#define Vval(ptr)	(0.877*(0.7011*((float)*(ptr+1))-   /* red   */	\
				0.5866*((float)*(ptr+2))-   /* green */	\
				0.1144*((float)*(ptr+3)) ) )/* blue  */

void do_mpeg( image )
     AVSfield_char* image;			/* the image to encode	    */
{
  FILE	*out;		/* the file to write. only 1 file open at any time  */
  int	x,y,
	err,			/* error value returned by file functions   */
	bufptr;				/* index into the buffer array	    */
  byte	*buffer,		/* array to store different images in	    */
	*ptr=NULL;				/* pointer into an image    */
  float	Y,U,V;						/* YUV values	    */
  char 	fname[PATH_MAX+10];	/* array to build up temp filename in	    */
  byte	ERROR=FALSE;

  AVSmodule_status("Adding Frame",0);

  if( mpegcount==0 ){
    /* This is the first image */
    width =EVEN( MAXX(image) );
    height=EVEN( MAXY(image) );
    printf("First Image recieved by Create MPEG dimensions %d %d\n",
	   width, height );
  } else {
    /* check that the image is the right size */
    if( width!=EVEN(MAXX(image)) || height!=EVEN(MAXY(image)) ){
      /* the size of the input has changed */
      AVSwarning("The size of the input image has been changed.\n"
		 "This image has been ignored.\n"
		 "Either return size of input image to %03d %03d,\n"
		 "or forget this movie and start a new one\n",width,height);
      return;
    }
  }

  buffer=(byte *) calloc( PIXELS(image), sizeof( byte ) );

  if( !buffer ){
    AVSerror("Failed to allocate memory for YUV buffer\n");
    return;
  }
      
  printf("Doing MPEG..for frame %d", mpegcount );

/* Y luminescence stored as one byte per pixel */
  /* make ptr point to alpha at top of image nad fill up buffer */
  ptr=I2DV(image,0,0);
  for( y=0 ; y < height ; y++ ){
    ptr=I2DV(image,0,y);
    for( x=0 ; x < width ; x++ ){
      buffer[x+y*width]=(byte) CORRECT(Yval(ptr));
      ptr+=4;
    }
  }
  /* WRITE the Y file */
  sprintf(fname,"%s_temp%d.Y",mpegfile,mpegcount);
  out=fopen(fname,"w");
  if( !out ){
    AVSerror("Could not open temporary MPEG file %s",fname);
    ERROR=TRUE;
  } else {
    err=fwrite(buffer , sizeof(byte) , PIXELS(image) , out);
    if( err != PIXELS(image) ){
      AVSerror("Failed to write whole Y file %s\nOnly wrote %d bytes",
	       fname,err);
      ERROR=TRUE;
    }
    fclose(out);
  }

  /* free the memory used by buffer that is no longer needed */
  realloc( buffer, PIXELS(image)/4 );

  AVSmodule_status("Adding Frame",33);
/* U blue chroma stored as one byte per 4 pixels */
  /* fill up buffer with U */
  bufptr=0;
  for( y=0 ; y < height ; y+=2 )
    for( x=0 ; x < width ; x+=2 ){
      ptr=I2DV( image, x  , y  );	U =Uval(ptr);	/* top left	*/
      ptr+=4;				U+=Uval(ptr);	/* top right	*/
      ptr=I2DV( image, x  , y+1  );	U+=Uval(ptr);	/* bottom left  */
      ptr+=4;				U+=Uval(ptr);	/* bottom right */
      buffer[bufptr++]=(byte) CORRECT(U/4.0+128.0);
    }

  /* save the U file */
  sprintf(fname,"%s_temp%d.U",mpegfile,mpegcount);
  out=fopen(fname,"w");
  if( !out ){
    AVSerror("Could not open temporary MPEG file %s",fname);
    ERROR=TRUE;
  } else {
    err=fwrite( buffer , sizeof(byte) , PIXELS(image)/4 , out);
    if( err != PIXELS(image)/4 ){
      AVSerror("Failed to write whole U file %s\nOnly wrote %d bytes",
	       fname,err);
      ERROR=TRUE;
    }
    fclose(out);
  }

  AVSmodule_status("Adding Frame",66);
/* V red chroma stored as one byte per 4 pixels */
  /* fill up buffer with V */
  bufptr=0;
  for( y=0 ; y < height ; y+=2 )
    for( x=0 ; x < width ; x+=2 ){
      ptr=I2DV( image, x  , y  );	V =Vval(ptr);	/* top left	*/
      ptr+=4;				V+=Vval(ptr);	/* top right	*/
      ptr=I2DV( image, x  , y+1  );	V+=Vval(ptr);	/* bottom left  */
      ptr+=4;				V+=Vval(ptr);	/* bottom right */
      buffer[bufptr++]=(byte) CORRECT(V/4.0+128.0);
    }

  /* save the V file */
  sprintf(fname,"%s_temp%d.V",mpegfile,mpegcount);
  out=fopen(fname,"w");
  if( !out ){
    AVSerror("Could not open temporary MPEG file %s",fname);
    ERROR=TRUE;
  } else {
    err=fwrite( buffer , sizeof(byte) , PIXELS(image)/4 , out );
    if( err != PIXELS(image)/4 ){
      AVSerror("Failed to write whole V file %s\nOnly wrote %d bytes",
	       fname,err);
      ERROR=TRUE;
    }
    fclose(out);
  }

  if( !ERROR )  mpegcount++;	/* only increse frame count if no errors    */

/* IAC CODE CHANGE :   free( buffer ); */
   free(buffer );
  printf("...Done   \r"); fflush(stdout);
}

/* repeat_mpeg makes a symbolic link to the last frame written and
 * increments mpegcount.
 */
void repeat_mpeg()
{
  char 	fname[PATH_MAX+10],
       	oldfname[PATH_MAX+10];
  int 	err=0;

  printf("Repeating MPEG..for frame %d", mpegcount );

  if( mpegcount==0 ){
    AVSwarning("Whoops. This is the first image in the sequence\n");
  } else {
    AVSmodule_status("Repeating Frame",0);
    sprintf(oldfname,"%s_temp%d.Y",mpegfile,mpegcount-1);
    sprintf(fname,   "%s_temp%d.Y",mpegfile,mpegcount);
    err+=symlink( oldfname, fname );

    AVSmodule_status("Repeating Frame",33);
    sprintf(oldfname,"%s_temp%d.U",mpegfile,mpegcount-1);
    sprintf(fname,   "%s_temp%d.U",mpegfile,mpegcount);
    err+=symlink( oldfname, fname );

    AVSmodule_status("Repeating Frame",66);
    sprintf(oldfname,"%s_temp%d.V",mpegfile,mpegcount-1);
    sprintf(fname,   "%s_temp%d.V",mpegfile,mpegcount);
    err+=symlink( oldfname, fname );

    if( err ){
      AVSerror("Failed when symlink-ing duplicate frame files in repeat_mpeg\nThis frame has been ignored\n");
      return;
    } else 
      mpegcount++;
  }
  printf("...Done   \r");	fflush(stdout);
}

/* takes all the mpeg temp files and makes a movie and then removes
   all the temp files
*/
void make_mpeg()
{
  char  lineout[PATH_MAX+PATH_MAX+PATH_MAX+1024];
  
  AVSmodule_status("Creating MPEG",50);
  sprintf(lineout,"mpeg %s -a 0 -b %d -h %d -v %d -s %s %s_temp >%s_message",
	  ( width%16==0 && height%16==0 ) ? "" : "-PF",
	  mpegcount-1, width, height, mpegfile, mpegfile, mpegfile );
  printf("\nEXECUTING : %s\n",lineout);
  /* system returns non null on error so... */
  if( system( lineout ) ){
    AVSerror("Something went wrong with MPEG movie.\n"
	     "look in %s_message for help",mpegfile);
  }
}

/* clear_mpegtemp -- removes all the temp files */
void clear_mpegtemp()
{
  int x,err;
  char fname[PATH_MAX+10];

  for( x=mpegcount-1 ; x>=0 ; x-- ){
    if( x%10==0 ){
      AVSmodule_status("Erasing Temp files",
		       (int) ((1.0-(float)((float)x/(float)mpegcount))*100.0));
    }
    sprintf(fname,"%s_temp%d.Y",mpegfile,x); err =remove(fname);
    sprintf(fname,"%s_temp%d.U",mpegfile,x); err+=remove(fname);
    sprintf(fname,"%s_temp%d.V",mpegfile,x); err+=remove(fname);
    if( err!=0 )
      AVSerror( "Failed when remove temporary files near\n%s\n",fname );
  }
  mpegcount=0;
  putchar('\n');
}
