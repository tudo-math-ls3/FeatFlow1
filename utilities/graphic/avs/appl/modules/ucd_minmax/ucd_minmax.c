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
/* ucd_minmax
 * override the data min/max values for
 * user supplied constants.
 *
 * Author:    Ian J. Curington, AVS Inc., UK
 *
 * Revision:
 *  27 May  92 ianc - Original
 *   4 June 92 ianc - added mesh_id check and override
 * Revision:
 *  14 April 97 cbek - added Disable
 */

/*-----------------------------------------------------*
 *                                                     *
 *        ****  ucd minmax module  ****                *
 *                                                     *
 * given a range of scalar values provided by the      *
 * the ucd input,         create a new ucd structure   *
 * which consists of the same ucd values and structure *
 * but with the minimun and maximum tagged to be       *
 * within a user-specified range.                      *
 *                                                     *
 * This module is a derivative of ucd_thresh.c         *
 *                                                     *
 *-----------------------------------------------------*/
/*
  Copyright (c) 1989 by
  Stardent Computer Inc.
  Copyright (c) 1992 by
  Advanced Visual Systems Inc.
  All Rights Reserved
  
  This software comprises unpublished confidential information of
  Stardent Computer Inc. and may not be used, copied or made
  available to anyone, except in accordance with the license
  under which it is furnished.
  
  This file is under sccs control at Stardent in:
  /user/ianc/sccs_4.0/avs/modules/ucd/filters/s.ucd_minmax.c
  
  */

/*********************************************/
/*  Include Definitions                      */
/*********************************************/

#include <stdio.h>
#include <avs/avs.h>
#include <avs/ucd_defs.h>
#include <avs/udata.h>

/*********************************************/
/*  Global Data Structures                   */
/*********************************************/

/* define seperate module formation, remove if mongo */
#define sep_exe 1

/* magic static number so mesh is assumed static */
#define FROZEN_MESH  -42    

extern char *AVSstatic;

typedef struct _ucd_stats {
	int *nd_comp_list, node_num_comp, node_offset;
	int *cl_comp_list, cell_num_comp;
} Ucd_Stats_Type;


/*  count the number of mid-edge nodes for a cell.  */

#define  count_me_nodes(FLAG) {\
	       if ((FLAG)) \
	       for (i = 0; i < 24; i++) \
	       if ((FLAG) & (0x1 << i)) node_csize++;\
       }

/*-----------------------------------------------------*
 *                                                     *
 *   ****  ucd_minmax Module Compute Func  ****        *
 *                                                     *
 *-----------------------------------------------------*/

ucd_minmax (ucd_input, ucd_output, range_min, range_max, freeze_mesh, disable ) 
    
    UCD_structure *ucd_input, **ucd_output;
    float *range_min, *range_max;
    int freeze_mesh,disable;
{
	
	char model_name[80], string[UCD_LABEL_LEN], data_type[80], ucell_name[80], 
	data_labels[UCD_LABEL_LEN], delim[2];
	
	float r, g, b, x, y, z, *elem_data, *disp, *cptr, value, minv, maxv,
	xmin, xmax, ymin, ymax, zmin, zmax, *cmap_ptr, surf_color[3], 
	scale, min_extent[3], max_extent[3], vcolors[12][3], iso_level,
	hi_val, lo_val, val, *min_node_data, *max_node_data, *node_data, *cell_data, 
	*new_cell_data,	*xc, *yc, *zc;
	
	int cell, node, new_cell, node_id, elem_id, i, j, found,  cell_type, me_flags, 
	*node_list, int_cell_name, mat_id, util_flag, num_cell_nodes, data_len, 
	cell_len, node_len, num_new_cells, out_nodes, num_nodes, num_cells, 
	*cell_list, cell_tsize, node_csize, ucd_flags, num_vcomp, *nd_comp_list, *cl_comp_list,
	cell_num_comp, node_num_comp, node_offset, ncells, *clist, old_cell_tsize, old_node_csize;
	
	int    cell_offset, new_offset, comp, veclen;
	float  *min_cell_data, *max_cell_data;
        int    mesh_id;

	Ucd_Stats_Type *ucd_stats;
	
	/**************
	 ***  body  ***
	 **************/
	

	if (!UCDstructure_get_header (ucd_input, model_name, &data_len, &ucd_flags,
				      &num_cells, &cell_len, &num_nodes, &node_len,
				      &util_flag)) {
		
		AVSerror ("Error in ucd_minmaxold: can't get header.\n");
		return (0);
	}
	
	if (node_len) { 
		if (!UCDstructure_get_node_data (ucd_input, &node_data)) {
			AVSerror ("Error in ucd_minmax: can't get node data.\n"); 
			return (0);
		}
	}
	else {
		AVSerror (" Error in ucd_minmax:  no nodal data.\n");
		return (0);
	}
	
	
	if (AVSinput_changed("Input", 0)) {
		UCDstructure_get_node_components (ucd_input, &nd_comp_list);
		
		if (node_len) {
			UCDstructure_get_node_labels (ucd_input, data_labels, delim);
			for (i = 0, node_num_comp = 1; data_labels[i] != '\0'; i++)
				if (data_labels[i] == delim[0]) node_num_comp++;
		}
		
		if (cell_len) {
			UCDstructure_get_cell_labels (ucd_input, data_labels, delim);
			for (i = 0, cell_num_comp = 1; data_labels[i] != '\0'; i++)
			if (data_labels[i] == delim[0]) cell_num_comp++;
		}
		if (AVSstatic == NULL) 
			ucd_stats = (Ucd_Stats_Type *)malloc(sizeof(Ucd_Stats_Type));
		else {
			ucd_stats = (Ucd_Stats_Type *)AVSstatic;
		} 
		
		ucd_stats->nd_comp_list = nd_comp_list;
		ucd_stats->node_num_comp = node_num_comp;
		ucd_stats->cl_comp_list = cl_comp_list;
		ucd_stats->cell_num_comp = cell_num_comp;
		
		AVSstatic = (char *)ucd_stats;
	}
	else {
		ucd_stats = (Ucd_Stats_Type *)AVSstatic;
		nd_comp_list = ucd_stats->nd_comp_list;
		node_num_comp = ucd_stats->node_num_comp;
		cl_comp_list = ucd_stats->cl_comp_list;
		cell_num_comp = ucd_stats->cell_num_comp;
	}
	
	cell_list = (int *)malloc(sizeof(int) * num_cells);
	num_new_cells = 0;
        node_offset = 0;
	
	/* determine if each cell is within the given range.  */
	gen_cell_list (ucd_input, num_cells, node_data, 
		       &num_new_cells, cell_list, &node_csize);
	
	
	/*  create a new ucd data set.  */
	util_flag = 0;
	cell_tsize = node_csize;
	
	UCDstructure_get_sizes (ucd_input, &old_cell_tsize, &old_node_csize);
	if (old_node_csize == 0)
		node_csize = 0;
	else
		node_csize = cell_tsize;
	if (num_new_cells == 0) {
		data_len = 0;
		cell_tsize = 0;
		cell_len = 0;
		num_nodes = 0;
		node_csize = 0;
		node_len = 0;
	}
	if (*ucd_output)
		UCDstructure_free (*ucd_output);
		
	*ucd_output = (UCD_structure *)UCDstructure_alloc (model_name, data_len, 
							   ucd_flags, num_new_cells, 
							   cell_tsize, cell_len,
							   num_nodes, node_csize, 
							   node_len, util_flag);
	if (num_new_cells) {
		
		UCDstructure_get_extent (ucd_input, min_extent, max_extent);
		UCDstructure_set_extent (*ucd_output, min_extent, max_extent);
		
		UCDstructure_get_node_positions (ucd_input, &xc, &yc, &zc);
		UCDstructure_set_node_positions (*ucd_output, xc, yc, zc);
		
		if (node_len) {
			UCDstructure_get_node_data (ucd_input, &node_data);
			UCDstructure_set_node_data (*ucd_output, node_data); 
			
			UCDstructure_get_node_components (ucd_input, &nd_comp_list);
			UCDstructure_set_node_components (*ucd_output, nd_comp_list, node_num_comp);
			
			UCDstructure_get_node_labels (ucd_input, data_labels, delim);
			UCDstructure_set_node_labels (*ucd_output, data_labels, delim);
			
			UCDstructure_get_node_units (ucd_input, data_labels, delim);
			UCDstructure_set_node_units (*ucd_output, data_labels, delim);

			min_node_data = (float *)malloc(sizeof(float) * node_len);
			max_node_data = (float *)malloc(sizeof(float) * node_len);

			/* Override the input data values here for nodes */
			UCDstructure_get_node_minmax (ucd_input, min_node_data, max_node_data);

                        if (disable==0)
				for ( j=0; j<node_len; j++)
				{
					min_node_data[j] = *range_min;
					max_node_data[j] = *range_max;
				}

			UCDstructure_set_node_minmax (*ucd_output, min_node_data, max_node_data);


/* IAC CODE CHANGE : 			free (min_node_data); */
			 free(min_node_data);

/* IAC CODE CHANGE : 			free (max_node_data); */
			 free(max_node_data);
		}

		if (cell_len) {
			UCDstructure_get_cell_components (ucd_input, &cl_comp_list);
			UCDstructure_set_cell_components (*ucd_output, cl_comp_list, cell_num_comp);

			UCDstructure_get_cell_data (ucd_input, &cell_data);
			
			new_cell_data = (float *)malloc(cell_len * num_new_cells *
							sizeof(float));
			UCDstructure_set_cell_data (*ucd_output, new_cell_data); 

			cell_offset = 0;
			new_offset = 0;
			for (comp = 0; comp < cell_num_comp; comp++) {
				veclen = cl_comp_list[comp];
				for (new_cell = 0; new_cell < num_new_cells; new_cell++) {
					cell = cell_list[new_cell];
					for (i = 0; i < veclen; i++) {
						new_cell_data[new_offset + 
							      new_cell * veclen + i] =
								      cell_data[cell_offset + cell * veclen + i];
					}
				}
				cell_offset += veclen * num_cells;
				new_offset += veclen * num_new_cells;
			}
		
			UCDstructure_get_cell_labels (ucd_input, data_labels, delim);
			UCDstructure_set_cell_labels (*ucd_output, data_labels, delim);
			
			UCDstructure_get_cell_units (ucd_input, data_labels, delim);
			UCDstructure_set_cell_units (*ucd_output, data_labels, delim);
			
			min_cell_data = (float *)malloc(sizeof(float) * cell_len);
			max_cell_data = (float *)malloc(sizeof(float) * cell_len);

			/* Override the input data values here for cells */
			UCDstructure_get_cell_minmax (ucd_input, min_cell_data, max_cell_data);

			if (disable==0)
				for ( j=0; j<node_len; j++)
				{
					min_cell_data[j] = *range_min;
					max_cell_data[j] = *range_max;
				}
			UCDstructure_set_cell_minmax (*ucd_output, min_cell_data, max_cell_data);


/* IAC CODE CHANGE : 			free (min_cell_data); */
			 free(min_cell_data);

/* IAC CODE CHANGE : 			free (max_cell_data); */
			 free(max_cell_data);

/* IAC CODE CHANGE : 			free (new_cell_data); */
			 free(new_cell_data);
		}

		if (data_len) {
			UCDstructure_get_data (ucd_input, &cell_data);
			UCDstructure_set_data (*ucd_output, cell_data); 
			
			UCDstructure_get_data_labels (ucd_input, data_labels, delim);
			UCDstructure_set_data_labels (*ucd_output, data_labels, delim);
			
			UCDstructure_get_data_units (ucd_input, data_labels, delim);
			UCDstructure_set_data_units (*ucd_output, data_labels, delim);
		}
		/*  set node id (assume node id is int).  */
		
		UTILucd_copy_cell_list(ucd_input, *ucd_output, num_new_cells, cell_list);
		
		UTILucd_copy_nodes(ucd_input, *ucd_output, 0);

		/*  this module does not modify the mesh topology or vertex list, */
		/*  so pass through the mesh_id flag so downstream modules can    */
		/*  cache the mesh. Force it to be a static value if requested.   */
		if ( freeze_mesh )
			mesh_id = FROZEN_MESH;
		else
			UCDstructure_get_mesh_id (ucd_input, &mesh_id);
		UCDstructure_set_mesh_id (*ucd_output, mesh_id);


		/*  if cell connectivity was set in the original
		    ucd struct then compute the connectivity.     */
		
		if (old_node_csize) {
			if (!UCDstructure_set_cell_connect(*ucd_output)) {
				AVSerror (" Error in read_ucd: node connectivity size too small.\n");
				return (0);
			}
		}
	}
	if (cell_list) 

/* IAC CODE CHANGE : 		free (cell_list); */
		 free(cell_list);
	
	return (1);
}


/*-----------------------------------------------------*
 *                                                     *
 *           ****  gen_cell_list  ****                 *
 *                                                     *
 *-----------------------------------------------------*/

gen_cell_list (ucd_input, num_cells, node_data,
	       num_new_cells, cell_list, 
	       pnode_csize)
    
    float *node_data;
    
    int num_cells, 
	    *num_new_cells, *cell_list, *pnode_csize;
    
    UCD_structure *ucd_input; 
{
	
	float val;
	
	int i, j, elem_id, int_cell_name, mat_id, cell_type, me_flags, *node_list,
	in, out, num_cell_nodes, node, new_cells, node_csize;
	
	/**************
	 ***  body  ***
	 **************/
	
	new_cells = 0;
	node_csize = 0;
	
	for (i = 0; i < num_cells; i++) {
		UCDcell_get_information (ucd_input, i, &elem_id, &int_cell_name, &mat_id, 
					 &cell_type, &me_flags, &node_list); 
		
		num_cell_nodes = UCD_num_nodes[cell_type];
		
		cell_list[new_cells++] = i;
		count_me_nodes(me_flags);
		node_csize += num_cell_nodes;
	}
	
	*num_new_cells = new_cells;
	*pnode_csize = node_csize;
}


/*-----------------------------------------------------*
 *                                                     *
 *           ****  ucd_minmax_init  ****               *
 *                                                     *
 *-----------------------------------------------------*/

ucd_minmax_init()
{
	AVSstatic = (char *)0;
}


/*-----------------------------------------------------*
 *                                                     *
 *           ****  ucd_minmax_finis  ****              *
 *                                                     *
 *-----------------------------------------------------*/

ucd_minmax_finis()
{
	Ucd_Stats_Type *ucd_stats;
	
	if (AVSstatic == NULL) return;
	
	ucd_stats = (Ucd_Stats_Type *)AVSstatic;
	

/* IAC CODE CHANGE : 	free (ucd_stats); */
	 free(ucd_stats);
}


/*-----------------------------------------------------*
 *                                                     *
 *           ****  ucd_minmax_desc  ****               *
 *                                                     *
 *-----------------------------------------------------*/

ucd_minmax_desc()
{
	char *string, tmp[20], *init, *sep;
	
	int ucd_minmax(), param;
	
	static char *choices = "<data 1>.<data 2>.<data 3>.<data 4>.<data 5>";
	
	/**************
	 ***  body  ***
	 **************/
	
	AVSset_module_name ("UCD minmax c", MODULE_FILTER);
	
	AVScreate_input_port ("Input", "ucd", REQUIRED);
	
	AVScreate_output_port ("Output", "ucd"); 
	
	param=AVSadd_float_parameter ("Min", 0., FLOAT_UNBOUND, FLOAT_UNBOUND);
	
	param=AVSadd_float_parameter ("Max", 1., FLOAT_UNBOUND, FLOAT_UNBOUND);

	param=AVSadd_parameter ("Freeze Mesh", "boolean", 0, 0, 1 );

	param=AVSadd_parameter ("Disable", "boolean", 0, 0, 1 );
	
	AVSset_init_proc (ucd_minmax_init);
	
	AVSset_destroy_proc (ucd_minmax_finis);
	
	AVSset_compute_proc (ucd_minmax);
	
#ifndef sep_exe
	AVSset_module_flags (REENTRANT|COOPERATIVE);
#endif
}

#ifdef sep_exe
AVSinit_modules()
{
	int ucd_minmax_desc();
	
	AVSmodule_from_desc(ucd_minmax_desc);
}
#endif
