/*   =S+H= 		Copyright (c) 1992 by
     =S+H= 		Advanced Visual Systems Inc.
     =S+H= 		All Rights Reserved
     =S+H= 
     =S+H= This software comprises unpublished confidential information of
     =S+H= Advanced Visual Systems Inc. and may not be used, copied or made
     =S+H= available to anyone, except in accordance with the license
     =S+H= under which it is furnished.
     =S+H= 
     =S+H= */
#  ifdef __STDC__ 						/* =S+H= */
#  ident "@(#)ucd_extract.c	8.1 AVS 92/10/11" 						/* =S+H= */
#  endif 							/* =S+H= */
#  ifdef SCCS 							/* =S+H= */
   static char sccsid[]="@(#)ucd_extract.c	8.1 AVS 92/10/11"; 				/* =S+H= */
#  endif 							/* =S+H= */
/*-----------------------------------------------------*
 *                                                     *
 *        ****  ucd extract module  ****               *
 *                                                     *
 *                                                     *
 * extract a single data component from a ucd data set.*
 * this only works for node data.                      *
 *-----------------------------------------------------*/

#include <stdio.h>
#include <avs/avs.h>
#include <avs/geom.h>
#include <avs/ucd_defs.h>
#include <avs/udata.h>

extern char *AVSstatic;

typedef struct _ucd_stats {
	float extent[6];
	int *nd_comp_list, nd_veclen, num_comp, offset;
} Ucd_Stats_Type;


#ifdef titan
#pragma OPT_LEVEL 0
#endif
ucd_extract (ucd_input, ucd_output, dummy, data_type) 
    
    char *data_type;
    
    int dummy;
    
    UCD_structure *ucd_input, **ucd_output;
{
	
	char model_name[80], string[40], ud_cell_type[80], data_labels[UCD_LABEL_LEN], 
	delim[2], old_data_type[UCD_LABEL_LEN];
	
	float r, g, b, x, y, z, *elem_data, *disp, *cptr, value, minv, maxv,
	xmin, xmax, ymin, ymax, zmin, zmax, 
	scale, min_extent[3], max_extent[3],
	*min_node_data, *max_node_data, *new_node_data;
	
	int cell, node, cell_id, i, j, found, comp_num, cell_type, me_flags, 
	*node_list, mat_id, util_flag, data_len, cell_len, node_len, ncells,
	num_new_cells, num_new_comp, cell_tsize, node_csize, node_id,
	input_node_len;
	
	float extent[6], vmax, vmin, *node_data, *xc, *yc, *zc;
	
	int num_nodes, num_cells, num_clist, *cell_list, ucd_flags, mesh_id;
	
	int num_vcomp, *nd_comp_list, num_comp, offset;
	
	Ucd_Stats_Type *ucd_stats;
	
	/**************
	 ***  body  ***
	 **************/
	
	if (!UCDstructure_get_header (ucd_input, model_name, &data_len, &ucd_flags,
				      &num_cells, &cell_len, &num_nodes, &node_len,
				      &util_flag)) {
		
		AVSerror ("Error in ucd_extractold: can't get header.\n");
		return (0);
	}
	
	input_node_len = node_len;
	if (!UCDstructure_get_node_data (ucd_input, &node_data)) {
		AVSerror ("Error in ucd_extract: can't get node data.\n"); 
		return (0);
	}
	
	UCDstructure_get_node_components (ucd_input, &nd_comp_list);
	
	
	if (AVSstatic == NULL) {
		ucd_stats = (Ucd_Stats_Type *)calloc(1,sizeof(Ucd_Stats_Type));
		ucd_stats->nd_comp_list = NULL;
		AVSstatic = (char *)ucd_stats;
	}
	else {
		ucd_stats = (Ucd_Stats_Type *)AVSstatic;
	} 
	
	if (AVSinput_changed("Input", 0)) {
		UCDstructure_get_node_labels (ucd_input, data_labels, delim);
		
		for(num_comp = 0,i = 0; i < node_len; num_comp++) 
			i += nd_comp_list[num_comp];
		
		/*  display a menu of the defined data types.  */
		strcpy(old_data_type, data_type);
		ucd_define_menu (ucd_input, node_len, cell_len, 1, &num_comp, &num_vcomp,
				 &nd_comp_list, 0, string);
		ucd_get_offset (ucd_input, num_nodes, num_comp, nd_comp_list, string, 
				&ucd_stats->offset);

		ucd_stats->num_comp = num_comp;
		ucd_stats->nd_comp_list = nd_comp_list;
		
		if (ucd_stats->num_comp > 1) {
			for (i = 0, found = 0; i < num_comp && !found; i++) {
				UCDstructure_get_node_label (ucd_input, i, string);
				found = !strcmp (string, old_data_type);
			}
			if (found) {
				AVSmodify_parameter ("node data",  AVS_VALUE, string,
						     data_labels, ".");
				ucd_get_offset (ucd_input, num_nodes, num_comp, nd_comp_list, string, 
						&ucd_stats->offset);
			}
			else {
				AVSmodify_parameter ("node data",  AVS_VALUE, "",
						     data_labels, ".");
				return (0);
			}
		}
	}
	else {
		ucd_stats = (Ucd_Stats_Type *)AVSstatic;
		nd_comp_list = ucd_stats->nd_comp_list;
		num_comp = ucd_stats->num_comp;
		ucd_get_offset (ucd_input, num_nodes, num_comp, nd_comp_list, data_type, 
				&ucd_stats->offset);
	}
	
	
	if (*ucd_output)
		UCDstructure_free (*ucd_output);
	/*  extract those cells.  */
	
	/*  get the component number.  */
	
	if (ucd_stats->num_comp == 1) {
		UCDstructure_get_node_label (ucd_input, 0, string);
		comp_num = 0;
	}
	else {
		for (i = 0, found = 0; i < num_comp && !found; i++) {
			UCDstructure_get_node_label (ucd_input, i, string);
			found = !strcmp (string, data_type);
		}
		comp_num = i - 1;
	}
	node_len = nd_comp_list[comp_num];
	util_flag = 0;
	
	
	/*  copy extracted node data.  */
	
	new_node_data = NULL;
	new_node_data = (float *)malloc(sizeof(float) * num_nodes * node_len);
	if (new_node_data == NULL) {
		AVSerror ("Error in ucd_extract: can't alloc node data.\n"); 
		return (0);
	}
	
	for (i = ucd_stats->offset, j = 0; i < node_len * num_nodes + 
	     ucd_stats->offset; i++, j++) 
		new_node_data[j] = node_data[i];
	
	UCDstructure_get_sizes (ucd_input, &cell_tsize, &node_csize);
	
	*ucd_output = (UCD_structure *)UCDstructure_alloc (model_name, data_len, 
							   ucd_flags, num_cells, cell_tsize, cell_len,
							   num_nodes, node_csize, node_len, util_flag);
	
	UCDstructure_get_extent (ucd_input, min_extent, max_extent);
	UCDstructure_set_extent (*ucd_output, min_extent, max_extent);
	
	UCDstructure_get_node_positions (ucd_input, &xc, &yc, &zc);
	UCDstructure_set_node_positions (*ucd_output, xc, yc, zc);
	
	UCDstructure_get_mesh_id (ucd_input, &mesh_id);
	UCDstructure_set_mesh_id (*ucd_output, mesh_id);

	UCDstructure_set_node_data (*ucd_output, new_node_data); 
	
	
	/*  copy component information.  */
	
	num_new_comp = 1;
	nd_comp_list[0] = node_len;
	
	UCDstructure_set_node_components (*ucd_output, nd_comp_list, num_new_comp);
	UCDstructure_set_node_labels (*ucd_output, string, ".");
	
	min_node_data = NULL;
	max_node_data = NULL;
	min_node_data = (float *)malloc(sizeof(float) * input_node_len);
	max_node_data = (float *)malloc(sizeof(float) * input_node_len);
	
	UCDstructure_get_node_minmax (ucd_input, min_node_data, max_node_data);
	
	max_node_data[0] = max_node_data[comp_num];
	min_node_data[0] = min_node_data[comp_num];
	
	UCDstructure_set_node_minmax (*ucd_output, min_node_data, max_node_data);
	

	UTILucd_copy_cells(ucd_input, *ucd_output);

	UTILucd_copy_nodes(ucd_input, *ucd_output, 1);

	if (new_node_data)
		free (new_node_data);
	if (min_node_data)
		free (min_node_data);
	if (max_node_data)
		free (max_node_data);
	return(1);
}


/*-----------------------------------------------------*
 *                                                     *
 *         ****  ucd_extract_init  ****                *
 *                                                     *
 *-----------------------------------------------------*/

ucd_extract_init()
{
	AVSstatic = (char *)0;
}


/*-----------------------------------------------------*
 *                                                     *
 *         ****  ucd_extract_finis  ****               *
 *                                                     *
 *-----------------------------------------------------*/

ucd_extract_finis()
{
	Ucd_Stats_Type *ucd_stats;
	
	if (AVSstatic == NULL) return;
	
	ucd_stats = (Ucd_Stats_Type *)AVSstatic;
	
	free (ucd_stats);
}


/*-----------------------------------------------------*
 *                                                     *
 *           ****  ucd_extract_desc  ****              *
 *                                                     *
 *-----------------------------------------------------*/

ucd_extract_desc()
{
	char *string, tmp[20], *init, *sep;
	
	int ucd_extract(), param;
	
	static char *choices = "<data 1>.<data 2>.<data 3>.<data 4>.<data 5>";
	
	/**************
	 ***  body  ***
	 **************/
	
	AVSset_module_name ("UCD extract c", MODULE_FILTER);
	
	AVScreate_input_port ("Input", "ucd", REQUIRED);
	
	AVScreate_output_port ("Output",  "ucd"); 
	
	param = AVSadd_parameter("Node Data", "string", "Node Data", "Node Data",
				 NULL);
	AVSconnect_widget (param, "text");
	
	AVSadd_parameter ("node data", "choice", "<data 1>", choices, ".");
	
	AVSset_init_proc (ucd_extract_init);
	
	AVSset_destroy_proc (ucd_extract_finis);
	
	AVSset_compute_proc (ucd_extract);
	
	AVSset_module_flags (COOPERATIVE);
}

AVSinit_modules()
{
	int ucd_extract_desc();
	
	AVSmodule_from_desc(ucd_extract_desc);
}
