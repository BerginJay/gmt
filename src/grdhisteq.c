/*--------------------------------------------------------------------
 *	$Id$
 *
 *	Copyright (c) 1991-2012 by P. Wessel, W. H. F. Smith, R. Scharroo, and J. Luis
 *	See LICENSE.TXT file for copying and redistribution conditions.
 *
 *	This program is free software; you can redistribute it and/or modify
 *	it under the terms of the GNU Lesser General Public License as published by
 *	the Free Software Foundation; version 3 or any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU Lesser General Public License for more details.
 *
 *	Contact info: gmt.soest.hawaii.edu
 *--------------------------------------------------------------------*/
/*
 * Brief synopsis: read a grid file and find the values which divide its range
 * into n_cell number of quantiles.
 *
 * Author:	W.H.F. Smith
 * Date: 	31 May 1990
 * Version:	5 API
 */
 
#define THIS_MODULE k_mod_grdhisteq /* I am grdhisteq */

#include "gmt.h"

struct GRDHISTEQ_CTRL {
	struct In {
		bool active;
		char *file;
	} In;
	struct C {	/* -C<n_cells>*/
		bool active;
		unsigned int value;
	} C;
	struct D {	/* -D[<file>] */
		bool active;
		char *file;
	} D;
	struct G {	/* -G<file> */
		bool active;
		char *file;
	} G;
	struct N {	/* -N[<norm>] */
		bool active;
		double norm;
	} N;
	struct Q {	/* -Q */
		bool active;
	} Q;
};

struct	INDEXED_DATA {
	float x;
	int i;
};

struct	CELL {
	float low;
	float high;
};

void *New_grdhisteq_Ctrl (struct GMT_CTRL *GMT) {	/* Allocate and initialize a new control structure */
	struct GRDHISTEQ_CTRL *C = NULL;
	
	C = GMT_memory (GMT, NULL, 1, struct GRDHISTEQ_CTRL);
	
	/* Initialize values whose defaults are not 0/false/NULL */
	return (C);
}

void Free_grdhisteq_Ctrl (struct GMT_CTRL *GMT, struct GRDHISTEQ_CTRL *C) {	/* Deallocate control structure */
	if (!C) return;
	if (C->In.file) free (C->In.file);	
	if (C->D.file) free (C->D.file);	
	if (C->G.file) free (C->G.file);	
	GMT_free (GMT, C);	
}

int GMT_grdhisteq_usage (struct GMTAPI_CTRL *C, int level)
{
	struct GMT_CTRL *GMT = C->GMT;

	gmt_module_show_name_and_purpose (THIS_MODULE);
	GMT_message (GMT, "usage: grdhisteq <ingrid> [-G<outgrid>] [-C<n_cells>] [-D[<table>]] [-N[<norm>]] [-Q]\n\t[%s] [%s]\n", GMT_Rgeo_OPT, GMT_V_OPT);
	
	if (level == GMTAPI_SYNOPSIS) return (EXIT_FAILURE);
	
	GMT_message (GMT, "\t<ingrid> is name of input grid file.\n");
	GMT_message (GMT, "\n\tOPTIONS:\n");
	GMT_message (GMT, "\t-C Set how many cells (divisions) of data range to make.\n");
	GMT_message (GMT, "\t-D Dump level information to <table> or stdout if not given.\n");
	GMT_message (GMT, "\t-G Create an equalized output grid file called <outgrid>.\n");
	GMT_message (GMT, "\t-N Use with -G to make an output grid file with standard normal scores.\n");
	GMT_message (GMT, "\t   Append <norm> to normalize the scores to <-1,+1>.\n");
	GMT_message (GMT, "\t-Q Use quadratic intensity scaling [Default is linear].\n");
	GMT_explain_options (GMT, "RV.");
	
	return (EXIT_FAILURE);
}

int GMT_grdhisteq_parse (struct GMTAPI_CTRL *C, struct GRDHISTEQ_CTRL *Ctrl, struct GMT_OPTION *options)
{
	/* This parses the options provided to grdhisteq and sets parameters in Ctrl.
	 * Note Ctrl has already been initialized and non-zero default values set.
	 * Any GMT common options will override values set previously by other commands.
	 * It also replaces any file names specified as input or output with the data ID
	 * returned when registering these sources/destinations with the API.
	 */

	unsigned int n_errors = 0, n_files = 0;
	int sval;
	struct GMT_OPTION *opt = NULL;
	struct GMT_CTRL *GMT = C->GMT;

	for (opt = options; opt; opt = opt->next) {	/* Process all the options given */

		switch (opt->option) {
			case '<':	/* Input file (only one is accepted) */
				Ctrl->In.active = true;
				if (n_files++ == 0) Ctrl->In.file = strdup (opt->arg);
				break;

			/* Processes program-specific parameters */

			case 'C':	/* Get # of cells */
				Ctrl->C.active = true;
				sval = atoi (opt->arg);
				n_errors += GMT_check_condition (GMT, sval <= 0, "Syntax error -C option: n_cells must be positive\n");
				Ctrl->C.value = sval;
				break;
			case 'D':	/* Dump info to file or stdout */
				Ctrl->D.active = true;
				if (opt->arg[0]) Ctrl->D.file = strdup (opt->arg);
				break;
			case 'G':	/* Output file for equalized grid */
				Ctrl->G.active = true;
				Ctrl->G.file = strdup (opt->arg);
				break;
			case 'N':	/* Get normalized scores */
				Ctrl->N.active = true;
				Ctrl->N.norm = atof (opt->arg);
				break;
			case 'Q':	/* Use quadratic scaling */
				Ctrl->Q.active = true;
				break;

			default:	/* Report bad options */
				n_errors += GMT_default_error (GMT, opt->option);
				break;
		}
	}

	n_errors += GMT_check_condition (GMT, n_files > 1, "Syntax error: Must specify a single input grid file\n");
	n_errors += GMT_check_condition (GMT, !Ctrl->In.file, "Syntax error: Must specify input grid file\n");
	n_errors += GMT_check_condition (GMT, Ctrl->N.active && !Ctrl->G.file, "Syntax error -N option: Must also specify output grid file with -G\n");
	n_errors += GMT_check_condition (GMT, !strcmp (Ctrl->In.file, "="), "Syntax error: Piping of input grid file not supported!\n");

	return (n_errors ? GMT_PARSE_ERROR : GMT_OK);
}

float get_cell (float x, struct CELL *cell, unsigned int n_cells_m1, unsigned int last_cell)
{
	unsigned int low, high, i;

	low = 0;
	high = n_cells_m1;
	i = last_cell;

	do {
		if (cell[i].low <= x && cell[i].high >= x) {
			last_cell = i;
			return ((float)i);
		}
		else if (cell[low].low <= x && cell[low].high >= x) {
			return ((float)low);
		}
		else if (cell[high].low <= x && cell[high].high >= x) {
			return ((float)high);
		}
		else if (cell[i].low > x) {
			high = i;
			i = (low + high) / 2;
		}
		else if (cell[i].high < x) {
			low = i;
			i = (low + high) / 2;
		}
	} while (true);
	return (0.0);	/* Cannot get here - just used to quiet compiler */
}

int do_hist_equalization (struct GMT_CTRL *GMT, struct GMT_GRID *Grid, char *outfile, unsigned int n_cells, bool quadratic, bool dump_intervals)
{	/* Do basic histogram equalization */
	uint64_t i, j, nxy;
	unsigned int last_cell, n_cells_m1 = 0, current_cell, pad[4];
	double delta_cell, target, out[3];
	struct CELL *cell = NULL;
	struct GMT_GRID *Orig = NULL;
	
	cell = GMT_memory (GMT, NULL, n_cells, struct CELL);

	/* Sort the data and find the division points */

	GMT_memcpy (pad, Grid->header->pad, 4, int);	/* Save the original pad */
	GMT_grd_pad_off (GMT, Grid);	/* Undo pad if one existed so we can sort the entire grid */
	if (outfile) Orig = GMT_duplicate_grid (GMT, Grid, true); /* Must keep original if readonly */
	GMT_sort_array (GMT, Grid->data, Grid->header->nm, GMTAPI_FLOAT);
	
	nxy = Grid->header->nm;
	while (nxy > 0 && GMT_is_fnan (Grid->data[nxy-1])) nxy--;	/* Only deal with real numbers */

	last_cell = n_cells / 2;
	n_cells_m1 = n_cells - 1;

	current_cell = i = 0;
	delta_cell = ((double)nxy) / ((double)n_cells);

	while (current_cell < n_cells) {

		if (current_cell == (n_cells - 1))
			j = nxy - 1;
		else if (quadratic) {	/* Use y = 2x - x**2 scaling  */
			target = (current_cell + 1.0) / n_cells;
			j = lrint (floor (nxy * (1.0 - sqrt (1.0 - target))));
		}
		else	/* Use simple linear scale  */
			j = lrint (floor ((current_cell + 1) * delta_cell)) - 1;

		cell[current_cell].low  = Grid->data[i];
		cell[current_cell].high = Grid->data[j];

		if (dump_intervals) {	/* Write records to file or stdout */
			out[GMT_X] = (double)Grid->data[i]; out[GMT_Y] = (double)Grid->data[j]; out[GMT_Z] = (double)current_cell;
			GMT_Put_Record (GMT->parent, GMT_WRITE_DOUBLE, out);
		}

		i = j;
		current_cell++;
	}
	if (dump_intervals && GMT_End_IO (GMT->parent, GMT_OUT, 0) != GMT_OK) {	/* Disables further data ioutput */
		return (GMT->parent->error);
	}
	
	if (outfile) {	/* Must re-read the grid and evaluate since it got sorted and trodden on... */
		for (i = 0; i < Grid->header->nm; i++) Grid->data[i] = (GMT_is_fnan (Orig->data[i])) ? GMT->session.f_NaN : get_cell (Orig->data[i], cell, n_cells_m1, last_cell);
		GMT_free_grid (GMT, &Orig, true);
	}

	GMT_grd_pad_on (GMT, Grid, pad);	/* Reinstate the original pad */
	GMT_free (GMT, cell);
	return (0);
}

int compare_indexed_floats (const void *point_1, const void *point_2)
{
	if (((struct INDEXED_DATA *)point_1)->x < ((struct INDEXED_DATA *)point_2)->x) return (-1);
	if (((struct INDEXED_DATA *)point_1)->x > ((struct INDEXED_DATA *)point_2)->x) return (+1);
	return (0);
}

int compare_indices (const void *point_1, const void *point_2)
{
	if (((struct INDEXED_DATA *)point_1)->i < ((struct INDEXED_DATA *)point_2)->i) return (-1);
	if (((struct INDEXED_DATA *)point_1)->i > ((struct INDEXED_DATA *)point_2)->i) return (1);
	return (0);
}

int do_gaussian_scores (struct GMT_CTRL *GMT, struct GMT_GRID *Grid, double norm)
{	/* Make an output grid file with standard normal scores */
	unsigned int row, col;
	uint64_t i = 0, j = 0, ij, nxy;
	double dnxy;
	struct INDEXED_DATA *indexed_data = NULL;

	indexed_data = GMT_memory (GMT, NULL, Grid->header->nm, struct INDEXED_DATA);

	nxy = Grid->header->nm;
	GMT_grd_loop (GMT, Grid, row, col, ij) {
		if (GMT_is_fnan (Grid->data[ij])) {	/* Put NaNs in the back */
			nxy--;
			indexed_data[nxy].i = ij;
			indexed_data[nxy].x = Grid->data[ij];
		}
		else {
			indexed_data[j].i = ij;
			indexed_data[j].x = Grid->data[ij];
			j++;
		}
	}

	/* Sort on data value  */

	qsort (indexed_data, nxy, sizeof (struct INDEXED_DATA), compare_indexed_floats);

	dnxy = 1.0 / (nxy + 1);

	if (norm != 0.0) norm /= fabs (GMT_zcrit (GMT, (double)dnxy));	/* Normalize by abs(max score) */

	for (i = 0; i < nxy; i++) {
		indexed_data[i].x = (float)GMT_zcrit (GMT, (i + 1.0) * dnxy);
		if (norm != 0.0) indexed_data[i].x *= (float)norm;
	}

	/* Sort on data index  */

	qsort (indexed_data, Grid->header->nm, sizeof (struct INDEXED_DATA), compare_indices);

	i = 0;
	GMT_grd_loop (GMT, Grid, row, col, ij) Grid->data[ij] = indexed_data[i++].x;	/* Load up the grid */

	GMT_free (GMT, indexed_data);
	return (0);
}

#define bailout(code) {GMT_Free_Options (mode); return (code);}
#define Return(code) {Free_grdhisteq_Ctrl (GMT, Ctrl); GMT_end_module (GMT, GMT_cpy); bailout (code);}

int GMT_grdhisteq (struct GMTAPI_CTRL *API, int mode, void *args)
{
	bool error = false;

	double wesn[4];
	
	struct GMT_GRID *Grid = NULL, *Out = NULL;
	struct GRDHISTEQ_CTRL *Ctrl = NULL;
	struct GMT_CTRL *GMT = NULL, *GMT_cpy = NULL;
	struct GMT_OPTION *options = NULL;

	/*----------------------- Standard module initialization and parsing ----------------------*/

	if (API == NULL) return (GMT_Report_Error (API, GMT_NOT_A_SESSION));
	options = GMT_Prep_Options (API, mode, args);	if (API->error) return (API->error);	/* Set or get option list */

	if (!options || options->option == GMTAPI_OPT_USAGE) bailout (GMT_grdhisteq_usage (API, GMTAPI_USAGE));	/* Return the usage message */
	if (options->option == GMTAPI_OPT_SYNOPSIS) bailout (GMT_grdhisteq_usage (API, GMTAPI_SYNOPSIS));	/* Return the synopsis */

	/* Parse the command-line arguments */

	GMT = GMT_begin_gmt_module (API, THIS_MODULE, &GMT_cpy); /* Save current state */
	if (GMT_Parse_Common (API, "-VR", "", options)) Return (API->error);
	Ctrl = New_grdhisteq_Ctrl (GMT);	/* Allocate and initialize a new control structure */
	if ((error = GMT_grdhisteq_parse (API, Ctrl, options))) Return (error);

	/*---------------------------- This is the grdhisteq main code ----------------------------*/

	GMT_memcpy (wesn, GMT->common.R.wesn, 4, double);	/* Current -R setting, if any */
	if ((Grid = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_HEADER, NULL, Ctrl->In.file, NULL)) == NULL) {
		Return (API->error);
	}
	if (GMT_is_subset (GMT, Grid->header, wesn)) GMT_err_fail (GMT, GMT_adjust_loose_wesn (GMT, wesn, Grid->header), "");	/* Subset requested; make sure wesn matches header spacing */
	if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_DATA, wesn, Ctrl->In.file, Grid) == NULL) {	/* Get subset */
		Return (API->error);
	}
	(void)GMT_set_outgrid (GMT, Grid, &Out);	/* true if input is a read-only array */
	GMT_grd_init (GMT, Out->header, options, true);

	if (Ctrl->N.active)
		error = do_gaussian_scores (GMT, Out, Ctrl->N.norm);
	else {
		if (Ctrl->D.active) {	/* Initialize file/stdout for table output */
			int out_ID;
			if (Ctrl->D.file && (out_ID = GMT_Register_IO (API, GMT_IS_DATASET, GMT_IS_FILE, GMT_IS_POINT, GMT_OUT, NULL, Ctrl->D.file)) == GMTAPI_NOTSET) {
				Return (EXIT_FAILURE);
			}
			if ((error = GMT_set_cols (GMT, GMT_OUT, 3)) != GMT_OK) {
				Return (error);
			}
			if (GMT_Init_IO (API, GMT_IS_DATASET, GMT_IS_POINT, GMT_OUT, GMT_REG_DEFAULT, 0, options) != GMT_OK) {	/* Registers default output destination, unless already set */
				Return (API->error);
			}
			if (GMT_Begin_IO (API, GMT_IS_DATASET, GMT_OUT) != GMT_OK) {	/* Enables data output and sets access mode */
				Return (API->error);
			}
		}
		if ((error = do_hist_equalization (GMT, Out, Ctrl->G.file, Ctrl->C.value, Ctrl->Q.active, Ctrl->D.active))) Return (EXIT_FAILURE);	/* Read error */
		/* do_hist_equalization will also call GMT_End_IO if Ctrl->D.active was true */
	}
	if (Ctrl->G.active && GMT_Write_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_GRID_ALL, NULL, Ctrl->G.file, Out) != GMT_OK) {
		Return (API->error);
	}

	Return (EXIT_SUCCESS);
}
