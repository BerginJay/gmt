#include "gmt_dev.h"

#define N_RHOI 9		/* Number of fixed infill densities */
#define N_RHOR 5		/* Number of fixed root densities */
#define N_TE 391		/* Number of elastic thicknesses */
#define N_X 1001		/* From x = 0 to x = 2000 in steps of 2 */

#define H 500			/* Half-height for y-box in search in km */
#define W 250			/* Half-width for y-box in search in km */
#define dTe 0.1			/* Te spacing in km */
#define TESTART	1.0		/* Start Te calculations */
#define INC 2			/* Grid spacing in km */
#define XSTART	0		/* Start point along x-axis for calculations */
#define XSTOP	2000	/* Start point along x-axis for calculations */

#define get_node(h,row,col) ((uint64_t)(((int64_t)(row))*((int64_t)h->n_columns)+(int64_t)(col)))

//#define CHECK

int main () {
	unsigned int ki, kr, kt, kx, col_start, col_stop, col, row, col_min, col_max, row_min, row_max, node, n, half_width, half_height;

	char gmodel[250] = {""};	/* Hold name of current FAA model */

	int rhoi[N_RHOI] = {2310, 2360, 2410, 2460, 2510, 2560, 2610, 2660, 2710};	/* List of infill densities */
	int rhor[N_RHOR] = {2730, 2780, 2830, 2880, 2930};	/* List of root densities */
	int x[N_X];	/* List of x-positions */
	double Te[N_TE];	/* List of Te values */
	double sum, sum2, d_faa, mean, std;

	struct GMT_GRID *G_obs = NULL, *G_model = NULL, *Mask = NULL;
	struct GMT_GRID_HEADER *h = NULL;

	struct GMTAPI_CTRL *API = NULL;

	for (kt = 0; kt < N_TE; kt++) Te[kt] = TESTART + kt * dTe;		/* Te array in km every 0.1 km */
	for (kx = 0; kx < N_X;  kx++) x[kx]  = XSTART  + kx * INC;		/* x array in km every 2 km */

	if ((API = GMT_Create_Session ("FLEX", 0, GMT_SESSION_NORMAL, NULL)) == NULL)
		exit (EXIT_FAILURE);

	/* Read FAA obs grid */
	if ((G_obs = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_CONTAINER_AND_DATA, NULL, "E_obl_full_area_faa_masked_ors_removed.grd", NULL)) == NULL)
		exit (EXIT_FAILURE);

	/* Read SID mask grid */
	if ((Mask = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_CONTAINER_AND_DATA, NULL, "E_obl_full_area_sid_mask.grd", NULL)) == NULL)
		exit (EXIT_FAILURE);

	h = G_obs->header;	/* Pointer to one of the grid headers */
	col_start = (unsigned int)rint ((XSTART - h->wesn[0]) / INC);	/* column in grid corresponding to x = XSTART */
	col_stop  = (unsigned int)rint ((XSTOP  - h->wesn[0]) / INC);	/* column in grid corresponding to x = XSTOP */
	row_min = (unsigned int)rint ((h->wesn[3] - (+H)) / INC);	/* row in grid corresponding to y = +H */
	row_max = (unsigned int)rint ((h->wesn[3] - (-H)) / INC);	/* row in grid corresponding to y = -H */
	half_width = W / INC;	half_height = H / INC;

	printf ("# xpos\tTe\trhoi\trhor\tmean\tstd\tnpoints\n");
#ifdef CHECK
	ki = 4;	/* Only do 1510 */
#else
	for (ki = 0; ki < N_RHOI; ki++)
#endif
	{
#ifdef CHECK
		kr = 2;	/* Only do 2830 */
#else
		for (kr = 0; kr < N_RHOR; kr++)
#endif
		{
			for (kt = 0; kt < N_TE; kt++) {
				sprintf (gmodel, "Egravdir/E_obl_FAA_R%d_I%d_%04.1fk_Total.grd", rhor[kr], rhoi[ki], Te[kt]);
				fprintf (stderr, "Processing %s\n", gmodel);
				/* Read model FAA (ki, kr, kt) grid */
				if ((G_model = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_CONTAINER_AND_DATA, NULL, gmodel, NULL)) == NULL)
					exit (-1);
				for (kx = 0; kx < N_X; kx++) {
					/* Determine col_min, col_max, row_min, row_max for grid subset */
					sum = sum2 = 0.0;	n = 0;	/* Reset counters */
					col_min = MAX (0, col_start + kx - half_width);	/* Start to the left but not off the grid */
					col_max = MIN (h->n_columns, col_start + kx + half_width);	/* Stop to the right but no off the grid */
					for (row = row_min; row < row_max; row++) {	/* Rows in subset */
						for (col = col_min; col < col_max; col++) {	/* Columns in subset */
							node = get_node (h, row, col);	/* Node in the grids */
							if (gmt_M_is_fnan (Mask->data[node])) continue;		/* No observation */
							if (gmt_M_is_fnan (G_obs->data[node])) continue;	/* No observation */
							d_faa = G_obs->data[node] - G_model->data[node];
							sum  += d_faa;
							sum2 += d_faa * d_faa;
							n++;
						}
					}
					if (n < 2) continue;	/* Not enough points */
					mean = sum / n;
					std = sqrt ((sum2 - mean * sum) / (n - 1.0));
					/* Report results */
					printf ("%d\t%4.1f\t%d\t%d\t%lg\t%lg\t%d\n", x[kx], Te[kt], rhoi[ki], rhor[kr], mean, std, n);
				}
				/* Free model gravity */
				if (GMT_Destroy_Data (API, &G_model) != GMT_NOERROR)
					exit (-1);
			}
		}
	}
	/* Free observed gratuity and mask */
	if (GMT_Destroy_Data (API, &G_obs) != GMT_NOERROR)
		exit (EXIT_FAILURE);
	if (GMT_Destroy_Data (API, &Mask) != GMT_NOERROR)
		exit (EXIT_FAILURE);
	if (GMT_Destroy_Session (API))
		exit (EXIT_FAILURE);
}
