/* Fake GMT test module to scan the large rms table for Emperor flexure and
 * output the desired results.  We read one file: SFILE below, and does a few things
 *
 * 1. We exclude solutions that do not match any given fixed paramters.
 *	    Use --xpos x to only keep solutions for position x
 *	    Use --rhoi x to only keep solutions for infill density x
 *	    Use --rhor x to only keep solutions for root density x
 *	    Use --Te x to only keep solutions for elastic thickness x
 *   Use --dump file to write the resulting subset to file
 *
 * 2. Given the subset we find the solution with lowers std per x
 *    This is written to stdout
 *
 * 3. Given the best Te(x), also write out a third file which only uses
 *    these Te values but writes out all the other solution for that x with
 *	  the same Te (e.g., other infill densities and rhor that gave the same)
 *    Use --best file to get this output.
 */

#include "gmt_dev.h"

#define MAX_SOLUTIONS	18063045
#define SFILE			"full_dump.txt"
#define N_X 1001		/* From x = 0 to x = 2000 in steps of 2 */
#define INC 2			/* Grid spacing in km */

struct SOLUTION {
	unsigned int x, rhoi, rhor, n;
	double Te, mean, std;
} solution[MAX_SOLUTIONS];

static int sort_on_x (const void *point_1v, const void *point_2v) {
	double x1, x2;
	const struct SOLUTION *point_1 = point_1v, *point_2 = point_2v;
	x1 = point_1->x;
	x2 = point_2->x;
	if (x1 < x2) return (-1);
	if (x1 > x2) return (+1);
	return (0);
}

static int sort_on_std (const void *point_1v, const void *point_2v) {
	double x1, x2;
	const struct SOLUTION *point_1 = point_1v, *point_2 = point_2v;
	x1 = point_1->std;
	x2 = point_2->std;
	if (x1 < x2) return (-1);
	if (x1 > x2) return (+1);
	return (0);
}

int main (int argc, char *argv[]) {
	unsigned int rhoi_limit = 0, rhor_limit = 0, Te_limit = 0, xpos_limit = 0;
	unsigned int ix, k, ns, last, next, n, n_total = 0;
	unsigned int rhoi, rhor, x, rhoi_fixed, rhor_fixed, xpos_fixed;
	double Te, Te_fixed, mean, std, best_Te[N_X];
	char line[BUFSIZ] = {""};
	FILE *fp = NULL, *fpo = NULL, *fpb = NULL;
	struct GMTAPI_CTRL *API = NULL;

	for (k = 1; k < argc; k++) {	/* March across options */
		if (!strcmp (argv[k], "--rhoi")) {	/* e.g., --rhoi 2510 */
			k++;	rhoi_limit = 1;
			rhoi_fixed = atoi (argv[k]);
			fprintf (stderr, "Excluding solutions with rhoi != %d\n", rhoi_fixed);
		}
		else if (!strcmp (argv[k], "--rhor")) {	/* e.g., --rhor 2830 */
			k++;	rhor_limit = 1;
			rhor_fixed = atoi (argv[k]);
			fprintf (stderr, "Excluding solutions with rhor != %d\n", rhor_fixed);
		}
		else if (!strcmp (argv[k], "--Te")) {	/* e.g., --Te 16 */
			k++;	Te_limit = 1;
			Te_fixed = atof (argv[k]);
			fprintf (stderr, "Excluding solutions with Te != %04.1lf\n", Te_fixed);
		}
		else if (!strcmp (argv[k], "--xpos")) {	/* e.g., --xpos 559 */
			k++;	xpos_limit = 1;
			xpos_fixed = atoi (argv[k]);
			fprintf (stderr, "Excluding solutions with xpos != %d\n", xpos_fixed);
		}
		else if (!strcmp (argv[k], "--dump")) {	/* e.g., --dump subset.txt */
			k++;	fpo = fopen (argv[k], "w");
			fprintf (stderr, "Write selected subset to file %s\n", argv[k]);
		}
		else if (!strcmp (argv[k], "--best")) {	/* e.g., --best answer.txt */
			k++;	fpb = fopen (argv[k], "w");
			fprintf (stderr, "Write best results to file %s\n", argv[k]);
		}
	}

	if ((API = GMT_Create_Session ("FLEX", 0, GMT_SESSION_NORMAL, NULL)) == NULL)
		exit (EXIT_FAILURE);

	if ((fp = fopen (SFILE, "r")) == NULL) {
			fprintf (stderr, "Unable to open file %s\n", SFILE);
		exit (EXIT_FAILURE);
	}

	fgets (line, BUFSIZ, fp);	/* Skip header */
	k = 0;
	while (fgets (line, BUFSIZ, fp)) {
		sscanf (line, "%d %lg %d %d %lg %lg %d\n", &x, &Te, &rhoi, &rhor, &mean, &std, &n);
		n_total++;	/* Solutions found */
		if (rhoi_limit && rhoi != rhoi_fixed) continue;			/* Not selected rhoi */
		if (rhor_limit && rhor != rhor_fixed) continue;			/* Not selected rhor */
		if (xpos_limit && x != xpos_fixed) continue;			/* Not selected x position */
		if (Te_limit && fabs (Te-Te_fixed) > 0.01) continue;	/* Not selected Te position */
		solution[k].x = x;
		solution[k].Te = Te;
		solution[k].rhoi = rhoi;
		solution[k].rhor = rhor;
		solution[k].mean = mean;
		solution[k].std = std;
		solution[k].n = n;
		k++;
	}
	ns = k;	/* We now have ns solutions in memory */
	fprintf (stderr, "Selected %d from a total of %d solutions\n", ns, n_total);

	/* Now sort these on x */
	qsort (solution, ns, sizeof (struct SOLUTION), sort_on_x);

	if (fpo) {	/* Want to dump those now */
		fprintf (fpo, "# xpos\tTe\trhoi\trhor\tmean\tstd\tnpoints\n");
		for (k = 0; k < ns; k++)
			fprintf (fpo, "%d\t%4.1f\t%d\t%d\t%lg\t%lg\t%d\n", solution[k].x, solution[k].Te, solution[k].rhoi, solution[k].rhor,
				solution[k].mean, solution[k].std, solution[k].n);
		fclose (fpo);
	}

	/* Now report best solution per xpos */
	last = k = 0;
	printf ("# xpos\tTe\trhoi\trhor\tmean\tstd\tnpoints\n");
	while (last != ns) {
		next = last;
		/* Find first x that differ */
		while (next < ns && solution[next].x == solution[last].x) next++;
		/* Sort the section from last to next-1 on std */
		qsort (&solution[last], next-last, sizeof (struct SOLUTION), sort_on_std);
		/* First entry in sorted list (Position last) is the best */
		printf ("%d\t%4.1f\t%d\t%d\t%lg\t%lg\t%d\n", solution[last].x, solution[last].Te, solution[last].rhoi, solution[last].rhor,
			solution[last].mean, solution[last].std, solution[last].n);
		best_Te[k++] = solution[last].Te;	/* Keep the best Te for each position */
		last = next;	/* Go to next set of xpos values */
	}

	if (fpb) {	/* Now report all density variations for the best Te per xpos */
		/* Re-sort these on x again */
		qsort (solution, ns, sizeof (struct SOLUTION), sort_on_x);
		fprintf (fpb, "# xpos\tTe\trhoi\trhor\tmean\tstd\tnpoints\n");
		for (k = 0; k < ns; k++) {
			ix = solution[k].x / INC;	/* x-index */
			if (fabs (solution[k].Te - best_Te[ix]) > 0.01) continue;
			fprintf (fpb, "%d\t%4.1f\t%d\t%d\t%lg\t%lg\t%d\n", solution[k].x, solution[k].Te, solution[k].rhoi, solution[k].rhor,
				solution[k].mean, solution[k].std, solution[k].n);
		}
		fclose (fpb);
	}

	if (GMT_Destroy_Session (API))
		exit (EXIT_FAILURE);
}
