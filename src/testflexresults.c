/* Fake GMT test module to scan the large rms table for Emperor flexure and
 * output the desired results. */

#include "gmt_dev.h"

#define MAX_SOLUTIONS	18063045

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
	unsigned int k, ns, last, next, n;
	unsigned int rhoi, rhor, x, rhoi_fixed, rhor_fixed, xpos_fixed;
	double Te, Te_fixed, mean, std;
	char line[BUFSIZ] = {""};
	FILE *fp = NULL, *fpo = NULL;
	struct GMTAPI_CTRL *API = NULL;

	for (k = 1; k < argc; k++) {	/* March across options */
		if (!strcmp (argv[k], "--rhoi")) {	/* e.g., --rhoi 2510 */
			k++;	rhoi_limit = 1;
			rhoi_fixed = atoi (argv[k]);
		}
		else if (!strcmp (argv[k], "--rhor")) {	/* e.g., --rhor 2830 */
			k++;	rhor_limit = 1;
			rhor_fixed = atoi (argv[k]);
		}
		else if (!strcmp (argv[k], "--Te")) {	/* e.g., --Te 16 */
			k++;	Te_limit = 1;
			Te_fixed = atof (argv[k]);
		}
		else if (!strcmp (argv[k], "--xpos")) {	/* e.g., --xpos 559 */
			k++;	xpos_limit = 1;
			xpos_fixed = atoi (argv[k]);
		}
		else if (!strcmp (argv[k], "--dump")) {	/* e.g., --dump subset.txt */
			k++;	fpo = fopen (argv[k], "w");
		}
	}

	if ((API = GMT_Create_Session ("FLEX", 0, GMT_SESSION_NORMAL, NULL)) == NULL)
		exit (EXIT_FAILURE);

	//if ((fp = fopen ("dump_all.txt", "r")) == NULL)
	if ((fp = fopen ("some_dump.txt", "r")) == NULL)
		exit (EXIT_FAILURE);

	fgets (line, BUFSIZ, fp);	/* Skip header */
	k = 0;
	while (fgets (line, BUFSIZ, fp)) {
		sscanf (line, "%d %lg %d %d %lg %lg %d\n", &x, &Te, &rhoi, &rhor, &mean, &std, &n);
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
	last = 0;
	while (last != ns) {
		next = last;
		/* Find first x that differ */
		while (next < ns && solution[next].x == solution[last].x) next++;
		/* Sort the section from last to next-1 on std */
		qsort (&solution[last], next-last, sizeof (struct SOLUTION), sort_on_std);
		k = last;	/* First entry in sorted list is the best */
		printf ("%d\t%4.1f\t%d\t%d\t%lg\t%lg\t%d\n", solution[k].x, solution[k].Te, solution[k].rhoi, solution[k].rhor,
			solution[k].mean, solution[k].std, solution[k].n);
		last = next;	/* Go to next set of xpos values */
	}

	if (GMT_Destroy_Session (API))
		exit (EXIT_FAILURE);
}
