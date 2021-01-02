/*--------------------------------------------------------------------
 *	Copyright (c) 1991-2020 by the GMT Team (https://www.generic-mapping-tools.org/team.html)
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
 *	Contact info: www.generic-mapping-tools.org
 *--------------------------------------------------------------------*/
/* grdslice reads a 2-D grid and detects isolated peaks by using
 * slice contouring.  We only examine interior, closed contours.
 * We contour from the maximum value and go downwards.  As the contour value
 * drops we will chop off the tops of smaller and smaller peaks.  The
 * first time we do this we initialize a new peak object and then the
 * contour rings further down are appended via pointers to the initial slice.
 * As the contours widen it is possible that more than one peak share
 * the same, broader base slices.  We may then traverse the list of peak
 * and their slices to do analysis.  Each slice stores its area and mean
 * location and fitted ellipsoidal parameters.
 * The usual input grid will be geographic lon/lat. However, we also accept
 * Mercator grids or Cartesian.  The analysis will be done using projected
 * Mercator (or original Cartesian) units but the results will be reported
 * in geographic (if input is geographic or Mercator).
 * The original application of this tool was to detect seamounts from grids
 * of gravity or vertical gradient anomalies [see ref list in docs].
 *
 * Author:	Paul Wessel and Seung-Sep Kim
 * Date:	1-DEC-2020 (original date 5-DEC-2006 and partly based on grdcontour)
 * Version:	2.0	Revised for GMT 6.2
 */

#include "gmt_dev.h"

#define THIS_MODULE_CLASSIC_NAME	"grdslice"
#define THIS_MODULE_MODERN_NAME	"grdslice"
#define THIS_MODULE_LIB		"potential"
#define THIS_MODULE_PURPOSE	"Detect isolated peaks by contour-slicing a grid"
#define THIS_MODULE_KEYS	"<G{,>D},DD),ED),FD)"
#define THIS_MODULE_NEEDS	"g"
#define THIS_MODULE_OPTIONS "-:RVfn"

struct GRDSLICE_CTRL {
	struct GMT_CONTOUR contour;
	struct GRDSLICE_Out {	/* -> */
		bool active;
		char *file;
	} Out;
	struct GRDSLICE_In {
		bool active;
		char *file;
	} In;
	struct GRDSLICE_A {	/* -A<area> */
		bool active;
		double cutoff;
	} A;
	struct GRDSLICE_C {	/* -C<cont_int> */
		bool active;
		double interval;
	} C;
	struct GRDSLICE_D {	/* -D<dist> */
		bool active;
		double cutoff;
	} D;
	struct GRDSLICE_E {	/* -D<slicefile> */
		bool active;
		char *file;
	} E;
	struct GRDSLICE_F {	/* -F<foundation> */
		bool active;
		char *file;
	} F;
	struct GRDSLICE_I {	/* -I<indexfile> */
		bool active;
		char *file;
	} I;
	struct GRDSLICE_L {	/* -L<Low/high> */
		bool active;
		double low, high;
	} L;
	struct GRDSLICE_Q {	/* -Q<factor> */
		bool active;
		double factor;
	} Q;
	struct GRDSLICE_S {	/* -S<smooth> */
		bool active;
		unsigned int value;
	} S;
	struct GRDSLICE_T {  /* -T<foundation>] */
		bool active;
		double cutoff;
	} T;
	struct GRDSLICE_Z {	/* -Z[+s<fact>][+o<shift>] */
		bool active;
		double scale, offset;
	} Z;
};

struct GRDSLICE_SLICE {	/* Hold each contour slice information */
	int n;				/* Number of points in this slice polygon */
	int id;				/* Unique ID number */
	int shared;			/* Number of peaks that share this slice */
	int F_id;			/* Id of foundation this slice belongs to */
	double cval;			/* Contour of this slice */
	double z;			/* Z-value of this slice, initially at peak (so may be > cval) */
	double *x, *y;			/* The array of Mercator or Cartesian (x,y) coordinates */
	double xmin, xmax, ymin, ymax;	/* Extreme Mercator or Cartesian coordinates of polygon */
	double x_mean, y_mean;	/* Mean geographic or Cartesian coordinate as approximation of center point */
	double area;			/* Area of polygon in km^2 (or whatever Cartesian unit^2) */
	double azimuth;			/* Azimuth of major axis of approximate ellipse with same area */
	double major, minor;	/* Length of axes (in km or Cartesian unit) of approximate ellipse */
	double fit;				/* 0-100% of how well an ellipse explains the shape of contour (100% is perfect) */
	struct GRDSLICE_SLICE *next;	/* Pointer to next slice in same contour level */
	struct GRDSLICE_SLICE *down;	/* Pointer to the next slice down in this stack of slices */
};

struct GRDSLICE_PEAK {	/* Hold start of peak and linked list of slices */
	int id;				/* Unique ID number */
	double x, y;		/* Mean Mercator or Cartesian coordinate of peak */
	double cval;		/* Z-value (contour) of this peak */
	double z;			/* Z-value (contour) of this peak */
	struct GRDSLICE_SLICE *start;		/* Pointer to top slice in stack */
};

static void *New_Ctrl (struct GMT_CTRL *GMT) {	/* Allocate and initialize a new control structure */
	struct GRDSLICE_CTRL *C;
	
	C = gmt_M_memory (GMT, NULL, 1, struct GRDSLICE_CTRL);
	
	/* Initialize values whose defaults are not 0/false/NULL */
	
	C->L.low = 0;
	C->L.high = DBL_MAX;
	C->Q.factor = 8.0;
	C->Z.scale = 1.0;
	
	return (C);
}

void Free_Ctrl (struct GMT_CTRL *GMT, struct GRDSLICE_CTRL *C) {	/* Deallocate control structure */
	gmt_M_str_free (C->E.file);
	gmt_M_str_free (C->F.file);
	gmt_M_str_free (C->I.file);
	gmt_M_free (GMT, C);
}

static void grdslice_reverse_polygon (double x[], double y[], uint64_t n) {
	/* Reverses a polygon if needed so all are CCW */
	uint64_t i, j;
	
	for (i = 0, j = n-1; i < n/2; i++, j--) {
		gmt_M_double_swap (x[i], x[j]);
		gmt_M_double_swap (y[i], y[j]);
	}
}

static int usage (struct GMTAPI_CTRL *API, int level) {
	const char *name = gmt_show_name_and_purpose (API, THIS_MODULE_LIB, THIS_MODULE_CLASSIC_NAME, THIS_MODULE_PURPOSE);
	if (level == GMT_MODULE_PURPOSE) return (GMT_NOERROR);

	GMT_Message (API, GMT_TIME_NONE, "usage: %s <grid> -C<cont_int> [-A<area>] [-D<dist>] [-E<slicefile>] [-F<foundationfile>] [-I<indexfile>] [-L<low/high>]\n", name);
	GMT_Message (API, GMT_TIME_NONE, "\t[%s] [-S<smooth>] [-Q<divisor>] [-T<foundation>] [%s] [-Z[+s<scale>][+o<offset>]] [%s]\n\n", GMT_Rgeo_OPT, GMT_V_OPT, GMT_n_OPT);

	if (level == GMT_SYNOPSIS) return (GMT_MODULE_SYNOPSIS);

	GMT_Message (API, GMT_TIME_NONE, "\t<grid> is the grid file to be contour sliced.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-C Set the contour slice interval.\n");
	GMT_Message (API, GMT_TIME_NONE, "\n\tOPTIONS:\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-A Eliminate slices whose areas are less than <area> in km^2.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-D Skip peaks that are within <dist> km from a larger peak [no skipping].\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-E Set filename for optional slice information [no slices].\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-F Set filename for optional foundation information [no foundations].\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-I Set filename for optional index information [no indeces].\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-L Limit contours to this range [Default is -L0/inf].\n");
	GMT_Option (API, "R");
	GMT_Message (API, GMT_TIME_NONE, "\t-Q Sub-pixel division for grid-search revising peak location [8].\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-S Smooth contours by interpolation at approximately <gridsize>/<smooth> intervals.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t-T Set the foundation level for reporting peaks.\n");
	GMT_Message (API, GMT_TIME_NONE, "\t   Note: With -L, <foundation> must be equal or larger than the <low> value.\n");
	GMT_Option (API, "V");
	GMT_Message (API, GMT_TIME_NONE, "\t-Z Subtract <shift> (via +o<shift> [0]) and multiply data by <fact> (via +s<fact> [1]).\n");
	GMT_Option (API, "f,n,.");

	return (GMT_MODULE_USAGE);
}

static int parse (struct GMT_CTRL *GMT, struct GRDSLICE_CTRL *Ctrl, struct GMT_OPTION *options) {
	unsigned int n_errors = 0, n_files_in = 0, n_files_out = 0, pos = 0;
	int j;
	char p[GMT_LEN64] = {""};
	struct GMT_OPTION *opt = NULL;
	struct GMTAPI_CTRL *API = GMT->parent;

	for (opt = options; opt; opt = opt->next) {	/* Process all the options given */

		switch (opt->option) {

			case '<':	/* Input file (only one is accepted) */
				if (n_files_in++ > 0) { n_errors++; continue; }
				Ctrl->In.active = true;
				if (opt->arg[0]) Ctrl->In.file = strdup (opt->arg);
				if (GMT_Get_FilePath (GMT->parent, GMT_IS_GRID, GMT_IN, GMT_FILE_REMOTE, &(Ctrl->In.file))) n_errors++;
				break;
			case '>':	/* Got named output file */
				if (n_files_out++ > 0) { n_errors++; continue; }
				Ctrl->Out.active = true;
				if (opt->arg[0]) Ctrl->Out.file = strdup (opt->arg);
				if (GMT_Get_FilePath (GMT->parent, GMT_IS_DATASET, GMT_OUT, GMT_FILE_LOCAL, &(Ctrl->Out.file))) n_errors++;
				break;

			/* Processes program-specific parameters */

			case 'A':	/* Ignore contour slices with area less tha cutoff */
				Ctrl->A.active = true;
				Ctrl->A.cutoff = atof (opt->arg);
				break;
			case 'C':	/* Contour slice interval */
				Ctrl->C.active = true;
				Ctrl->C.interval = atof (opt->arg);
				break;
			case 'D':	/* Ignore small peaks too close to larger peaks */
				Ctrl->D.active = true;
				Ctrl->D.cutoff = atof (opt->arg);
				break;
			case 'E':	/* Slice file */
				Ctrl->E.active = true;
				gmt_M_str_free (Ctrl->E.file);
				Ctrl->E.file = strdup (opt->arg);
				break;
			case 'F':	/* Foundation file */
				Ctrl->F.active = true;
				gmt_M_str_free (Ctrl->F.file);
				Ctrl->F.file = strdup (opt->arg);
				break;
			case 'I':	/* Index file */
				Ctrl->I.active = true;
				gmt_M_str_free (Ctrl->I.file);
				Ctrl->I.file = strdup (opt->arg);
				break;
			case 'L':	/* Limit the contour range */
				Ctrl->L.active = true;
				sscanf (opt->arg, "%lf/%lf", &Ctrl->L.low, &Ctrl->L.high);
				break;
			case 'Q':	/* Sub-pixel factor */
				Ctrl->Q.active = true;
				Ctrl->Q.factor = atof (opt->arg);
				break;
			case 'S':	/* Smooth the contours */
				Ctrl->S.active = true;
				j = atoi (opt->arg);
				n_errors += gmt_M_check_condition (GMT, j < 0, "Option -S: Smooth_factor must be > 0\n");
				Ctrl->S.value = j;
				break;
			case 'T':	/* Set foundation contour and area limits */
				Ctrl->T.active = true;
				Ctrl->T.cutoff  = atof (opt->arg);
				break;
			case 'Z':	/* Sake/offset grid before analysis */
				Ctrl->Z.active = true;
				while (gmt_getmodopt (GMT, 'Z', opt->arg, "so", &pos, p, &n_errors) && n_errors == 0) {
					switch (p[0]) {
						case 's':	Ctrl->Z.scale  = atof (&p[1]);	break;
						case 'o':	Ctrl->Z.offset = atof (&p[1]);	break;
						default: 	/* These are caught in gmt_getmodopt so break is just for Coverity */
							break;
					}
				}
				break;
			default:	/* Report bad options */
				n_errors += gmt_default_error (GMT, opt->option);
				break;
		}
	}

	n_errors += gmt_M_check_condition (GMT, n_files_in != 1 || Ctrl->In.file == NULL, "Must specify a single grid file\n");
	n_errors += gmt_M_check_condition (GMT, Ctrl->C.interval <= 0.0, "Option -C: Must specify positive contour interval\n");
	n_errors += gmt_M_check_condition (GMT, Ctrl->L.low > Ctrl->L.high, "Option -L: lower limit > upper!\n");
	n_errors += gmt_M_check_condition (GMT, Ctrl->T.active && Ctrl->T.cutoff < Ctrl->L.low, "Option -T: Bottom level < lower limit set in -L!\n");
	n_errors += gmt_M_check_condition (GMT, Ctrl->Z.scale == 0.0, "Option -Z: factor must be nonzero\n");

	return (n_errors ? GMT_PARSE_ERROR : GMT_NOERROR);

}

GMT_LOCAL double grdslice_centroid_area (struct GMT_CTRL *GMT, double *x, double *y, unsigned int n, struct GRDSLICE_SLICE *p, struct GMT_GRID *G, unsigned int geo, bool refine, double factor, double pos[]) {
	unsigned int i, j, nx, ny;
	double area, xmin, xmax, ymin, ymax, zz, *xx = NULL, *yy = NULL;
	struct GMT_DATASEGMENT *P = NULL;
	struct GMT_GRID_HEADER_HIDDEN *HH = NULL;

	if (!refine) {	/* Just get area and centroid */
		area = gmt_centroid_area (GMT, x, y, n, geo, pos);	/* First get the area only */

		if (area > 0.0)	/* Enforce CCW polygons only */
			grdslice_reverse_polygon (x, y, n);
		else
			area = fabs (area);

		return area;
	}

	/* Get here if we need to refine the peak location. (x,y) are in same coordinate system as grid G */

	HH = gmt_get_H_hidden (G->header);
	P = GMT_Alloc_Segment (GMT->parent, GMT_NO_STRINGS, 0, 2, NULL, NULL);
	P->data[GMT_X] = x;
	P->data[GMT_Y] = y;
	pos[GMT_X] = p->x_mean;	/* Initial location is mean of contour perimeter */
	pos[GMT_Z] = p->y_mean;
	pos[GMT_Z] = p->cval;	/* Initial peak value is contour value */
	P->n_rows = n;
	gmt_set_seg_minmax (GMT, GMT_IS_POLY, 2, P);	/* Update min/max of x/y only */

	/* Round the min/max coordinates outwards to grid steps */
	xmin = floor (P->min[GMT_X] / G->header->inc[GMT_X]) * G->header->inc[GMT_X];
	xmax = ceil  (P->max[GMT_X] / G->header->inc[GMT_X]) * G->header->inc[GMT_X];
	ymin = floor (P->min[GMT_Y] / G->header->inc[GMT_Y]) * G->header->inc[GMT_Y];
	ymax = ceil  (P->max[GMT_Y] / G->header->inc[GMT_Y]) * G->header->inc[GMT_Y];
	nx = gmt_make_equidistant_array (GMT, xmin, xmax, G->header->inc[GMT_X] / factor, &xx);
	ny = gmt_make_equidistant_array (GMT, ymin, ymax, G->header->inc[GMT_Y] / factor, &yy);
	for (j = 0; j < ny; j++) {
		for (i = 0; i < nx; i++) {
			if (gmt_inonout (GMT, xx[i], yy[j], P) == GMT_OUTSIDE) continue;
			zz = (HH->has_NaNs == GMT_GRID_NO_NANS) ? gmt_bcr_get_z_fast (GMT, G, xx[i], yy[j]) : gmt_bcr_get_z (GMT, G, xx[i], yy[j]);
			if (zz > pos[GMT_Z]) {	/* Update location and value of largest value inside the contour */
				pos[GMT_X] = xx[i];
				pos[GMT_Y] = yy[j];
				pos[GMT_Z] = zz;
			}
		}
	}
	/* Blank the P arrays and wipe memory */
	P->data[GMT_X] = P->data[GMT_Y] = NULL;
	gmt_free_segment (GMT, &P);
	gmt_M_free (GMT, xx);
	gmt_M_free (GMT, yy);
	return (0.0);
}

GMT_LOCAL void grdslice_fit_ellipse (struct GMT_CTRL *GMT, double *dx, double *dy, unsigned int n, struct GRDSLICE_SLICE *p, double *pos, double area) {
	/* Find orientation of major/minor axes and aspect ratio from reduced, projected x,y coordinates */
	double A[4] = {0.0, 0.0, 0.0, 0.0}, EigenValue[2], EigenVector[4], angle, aspect, f, g, cos_g, sin_g, xr, yr;
	double minor, major, Sx = 0.0, Sy = 0.0, r_fit, r, dr, rms, work1[2], work2[2];
	double *b = gmt_M_memory (GMT, NULL, n, double);
	unsigned int M = 2, nrots, i;
	struct GMTAPI_CTRL *API = GMT->parent;

	/* Object is to find x, y, azimuth, major, minor axes based on counter in x, y.  This is 5 parameters in a nonlinear inversion.
	 * Here w cheat a bit and do this:
	 * 1) Use the centroid of x,y as the ellipse center - these are already computed in pos[].
	 * 2) Use eigenvalue analysis to determine the azimuth of the ellipse
	 * 3) Do a simple L2 fit to the contour points to determine the aspect ration major/minor
	 * 4) Do another simple L2 fit to determine the radial scaling that gives the final major/minor values in the x,y system
	 */
	for (i = 0; i < n; i++) {	/* Build dot-product 2x2 matrix using Mercator or Cartesian coordinates */
		dx[i] -= pos[GMT_X];	/* Compute deviations from the mean location */
		dy[i] -= pos[GMT_Y];
		A[0] += dx[i] * dx[i];	/* Add up the covariances */
		A[1] += dx[i] * dy[i];
		A[3] += dy[i] * dy[i];
		b[i] = atan2 (dy[i], dx[i]) * R2D;	/* Compute angle from origin to each point */
	} 
	A[2] = A[1];	/* Because of symmetry */
	if (gmt_jacobi (GMT, A, M, M, EigenValue, EigenVector, work1, work2, &nrots)) {	/* Solve eigen-system A = EigenVector * EigenValue * EigenVector^T */
		GMT_Report (API, GMT_MSG_WARNING, "Eigenvalue routine failed to converge in 50 sweeps.\n");
	}
	angle = atan2 (EigenVector[1], EigenVector[0]) * R2D;	/* Actual, not azimuth yet - just angle CCW from horizontal */
	p->azimuth = 90.0 - angle;		/* Now it is a proper azimuth */
	if (p->azimuth <= -180.0) p->azimuth += 360.0;
	if (p->azimuth > 180.0) p->azimuth -= 360.0;
	
	/* Fit model to dx, dy by minimizing dr^2 and solve for aspect ratio. We use the form
	 * x' = major * cos(g) and y' = minor * sin(g), form E = (dx - x')^2 + (dy - y')^2 and solve dE/dmajor = 0 etc*/

	for (i = 0; i < n; i++) {
		b[i] -= angle;	/* Now b has angles from major axis */
		sincosd (b[i], &sin_g, &cos_g);
		Sx += dx[i] * cos_g;
		Sy += dy[i] * sin_g;
	}
	aspect = Sx / Sy;

	/* Given final aspect ratio and actual area we can solve for axes in km via A = major * minor * pi */
	p->minor = sqrt (p->area / (aspect * M_PI));	/* Ellipse axes in km */
	p->major = p->minor * aspect;
	minor = sqrt (area / (aspect * M_PI));		/* Axes in map units, to be scaled below */
	major = minor * aspect;

	/* Find best L2 radial scale f for fitting dx, dy using E = (r_obs - f*r_fit)^2 and dE/df = 0 */
	for (i = 0, Sx = Sy = 0.0; i < n; i++) {
		sincosd (b[i], &sin_g, &cos_g);
		r = hypot (dx[i], dy[i]);	/* radius to observed point  */
		r_fit = major * minor / hypot (major * sin_g, minor * cos_g);	/* Predicted radius, but off by unknown factor f */
		Sx += r * r_fit;
		Sy += r_fit * r_fit;
	}
	f = Sx / Sy;	/* Scale that turns major and minor into the best-fit values for dx, dy */
	major *= f;	minor *= f;

	/* Finally compute a misfit value between observed contour and best-fit ellipse */
	for (i = 0; i < n; i++) {
		sincosd (b[i], &sin_g, &cos_g);
		r = hypot (dx[i], dy[i]);		/* radius to observed point  */
		r_fit = major * minor / hypot (major * sin_g, minor * cos_g);	/* Predicted radius, but off by factor f */
		dr = r - r_fit;		/* Radial misfit */
		rms += dr * dr;		/* Sum up the rms w.r.t. model */
	}
	p->fit = 100.0 * (1 - sqrt (rms / n) / major);	/* Our fit parameter */
	gmt_M_free (GMT, b);
}

#define bailout(code) {gmt_M_free_options (mode); return (code);}
#define Return(code) {Free_Ctrl (GMT, Ctrl); gmt_end_module (GMT, GMT_cpy); bailout (code);}

EXTERN_MSC int GMT_grdslice (void *V_API, int mode, void *args) {
	bool begin = true, inside, is_mercator = false, dump, has_foundation, *skip = NULL;

	int error, cs;

	unsigned int c, n_contours, n_edges, slc, n_inside, n_peaks = 0, n_slices = 0, n_foundations, index, n_skipped = 0, F_id;
	unsigned int geo = 0, *np = NULL, *edge = NULL;

	uint64_t ij, n, i, n_alloc = GMT_CHUNK, row, seg, dim[4] = {1, 1, 0, 0};
	int64_t ns;

	char header[GMT_LEN256] = {""}, *xname[2] = {"x", "lon"}, *yname[2] = {"y", "lat"};

	double aspect, cval, min, max, small, scale = 1.0, area, pos[3], *x_orig = NULL, *y_orig = NULL;
	double small_x, small_y, lon, lat = 0.0, min_area, max_area, merc_x0 = 0.0, merc_y0 = 0.0;
	double wesn[4], A[4], EigenValue[2], EigenVector[4], out[9], *x = NULL, *y = NULL, *contour = NULL;
	double wesn_m[4] = {GMT_IMG_MINLON, GMT_IMG_MAXLON, GMT_IMG_MINLAT_80, GMT_IMG_MAXLAT_80};

	struct GMT_GRID *G = NULL, *G_orig = NULL;
	struct GMT_DATASET *Center = NULL, *Slice = NULL, *Foundation = NULL, *Index = NULL;
	struct GMT_DATASEGMENT *S = NULL, *SI = NULL;
	struct GRDSLICE_SLICE *this_slice = NULL, *poly = NULL, **slice = NULL, **last = NULL;
	struct GRDSLICE_PEAK **peak = NULL, *this_peak = NULL;
	struct GRDSLICE_CTRL *Ctrl = NULL;
	struct GMT_CTRL *GMT = NULL, *GMT_cpy = NULL;
	struct GMT_OPTION *options = NULL;
	struct GMTAPI_CTRL *API = gmt_get_api_ptr (V_API);	/* Cast from void to GMTAPI_CTRL pointer */

	/*----------------------- Standard module initialization and parsing ----------------------*/

	if (API == NULL) return (GMT_NOT_A_SESSION);
	if (mode == GMT_MODULE_PURPOSE) return (usage (API, GMT_MODULE_PURPOSE));	/* Return the purpose of program */
	options = GMT_Create_Options (API, mode, args);	if (API->error) return (API->error);	/* Set or get option list */

	if ((error = gmt_report_usage (API, options, 0, usage)) != GMT_NOERROR) bailout (error);	/* Give usage if requested */

	/* Parse the command-line arguments */

	if ((GMT = gmt_init_module (API, THIS_MODULE_LIB, THIS_MODULE_CLASSIC_NAME, THIS_MODULE_KEYS, THIS_MODULE_NEEDS, NULL, &options, &GMT_cpy)) == NULL) bailout (API->error); /* Save current state */
	if (GMT_Parse_Common (API, THIS_MODULE_OPTIONS, options)) Return (API->error);
	Ctrl = New_Ctrl (GMT);	/* Allocate and initialize a new control structure */
	if ((error = parse (GMT, Ctrl, options)) != GMT_NOERROR) Return (error);

	/*---------------------------- This is the grdslice main code ----------------------------*/

	if (!strcmp (Ctrl->In.file,  "=")) {
		GMT_Report (API, GMT_MSG_ERROR, "Piping of grids not supported!\n");
		Return (EXIT_FAILURE);
	}

	GMT_Report (API, GMT_MSG_INFORMATION, "Allocate memory and read data file\n");

	gmt_M_memcpy (wesn, GMT->common.R.wesn, 4, double);	/* Current -R setting, if any (else it is 0/0/0/0 and ignored) */

	if ((G = GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_CONTAINER_ONLY, NULL, Ctrl->In.file, NULL)) == NULL) {
		Return (API->error);
	}

	if (gmt_M_is_subset (GMT, G->header, wesn)) {	/* Subset requested; make sure wesn matches header spacing since may have rough edges */
		if ((error = gmt_M_err_fail (GMT, gmt_adjust_loose_wesn (GMT, wesn, G->header), "")))
			Return (GMT_RUNTIME_ERROR);
	}

	/* Read in the grid, subset or not */
	if (GMT_Read_Data (API, GMT_IS_GRID, GMT_IS_FILE, GMT_IS_SURFACE, GMT_DATA_ONLY, wesn, Ctrl->In.file, G) == NULL) {
		Return (API->error);
	}

	if (gmt_M_is_geographic (GMT, GMT_IN)) geo = 1;	/* Geo means input grid is in geographic degrees */
	is_mercator = (strstr (G->header->remark, "Spherical Mercator Projected with -Jm1") != NULL);	/* Do we detect a Mercator grid? */

	if (is_mercator || geo) {	/* One or the other flavor of geographic data */
		/* Indirectly set up a plain Mercator projection for a spherical Earth with -Jm1i -R0/360/-lat/+lat, unless Cartesian data */
		double clon;	/* Determine a reasonable central meridian in steps of 90 */
		clon = 90.0 * ((GMT->common.R.active[RSET]) ? rint (0.5 * (wesn[XLO] + wesn[XHI] / 90.0)) : rint (0.5 * (G->header->wesn[XLO] + G->header->wesn[XHI] / 90.0)));
		wesn_m[XLO] = clon - 180.0;	wesn_m[XHI] = clon + 180.0;
		GMT->current.setting.proj_ellipsoid = gmt_get_ellipsoid (GMT, "Sphere");
		GMT->current.proj.units_pr_degree = true;
		GMT->current.proj.pars[0] = clon;
		GMT->current.proj.pars[1] = 0.0;
		GMT->current.proj.pars[2] = 1.0;
		GMT->current.proj.projection = GMT->current.proj.projection_GMT = GMT_MERCATOR;
		GMT->common.J.active = true;

		if (gmt_map_setup (GMT, wesn_m)) Return (GMT_PROJECTION_ERROR);

		if (is_mercator)	/* Get Mercator coordinates of our origin at (0,0) */
			gmt_geo_to_xy (GMT, 0.0, 0.0, &merc_x0, &merc_y0);
	}

	if (!(Ctrl->Z.scale == 1.0 && Ctrl->Z.offset == 0.0)) {	/* Must transform z grid */
		GMT_Report (API, GMT_MSG_INFORMATION, "Subtracting %g and multiplying grid by %g\n", Ctrl->Z.offset, Ctrl->Z.scale);
		G->header->z_min = (G->header->z_min - Ctrl->Z.offset) * Ctrl->Z.scale;
		G->header->z_max = (G->header->z_max - Ctrl->Z.offset) * Ctrl->Z.scale;
		if (Ctrl->Z.scale < 0.0) gmt_M_double_swap (G->header->z_min, G->header->z_max);
		/* Since gmt_scale_and_offset_f applies z' = z * scale + offset we must adjust Z.offset first: */
		Ctrl->Z.offset *= Ctrl->Z.scale;
		gmt_scale_and_offset_f (GMT, G->data, G->header->size, Ctrl->Z.scale, -Ctrl->Z.offset);
	}
	if (Ctrl->L.low  > G->header->z_min) G->header->z_min = Ctrl->L.low;	/* Possibly clip the z range */
	if (Ctrl->L.high < G->header->z_max) G->header->z_max = Ctrl->L.high;

	small = Ctrl->C.interval * 1.0e-6;	/* Noise threshold for contouring */

	small_x = 0.01 * G->header->inc[GMT_X];	small_y = 0.01 * G->header->inc[GMT_Y];	/* Noise in detecting closed contours */
	min_area = DBL_MAX;	max_area = -DBL_MAX;

	n_alloc = GMT_CHUNK;
	contour = gmt_M_memory (GMT, NULL, n_alloc, double);

	/* Set up contour intervals automatically from Ctrl->C.interval */

	min = floor (G->header->z_min / Ctrl->C.interval) * Ctrl->C.interval; if (min < G->header->z_min) min += Ctrl->C.interval;
	max = ceil  (G->header->z_max / Ctrl->C.interval) * Ctrl->C.interval; if (max > G->header->z_max) max -= Ctrl->C.interval;
	for (cs = irint (min/Ctrl->C.interval), n_contours = 0; cs <= irint (max/Ctrl->C.interval); cs++, n_contours++) {
		if (n_contours == n_alloc) {
			n_alloc += GMT_CHUNK;
			contour = gmt_M_memory (GMT, contour, n_alloc, double);
		}
		contour[n_contours] = cs * Ctrl->C.interval;
	}
	contour = gmt_M_memory (GMT, contour, n_contours, double);

	/* Because we are doing single-precision, we cannot subtract incrementally but must start with the
	 * original grid values and subtract the current contour value. */

	if ((G_orig = GMT_Duplicate_Data (API, GMT_IS_GRID, GMT_DUPLICATE_DATA, G)) == NULL) {
		gmt_M_free (GMT, contour);
		Return (GMT_RUNTIME_ERROR);
	}

	/* Get initial memory allocations for slices and peaks */
	slice = gmt_M_memory (GMT, NULL, n_contours, struct GRDSLICE_SLICE *);
	last  = gmt_M_memory (GMT, NULL, n_contours, struct GRDSLICE_SLICE *);
	peak  = gmt_M_memory (GMT, NULL, n_alloc,    struct GRDSLICE_PEAK *);

	/* Allocate edge array for contouring */
	edge = gmt_contour_edge_init (GMT, G->header, &n_edges);

	for (cs = (int)(n_contours-1); cs >= 0; cs--) {	/* For each contour value cval but starting from the top instead of base */
		c = (unsigned int)cs;

		/* Reset markers and set up new zero-contour */

		cval = contour[c];	/* Current slice contour */

		for (ij = 0; ij < G->header->size; ij++) {	/* Ensure G->data[ij] != 0.0 */
			G->data[ij] = G_orig->data[ij] - (gmt_grdfloat)cval;		/* If there are NaNs they will remain NaNs */
			if (G->data[ij] == 0.0) G->data[ij] += (gmt_grdfloat)small;	  /* There will be no actual zero-values, just -ve and +ve values */
		}

		/* Trace contours and skip all open contours */
		begin = true;
		while ((ns = gmt_contours (GMT, G, Ctrl->S.value, GMT->current.setting.interpolant, 0, edge, &begin, &x, &y)) > 0) {
			n = (uint64_t)ns;
			if (!(fabs (x[0] - x[n-1]) < small_x && fabs (y[0] - y[n-1]) < small_y)) {	/* Not a closed contour, skip */
				gmt_M_free (GMT, x);
				gmt_M_free (GMT, y);
				continue;
			}
			x[n-1] = x[0];	y[n-1] = y[0];	/* Close the contour exactly, so first and last point are identical */
			
			/* Allocate structure for the new contour slice */
			
			this_slice    = gmt_M_memory (GMT, NULL, 1, struct GRDSLICE_SLICE);
			this_slice->x = gmt_M_memory (GMT, NULL, n, double);
			this_slice->y = gmt_M_memory (GMT, NULL, n, double);

			
			/* Get signed area and point of max grid value inside contour from geographic, Mercator, or Cartesian coordinates */
			area = grdslice_centroid_area (GMT, x, y, n, NULL, NULL, geo, false, 0.0, pos);

			/* Below we want all coordinates to be projected Mercator (or the original Cartesian) except x_mean/y_mean which will be in degrees (unless Cartesian) */

			x_orig = gmt_M_memory (GMT, NULL, n, double);
			y_orig = gmt_M_memory (GMT, NULL, n, double);
			gmt_M_memcpy (x_orig, x, n, double);	/* Copy original grid coordinates */
			gmt_M_memcpy (y_orig, y, n, double);
			if (geo) {	/* Contour coordinates (x, y) are geographic (lon, lat), convert to Mercator */
				for (i = 0; i < n; i++)
					gmt_geo_to_xy (GMT, x[i], y[i], &this_slice->x[i], &this_slice->y[i]);
				gmt_M_memcpy (x, this_slice->x, n, double);	/* Copy polygon projected coordinates back to this array */
				gmt_M_memcpy (y, this_slice->y, n, double);
				this_slice->x_mean = pos[GMT_X]; this_slice->y_mean = pos[GMT_Y];
				gmt_geo_to_xy (GMT, this_slice->x_mean, this_slice->y_mean, &pos[GMT_X], &pos[GMT_Y]);
			}
			else if (is_mercator) {	/* All values are in Mercator units */
				gmt_M_memcpy (this_slice->x, x, n, double);	/* Copy polygon projected coordinates to this array */
				gmt_M_memcpy (this_slice->y, y, n, double);
				gmt_xy_to_geo (GMT, &this_slice->x_mean, &this_slice->y_mean, pos[GMT_X] + merc_x0, pos[GMT_Y] + merc_y0);
				scale = GMT->current.proj.DIST_KM_PR_DEG * cosd (this_slice->y_mean);	/* Needed to convert Mercator areas to km below */
			}
			else {	/* Cartesian: no conversions */
				gmt_M_memcpy (this_slice->x, x, n, double);	/* Copy polygon Cartesian coordinates to this array as well */
				gmt_M_memcpy (this_slice->y, y, n, double);
				this_slice->x_mean = pos[GMT_X]; this_slice->y_mean = pos[GMT_Y];
			}

			/* Here, everything is Mercator or Cartesian coordinates */
			/* Get information about this slice and store it in structure as projected (or Cartesian) coordinates */

			this_slice->xmin = this_slice->ymin = +DBL_MAX;
			this_slice->xmax = this_slice->ymax = -DBL_MAX;
			this_slice->n = n;
			this_slice->z = this_slice->cval = cval;
			this_slice->area = area * (scale * scale);	/* Now in km^2 unless for Cartesian grids */

			GMT_Report (API, GMT_MSG_INFORMATION, "Area = %g with scale = %g for z = %g [lat = %g]\n", area, scale, this_slice->z, lat);

			for (i = 0; i < n; i++) {	/* Determine extreme projected coordinate values */
				if (x[i] < this_slice->xmin) this_slice->xmin = this_slice->x[i];
				if (x[i] > this_slice->xmax) this_slice->xmax = this_slice->x[i];
				if (y[i] < this_slice->ymin) this_slice->ymin = this_slice->y[i];
				if (y[i] > this_slice->ymax) this_slice->ymax = this_slice->y[i];
			}

			/* Find orientation of major/minor axes and aspect ratio from reduced, projected x,y coordinates */

			grdslice_fit_ellipse (GMT, x, y, n, this_slice, pos, area);
			
			/* Update information of min/max area */
			if (this_slice->area > max_area) max_area = this_slice->area;
			if (this_slice->area < min_area) min_area = this_slice->area;
			
			n_slices++;
			
			/* Append this slice to end of slice array */
			
			if (last[c]) {	/* We already had at least one slice at this level - append it at end of the list */
				last[c]->next = this_slice;
				last[c] = last[c]->next;	/* Update the pointer to the last slice in this list */
			}
			else	/* First time - initialize slice and last to point to this slice */
				slice[c] = last[c] = this_slice;

			/* Determine if any previous (higher) level contours are inside the new (lower) contour */
			
			n_inside = 0;	/* Not inside anything so far */
			if (c < (n_contours-1) && slice[c+1]) {	/* There is a previous contour level and it has contours - must check for nesting */
				poly = slice[c+1];		/* First contour polygon for the previous contour level */
				while (poly) {			/* As long as there are more polygons */
					inside = true;		/* Set to true initially but this is usually reversed by one of two tests below: */
					if (poly->x[0] < this_slice->xmin || poly->x[0] > this_slice->xmax || poly->y[0] < this_slice->ymin || poly->y[0] > this_slice->ymax) {
						inside = false;		/* Outside polygon's extreme x/y-range */
					}
					else if (gmt_non_zero_winding (GMT, poly->x[0], poly->y[0], this_slice->x, this_slice->y, this_slice->n) < 2) {
						inside = false;		/* Inside rectangle, but still outside the polygon (returns 0 = outside, 1 on line (should not happen), 2 inside) */
					}
					if (inside) {	/* Here we know that the previous contour slice is inside the new one so we point to the new slice as "downstream" from the old one */
						poly->down = this_slice;
						poly->shared++;	/* Number of peaks sharing this slice */
						n_inside++;
					}
					poly = poly->next;	/* Go to next polygon in the list of contours at the previous level */
				}
			}
			if (n_inside == 0) {	/* No previous contours contained by this contour - initialize a new peak location at the center of this slice */
				/* Must revise the peak location and amplitude via BCR brute force search */
				(void)grdslice_centroid_area (GMT, x_orig, y_orig, n, this_slice, G_orig, geo, true, Ctrl->Q.factor, pos);
				if (is_mercator)
					gmt_xy_to_geo (GMT, &this_slice->x_mean, &this_slice->y_mean, pos[GMT_X] + merc_x0, pos[GMT_Y] + merc_y0);
				else
					this_slice->x_mean = pos[GMT_X], this_slice->y_mean = pos[GMT_Y];
				this_slice->z = pos[GMT_Z];
				this_peak = gmt_M_memory (GMT, NULL, 1, struct GRDSLICE_PEAK);
				this_peak->x = this_slice->x_mean;	/* Just use mean location for now - perhaps later choose the actual grid maximum */
				this_peak->y = this_slice->y_mean;
				this_peak->cval = this_slice->cval;		/* Contour value */
				this_peak->z = this_slice->z;		/* The local high value in the grid */
				this_peak->start = this_slice;		/* This is the top slice in the stack below this point */
				this_peak->id = n_peaks;			/* This is the unique peak ID */
				peak[n_peaks] = this_peak;			/* Add peak to peak array */
				n_peaks++;
				if (n_peaks == n_alloc) {	/* Get more memory */
					n_alloc += GMT_CHUNK;
					peak = gmt_M_memory (GMT, peak, n_alloc, struct GRDSLICE_PEAK *);
				}
			}
			gmt_M_free (GMT, x);	gmt_M_free (GMT, y);	/* Free original memory returned by gmt_contours */
			gmt_M_free (GMT, x_orig);	gmt_M_free (GMT, y_orig);	/* Free original x/y coordinates */

		}
		GMT_Report (API, GMT_MSG_INFORMATION, "Tracing the %8.2f contour: # of slices: %6d # of peaks: %6d\n", cval, n_slices, n_peaks);
	}
	gmt_M_free (GMT, edge);
	gmt_M_free (GMT, contour);

	if (is_mercator) geo = 1;	/* Since all coordinates have been converted by now */

	peak = gmt_M_memory (GMT, peak, n_peaks, struct GRDSLICE_PEAK *);

	if (Ctrl->D.active) {	/* Calculate distances between peaks in km and reject smaller peaks within set distance from a larger peak */
		int which;
		double dist;
		skip = gmt_M_memory (GMT, NULL, n_peaks, bool);
		for (c = 0; c < n_peaks; c++) {	/* For each peak */
			if (skip[c]) continue;	/* Already flagged to be skipped */
			for (i = c+1; i < n_peaks; i++) {	/* For each other peak */
				if (skip[i]) continue;	/* Already flagged to be skipped */
				if (geo)	/* Compute distances between peaks in km */
					dist = 0.001 * gmt_great_circle_dist_meter (GMT, peak[c]->start->x_mean, peak[c]->start->y_mean, peak[i]->start->x_mean, peak[i]->start->y_mean);
				else 	/* Cartesian distances */
					dist = hypot (peak[c]->start->x_mean - peak[c]->start->y_mean, peak[i]->start->x_mean - peak[i]->start->y_mean);
				if (dist < Ctrl->D.cutoff) {	/* These two peaks are too close, eliminate the smaller one */
					which = (peak[c]->start->z >= peak[i]->start->z) ? i : c;
					skip[which] = true;
				}
			}
		}
		for (c = n_skipped = 0; c < n_peaks; c++)	/* Count skipped peaks */
			if (skip[c]) n_skipped++;
	}

	/* Count slices and foundations */
	np = gmt_M_memory (GMT, NULL, n_contours, unsigned int);	/* Store number of slices per peak */
	for (c = n_slices = n_foundations = 0; c < n_contours; c++) {	/* For each peak */
		poly = slice[c];		/* First contour polygon at this contour level */
		np[c] = 0;	/* Number of slices per level c */
		while (poly) {			/* As long as there are more polygons at this level */
			if (Ctrl->T.active && doubleAlmostEqualZero (poly->z, Ctrl->T.cutoff) && poly->area >= Ctrl->A.cutoff) {
				poly->F_id = n_foundations;
				n_foundations++;
			}
			np[c]++;
			poly = poly->next;	/* Go to next polygon in this contour list */
		}
		n_slices += np[c];		/* Total number of slices found so far */
	}

	GMT_Report (API, GMT_MSG_INFORMATION, "Peaks found: %d\n", n_peaks);
	if (n_skipped) GMT_Report (API, GMT_MSG_INFORMATION, "Peaks eliminated due to -D minimum separation: %d\n", n_skipped);
	GMT_Report (API, GMT_MSG_INFORMATION, "Slices: %d\n", n_slices);
	if (n_foundations == 0) {	/* Found nothing to write home about */
		if (Ctrl->F.active)
			GMT_Report (API, GMT_MSG_WARNING, "Option -F: No peak foundations detected despite -T being set\n");
		if (Ctrl->I.active)
			GMT_Report (API, GMT_MSG_WARNING, "Option -I: No peak indices detected despite -T being set\n");
		Ctrl->F.active = Ctrl->I.active = false;	/* To prevent any attempt to write nothing below */
	}
	else
		GMT_Report (API, GMT_MSG_INFORMATION, "Indices/Foundations: %d\n", n_foundations);
	if (geo) gmt_set_geographic (GMT, GMT_OUT);	/* We wish to write geographic coordinates */
	gmt_set_tableheader (GMT, GMT_OUT, true);

	dim[GMT_SEG] = n_peaks - n_skipped;	dim[GMT_COL] = 9;	/* Default output format for centers */

	if ((Center = GMT_Create_Data (API, GMT_IS_DATASET, GMT_IS_LINE, 0, dim, NULL, NULL, 0, 0, NULL)) == NULL) {
		Return (API->error);
	}
	for (c = seg = 0; c < n_peaks; c++) {	/* For each peak */
		if (Ctrl->D.active && skip[c])	continue;	/* Skipped peak */
		poly = peak[c]->start;		/* First contour polygon originating form this peak */
		row = 0;
		has_foundation = false;
		while (poly) {			/* As long as there are more polygons at this level */
			if (Ctrl->T.active && doubleAlmostEqualZero (poly->z, Ctrl->T.cutoff) && poly->area >= Ctrl->A.cutoff) {
				F_id = poly->F_id;
				has_foundation = true;
			}
			row++;
			poly = poly->down;		/* Go to next polygon in this contour list */
		}
		if (has_foundation) {
			poly = peak[c]->start;		/* First contour polygon originating form this peak */
			while (poly) {			/* As long as there are more polygons at this level */
				poly->F_id = F_id;
				poly = poly->down;	/* Go to next polygon above in this contour list */
			}
		}
		sprintf (header, "-Z%g -L%d", peak[c]->cval, peak[c]->id);
		S = GMT_Alloc_Segment (API, GMT_NO_STRINGS, row, Center->n_columns, header, Center->table[0]->segment[seg]);
		poly = peak[c]->start;		/* Back to first contour again */
		row = 0;
		while (poly) {			/* As long as there are more polygons at this level. Format: x y z id area major minor azimuth fit */
			S->data[GMT_X][row] = poly->x_mean;	S->data[GMT_Y][row] = poly->y_mean;	S->data[GMT_Z][row] = poly->z;
			S->data[3][row] = peak[c]->id; S->data[4][row] = poly->area; S->data[5][row] = poly->major;
			S->data[6][row] = poly->minor; S->data[7][row] = poly->azimuth; S->data[8][row] = poly->fit;
			poly = poly->down;		/* Go to next polygon in this contour list */
			row++;
		}
		seg++;
	}
	sprintf (header, "%s\t%s\tz\tid\tarea\tmajor\tminor\tazimuth\tfit", xname[is_mercator], yname[is_mercator]);
	if (GMT_Set_Comment (API, GMT_IS_DATASET, GMT_COMMENT_IS_COLNAMES, header, Center)) {	/* Not working yet, gives command line instead */
		Return (API->error);
	}
	if (GMT_Write_Data (API, GMT_IS_DATASET, GMT_IS_FILE, GMT_IS_LINE, GMT_WRITE_NORMAL, NULL, Ctrl->Out.file, Center) != GMT_NOERROR) {
		Return (API->error);
	}
	if (Ctrl->Out.file)
		GMT_Report (API, GMT_MSG_INFORMATION, "%d centers written to %s\n", seg, Ctrl->Out.file);
	else
		GMT_Report (API, GMT_MSG_INFORMATION, "%d centers written to stdout\n", seg);

	/* Optional outputs, if any */
	GMT_Report (API, GMT_MSG_INFORMATION, "Slices: %d\n", n_slices);

	if (Ctrl->E.active || Ctrl->F.active || Ctrl->I.active) {
		if (Ctrl->E.active) {	/* Must write all slices to one table */
			dim[GMT_SEG] = n_slices;	dim[GMT_COL] = 3;
			if ((Slice = GMT_Create_Data (API, GMT_IS_DATASET, GMT_IS_POLY, 0, dim, NULL, NULL, 0, 0, NULL)) == NULL) {
				Return (API->error);
			}
		}
		if (Ctrl->F.active) {	/* Must write foundation polygons, one segment per foundation */
			dim[GMT_SEG] = n_foundations;	dim[GMT_COL] = 3;
			if ((Foundation = GMT_Create_Data (API, GMT_IS_DATASET, GMT_IS_POLY, 0, dim, NULL, NULL, 0, 0, NULL)) == NULL) {
				Return (API->error);
			}
		}
		if (Ctrl->I.active) {	/* Must write index points to a single segment in the table. Format: lon lat x y id */
			dim[GMT_SEG] = 1;	dim[GMT_ROW] = n_foundations;	dim[GMT_COL] = 5;
			if ((Index = GMT_Create_Data (API, GMT_IS_DATASET, GMT_IS_POINT, 0, dim, NULL, NULL, 0, 0, NULL)) == NULL) {
				Return (API->error);
			}
			SI = Index->table[0]->segment[0];	/* Only a single segment of indices */
		}

		for (c = seg = slc = index = 0; c < n_contours; c++) {	/* For each slice */
			poly = slice[c];		/* First contour polygon at this contour level */
			while (poly) {			/* As long as there are more polygons at this level */
				if (geo)
					gmt_geo_to_xy (GMT, poly->x_mean, poly->y_mean, &out[GMT_X], &out[GMT_Y]);
				else
					out[GMT_X] = poly->x_mean, out[GMT_Y] = poly->y_mean;
				if (Ctrl->T.active && doubleAlmostEqualZero (poly->z, Ctrl->T.cutoff) && poly->area >= Ctrl->A.cutoff) {
					if (Ctrl->I.active) {	/* Fill out one index record with format: lon lat x y id */
						SI->data[GMT_X][index] = poly->x_mean;
						SI->data[GMT_Y][index] = poly->y_mean;
						SI->data[2][index] = out[GMT_X];
						SI->data[3][index] = out[GMT_Y];
						SI->data[4][index] = index;
					}
					if (Ctrl->F.active) {	/* Write out the foundation contour. Format: x y z */
						sprintf (header, "%d -Z%g -L%g -N%d -S%g/%g/%g/%g/%g/%g", index, poly->area, poly->cval, poly->shared, out[GMT_X], out[GMT_Y], poly->azimuth, poly->major, poly->minor, poly->fit);
						S = GMT_Alloc_Segment (API, GMT_NO_STRINGS, poly->n, Foundation->n_columns, header, Foundation->table[0]->segment[index]);
						for (row = 0; row < poly->n; row++) {
							if (geo)	/* Get lon, lat */
								gmt_xy_to_geo (GMT, &S->data[GMT_X][row], &S->data[GMT_Y][row], poly->x[row] + merc_x0, poly->y[row] + merc_y0);
							else
								S->data[GMT_X][row] = poly->x[row], S->data[GMT_Y][row] = poly->y[row];
							S->data[GMT_Z][row] = poly->z;
						}
					}
					index++;
				}
				if (Ctrl->E.active) {	/* Write out all the polygon perimeters. Format: x y z  */
					sprintf (header, "-Z%g -L%g -N%d -S%g/%g/%g/%g/%g/%g -I%d", poly->area, poly->cval, poly->shared, out[GMT_X], out[GMT_Y], poly->azimuth, poly->major, poly->minor, poly->fit, poly->F_id);
					S = GMT_Alloc_Segment (API, GMT_NO_STRINGS, poly->n, Slice->n_columns, header, Slice->table[0]->segment[slc]);
					for (row = 0; row < poly->n; row++) {
						if (geo)	/* Get lon, lat */
							gmt_xy_to_geo (GMT, &S->data[GMT_X][row], &S->data[GMT_Y][row], poly->x[row] + merc_x0, poly->y[row] + merc_y0);
						else
							S->data[GMT_X][row] = poly->x[row], S->data[GMT_Y][row] = poly->y[row];
						S->data[GMT_Z][row] = poly->z;
					}
					slc++;
				}
				poly = poly->next;		/* Go to next polygon in this contour list */
			}
		}
		if (Ctrl->E.active) {
			sprintf (header, "%s\t%s\tz", xname[is_mercator], yname[is_mercator]);
			if (GMT_Set_Comment (API, GMT_IS_DATASET, GMT_COMMENT_IS_COLNAMES, header, Slice)) {	/* Not working yet, gives command line instead */
				Return (API->error);
			}
			if (GMT_Write_Data (API, GMT_IS_DATASET, GMT_IS_FILE, GMT_IS_POLY, GMT_WRITE_NORMAL, NULL, Ctrl->E.file, Slice) != GMT_NOERROR) {
				Return (API->error);
			}
			GMT_Report (API, GMT_MSG_INFORMATION, "%d slices written to %s\n", slc, Ctrl->E.file);
		}
		if (Ctrl->F.active) {
			sprintf (header, "%s\t%s\tz", xname[is_mercator], yname[is_mercator]);
			if (GMT_Set_Comment (API, GMT_IS_DATASET, GMT_COMMENT_IS_COLNAMES, header, Foundation)) {	/* Not working yet, gives command line instead */
				Return (API->error);
			}
			if (GMT_Write_Data (API, GMT_IS_DATASET, GMT_IS_FILE, GMT_IS_POLY, GMT_WRITE_NORMAL, NULL, Ctrl->F.file, Foundation) != GMT_NOERROR) {
				Return (API->error);
			}
			GMT_Report (API, GMT_MSG_INFORMATION, "%d foundations written to %s\n", index, Ctrl->F.file);
		}
		if (Ctrl->I.active) {
			sprintf (header, "%s\t%s\tx\ty\tID", xname[is_mercator], yname[is_mercator]);
			if (GMT_Set_Comment (API, GMT_IS_DATASET, GMT_COMMENT_IS_COLNAMES, header, Index)) {	/* Not working yet, gives command line instead */
				Return (API->error);
			}
			if (GMT_Write_Data (API, GMT_IS_DATASET, GMT_IS_FILE, GMT_IS_POINT, GMT_WRITE_NORMAL, NULL, Ctrl->I.file, Index) != GMT_NOERROR) {
				Return (API->error);
			}
			GMT_Report (API, GMT_MSG_INFORMATION, "%d indices written to %s\n", index, Ctrl->I.file);
		}
	}

	/* Free memory */

	for (c = 0; c < n_peaks; c++) {	/* For each peak */
		if (peak[c]) gmt_M_free (GMT, peak[c]);
	}
	gmt_M_free (GMT, peak);

	for (c = 0; c < n_contours; c++) {	/* For each slice */
		poly = slice[c];
		while (poly) {	/* As long as there are more polygons at this level */
			this_slice = poly;
			poly = poly->next;	/* Go to next polygon in this contour list */
			gmt_M_free (GMT, this_slice->x);
			gmt_M_free (GMT, this_slice->y);
			gmt_M_free (GMT, this_slice);
		}
	}
	gmt_M_free (GMT, slice);
	gmt_M_free (GMT, last);
	gmt_M_free (GMT, np);

	GMT_Report (API, GMT_MSG_INFORMATION, "Done, min/max areas: %g %g\n", min_area, max_area);

	Return (GMT_NOERROR);
}
