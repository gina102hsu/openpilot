/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_392995883448265057);
void inv_err_fun(double *nom_x, double *true_x, double *out_5772860993230091607);
void H_mod_fun(double *state, double *out_9187240664393523899);
void f_fun(double *state, double dt, double *out_3582577477181916368);
void F_fun(double *state, double dt, double *out_2285739652126540264);
void h_3(double *state, double *unused, double *out_2530045924594199149);
void H_3(double *state, double *unused, double *out_6382007605522758295);
void h_4(double *state, double *unused, double *out_8996411428999236214);
void H_4(double *state, double *unused, double *out_4387151028043857657);
void h_9(double *state, double *unused, double *out_3337231452119161023);
void H_9(double *state, double *unused, double *out_2986822082298273303);
void h_10(double *state, double *unused, double *out_8408478419395825874);
void H_10(double *state, double *unused, double *out_1162930201856238328);
void h_12(double *state, double *unused, double *out_8258375951659881901);
void H_12(double *state, double *unused, double *out_7534016994323241433);
void h_13(double *state, double *unused, double *out_614960214416859041);
void H_13(double *state, double *unused, double *out_5080997235599590610);
void h_14(double *state, double *unused, double *out_3337231452119161023);
void H_14(double *state, double *unused, double *out_2986822082298273303);
void h_19(double *state, double *unused, double *out_1038943756963724727);
void H_19(double *state, double *unused, double *out_7089476489701242269);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);