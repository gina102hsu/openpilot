/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_5248935579131712941);
void inv_err_fun(double *nom_x, double *true_x, double *out_2549868245409028765);
void H_mod_fun(double *state, double *out_6426779444371049680);
void f_fun(double *state, double dt, double *out_8577481857843264736);
void F_fun(double *state, double dt, double *out_2862380066592312672);
void h_25(double *state, double *unused, double *out_7097572299460496017);
void H_25(double *state, double *unused, double *out_1576780616160605272);
void h_24(double *state, double *unused, double *out_4763024709455486134);
void H_24(double *state, double *unused, double *out_8247603600142628412);
void h_26(double *state, double *unused, double *out_6131103443843024056);
void H_26(double *state, double *unused, double *out_220313274239780944);
void h_27(double *state, double *unused, double *out_3832451680646436156);
void H_27(double *state, double *unused, double *out_7050806772081510156);
void h_29(double *state, double *unused, double *out_6988761827225997278);
void H_29(double *state, double *unused, double *out_3791518996948079736);
void h_28(double *state, double *unused, double *out_6180597138765863042);
void H_28(double *state, double *unused, double *out_7686849636252634942);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
