
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_5248935579131712941) {
   out_5248935579131712941[0] = delta_x[0] + nom_x[0];
   out_5248935579131712941[1] = delta_x[1] + nom_x[1];
   out_5248935579131712941[2] = delta_x[2] + nom_x[2];
   out_5248935579131712941[3] = delta_x[3] + nom_x[3];
   out_5248935579131712941[4] = delta_x[4] + nom_x[4];
   out_5248935579131712941[5] = delta_x[5] + nom_x[5];
   out_5248935579131712941[6] = delta_x[6] + nom_x[6];
   out_5248935579131712941[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_2549868245409028765) {
   out_2549868245409028765[0] = -nom_x[0] + true_x[0];
   out_2549868245409028765[1] = -nom_x[1] + true_x[1];
   out_2549868245409028765[2] = -nom_x[2] + true_x[2];
   out_2549868245409028765[3] = -nom_x[3] + true_x[3];
   out_2549868245409028765[4] = -nom_x[4] + true_x[4];
   out_2549868245409028765[5] = -nom_x[5] + true_x[5];
   out_2549868245409028765[6] = -nom_x[6] + true_x[6];
   out_2549868245409028765[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_6426779444371049680) {
   out_6426779444371049680[0] = 1.0;
   out_6426779444371049680[1] = 0.0;
   out_6426779444371049680[2] = 0.0;
   out_6426779444371049680[3] = 0.0;
   out_6426779444371049680[4] = 0.0;
   out_6426779444371049680[5] = 0.0;
   out_6426779444371049680[6] = 0.0;
   out_6426779444371049680[7] = 0.0;
   out_6426779444371049680[8] = 0.0;
   out_6426779444371049680[9] = 1.0;
   out_6426779444371049680[10] = 0.0;
   out_6426779444371049680[11] = 0.0;
   out_6426779444371049680[12] = 0.0;
   out_6426779444371049680[13] = 0.0;
   out_6426779444371049680[14] = 0.0;
   out_6426779444371049680[15] = 0.0;
   out_6426779444371049680[16] = 0.0;
   out_6426779444371049680[17] = 0.0;
   out_6426779444371049680[18] = 1.0;
   out_6426779444371049680[19] = 0.0;
   out_6426779444371049680[20] = 0.0;
   out_6426779444371049680[21] = 0.0;
   out_6426779444371049680[22] = 0.0;
   out_6426779444371049680[23] = 0.0;
   out_6426779444371049680[24] = 0.0;
   out_6426779444371049680[25] = 0.0;
   out_6426779444371049680[26] = 0.0;
   out_6426779444371049680[27] = 1.0;
   out_6426779444371049680[28] = 0.0;
   out_6426779444371049680[29] = 0.0;
   out_6426779444371049680[30] = 0.0;
   out_6426779444371049680[31] = 0.0;
   out_6426779444371049680[32] = 0.0;
   out_6426779444371049680[33] = 0.0;
   out_6426779444371049680[34] = 0.0;
   out_6426779444371049680[35] = 0.0;
   out_6426779444371049680[36] = 1.0;
   out_6426779444371049680[37] = 0.0;
   out_6426779444371049680[38] = 0.0;
   out_6426779444371049680[39] = 0.0;
   out_6426779444371049680[40] = 0.0;
   out_6426779444371049680[41] = 0.0;
   out_6426779444371049680[42] = 0.0;
   out_6426779444371049680[43] = 0.0;
   out_6426779444371049680[44] = 0.0;
   out_6426779444371049680[45] = 1.0;
   out_6426779444371049680[46] = 0.0;
   out_6426779444371049680[47] = 0.0;
   out_6426779444371049680[48] = 0.0;
   out_6426779444371049680[49] = 0.0;
   out_6426779444371049680[50] = 0.0;
   out_6426779444371049680[51] = 0.0;
   out_6426779444371049680[52] = 0.0;
   out_6426779444371049680[53] = 0.0;
   out_6426779444371049680[54] = 1.0;
   out_6426779444371049680[55] = 0.0;
   out_6426779444371049680[56] = 0.0;
   out_6426779444371049680[57] = 0.0;
   out_6426779444371049680[58] = 0.0;
   out_6426779444371049680[59] = 0.0;
   out_6426779444371049680[60] = 0.0;
   out_6426779444371049680[61] = 0.0;
   out_6426779444371049680[62] = 0.0;
   out_6426779444371049680[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_8577481857843264736) {
   out_8577481857843264736[0] = state[0];
   out_8577481857843264736[1] = state[1];
   out_8577481857843264736[2] = state[2];
   out_8577481857843264736[3] = state[3];
   out_8577481857843264736[4] = state[4];
   out_8577481857843264736[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_8577481857843264736[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_8577481857843264736[7] = state[7];
}
void F_fun(double *state, double dt, double *out_2862380066592312672) {
   out_2862380066592312672[0] = 1;
   out_2862380066592312672[1] = 0;
   out_2862380066592312672[2] = 0;
   out_2862380066592312672[3] = 0;
   out_2862380066592312672[4] = 0;
   out_2862380066592312672[5] = 0;
   out_2862380066592312672[6] = 0;
   out_2862380066592312672[7] = 0;
   out_2862380066592312672[8] = 0;
   out_2862380066592312672[9] = 1;
   out_2862380066592312672[10] = 0;
   out_2862380066592312672[11] = 0;
   out_2862380066592312672[12] = 0;
   out_2862380066592312672[13] = 0;
   out_2862380066592312672[14] = 0;
   out_2862380066592312672[15] = 0;
   out_2862380066592312672[16] = 0;
   out_2862380066592312672[17] = 0;
   out_2862380066592312672[18] = 1;
   out_2862380066592312672[19] = 0;
   out_2862380066592312672[20] = 0;
   out_2862380066592312672[21] = 0;
   out_2862380066592312672[22] = 0;
   out_2862380066592312672[23] = 0;
   out_2862380066592312672[24] = 0;
   out_2862380066592312672[25] = 0;
   out_2862380066592312672[26] = 0;
   out_2862380066592312672[27] = 1;
   out_2862380066592312672[28] = 0;
   out_2862380066592312672[29] = 0;
   out_2862380066592312672[30] = 0;
   out_2862380066592312672[31] = 0;
   out_2862380066592312672[32] = 0;
   out_2862380066592312672[33] = 0;
   out_2862380066592312672[34] = 0;
   out_2862380066592312672[35] = 0;
   out_2862380066592312672[36] = 1;
   out_2862380066592312672[37] = 0;
   out_2862380066592312672[38] = 0;
   out_2862380066592312672[39] = 0;
   out_2862380066592312672[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_2862380066592312672[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_2862380066592312672[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2862380066592312672[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2862380066592312672[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_2862380066592312672[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_2862380066592312672[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_2862380066592312672[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_2862380066592312672[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_2862380066592312672[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_2862380066592312672[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2862380066592312672[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2862380066592312672[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_2862380066592312672[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_2862380066592312672[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_2862380066592312672[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2862380066592312672[56] = 0;
   out_2862380066592312672[57] = 0;
   out_2862380066592312672[58] = 0;
   out_2862380066592312672[59] = 0;
   out_2862380066592312672[60] = 0;
   out_2862380066592312672[61] = 0;
   out_2862380066592312672[62] = 0;
   out_2862380066592312672[63] = 1;
}
void h_25(double *state, double *unused, double *out_7097572299460496017) {
   out_7097572299460496017[0] = state[6];
}
void H_25(double *state, double *unused, double *out_1576780616160605272) {
   out_1576780616160605272[0] = 0;
   out_1576780616160605272[1] = 0;
   out_1576780616160605272[2] = 0;
   out_1576780616160605272[3] = 0;
   out_1576780616160605272[4] = 0;
   out_1576780616160605272[5] = 0;
   out_1576780616160605272[6] = 1;
   out_1576780616160605272[7] = 0;
}
void h_24(double *state, double *unused, double *out_4763024709455486134) {
   out_4763024709455486134[0] = state[4];
   out_4763024709455486134[1] = state[5];
}
void H_24(double *state, double *unused, double *out_8247603600142628412) {
   out_8247603600142628412[0] = 0;
   out_8247603600142628412[1] = 0;
   out_8247603600142628412[2] = 0;
   out_8247603600142628412[3] = 0;
   out_8247603600142628412[4] = 1;
   out_8247603600142628412[5] = 0;
   out_8247603600142628412[6] = 0;
   out_8247603600142628412[7] = 0;
   out_8247603600142628412[8] = 0;
   out_8247603600142628412[9] = 0;
   out_8247603600142628412[10] = 0;
   out_8247603600142628412[11] = 0;
   out_8247603600142628412[12] = 0;
   out_8247603600142628412[13] = 1;
   out_8247603600142628412[14] = 0;
   out_8247603600142628412[15] = 0;
}
void h_26(double *state, double *unused, double *out_6131103443843024056) {
   out_6131103443843024056[0] = state[7];
}
void H_26(double *state, double *unused, double *out_220313274239780944) {
   out_220313274239780944[0] = 0;
   out_220313274239780944[1] = 0;
   out_220313274239780944[2] = 0;
   out_220313274239780944[3] = 0;
   out_220313274239780944[4] = 0;
   out_220313274239780944[5] = 0;
   out_220313274239780944[6] = 0;
   out_220313274239780944[7] = 1;
}
void h_27(double *state, double *unused, double *out_3832451680646436156) {
   out_3832451680646436156[0] = state[3];
}
void H_27(double *state, double *unused, double *out_7050806772081510156) {
   out_7050806772081510156[0] = 0;
   out_7050806772081510156[1] = 0;
   out_7050806772081510156[2] = 0;
   out_7050806772081510156[3] = 1;
   out_7050806772081510156[4] = 0;
   out_7050806772081510156[5] = 0;
   out_7050806772081510156[6] = 0;
   out_7050806772081510156[7] = 0;
}
void h_29(double *state, double *unused, double *out_6988761827225997278) {
   out_6988761827225997278[0] = state[1];
}
void H_29(double *state, double *unused, double *out_3791518996948079736) {
   out_3791518996948079736[0] = 0;
   out_3791518996948079736[1] = 1;
   out_3791518996948079736[2] = 0;
   out_3791518996948079736[3] = 0;
   out_3791518996948079736[4] = 0;
   out_3791518996948079736[5] = 0;
   out_3791518996948079736[6] = 0;
   out_3791518996948079736[7] = 0;
}
void h_28(double *state, double *unused, double *out_6180597138765863042) {
   out_6180597138765863042[0] = state[5];
   out_6180597138765863042[1] = state[6];
}
void H_28(double *state, double *unused, double *out_7686849636252634942) {
   out_7686849636252634942[0] = 0;
   out_7686849636252634942[1] = 0;
   out_7686849636252634942[2] = 0;
   out_7686849636252634942[3] = 0;
   out_7686849636252634942[4] = 0;
   out_7686849636252634942[5] = 1;
   out_7686849636252634942[6] = 0;
   out_7686849636252634942[7] = 0;
   out_7686849636252634942[8] = 0;
   out_7686849636252634942[9] = 0;
   out_7686849636252634942[10] = 0;
   out_7686849636252634942[11] = 0;
   out_7686849636252634942[12] = 0;
   out_7686849636252634942[13] = 0;
   out_7686849636252634942[14] = 1;
   out_7686849636252634942[15] = 0;
}
}

extern "C"{
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
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;
  
  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);
  
  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H); 
  
  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();
   

    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;
  
  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);
 
  // update cov 
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
