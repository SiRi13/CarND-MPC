#ifndef MPC_H
#define MPC_H

#define VARS_BOUND 1.0e19;
#define THROTTLE_BOUND 1.0;

#include <vector>
#include "Eigen-3.3/Eigen/Core"

using namespace std;

class MPC {
 public:
  // TODO: Set the timestep length and duration
  static constexpr unsigned N = 20;
  static constexpr double dt = 0.1;

  // This value assumes the model presented in the classroom is used.
  //
  // It was obtained by measuring the radius formed by running the vehicle in
  // the simulator around in a circle with a constant steering angle and
  // velocity on a flat terrain.
  //
  // Lf was tuned until the the radius formed by the simulating the model
  // presented in the classroom matched the previous radius.
  //
  // This is the length from front to CoG that has a similar radius.
  static constexpr double Lf = 2.67;

  static constexpr double REF_V = 40;

  static constexpr double FAC_CTE = 2;
  static constexpr double FAC_EPSI = 20;
  static constexpr double FAC_V = 1;
  static constexpr double FAC_A = 20;
  static constexpr double FAC_DELTA = 100000;
  static constexpr double FAC_A_2 = 1;
  static constexpr double FAC_DELTA_2 = 1;

  // The solver takes all the state variables and actuator
  // variables in a singular vector. Thus, we should to establish
  // when one variable starts and another ends to make our lifes easier.
  size_t x_start = 0;
  size_t y_start = x_start + N;
  size_t psi_start = y_start + N;
  size_t v_start = psi_start + N;
  size_t cte_start = v_start + N;
  size_t epsi_start = cte_start + N;
  size_t delta_start = epsi_start + N;
  size_t a_start = delta_start + N - 1;

  vector<double> mpc_x, mpc_y;

  MPC();

  virtual ~MPC();

  // Solve the model given an initial state and polynomial coefficients.
  // Return the first actuatotions.
  vector<double> Solve(Eigen::VectorXd state, Eigen::VectorXd coeffs);
};

#endif /* MPC_H */
