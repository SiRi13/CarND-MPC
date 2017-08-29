#include <cppad/cppad.hpp>
#include <cppad/ipopt/solve.hpp>
#include "Eigen-3.3/Eigen/Core"
#include "MPC.h"

using CppAD::AD;

class FG_eval {
 public:
  // Fitted polynomial coefficients
  Eigen::VectorXd coeffs;

  MPC mpc;

  FG_eval(MPC &mpc, Eigen::VectorXd coeffs) {
    this->coeffs = coeffs;
    this->mpc = mpc;
  }

  typedef CPPAD_TESTVECTOR(AD<double>) ADvector;
  void operator()(ADvector &fg, const ADvector &vars) {
    // implement MPC
    // `fg` a vector of the cost constraints, `vars` is a vector of variable
    // values (state & actuators)
    // NOTE: You'll probably go back and forth between this function and
    // the Solver function below.

    // The cost stored is the first element of `fg`.
    // Any additions to the cost should be added to `fg[0]`.
    fg[0] = 0;

    // Reference State Cost
    // Define the cost related the reference state and
    // any anything you think may be beneficial.
    unsigned t = 0;
    for (t = 0; t < MPC::N; ++t) {
      // TODO: tune ref state
      fg[0] +=
          MPC::FAC_CTE * CppAD::pow(vars[mpc.cte_start + t] - MPC::REF_CTE, 2);
      fg[0] += MPC::FAC_EPSI *
               CppAD::pow(vars[mpc.epsi_start + t] - MPC::REF_EPSI, 2);
      fg[0] += MPC::FAC_V * CppAD::pow(vars[mpc.v_start + t] - MPC::REF_V, 2);

      if (t < MPC::N - 1) {
        fg[0] += MPC::FAC_DELTA * CppAD::pow(vars[mpc.delta_start + t], 2);
        fg[0] += MPC::FAC_A * CppAD::pow(vars[mpc.a_start + t], 2);
      }

      if (t < MPC::N - 2) {
        fg[0] += MPC::FAC_DELTA_2 * CppAD::pow(vars[mpc.delta_start + t + 1] -
                                                   vars[mpc.delta_start + t],
                                               2);
        fg[0] +=
            MPC::FAC_A_2 *
            CppAD::pow(vars[mpc.a_start + t + 1] - vars[mpc.a_start + t], 2);
      }
    }

    //
    // Setup Constraints
    //
    // NOTE: In this section you'll setup the model constraints.

    // Initial constraints
    //
    // We add 1 to each of the starting indices due to cost being located at
    // index 0 of `fg`.
    // This bumps up the position of all the other values.
    fg[1 + mpc.x_start] = vars[mpc.x_start];
    fg[1 + mpc.y_start] = vars[mpc.y_start];
    fg[1 + mpc.psi_start] = vars[mpc.psi_start];
    fg[1 + mpc.v_start] = vars[mpc.v_start];
    fg[1 + mpc.cte_start] = vars[mpc.cte_start];
    fg[1 + mpc.epsi_start] = vars[mpc.epsi_start];

    // The rest of the constraints
    for (t = 1; t < MPC::N; ++t) {
      // at time t
      AD<double> x0 = vars[mpc.x_start + t - 1];
      AD<double> y0 = vars[mpc.y_start + t - 1];
      AD<double> psi0 = vars[mpc.psi_start + t - 1];
      AD<double> v0 = vars[mpc.v_start + t - 1];
      AD<double> cte0 = vars[mpc.cte_start + t - 1];
      AD<double> ePsi0 = vars[mpc.epsi_start + t - 1];

      AD<double> delta0 = vars[mpc.delta_start + t - 1];
      AD<double> a0 = vars[mpc.a_start + t - 1];

      // at time t + 1
      AD<double> x1 = vars[mpc.x_start + t];
      AD<double> y1 = vars[mpc.y_start + t];
      AD<double> psi1 = vars[mpc.psi_start + t];
      AD<double> v1 = vars[mpc.v_start + t];
      AD<double> cte1 = vars[mpc.cte_start + t];
      AD<double> ePsi1 = vars[mpc.epsi_start + t];

      AD<double> f0 = coeffs[0] + coeffs[1] * x0 +
                      coeffs[2] * CppAD::pow(x0, 2) +
                      coeffs[3] * CppAD::pow(x0, 3);
      AD<double> psiDes0 = CppAD::atan(coeffs[1] + 2 * coeffs[2] * x0 +
                                       3 * coeffs[3] * CppAD::pow(x0, 2));

      // Setup the rest of the model constraints
      fg[1 + mpc.x_start + t] = x1 - (x0 + v0 * CppAD::cos(psi0) * mpc.dt);
      fg[1 + mpc.y_start + t] = y1 - (y0 + v0 * CppAD::sin(psi0) * mpc.dt);
      fg[1 + mpc.psi_start + t] =
          psi1 - (psi0 + v0 * delta0 / MPC::Lf * mpc.dt);
      fg[1 + mpc.v_start + t] = v1 - (v0 + a0 * mpc.dt);
      fg[1 + mpc.cte_start + t] =
          cte1 - ((f0 - y0) + (v0 * CppAD::sin(ePsi0) * mpc.dt));
      fg[1 + mpc.epsi_start + t] =
          ePsi1 - ((psi0 - psiDes0) + v0 * delta0 / MPC::Lf * mpc.dt);
    }
  }
};
