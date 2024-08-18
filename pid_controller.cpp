#include <Rcpp.h>

using namespace Rcpp;

#include "pid_controller.hpp"

typedef pid::controller<double> PIDController;

RCPP_MODULE(mod_pid) {
  class_<PIDController>("PIDController")
    .constructor<PIDController::scalar_type, PIDController::scalar_type, PIDController::scalar_type>()
    .property("measure", &PIDController::measure, &PIDController::set_measure)
    .property("control", &PIDController::control, &PIDController::set_control)
    .method("reset", &PIDController::reset)
    .method("monotonic", &PIDController::monotonic)
    .property("out", &PIDController::out)
    .property("proportional", &PIDController::proportional)
    .property("integral", &PIDController::integral)
    .property("derivative", &PIDController::derivative);
}

