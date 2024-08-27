/*!
 * \file pid_controller.h
 * \copyright (c) 2024, Roy Ratcliffe, Northumberland, United Kingdom
 * SPDX-License-Identifier: MIT
 *
 * Permission is hereby granted, free of charge,  to any person obtaining a
 * copy  of  this  software  and    associated   documentation  files  (the
 * "Software"), to deal in  the   Software  without  restriction, including
 * without limitation the rights to  use,   copy,  modify,  merge, publish,
 * distribute, sublicense, and/or sell  copies  of   the  Software,  and to
 * permit persons to whom the Software is   furnished  to do so, subject to
 * the following conditions:
 *
 *     The above copyright notice and this permission notice shall be
 *     included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT  WARRANTY OF ANY KIND, EXPRESS
 * OR  IMPLIED,  INCLUDING  BUT  NOT   LIMITED    TO   THE   WARRANTIES  OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR   PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS  OR   COPYRIGHT  HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY,  WHETHER   IN  AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM,  OUT  OF   OR  IN  CONNECTION  WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#pragma once

#ifdef __cplusplus
namespace pid {
template <typename T, typename U = T> class controller {
public:
  typedef T scalar_type;
  typedef U factor_type;

  controller(factor_type kp, factor_type ki, factor_type kd)
      : kp_(kp), ki_(ki), kd_(kd),
        // proportional, integral and derivative terms
        p_(0), i_(0), d_(0),
        // control point, measure point, output of PID controller
        control_(0), measure_(0), out_(0) {}

  void set_measure(scalar_type measure) { measure_ = measure; }
  scalar_type measure() const { return measure_; }

  void set_control(scalar_type control) { control_ = control; }
  scalar_type control() const { return control_; }

  //! \brief Resets the controller.
  //! \details Resetting the PID controller restarts the proportional error as
  //! the difference between the new control and the existing measure. The
  //! integral and derivative reset to zero. Always reset after setting a new
  //! control.
  void reset() {
    p_ = control_ - measure_;
    i_ = 0;
    d_ = 0;
  }

  void monotonic() {
    const scalar_type p = control_ - measure_;
    const scalar_type i = i_ + p;
    const scalar_type d = p - p_;
    out_ = kp_ * p + ki_ * i + kd_ * d;
    p_ = p;
    i_ = i;
    d_ = d;
  }

  scalar_type out() const { return out_; }

  // general accessors

  void set_proportional_factor(factor_type kp) { kp_ = kp; }
  factor_type proportional_factor() const { return kp_; }
  scalar_type proportional() const { return p_; }

  void set_integral_factor(factor_type ki) { ki_ = ki; }
  factor_type integral_factor() const { return ki_; }
  scalar_type integral() const { return i_; }

  void set_derivative_factor(factor_type kd) { kd_ = kd; }
  factor_type derivative_factor() const { return kd_; }
  scalar_type derivative() const { return d_; }

private:
  factor_type kp_, ki_, kd_;
  scalar_type p_, i_, d_;
  scalar_type control_, measure_, out_;
};
} // namespace pid
#endif

struct pid_float_s {
  float kp_, ki_, kd_;
  float p_, i_, d_;
  float control_, measure_, out_;
};

/*!
 * \brief Applies Proportion, Integral and Derivative monotonically.
 * \details First, set up the control point. Feed in the measurement samples
 * while applying the monotonic method until the output matches the control. The
 * name implies that the control hardware runs it periodically at a real-time
 * fixed rate---neither faster nor slower. The application is real-time.
 */
static inline void pid_float_monotonic(struct pid_float_s *pid_float) {
  const float p = pid_float->control_ - pid_float->measure_;
  const float i = pid_float->i_ + p;
  const float d = p - pid_float->p_;
  pid_float->out_ = pid_float->kp_ * p + pid_float->ki_ * i + pid_float->kd_ * d;
  pid_float->p_ = p;
  pid_float->i_ = i;
  pid_float->d_ = d;
}
