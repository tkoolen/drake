#pragma once

#include <map>
#include <string>

// TODO(tkoolen): this is probably general enough to move to a different
// namespace/folder.
namespace drake {
namespace examples {
namespace valkyrie {

struct JointGains {
  // TODO(tkoolen) perhaps rename these fields.

  double k_q_p_; // corresponds to kp_position in drcsim API
  double k_q_i_; // corresponds to ki_position in drcsim API
  double k_qd_p_; // corresponds to kp_velocity in drcsim API
  double k_f_p_;
  double ff_qd_; // maps to kd_position in drcsim API (there isnt an equivalent gain in the bdi api)
  double ff_qd_d_;
  double ff_f_d_;
  double ff_const_;
};

struct JointCommand {
  double position_;
  double velocity_;
  double effort_;

  JointGains gains_;
};

typedef std::map<std::string, JointCommand> RobotCommand;

} // valkyrie
} // examples
} // drake
