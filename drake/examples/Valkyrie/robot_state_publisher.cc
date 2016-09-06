#include "robot_state_publisher.h"

namespace drake {
namespace examples {

RobotStatePublisher::RobotStatePublisher(const std::string &channel,
                                         lcm::LCM *lcm) : channel_(channel),
                                                          lcm_(lcm){
  // TODO: declare input ports
}

RobotStatePublisher::~RobotStatePublisher() { }

std::string RobotStatePublisher::get_name() const {
  return "RobotStatePublisher::" + channel_;
}

void RobotStatePublisher::DoPublish(const systems::ContextBase<double> &context) const {

  const int lcm_message_length = message_.getEncodedSize();
  message_bytes_.resize(static_cast<size_t>(lcm_message_length));
  message_.encode(message_bytes_.data(), 0, lcm_message_length);
  lcm_->publish(channel_, message_bytes_.data(),
                static_cast<unsigned int>(message_bytes_.size()));
}

} // examples
} // drake
