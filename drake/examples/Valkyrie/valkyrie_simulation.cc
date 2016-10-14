#include <memory>

#include <lcm/lcm-cpp.hpp>
#include <drake/lcm/drake_lcm.h>

#include "lcmtypes/bot_core/atlas_command_t.hpp"
#include "lcmtypes/bot_core/robot_state_t.hpp"

#include "drake/systems/framework/primitives/pass_through.h"
#include "drake/common/drake_path.h"
#include "drake/common/text_logging.h"
#include "drake/systems/analysis/simulator.h"
#include "drake/systems/framework/diagram_builder.h"
#include "drake/systems/framework/primitives/constant_vector_source.h"
#include "drake/systems/lcm/lcm_publisher_system.h"
#include "drake/systems/lcm/lcm_subscriber_system.h"
#include "drake/systems/lcm/lcmt_drake_signal_translator.h"
#include "drake/systems/plants/parser_urdf.h"
#include "drake/systems/plants/rigid_body_plant/rigid_body_plant.h"
#include "drake/examples/Valkyrie/robot_state_encoder.h"
#include "simple_robot_command_to_desired_effort_converter.h"

// TODO: use syntactic sugar for connections from cars-simulator2 branch

namespace drake {
namespace examples {
namespace valkyrie {

using std::move;
using std::make_unique;
using systems::RigidBodyPlant;
using systems::lcm::LcmPublisherSystem;
using systems::lcm::LcmSubscriberSystem;
using systems::lcm::LcmtDrakeSignalTranslator;
using systems::DiagramContext;
using lcm::DrakeLcm;
using Eigen::VectorXd;
using drake::systems::plants::joints::kRollPitchYaw;
using bot_core::atlas_command_t;
using bot_core::robot_state_t;

int main(int argc, const char **argv) {
  drake::log()->set_level(spdlog::level::trace);
  auto builder = std::make_unique<systems::DiagramBuilder<double>>();

  // Create RigidBodyTree.
  auto tree_ptr = make_unique<RigidBodyTree>();
  drake::parsers::urdf::AddModelInstanceFromUrdfFile(
      drake::GetDrakePath() +
          "/examples/Valkyrie/urdf/urdf/"
          "valkyrie_A_sim_drake_one_neck_dof_wide_ankle_rom.urdf",
      kRollPitchYaw, nullptr /* weld to frame */, tree_ptr.get());

  // Instantiate a RigidBodyPlant from the RigidBodyTree.
  RigidBodyPlant<double> plant(move(tree_ptr));
  const auto& tree = plant.get_rigid_body_tree();

  // LCM communication.
  DrakeLcm lcm;

  // LCM inputs.
  auto atlas_command_subscriber = LcmSubscriberSystem::Make<atlas_command_t>("ATLAS_COMMAND", &lcm);
  std::vector<const RigidBodyActuator*> rigid_body_actuator_ptrs;
  for (const auto& actuator : tree.actuators) {
    rigid_body_actuator_ptrs.push_back(&actuator);
  }
  SimpleRobotCommandToDesiredEffortConverter robot_command_to_desired_effort_converter(rigid_body_actuator_ptrs);

  // LCM outputs.
  RobotStateEncoder robot_state_encoder(plant.get_rigid_body_tree());
  auto robot_state_publisher = LcmPublisherSystem::Make<robot_state_t>("EST_ROBOT_STATE", &lcm);

  // Placeholder for actuator dynamics.
  systems::PassThrough<double> actuators(plant.get_input_size());

  // Placeholder for effort sensors
  systems::PassThrough<double> effort_sensors(plant.get_input_size());

  // TODO: Joint sensors.

  // TODO: Force/torque sensors.

  // Connections
  builder->Connect(*atlas_command_subscriber, robot_command_to_desired_effort_converter);
  builder->Connect(robot_command_to_desired_effort_converter, actuators);
  builder->Connect(actuators, plant);
  builder->Connect(actuators, effort_sensors);
  builder->Connect(plant.get_output_port(0), joint_sensors.get_kinematics_results_port());
  builder->Connect(joint_sensors.get_configuration_port(), robot_state_encoder.get_configuration_port());
  builder->Connect(joint_sensors.get_velocity_port(), robot_state_encoder.get_velocity_port());
  builder->Connect(effort_sensors.get_output_port(0), robot_state_encoder.get_effort_port());
  // TODO: design and connect wrench ports.
  builder->Connect(robot_state_encoder, *robot_state_publisher);

  auto diagram = builder->Build();

  // Create simulator.
  auto simulator = std::make_unique<systems::Simulator<double>>(*diagram);
  simulator->Initialize();
  auto context = simulator->get_mutable_context();

  // Set initial state.
  auto plant_context = diagram->GetMutableSubsystemContext(context, &plant);
  plant.SetZeroConfiguration(plant_context);

  while (true) {
    const double time = context->get_time();
    SPDLOG_TRACE(drake::log(), "Time is now {}", time);
    simulator->StepTo(time + 0.01);
  }
}

}  // valkyrie
}  // examples
}  // drake

int main(int argc, const char* argv[]) {
  return drake::examples::valkyrie::main(argc, argv);
}
