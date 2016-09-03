#include <memory>

#include <drake/lcm/drake_lcm.h>

#include "lcmtypes/bot_core/atlas_command_t.hpp"
#include "lcmtypes/bot_core/robot_state_t.hpp"

#include "actuator_effort_to_rigid_body_plant_input_converter.h"
#include "drake/common/drake_path.h"
#include "drake/common/text_logging.h"
#include "drake/examples/Valkyrie/robot_state_encoder.h"
#include "drake/systems/analysis/simulator.h"
#include "drake/systems/framework/diagram_builder.h"
#include "drake/systems/framework/primitives/constant_vector_source.h"
#include "drake/systems/framework/primitives/pass_through.h"
#include "drake/systems/lcm/lcm_publisher_system.h"
#include "drake/systems/lcm/lcm_subscriber_system.h"
#include "drake/systems/lcm/lcmt_drake_signal_translator.h"
#include "drake/systems/plants/parser_urdf.h"
#include "drake/systems/plants/rigid_body_plant/rigid_body_plant.h"
#include "robot_command_to_desired_effort_converter.h"

// TODO: use syntactic sugar for connections from cars-simulator2 branch

namespace drake {
namespace systems {

using std::move;
using std::unique_ptr;
using std::make_unique;
using std::map;

using Eigen::VectorXd;
using bot_core::atlas_command_t;
using bot_core::robot_state_t;

using drake::lcm::DrakeLcm;
using lcm::LcmSubscriberSystem;
using lcm::LcmPublisherSystem;
using plants::joints::kRollPitchYaw;

int main(int argc, const char** argv) {
  drake::log()->set_level(spdlog::level::trace);
  auto builder = std::make_unique<DiagramBuilder<double>>();

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

  // RigidBodyActuators.
  std::vector<const RigidBodyActuator*> actuators;
  for (const auto& actuator : tree.actuators) {
    actuators.push_back(&actuator);
  }
  // Currently, all RigidBodyActuators are assumed to be one-dimensional.
  const int actuator_effort_length = 1;

  // LCM communication.
  DrakeLcm lcm;

  // LCM inputs.
  auto atlas_command_subscriber =
      LcmSubscriberSystem::Make<atlas_command_t>("ROBOT_COMMAND", &lcm);
  RobotCommandToDesiredEffortConverter
      robot_command_to_desired_effort_converter(actuators);

  // Placeholder for actuator dynamics.
  map<const RigidBodyActuator*, unique_ptr<System<double>>> actuator_dynamics;
  for (const auto& actuator : actuators) {
    actuator_dynamics.emplace(std::make_pair(
        actuator, make_unique<PassThrough<double>>(actuator_effort_length)));
  }

  // Conversion from desired efforts to RigidBodyPlant input vector.
  ActuatorEffortToRigidBodyPlantInputConverter<double>
      actuator_effort_to_rigid_body_plant_input_converter(actuators);

  // Placeholder for effort sensors.
  map<const RigidBodyActuator*, unique_ptr<System<double>>> effort_sensors;
  for (const auto& actuator : actuators) {
    effort_sensors.emplace(std::make_pair(
        actuator, make_unique<PassThrough<double>>(actuator_effort_length)));
  }

  // LCM outputs.
  RobotStateEncoder robot_state_encoder(plant.get_rigid_body_tree());
  auto robot_state_publisher =
      LcmPublisherSystem::Make<robot_state_t>("EST_ROBOT_STATE", &lcm);

  // TODO: Force/torque sensors.


  // Connections.
  // LCM message to desired effort conversion.
  builder->Connect(*atlas_command_subscriber,
                   robot_command_to_desired_effort_converter);

  for (const auto& actuator : actuators) {
    // Desired effort inputs to actuator dynamics.
    const auto& desired_effort_output =
        robot_command_to_desired_effort_converter
            .get_desired_effort_output_port(*actuator);
    const auto& desired_effort_input =
        actuator_dynamics.at(actuator)->get_input_port(0);
    builder->Connect(desired_effort_output, desired_effort_input);

    // Efforts to effort sensors.
    const auto& effort_output_port =
        actuator_dynamics.at(actuator)->get_output_port(0);
    const auto& measured_effort_input_port =
        effort_sensors.at(actuator)->get_input_port(0);
    builder->Connect(effort_output_port, measured_effort_input_port);

    // Efforts to rigid body plant input
    builder->Connect(effort_output_port,
                     actuator_effort_to_rigid_body_plant_input_converter
                         .get_effort_input_port(*actuator));
    // TODO: can you connect one output port to multiple input ports?

    // Effort sensors to robot state encoder.
    const auto& measured_effort_output_port = effort_sensors.at(actuator)->get_output_port(0);
    const auto& state_encoder_effort_input_port = robot_state_encoder.get_effort_port(*actuator);
    builder->Connect(measured_effort_output_port, state_encoder_effort_input_port);\
  }

  // Plant input to plant.
  builder->Connect(actuator_effort_to_rigid_body_plant_input_converter, plant);

  // Kinematics results to robot state encoder.
  builder->Connect(plant.get_output_port(1), robot_state_encoder.get_kinematics_results_port()); // FIXME: magic number

  // TODO: design and connect wrench ports.

  // Robot state encoder to robot state publisher.
  builder->Connect(robot_state_encoder, *robot_state_publisher);

  auto diagram = builder->Build();

  // Create simulator.
  auto simulator = std::make_unique<Simulator<double>>(*diagram);
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

}  // systems
}  // drake

int main(int argc, const char* argv[]) {
  return drake::systems::main(argc, argv);
}
