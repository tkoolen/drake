#pragma once

// TODO(tkoolen): namespace
namespace drake {

enum class Behavior {
  kFreeze,
  kPositionControl,
  kNormal
};

struct BehaviorTransition {
  Behavior behavior_to_transition_to_;
  double duration_;
};

} // drake
