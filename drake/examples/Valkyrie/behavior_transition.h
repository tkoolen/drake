#pragma once

// TODO(tkoolen): this is probably general enough to move to a different
// namespace/folder.
namespace drake {
namespace examples {
namespace valkyrie {

enum class Behavior {
  kFreeze,
  kPositionControl,
  kNormal
};

struct BehaviorTransition {
  Behavior behavior_to_transition_to_;
  double duration_;
};

} // valkyrie
} // examples
} // drake
