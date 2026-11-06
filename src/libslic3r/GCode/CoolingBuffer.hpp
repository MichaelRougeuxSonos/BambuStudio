#ifndef slic3r_CoolingProcess_hpp_
#define slic3r_CoolingProcess_hpp_
#include "../libslic3r.h"
#include "GCodeEditor.hpp"

namespace Slic3r {

class CoolingBuffer
{
public:
    CoolingBuffer(){};

    float calculate_layer_slowdown(std::vector<PerExtruderAdjustments> &per_extruder_adjustments);

    // Slow down an already cooled layer further, so that it takes target_time (bounded by the minimum print speed).
    // Returns the layer time after the stretch, which is never below the time before it.
    float stretch_layer_time(std::vector<PerExtruderAdjustments> &per_extruder_adjustments, float target_time);

private:
    // Old logic: proportional.
    bool m_cooling_logic_proportional = false;
};

}
#endif