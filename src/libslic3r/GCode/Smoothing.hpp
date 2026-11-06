#ifndef slic3r_Smoothing_hpp_
#define slic3r_Smoothing_hpp_
#include "../libslic3r.h"
#include <libslic3r/GCode/GCodeEditor.hpp>
#include <math.h>
#include <string>
#include <vector>

namespace Slic3r {

static const int   guassian_window_size               = 11;
static const int   guassian_r                         = 2;
static const int   guassian_stop_threshold            = 5;
static const float guassian_layer_time_stop_threshold = 3.0;
static const int   max_steps_count             = 1000;

struct CoolingNode
{
    // extruder pos, line pos;
    std::vector<std::pair<int, int>> outwall_line;
    float                            max_feedrate    = 0;
    float                            filter_feedrate = 0;
    double                           rate            = 1;
};

struct OutwallCollection
{
    int                        object_id;
    std::map<int, CoolingNode> cooling_nodes;
};

class SmoothCalculator
{

public:
    std::vector<std::map<int, std::pair<int, int>>> objects_node_range;
    std::vector<std::vector<OutwallCollection>> layers_wall_collection;
    std::vector<float>                          layers_cooling_time;

    SmoothCalculator(const int objects_size, const double gap_limit) : layer_time_smoothing_threshold(gap_limit)
    {
        guassian_filter_generator();
        objects_node_range.resize(objects_size);
    }

    SmoothCalculator(const int objects_size)
    {
        guassian_filter_generator();
        objects_node_range.resize(objects_size);
    }

    void append_data(const std::vector<OutwallCollection> &wall_collection, float cooling_time)
    {
        layers_wall_collection.push_back(wall_collection);
        layers_cooling_time.push_back(cooling_time);
    }

    void append_data(const std::vector<OutwallCollection> &wall_collection)
    {
        layers_wall_collection.push_back(wall_collection);
    }

    void build_node(std::vector<OutwallCollection> &wall_collection, const std::vector<int> &object_label, const std::vector<PerExtruderAdjustments> &per_extruder_adjustments);

    float recaculate_layer_time(int layer_id, std::vector<PerExtruderAdjustments> &extruder_adjustments);

    void smooth_layer_speed();

private:
    // guassian filter
    double guassian_function(double x, double r) {
        return exp(-x * x / (2 * r * r)) / (r * sqrt(2 * PI));
    }

    void guassian_filter_generator() {
        double r = guassian_r;
        int    half_win_size = guassian_window_size / 2;
        for (int start = -half_win_size; start <= half_win_size; ++start) {
            double y = guassian_function(start, r);
            filter_sum += y;
            guassian_filter.push_back(y);
        }
    }

    void init_object_node_range();

    // filter the data
    void layer_speed_filter(const int object_id, const int node_id);

    bool speed_filter_continue(const int object_id, const int node_id);

    // filter the data
    void filter_layer_time();

    bool layer_time_filter_continue();
    void smooth_layer_time();

    std::vector<double>                             guassian_filter;
    double                                          filter_sum                     = .0f;
    float                                           layer_time_smoothing_threshold = 30.0f;
};

// Raise layer times so that no two adjacent layers differ by more than max_delta seconds.
// Times are only ever increased, and a layer constrains its neighbours by its real time (no cap).
// Layers before first_layer are neither raised nor used as a constraint.
void limit_layer_time_delta(std::vector<float> &times, float max_delta, size_t first_layer = 1);

// How much time each layer still lacks to satisfy the delta limit (zero for layers that are fine).
std::vector<float> layer_time_shortfalls(const std::vector<float> &times, float max_delta, size_t first_layer = 1);

// Parameters of a pause that is appended to a layer to make up for layer time that slowing down could not add.
struct LayerTimePause
{
    // Time the pause should add to the layer, including the moves it needs.
    float  seconds = 0.f;
    // Z of the layer that has just been printed, and the highest Z the printer can reach.
    double z                       = 0.;
    double max_z                   = 0.;
    double retract_length          = 0.;
    // mm/s; the deretraction speed falls back to the retraction speed when it is zero.
    double retract_speed           = 0.;
    double deretract_speed         = 0.;
    bool   use_firmware_retraction = false;
    // The filament is only retracted with relative E distances, or with firmware retraction.
    bool   relative_e              = true;
};

// Retract, lift the nozzle away from the print, wait and come back. Empty if the pause is too short to matter.
std::string layer_time_pause_gcode(const LayerTimePause &pause);

// Z of the last "; Z_HEIGHT: " comment in the G-code of a layer.
bool layer_gcode_z(const std::string &gcode, double &z);

} // namespace Slic3r

#endif