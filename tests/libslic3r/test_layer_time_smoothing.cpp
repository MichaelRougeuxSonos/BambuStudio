#include <catch2/catch.hpp>

#include "libslic3r/GCode/CoolingBuffer.hpp"
#include "libslic3r/GCode/Smoothing.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <random>
#include <string>

using namespace Slic3r;

namespace {

// Reference implementation: t'[i] = max over j of (t[j] - delta * |i - j|), for the layers from first_layer on.
std::vector<float> brute_force_limit(const std::vector<float> &times, float delta, size_t first_layer)
{
    std::vector<float> out = times;
    for (size_t i = first_layer; i < times.size(); ++i)
        for (size_t j = first_layer; j < times.size(); ++j) {
            const float dist = float(i > j ? i - j : j - i);
            out[i]           = std::max(out[i], times[j] - delta * dist);
        }
    return out;
}

} // namespace

TEST_CASE("limit_layer_time_delta ramps short layers up to a long one without a cap", "[LayerTimeSmoothing]")
{
    // The long layer constrains its neighbours by its real time (360 s), not by some cooling threshold.
    std::vector<float> times{100.f, 30.f, 30.f, 360.f, 30.f, 30.f};
    limit_layer_time_delta(times, 100.f);

    const std::vector<float> expected{100.f, 160.f, 260.f, 360.f, 260.f, 160.f};
    REQUIRE(times.size() == expected.size());
    for (size_t i = 0; i < times.size(); ++i)
        CHECK(times[i] == Approx(expected[i]));
}

TEST_CASE("limit_layer_time_delta never lowers a layer and bounds adjacent differences", "[LayerTimeSmoothing]")
{
    std::mt19937                          rng(12345);
    std::uniform_real_distribution<float> time_dist(1.f, 400.f);
    std::uniform_real_distribution<float> delta_dist(1.f, 60.f);

    for (int trial = 0; trial < 200; ++trial) {
        const size_t       n     = 2 + size_t(rng() % 60);
        const float        delta = delta_dist(rng);
        std::vector<float> original(n);
        for (float &t : original)
            t = time_dist(rng);

        std::vector<float> times = original;
        limit_layer_time_delta(times, delta);

        // The first layer is left alone and does not constrain the layers above it.
        CHECK(times[0] == original[0]);
        const std::vector<float> expected = brute_force_limit(original, delta, 1);
        for (size_t i = 0; i < n; ++i) {
            CHECK(times[i] >= original[i]);
            CHECK(times[i] == Approx(expected[i]).margin(1e-2));
            if (i >= 2)
                CHECK(std::abs(times[i] - times[i - 1]) <= delta + 1e-2f);
        }
    }
}

TEST_CASE("limit_layer_time_delta leaves the layers alone when it should", "[LayerTimeSmoothing]")
{
    SECTION("disabled by a zero or negative delta")
    {
        std::vector<float> times{10.f, 300.f, 5.f};
        const auto         original = times;
        limit_layer_time_delta(times, 0.f);
        CHECK(times == original);
        limit_layer_time_delta(times, -3.f);
        CHECK(times == original);
    }
    SECTION("fewer than two layers to smooth")
    {
        std::vector<float> empty;
        limit_layer_time_delta(empty, 10.f);
        CHECK(empty.empty());

        std::vector<float> one{50.f};
        limit_layer_time_delta(one, 10.f);
        CHECK(one == std::vector<float>{50.f});

        std::vector<float> first_and_one{500.f, 5.f};
        limit_layer_time_delta(first_and_one, 10.f);
        CHECK(first_and_one == std::vector<float>{500.f, 5.f});
    }
    SECTION("a long first layer does not drag the layers above it")
    {
        std::vector<float> times{1000.f, 10.f, 10.f, 10.f};
        limit_layer_time_delta(times, 5.f);
        CHECK(times == std::vector<float>{1000.f, 10.f, 10.f, 10.f});
    }
    SECTION("first_layer = 0 includes the first layer")
    {
        std::vector<float> times{1000.f, 10.f, 10.f, 10.f};
        limit_layer_time_delta(times, 5.f, 0);
        CHECK(times == std::vector<float>{1000.f, 995.f, 990.f, 985.f});
    }
}

TEST_CASE("layer_time_shortfalls is what limit_layer_time_delta adds", "[LayerTimeSmoothing]")
{
    const std::vector<float> times{100.f, 30.f, 30.f, 360.f, 30.f, 30.f};
    const std::vector<float> shortfalls = layer_time_shortfalls(times, 100.f);
    const std::vector<float> expected{0.f, 130.f, 230.f, 0.f, 230.f, 130.f};
    REQUIRE(shortfalls.size() == expected.size());
    for (size_t i = 0; i < expected.size(); ++i)
        CHECK(shortfalls[i] == Approx(expected[i]));

    CHECK(layer_time_shortfalls(times, 0.f) == std::vector<float>(times.size(), 0.f));
}

TEST_CASE("layer_gcode_z reads the layer height of the layer", "[LayerTimeSmoothing]")
{
    double z = -1.;
    CHECK_FALSE(layer_gcode_z("G1 X1 Y1 E.1\n", z));
    CHECK(z == -1.);

    CHECK(layer_gcode_z("; CHANGE_LAYER\n; Z_HEIGHT: 43.04\n; LAYER_HEIGHT: 0.12\nG1 X1\n", z));
    CHECK(z == Approx(43.04));

    // The last one wins: a layer can carry the preceding support layers.
    CHECK(layer_gcode_z("; Z_HEIGHT: 1.2\nG1 X1\n; Z_HEIGHT: 1.4\nG1 X2\n", z));
    CHECK(z == Approx(1.4));
}

namespace {

LayerTimePause make_pause(float seconds)
{
    LayerTimePause pause;
    pause.seconds        = seconds;
    pause.z              = 10.;
    pause.max_z          = 325.;
    pause.retract_length = 0.8;
    pause.retract_speed  = 30.;
    return pause;
}

double wait_seconds(const std::string &gcode)
{
    const size_t pos = gcode.find("G4 P");
    return pos == std::string::npos ? -1. : std::atof(gcode.c_str() + pos + 4) / 1000.;
}

} // namespace

TEST_CASE("layer_time_pause_gcode retracts, lifts, waits and comes back", "[LayerTimeSmoothing]")
{
    SECTION("moves are in order and the whole pause takes the requested time")
    {
        const std::string gcode = layer_time_pause_gcode(make_pause(30.f));
        const size_t retract = gcode.find("G1 E-0.8 F1800");
        const size_t lift    = gcode.find("G1 Z15.000 F600");
        const size_t wait    = gcode.find("G4 P");
        const size_t lower   = gcode.find("G1 Z10.000 F600");
        const size_t restore = gcode.find("G1 E0.8 F1800");
        REQUIRE(retract != std::string::npos);
        REQUIRE(lift != std::string::npos);
        REQUIRE(wait != std::string::npos);
        REQUIRE(lower != std::string::npos);
        REQUIRE(restore != std::string::npos);
        CHECK(retract < lift);
        CHECK(lift < wait);
        CHECK(wait < lower);
        CHECK(lower < restore);
        // 2 * 5 mm at 10 mm/s and 2 * 0.8 mm at 30 mm/s are spent moving, the rest is waiting.
        CHECK(wait_seconds(gcode) == Approx(30. - 1. - 2. * 0.8 / 30.).margin(0.002));
    }
    SECTION("a pause that is too short to matter is skipped")
    {
        CHECK(layer_time_pause_gcode(make_pause(0.4f)).empty());
        CHECK(layer_time_pause_gcode(make_pause(0.f)).empty());
    }
    SECTION("a short pause is a plain wait, without retracting or lifting")
    {
        const std::string gcode = layer_time_pause_gcode(make_pause(2.f));
        CHECK(wait_seconds(gcode) == Approx(2.).margin(0.001));
        CHECK(gcode.find("G1 ") == std::string::npos);
        // At 3 s the nozzle moves away again.
        CHECK(layer_time_pause_gcode(make_pause(3.f)).find("G1 Z") != std::string::npos);
    }
    SECTION("the deretraction speed is used when set")
    {
        LayerTimePause pause = make_pause(30.f);
        pause.deretract_speed = 20.;
        CHECK(layer_time_pause_gcode(pause).find("G1 E0.8 F1200") != std::string::npos);
    }
    SECTION("no lift when the nozzle is already at the top of the build volume")
    {
        LayerTimePause pause = make_pause(30.f);
        pause.z = 324.9;
        const std::string gcode = layer_time_pause_gcode(pause);
        CHECK(gcode.find("G1 Z") == std::string::npos);
        CHECK(wait_seconds(gcode) == Approx(30. - 2. * 0.8 / 30.).margin(0.002));
    }
    SECTION("the lift is limited by the printable height")
    {
        LayerTimePause pause = make_pause(30.f);
        pause.z = 323.;
        CHECK(layer_time_pause_gcode(pause).find("G1 Z325.000") != std::string::npos);
    }
    SECTION("firmware retraction uses G10 and G11")
    {
        LayerTimePause pause = make_pause(30.f);
        pause.use_firmware_retraction = true;
        const std::string gcode       = layer_time_pause_gcode(pause);
        CHECK(gcode.find("G10") != std::string::npos);
        CHECK(gcode.find("G11") != std::string::npos);
        CHECK(gcode.find("G1 E") == std::string::npos);
    }
    SECTION("no retraction with absolute E distances, or without a retraction length")
    {
        LayerTimePause pause = make_pause(30.f);
        pause.relative_e     = false;
        CHECK(layer_time_pause_gcode(pause).find("G1 E") == std::string::npos);
        CHECK(layer_time_pause_gcode(pause).find("G4 P") != std::string::npos);

        pause                = make_pause(30.f);
        pause.retract_length = 0.;
        CHECK(layer_time_pause_gcode(pause).find("G1 E") == std::string::npos);
    }
}

namespace {

constexpr float min_speed = 10.f;

CoolingLine make_line(unsigned int type, float length, float feedrate)
{
    CoolingLine line(type, 0, 0);
    line.length   = length;
    line.feedrate = feedrate;
    line.time     = length / feedrate;
    line.time_max = (type & CoolingLine::TYPE_ADJUSTABLE) ? std::max(line.time, length / min_speed) : line.time;
    return line;
}

constexpr unsigned int infill_type    = CoolingLine::TYPE_ADJUSTABLE | CoolingLine::TYPE_G1;
constexpr unsigned int perimeter_type = CoolingLine::TYPE_ADJUSTABLE | CoolingLine::TYPE_EXTERNAL_PERIMETER | CoolingLine::TYPE_G1;
constexpr unsigned int travel_type    = CoolingLine::TYPE_G0;

float total_time(const std::vector<PerExtruderAdjustments> &adjustments)
{
    float total = 0.f;
    for (const PerExtruderAdjustments &adj : adjustments)
        total += adj.elapsed_time_total();
    return total;
}

// 1 s of infill (100 mm at 100 mm/s), 2 s of external perimeter (100 mm at 50 mm/s) and 1 s of travel.
std::vector<PerExtruderAdjustments> make_layer()
{
    std::vector<PerExtruderAdjustments> adjustments(1);
    adjustments[0].lines.push_back(make_line(infill_type, 100.f, 100.f));
    adjustments[0].lines.push_back(make_line(perimeter_type, 100.f, 50.f));
    adjustments[0].lines.push_back(make_line(travel_type, 100.f, 100.f));
    return adjustments;
}

} // namespace

TEST_CASE("CoolingBuffer::stretch_layer_time stretches a cooled layer to a longer time", "[LayerTimeSmoothing]")
{
    CoolingBuffer cooling;

    SECTION("a layer that is long enough already is not touched")
    {
        auto layer = make_layer();
        CHECK(cooling.stretch_layer_time(layer, 3.f) == Approx(4.f));
        for (const CoolingLine &line : layer[0].lines)
            CHECK_FALSE(line.slowdown);
    }
    SECTION("the external perimeters are left alone while other lines can still absorb the stretch")
    {
        auto layer = make_layer();
        CHECK(cooling.stretch_layer_time(layer, 8.f) == Approx(8.f));
        CHECK(total_time(layer) == Approx(8.f));
        CHECK(layer[0].lines[0].feedrate == Approx(20.f));
        CHECK(layer[0].lines[1].feedrate == Approx(50.f));
        CHECK_FALSE(layer[0].lines[1].slowdown);
        CHECK_FALSE(layer[0].lines[2].slowdown);
    }
    SECTION("the external perimeters take the rest once the other lines are at the minimum speed")
    {
        auto layer = make_layer();
        CHECK(cooling.stretch_layer_time(layer, 20.f) == Approx(20.f));
        CHECK(layer[0].lines[0].feedrate == Approx(min_speed));
        CHECK(layer[0].lines[1].time == Approx(9.f));
        CHECK(layer[0].lines[1].feedrate == Approx(100.f / 9.f));
    }
    SECTION("an unreachable target stops at the minimum print speed")
    {
        auto layer = make_layer();
        CHECK(cooling.stretch_layer_time(layer, 100.f) == Approx(21.f));
        CHECK(layer[0].lines[0].feedrate == Approx(min_speed));
        CHECK(layer[0].lines[1].feedrate == Approx(min_speed));
        CHECK(layer[0].lines[2].time == Approx(1.f));
    }
    SECTION("the stretch is shared by all the extruders of the layer")
    {
        std::vector<PerExtruderAdjustments> layer(2);
        layer[0].lines.push_back(make_line(infill_type, 100.f, 100.f));
        layer[1].lines.push_back(make_line(infill_type, 300.f, 100.f));
        CHECK(cooling.stretch_layer_time(layer, 8.f) == Approx(8.f));
        // Both extruders are slowed down by the same factor.
        CHECK(layer[0].lines[0].time == Approx(2.f));
        CHECK(layer[1].lines[0].time == Approx(6.f));
    }
    SECTION("zero length and zero time lines do not produce a bad feedrate")
    {
        auto layer = make_layer();
        layer[0].lines.push_back(make_line(infill_type, 0.f, 100.f));
        cooling.stretch_layer_time(layer, 8.f);
        for (const CoolingLine &line : layer[0].lines)
            CHECK(std::isfinite(line.feedrate));
    }
}
