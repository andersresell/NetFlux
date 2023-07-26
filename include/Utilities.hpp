#pragma once
#include "Includes.hpp"
// #include "containers/StaticContainer.hpp"

constexpr ShortIndex N_EQS_EULER{N_DIM + 2};

// --------------------------------------------------------------------
// Some enums that specify solver behaviour and corresponding string -> enum maps for some
// --------------------------------------------------------------------
enum class MainSolverType
{
    Euler
};

const map<string, MainSolverType> main_solver_from_string{
    {"Euler", MainSolverType::Euler}};

enum class SolverType
{
    Euler
}; // The main solver can (but need not) be comprised of multiple sub solvers

const map<SolverType, string> string_from_solver_type{
    {SolverType::Euler, "Euler"}};

enum class TimeIntegrationType
{
    Explicit,
    Implicit
};

enum class TimeScheme
{
    ExplicitEuler,
    TVD_RK3
};

const map<string, TimeScheme> time_scheme_from_string{
    {"ExplicitEuler", TimeScheme::ExplicitEuler},
    {"TVD_RK3", TimeScheme::TVD_RK3}};

enum class SpatialOrder
{
    First,
    Second
};

const map<string, SpatialOrder> spatial_order_from_string{
    {"First", SpatialOrder::First},
    {"Second", SpatialOrder::Second}};

enum class GradientScheme
{
    GreenGauss
};

const map<string, GradientScheme> gradient_scheme_from_string{
    {"GreenGauss", GradientScheme::GreenGauss}};

enum class InviscidFluxScheme
{
    Rusanov,
    HLLC
};

const map<string, InviscidFluxScheme> inviscid_flux_scheme_from_string{
    {"Rusanov", InviscidFluxScheme::Rusanov},
    {"HLLC", InviscidFluxScheme::HLLC}};

enum class Limiter
{
    NONE,
    Barth
};

const map<string, Limiter> limiter_from_string{
    {"NONE", Limiter::NONE},
    {"Barth", Limiter::Barth}};

enum class BoundaryType
{
    NoSlipWall,
    SlipWall,
    FarField
};

const map<string, BoundaryType> boundary_type_from_string{
    {"NoSlipWall", BoundaryType::NoSlipWall},
    {"SlipWall", BoundaryType::SlipWall},
    {"FarField", BoundaryType::FarField}};

enum class InitialConditionOption
{
    Freestream
};

const map<string, InitialConditionOption> initial_condition_option_from_string{
    {"Freestream", InitialConditionOption::Freestream}};

namespace standard_air
{
    constexpr Scalar gas_constant{287.058};
    constexpr Scalar gamma{1.4};
    constexpr Scalar density{1.225};
    constexpr Scalar pressure{101325.0};
    constexpr Scalar temperature = pressure / (density * gas_constant);
}

class Timer
{
    using Clock = std::chrono::high_resolution_clock;
    using Time = std::chrono::_V2::system_clock::time_point;
    using Milliseconds = std::chrono::milliseconds;

    Time start_time; // Used for measuring time

public:
    void start_counter() { start_time = Clock::now(); }
    string get_elapsed_time() const;
};