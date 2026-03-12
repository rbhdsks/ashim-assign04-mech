#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

// ============================================================================
// Question 1(c)
// ----------------------------------------------------------------------------
// This file solves only Question 1(c) from the assignment.
// The requested method here is the hybrid approach:
//   1. Perform 10 bisection iterations.
//   2. Use the latest bisection value as the initial guess for Newton's method.
//   3. Continue counting the iteration numbers after the switch.
//
// The output stores both phases so the plots can show the switch point and the
// phase-specific convergence behavior.
// ============================================================================

constexpr double kPi = 3.14159265358979323846;
constexpr int kTrapGridPoints = 51;

// Physical parameters given in the assignment.
struct ExperimentData {
    double initial_temperature_c = 25.0;
    double ambient_temperature_c = 300.0;
    double melt_temperature_c = 53.5;
    double thermal_diffusivity = 1.0e-4;
    double thermal_conductivity = 400.0;
    double time_seconds = 400.0;
};

// One history row from either the bisection phase or the Newton phase.
struct ScalarIteration {
    int iteration = 0;
    double estimate = 0.0;
    double residual = 0.0;
    double error = std::numeric_limits<double>::quiet_NaN();
    double rate = std::numeric_limits<double>::quiet_NaN();
    std::string phase;
};

// Holds the full hybrid-method output.
struct ScalarMethodResult {
    std::vector<ScalarIteration> history;
    bool converged = false;
    std::string status;
};

// Computes erf(beta) with the trapezoidal rule using exactly 51 grid points.
double erfTrapezoidal(double beta) {
    if (beta == 0.0) {
        return 0.0;
    }

    double sign = 1.0;
    if (beta < 0.0) {
        sign = -1.0;
        beta = -beta;
    }

    const double dx = beta / static_cast<double>(kTrapGridPoints - 1);
    double sum = std::exp(-0.0) + std::exp(-(beta * beta));
    for (int i = 1; i < kTrapGridPoints - 1; ++i) {
        const double x = static_cast<double>(i) * dx;
        sum += 2.0 * std::exp(-(x * x));
    }

    const double integral = 0.5 * dx * sum;
    return sign * (2.0 / std::sqrt(kPi)) * integral;
}

// Complementary error function required by the analytical temperature formula.
double erfcTrapezoidal(double beta) {
    return 1.0 - erfTrapezoidal(beta);
}

// Represents the scalar nonlinear equation for the unknown convection coefficient.
class HeatTransferProblem {
  public:
    explicit HeatTransferProblem(const ExperimentData &data) : data_(data) {
        ratio_ = (data_.melt_temperature_c - data_.initial_temperature_c) /
                 (data_.ambient_temperature_c - data_.initial_temperature_c);
        beta_scale_ = std::sqrt(data_.thermal_diffusivity * data_.time_seconds) /
                      data_.thermal_conductivity;
    }

    double ratio() const { return ratio_; }

    // Nonlinear residual f(h).
    double residual(double h) const {
        const double beta = beta_scale_ * h;
        return 1.0 - std::exp(beta * beta) * erfcTrapezoidal(beta) - ratio_;
    }

    // Analytical derivative f'(h) needed once the method switches to Newton.
    double derivative(double h) const {
        const double beta = beta_scale_ * h;
        const double exp_term = std::exp(beta * beta);
        const double erfc_term = erfcTrapezoidal(beta);
        return 2.0 * beta_scale_ *
               ((1.0 / std::sqrt(kPi)) - beta * exp_term * erfc_term);
    }

  private:
    ExperimentData data_;
    double ratio_ = 0.0;
    double beta_scale_ = 0.0;
};

// Computes per-iteration error and the correct convergence-rate convention
// for each phase of the hybrid algorithm.
void finalizeHistory(ScalarMethodResult &result) {
    if (result.history.empty()) {
        return;
    }

    const double reference = result.history.back().estimate;
    for (std::size_t i = 0; i < result.history.size(); ++i) {
        result.history[i].error = std::fabs(result.history[i].estimate - reference);

        if (i == 0) {
            result.history[i].rate = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        const bool same_phase = result.history[i].phase == result.history[i - 1].phase;
        const double previous_error = result.history[i - 1].error;
        if (!same_phase || previous_error == 0.0) {
            result.history[i].rate = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        // Apply the convention that matches the current phase.
        if (result.history[i].phase == "bisection") {
            result.history[i].rate = result.history[i].error / previous_error;
        } else {
            result.history[i].rate = result.history[i].error / (previous_error * previous_error);
        }
    }
}

// Performs the Question 1(c) hybrid strategy.
ScalarMethodResult solveHybrid(const HeatTransferProblem &problem,
                               double left,
                               double right,
                               int bisection_iterations,
                               double tolerance,
                               int max_newton_iterations) {
    ScalarMethodResult result;

    // --------------------------
    // Phase 1: Bisection startup
    // --------------------------
    double a = left;
    double b = right;
    double fa = problem.residual(a);
    const double fb = problem.residual(b);

    if (fa * fb > 0.0) {
        result.status = "Initial interval does not bracket the root.";
        return result;
    }

    for (int iteration = 1; iteration <= bisection_iterations; ++iteration) {
        const double midpoint = 0.5 * (a + b);
        const double fm = problem.residual(midpoint);
        result.history.push_back({iteration, midpoint, fm, 0.0, 0.0, "bisection"});

        if (fa * fm < 0.0) {
            b = midpoint;
        } else {
            a = midpoint;
            fa = fm;
        }
    }

    // ---------------------------------
    // Phase 2: Newton refinement phase
    // ---------------------------------
    double current = result.history.back().estimate;
    for (int iteration = 1; iteration <= max_newton_iterations; ++iteration) {
        const double residual = problem.residual(current);
        const double derivative = problem.derivative(current);

        if (!std::isfinite(residual) || !std::isfinite(derivative) || std::fabs(derivative) < 1.0e-14) {
            result.status = "Hybrid method encountered an invalid or near-zero Newton derivative.";
            break;
        }

        const double next = current - residual / derivative;
        const double next_residual = problem.residual(next);

        // Continue the numbering from iteration 11 onward, as required.
        result.history.push_back({bisection_iterations + iteration, next, next_residual, 0.0, 0.0, "newton"});

        if (!std::isfinite(next) || !std::isfinite(next_residual)) {
            result.status = "Hybrid method produced a non-finite Newton iterate.";
            break;
        }

        if (std::fabs(next - current) < tolerance && std::fabs(next_residual) < tolerance) {
            result.converged = true;
            result.status = "Converged after 10 bisection steps followed by Newton iterations.";
            break;
        }

        current = next;
    }

    if (!result.converged && result.status.empty()) {
        result.status = "Maximum Newton iterations reached in the hybrid method.";
    }

    finalizeHistory(result);
    return result;
}

// Writes the hybrid iteration history to CSV.
void writeCsv(const fs::path &path, const ScalarMethodResult &result) {
    std::ofstream out(path);
    out << "iteration,estimate,residual,error,rate,phase\n";
    out << std::setprecision(15);

    for (const auto &row : result.history) {
        out << row.iteration << ','
            << row.estimate << ','
            << row.residual << ','
            << row.error << ',';

        if (std::isnan(row.rate)) {
            out << ',';
        } else {
            out << row.rate << ',';
        }

        out << row.phase << '\n';
    }
}

// Writes the text summary for Question 1(c).
void writeSummary(const fs::path &path,
                  const ExperimentData &data,
                  const HeatTransferProblem &problem,
                  const ScalarMethodResult &result) {
    std::ofstream out(path);
    out << std::fixed << std::setprecision(10);

    out << "Question 1(c) Summary\n";
    out << "=====================\n\n";
    out << "Method   = Hybrid (10 bisection steps + Newton)\n";
    out << "Ti       = " << data.initial_temperature_c << " C\n";
    out << "Tinf     = " << data.ambient_temperature_c << " C\n";
    out << "Tmelt    = " << data.melt_temperature_c << " C\n";
    out << "alpha    = " << data.thermal_diffusivity << " m^2/s\n";
    out << "k        = " << data.thermal_conductivity << " W/m-K\n";
    out << "time     = " << data.time_seconds << " s\n";
    out << "ratio    = " << problem.ratio() << "\n";
    out << "status   = " << result.status << "\n";

    if (!result.history.empty()) {
        out << "iterations = " << result.history.size() << "\n";
        out << "h          = " << result.history.back().estimate << " W/m^2-K\n";
        out << "residual   = " << result.history.back().residual << "\n";
    }
}

// Prints the final hybrid answer to the terminal.
void printSummary(const ScalarMethodResult &result) {
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Question 1(c)\n";
    if (!result.history.empty()) {
        std::cout << "  Method: Hybrid (10 bisection steps + Newton)\n";
        std::cout << "  h = " << result.history.back().estimate
                  << " W/m^2-K, iterations = " << result.history.size()
                  << ", status = " << result.status << '\n';
    }
}

} // namespace

int main() {
    // Ensure the output directory exists.
    fs::create_directories("output");

    // Construct the scalar nonlinear equation from the assignment parameters.
    const ExperimentData data;
    const HeatTransferProblem problem(data);

    // Run only Question 1(c): hybrid bisection + Newton.
    const auto result = solveHybrid(problem, 1.0, 10000.0, 10, 1.0e-5, 30);

    // Save the Question 1(c) outputs.
    writeCsv("output/q1c_hybrid.csv", result);
    writeSummary("output/q1c_summary.txt", data, problem, result);

    // Print the final result to the terminal.
    printSummary(result);

    return 0;
}
