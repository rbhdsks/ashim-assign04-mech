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
// Question 1(a)
// ----------------------------------------------------------------------------
// This file solves only Question 1(a) from the assignment.
// The requested method here is the bisection method.
//
// Even though this program is only for 1(a), it still stores error and rate
// information so the plotting/report pipeline remains complete and consistent.
// ============================================================================

// Numerical constant used in the Gauss error function definition.
constexpr double kPi = 3.14159265358979323846;

// The assignment explicitly states that the trapezoidal rule must use 51 points.
constexpr int kTrapGridPoints = 51;

// All physical data from the assignment statement is collected here.
struct ExperimentData {
    double initial_temperature_c = 25.0;
    double ambient_temperature_c = 300.0;
    double melt_temperature_c = 53.5;
    double thermal_diffusivity = 1.0e-4;
    double thermal_conductivity = 400.0;
    double time_seconds = 400.0;
};

// Stores one iteration row for the scalar root-finding history.
struct ScalarIteration {
    int iteration = 0;
    double estimate = 0.0;
    double residual = 0.0;
    double error = std::numeric_limits<double>::quiet_NaN();
    double rate = std::numeric_limits<double>::quiet_NaN();
    std::string phase;
};

// Stores the complete output of the bisection run.
struct ScalarMethodResult {
    std::vector<ScalarIteration> history;
    bool converged = false;
    std::string status;
};

// Evaluates erf(beta) numerically using the trapezoidal rule requested
// in the assignment prompt.
double erfTrapezoidal(double beta) {
    // erf(0) is known exactly.
    if (beta == 0.0) {
        return 0.0;
    }

    // Use the odd symmetry erf(-x) = -erf(x) to handle negative inputs cleanly.
    double sign = 1.0;
    if (beta < 0.0) {
        sign = -1.0;
        beta = -beta;
    }

    // Uniform step size over the interval [0, beta].
    const double dx = beta / static_cast<double>(kTrapGridPoints - 1);

    // Start with the endpoint contributions.
    double sum = std::exp(-0.0) + std::exp(-(beta * beta));

    // Add the interior points, each carrying weight 2 in the trapezoidal rule.
    for (int i = 1; i < kTrapGridPoints - 1; ++i) {
        const double x = static_cast<double>(i) * dx;
        sum += 2.0 * std::exp(-(x * x));
    }

    // Convert the integral into erf(beta).
    const double integral = 0.5 * dx * sum;
    return sign * (2.0 / std::sqrt(kPi)) * integral;
}

// Complementary error function required by the transient heat-transfer formula.
double erfcTrapezoidal(double beta) {
    return 1.0 - erfTrapezoidal(beta);
}

// Encapsulates the scalar nonlinear equation f(h) = 0 for Question 1.
class HeatTransferProblem {
  public:
    explicit HeatTransferProblem(const ExperimentData &data) : data_(data) {
        // Temperature ratio taken directly from the measured temperatures.
        ratio_ = (data_.melt_temperature_c - data_.initial_temperature_c) /
                 (data_.ambient_temperature_c - data_.initial_temperature_c);

        // Because the coating is negligibly thin, z = 0 and the equation reduces
        // to a function of beta = h * sqrt(alpha * t) / k.
        beta_scale_ = std::sqrt(data_.thermal_diffusivity * data_.time_seconds) /
                      data_.thermal_conductivity;
    }

    // Expose the ratio for summaries and explanations.
    double ratio() const { return ratio_; }

    // Residual of the nonlinear equation at a trial heat-transfer coefficient h.
    double residual(double h) const {
        const double beta = beta_scale_ * h;
        return 1.0 - std::exp(beta * beta) * erfcTrapezoidal(beta) - ratio_;
    }

  private:
    ExperimentData data_;
    double ratio_ = 0.0;
    double beta_scale_ = 0.0;
};

// After the method finishes, compute the error at every iteration by treating
// the latest iterate as the "exact" answer, exactly as the assignment states.
void finalizeHistory(ScalarMethodResult &result) {
    if (result.history.empty()) {
        return;
    }

    const double reference = result.history.back().estimate;
    for (std::size_t i = 0; i < result.history.size(); ++i) {
        result.history[i].error = std::fabs(result.history[i].estimate - reference);

        // The first iteration has no previous error, so the rate is undefined.
        if (i == 0) {
            result.history[i].rate = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        const double previous_error = result.history[i - 1].error;
        if (previous_error == 0.0) {
            result.history[i].rate = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        // For bisection, the assignment asks for e_(k+1) / e_k.
        result.history[i].rate = result.history[i].error / previous_error;
    }
}

// Performs the bisection algorithm on the bracket [left, right].
ScalarMethodResult solveBisection(const HeatTransferProblem &problem,
                                  double left,
                                  double right,
                                  double tolerance,
                                  int max_iterations) {
    ScalarMethodResult result;

    // Copy the input interval into working variables.
    double a = left;
    double b = right;
    double fa = problem.residual(a);
    const double fb = problem.residual(b);

    // Bisection requires a sign change across the initial interval.
    if (fa * fb > 0.0) {
        result.status = "Initial interval does not bracket the root.";
        return result;
    }

    for (int iteration = 1; iteration <= max_iterations; ++iteration) {
        // Midpoint of the current interval becomes the new approximation.
        const double midpoint = 0.5 * (a + b);
        const double fm = problem.residual(midpoint);

        // Save the current iteration for plotting and reporting.
        result.history.push_back({iteration, midpoint, fm, 0.0, 0.0, "bisection"});

        // Stop when the active interval is sufficiently small.
        if (0.5 * (b - a) < tolerance) {
            result.converged = true;
            result.status = "Converged using the bisection stopping criterion.";
            break;
        }

        // Keep the half interval that still contains the sign change.
        if (fa * fm < 0.0) {
            b = midpoint;
        } else {
            a = midpoint;
            fa = fm;
        }
    }

    if (!result.converged && !result.history.empty()) {
        result.status = "Maximum iterations reached before satisfying tolerance.";
    }

    finalizeHistory(result);
    return result;
}

// Writes the iteration table to CSV so the plotting script can read it.
void writeCsv(const fs::path &path, const ScalarMethodResult &result) {
    std::ofstream out(path);
    out << "iteration,estimate,residual,error,rate,phase\n";
    out << std::setprecision(15);

    for (const auto &row : result.history) {
        out << row.iteration << ','
            << row.estimate << ','
            << row.residual << ','
            << row.error << ',';

        // Leave the rate cell blank if the value is undefined.
        if (std::isnan(row.rate)) {
            out << ',';
        } else {
            out << row.rate << ',';
        }

        out << row.phase << '\n';
    }
}

// Writes a concise human-readable text summary for Question 1(a).
void writeSummary(const fs::path &path,
                  const ExperimentData &data,
                  const HeatTransferProblem &problem,
                  const ScalarMethodResult &result) {
    std::ofstream out(path);
    out << std::fixed << std::setprecision(10);

    out << "Question 1(a) Summary\n";
    out << "=====================\n\n";
    out << "Method   = Bisection\n";
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

// Prints the same final information to the terminal.
void printSummary(const ScalarMethodResult &result) {
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Question 1(a)\n";
    if (!result.history.empty()) {
        std::cout << "  Method: Bisection\n";
        std::cout << "  h = " << result.history.back().estimate
                  << " W/m^2-K, iterations = " << result.history.size()
                  << ", status = " << result.status << '\n';
    }
}

} // namespace

int main() {
    // Ensure the output directory is present before writing files.
    fs::create_directories("output");

    // Build the scalar nonlinear problem from the assignment data.
    const ExperimentData data;
    const HeatTransferProblem problem(data);

    // Run only the method requested in Question 1(a).
    const auto result = solveBisection(problem, 1.0, 10000.0, 1.0e-5, 100);

    // Save the numeric history and the summary file.
    writeCsv("output/q1a_bisection.csv", result);
    writeSummary("output/q1a_summary.txt", data, problem, result);

    // Print the final answer for quick verification.
    printSummary(result);

    return 0;
}
