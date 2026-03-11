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

// Numerical constant used in the Gauss error function definition.
constexpr double kPi = 3.14159265358979323846;

// The assignment explicitly asks us to use 51 grid points in the trapezoidal rule.
constexpr int kTrapGridPoints = 51;

// All physical input data for Question 1 is grouped here so the solver code stays clean.
struct ExperimentData {
    double initial_temperature_c = 25.0;
    double ambient_temperature_c = 300.0;
    double melt_temperature_c = 53.5;
    double thermal_diffusivity = 1.0e-4;
    double thermal_conductivity = 400.0;
    double time_seconds = 400.0;
};

// One row of scalar iteration history for bisection, Newton, or the hybrid method.
struct ScalarIteration {
    int iteration = 0;
    double estimate = 0.0;
    double residual = 0.0;
    double error = std::numeric_limits<double>::quiet_NaN();
    double rate = std::numeric_limits<double>::quiet_NaN();
    std::string phase;
};

// Stores the full result of one scalar root-finding method.
struct ScalarMethodResult {
    std::string name;
    std::vector<ScalarIteration> history;
    bool converged = false;
    std::string status;
};

// Numerically evaluates erf(beta) using the trapezoidal rule requested in the assignment.
double erfTrapezoidal(double beta) {
    // erf(0) is exactly zero, so we can return immediately.
    if (beta == 0.0) {
        return 0.0;
    }

    // Handle negative beta values using erf(-x) = -erf(x).
    double sign = 1.0;
    if (beta < 0.0) {
        sign = -1.0;
        beta = -beta;
    }

    // Uniform step size over [0, beta].
    const double dx = beta / static_cast<double>(kTrapGridPoints - 1);

    // Trapezoidal rule endpoint terms.
    double sum = std::exp(-0.0) + std::exp(-(beta * beta));

    // Add the interior terms with weight 2.
    for (int i = 1; i < kTrapGridPoints - 1; ++i) {
        const double x = static_cast<double>(i) * dx;
        sum += 2.0 * std::exp(-(x * x));
    }

    // Convert the numerical integral into erf(beta).
    const double integral = 0.5 * dx * sum;
    return sign * (2.0 / std::sqrt(kPi)) * integral;
}

// erfc(beta) is defined directly from erf(beta).
double erfcTrapezoidal(double beta) {
    return 1.0 - erfTrapezoidal(beta);
}

// Encapsulates the reduced surface-temperature equation for the unknown coefficient h.
class HeatTransferProblem {
  public:
    explicit HeatTransferProblem(ExperimentData data) : data_(data) {
        // This is the left-hand-side temperature ratio from the assignment data.
        ratio_ = (data_.melt_temperature_c - data_.initial_temperature_c) /
                 (data_.ambient_temperature_c - data_.initial_temperature_c);

        // Because z = 0 at the surface, beta becomes h * sqrt(alpha * t) / k.
        beta_scale_ = std::sqrt(data_.thermal_diffusivity * data_.time_seconds) /
                      data_.thermal_conductivity;
    }

    // Exposes the target ratio so we can print it in the report/summary.
    double ratio() const { return ratio_; }

    // f(h) = 0 is the scalar nonlinear equation we need to solve.
    double residual(double h) const {
        const double beta = beta_scale_ * h;
        return 1.0 - std::exp(beta * beta) * erfcTrapezoidal(beta) - ratio_;
    }

    // Analytical derivative of f(h), used by Newton's method.
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

// Post-processes iteration history to compute the per-iteration error and convergence rate.
void finalizeHistory(ScalarMethodResult &result) {
    if (result.history.empty()) {
        return;
    }

    // The assignment asks us to treat the latest value as the "exact" solution.
    const double reference = result.history.back().estimate;
    for (std::size_t i = 0; i < result.history.size(); ++i) {
        result.history[i].error = std::fabs(result.history[i].estimate - reference);

        // There is no previous iteration for k = 0, so the rate is undefined here.
        if (i == 0) {
            result.history[i].rate = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        // Only compare rates within the same method phase.
        const bool same_phase = result.history[i].phase == result.history[i - 1].phase;
        const double previous_error = result.history[i - 1].error;
        if (!same_phase || previous_error == 0.0) {
            result.history[i].rate = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        // The assignment gives different convergence-rate formulas for bisection and Newton.
        if (result.history[i].phase == "bisection") {
            result.history[i].rate = result.history[i].error / previous_error;
        } else {
            result.history[i].rate = result.history[i].error / (previous_error * previous_error);
        }
    }
}

ScalarMethodResult solveBisection(const HeatTransferProblem &problem,
                                  double left,
                                  double right,
                                  double tolerance,
                                  int max_iterations) {
    // This function implements Question 1(a).
    ScalarMethodResult result;
    result.name = "bisection";

    // Start from the interval given in the assignment.
    double a = left;
    double b = right;
    double fa = problem.residual(a);
    double fb = problem.residual(b);

    // Bisection only works if the initial interval actually brackets a root.
    if (fa * fb > 0.0) {
        result.status = "Initial interval does not bracket the root.";
        return result;
    }

    for (int iteration = 1; iteration <= max_iterations; ++iteration) {
        // The current approximation is the midpoint of the active interval.
        const double midpoint = 0.5 * (a + b);
        const double fm = problem.residual(midpoint);
        result.history.push_back({iteration, midpoint, fm, 0.0, 0.0, "bisection"});

        // Stop when the interval width is small enough.
        if (0.5 * (b - a) < tolerance) {
            result.converged = true;
            result.status = "Converged using the bisection stopping criterion.";
            break;
        }

        // Keep the half-interval that still contains the sign change.
        if (fa * fm < 0.0) {
            b = midpoint;
            fb = fm;
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

ScalarMethodResult solveNewton(const HeatTransferProblem &problem,
                               double initial_guess,
                               double tolerance,
                               int max_iterations) {
    // This function implements Question 1(b).
    ScalarMethodResult result;
    result.name = "newton";

    // Use the assignment's required initial guess.
    double current = initial_guess;
    for (int iteration = 1; iteration <= max_iterations; ++iteration) {
        // Evaluate f(h) and f'(h) at the current iterate.
        const double residual = problem.residual(current);
        const double derivative = problem.derivative(current);

        // Abort if Newton's formula would be numerically unsafe.
        if (!std::isfinite(residual) || !std::isfinite(derivative) || std::fabs(derivative) < 1.0e-14) {
            result.status = "Newton's method encountered an invalid or near-zero derivative.";
            break;
        }

        // Standard Newton update: h_(k+1) = h_k - f(h_k) / f'(h_k).
        const double next = current - residual / derivative;
        const double next_residual = problem.residual(next);
        result.history.push_back({iteration, next, next_residual, 0.0, 0.0, "newton"});

        // Guard against overflow or NaN values.
        if (!std::isfinite(next) || !std::isfinite(next_residual)) {
            result.status = "Newton's method produced a non-finite iterate.";
            break;
        }

        // Stop only when both the step size and the residual are small.
        if (std::fabs(next - current) < tolerance && std::fabs(next_residual) < tolerance) {
            result.converged = true;
            result.status = "Converged using Newton's method.";
            break;
        }

        // Continue Newton iteration from the new value.
        current = next;
    }

    if (!result.converged && result.status.empty()) {
        result.status = "Maximum iterations reached. The method did not converge from h = 5000.";
    }

    finalizeHistory(result);
    return result;
}

ScalarMethodResult solveHybrid(const HeatTransferProblem &problem,
                               double left,
                               double right,
                               int bisection_iterations,
                               double tolerance,
                               int max_newton_iterations) {
    // This function implements Question 1(c): 10 bisection steps, then Newton.
    ScalarMethodResult result;
    result.name = "hybrid";

    // Phase 1: safely move near the root using bisection.
    double a = left;
    double b = right;
    double fa = problem.residual(a);
    for (int iteration = 1; iteration <= bisection_iterations; ++iteration) {
        const double midpoint = 0.5 * (a + b);
        const double fm = problem.residual(midpoint);
        result.history.push_back({iteration, midpoint, fm, 0.0, 0.0, "bisection"});

        // Keep the subinterval that contains the root.
        if (fa * fm < 0.0) {
            b = midpoint;
        } else {
            a = midpoint;
            fa = fm;
        }
    }

    // Phase 2: switch to Newton starting from the 10th bisection estimate.
    double current = result.history.back().estimate;
    for (int iteration = 1; iteration <= max_newton_iterations; ++iteration) {
        const double residual = problem.residual(current);
        const double derivative = problem.derivative(current);

        // Stop if the Newton correction cannot be trusted.
        if (!std::isfinite(residual) || !std::isfinite(derivative) || std::fabs(derivative) < 1.0e-14) {
            result.status = "Hybrid method encountered an invalid or near-zero Newton derivative.";
            break;
        }

        // Continue the iteration count from where bisection left off.
        const double next = current - residual / derivative;
        const double next_residual = problem.residual(next);
        result.history.push_back({bisection_iterations + iteration, next, next_residual, 0.0, 0.0, "newton"});

        // Same stopping test as the standalone Newton method.
        if (std::fabs(next - current) < tolerance && std::fabs(next_residual) < tolerance) {
            result.converged = true;
            result.status = "Converged after 10 bisection steps followed by Newton iterations.";
            break;
        }

        // Update the working value and keep iterating.
        current = next;
    }

    if (!result.converged && result.status.empty()) {
        result.status = "Maximum Newton iterations reached in the hybrid method.";
    }

    finalizeHistory(result);
    return result;
}

// Writes one CSV file per method so the plotting script can build the required graphs.
void writeCsv(const fs::path &path, const ScalarMethodResult &result) {
    std::ofstream out(path);
    out << "iteration,estimate,residual,error,rate,phase\n";
    out << std::setprecision(15);

    // Each row stores the iteration count, current estimate, residual, error, rate, and phase.
    for (const auto &row : result.history) {
        out << row.iteration << ','
            << row.estimate << ','
            << row.residual << ','
            << row.error << ',';

        // The first row or phase transition can have an undefined rate.
        if (std::isnan(row.rate)) {
            out << ',';
        } else {
            out << row.rate << ',';
        }
        out << row.phase << '\n';
    }
}

void writeSummary(const fs::path &path,
                  const ExperimentData &data,
                  const HeatTransferProblem &problem,
                  const ScalarMethodResult &bisection,
                  const ScalarMethodResult &newton,
                  const ScalarMethodResult &hybrid) {
    // Write a short human-readable summary file for Question 1.
    std::ofstream out(path);
    out << std::fixed << std::setprecision(10);
    out << "Question 1 Summary\n";
    out << "==================\n\n";
    out << "Ti   = " << data.initial_temperature_c << " C\n";
    out << "Tinf = " << data.ambient_temperature_c << " C\n";
    out << "Tm   = " << data.melt_temperature_c << " C\n";
    out << "alpha= " << data.thermal_diffusivity << " m^2/s\n";
    out << "k    = " << data.thermal_conductivity << " W/m-K\n";
    out << "t    = " << data.time_seconds << " s\n";
    out << "Target ratio = " << problem.ratio() << "\n\n";

    // Print the final result for each scalar method.
    out << "Bisection: " << bisection.status << "\n";
    if (!bisection.history.empty()) {
        out << "  iterations = " << bisection.history.size() << "\n";
        out << "  h          = " << bisection.history.back().estimate << " W/m^2-K\n";
        out << "  residual   = " << bisection.history.back().residual << "\n";
    }

    out << "Newton: " << newton.status << "\n";
    if (!newton.history.empty()) {
        out << "  iterations = " << newton.history.size() << "\n";
        out << "  h          = " << newton.history.back().estimate << " W/m^2-K\n";
        out << "  residual   = " << newton.history.back().residual << "\n";
    }

    out << "Hybrid: " << hybrid.status << "\n";
    if (!hybrid.history.empty()) {
        out << "  iterations = " << hybrid.history.size() << "\n";
        out << "  h          = " << hybrid.history.back().estimate << " W/m^2-K\n";
        out << "  residual   = " << hybrid.history.back().residual << "\n";
    }
}

void printSummary(const ScalarMethodResult &bisection,
                  const ScalarMethodResult &newton,
                  const ScalarMethodResult &hybrid) {
    // Mirror the most important results to the terminal for quick checking.
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Question 1\n";
    if (!bisection.history.empty()) {
        std::cout << "  Bisection: h = " << bisection.history.back().estimate
                  << " W/m^2-K, iterations = " << bisection.history.size()
                  << ", status = " << bisection.status << '\n';
    }
    if (!newton.history.empty()) {
        std::cout << "  Newton:    h = " << newton.history.back().estimate
                  << " W/m^2-K, iterations = " << newton.history.size()
                  << ", status = " << newton.status << '\n';
    }
    if (!hybrid.history.empty()) {
        std::cout << "  Hybrid:    h = " << hybrid.history.back().estimate
                  << " W/m^2-K, iterations = " << hybrid.history.size()
                  << ", status = " << hybrid.status << '\n';
    }
}

} // namespace

int main() {
    // Make sure the output directory exists before writing CSV and summary files.
    fs::create_directories("output");

    // Load the assignment data and build the scalar residual function.
    const ExperimentData data;
    const HeatTransferProblem problem(data);

    // Run all three methods requested in Question 1.
    const auto bisection = solveBisection(problem, 1.0, 10000.0, 1.0e-5, 100);
    const auto newton = solveNewton(problem, 5000.0, 1.0e-5, 30);
    const auto hybrid = solveHybrid(problem, 1.0, 10000.0, 10, 1.0e-5, 30);

    // Save the iteration histories for plotting and grading.
    writeCsv("output/part1_bisection.csv", bisection);
    writeCsv("output/part1_newton.csv", newton);
    writeCsv("output/part1_hybrid.csv", hybrid);

    // Save and print the final human-readable summary.
    writeSummary("output/question1_summary.txt", data, problem, bisection, newton, hybrid);
    printSummary(bisection, newton, hybrid);

    // Normal successful program exit.
    return 0;
}
