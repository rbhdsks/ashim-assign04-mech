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
// Question 1(b)
// ----------------------------------------------------------------------------
// This file solves only Question 1(b) from the assignment.
// The requested method here is Newton's method with initial guess h = 5000.
//
// The code intentionally keeps all iteration data, including divergence, because
// that behavior is part of the actual numerical result for this starting guess.
// ============================================================================

constexpr double kPi = 3.14159265358979323846;
constexpr int kTrapGridPoints = 51;

// Physical data from the assignment prompt.
struct ExperimentData {
    double initial_temperature_c = 25.0;
    double ambient_temperature_c = 300.0;
    double melt_temperature_c = 53.5;
    double thermal_diffusivity = 1.0e-4;
    double thermal_conductivity = 400.0;
    double time_seconds = 400.0;
};

// One row of scalar iteration history for Newton's method.
struct ScalarIteration {
    int iteration = 0;
    double estimate = 0.0;
    double residual = 0.0;
    double error = std::numeric_limits<double>::quiet_NaN();
    double rate = std::numeric_limits<double>::quiet_NaN();
    std::string phase;
};

// Full output of the Newton solve.
struct ScalarMethodResult {
    std::vector<ScalarIteration> history;
    bool converged = false;
    std::string status;
};

// Numerical evaluation of erf(beta) using the required trapezoidal rule.
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

// Complementary error function appearing in the transient temperature equation.
double erfcTrapezoidal(double beta) {
    return 1.0 - erfTrapezoidal(beta);
}

// Encodes the scalar equation f(h) = 0 and its derivative f'(h).
class HeatTransferProblem {
  public:
    explicit HeatTransferProblem(const ExperimentData &data) : data_(data) {
        ratio_ = (data_.melt_temperature_c - data_.initial_temperature_c) /
                 (data_.ambient_temperature_c - data_.initial_temperature_c);
        beta_scale_ = std::sqrt(data_.thermal_diffusivity * data_.time_seconds) /
                      data_.thermal_conductivity;
    }

    double ratio() const { return ratio_; }

    // Scalar nonlinear residual for the unknown h.
    double residual(double h) const {
        const double beta = beta_scale_ * h;
        return 1.0 - std::exp(beta * beta) * erfcTrapezoidal(beta) - ratio_;
    }

    // Analytical derivative needed by Newton's update formula.
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

// Computes the assignment's error and convergence-rate quantities after the run.
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

        const double previous_error = result.history[i - 1].error;
        if (previous_error == 0.0) {
            result.history[i].rate = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        // For Newton's method, the assignment asks for e_(k+1) / e_k^2.
        result.history[i].rate = result.history[i].error / (previous_error * previous_error);
    }
}

// Performs Newton's method starting from the assignment's required guess h = 5000.
ScalarMethodResult solveNewton(const HeatTransferProblem &problem,
                               double initial_guess,
                               double tolerance,
                               int max_iterations) {
    ScalarMethodResult result;

    double current = initial_guess;
    for (int iteration = 1; iteration <= max_iterations; ++iteration) {
        // Evaluate the residual and derivative at the current guess.
        const double residual = problem.residual(current);
        const double derivative = problem.derivative(current);

        // Newton's method becomes unsafe if the derivative is invalid or tiny.
        if (!std::isfinite(residual) || !std::isfinite(derivative) || std::fabs(derivative) < 1.0e-14) {
            result.status = "Newton's method encountered an invalid or near-zero derivative.";
            break;
        }

        // Standard Newton correction.
        const double next = current - residual / derivative;
        const double next_residual = problem.residual(next);

        result.history.push_back({iteration, next, next_residual, 0.0, 0.0, "newton"});

        // Guard against overflow and NaN propagation.
        if (!std::isfinite(next) || !std::isfinite(next_residual)) {
            result.status = "Newton's method produced a non-finite iterate.";
            break;
        }

        // Use both step size and residual as the stopping condition.
        if (std::fabs(next - current) < tolerance && std::fabs(next_residual) < tolerance) {
            result.converged = true;
            result.status = "Converged using Newton's method.";
            break;
        }

        // Continue iterating from the updated value.
        current = next;
    }

    if (!result.converged && result.status.empty()) {
        result.status = "Maximum iterations reached. The method did not converge from h = 5000.";
    }

    finalizeHistory(result);
    return result;
}

// Writes the Newton iteration history to CSV for plotting.
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

// Writes a readable summary file specifically for Question 1(b).
void writeSummary(const fs::path &path,
                  const ExperimentData &data,
                  const HeatTransferProblem &problem,
                  const ScalarMethodResult &result) {
    std::ofstream out(path);
    out << std::fixed << std::setprecision(10);

    out << "Question 1(b) Summary\n";
    out << "=====================\n\n";
    out << "Method   = Newton\n";
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

// Prints the final terminal summary for the Newton run.
void printSummary(const ScalarMethodResult &result) {
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Question 1(b)\n";
    if (!result.history.empty()) {
        std::cout << "  Method: Newton\n";
        std::cout << "  h = " << result.history.back().estimate
                  << " W/m^2-K, iterations = " << result.history.size()
                  << ", status = " << result.status << '\n';
    }
}

} // namespace

int main() {
    // Ensure output/ exists before writing CSV and summary files.
    fs::create_directories("output");

    // Build the heat-transfer problem from the assignment data.
    const ExperimentData data;
    const HeatTransferProblem problem(data);

    // Run only Question 1(b): Newton's method.
    const auto result = solveNewton(problem, 5000.0, 1.0e-5, 30);

    // Save the results in separate Question 1(b) files.
    writeCsv("output/q1b_newton.csv", result);
    writeSummary("output/q1b_summary.txt", data, problem, result);

    // Print the final information to the terminal.
    printSummary(result);

    return 0;
}
