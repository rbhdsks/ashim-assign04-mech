#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr int kTrapGridPoints = 51;

struct ExperimentData {
    double initial_temperature_c = 25.0;
    double ambient_temperature_c = 300.0;
    double melt_temperature_c = 53.5;
    double thermal_diffusivity = 1.0e-4;
    double thermal_conductivity = 400.0;
    double time_seconds = 400.0;
};

struct ScalarIteration {
    int iteration = 0;
    double estimate = 0.0;
    double residual = 0.0;
    double error = std::numeric_limits<double>::quiet_NaN();
    double rate = std::numeric_limits<double>::quiet_NaN();
    std::string phase;
};

struct ScalarMethodResult {
    std::string name;
    std::vector<ScalarIteration> history;
    bool converged = false;
    std::string status;
};

struct SystemIteration {
    int iteration = 0;
    std::array<double, 3> state{};
    std::array<double, 3> residual{};
    double delta_norm = 0.0;
    double residual_norm = 0.0;
};

struct SystemSolveResult {
    std::array<double, 3> initial_guess{};
    std::vector<SystemIteration> history;
    bool converged = false;
    std::string status;
};

double maxAbs(const std::array<double, 3> &values) {
    return std::max({std::fabs(values[0]), std::fabs(values[1]), std::fabs(values[2])});
}

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

double erfcTrapezoidal(double beta) {
    return 1.0 - erfTrapezoidal(beta);
}

class HeatTransferProblem {
  public:
    explicit HeatTransferProblem(ExperimentData data) : data_(data) {
        ratio_ = (data_.melt_temperature_c - data_.initial_temperature_c) /
                 (data_.ambient_temperature_c - data_.initial_temperature_c);
        beta_scale_ = std::sqrt(data_.thermal_diffusivity * data_.time_seconds) /
                      data_.thermal_conductivity;
    }

    double ratio() const { return ratio_; }

    double residual(double h) const {
        const double beta = beta_scale_ * h;
        return 1.0 - std::exp(beta * beta) * erfcTrapezoidal(beta) - ratio_;
    }

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

void finalizeScalarHistory(ScalarMethodResult &result) {
    if (result.history.empty()) {
        return;
    }

    const double reference = result.history.back().estimate;
    result.history.front().error = std::fabs(result.history.front().estimate - reference);
    result.history.front().rate = std::numeric_limits<double>::quiet_NaN();

    for (std::size_t i = 0; i < result.history.size(); ++i) {
        result.history[i].error = std::fabs(result.history[i].estimate - reference);
        if (i == 0) {
            result.history[i].rate = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

        const double previous_error = result.history[i - 1].error;
        const bool same_phase = result.history[i].phase == result.history[i - 1].phase;
        if (!same_phase || previous_error == 0.0) {
            result.history[i].rate = std::numeric_limits<double>::quiet_NaN();
            continue;
        }

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
    ScalarMethodResult result;
    result.name = "bisection";

    const double left_residual = problem.residual(left);
    const double right_residual = problem.residual(right);
    if (left_residual * right_residual > 0.0) {
        result.status = "Initial interval does not bracket the root.";
        return result;
    }

    double a = left;
    double b = right;
    for (int iteration = 1; iteration <= max_iterations; ++iteration) {
        const double midpoint = 0.5 * (a + b);
        const double residual = problem.residual(midpoint);
        result.history.push_back({iteration, midpoint, residual, 0.0, 0.0, "bisection"});

        if (0.5 * (b - a) < tolerance) {
            result.converged = true;
            result.status = "Converged using the bisection stopping criterion.";
            break;
        }

        if (problem.residual(a) * residual < 0.0) {
            b = midpoint;
        } else {
            a = midpoint;
        }
    }

    if (!result.converged && !result.history.empty()) {
        result.status = "Maximum iterations reached before satisfying tolerance.";
    }

    finalizeScalarHistory(result);
    return result;
}

ScalarMethodResult solveNewton(const HeatTransferProblem &problem,
                               double initial_guess,
                               double tolerance,
                               int max_iterations) {
    ScalarMethodResult result;
    result.name = "newton";

    double current = initial_guess;
    for (int iteration = 1; iteration <= max_iterations; ++iteration) {
        const double residual = problem.residual(current);
        const double derivative = problem.derivative(current);
        if (!std::isfinite(residual) || !std::isfinite(derivative) || std::fabs(derivative) < 1.0e-14) {
            result.status = "Newton's method encountered an invalid or near-zero derivative.";
            break;
        }

        const double next = current - residual / derivative;
        const double next_residual = problem.residual(next);
        result.history.push_back({iteration, next, next_residual, 0.0, 0.0, "newton"});

        if (!std::isfinite(next) || !std::isfinite(next_residual)) {
            result.status = "Newton's method produced a non-finite iterate.";
            break;
        }

        if (std::fabs(next - current) < tolerance && std::fabs(next_residual) < tolerance) {
            result.converged = true;
            result.status = "Converged using Newton's method.";
            break;
        }

        current = next;
    }

    if (!result.converged && result.status.empty()) {
        result.status = "Maximum iterations reached. The method did not converge from h = 5000.";
    }

    finalizeScalarHistory(result);
    return result;
}

ScalarMethodResult solveHybrid(const HeatTransferProblem &problem,
                               double left,
                               double right,
                               int bisection_iterations,
                               double tolerance,
                               int max_newton_iterations) {
    ScalarMethodResult result;
    result.name = "hybrid";

    double a = left;
    double b = right;
    for (int iteration = 1; iteration <= bisection_iterations; ++iteration) {
        const double midpoint = 0.5 * (a + b);
        const double residual = problem.residual(midpoint);
        result.history.push_back({iteration, midpoint, residual, 0.0, 0.0, "bisection"});

        if (problem.residual(a) * residual < 0.0) {
            b = midpoint;
        } else {
            a = midpoint;
        }
    }

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
        const int output_iteration = bisection_iterations + iteration;
        result.history.push_back({output_iteration, next, next_residual, 0.0, 0.0, "newton"});

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

    finalizeScalarHistory(result);
    return result;
}

std::array<double, 3> nonlinearSystem(const std::array<double, 3> &state) {
    const double x = state[0];
    const double y = state[1];
    const double z = state[2];
    return {
        x + y + z - 3.0,
        x * x + y * y + z * z - 5.0,
        std::exp(x) + x * y - x * z - 1.0,
    };
}

std::array<std::array<double, 3>, 3> jacobian(const std::array<double, 3> &state) {
    const double x = state[0];
    const double y = state[1];
    const double z = state[2];
    return {{
        {{1.0, 1.0, 1.0}},
        {{2.0 * x, 2.0 * y, 2.0 * z}},
        {{std::exp(x) + y - z, x, -x}},
    }};
}

std::array<double, 3> solveLinear3x3(std::array<std::array<double, 3>, 3> matrix,
                                     std::array<double, 3> rhs) {
    for (int pivot = 0; pivot < 3; ++pivot) {
        int best_row = pivot;
        for (int row = pivot + 1; row < 3; ++row) {
            if (std::fabs(matrix[row][pivot]) > std::fabs(matrix[best_row][pivot])) {
                best_row = row;
            }
        }

        if (std::fabs(matrix[best_row][pivot]) < 1.0e-14) {
            throw std::runtime_error("Jacobian is singular.");
        }

        if (best_row != pivot) {
            std::swap(matrix[pivot], matrix[best_row]);
            std::swap(rhs[pivot], rhs[best_row]);
        }

        for (int row = pivot + 1; row < 3; ++row) {
            const double factor = matrix[row][pivot] / matrix[pivot][pivot];
            for (int col = pivot; col < 3; ++col) {
                matrix[row][col] -= factor * matrix[pivot][col];
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }

    std::array<double, 3> solution{};
    for (int row = 2; row >= 0; --row) {
        double value = rhs[row];
        for (int col = row + 1; col < 3; ++col) {
            value -= matrix[row][col] * solution[col];
        }
        solution[row] = value / matrix[row][row];
    }
    return solution;
}

SystemSolveResult solveSystemWithNewton(const std::array<double, 3> &initial_guess,
                                        double tolerance,
                                        int max_iterations) {
    SystemSolveResult result;
    result.initial_guess = initial_guess;

    std::array<double, 3> current = initial_guess;
    for (int iteration = 1; iteration <= max_iterations; ++iteration) {
        const auto residual = nonlinearSystem(current);
        const auto jac = jacobian(current);

        std::array<double, 3> rhs = {-residual[0], -residual[1], -residual[2]};
        std::array<double, 3> delta{};
        try {
            delta = solveLinear3x3(jac, rhs);
        } catch (const std::exception &error) {
            result.status = error.what();
            break;
        }

        std::array<double, 3> next = {
            current[0] + delta[0],
            current[1] + delta[1],
            current[2] + delta[2],
        };

        const auto next_residual = nonlinearSystem(next);
        const double delta_norm = maxAbs(delta);
        const double residual_norm = maxAbs(next_residual);
        result.history.push_back({iteration, next, next_residual, delta_norm, residual_norm});

        if (delta_norm < tolerance && residual_norm < tolerance) {
            result.converged = true;
            result.status = "Converged using multivariable Newton-Raphson.";
            break;
        }

        current = next;
    }

    if (!result.converged && result.status.empty()) {
        result.status = "Maximum iterations reached before convergence.";
    }

    return result;
}

void writeScalarCsv(const fs::path &path, const ScalarMethodResult &result) {
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

void writeSystemCsv(const fs::path &path, const SystemSolveResult &result) {
    std::ofstream out(path);
    out << "iteration,x,y,z,f1,f2,f3,delta_norm,residual_norm\n";
    out << std::setprecision(15);
    for (const auto &row : result.history) {
        out << row.iteration << ','
            << row.state[0] << ','
            << row.state[1] << ','
            << row.state[2] << ','
            << row.residual[0] << ','
            << row.residual[1] << ','
            << row.residual[2] << ','
            << row.delta_norm << ','
            << row.residual_norm << '\n';
    }
}

void writeSummary(const fs::path &path,
                  const ExperimentData &data,
                  const HeatTransferProblem &problem,
                  const ScalarMethodResult &bisection,
                  const ScalarMethodResult &newton,
                  const ScalarMethodResult &hybrid,
                  const SystemSolveResult &system_one,
                  const SystemSolveResult &system_two) {
    std::ofstream out(path);
    out << std::fixed << std::setprecision(10);
    out << "ME 5107 Assignment 04 Summary\n";
    out << "========================================\n\n";
    out << "Problem 1 data\n";
    out << "Ti   = " << data.initial_temperature_c << " C\n";
    out << "Tinf = " << data.ambient_temperature_c << " C\n";
    out << "Tm   = " << data.melt_temperature_c << " C\n";
    out << "alpha= " << data.thermal_diffusivity << " m^2/s\n";
    out << "k    = " << data.thermal_conductivity << " W/m-K\n";
    out << "t    = " << data.time_seconds << " s\n";
    out << "Target ratio = " << problem.ratio() << "\n\n";

    out << "Problem 1 results\n";
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

    out << "\nProblem 2 results\n";
    out << "Initial guess 1 = [0.1, 1.2, 2.5]\n";
    out << "  status        = " << system_one.status << "\n";
    if (!system_one.history.empty()) {
        const auto &root = system_one.history.back().state;
        out << "  iterations    = " << system_one.history.size() << "\n";
        out << "  root          = [" << root[0] << ", " << root[1] << ", " << root[2] << "]\n";
    }

    out << "Initial guess 2 = [1, 0, 1]\n";
    out << "  status        = " << system_two.status << "\n";
    if (!system_two.history.empty()) {
        const auto &root = system_two.history.back().state;
        out << "  iterations    = " << system_two.history.size() << "\n";
        out << "  root          = [" << root[0] << ", " << root[1] << ", " << root[2] << "]\n";
    }

    if (!system_one.history.empty() && !system_two.history.empty()) {
        const auto &r1 = system_one.history.back().state;
        const auto &r2 = system_two.history.back().state;
        const double distance = std::max({std::fabs(r1[0] - r2[0]),
                                          std::fabs(r1[1] - r2[1]),
                                          std::fabs(r1[2] - r2[2])});
        out << "  same root?    = " << (distance < 1.0e-6 ? "Yes" : "No") << "\n";
    }
}

void printConsoleSummary(const ScalarMethodResult &bisection,
                         const ScalarMethodResult &newton,
                         const ScalarMethodResult &hybrid,
                         const SystemSolveResult &system_one,
                         const SystemSolveResult &system_two) {
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Problem 1\n";
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

    std::cout << "\nProblem 2\n";
    if (!system_one.history.empty()) {
        const auto &root = system_one.history.back().state;
        std::cout << "  Guess [0.1, 1.2, 2.5] -> ["
                  << root[0] << ", " << root[1] << ", " << root[2] << "]"
                  << ", iterations = " << system_one.history.size() << '\n';
    }
    if (!system_two.history.empty()) {
        const auto &root = system_two.history.back().state;
        std::cout << "  Guess [1, 0, 1] -> ["
                  << root[0] << ", " << root[1] << ", " << root[2] << "]"
                  << ", iterations = " << system_two.history.size() << '\n';
    }
}

} // namespace

int main() {
    const ExperimentData data;
    const HeatTransferProblem problem(data);

    fs::create_directories("output");

    const auto bisection = solveBisection(problem, 1.0, 10000.0, 1.0e-5, 100);
    const auto newton = solveNewton(problem, 5000.0, 1.0e-5, 30);
    const auto hybrid = solveHybrid(problem, 1.0, 10000.0, 10, 1.0e-5, 30);

    const auto system_one = solveSystemWithNewton({0.1, 1.2, 2.5}, 1.0e-10, 100);
    const auto system_two = solveSystemWithNewton({1.0, 0.0, 1.0}, 1.0e-10, 100);

    writeScalarCsv("output/part1_bisection.csv", bisection);
    writeScalarCsv("output/part1_newton.csv", newton);
    writeScalarCsv("output/part1_hybrid.csv", hybrid);
    writeSystemCsv("output/part2_guess1.csv", system_one);
    writeSystemCsv("output/part2_guess2.csv", system_two);
    writeSummary("output/summary.txt", data, problem, bisection, newton, hybrid, system_one, system_two);

    printConsoleSummary(bisection, newton, hybrid, system_one, system_two);
    return 0;
}
