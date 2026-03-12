#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

// ============================================================================
// Question 2
// ----------------------------------------------------------------------------
// This file solves Question 2 from the assignment.
// Unlike Question 1, Question 2 is not split into subparts in the prompt, so
// it remains as a single standalone program.
// ============================================================================

// One row of iteration history for the 3-variable Newton-Raphson solver.
struct SystemIteration {
    int iteration = 0;
    std::array<double, 3> state{};
    std::array<double, 3> residual{};
    double delta_norm = 0.0;
    double residual_norm = 0.0;
};

// Stores the complete result for one initial guess.
struct SystemSolveResult {
    std::array<double, 3> initial_guess{};
    std::vector<SystemIteration> history;
    bool converged = false;
    std::string status;
};

// Uses the infinity norm because it is simple and works well for stopping tests.
double maxAbs(const std::array<double, 3> &values) {
    return std::max({std::fabs(values[0]), std::fabs(values[1]), std::fabs(values[2])});
}

// Returns the vector F(x, y, z) for the nonlinear system in Question 2.
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

// Returns the Jacobian matrix dF/d[x y z] needed by Newton-Raphson.
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

// Solves a 3x3 linear system using Gaussian elimination with partial pivoting.
std::array<double, 3> solveLinear3x3(std::array<std::array<double, 3>, 3> matrix,
                                     std::array<double, 3> rhs) {
    for (int pivot = 0; pivot < 3; ++pivot) {
        // Choose the row with the largest pivot magnitude for numerical stability.
        int best_row = pivot;
        for (int row = pivot + 1; row < 3; ++row) {
            if (std::fabs(matrix[row][pivot]) > std::fabs(matrix[best_row][pivot])) {
                best_row = row;
            }
        }

        // If the pivot is effectively zero, the Jacobian is singular.
        if (std::fabs(matrix[best_row][pivot]) < 1.0e-14) {
            throw std::runtime_error("Jacobian is singular.");
        }

        // Move the best pivot row into place.
        if (best_row != pivot) {
            std::swap(matrix[pivot], matrix[best_row]);
            std::swap(rhs[pivot], rhs[best_row]);
        }

        // Eliminate the entries below the pivot.
        for (int row = pivot + 1; row < 3; ++row) {
            const double factor = matrix[row][pivot] / matrix[pivot][pivot];
            for (int col = pivot; col < 3; ++col) {
                matrix[row][col] -= factor * matrix[pivot][col];
            }
            rhs[row] -= factor * rhs[pivot];
        }
    }

    // Recover the solution by back-substitution.
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
    // This function performs the multivariable Newton-Raphson iteration.
    SystemSolveResult result;
    result.initial_guess = initial_guess;

    // Start from the guess provided by the assignment.
    std::array<double, 3> current = initial_guess;
    for (int iteration = 1; iteration <= max_iterations; ++iteration) {
        // Evaluate the nonlinear system and its Jacobian at the current point.
        const auto residual = nonlinearSystem(current);
        const auto jac = jacobian(current);

        // Newton's method solves J * delta = -F to get the correction delta.
        std::array<double, 3> rhs = {-residual[0], -residual[1], -residual[2]};

        std::array<double, 3> delta{};
        try {
            delta = solveLinear3x3(jac, rhs);
        } catch (const std::exception &error) {
            result.status = error.what();
            break;
        }

        // Update the state vector with the Newton correction.
        std::array<double, 3> next = {
            current[0] + delta[0],
            current[1] + delta[1],
            current[2] + delta[2],
        };

        // Recompute the residual at the updated point for convergence monitoring.
        const auto next_residual = nonlinearSystem(next);
        const double delta_norm = maxAbs(delta);
        const double residual_norm = maxAbs(next_residual);
        result.history.push_back({iteration, next, next_residual, delta_norm, residual_norm});

        // Stop only when both the correction and the residual are small.
        if (delta_norm < tolerance && residual_norm < tolerance) {
            result.converged = true;
            result.status = "Converged using multivariable Newton-Raphson.";
            break;
        }

        // Continue iterating from the new point.
        current = next;
    }

    if (!result.converged && result.status.empty()) {
        result.status = "Maximum iterations reached before convergence.";
    }

    return result;
}

// Writes the full iteration history for one starting guess to CSV.
void writeCsv(const fs::path &path, const SystemSolveResult &result) {
    std::ofstream out(path);
    out << "iteration,x,y,z,f1,f2,f3,delta_norm,residual_norm\n";
    out << std::setprecision(15);

    // Each row records the current state, residual, and convergence norms.
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
                  const SystemSolveResult &guess_one,
                  const SystemSolveResult &guess_two) {
    // Write a short report for Question 2 so the root comparison is easy to read.
    std::ofstream out(path);
    out << std::fixed << std::setprecision(10);
    out << "Question 2 Summary\n";
    out << "==================\n\n";

    out << "Initial guess 1 = [0.1, 1.2, 2.5]\n";
    out << "  status        = " << guess_one.status << "\n";
    if (!guess_one.history.empty()) {
        const auto &root = guess_one.history.back().state;
        out << "  iterations    = " << guess_one.history.size() << "\n";
        out << "  root          = [" << root[0] << ", " << root[1] << ", " << root[2] << "]\n";
    }

    out << "Initial guess 2 = [1, 0, 1]\n";
    out << "  status        = " << guess_two.status << "\n";
    if (!guess_two.history.empty()) {
        const auto &root = guess_two.history.back().state;
        out << "  iterations    = " << guess_two.history.size() << "\n";
        out << "  root          = [" << root[0] << ", " << root[1] << ", " << root[2] << "]\n";
    }

    // Compare the two final roots to answer the question directly.
    if (!guess_one.history.empty() && !guess_two.history.empty()) {
        const auto &r1 = guess_one.history.back().state;
        const auto &r2 = guess_two.history.back().state;
        const double distance = std::max({std::fabs(r1[0] - r2[0]),
                                          std::fabs(r1[1] - r2[1]),
                                          std::fabs(r1[2] - r2[2])});
        out << "  same root?    = " << (distance < 1.0e-6 ? "Yes" : "No") << "\n";
    }
}

void printSummary(const SystemSolveResult &guess_one, const SystemSolveResult &guess_two) {
    // Print the final roots to the console for quick verification.
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Question 2\n";
    if (!guess_one.history.empty()) {
        const auto &root = guess_one.history.back().state;
        std::cout << "  Guess [0.1, 1.2, 2.5] -> ["
                  << root[0] << ", " << root[1] << ", " << root[2] << "]"
                  << ", iterations = " << guess_one.history.size() << '\n';
    }
    if (!guess_two.history.empty()) {
        const auto &root = guess_two.history.back().state;
        std::cout << "  Guess [1, 0, 1] -> ["
                  << root[0] << ", " << root[1] << ", " << root[2] << "]"
                  << ", iterations = " << guess_two.history.size() << '\n';
    }
}

} // namespace

int main() {
    // Ensure the output directory exists before writing files.
    fs::create_directories("output");

    // Run Newton-Raphson for the two initial guesses given in the assignment.
    const auto guess_one = solveSystemWithNewton({0.1, 1.2, 2.5}, 1.0e-10, 100);
    const auto guess_two = solveSystemWithNewton({1.0, 0.0, 1.0}, 1.0e-10, 100);

    // Save the raw iteration history for both runs.
    writeCsv("output/q2_guess1.csv", guess_one);
    writeCsv("output/q2_guess2.csv", guess_two);

    // Save and print the concise summary required for the report.
    writeSummary("output/q2_summary.txt", guess_one, guess_two);
    printSummary(guess_one, guess_two);

    // Normal successful program exit.
    return 0;
}
