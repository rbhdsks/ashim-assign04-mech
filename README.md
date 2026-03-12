# ME 5107 Assignment 04 Solution

This workspace now contains a complete C++ solution for the scanned assignment, split by the exact assignment numbering: `1(a)`, `1(b)`, `1(c)`, and `2`, plus the cleaned question transcription, generated result files, and a short report.

## Files

- `src/q1a.cpp`: C++ solution for Question 1(a)
- `src/q1b.cpp`: C++ solution for Question 1(b)
- `src/q1c.cpp`: C++ solution for Question 1(c)
- `src/question2.cpp`: C++ solution for Question 2
- `Makefile`: build and run commands
- `question.txt`: cleaned text version of the scanned assignment
- `report.md`: explanation report with numerical conclusions
- `plot_results.py`: creates the required plots from the generated CSV files
- `output/`: CSV histories, summary text, and plot images

## How to build and run

```bash
make
./q1a
./q1b
./q1c
./q2
python3 plot_results.py
```

You can also use:

```bash
make run-q1a
make run-q1b
make run-q1c
make run-q2
```

## What the code does

For Question 1(a), 1(b), and 1(c):

- evaluates the error function numerically with the trapezoidal rule using `N = 51`
- solves for the convection coefficient `h` separately using
  - bisection in `q1a.cpp`
  - Newton's method in `q1b.cpp`
  - the hybrid method in `q1c.cpp`
- stores each method history in its own CSV file

For Question 2:

- solves the 3x3 nonlinear system using multivariable Newton-Raphson
- runs the solver from both initial guesses given in the question

## Generated outputs

After running the executables:

- `output/q1a_summary.txt` contains the Question 1(a) answer
- `output/q1b_summary.txt` contains the Question 1(b) answer
- `output/q1c_summary.txt` contains the Question 1(c) answer
- `output/q2_summary.txt` contains the Question 2 answer
- `output/q1a_bisection.csv`, `output/q1b_newton.csv`, and `output/q1c_hybrid.csv` contain the scalar iteration histories
- `output/q2_guess1.csv` and `output/q2_guess2.csv` contain the Question 2 iteration histories

After running the plot script:

- `output/plots/q1a_bisection_history.png`
- `output/plots/q1b_newton_history.png`
- `output/plots/q1c_hybrid_history.png`

## Main results

- Question 1(a) bisection converges to `h ~= 200.196878 W/m^2-K`
- Question 1(b) Newton with initial guess `5000` does not converge
- Question 1(c) hybrid converges to `h ~= 200.196877 W/m^2-K`
- The two starting guesses in Problem 2 converge to different roots

See `report.md` for the explanation.
