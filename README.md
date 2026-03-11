# ME 5107 Assignment 04 Solution

This workspace now contains a complete C++ solution for the scanned assignment, split into separate C++ files for Question 1 and Question 2, plus the cleaned question transcription, generated result files, and a short report.

## Files

- `src/question1.cpp`: C++ solution for Question 1
- `src/question2.cpp`: C++ solution for Question 2
- `Makefile`: build and run commands
- `question.txt`: cleaned text version of the scanned assignment
- `report.md`: explanation report with numerical conclusions
- `plot_results.py`: creates the required plots from the generated CSV files
- `output/`: CSV histories, summary text, and plot images

## How to build and run

```bash
make
./question1
./question2
python3 plot_results.py
```

You can also use:

```bash
make run-q1
make run-q2
```

## What the code does

For Problem 1:

- evaluates the error function numerically with the trapezoidal rule using `N = 51`
- solves for the convection coefficient `h` using
  - bisection
  - Newton's method
  - a hybrid of 10 bisection steps followed by Newton's method
- stores iteration histories in CSV files

For Problem 2:

- solves the 3x3 nonlinear system using multivariable Newton-Raphson
- runs the solver from both initial guesses given in the question

## Generated outputs

After running the executables:

- `output/question1_summary.txt` contains the Question 1 answers
- `output/question2_summary.txt` contains the Question 2 answers
- `output/part1_*.csv` contain iteration-by-iteration histories for Problem 1
- `output/part2_*.csv` contain iteration-by-iteration histories for Problem 2

After running the plot script:

- `output/plots/q1a_bisection_history.png`
- `output/plots/q1b_newton_history.png`
- `output/plots/q1c_hybrid_history.png`

## Main results

- Bisection converges to `h ~= 200.196878 W/m^2-K`
- Newton's method with initial guess `5000` does not converge
- The hybrid method converges to `h ~= 200.196877 W/m^2-K`
- The two starting guesses in Problem 2 converge to different roots

See `report.md` for the explanation.
