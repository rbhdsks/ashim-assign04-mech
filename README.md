# ME 5107 Assignment 04 Solution

This workspace now contains a complete C++ solution for the scanned assignment, plus the cleaned question transcription, generated result files, and a short report.

## Files

- `src/assignment04.cpp`: main C++ solution
- `Makefile`: build and run commands
- `question.txt`: cleaned text version of the scanned assignment
- `report.md`: explanation report with numerical conclusions
- `plot_results.py`: creates the required plots from the generated CSV files
- `output/`: CSV histories, summary text, and plot images

## How to build and run

```bash
make
./assignment04
python3 plot_results.py
```

You can also use:

```bash
make run
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

After running the executable:

- `output/summary.txt` contains the main numerical answers
- `output/part1_*.csv` contain iteration-by-iteration histories for Problem 1
- `output/part2_*.csv` contain iteration-by-iteration histories for Problem 2

After running the plot script:

- `output/plots/bisection_history.png`
- `output/plots/newton_history.png`
- `output/plots/hybrid_history.png`

## Main results

- Bisection converges to `h ~= 200.196878 W/m^2-K`
- Newton's method with initial guess `5000` does not converge
- The hybrid method converges to `h ~= 200.196877 W/m^2-K`
- The two starting guesses in Problem 2 converge to different roots

See `report.md` for the explanation.
