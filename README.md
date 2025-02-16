# Tullock Contest Solver

This project implements a solver for **Tullock contests**, computing **Pure Nash Equilibrium (PNE)** and **ε-Approximate Nash Equilibrium (ε-NE)** based on contest instances. The solver applies different algorithms based on the elasticity parameters of the contestants.

## **Features**

- Computes **PNE** and **ε-NE** in Tullock contests.
- Handles different cases based on elasticity parameters (`r` values).
- Uses **binary search and subset sum approximation** techniques.
- Reads input from files and writes output to specified files.

## **Installation & Environment**

### **Requirements**

Ensure you have the following dependencies installed:

- Python **3.8+**
- Required libraries:
  - `numpy`
  - `scipy`

### **Installation**

Clone the repository and install dependencies:

```
git clone https://github.com/1653133307/Tullock.git
cd Tullock
pip install -r requirements.txt
```

Alternatively, install required libraries manually:

```
pip install numpy scipy
```

## **Usage**

### **Run the Solver**

You can run the solver by executing:

```
python Solver.py
```

Make sure your input file is in the correct format (see below).

## **Input Format**

The input file should be a **text file (**`**.txt**`**)** with the following structure:

```
n
R
a_1 a_2 ... a_n
r_1 r_2 ... r_n
```

Where:

- `n` → **Number of contestants** (integer)
- `R` → **Reward amount** (float)
- `a_1, a_2, ..., a_n` → **Efficiency parameters** (list of floats)
- `r_1, r_2, ..., r_n` → **Elasticity parameters** (list of floats)

### **Example Input (**input_sample_3.txt)

```
5
10
4.332009257737827 4.122455563544566 2.7494592222802523 4.262567021416704 4.053124274123181
1.9957465259789928 1.5590390196424875 1.8678656240993055 1.9391436075978281 1.9265985561604133
```

## **Output Format**

The solver writes results to an output file. Each equilibrium solution is printed in the format:

```
Solution A_star: <A_star_value>
Player 1: Share = <share_1>, Effort = <effort_1>
Player 2: Share = <share_2>, Effort = <effort_2>
...
```

If no PNE exists, it returns:

```
No PNE
```

### **Example Output (output_sample_3.txt)**

```
Processing file: input_sample_3.txt

No PNE exists
```

## **File Structure**

```
Tullock/
│── Solver.py                  # Main script to run the solver
│── TullockContestSolver.py     # Core logic for computing equilibria
│── input_sample_1.txt          # Example input file
│── output_sample_1.txt         # Example output file
│── README.md                   # Documentation
│── requirements.txt            # Python dependencies
```

## **How It Works**

1. The program **reads input** from a text file.
2. It initializes the **TullockContestSolver** with:
   - `n`: number of contestants
   - `R`: reward amount
   - `a`: efficiency parameters
   - `r`: elasticity parameters
3. Based on `r` values, it applies:
   - **Algorithm 1** if all `r ≤ 1`
   - **Algorithm 2** if all`r ≤ 1` or `r > 2`
   - **Algorithm 3** otherwise (ε-approximate Nash Equilibrium)
4. The solver **computes the equilibrium solution** and writes it to the output file.

## **License**

This project is open-source under the **MIT License**.

## **Contact & Contributions**

For questions or contributions, feel free to submit an issue or pull request on [**GitHub**](https://github.com/1653133307/Tullock).

🚀 **Happy Computing!**