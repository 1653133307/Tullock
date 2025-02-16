from TullockContestSolver import TullockContestSolver

def read_input_from_file(file_name):
    """
    Read input data from a file.
    file_name: str, the name of the input file.
    Returns: tuple (n, R, a, r)
    """
    with open(file_name, 'r') as file:
        n = int(file.readline().strip())  # Number of players
        R = float(file.readline().strip())  # Reward
        a = list(map(float, file.readline().strip().split()))  # Efficiency parameters
        r = list(map(float, file.readline().strip().split()))  # Elasticity parameters
    return n, R, a, r
# Example Usage
import time

if __name__ == "__main__":
    file_name = "input_sample_3.txt"
    output_file = "output_sample_3.txt"  # Output file to store results
    with open(output_file, "w") as f_out:  # Open file in write mode
        f_out.write(f"Processing file: {file_name}\n\n")
        # Read input data from file
        n, R, a, r = read_input_from_file(file_name)
        print(f"Processing file: {file_name}")
        solver = TullockContestSolver(n, R, a, r)
        A_solutions = solver.solve()
        if isinstance(A_solutions, list):
            for A_star, s in A_solutions:
                f_out.write(f"Solution A_star: {A_star}\n")
                shares_and_efforts = solver.compute_shares_and_efforts(A_star, s)
                if shares_and_efforts is []:
                    continue
                for i, (share, effort) in enumerate(shares_and_efforts):
                    f_out.write(f"Player {i + 1}: Share = {share}, Effort = {effort}\n")
        else:
            f_out.write(f"{A_solutions}\n")
        f_out.write(f"\n\n")
    print(f"All results have been saved to {output_file}")