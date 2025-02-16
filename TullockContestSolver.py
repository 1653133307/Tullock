class TullockContestSolver:
    """
    A module for solving Tullock contest problems to compute Pure Nash Equilibrium (PNE)
    or epsilon-Approximate Nash Equilibrium (ε-NE) based on the contest instance.
    """
    def __init__(self, n, R, a, r):
        """
        Initialize the Tullock contest solver.
        n: int, number of contestants
        R: float, reward amount
        a: list of floats, efficiency parameters
        r: list of floats, elasticity parameters
        """
        self.n = n
        self.R = R
        self.a = a
        self.r = r

        self.r1 = min(self.r)  # Smallest r in the group
        self.r2 = min([ri for ri in self.r if ri > 1], default=float('inf'))  # Smallest r > 1
        self.r3 = min(self.r1 ** 2, self.r2 * (self.r2 - 1)) if self.r2 != float('inf') else self.r1 ** 2
        self.a1 = min(self.a)
        self.mmax = sum(ai * (self.R ** ri) for ai, ri in zip(self.a, self.r))
        self.eps = 1.0

    def solve(self):
        """
        Solve the contest based on the elasticity parameter ranges.
        Returns: The result based on the selected algorithm.
        """
        # Case 1: All players' r <= 1, use Algorithm 1
        if all(ri <= 1 for ri in self.r):
            return self._algorithm1()

        # Case 2: All players' r > 2, return "No PNE"
        if all(ri > 2 for ri in self.r):
            return "No PNE exists"

        # Case 3: Some players' r <= 1, others' r > 2, and no r in (1, 2]
        if all(ri <= 1 or ri > 2 for ri in self.r):
            return self._algorithm2()

        # Case 4: Remaining cases, use Algorithm 3
        return self._algorithm3()

    def _algorithm1(self):
        """
        Implementation of Algorithm 1.
        Returns: Result of Algorithm 1.
        """
        Y = []  # Initialize an empty list to store outputs
        A_max = self.mmax # Maximum value for binary search
        low, high = 0, A_max

        while high - low > 1e-6:
            A_star = (low + high) / 2
            S_A = self._compute_S(A_star)
            if abs(S_A - 1) < 1e-6:  # Check if S(A*) equals 1
                Y.append((A_star, self._compute_sigma(A_star))) # Append the PNE solution to Y
                break
            elif S_A < 1:
                high = A_star
            else:
                low = A_star

        return Y if Y else "Error"

    def _algorithm2(self):
        """
        Implementation of Algorithm 2.
        Returns: Result of Algorithm 2.
        """
        Y = []  # Initialize an empty list to store outputs
        I1 = [i for i in range(self.n) if self.r[i] <= 1]  # Players with r_i <= 1
        I2 = [i for i in range(self.n) if self.r[i] > 1]  # Players with r_i > 1

        # Case 1: Only players in I1 are active
        IA = I1
        A_max = self.n * self.R  # Maximum value for binary search
        low, high = max([self._compute_A_bar(i) for i in I1]), A_max

        while high - low > 1e-7:
            A_star = (low + high) / 2
            S_A = sum(self._compute_sigma_less_than_1(i, A_star) for i in IA)

            if abs(S_A - 1) < 1e-7:
                sig1 = self._compute_sigma(A_star)
                sig1 = [0 if self.r[i] > 2 else sig1[i] for i in range(self.n)]
                Y.append((A_star, sig1))
                break
            elif S_A < 1:
                low = A_star
            else:
                high = A_star

        # Case 2: Adding one player from I2
        for j in I2:
            IA = I1 + [j]
            low, high = self._compute_A_bar(j), self._compute_A_bar_upper(j)

            while high - low > 1e-7:
                A_star = (low + high) / 2
                S_A = sum(
                    self._compute_sigma_greater_than_1(i, A_star) if i == j else self._compute_sigma_less_than_1(i,
                                                                                                                 A_star)
                    for i in IA)
                if abs(S_A - 1) < 1e-5:
                    sig2 = self._compute_sigma(S_A)
                    sig2 = [sig2[i] if i in IA else 0 for i in range(self.n)]
                    Y.append((A_star, sig2))
                    break
                elif S_A < 1:
                    low = A_star
                else:
                    high = A_star

        return Y if Y else "No PNE exists"

    def _algorithm3(self, epsilon=1e-2):
        """
        Implementation of Algorithm 3 to search for ε-approximate solutions.
        Returns: A list of ε-approximate solutions.
        """
        Y = []  # Initialize an empty list to store outputs
        temp = []
        I2 = [i for i in range(self.n) if self.r[i] > 1]  # Players with r_i > 1

        # Step 1: Compute A_i_bar and A_i_bar_upper for all players in I2
        A1 = set()
        for i in I2:
            A1.add(self._compute_A_bar(i))
            A1.add(self._compute_A_bar_upper(i))

        # Step 2: Sort A1 in ascending order
        A1_sorted = sorted(A1)
        A2 = set()
        A_min, A_max = A1_sorted[0], A1_sorted[-1]
        L_sigma = A_max/A_min * self.r3
        L_p = ((A_max/self.a1) ** (1/self.r1 - 1)) * (1/self.r1)
        L = L_sigma*self.R + L_p * self.a1*(1+A_max*L_sigma)
        #epsilon = epsilon / (self.n * L)
        epsilon = epsilon / self.n
        current = A_min

        # Step 3: Binary search for ε-approximate solutions in boundary regions
        boundary_intervals = [(0, min(A_min, self.mmax)), (A_max, self.mmax)]
        for low, high in boundary_intervals:
            solution = self._binary_search_for_approx_solution(low, high, epsilon)
            if solution is not None:
                if low == 0:
                    Y.append((solution, self._compute_sigma(solution)))
                else:
                    sig = self._compute_sigma(solution)
                    sig = [0 if self.r[i] > 1 else sig[i] for i in range(self.n)]
                    Y.append((solution, sig))
        #("Y1",Y)
        if A_min > (self.mmax):
            return Y if Y else "No PNE exists"

        # Step 4: Generate intermediate points with δ spacing
        while A_max - (current * (1 + self.r3 * epsilon)) > 1e-5:
            current *= (1 + self.r3 * epsilon)
            A2.add(current)

        # Step 5: Combine A1 and A2
        A = sorted(A1.union(A2))
        #print(f"Size of A: {len(A)}")

        # Step 6: Verify each node in A using APPROX-SUBSET-SUM
        for A_star in A:
            if self._approx_subset_sum(A_star, epsilon):
                temp.append(A_star)

        for A_star in temp:
            res = self._approx_subset_sum2(A_star, epsilon)
            for t, A_solu, s_r in res:
                print(s_r)
                if t != 0 and (abs(sum(i for i in s_r)-1.0) < 0.1):
                    Y.append((t, s_r))

        return Y if Y else "No PNE exists"

    def _binary_search_for_approx_solution(self, low, high, epsilon=1e-6):
        """
        Perform binary search to find an ε-approximate solution.
        Returns: An approximate solution if found, else None.
        """
        while high - low > epsilon:
            mid = (low + high) / 2
            S_A = self._compute_S(mid)

            if abs(S_A - 1) < epsilon:
                return mid
            elif S_A < 1:
                high = mid
            else:
                low = mid

        return None

    def _compute_active_player_set(self, A):
        """
        Compute the active player set S0 based on the current value of A.
        A: float, the current aggregate action value.
        Returns: S0, a list of indices for active players.
        """
        S0 = []
        for i in range(self.n):
            if self.r[i] < 1:
                S0.append(i)
            elif self.r[i] == 1:
                if A < self.a[i] * self.R:
                    S0.append(i)
            elif self.r[i] > 1:
                if A < self._compute_A_bar(i):
                    S0.append(i)
        return S0

    def _compute_sigma(self, A):
        """
        Compute sigma values for all players based on the current value of A.
        A: float, the current aggregate action value.
        Returns: A list of sigma values for all players.
        """
        sigma = []
        for i in range(self.n):
            if self.r[i] < 1:
                sigma.append(self._compute_sigma_less_than_1(i, A))
            elif self.r[i] == 1:
                sigma.append(self._compute_sigma_equal_1(i, A))
            elif self.r[i] > 1:
                sigma.append(self._compute_sigma_greater_than_1(i, A))
        return sigma

    def _approx_subset_sum2(self, A, epsilon):
        """
        Verify if a subset of sigma combines with S0 to form an ε-approximate sum in (1 - ε, 1 + ε).
        A: float, the current aggregate action value.
        epsilon: float, the approximation parameter.
        Returns: A list of tuples (target value, ans[x]) where ans[x] is a list of sigma values per player.
        """
        S = self._compute_sigma(A)  # Compute sigma values
        ans_template = [0] * self.n  # Initialize ans list
        S0_indices = self._compute_active_player_set(A)

        # Step 1: Calculate S0 and modify S and ans_template
        S0 = 0
        ans = ans_template[:]  # Create a copy of ans_template
        for idx in S0_indices:
            S0 += S[idx]
            ans[idx] = S[idx]  # Store corresponding sigma value
            S[idx] = 0  # Remove used values from S

        # Step 2: Construct S_filtered and track original indices
        S_filtered_with_indices = [(S[i], i) for i in range(self.n) if S[i] != 0]
        S_filtered = [x[0] for x in S_filtered_with_indices]  # Extract values
        S_filtered_indices = [x[1] for x in S_filtered_with_indices]  # Extract original indices

        n = len(S_filtered)

        # Initialize tracking dictionaries
        X = [S0]
        ans_dict = {S0: ans[:]}

        Y = [S0]
        ans_dict_Y = {S0: ans[:]}

        for i in range(n):
            new_X = [x + S_filtered[i] for x in X]
            new_Y = [y + S_filtered[i] for y in Y]

            X = self._merge_list(X, new_X)
            Y = self._merge_list(Y, new_Y)

            # Update ans_dict

            for j in range(len(new_X)):
                new_x = new_X[j]
                ans_dict[new_x] = ans_dict[X[j]][:]  # **确保继承完整 subset**
                ans_dict[new_x][S_filtered_indices[i]] = S_filtered[i]  # **正确填充 sigma**

            for j in range(len(new_Y)):
                new_y = new_Y[j]
                ans_dict_Y[new_y] = ans_dict_Y[Y[j]][:]
                ans_dict_Y[new_y][S_filtered_indices[i]] = S_filtered[i]

                # Trim results
            X = self._trim_from_below(X, epsilon / (2 * n))
            Y = self._trim_from_above(Y, epsilon / (2 * n))

        valid_subsets = list(set(
            [(A, x, tuple(ans_dict[x])) for x in X if 1 - epsilon <= x <= 1 + epsilon] +
            [(A, y, tuple(ans_dict_Y[y])) for y in Y if 1 - epsilon <= y <= 1 + epsilon]
        ))

        return valid_subsets if valid_subsets else [[0,0,0]]  # Return target values and corresponding full ans arrays

    def _approx_subset_sum(self, A, epsilon):
        """
        Verify if a subset of sigma combines with S0 to form an ε-approximate sum in (1 - ε, 1 + ε).
        A: float, the current aggregate action value.
        epsilon: float, the approximation parameter.
        Returns: A list of tuples (target value, ans[x]) where ans[x] is a list of sigma values per player.
        """

        S0_indices = self._compute_active_player_set(A)
        S0 = sum(self._compute_sigma(A)[i] for i in S0_indices)
        sigma = [self._compute_sigma(A)[i] for i in range(self.n) if i not in S0_indices]

        n = len(sigma)
        X = [S0]
        Y = [S0]

        for i in range(n):
            X = self._merge_list(X, [x + sigma[i] for x in X])
            X = self._trim_from_below(X, epsilon / (2 * n))
            Y = self._merge_list(Y, [y + sigma[i] for y in Y])
            Y = self._trim_from_above(Y, epsilon / (2 * n))

        return list(set(x for x in X + Y if 1 - epsilon <= x <= 1 + epsilon))

    def _merge_list(self, L1, L2):
        """
        Merge two sorted lists without duplicates.
        Returns: Merged list.
        """
        return sorted(set(L1 + L2))

    def _trim_from_below(self, L, delta):
        """
        Trim a sorted list from below based on delta.
        Returns: Trimmed list.
        """
        L = sorted(L)
        last = L[0]
        trimmed = [last]

        for y in L[1:]:
            if y > (1 + delta) * last:
                trimmed.append(y)
                last = y

        return trimmed

    def _trim_from_above(self, L, delta):
        """
        Trim a sorted list from above based on delta.
        Returns: Trimmed list.
        """
        L = sorted(L, reverse=True)
        last = L[0]
        trimmed = [last]

        for y in L[1:]:
            if y < last * (1 - delta):
                trimmed.append(y)
                last = y

        return list(reversed(trimmed))

    def _compute_A_bar(self, i):
        """
        Compute the lower bound of A for player i when r_i > 1.
        """
        r_i = self.r[i]
        return self.a[i] * (self.R ** r_i) * ((r_i - 1) ** (r_i - 1)) / (r_i ** r_i)

    def _compute_A_bar_upper(self, i):
        """
        Compute the upper bound of A for player i when r_i > 1.
        """
        r_i = self.r[i]
        return self.a[i] * (self.R ** r_i) * ((r_i - 1) / r_i) ** (r_i - 1)

    def _compute_S(self, A_star):
        """
        Compute S(A*).
        A_star: float, the aggregate action being tested.
        Returns: Computed S(A*).
        """
        sigma_sum = 0

        for i in range(self.n):
            if self.r[i] < 1:
                sigma_sum += self._compute_sigma_less_than_1(i, A_star)
            elif self.r[i] == 1:
                sigma_sum += self._compute_sigma_equal_1(i, A_star)
            elif self.r[i] > 1:
                sigma_sum += self._compute_sigma_greater_than_1(i, A_star)

        return sigma_sum

    def _computer_threshold(self,i):
        from scipy.optimize import fsolve
        def equation(A_threshold):
            value = ((1e-3) * A_threshold) / self.a[i]
            if value <= 0:
                print("error")
                return float('inf')  # Avoid invalid power operation
            g_prime = value ** (1 / self.r[i] - 1)
            return (1 - (1e-3)) * self.R - g_prime * (A_threshold / (self.a[i] * self.r[i]))

        A_threshold_i_initial_guess = self.R / 2
        try:
            # Solve using fsolve
            A_threshold_i_solution = fsolve(equation, A_threshold_i_initial_guess, xtol=1e-3)[0]

            # Validate solution
            #if A_threshold_i_solution <= 0 or A_threshold_i_solution > R:
            #    raise ValueError("Solution out of range.")
        except Exception as e:
            #print(f"fsolve failed for player {i} with A_star={A_star}: {e}")
            A_threshold_i_solution = self.R / 2 # Default to midpoint value

        return A_threshold_i_solution
    def _compute_sigma_less_than_1(self, i, A_star):
        from scipy.optimize import newton
        """
        Compute sigma_i when r_i < 1.
        This assumes A_star > 0, ensuring sigma_i is always in (0, 1].
        """
        A_threshold = self._computer_threshold(i)

        from scipy.optimize import fsolve
        if abs(A_star - (self.mmax)) <= 1e-5 or abs(A_star - A_threshold) <= 1e-5:
            return 0
        
        def equation(sigma_i):
            value = sigma_i * A_star / self.a[i]
            if value <= 0:
                return float('inf')  # Avoid invalid power operation
            g_prime = value ** (1 / self.r[i] - 1)
            return (1 - sigma_i) * self.R - g_prime * (A_star/(self.a[i]*self.r[i]))
        
        # Iteratively adjust initial guesses
        sigma_i_solution = 0
        initial_guesses = [1e-4] + [0.1 * k for k in range(1, 10)] + [1 - 1e-4]

        for guess in initial_guesses:
            try:
                def wrapped_equation(sigma_i):
                    value = sigma_i * A_star / self.a[i]
                    if value <= 0:
                        raise ValueError("Negative value encountered")
                    g_prime = value ** (1 / self.r[i] - 1)
                    return (1 - sigma_i) * self.R - g_prime * (A_star/(self.a[i]*self.r[i]))

                sigma_i_solution = fsolve(wrapped_equation, guess, xtol=1e-4)[0]
                # Validate solution
                if 0 < sigma_i_solution <= 1:
                    return sigma_i_solution
            except ValueError:
                #print(f"Switching initial guess due to invalid value for player {i} with A_star={A_star}")
                continue  # Try the next initial guess
            except Exception as e:
                print(f"fsolve failed for player {i} with initial guess {guess} and A_star={A_star}: {e}")

        # If all attempts fail, return 0
        return 0

    def _compute_sigma_equal_1(self, i, A_star):
        """
        Compute sigma_i when r_i == 1.
        i: int, player index
        A_star: float, the aggregate action being tested.
        Returns: Computed sigma_i.
        """
        return max(0, 1 - (A_star / (self.a[i] * self.R)))

    def _compute_sigma_greater_than_1(self, i, A_star,epsilon = 1e-4):
        """
        Compute sigma_i when r_i > 1.
        """
        from scipy.optimize import minimize_scalar
        from scipy.optimize import brentq
        eps = epsilon / self.n
        r_i = self.r[i]
        A_i_bar_upper = self._compute_A_bar_upper(i)

        # Special cases
        if abs(A_star - A_i_bar_upper) < 1e-6:
            return (r_i - 1) / r_i  # Precision issues when A_star equals A_i_bar

        if A_star > A_i_bar_upper:
            return 0  # If A_star exceeds A_i_bar, sigma_i is 0

        def equation2(sigma_i):
            return (
                self.a[i] * (self.R ** r_i) * (r_i ** r_i) * ((1 - sigma_i) ** r_i) * (sigma_i ** (r_i - 1)) - A_star
            )

        # Solve using Brent’s method
        try:
            sigma_i_solution = brentq(equation2, (r_i - 1) / r_i, 1,xtol=1e-5, rtol=1e-5)
            #print(f"Brentq succeeded! Solution: {sigma_i_solution}")
            x2 = sigma_i_solution
        except ValueError:
            raise ValueError(f"Brent’s method failed for player {i} with A_star={A_star}")
        #print(x1,x2)
        return x2


    def _compute_best_effort(self, i, A_star, effort_i):
        from scipy.optimize import minimize_scalar
        def equation(x_i):
            return (
                    (self.a[i] * (x_i ** self.r[i])) /
                    (A_star - effort_i + self.a[i] * (x_i ** self.r[i])) * self.r[i] - x_i
            )

        # Maximizing equation → minimize the negative
        result = minimize_scalar(lambda x: -equation(x), bounds=(0, A_star), method='bounded', options={'xatol': 1e-6})

        if result.success:
            return result.x  # Return the optimal x_i
        else:
            raise ValueError("Optimization failed to find the best x_i")

    def compute_shares_and_efforts(self, A_star,s):
        """
        Compute the share and effort for each player given A_star.
        A_star: float, the aggregate action value.
        Returns: A list of (share, effort) tuples for each player.
        """
        results = []
        error = 0.0
        #print("s ",s)
        for i in range(self.n):
            action_i = s[i] * A_star  # Total action contributed
            effort_i = (action_i / self.a[i]) ** (1 / self.r[i])  # Effort derived from action and efficiency
            results.append((s[i], effort_i))
        for i in range(len(results)):
            if results[i][1] == 0:
                continue
            best_x_i = self._compute_best_effort(i, A_star,results[i][1])
            error = max(error, abs(best_x_i - results[i][0]))
        if (error < self.eps):
            return results
        else:
            return [],0
