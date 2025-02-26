# Genetic Algorithm Pipeline
## Based on Lam et. al. 2024.

Algorithm to figure out the optimal experimental sampling schedule.
Looks at the confidence interval for a parameter of interest across a number of randomly generated samples. The interval widths are calculated and used as a proxy for "fitness". The best ones are then passed on to the next "generation", with some noise added in all but the best case. Repeating in this way, the sampling schedule with the smallest parameter confidence bounds can be identified.

1. Generate M' sample schedules
2. Find PLB (profile likelihood bounds) for each
3. Save m' in M' for which PLB(m')>PLB(t') for all t' in M'
4. Order the set so: if l' >= m' in M' means PLB(l')>=PLB(m').
5. Carry forward the best through some weighted cloning process (now called M'')
6. Mutate all but the best with some normally distrivuted noise.
7. Repeat 2-6.

Constraint 1: m,n in M' means abs(m-n) > \delta{t}.
