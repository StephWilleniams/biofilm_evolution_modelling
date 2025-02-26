
1. Generate M' sample schedules
2. Find PLB (profile likelihood bounds) for each
3. Save m' in M' for which PLB(m')>PLB(t') for all t' in M'
4. Order the set so: if l' >= m' in M' means PLB(l')>=PLB(m').
5. Carry forward the best through some weighted cloning process (now called M'')
6. Mutate all but the best with some normally distrivuted noise.
7. Repeat 2-6.

Constraint 1: m,n in M' means abs(m-n) > \delta{t}.

To do:
- Make a function to carry out (2).
- Make a function to order the set for 3/4
- Determine ways in which the cloning can be done.
- Figure out the best number of repeats for (7), or decide a convergence condition.
- Figure out how to efficiently check constraint 1 in mutated schedules.
