
What could this project be?

- Many bacteria must exist in multiple distinct settings.
- Each of these settings contain different sets of challenges.
- It is not well known how optimising for setting A, can improve or impede viability of organisms in setting B.
- In order to explore this idea, we describe an e

- An example of this kind of setting optimisation in Seiler 
- This was a population which were encouraged by the presence of predators to form opposing (non-edible) phenotypes
- To extend this, we have begun to think about experiments in which these experiments are serially continued along some large time-frame to reinforce these behaviours
- (1) By fitting the data to mathematical models, which are designed to 

- re-wrote a bunch of the code in python to get a better feel for what it's doing.
- Realised there was some issue in the way I was trying to do the parameter estimation for the synthetic data. 
- In particular, I was fitting to a normal noise, not a log normal. 
- As well as this, I was using the residual of the data, not the ratio. I'm not entirely sure why the ratio is the right thing to do, I will look at the paper from Ruth/Simp to see what they say about this.
- The ratio of the data to the original is now IID, which is good.
- However, the % difference between the predicted params vs. the originals doesn't seem to be reducible with increasing data in many cases. This could be an indication that they are not practically identifiable? 
