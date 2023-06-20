# landmulti

Executable codes of the statistical manuscript:

"Enhancing Long-term Survival Prediction with Multiple Short-term Events: Landmarking with A Flexible Varying Coefficient Model"

1. main.R: R file contains wrapper functions used to implement the proposed method. It shows how to analyze the example dataset "mydata.csv".

2. HelperFunctions.R: helper functions that are used in main.R.

3. mydata.csv: example data file that contains columns:
   
                                        (1) time: long-term survival time, Y;
   
                                        (2) outcome: binary censoring indicator with 0 indicates censoring and 1 indicates event, \delta;
   
                                        (3) st1: event time for short-term outcome 1, S1;
   
                                        (4) st2: event time for short-term outcome 2, S2;
   
                                        (5) age: continuous covariate, age, X;
   
                                        (6) ID: subject id.

