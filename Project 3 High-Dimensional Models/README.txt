README

For the third project on multidimensional models, the students were provided with data on the car market for N = 5 major European countries over T = 30 years, covering a total of J = 356 unique car types and number of purchases per car type.

The overarching research question of this empirical project was: Should a car manufacturer charge a higher price in the home market versus abroad. 

Students had to:
1) Choose a sample and discuss identification of the coefficient on price and the home dummy 
2) Motivate an OLS estimator that allows you to consistently estimate a discrete choice model of car purchases
3) Calculate the price-elasticity of demand and discuss home-market advantage

--------
Main file: P3_main_file_16.12.mlx 

Line 252 onwards display code produced by students as answers to questions in problem set.
Before line 252 most code was to some degree pre-written and students had to make choices of what variables/controls to keep or to drop. 
--------
Other files partly produced by students: 
ols.m, estimation_clog.m, estimation_clog_search.m, clogit_w8.m

Students had to complete the critical functions within each model: Among the tasks was to select the choice of optimizer (fminsearch vs. fminunc), specify the correct  variance matrices, specify the correct equations for estimation of parameters and standard errors.
--------
Auxiliare files provided by professor: read_data.m, print_par.m and others

