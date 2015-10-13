## **To Do:** ##

* Add better handling for the likelihood and data generation as shape -> 0

* Add pgev function for the univariate case and export qgev as a function

* Make sure outputs are dataframes, not matrices

* Need to check MPS estimators for gpdfit function (gives weird results when shape < 0)

* Need to fix gevspatscore code (may crash sometimes since it relies on gev.pbscore)

* Add pvalue adjustment function (give example from code, and cite paper)

* Combine the ed and score sequential tests for GEV into one function

* Change fortmax data to sealevel data (just create the top ten, with tau=60) to demonstrate the difference in CIs from r=1 to 10

* Add gpd return levels (profile likelihood and delta method)

* Incorporate non-stationarity into the fitting and PB/ED tests