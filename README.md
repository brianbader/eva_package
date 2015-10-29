## **To Do:** ##

* ~~Add better handling for the likelihood and data generation as shape -> 0~~

* Add interpolation option for AD and CVM tests

* ~~Add pgev function for the univariate case and export qgev as a function~~

* ~~Need to check MPS estimators for gpdfit function (gives weird results when shape < 0)~~

* Need to fix gevspatscore code (may crash sometimes since it relies on gev.pbscore)

* ~~Combine the ed and score sequential tests for GEV into one function~~

* Add gpd return levels (profile likelihood and delta method)

* Make sure likelihood density function is being referred to in all function

* Add more examples and descriptions

* Incorporate non-stationarity into the fitting and PB/ED tests