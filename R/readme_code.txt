Folder contains R scripts used in this study to replicate results.

#Step1:Decomposition of NDVI time-series of random MODIS forested pixels over Spain downloaded and pre-processed throught Google earth engine.
#Step2:Extraction of SPEI time-series from SPEI-scale dataset in .nc format and mapping SPEI change over spain and NDVI time-series.
#Step3: Multiple linear regression with NDVI, SPEI and TIME computing short-term (CGW) and long-term climate change effect (CCI) in the continuous analysis of forest canopy resilience. 
#Step3.5: Determining average water balance (WB) of Spain computed from GLDAS 2.1 environmental data (WB=precipitation-potential evapotranpiration). Repetition of step 3 using WB instead of SPEI (sensitivity analysis).
#Step4: Decomposition of NDVI time-series running DBEST package to extract short-term resilience metrics: sens=ndvi,lai loss; recov=ndvi,lai gain; resil=ndvi,lai gain/loss.
#Step5: Distinguishing protected areas from non-protected areas. Extracting landcover types from Corine dataset 1990, 2000. Looking at species groups.
#Step6: Resilience metrics analysis alongside environmental data using spatial auto-correlation regression models.

For more information on methods and for citing the article.
Khoury S, Coomes DA. Resilience of Spanish forests to recent droughts and climate change. Glob Change Biol. 2020;00:1–20. https://doi.org/10.1111/gcb.1526817