# Create and compare GIA models over Antarctica; 
:zap: **repo for the Antarctic GIA project, updated August 2016: watch this space  - a new publication to appear in JGR:Solid Earth is coming up soon!** :zap:

# What is GIA
For those who may not be familiar with the term (I wont blame you!), GIA stands for Glacio-Isostatic Adjustment. This is a response of the Earth's mantle to changes in the surface load. 
The Earth's mantle is a viscous-elastic material which responds very slowly to changes in the surface. For instance, after the last glaciation (~20.000 years ago), the large ice-sheets lying over northern Europe and north America melted away, causing a "rebound" (or an uplift trend) of the solid Earth, which, in fact, is still going on!.
We can still feel these changes nowadays, but GIA is a really difficult process to be measured, specially in Antarctica! 

There are a few attempts to either model or solve for GIA processes in Antarctica. This project attempts to find out how these solutions actually work by comparing the modelled uplift rates with a network of GPS data gathered in Antarctica. 
We also bring along a new solution for GIA as a result of our RATES project (https://sites.google.com/site/wwwratesantarcticanet/)


# Code
The routine takes in a raster dataset (GIA models are not provided in this repo), changes the values to make it binary (i.e. 1/0) and then out pops a shapefile which adopts the values of the input raster.

The mask passed in is a transformed version of the GIMP DEM as modified following [Sasgen et al. 2013] http://www.the-cryosphere.net/7/1499/2013/tc-7-1499-2013.pdf. 
The code framework will work for any mask, just change the path and the values which you need to alter to make it 1/0.
The code is taken largely from this post: [https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/](https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/)

Then it loads the uplift rates from the GPS network and performs the statistical assesment.
If running on windows, try using the polygonizer function as opposed to gdal_polygonizeR.

Also, this will only work if you have [gdal](http://www.gdal.org/) downloaded locally!

![GIA models in Antarctica](Fig2.jpg?raw=true "GIA models for Antarctica")
![How do GIA models in Antarctica actually perform - An assessment using GPS data](Fig3.jpg?raw=true "How do GIA models in Antarctica actually perform - An assessment using GPS data")