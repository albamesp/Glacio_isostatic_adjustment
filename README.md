# Create and compare GIA models over Antarctica; 
:zap: **repo for the Antarctic GIA project, updated August 2016: watch this space  - a new publication to appear in JGR:Solid Earth is coming up soon!** :zap:

# Code
The routine takes in a raster dataset (GIA models are not provided in this repo), changes the values to make it binary (i.e. 1/0) and then out pops a shapefile which adopts the values of the input raster.

The mask passed in is a transformed version of the GIMP DEM as modified following [Sasgen et al. 2013] http://www.the-cryosphere.net/7/1499/2013/tc-7-1499-2013.pdf. 
The code framework will work for any mask, just change the path and the values which you need to alter to make it 1/0.

The code is taken largely from this post: [https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/](https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/)

If running on windows, try using the polygonizer function as opposed to gdal_polygonizeR.

Also, this will only work if you have [gdal](http://www.gdal.org/) downloaded locally!

![GIA models in Antarctica](Fig2.jpg?raw=true "GIA models for Antarctica")
