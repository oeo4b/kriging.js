kriging.js
==========

**kriging.js** is a Javascript library providing spatial prediction and mapping capabilities via the ordinary kriging algorithm. 

Kriging is a type of gaussian process where 2-dimensional coordinates are mapped to some target (*f(x,y) => t*) using kernel regression. This algorithm has been specifically designed to accurately model smaller data sets by assmuming a set of bayesian priors on both the variogram parameters and the gaussian process.

Fitting a model
---------------

The first step is to link **kriging.js** to your html code and assign your coordinates and target variable to 3 separate arrays.

``` html
<script src="kriging.js" type="text/javascript"></script>
<script type="text/javascript">
	var t = [ /* Target variable */ ];
	var x = [ /* X-axis coordinates */; ];
	var y = [ /* Y-axis coordinates */; ];
	var model = kriging.train(t, x, y, "gaussian", 100, 100);
</script>
```


Creating a map
--------------
