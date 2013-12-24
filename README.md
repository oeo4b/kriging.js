kriging.js
==========

**kriging.js** is a Javascript library providing spatial prediction and mapping capabilities via the ordinary kriging algorithm. 

Kriging is a type of gaussian process where 2-dimensional coordinates are mapped to some target variable using kernel regression. This algorithm has been specifically designed to accurately model smaller data sets by assmuming a set of bayesian priors on both the variogram parameters and the gaussian process.

Fitting a model
---------------

The first step is to link **kriging.js** to your html code and assign your coordinate and target variables to 3 separate arrays.

``` html
<script src="kriging.js" type="text/javascript"></script>
<script type="text/javascript">
	var t = [ /* Target variable */ ];
	var x = [ /* X-axis coordinates */; ];
	var y = [ /* Y-axis coordinates */; ];
	var model = "exponential";
	var alpha = 100, beta = 1;
	var variogram = kriging.train(t, x, y, model, alpha, beta);
</script>
```

The train method in the kriging object fits your input to whatever variogram model you specify - gaussian, exponential or spherical - and returns a variogram object. 

Bayesian priors
---------------

Notice the alpha and beta variables, these correspond to the variance parameters of the gaussian priors of the variogram model and the gaussian process, respectively. A diffuse alpha (α) prior (~100) and standard normal variance beta (β) prior (1) is recommended; a formal mathematical explanation of the model is provided below.

Predicting new values
---------------------

Values can be predicted for new coordinate pairs by using the predict method in the kriging object.

``` javascript
  var xnew, ynew /* Pair of new coordinates to predict */;
  var tpredicted = kriging.predict(xnew, ynew, variogram);
  
```

Creating a map
--------------


Variogram and Probability Model
-------------------------------

The various variogram models can be interpreted as kernel functions for 2-dimensional coordinates **a**, **b** and parameters nugget, range, sill and A. Reparameterized as a linear function, with w = [nugget, (sill-nugget)/range], this becomes:
  
- Gaussian: k(**a**,**b**) = w[0] + w[1] * ( 1 - exp(-(||**a**-**b**||/range)^2 / A) )
- Exponential: k(**a**,**b**) = w[0] + w[1] * ( 1 - exp(-(||**a**-**b**||/range) / A) )
- Spherical: k(**a**,**b**) = w[0] + w[1] * ( 1.5 * (||**a**-**b**||/range) - 0.5 * (||**a**-**b**||/range)^3 )

The variance parameter (α) of the prior distribution for w should be manually set, according to:

- w ~ N(w|**0**, α**I**)

Using the fitted kernel function hyperparameters and setting K as the Gram matrix, the prior and likelihood for the gaussian process become:

- **y**       ~ N(**y**|**0**, K)
- **t**|**y** ~ N(**t**|**y**, β**I**)

The variance parameter (β) of the likelihood reflects the error in the gaussian process and should be manually set. 




