kriging.js
==========

kriging.js is an HTML5 implementation of an ordinary kriging algorithm using the canvas element. This library is useful for spatial predictions and visualization of geographical data.

Usage
-----

This library can be used by first instantiating an object of class kriging and supplying the following arguments to the constructor and other methods:

``` javascript
var x = new kriging("mycanvas");
x.krig(longitude, latitude, response, polygons);
x.map(center, zoom);
```

Make sure you've created a canvas element somewhere in your page before you start the map:

``` html
<canvas id="mycanvas"></canvas>
```

All the supplied arguments should be 1-dimensional Array objects and `zoom` should be a positive value. The `polygons` argument should be a multi-dimensional Array object. For polygons 1 through n, each x and y Array object refers to the longitude and latitude, respectively, coordinate vertices for that particular polygon:

``` javascript
var polygons = [[x1, y1], [x2, y2], ..., [xn, yn]];
```

License
-------

Copyright 2012

