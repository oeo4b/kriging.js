kriging.js
==========

kriging.js is an HTML5 implementation of an ordinary kriging algorithm using the canvas element. This library is useful for spatial predictions and visualization of geographical data.

Usage
-----

This library can be used by instantiating an object of class kriging, and supplying the following arguments to the contructor and other methods.

``` js
var x = new kriging("mycanvas");
x.krig(longitude, latitude, response, polygons);
x.map(center, zoom);
```

That's it!

License
-------

Copyright 2012

