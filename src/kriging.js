/**
 *
 * kriging.js
 *
 * Copyright 2012
 */


/* Extend the Array class */
Array.prototype.max = function() {
  return Math.max.apply(null, this)
}

Array.prototype.min = function() {
  return Math.min.apply(null, this)
}

Array.prototype.mean = function() {
  for(var i=0, sum=0; i<this.length;i++) {
    sum += this[i];
  }
  return sum / this.length;
}

/**
 * Ported R functions
 */
/* Repeat a value */
var R_rep = function(x, times) {
  var i = new Array(times);
  for(var j=0; j<i.length; j++) {
    i[j] = x;
  }
  return i;
}

/* Matrix transpose */
var R_t = function(x) {
  /* Must be a 2-dimensional matrix */
  var i, j, n, m;
  n = x.length;
  m = x[0].length;
  
  var y = new Array(m);
  for(i=0;i<m;i++) {
    y[i] = new Array(n);
    for(j=0;j<n;j++) {
	y[i][j] = x[j][i];
    }
  }
  return y;
}


/* Determinant */
var R_det = function(x, n) {
  var i, j, k, l;
  var det = 0;
  var m = new Array(n-1);
  for(i=0;i<(n-1);i++) {
    m[i] = new Array(n-1);
  }

  if(n<1) return;
  else {
    if(n==1) det = x[0][0];
    else {
      if(n==2) det = x[0][0]*x[1][1] - x[1][0]*x[0][1];
      else {
        det = 0;
        for(i=0;i<n;i++) {
          for(j=1;j<n;j++) {
	    k=0;
            for(l=0;l<n;l++) {
	      if(l==i) continue;
              m[j-1][k] = x[j][l];
              k++;
            }
          }
          det += Math.pow(-1, i+2) * x[0][i] * R_det(m, n-1);
        }
      }
    }
    return det;
  }
}

/* Non-R function -- essential for R_solve */
var cofactor = function(x, n) {
  var i, j, k, l, m, o;
  var det;
  var c = new Array(n-1);
  var y = new Array(n);

  for(i=0;i<n;i++) y[i] = new Array(n);
  for(i=0;i<(n-1);i++) c[i] = new Array(n-1);
  for(i=0;i<n;i++) {
    for(j=0;j<n;j++) {
      k=0;
      for(l=0;l<n;l++) {
	if(l==j) continue;
        m = 0;
        for(o=0;o<n;o++) {
	  if(o==i) continue;
          c[k][m] = x[l][o];
          m++;
        }
        k++;
      }
      det = R_det(c, n-1);
      y[j][i] = Math.pow(-1, j+i+2) * det;
    }
  }
  return y;
}

/* Matrix inversion */
var R_solve = function(x) {
  /* Solve to determine the adjunct matrix */
  var i, j;
  var adj = R_t(cofactor(x, x.length));
  var inv_det_a = 1 / R_det(x, x.length);
  var y = new Array(x.length);

  for(i=0;i<x.length;i++) {
    y[i] = new Array(x.length);
    for(j=0;j<x.length;j++) {
      y[i][j] = inv_det_a * adj[i][j];
    }
  }

  return y;
}

/* Fit a linear model */
var R_lm = function(y, x) {
  var n = y.length;

  /* Add an intercept term to the design matrix */
  x = [R_rep(1, n), x];
  y = [y];

  /* OLS estimate */
  return matrixmult(matrixmult(R_solve(matrixmult(R_t(x), x)), R_t(x)), y);
}

/* Cluster analysis */
var R_kmeans = function(x, y, centers) {

}

/**
 * Matrix multiplication 
 */
var matrixmult = function(y, x) {
  var i, j, k;
  var n = x.length;
  var m = x[0].length;
  if(m!=y.length) return false;
  var p = y[0].length;
  var z = new Array(n);

  for(i=0;i<n;i++) {
    z[i] = new Array(p);
    for(j=0;j<p;j++) {
      z[i][j] = 0;
      for(k=0;k<m;k++) {
        z[i][j] += x[i][k] * y[k][j]; 
      }
    }
  }
  return z;
}

/* Point-in-polygon */
var pip = function(X, Y, x, y) {
  var i, j;
  var c = false;
  for(i=0, j=X.length-1; i<X.length; j = i++) {
    if( ((Y[i]>y) != (Y[j]>y)) && (x<(X[j]-X[i]) * (y-Y[i]) / (Y[j]-Y[i]) + X[i]) ) {
      c = !c;
    }
  }
  return c;
}


/**
 * Defines the kriging class
 * 
 */
function kriging(id) {
  /* Output testing */
  var o = document.getElementById("output");

  /* Global vars */
  var canvaspad = 50;
  var pixelsize = 8;
  var yxratio = 1;

  /* Canvas element */
  var canvasobj = document.getElementById(id);
  this.canvas = document.getElementById(id);
  this.canvas.ctx = this.canvas.getContext("2d");

  /* New objects */
  this.canvas.model = new Object();

  /* Kriging method 
   * Usage: kriging(longitude, latitude, response, polygons)
   */
  this.krig = function(x, y, response, polygons) {
    /* Bring the polygons and frame properties into the DOM */
    this.canvas.model.x = x;
    this.canvas.model.y = y;
    this.canvas.model.response = response;
    this.canvas.model.response_min = response.min();
    this.canvas.model.response_max = response.max();
    this.canvas.model.response_range = response.max() - response.min();
    this.canvas.polygons = polygons;
   
    /**
     * Calculate the euclidean distance matrix for the coordinates 
     * and the outcome variable
     */
    this.canvas.model.n = response.length;
    var i, j, k;
    var D = new Array(this.canvas.model.n);
    var V = new Array(this.canvas.model.n);
    for(i=0; i<this.canvas.model.n; i++) {
      D[i] = new Array(this.canvas.model.n);
      V[i] = new Array(this.canvas.model.n);
      for(j=0; j<this.canvas.model.n; j++) { 
        D[i][j] = Math.sqrt(Math.pow(x[i]-x[j], 2) + Math.pow(y[i]-y[j], 2));
        V[i][j] = Math.abs(response[i] - response[j]);
      }
    }

    /* Fit the observations to the variogram */
    var lags = 10;
    var sum_z, n_h;
    this.canvas.model.semivariance = new Array();
    this.canvas.model.distance = new Array();
    var cutoff = Math.sqrt(Math.pow(x.max() - x.min(), 2) + Math.pow(y.max() - y.min(), 2)) / 3;
    for(i=0; i<lags; i++) {
      sum_z = 0;
      n_h = 0;
      for(j=0; j<this.canvas.model.n; j++) {
	  for(k=j+1; k<this.canvas.model.n; k++) {
	    if(D[j][k] <= ((i+1)*cutoff/lags)) {
	      sum_z += Math.pow(V[j][k], 2);
            n_h++;
	    }
	  }
      }
      if(!isNaN(sum_z / n_h)) {
        this.canvas.model.semivariance.push(sum_z / n_h);
        this.canvas.model.distance.push((i+1)*cutoff/lags);
      }
    }

    /* Check for enough points in the lag model */
    if(this.canvas.model.semivariance.length<3) {
      /* ERROR -- quit app */
    }

    /* Estimate the model parameters */
    var coef = R_lm(this.canvas.model.semivariance, this.canvas.model.distance);
    this.canvas.model.nugget = coef[0][0]; /* Intercept */
    this.canvas.model.range = this.canvas.model.distance.max();
    this.canvas.model.sill = coef[0][1] * this.canvas.model.range;

    /**
      * Calculate the inverted (n+1) x (n+1) matrix
      * Used to calculate weights
      */
    var X = new Array(this.canvas.model.n+1);
    for(i=0;i<=this.canvas.model.n;i++) {
      X[i] = new Array(this.canvas.model.n+1);
      for(j=0;j<=this.canvas.model.n;j++) {
        if(i==this.canvas.model.n && j!=this.canvas.model.n) X[i][j] = 1;
        else {
          if(i!=this.canvas.model.n && j==this.canvas.model.n) X[i][j] = 1;
          else {
            if(i==this.canvas.model.n && j==this.canvas.model.n) X[i][j] = 0;
            else {
              X[i][j] = this.canvas.model.spherical(D[i][j]);
            }
          }
        }
      }
    }
    
    /* Invert the matrix */
    this.canvas.model.X_inv = R_solve(X);
  }

  /* Variogram models */
  this.canvas.model.exponential = function(h) {
    if(h==0) return 0;
    else {
      return this.nugget + (this.sill - this.nugget) * (1 - Math.exp((-3*Math.abs(h)) / this.range));
    }
  }

  this.canvas.model.spherical = function(h) {
    if(h>this.range) return this.sill;
    if(h<=this.range && h>0) {
      return this.nugget + (this.sill - this.nugget) * ((3*h)/(2*this.range) - Math.pow(h, 3) / (2*Math.pow(this.range, 3)));
    }
    else return 0;
  }

  /* Model prediction method */
  this.canvas.model.pred = function(x, y) {
    var i;
    var L = R_rep(1, this.n+1);
    for(i=0;i<this.n;i++) {
      L[i] = Math.sqrt(Math.pow(this.x[i] - x, 2) + Math.pow(this.y[i] - y, 2))
    }
    var R = matrixmult(this.X_inv, [L])[0];
    R.pop();
    return matrixmult(R_t([R]), [this.response])[0][0];
  }

  /**
   * Set up the map properties, event handlers and initialize the map.
   */
  this.map = function(center, zoom) {
    /* Set up the canvas frame */
    this.canvas.height = window.innerHeight - this.canvas.offsetTop - canvaspad;
    this.canvas.style.border = "";


    /**
     * Loop through the polygons to determine the limits based on the 
     * area of each of the polygons.
     * AND
     * Create an Array containing the center coordinates for each polygon 
     * to be used during the sorting algorithm. 
     */
    var i, j;
    this.canvas.polygoncenters = new Array(this.canvas.polygons.length);
    this.canvas.polygonsorted = new Array(this.canvas.polygons.length);

    for(i=0;i<this.canvas.polygons.length;i++) {
      if(i==0) {
	this.canvas.xlim = [this.canvas.polygons[i][0].min(), this.canvas.polygons[i][0].max()];
      }
      else {
        if(this.canvas.polygons[i][0].min()<this.canvas.xlim[0]) this.canvas.xlim[0] = this.canvas.polygons[i][0].min();
        if(this.canvas.polygons[i][0].max()>this.canvas.xlim[1]) this.canvas.xlim[1] = this.canvas.polygons[i][0].max();
      }
      this.canvas.polygoncenters[i] = [this.canvas.polygons[i][0].mean(), this.canvas.polygons[i][1].mean()];
      this.canvas.polygonsorted[i] = 0;
    }


    /**
     * Calculate the ratio and pixel size for conversion 
     * between units.
     */
    this.canvas.xratio = (this.canvas.xlim[1]-this.canvas.xlim[0]) / this.canvas.width;
    this.canvas.yratio = this.canvas.xratio * yxratio;

    this.canvas.xlim = [center[0] - 0.5 * this.canvas.width * this.canvas.xratio, center[0] + 0.5 * this.canvas.width * this.canvas.xratio];
    this.canvas.ylim = [center[1] - 0.5 * this.canvas.height * this.canvas.yratio, center[1] + 0.5 * this.canvas.height * this.canvas.yratio];

    this.canvas.xpixel = pixelsize * this.canvas.xratio;
    this.canvas.ypixel = pixelsize * this.canvas.yratio;


    /* Mouse event properties */
    this.canvas.mousex = 0;
    this.canvas.mousey = 0;
    this.canvas.startx = 0;
    this.canvas.starty = 0;
    this.canvas.mousedown = false;

    /* Resize event handlers */
    window.onresize = function(e) {
      canvasobj.height = window.innerHeight - canvasobj.offsetTop - canvaspad;
      canvasobj.ylim[0] = canvasobj.ylim[1] - canvasobj.height * canvasobj.yratio;
      canvasobj.render();
    }


    /**
     * Mouse event handlers  
     */
    document.addEventListener("mousemove", function(e) {
      /* Reset mouse coordinates */
      canvasobj.mousex = e.pageX - canvasobj.offsetLeft;
      canvasobj.mousey = e.pageY - canvasobj.offsetTop;
      var predicted = canvasobj.model.pred(canvasobj.xlim[0] + canvasobj.xratio*canvasobj.mousex, canvasobj.ylim[0] - canvasobj.yratio*canvasobj.mousey);
      document.getElementById("output").innerHTML = "Predicted value: " + predicted;

     /* Drag the map if mouse is clicked */
      if(canvasobj.mousedown) {
	//document.onselectstart = function() { return false; }

        canvasobj.xlim = [canvasobj.xlim[0] + canvasobj.xratio*(canvasobj.startx-canvasobj.mousex), canvasobj.xlim[1] + canvasobj.xratio*(canvasobj.startx-canvasobj.mousex)]; 
        canvasobj.ylim = [canvasobj.ylim[0] - canvasobj.yratio*(canvasobj.starty-canvasobj.mousey), canvasobj.ylim[1] - canvasobj.yratio*(canvasobj.starty-canvasobj.mousey)]; 
        canvasobj.startx = canvasobj.mousex;
        canvasobj.starty = canvasobj.mousey;
        canvasobj.render();
      }
    });

    this.canvas.addEventListener("mousedown", function(e) {
      this.startx = this.mousex;
      this.starty = this.mousey;
      this.mousedown = true;
    });

    document.addEventListener("mouseup", function(e) {
      canvasobj.mousedown = false;
    });


    /**
     * Navigation event handlers (zoom-in / zoom-out buttons)
     */
    var zoomin = document.createElement("input");
    var zoomout = document.createElement("input");
    zoomin.style.position = "absolute";
    zoomout.style.position = "absolute";
    zoomout.style.top = zoomin.size*2;
    zoomin.type = "button";
    zoomout.type = "button";
    zoomin.value = "+";
    zoomout.value = "~";
    zoomin.onclick = function() { canvasobj.zoom(0.8, yxratio, pixelsize);}
    zoomout.onclick = function() { canvasobj.zoom(1.2, yxratio, pixelsize);}
    this.canvas.parentNode.insertBefore(zoomin, this.canvas);
    this.canvas.parentNode.insertBefore(zoomout, this.canvas);
    
    /* Start the map */
    this.canvas.zoom(zoom, yxratio, pixelsize);
  }


  /**
   * Navigation
   */
  this.canvas.zoom = function(zoom, yxratio, pixelsize) {
    /* Re-size the limits */
    var newlen = [zoom * (this.xlim[1]-this.xlim[0])/2,
                  zoom * (this.ylim[1]-this.ylim[0])/2];
    var center = [(this.xlim[1]-this.xlim[0])/2 + this.xlim[0], 
                  (this.ylim[1]-this.ylim[0])/2 + this.ylim[0]];

    /* Reset the properties */
    this.xlim = [center[0]-newlen[0], center[0]+newlen[0]];
    this.ylim = [center[1]-newlen[1], center[1]+newlen[1]];
    this.xratio = (this.xlim[1]-this.xlim[0]) / this.width;
    this.yratio = this.xratio * yxratio;
    this.xpixel = pixelsize * this.xratio;
    this.ypixel = pixelsize * this.yratio;

    /* Render the map */
    this.render();
  }


  /** 
   * Methods for drawing onto the canvas
   */

  /* Color spectrums */
  this.canvas.colorspectrum = new Object();
  this.canvas.colorspectrum.heatcolors = ["#FF0000","#FF0E00","#FF1C00","#FF2A00","#FF3900","#FF4700","#FF5500","#FF6300","#FF7100","#FF8000","#FF8E00","#FF9C00","#FFAA00","#FFB800","#FFC600","#FFD500","#FFE300","#FFF100","#FFFF00","#FFFF15","#FFFF40","#FFFF6A","#FFFF95","#FFFFBF","#FFFFEA"];
  this.canvas.colorspectrum.terraincolors = ["#00A600", "#10AC00", "#20B100", "#32B700", "#45BD00", "#59C300", "#6DC900", "#83CE00", "#9AD400", "#B2DA00", "#CBE000", "#E6E600", "#E6D612", "#E7C924", "#E8BF36", "#E9B848", "#EAB35A", "#EBB16D", "#ECB27F", "#EDB592", "#EEBCA5", "#EFC5B8", "#F0D1CB", "#F1E0DF", "#F2F2F2"];
  this.canvas.colorspectrum.topocolors = ["#4C00FF", "#2600FF", "#0000FF", "#0026FF", "#004CFF", "#0073FF", "#0099FF", "#00BFFF", "#00E5FF", "#00FF4D", "#00FF21", "#0BFF00", "#37FF00", "#62FF00", "#8EFF00", "#BAFF00", "#E6FF00", "#FFFF00", "#FFF219", "#FFE833", "#FFE04D", "#FFDC66", "#FFDB80", "#FFDC99", "#FFE0B3"];
  this.canvas.colorspectrum.cmcolors = ["#80FFFF", "#8AFFFF", "#95FFFF", "#9FFFFF", "#AAFFFF", "#B5FFFF", "#BFFFFF", "#CAFFFF", "#D4FFFF", "#DFFFFF", "#EAFFFF", "#F4FFFF", "#FFFFFF", "#FFF4FF", "#FFEAFF", "#FFDFFF", "#FFD5FF", "#FFCAFF", "#FFBFFF", "#FFB5FF", "#FFAAFF", "#FF9FFF", "#FF95FF", "#FF8AFF", "#FF80FF"];

  this.canvas.render = function() {
    this.clear();
    this.background();
  }  

  this.canvas.pixel = function(x, y, col) {
    this.ctx.fillStyle = col;
    this.ctx.fillRect((x-this.xlim[0])/this.xratio - pixelsize/2 + 1, this.height - (y-this.ylim[0])/this.yratio - pixelsize/2 + 1, pixelsize - 2, pixelsize - 2);
  }
  
  this.canvas.clear = function() {
    this.ctx.clearRect(0, 0, this.width, this.height);
  }


  /* Fills the background with the polygons */
  this.canvas.background = function() {
    /**
     * 1) Nearest-neighbor
     * 2) Loop through sorted polygons 
     * 3) Point-in-polygon for eligible points
     * 4) Break if no points in polygon
     */

    var i, j, k;
    for(i=0; i<this.polygoncenters.length; i++) {
      this.polygonsorted[i] = Math.sqrt(Math.pow(this.polygoncenters[i][0] - this.xlim.mean(), 2) + Math.pow(this.polygoncenters[i][1] - this.ylim.mean(), 2));
    }
    var maxVal = this.polygonsorted.max();
    var nearest = this.polygonsorted.indexOf(this.polygonsorted.min());
    var xbox = [0,0];
    var ybox = [0,0];
    var color;

    for(i=0;i<this.polygons.length;i++) {
      /* Calculate the intersecting box */
      if(this.xlim[0]>this.polygons[nearest][0].min()) xbox[0] = this.xpixel * Math.floor(this.xlim[0]/this.xpixel);
      else xbox[0] = this.xpixel * Math.floor(this.polygons[nearest][0].min()/this.xpixel);

      if(this.xlim[1]<this.polygons[nearest][0].max()) xbox[1] = this.xpixel * Math.ceil(this.xlim[1]/this.xpixel);
      else xbox[1] = this.xpixel * Math.ceil(this.polygons[nearest][0].max()/this.xpixel);

      if(this.ylim[0]>this.polygons[nearest][1].min()) ybox[0] = this.ypixel * Math.floor(this.ylim[0]/this.ypixel);
      else ybox[0] = this.ypixel * Math.floor(this.polygons[nearest][1].min()/this.ypixel);

      if(this.ylim[0]<this.polygons[nearest][1].max()) ybox[1] = this.ypixel * Math.ceil(this.ylim[1]/this.ypixel);
      else ybox[1] = this.ypixel * Math.ceil(this.polygons[nearest][1].max()/this.ypixel);

      for(j = xbox[0]; j <= xbox[1]; j += this.xpixel) {
        for(k = ybox[0]; k <= ybox[1]; k += this.ypixel) {
	  if(pip(this.polygons[nearest][0], this.polygons[nearest][1], j, k)) {
	    color = Math.round(24 * (this.model.pred(j, k) - this.model.response_min) / (this.model.response_range));
            if(color<0) color = 0;
            else if(color>24) color = 24;
            this.pixel(j, k, this.colorspectrum.terraincolors[color])
	  }
        }
      }

      this.polygonsorted[nearest] = maxVal;
      nearest = this.polygonsorted.indexOf(this.polygonsorted.min());

    }
  }

  return true;
}