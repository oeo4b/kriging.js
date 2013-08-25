/**
 *
 * kriging.js
 *
 * Copyright 2012-2013
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

/* Non-R function -- essential for R_solve_ */
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

/* Matrix inversion -- Gauss-jordan elimination */
var R_solve = function(a) {
  var n = a.length;
  var m = n;
  var b = new Array(n);
  var indxc = new Array(n);
  var indxr = new Array(n);
  var ipiv = new Array(n);

  var i, icol, irow, j, k, l, ll;
  var big, dum, pivinv, temp;

  for(i=0;i<n;i++) {
    b[i] = new Array(n);
    for(j=0;j<n;j++) {
      if(i==j) b[i][j] = 1;
      else b[i][j] = 0;
    }
  }
  for(j=0;j<n;j++) ipiv[j] = 0;
  for(i=0;i<n;i++) {
    big = 0;
    for(j=0;j<n;j++) {
      if(ipiv[j]!=1) {
	for(k=0;k<n;k++) {
	  if(ipiv[k]==0) {
	    if(Math.abs(a[j][k])>=big) {
	      big = Math.abs(a[j][k]);
              irow = j;
              icol = k;
            }
          }
        }
      }
    }
    ++(ipiv[icol]);

    if(irow!=icol) {
      for(l=0;l<n;l++) {
        temp = a[irow][l];
        a[irow][l] = a[icol][l];
        a[icol][l] = temp;
      }
      for(l=0;l<m;l++) {
        temp = b[irow][l];
        b[irow][l] = b[icol][l];
        b[icol][l] = temp;
      }
    }

    indxr[i] = irow;
    indxc[i] = icol;

    if(a[icol][icol]==0) { /* Singular matrix */
      return false;
    }

    pivinv = 1 / a[icol][icol];
    a[icol][icol] = 1;
    for(l=0;l<n;l++) a[icol][l] *= pivinv;
    for(l=0;l<m;l++) b[icol][l] *= pivinv;

    for(ll=0;ll<n;ll++) {
      if(ll!=icol) {
        dum = a[ll][icol];
        a[ll][icol] = 0;
        for(l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
        for(l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
      }
    }
  }

  for(l=(n-1);l>=0;l--) {
    if(indxr[l]!=indxc[l]) {
      for(k=0;k<n;k++) {
        temp = a[k][indxr[l]];
        a[k][indxr[l]] = a[k][indxc[l]];
        a[k][indxc[l]] = temp;
      }
    }
  }

  return a;
}

var R_solve_cramers_rule = function(x) {
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
  var pixelsize = 1;
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
    zoomin.value = "zoom in";
    zoomout.value = "zoom out";
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
  this.canvas.colorspectrum.heatcolors = ["#FF0000", "#FF0700", "#FF0E00", "#FF1500", "#FF1C00", "#FF2200", "#FF2900", "#FF3000", "#FF3700", "#FF3E00", "#FF4500", "#FF4C00", "#FF5300", "#FF5A00", "#FF6000", "#FF6700", "#FF6E00", "#FF7500", "#FF7C00", "#FF8300", "#FF8A00", "#FF9100", "#FF9800", "#FF9F00", "#FFA500", "#FFAC00", "#FFB300", "#FFBA00", "#FFC100", "#FFC800", "#FFCF00", "#FFD600", "#FFDD00", "#FFE300", "#FFEA00", "#FFF100", "#FFF800", "#FFFF00", "#FFFF0B", "#FFFF20", "#FFFF35", "#FFFF4A", "#FFFF60", "#FFFF75", "#FFFF8A", "#FFFF9F", "#FFFFB5", "#FFFFCA", "#FFFFDF", "#FFFFF4"];
  this.canvas.colorspectrum.terraincolors = ["#00A600", "#07A800", "#0EAB00", "#16AE00", "#1DB000", "#25B300", "#2DB600", "#36B800", "#3EBB00", "#47BE00", "#50C000", "#59C300", "#63C600", "#6CC800", "#76CB00", "#80CE00", "#8BD000", "#95D300", "#A0D600", "#ABD800", "#B6DB00", "#C2DE00", "#CEE000", "#D9E300", "#E6E600", "#E6DD09", "#E7D612", "#E7CF1C", "#E8C825", "#E8C32E", "#E9BE38", "#E9BA41", "#EAB74B", "#EAB454", "#EBB25E", "#EBB167", "#ECB171", "#ECB17B", "#EDB285", "#EDB48E", "#EEB798", "#EEBAA2", "#EFBFAC", "#EFC4B6", "#F0C9C0", "#F0D0CA", "#F1D7D4", "#F1DFDE", "#F2E8E8", "#F2F2F2"];
  this.canvas.colorspectrum.topocolors = ["#4C00FF", "#3B00FF", "#2800FF", "#1600FF", "#0400FF", "#000DFF", "#001FFF", "#0032FF", "#0043FF", "#0055FF", "#0068FF", "#007AFF", "#008BFF", "#009EFF", "#00AFFF", "#00C1FF", "#00D3FF", "#00E5FF", "#00FF4D", "#00FF38", "#00FF24", "#00FF0F", "#05FF00", "#1AFF00", "#2EFF00", "#42FF00", "#57FF00", "#6BFF00", "#80FF00", "#94FF00", "#A8FF00", "#BDFF00", "#D1FF00", "#E6FF00", "#FFFF00", "#FFF90C", "#FFF318", "#FFED24", "#FFE930", "#FFE53B", "#FFE247", "#FFDF53", "#FFDD5F", "#FFDC6B", "#FFDB77", "#FFDB83", "#FFDB8F", "#FFDC9B", "#FFDEA7", "#FFE0B3"];
  this.canvas.colorspectrum.cmcolors = ["#80FFFF", "#85FFFF", "#8AFFFF", "#8FFFFF", "#94FFFF", "#99FFFF", "#9EFFFF", "#A3FFFF", "#A8FFFF", "#ADFFFF", "#B3FFFF", "#B8FFFF", "#BDFFFF", "#C2FFFF", "#C7FFFF", "#CCFFFF", "#D1FFFF", "#D6FFFF", "#DBFFFF", "#E0FFFF", "#E6FFFF", "#EBFFFF", "#F0FFFF", "#F5FFFF", "#FAFFFF", "#FFFAFF", "#FFF5FF", "#FFF0FF", "#FFEBFF", "#FFE6FF", "#FFE0FF", "#FFDBFF", "#FFD6FF", "#FFD1FF", "#FFCCFF", "#FFC7FF", "#FFC2FF", "#FFBDFF", "#FFB8FF", "#FFB3FF", "#FFADFF", "#FFA8FF", "#FFA3FF", "#FF9EFF", "#FF99FF", "#FF94FF", "#FF8FFF", "#FF8AFF", "#FF85FF", "#FF80FF"];

  this.canvas.render = function() {
    this.clear();
    this.background();
    this.points();
  }  

  this.canvas.pixel = function(x, y, col) {
    this.ctx.fillStyle = col;
    
    /* Spaced-out pixels */
    //this.ctx.fillRect((x-this.xlim[0])/this.xratio - pixelsize/2 + 1, this.height - (y-this.ylim[0])/this.yratio - pixelsize/2 + 1, pixelsize - 2, pixelsize - 2);

    /* Solid map */
    this.ctx.fillRect((x-this.xlim[0])/this.xratio - pixelsize/2, this.height - (y-this.ylim[0])/this.yratio - pixelsize/2, pixelsize, pixelsize);
  }

  this.canvas.focus = function(x, y, col) {
    this.ctx.beginPath();
    this.ctx.arc((x-this.xlim[0])/this.xratio - (2*pixelsize)/2, this.height - (y-this.ylim[0])/this.yratio - (2*pixelsize)/2, 2*pixelsize, 0, 2 * Math.PI, false);
    this.ctx.fillStyle = col;
    this.ctx.fill();
  }
  
  this.canvas.clear = function() {
    this.ctx.clearRect(0, 0, this.width, this.height);
  }

  /* Plot observed points */
  this.canvas.points = function() {
    var i;
    for(i=0;i<this.model.n;i++) {
      this.focus(this.model.x[i], this.model.y[i], "black");
    }
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
	    color = Math.round(49 * (this.model.pred(j, k) - this.model.response_min) / (this.model.response_range));
            if(color<0) color = 0;
            else if(color>49) color = 49;
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