/*
 Definition of map class

 */

Array.prototype.max = function() {
  return Math.max.apply(null, this)
}

Array.prototype.min = function() {
  return Math.min.apply(null, this)
}

/* Constructor */
function map(id, polygons, width, xlim, ylim) {
  /* Output testing */
  var o = document.getElementById("output");

  /* Global vars */
  var canvaspad = 50;
  var pixelsize = 4;

  /* Map switches */

  /* Canvas element */
  var canvasobj = document.getElementById(id);
  this.canvas = document.getElementById(id);
  this.canvas.ctx = this.canvas.getContext("2d");

  /* Generic methods */
  this.canvas.pip = function(X, Y, x, y) {
    var i, j;
    var c = false;

    for(i=0, j=X.length-1; i<X.length; j = i++) {
      if( ((Y[i]>y) != (Y[j]>y)) && (x<(X[j]-X[i]) * (y-Y[i]) / (Y[j]-Y[i]) + X[i]) ) {
	c = !c;
      }
    }
    return c;
  }

  /* Set up the canvas frame */
  this.canvas.width = width;
  this.canvas.height = window.innerHeight - this.canvas.offsetTop - canvaspad;
  this.canvas.xlim = xlim;
  this.canvas.ylim = ylim;
  this.canvas.xratio = (this.canvas.xlim[1]-this.canvas.xlim[0]) / this.canvas.width;
  this.canvas.yratio = (this.canvas.ylim[1]-this.canvas.ylim[0]) / this.canvas.height;
  this.canvas.xpixel = pixelsize * this.canvas.xratio;
  this.canvas.ypixel = pixelsize * this.canvas.yratio;
  this.canvas.style.border = "1px solid black";

  /* Mouse event properties */
  this.canvas.mousex = 0;
  this.canvas.mousey = 0;
  this.canvas.startx = 0;
  this.canvas.starty = 0;
  this.canvas.mousedown = false;

  /* Resize event handlers */
  window.onresize = function(e) {
    canvasobj.height = window.innerHeight - canvasobj.offsetTop - canvaspad;
    canvasobj.render();
  }

  /* Mouse event handlers */
  document.addEventListener("mousemove", function(e) {
    /* Reset mouse coordinates */
    canvasobj.mousex = e.pageX - canvasobj.offsetLeft;
    canvasobj.mousey = e.pageY - canvasobj.offsetTop;

   /* Drag the map if mouse is clicked */
    if(canvasobj.mousedown) {
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

  document.onselectstart = function() {
    if(canvasobj.mousedown) {
	this.style.cursor = "move";
	o.innerHTML = 1;
      return false;
    }
    else {
	this.style.cursor = "default";
	o.innerHTML = 2;
      return true;
    }
  }


  /* Rendering methods */
  this.canvas.render = function() {
    this.clear();
    this.background();
  }  

  this.canvas.pixel = function(x, y, col) {
    this.ctx.fillStyle = col;
    this.ctx.fillRect((x-this.xlim[0])/this.xratio - pixelsize/2, this.height - (y-this.ylim[0])/this.yratio - pixelsize/2, pixelsize, pixelsize);
  }
  
  this.canvas.clear = function() {
    this.ctx.clearRect(0, 0, this.width, this.height);
  }

  /* Fills the background with the polygons */
  this.canvas.background = function() {
    var i, j, k;

    for(i=0;i<polygons.length;i++) {
      for(j = this.xpixel * Math.floor(polygons[i][0].min()/this.xpixel); j <= (this.xpixel * Math.ceil(polygons[i][0].max()/this.xpixel)); j+=this.xpixel) {
	for(k = this.ypixel * Math.floor(polygons[i][1].min()/this.ypixel); k <= (this.ypixel * Math.ceil(polygons[i][1].max()/this.ypixel)); k+=this.ypixel) { 
	  if(this.pip(polygons[i][0], polygons[i][1], j, k)) {
	    this.pixel(j, k, "black");
	  }
	}
      }   
    }
  }

  this.render = function() {
    this.canvas.render();
  }

  return true;
}
