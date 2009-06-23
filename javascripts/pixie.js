/*global jQuery */
(function($) {
  //
  // PNG drawing library for JavaScript.
  // Copyright (C) 1999 by Roger E Critchlow Jr,
  // Santa Fe, New Mexico, USA.
  //
  // Licensed under the Academic Free License version 2.1
  //
  // The home page for Pnglets is http://www.elf.org/pnglets,
  // a copy of the AFL may be found at http://www.opensource.org/licenses/afl-2.1.php,
  // Pnglets were inspired by and copied from gd1.3, http://www.boutell.com/gd,
  // other parts were inspired by or copied from Tcl/Tk, http://www.scriptics.com,
  // and some algorithms were taken from Foley & van Dam 2nd Edition.
  //
  // Thanks to Alex Vincent for pointing out the advantages of eliminating strict
  // javascript warnings.
  //

  // create a new Pnglet of specified width, height, and depth
  // width and height are specified in pixels
  // depth is really the number of palette entries
  function Pnglet(width,height,depth) {
    this.width = width || 16;
    this.height = height || 16;
    this.depth = Math.min(256, depth || 16);

    // pixel data and row filter identifier size
    this.pix_size = height*(width+1);

    // deflate header, pix_size, block headers, adler32 checksum
    this.data_size = 2 + this.pix_size + 5*Math.floor((this.pix_size+0xffff-1)/0xffff) + 4;

    // offsets and sizes of Png chunks
    this.ihdr_offs = 0;					// IHDR offset and size
    this.ihdr_size = 4+4+13+4;
    this.plte_offs = this.ihdr_offs+this.ihdr_size;	// PLTE offset and size
    this.plte_size = 4+4+3*depth+4;
    this.trns_offs = this.plte_offs+this.plte_size;	// tRNS offset and size
    this.trns_size = 4+4+depth+4;
    this.idat_offs = this.trns_offs+this.trns_size;	// IDAT offset and size
    this.idat_size = 4+4+this.data_size+4;
    this.iend_offs = this.idat_offs+this.idat_size;	// IEND offset and size
    this.iend_size = 4+4+4;
    this.png_size = this.iend_offs+this.iend_size;	// total PNG size

    // array of one byte strings
    this.png = new Array(this.png_size);

    // functions for initializing data
    function initialize(png, offs, str) {
      for (var i = 1; i < arguments.length; i += 1)
        if (typeof arguments[i].length != "undefined")
          for (var j = 0; j < arguments[i].length; j += 1)
            png[offs++] = arguments[i].charAt(j);
    };
    function byte2(w) { return String.fromCharCode((w>>8)&255, w&255); };
    function byte4(w) { return String.fromCharCode((w>>24)&255, (w>>16)&255, (w>>8)&255, w&255); };
    function byte2lsb(w) { return String.fromCharCode(w&255, (w>>8)&255); };

    // initialize everything to zero byte
    for (var i = 0; i < this.png_size; i += 1)
      this.png[i] = String.fromCharCode(0);

    // initialize non-zero elements
    initialize(this.png, this.ihdr_offs, byte4(this.ihdr_size-12), 'IHDR',
               byte4(width), byte4(height), String.fromCharCode(8, 3));
    initialize(this.png, this.plte_offs, byte4(this.plte_size-12), 'PLTE');
    initialize(this.png, this.trns_offs, byte4(this.trns_size-12), 'tRNS');
    initialize(this.png, this.idat_offs, byte4(this.idat_size-12), 'IDAT');
    initialize(this.png, this.iend_offs, byte4(this.iend_size-12), 'IEND');

    // initialize deflate header
    var header = ((8 + (7<<4)) << 8) | (3 << 6);
    header += 31 - (header % 31);
    initialize(this.png, this.idat_offs+8, byte2(header));

    // initialize deflate block headers
    for (i = 0; i*0xffff < this.pix_size; i += 1) {
    var size, bits;
    if (i + 0xffff < this.pix_size) {
        size = 0xffff;
        bits = String.fromCharCode(0);
    } else {
        size = this.pix_size - i*0xffff;
        bits = String.fromCharCode(1);
    }
    initialize(this.png, this.idat_offs+8+2+i*(5+0xffff), bits, byte2lsb(size), byte2lsb(~size));
    }

    // initialize palette hash
    this.palette = new Object();
    this.pindex = 0;
  }

  // version string/number
  Pnglet.version = "19990427.0";

  // test if coordinates are within bounds
  Pnglet.prototype.inBounds = function(x,y) { return x >= 0 && x < this.width && y >= 0 && y < this.height; }

  // clip an x value to the window width
  Pnglet.prototype.clipX = function(x) { return (x < 0) ? 0 : (x >= this.width) ? this.width-1 : x ; }

  // clip a y value to the window height
  Pnglet.prototype.clipY = function(y) { return (y < 0) ? 0 : (y >= this.height) ? this.height-1 : y ; }

  // compute the index into a png for a given pixel
  Pnglet.prototype.index = function(x,y) {
    var i = y*(this.width+1)+x+1;
    var j = this.idat_offs+8+2+Math.floor((i/0xffff)+1)*5+i;
    return j;
  }

  // make a color in a Pnglet
  Pnglet.prototype.color = function(red, green, blue, alpha) {
    alpha = alpha >= 0 ? alpha : 255;
    var rgba = (((((alpha<<8)+red)<<8)+green)<<8)+blue;
    if ( typeof this.palette[rgba] == "undefined") {
    if (this.pindex == this.depth) return String.fromCharCode(0);
    this.palette[rgba] = String.fromCharCode(this.pindex);
    this.png[this.plte_offs+8+this.pindex*3+0] = String.fromCharCode(red);
    this.png[this.plte_offs+8+this.pindex*3+1] = String.fromCharCode(green);
    this.png[this.plte_offs+8+this.pindex*3+2] = String.fromCharCode(blue);
    this.png[this.trns_offs+8+this.pindex] = String.fromCharCode(alpha);
    this.pindex += 1;
    }
    return this.palette[rgba];
  }

  // return true if this is a color
  Pnglet.prototype.isColor = function(color) {
    return typeof(color) == 'string' &&
    color.length == 1 &&
    color.charCodeAt(0) >= 0 &&
    color.charCodeAt(0) < this.depth;
  }

  // find the red, green, blue, or alpha value of a Pnglet color
  Pnglet.prototype.red = function(color) { return this.png[this.plte_offs+8+color.charCodeAt(0)*3+0].charCodeAt(0); }
  Pnglet.prototype.green = function(color) { return this.png[this.plte_offs+8+color.charCodeAt(0)*3+1].charCodeAt(0); }
  Pnglet.prototype.blue = function(color) { return this.png[this.plte_offs+8+color.charCodeAt(0)*3+2].charCodeAt(0); }
  Pnglet.prototype.alpha = function(color) { return this.png[this.trns_offs+8+color.charCodeAt(0)].charCodeAt(0); }

  // draw a point or points
  Pnglet.prototype.point = function(pointColor, x0, y0) {
    var a = arguments;
    this.pointNXY(pointColor, (a.length-1)/2, function(i) { return a[2*i+1]; }, function(i) { return a[2*i+2]; });
  }

  Pnglet.prototype.pointNXY = function(pointColor, n, x, y) {
    if ( ! this.isColor(pointColor))
    return;
    for (var i = 0; i < n; i += 1) {
    var x1 = x(i), y1 = y(i);
    if (this.inBounds(x1,y1))
        this.png[this.index(x1,y1)] = pointColor;
    }
  }

  // read a pixel
  Pnglet.prototype.getPoint = function(x,y) { return this.inBounds(x,y) ? this.png[this.index(x,y)] : String.fromCharCode(0); }

  // draw a horizontal line
  Pnglet.prototype.horizontalLine = function(lineColor, x1, x2, y) {
    if ( ! this.isColor(lineColor))
    return;
    x1 = this.clipX(x1);
    x2 = this.clipX(x2);
    var x;
    if (x1 < x2)
    for (x = x1; x <= x2; x += 1)
        this.png[this.index(x,y)] = lineColor;
    else
    for (x = x2; x <= x1; x += 1)
        this.png[this.index(x,y)] = lineColor;
  }

  // draw a vertical line
  Pnglet.prototype.verticalLine = function(lineColor, x, y1, y2) {
    if ( ! this.isColor(lineColor))
    return;
    y1 = this.clipY(y1);
    y2 = this.clipY(y2);
    var y;
    if (y1 < y2)
    for (y = y1; y <= y2; y += 1)
        this.png[this.index(x,y)] = lineColor;
    else
    for (y = y2; y <= y1; y += 1)
        this.png[this.index(x,y)] = lineColor;
  }

  // draw a general line
  Pnglet.prototype.generalLine = function(lineColor, x1, y1, x2, y2) {
    if ( ! this.isColor(lineColor))
    return;
    var dx = Math.abs(x2-x1), dy = Math.abs(y2-y1);
    var incr1, incr2, d, x, y, xend, yend, xdirflag, ydirflag, xinc, yinc;
    if (dy <= dx) {
    d = 2*dy - dx;
    incr1 = 2*dy;
    incr2 = 2 * (dy - dx);
    if (x1 > x2) {
        x = x2;
        y = y2;
        ydirflag = -1;
        xend = x1;
    } else {
        x = x1;
        y = y1;
        ydirflag = 1;
        xend = x2;
    }
    yinc = (((y2 - y1) * ydirflag) > 0) ? 1 : -1;
    this.point(lineColor, x, y);
    while (x++ < xend) {
        if (d < 0) {
      d += incr1;
        } else {
      y += yinc;
      d += incr2;
        }
        this.point(lineColor, x, y);
    }
    } else {			/* dy > dx */
    d = 2*dx - dy;
    incr1 = 2*dx;
    incr2 = 2 * (dx - dy);
    if (y1 > y2) {
        y = y2;
        x = x2;
        yend = y1;
        xdirflag = -1;
    } else {
        y = y1;
        x = x1;
        yend = y2;
        xdirflag = 1;
    }
    xinc = (((x2 - x1) * xdirflag) > 0) ? 1 : -1;
    this.point(lineColor, x, y);
    while (y++ < yend) {
        if (d < 0) {
      d += incr1;
        } else {
      x += xinc;
      d += incr2;
        }
        this.point(lineColor, x, y);
    }
    }
  }

  // draw a line
  Pnglet.prototype.line = function(lineColor, x0, y0) {
    var a = arguments;
    this.lineNXY(lineColor, (a.length-1)/2, function(i) { return a[2*i+1]; }, function(i) { return a[2*i+2]; });
  }

  Pnglet.prototype.lineNXY = function(lineColor, n, x, y) {
    if ( ! this.isColor(lineColor))
    return;
    var x1 = x(0), y1 = y(0);
    for (var i = 1; i < n; i += 1) {
    var x2 = x(i), y2 = y(i);
    if (x1 == x2)
        this.verticalLine(lineColor, x1, y1, y2);
    else if (y1 == y2)
        this.horizontalLine(lineColor, x1, x2, y1);
    else
        this.generalLine(lineColor, x1, y1, x2, y2);
    x1 = x2;
    y1 = y2;
    }
  }

  // draw a polygon
  Pnglet.prototype.polygon = function(outlineColor, fillColor, x1, y1) {
    var a = arguments;
    this.polygonNXY(outlineColor, fillColor, (a.length-2)/2, function(i) {return a[2*i+2];}, function(i) {return a[2*i+3];});
  }

  Pnglet.prototype.polygonNXY = function(outlineColor, fillColor, n, x, y) {
    if (n <= 0)
    return;
    if (this.isColor(fillColor))
    this.concaveNXY(fillColor, n, x, y);
    if (this.isColor(outlineColor))
    this.lineNXY(outlineColor, n+1, function(i) { return x(i%n); }, function(i) { return y(i%n); });
  }

  /*
   * Concave Polygon Scan Conversion
   * by Paul Heckbert
   * from "Graphics Gems", Academic Press, 1990
   */
  Pnglet.prototype.concaveNXY = function(fillColor, n, ptx, pty) {
    function Edge(ex, edx, ei) {	/* a polygon edge */
    this.x = ex;	/* x coordinate of edge's intersection with current scanline */
    this.dx = edx;	/* change in x with respect to y */
    this.i = ei;	/* edge number: edge i goes from pt[i] to pt[i+1] */
    };
    function cdelete(di) {	/* remove edge i from active list */
      for (var j = 0; j < active.length; j += 1)
        if (active[j].i == di)
          active.splice(j, 1);
    };
    function cinsert(ii, iy) {	/* append edge i to end of active list */
      var ij = ii<n-1 ? ii+1 : 0;
      var px, py, qx, qy;
      if (pty(ii) < pty(ij)) {
        px = ptx(ii); py = pty(ii);
        qx = ptx(ij); qy = pty(ij);
      } else {
        px = ptx(ij); py = pty(ij);
        qx = ptx(ii); qy = pty(ii);
      }
      /* initialize x position at intersection of edge with scanline y */
      var dx = (qx-px)/(qy-py);
      active.push(new Edge(dx*(iy+.5-py)+px, dx, ii));
    };

    var ind = new Array(n);		/* list of vertex indices, sorted by pt[ind[j]].y */
    var active = new Array(0);		/* start with an empty active list */

    /* create y-sorted array of indices ind[k] into vertex list */
    for (var k = 0; k < n; k += 1) ind[k] = k;
    ind.sort(function(i1, i2) { return pty(i1) <= pty(i2) ? -1 : 1; });
    k = 0;                        /* ind[k] is next vertex to process */
    var y0 = Math.max(0, Math.ceil(pty(ind[0])+.5));			/* ymin of polygon */
    var y1 = Math.min(this.height, Math.floor(pty(ind[n-1])-.5));	/* ymax of polygon */

    for (var y = y0; y <= y1; y += 1) {		/* step through scanlines */
    /* scanline y is at y+.5 in continuous coordinates */

    /* check vertices between previous scanline and current one, if any */
    for (; k<n && pty(ind[k]) <= y+.5; k += 1) {
        /* to simplify, if pt.y=y+.5, pretend it's above */
        /* invariant: y-.5 < pt[i].y <= y+.5 */
        var i = ind[k];

        /*
         * insert or delete edges before and after vertex i (i-1 to i,
         * and i to i+1) from active list if they cross scanline y
         */
        var j = (i-1+n)%n;		/* vertex previous to i */
        if (pty(j) <= y-.5)	{	/* old edge, remove from active list */
      cdelete(j);
        } else if (pty(j) > y+.5) {	/* new edge, add to active list */
      cinsert(j, y);
        }
        if (i != ind[k]) {
      alert("Your browser's implementation of JavaScript is seriously broken,\n"+
            "as in variables are changing value of their own volition.\n"+
            "You should upgrade to a newer version browser.");
      return;
        }
        j = (i+1)%n;		/* vertex next after i */
        if (pty(j) <= y-.5)	{	/* old edge, remove from active list */
      cdelete(i);
        } else if (pty(j) > y+.5) {	/* new edge, add to active list */
      cinsert(i, y);
        }
    }

    /* sort active edge list by active[j].x */
    active.sort(function(u,v) { return u.x <= v.x ? -1 : 1; });

    /* draw horizontal segments for scanline y */
    for (j = 0; j < active.length; j += 2) {	/* draw horizontal segments */
        /* span 'tween j & j+1 is inside, span tween j+1 & j+2 is outside */
        var xl = Math.ceil(active[j].x+.5);		/* left end of span */
        if (xl<0) xl = 0;
        var xr = Math.floor(active[j+1].x-.5);	/* right end of span */
        if (xr>this.width-1) xr = this.width-1;
        if (xl<=xr)
      this.horizontalLine(fillColor, xl, xr, y);	/* draw pixels in span */
        active[j].x += active[j].dx;	/* increment edge coords */
        active[j+1].x += active[j+1].dx;
    }
    }
  }

  // draw a rectangle
  Pnglet.prototype.rectangle = function(outlineColor, fillColor, x0,y0,x1,y1) {
    if (this.isColor(fillColor))
    for (var y = y0; y < y1; y += 1)
        this.horizontalLine(fillColor, x0+1, x1-2, y);
    if (this.isColor(outlineColor)) {
    this.horizontalLine(outlineColor, x0, x1-1, y0);
    this.horizontalLine(outlineColor, x0, x1-1, y1-1);
    this.verticalLine(outlineColor, x0, y0, y1-1);
    this.verticalLine(outlineColor, x1-1, y0, y1-1);
    }
  }

  // draw an arc
  Pnglet.prototype.arc = function(outlineColor, cx,cy,w,h, s,e) {
    var p = this.midpointEllipse(cx,cy, w,h, s,e);
    function x(i) { return p[i*2]; };
    function y(i) { return p[i*2+1]; };
    this.lineNXY(outlineColor, p.length/2, x, y);
  }

  // draw an oval
  Pnglet.prototype.oval = function(outlineColor, fillColor, cx,cy,w,h) {
    var p = this.midpointEllipse(cx,cy, w,h, 0,359);
    function x(i) { return p[i*2]; };
    function y(i) { return p[i*2+1]; };
    this.polygonNXY(outlineColor, fillColor, p.length/2, x, y);
  }

  // draw an arc with chord
  Pnglet.prototype.chord = function(outlineColor, fillColor, cx,cy,w,h, s,e) {
    var p = this.midpointEllipse(cx,cy, w,h, s,e);
    function x(i) { return p[i*2]; };
    function y(i) { return p[i*2+1]; };
    this.polygonNXY(outlineColor, fillColor, p.length/2, x, y);
  }

  // draw an arc with pieslice
  Pnglet.prototype.pieslice = function(outlineColor, fillColor, cx,cy,w,h, s,e) {
    var p = this.midpointEllipse(cx,cy, w,h, s,e);
    p[p.length] = cx;
    p[p.length] = cy;
    function x(i) { return p[i*2]; };
    function y(i) { return p[i*2+1]; };
    this.polygonNXY(outlineColor, fillColor, p.length/2, x, y);
  }

  // oval arcs
  // generate points of oval circumference
  // midpoint ellipse, Foley & van Dam, 2nd Edition, p. 90, 1990
  Pnglet.prototype.midpointEllipse = function(cx,cy, w,h, s,e) {
    var a = Math.floor(w/2), b = Math.floor(h/2);
    var a2 = a*a, b2 = b*b, x = 0, y = b;
    var d1 = b2 - a2*b + a2/4;
    cx = Math.floor(cx);
    cy = Math.floor(cy);
    var p = new Array();

    // quadrant I, anticlockwise
    p.push(x,-y);
    while (a2*(y-1/2) > b2*(x+1)) {
    if (d1 < 0) {
        d1 += b2*(2*x+3);
    } else {
        d1 += b2*(2*x+3) + a2*(-2*y + 2);
        y -= 1;
    }
    x += 1;
    p.unshift(x,-y);
    }
    var d2 = b2*(x+1/2)*(x+1/2) + a2*(y-1)*(y-1) - a2*b2;
    while (y > 0) {
    if (d2 < 0) {
        d2 += b2*(2*x+2) + a2*(-2*y+3);
        x += 1;
    } else {
        d2 += a2*(-2*y+3);
    }
    y -= 1;
    p.unshift(x,-y);
    }
    // quadrant II, anticlockwise
    var n4 = p.length;
    for (var i = n4-4; i >= 0; i -= 2)
      p.push(-p[i], p[i+1]);
    // quadrants III and IV, anticlockwise
    var n2 = p.length;
    for (i = n2-4; i > 0; i -= 2)
      p.push(p[i], -p[i+1]);

    // compute start and end indexes from start and extent
    e %= 360;
    if (e < 0) {
    s += e;
    e = -e;
    }
    s %= 360;
    if (s < 0)
    s += 360;
    var is = Math.floor(s/359 * p.length/2);
    var ie = Math.floor(e/359 * p.length/2)+1;
    p = p.slice(is*2).concat(p.slice(0, is*2)).slice(0, ie*2);

    // displace to center
    for (i = 0; i < p.length; i += 2) {
    p[i] += cx;
    p[i+1] += cy;
    }
    return p;
  }

  // fill a region
  // from gd1.3 with modifications
  Pnglet.prototype.fill = function(outlineColor,fillColor,x,y) {
    if (outlineColor) {			// fill to outline color
    /* Seek left */
    var leftLimit = -1;
    for (var i = x; i >= 0 && this.getPoint(i, y) != outlineColor; i -= 1)
        leftLimit = i;

    if (leftLimit == -1)
        return;

    /* Seek right */
    var rightLimit = x;
    for (i = (x+1); i < this.width && this.getPoint(i, y) != outlineColor; i += 1)
        rightLimit = i;

    /* fill extent found */
    this.horizontalLine(fillColor, leftLimit, rightLimit, y);

    /* Seek above and below */
    for (var dy = -1; dy <= 1; dy += 2) {
        if (this.inBounds(x,y+dy)) {
      var lastBorder = 1;
      for (i = leftLimit; i <= rightLimit; i++) {
            var c = this.getPoint(i, y+dy);
            if (lastBorder) {
        if ((c != outlineColor) && (c != fillColor)) {
                this.fill(outlineColor, fillColor, i, y+dy);
                lastBorder = 0;
        }
            } else if ((c == outlineColor) || (c == fillColor)) {
        lastBorder = 1;
            }
      }
        }
    }

    } else {			// flood fill color at x, y
    /* Test for completion */
    var oldColor = this.getPoint(x, y);
    if (oldColor == fillColor)
        return;

    /* Seek left */
    leftLimit = (-1);
    for (i = x; i >= 0 && this.getPoint(i, y) == oldColor; i--)
        leftLimit = i;

    if (leftLimit == -1)
        return;

    /* Seek right */
    rightLimit = x;
    for (i = (x+1); i < this.width && this.getPoint(i, y) == oldColor; i++)
        rightLimit = i;

    /* Fill extent found */
    this.horizontalLine(fillColor, leftLimit, rightLimit, y);

    /* Seek above and below */
    for (dy = -1; dy <= 1; dy += 2) {
        if (this.inBounds(x,y+dy)) {
      lastBorder = 1;
      for (i = leftLimit; i <= rightLimit; i++) {
            c = this.getPoint(i, y+dy);
            if (lastBorder) {
        if (c == oldColor) {
                this.fill(null, fillColor, i, y+dy);
                lastBorder = 0;
        }
            } else if (c != oldColor) {
        lastBorder = 1;
            }
      }
        }
    }
    }
  }

  // smoothed points
  Pnglet.prototype.smoothPoint = function(smoothSteps, pointColor, x0, y0) {
    var a = arguments, self = this, n = (a.length-2)/2;
    this.smooth(smoothSteps,
                function(n, x, y) { self.pointNXY(pointColor, n, x, y); },
                n,
                function(i) { return a[2*i+2]; },
                function(i) { return a[2*i+3]; });
  }

  // smoothed polyline
  Pnglet.prototype.smoothLine = function(smoothSteps, lineColor, x0, y0) {
    var a = arguments, self = this, n = (a.length-2)/2;
    this.smooth(smoothSteps,
                function(n, x, y) { self.lineNXY(lineColor, n, x, y); },
                n,
                function(i) { return a[2*i+2]; },
                function(i) { return a[2*i+3]; });
  }

  // smoothed polygon
  Pnglet.prototype.smoothPolygon = function(smoothSteps, outlineColor, fillColor, x0, y0) {
    var a = arguments, self = this, n = (a.length-3)/2 + 1;
    this.smooth(smoothSteps,
                function(n, x, y) { self.polygonNXY(outlineColor, fillColor, n, x, y); },
                n,
                function(i) { return a[2*(i%(n-1))+3]; },
                function(i) { return a[2*(i%(n-1))+4]; });
  }

  // generate smoothSteps points for the line segment connecting
  // each consecutive pair of points in x(i), y(i).
  // adapted from the source for tk8.1b3, http://www.scriptics.com
  Pnglet.prototype.smooth = function(smoothSteps, fNXY, n, x, y) {
    var control = new Array(8);
    var outputPoints = 0;
    var dblPoints = new Array();

    // compute numSteps of smoothed points
    // according to the basis in control[]
    // placing points into coordPtr[coordOff]
    function smoothPoints(control, numSteps, coordPtr, coordOff) {
    for (var i = 1; i <= numSteps; i++, coordOff += 2) {
        var t = i/numSteps, t2 = t*t, t3 = t2*t,
      u = 1.0 - t, u2 = u*u, u3 = u2*u;
        coordPtr[coordOff+0] = control[0]*u3 + 3.0 * (control[2]*t*u2 + control[4]*t2*u) + control[6]*t3;
        coordPtr[coordOff+1] = control[1]*u3 + 3.0 * (control[3]*t*u2 + control[5]*t2*u) + control[7]*t3;
    }
    };

    /*
     * If the curve is a closed one then generate a special spline
     * that spans the last points and the first ones.  Otherwise
     * just put the first point into the output.
     */

    var closed = (x(0) == x(n-1)) && (y(0) == y(n-1));
    if (closed) {
    control[0] = 0.5*x(n-2) + 0.5*x(0);
    control[1] = 0.5*y(n-2) + 0.5*y(0);
    control[2] = 0.167*x(n-2) + 0.833*x(0);
    control[3] = 0.167*y(n-2) + 0.833*y(0);
    control[4] = 0.833*x(0) + 0.167*x(1);
    control[5] = 0.833*y(0) + 0.167*y(1);
    control[6] = 0.5*x(0) + 0.5*x(1);
    control[7] = 0.5*y(0) + 0.5*y(1);
    dblPoints[2*outputPoints+0] = control[0];
    dblPoints[2*outputPoints+1] = control[1];
    outputPoints += 1;
    smoothPoints(control, smoothSteps, dblPoints, 2*outputPoints);
    outputPoints += smoothSteps;
    } else {
    dblPoints[2*outputPoints+0] = x(0);
    dblPoints[2*outputPoints+1] = y(0);
    outputPoints += 1;
    }

    for (var i = 2; i < n; i += 1) {
    var j = i - 2;
    /*
     * Set up the first two control points.  This is done
     * differently for the first spline of an open curve
     * than for other cases.
     */
    if ((i == 2) && !closed) {
        control[0] = x(j);
        control[1] = y(j);
        control[2] = 0.333*x(j) + 0.667*x(j+1);
        control[3] = 0.333*y(j) + 0.667*y(j+1);
    } else {
        control[0] = 0.5*x(j) + 0.5*x(j+1);
        control[1] = 0.5*y(j) + 0.5*y(j+1);
        control[2] = 0.167*x(j) + 0.833*x(j+1);
        control[3] = 0.167*y(j) + 0.833*y(j+1);
    }

    /*
     * Set up the last two control points.  This is done
     * differently for the last spline of an open curve
     * than for other cases.
     */

    if ((i == (n-1)) && !closed) {
        control[4] = .667*x(j+1) + .333*x(j+2);
        control[5] = .667*y(j+1) + .333*y(j+2);
        control[6] = x(j+2);
        control[7] = y(j+2);
    } else {
        control[4] = .833*x(j+1) + .167*x(j+2);
        control[5] = .833*y(j+1) + .167*y(j+2);
        control[6] = 0.5*x(j+1) + 0.5*x(j+2);
        control[7] = 0.5*y(j+1) + 0.5*y(j+2);
    }

    /*
     * If the first two points coincide, or if the last
     * two points coincide, then generate a single
     * straight-line segment by outputting the last control
     * point.
     */

    if (((x(j) == x(j+1)) && (y(j) == y(j+1)))
          || ((x(j+1) == x(j+2)) && (y(j+1) == y(j+2)))) {
      dblPoints[2*outputPoints+0] = control[6];
      dblPoints[2*outputPoints+1] = control[7];
      outputPoints += 1;
      continue;
    }

    /*
     * Generate a Bezier spline using the control points.
     */
    smoothPoints(control, smoothSteps, dblPoints, 2*outputPoints);
    outputPoints += smoothSteps;
    }

    // draw the points
    // anonymous functions don't work here
    // they result in "undefined" point values
    function myx(i) { return Math.round(dblPoints[2*i]); }
      function myy(i) { return Math.round(dblPoints[2*i+1]); }
      fNXY(outputPoints, myx, myy);
  }

  // output a PNG string
  Pnglet.prototype.output = function() {
    // output translations
    function initialize(png, offs, str) {
    for (var i = 1; i < arguments.length; i += 1)
        if (typeof arguments[i].length != "undefined")
          for (var j = 0; j < arguments[i].length; j += 1)
            png[offs++] = arguments[i].charAt(j);
    }
      function byte4(w) { return String.fromCharCode((w>>24)&255, (w>>16)&255, (w>>8)&255, w&255); }

      // compute adler32 of output pixels + row filter bytes
      var BASE = 65521; /* largest prime smaller than 65536 */
    var NMAX = 5552;  /* NMAX is the largest n such that 255n(n+1)/2 + (n+1)(BASE-1) <= 2^32-1 */
    var s1 = 1;
    var s2 = 0;
    var n = NMAX;
    for (var y = 0; y < this.height; y += 1)
    for (var x = -1; x < this.width; x += 1) {
        s1 += this.png[this.index(x,y)].charCodeAt(0);
        s2 += s1;
        if ((n -= 1) == 0) {
      s1 %= BASE;
      s2 %= BASE;
      n = NMAX;
        }
    }
    s1 %= BASE;
    s2 %= BASE;
    initialize(this.png, this.idat_offs+this.idat_size-8, byte4((s2 << 16) | s1));

    // compute crc32 of the PNG chunks
    function crc32(png, offs, size) {
    var crc = -1;		// initialize crc
    for (var i = 4; i < size-4; i += 1)
        crc = Pnglet.crc32_table[(crc ^ png[offs+i].charCodeAt(0)) & 0xff] ^ ((crc >> 8) & 0x00ffffff);
    initialize(png, offs+size-4, byte4(crc ^ -1));
    }

      crc32(this.png, this.ihdr_offs, this.ihdr_size);
    crc32(this.png, this.plte_offs, this.plte_size);
    crc32(this.png, this.trns_offs, this.trns_size);
    crc32(this.png, this.idat_offs, this.idat_size);
    crc32(this.png, this.iend_offs, this.iend_size);

    // convert PNG to string
    return "\211PNG\r\n\032\n"+this.png.join('');
  }

  /* Table of CRCs of all 8-bit messages. */
  Pnglet.crc32_table = new Array(256);
  for (var n = 0; n < 256; n++) {
    var c = n;
    for (var k = 0; k < 8; k++) {
      if (c & 1)
        c = -306674912 ^ ((c >> 1) & 0x7fffffff);
      else
        c = (c >> 1) & 0x7fffffff;
    }
    Pnglet.crc32_table[n] = c;
  }

  // Global lookup arrays for base64 conversions
  var enc64List, dec64List;
  // Load the lookup arrays once
  (function() {
      enc64List = new Array();
      dec64List = new Array();
      var i;
      for (i = 0; i < 26; i++) {
          enc64List[enc64List.length] = String.fromCharCode(65 + i);
      }
      for (i = 0; i < 26; i++) {
          enc64List[enc64List.length] = String.fromCharCode(97 + i);
      }
      for (i = 0; i < 10; i++) {
          enc64List[enc64List.length] = String.fromCharCode(48 + i);
      }
      enc64List[enc64List.length] = "+";
      enc64List[enc64List.length] = "/";
      for (i = 0; i < 128; i++) {
          dec64List[dec64List.length] = -1;
      }
      for (i = 0; i < 64; i++) {
          dec64List[enc64List[i].charCodeAt(0)] = i;
      }
  })();

  var base64Encode = function(str) {
      var c, d, e, end = 0;
      var u, v, w, x;
      var ptr = -1;
      var input = str.split("");
      var output = "";
      while(end == 0) {
          c = (typeof input[++ptr] != "undefined") ? input[ptr].charCodeAt(0) :
              ((end = 1) ? 0 : 0);
          d = (typeof input[++ptr] != "undefined") ? input[ptr].charCodeAt(0) :
              ((end += 1) ? 0 : 0);
          e = (typeof input[++ptr] != "undefined") ? input[ptr].charCodeAt(0) :
              ((end += 1) ? 0 : 0);
          u = enc64List[c >> 2];
          v = enc64List[(0x00000003 & c) << 4 | d >> 4];
          w = enc64List[(0x0000000F & d) << 2 | e >> 6];
          x = enc64List[e & 0x0000003F];

          // handle padding to even out unevenly divisible string lengths
          if (end >= 1) {x = "=";}
          if (end == 2) {w = "=";}

          if (end < 3) {output += u + v + w + x;}
      }
      return output;
  }

  var base64Decode = function(str) {
      var c=0, d=0, e=0, f=0, i=0, n=0;
      var input = str.split("");
      var output = "";
      var ptr = 0;
      do {
          f = input[ptr++].charCodeAt(0);
          i = dec64List[f];
          if ( f >= 0 && f < 128 && i != -1 ) {
              if ( n % 4 == 0 ) {
                  c = i << 2;
              } else if ( n % 4 == 1 ) {
                  c = c | ( i >> 4 );
                  d = ( i & 0x0000000F ) << 4;
              } else if ( n % 4 == 2 ) {
                  d = d | ( i >> 2 );
                  e = ( i & 0x00000003 ) << 6;
              } else {
                  e = e | i;
              }
              n++;
              if ( n % 4 == 0 ) {
                  output += String.fromCharCode(c) +
                            String.fromCharCode(d) +
                            String.fromCharCode(e);
              }
          }
      }
      while (typeof input[ptr] != "undefined");
      output += (n % 4 == 3) ? String.fromCharCode(c) + String.fromCharCode(d) :
                ((n % 4 == 2) ? String.fromCharCode(c) : "");
      return output;
  }

  /************************************************************************
   *
   * pixie.js - The jQuery pixel editor.
   *
   ************************************************************************/

  var imageDir = "../images/pixie/";

  var actions = {
    undo: {
      name: "Undo",
      undoable: false,
      hotkeys: ["ctrl+z"],
      perform: function(canvas) {
        canvas.undo();
      }
    },

    redo: {
      name: "Redo",
      undoable: false,
      hotkeys: ["ctrl+y"],
      perform: function(canvas) {
        canvas.redo();
      }
    },

    clear: {
      name: "Clear Layer",
      perform: function(canvas) {
        canvas.eachPixel(function(pixel) {
          pixel.color(null);
        });
      }
    },

    preview: {
      name: "Preview",
      perform: function(canvas) {
        canvas.showPreview();
      }
    },

    save: {
      name: "Download Image",
      hotkeys: ["ctrl+s"],
      perform: function(canvas) {
        document.location.href = 'data:image/octet-stream;base64,' + canvas.toBase64();
      }
    },

    left: {
      name: "Move Left",
      menu: false,
      hotkeys: ['left'],
      perform: function(canvas) {
        canvas.eachPixel(function(pixel, x, y) {
          var rightPixel = canvas.getPixel(x + 1, y);

          if(rightPixel) {
            pixel.color(rightPixel.color());
          } else {
            pixel.color(null);
          }
        });
      }
    },

    right: {
      name: "Move Right",
      menu: false,
      hotkeys: ['right'],
      perform: function(canvas) {
        var width = canvas.width;
        var height = canvas.height;
        for(var x = width-1; x >= 0; x--) {
          for(var y = 0; y < height; y++) {
            var currentPixel = canvas.getPixel(x, y);
            var leftPixel = canvas.getPixel(x-1, y);

            if(leftPixel) {
              currentPixel.color(leftPixel.color());
            } else {
              currentPixel.color(null);
            }
          }
        }
      }
    },

    up: {
      name: "Move Up",
      menu: false,
      hotkeys: ['up'],
      perform: function(canvas) {
        canvas.eachPixel(function(pixel, x, y) {
          var lowerPixel = canvas.getPixel(x, y + 1);

          if(lowerPixel) {
            pixel.color(lowerPixel.color());
          } else {
            pixel.color(null);
          }
        });
      }
    },

    down: {
      name: "Move Down",
      menu: false,
      hotkeys: ['down'],
      perform: function(canvas) {
        var width = canvas.width;
        var height = canvas.height;
        for(var x = 0; x < width; x++) {
          for(var y = height-1; y >= 0; y--) {
            var currentPixel = canvas.getPixel(x, y);
            var upperPixel = canvas.getPixel(x, y-1);

            if(upperPixel) {
              currentPixel.color(upperPixel.color());
            } else {
              currentPixel.color(null);
            }
          }
        }
      }
    }
  };

  var CloneTool = function() {
    var cloneX, cloneY, targetX, targetY;

    return {
      name: "Clone",
      hotkeys: ['C'],
      icon: imageDir + "clone.png",
      cursor: "url("+ imageDir +"clone.png) 0 0, default",
      mousedown: function(e) {
        if(e.shiftKey) {
          cloneX = this.x;
          cloneY = this.y;
        } else {
          targetX = this.x;
          targetY = this.y;
          var selection = this.canvas.getPixel(cloneX, cloneY);

          if(selection) {
            this.color(selection.color());
          }
        }
      },
      mouseenter: function(e) {
        var deltaX = this.x - targetX;
        var deltaY = this.y - targetY;
        var selection = this.canvas.getPixel(cloneX + deltaX, cloneY + deltaY);

        if(selection) {
          this.color(selection.color());
        }
      }
    };
  };

  var tools = {
    pencil: {
      name: "Pencil",
      hotkeys: ['P'],
      icon: imageDir + "pencil.png",
      cursor: "url(" + imageDir + "pencil.png) 4 14, default",
      mousedown: function(e, color) {
        this.color(color);
      },
      mouseenter: function(e, color) {
        this.color(color);
      }
    },
    
    brush: {
      name: "Brush",
      hotkeys: ['B'],
      icon: imageDir + "paintbrush.png",
      cursor: "url(" + imageDir + "paintbrush.png) 4 14, default",
      mousedown: function(e, color) {
        this.color(color);

        $.each(this.canvas.getNeighbors(this.x, this.y), function(i, neighbor) {
          if(neighbor) {
            neighbor.color(color);
          }
        });
      },
      mouseenter: function(e, color) {
        this.color(color);

        $.each(this.canvas.getNeighbors(this.x, this.y), function(i, neighbor) {
          if(neighbor) {
            neighbor.color(color);
          }
        });
      }
    },

    dropper: {
      name: "Dropper",
      hotkeys: ['I'],
      icon: imageDir + "dropper.png",
      cursor: "url(" + imageDir + "dropper.png) 13 13, default",
      mousedown: function() {
        this.canvas.color(this.color());
        this.canvas.setTool(tools.pencil);
      }
    },

    eraser: {
      name: "Eraser",
      hotkeys: ['E'],
      icon: imageDir + "eraser.png",
      cursor: "url(" + imageDir + "eraser.png) 4 11, default",
      mousedown: function() {
        this.color(null);
      },
      mouseenter: function() {
        this.color(null);
      }
    },

    fill: {
      name: "Fill",
      hotkeys: ['F'],
      icon: imageDir + "fill.png",
      cursor: "url(" + imageDir + "fill.png) 12 13, default",
      mousedown: function(e, newColor, pixel) {
        // Store original pixel's color here
        var originalColor = this.color();

        // Return if original color is same as currentColor
        if(newColor === originalColor) {
          return;
        }

        var q = new Array();
        pixel.color(newColor);
        q.push(pixel);

        while(q.length > 0) {
          pixel = q.pop();

          // Add neighboring pixels to the queue
          var neighbors = this.canvas.getNeighbors(pixel.x, pixel.y);

          $.each(neighbors, function(index, neighbor) {
            if(neighbor && neighbor.css("backgroundColor") === originalColor) {
              neighbor.color(newColor);
              q.push(neighbor);
            }
          });
        }
      }
    },

    clone: CloneTool()
  };

  var falseFn = function() {return false};
  var div = '<div></div>';
  var clear = '<div class="clear"></div>';
  var ColorPicker = function() {
    return $('<input></input>').addClass('color').colorPicker();
  };

  var rgbParser = /^rgb\((\d{1,3}),\s*(\d{1,3}),\s*(\d{1,3})\)$/;

  var UndoStack = function() {
    var undos = [];
    var redos = [];
    var empty = true;

    return {
      popUndo: function() {
        var undo = undos.pop();

        if(undo) {
          redos.push(undo);
        }

        return undo;
      },

      popRedo: function() {
        var undo = redos.pop();

        if(undo) {
          undos.push(undo);
        }

        return undo;
      },

      next: function() {
        var last = undos[undos.length - 1];
        if(!last || !empty) {
          undos.push({});
          empty = true;
          // New future, no more redos
          redos = [];
        }
      },

      add: function(object, data) {
        var last = undos[undos.length - 1];

        // Only store this objects data if it is not already present.
        if(!last[object]) {
          last[object] = data;
          empty = false;
        }

        return this;
      }
    };
  };
  
  $.fn.pixie = function(options) {
    
    options = options || {};
    var width = options.width || 16;
    var height = options.height || 16;
    var pixelWidth = 16;
    var pixelHeight = 16;
    var initializer = options.initializer;
    var layers = options.layers || 2;

    return this.each(function() {
      var pixie = $(div).addClass('pixie');
      var actionsMenu = $(div).addClass('actions');
      var canvas = $(div).addClass('canvas').css({width: pixelWidth*width, height: pixelHeight*height});
      var toolbar = $(div).addClass('toolbar');
      var colorbar = $(div).addClass('toolbar');
      var preview = $(div).addClass('preview').css({width: width, height: height});
      var layerMenu = $(div).addClass('actions');

      var undoStack = UndoStack();

      var currentTool = undefined;
      var active = false;
      var mode = undefined;
      var primaryColorPicker = ColorPicker();
      var secondaryColorPicker = ColorPicker();

      colorbar.append(
        primaryColorPicker
      ).append(
        secondaryColorPicker
      );

      pixie
        .bind('contextmenu', falseFn)
        .bind("mousedown", function(e){
          var target = $(e.target);
          
          if(target.is(".swatch")) {
            canvas.color(target.css('backgroundColor'), e.button !== 0);
          }
        })
        .bind("mouseup", function(e) {
          active = false;
          mode = undefined;
        });

      var pixels = [];

      for(var layer = 0; layer < layers; layer++) {
        var layerDiv = $(div).addClass("layer").css("zIndex", layer);
        pixels[layer] = [];

        (function(currentLayer) {
          var layerSelection = $("<a href='#' title='Layer "+ currentLayer +"'>"+ currentLayer +"</a>")
            .addClass('tool')
            .bind("mousedown", function(e) {
              layer = currentLayer;
              layerMenu.children().removeClass("active");
              $(this).addClass("active");
            })
            .click(falseFn)

          if(currentLayer === 0) {
            layerSelection.addClass("active");
          }

          layerMenu.append(layerSelection);
        })(layer);

        for(var row = 0; row < height; row++) {
          pixels[layer][row] = [];

          for(var col = 0; col < width; col++) {
            var pixel = $(div).addClass('pixel');
            pixels[layer][row][col] = pixel;

            $.extend(pixel, {
              x: col,
              y: row,
              z: layer,
              canvas: canvas,
              toString: function() {
                return "[Pixel: " + this.x + ", " + this.y + ", " + this.z + "]";
              },
              color: function(color) {
                if(arguments.length >= 1) {
                  undoStack.add(this, {pixel: this, color: this.css("backgroundColor")});
                  this.css("backgroundColor", color);
                  return this;
                } else {
                  return this.css("backgroundColor");
                }
              }
            });


            // Only the top layer should be sensitive to events
            if(layer == layers - 1) {
              (function(pixel, row, col){
                pixel
                  .bind("mousedown", function(e){
                    undoStack.next();
                    active = true;
                    if(e.button === 0) {
                      mode = "P";
                    } else {
                      mode = "S";
                    }

                    e.preventDefault();
                  })
                  .bind("mousedown mouseup mouseenter", function(e) {
                    var p = pixels[layer][row][col];
                    if(active && currentTool && currentTool[e.type]) {
                      currentTool[e.type].call(p, e, canvas.color(), p);
                    }
                  });
              })(pixel, row, col);
            }

            layerDiv.append(pixel);

          }

          layerDiv.append(clear);
        }

        canvas.append(layerDiv);
      }

      canvas.append(clear);

      layer = 0;

      $.extend(canvas, {
        eachPixel: function(fn, z) {
          if(z === undefined) {
            z = layer;
          }

          for(row = 0; row < height; row++) {
            for(col = 0; col < width; col++) {
              var pixel = pixels[z][row][col];
              fn.call(pixel, pixel, col, row);
            }
          }

          return canvas;
        },

        getPixel: function(x, y, z) {
          if(z === undefined) {
            z = layer;
          }

          if(y >= 0 && y < height) {
            if(x >= 0 && x < width) {
              return pixels[z][y][x];
            }
          }

          return undefined;
        },
        
        getNeighbors: function(x, y, z) {
          if(z === undefined) {
            z = layer;
          }

          return [
            this.getPixel(x+1, y, z),
            this.getPixel(x, y+1, z),
            this.getPixel(x-1, y, z),
            this.getPixel(x, y-1, z)
          ];
        },
        
        toHex: function(bits) {
          var s = parseInt(bits).toString(16);
          if(s.length == 1) {
            s = '0' + s
          }
          return s;
        },

        parseColor: function(colorString) {
          if(!colorString || colorString == "transparent") {
            return false;
          }

          var bits = rgbParser.exec(colorString);
          return [
            this.toHex(bits[1]),
            this.toHex(bits[2]),
            this.toHex(bits[3])
          ].join('').toUpperCase();
        },

        color: function(color, alternate) {
          // Handle cases where nothing, or only true or false is passed in
          // i.e. when getting the alternate color `canvas.color(true)`
          if(arguments.length === 0 || color === false) {
            return mode == "S" ? 
              secondaryColorPicker.css('backgroundColor') :
              primaryColorPicker.css('backgroundColor');
          } else if(color === true) {
            // Switch color choice when alterate is true
            return mode == "S" ?
              primaryColorPicker.css('backgroundColor') :
              secondaryColorPicker.css('backgroundColor');
          }

          var parsedColor;
          if(color[0] != "#") {
            parsedColor = "#" + (this.parseColor(color) || "FFFFFF");
          } else {
            parsedColor = color;
          }

          if((mode == "S") ^ alternate) {
            secondaryColorPicker.val(parsedColor);
            secondaryColorPicker[0].onblur();
          } else {
            primaryColorPicker.val(parsedColor);
            primaryColorPicker[0].onblur();
          }

          return this;
        },

        addSwatch: function(color) {
          colorbar.append(
            $(div)
              .addClass('swatch')
              .css({backgroundColor: color})
          );
        },
        
        addAction: function(action) {
          var titleText = action.name;
          var undoable = action.undoable;

          function doIt() {
            if(undoable !== false) {
              undoStack.next();
            }
            action.perform(canvas);
          }

          if(action.hotkeys) {
            titleText += " ("+ action.hotkeys +")";

            $.each(action.hotkeys, function(i, hotkey) {
              $(document).bind('keydown', hotkey, function(e) {
                doIt();
                e.preventDefault();
              });
            });
          }

          if(action.menu !== false) {
            actionsMenu.append(
              $("<a href='#' title='"+ titleText +"'>"+ action.name +"</a>")
                .addClass('tool')
                .bind("mousedown", function(e) {
                  doIt();
                })
                .click(falseFn)
            );
          }
        },

        addTool: function(tool) {
          var alt = tool.name;

          function setMe() {
            canvas.setTool(tool);
            toolbar.children().removeClass("active");
            toolDiv.addClass("active");
          }

          if(tool.hotkeys) {
            alt += " ("+ tool.hotkeys +")";

            $.each(tool.hotkeys, function(i, hotkey) {
              $(document).bind('keydown', hotkey, function(e) {
                setMe();
                e.preventDefault();
              });
            });
          }

          var toolDiv = $("<img src='"+ tool.icon +"' alt='"+ alt +"' title='"+ alt +"'></img>")
            .addClass('tool')
            .bind('mousedown', function(e) {
              setMe();
            });

          toolbar.append(toolDiv);
        },

        setTool: function(tool) {
          currentTool = tool;
          if(tool.cursor) {
            pixie.css({cursor: tool.cursor});
          } else {
            pixie.css({cursor: "pointer"});
          }
        },

        toPNG: function() {
          var p = new Pnglet(width, height, 256);

          for(var z = 0; z < layers; z++) {
            this.eachPixel(function(pixel, x, y) {
              if(pixel.css('backgroundColor') == "transparent") {
                if(z === 0) {
                  p.point(p.color(0, 0, 0, 0), x, y);
                }
              } else {
                var rgb = rgbParser.exec(pixel.css('backgroundColor'));
                var alpha = 255;
                p.point(p.color(rgb[1], rgb[2], rgb[3], alpha), x, y);
              }
            }, z);
          }

          return p.output();
        },

        toBase64: function() {
          return base64Encode(this.toPNG());
        },

        toDataURL: function() {
          return 'url(data:image/png;base64,' + this.toBase64() + ')';
        },

        showPreview: function() {
          preview.css('backgroundImage', this.toDataURL());
        },

        undo: function() {
          var data = undoStack.popUndo();
          var swap;

          if(data) {
            $.each(data, function() {
              swap = this.color;
              this.color = this.pixel.css('backgroundColor');
              this.pixel.css('backgroundColor', (swap));
            });
          }
        },

        redo: function() {
          var data = undoStack.popRedo();
          var swap;

          if(data) {
            $.each(data, function() {
              swap = this.color;
              this.color = this.pixel.css('backgroundColor');
              this.pixel.css('backgroundColor', (swap));
            });
          }
        },

        width: width,
        height: height
      });

      $.each(actions, function(key, action) {
        canvas.addAction(action);
      });

      $.each(["#000", "#FFF", "#666", "#CCC", "#800", "#080", "#008", "#880", "#808", "#088"], function(i, color) {
        canvas.addSwatch(color);
      });

      $.each(tools, function(key, tool) {
        canvas.addTool(tool);
      });

      canvas.setTool(tools.pencil);
      toolbar.children().eq(0).addClass("active");

      if(initializer) {
        initializer(canvas);
      }

      pixie
        .append(actionsMenu)
        .append(toolbar)
        .append(canvas)
        .append(colorbar)
        .append(preview)
        .append(clear)
        .append(layerMenu);

      $(this).append(pixie);
    });
  };
})(jQuery);