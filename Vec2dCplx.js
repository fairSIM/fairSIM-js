/*
This file is part of Free Analysis and Interactive Reconstruction
for Structured Illumination Microscopy in JavaScript (fairSIM-js).

fairSIM-js is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

fairSIM-js is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with fairSIM-js.  If not, see <http://www.gnu.org/licenses/>
*/


function valOtf( dist ) {
    if ((dist<0)||(dist>=1))
	return 0.0;
    if (dist == 1 )
	return 1.0;
    return (2/Math.PI)*(Math.acos(dist) - dist*Math.sqrt(1-dist*dist));
}

function OtfVals() {
    this.objNA = 1.4;
    this.emLambda = 525;
    this.pxlSize = 0.08; 
    this.attFactor = 0.4;
}

class Vec2dCplx {

    constructor( size ) {
	this.size = size;
	this.length = size*size;
        this.data = new Float32Array(this.length*2);

        this.fftForward   = new FFT.complex(size, false);  
	this.fftBackwards = new FFT.complex(size, true);  
    }
   
    // duplicate this vector
    duplicate() {
	var ret = new Vec2dCplx( this.size );
	for ( var i=0; i<this.length*2; i++)
	    ret.data[i] = this.data[i];
	return ret;
    }

    // Fourier-transform this vector (forward or backwards)
    fft2d( dir)  {
	
	if (dir) {
	    // horz.
	    for ( var i = 0 ; i< this.size; i++ ) {
		this.fftForward.process( this.data, i*this.size, 1, this.data, i*this.size, 1, 'complex');
	    }
	    // vert.
	    for ( var i = 0 ; i< this.size; i++ ) {
		this.fftForward.process( this.data, i, this.size, this.data, i, this.size, 'complex');
	    }
	} else {
	    // horz.
	    for ( var i = 0 ; i< this.size; i++ ) {
		this.fftBackwards.process( this.data, i*this.size, 1, this.data, i*this.size, 1, 'complex');
	    }
	    for ( var i = 0 ; i< 2*this.length; i++) {
		this.data[i] /= 1.*this.length;
	    }
	    // vert.
	    for ( var i = 0 ; i< this.size; i++ ) {
		this.fftBackwards.process( this.data, i, this.size, this.data, i, this.size, 'complex');
	    }


	}

    }

    // set real-valued entries of this vector from 'input'
    // sets imaginary-valued entries to 0
    setFromTiff(input) {

	for ( var i=0 ; i<this.length ; i++ ) {
	    this.data [ i*2 ] = 1.*input[i];
	    this.data [ i*2+1 ] =  0.;
	}
    }   

    // set real-valued entries of this vector from 'input'
    // sets imaginary-valued entries to 0
    setFromRGB(input) {

	for ( var i=0 ; i<this.length ; i++ ) {
	    this.data [ i*2 ] = input[i*4];
	    this.data [ i*2+1 ] =  0.;
	}
    }   

    // set each imaginary component to zero
    clearImag() {
	for (var i=0; i<this.length; i++) {
	    this.data[2*i+1] = 0.0;
	}
    }

    // set entry at x,y to val, imag. component to zero
    set( x, y, val ) {
	this.data[ (this.size * y + x)*2 + 0 ] = val;
	this.data[ (this.size * y + x)*2 + 1 ] = 0;
    }

    // set entry at x,y, to (re,im)
    set( x, y, re, im ) {
	this.data[ (this.size * y + x)*2 + 0 ] = re;
	this.data[ (this.size * y + x)*2 + 1 ] = im;
    }


    // point-wise multiplication with other vector
    times( vec ) {
	if ( vec.length != this.length ) {
	    throw "times: vector size mismatch";
	}
	for ( var i=0 ; i<this.length ; i++ ) {
	    var r1 = this.data[ i*2 + 0 ];	    
	    var i1 = this.data[ i*2 + 1 ];	    
	    var r2 = vec.data[ i*2 + 0 ];	    
	    var i2 = vec.data[ i*2 + 1 ];	    

	    var rr = r1 * r2 - i1 * i2 ;
	    var ir = r1 * i2 + r2 * i1 ; 
	    
	    this.data[ i * 2 + 0 ]  = rr;
	    this.data[ i * 2 + 1 ]  = ir;
	}
    }

    // multiply vector with complex scalar
    mult( re, im ) {
	for ( var i=0 ; i<this.length ; i++ ) {
	    var r1 = this.data[ i*2 + 0 ];	    
	    var i1 = this.data[ i*2 + 1 ];	    
	    var rr = r1 * re - i1 * im ;
	    var ir = r1 * im + re * i1 ; 
	    this.data[ i * 2 + 0 ]  = rr;
	    this.data[ i * 2 + 1 ]  = ir;
	}
    }

    // add a Vector to this one
    add( vec ) {
	if ( vec.length != this.length ) {
	    throw "add: vector size mismatch";
	}
	for ( var i=0 ; i<this.length*2 ; i++ ) {
	    this.data[i]  += vec.data[i];
	}
    }

    // conjugate each entry
    conj() {
	for (var i=0; i<this.length; i++) {
	    this.data[2*i+1] *= -1.0;
	}
    }
    
    // compare two vectors, return sum of abs. point-wise difference
    comp( vec ) {
	if ( vec.length != this.length ) {
	    throw "add: vector size mismatch";
	}
	var res = 0.0;
	for ( var i=0 ; i<this.length ; i++ ) {
	    res += Math.abs( this.data[2*i] - vec.data[2*i] ); 
	    res += Math.abs( this.data[2*i+1] - vec.data[2*i+1] ); 
	}
	return res;
    }


    // find the maximum, excluding regions close to origin
    findMax( minDist = 100 ) {

	var max = 0, maxRe, maxIm;
	var posX = 0;
	var posY = 0; 

	for ( var y=0; y<this.size; y++ ) 
	for ( var x=0; x<this.size; x++ ) {
	    var xr = x;
	    var yr = y;
	    if (xr>this.size/2) xr -= this.size;	    
	    if (yr>this.size/2) yr -= this.size;
	    if ( Math.sqrt( xr*xr + yr*yr ) > minDist ) {
		var data = this.data;
		var i = x + this.size*y;
		var abs = Math.sqrt( data[2*i]*data[2*i] + data[2*i+1]*data[2*i+1]);
		if ( abs > max ) {
		    max = abs;
		    maxRe = data[2*i];
		    maxIm = data[2*i+1];
		    posX = x;
		    posY = y;
		}
	    }
	} 
    
	return [ posX, posY, max, maxRe, maxIm ];

    }

    // set every entry to zero
    zero() {
	for (var i=0; i<this.length*2; i++) {
	    this.data[i]=0.0;
	}
    }


    // apply a window function do make the FFT nicer
    applyWindow( pxl = 10 ) {

	const w = this.size;

	for ( var i=0 ; i< pxl ; i++ ) {
	    
	    const v = .5*(1-Math.cos(Math.PI*i/(pxl-1)));
	    
	    for ( var j=0 ; j< this.size ; j++ ) {
		this.data [ 2 * ( i + w*j ) + 0 ]  *= v;
		this.data [ 2 * ( i + w*j ) + 1 ]  *= v;
		this.data [ 2 * ( j + w*i ) + 0 ]  *= v;
		this.data [ 2 * ( j + w*i ) + 1 ]  *= v;
		this.data [ 2 * ( (w-1-i) + w*j ) + 0 ]  *= v;
		this.data [ 2 * ( (w-1-i) + w*j ) + 1 ]  *= v;
		this.data [ 2 * ( (w-1-j) + w*(w-i-1) ) + 0 ]  *= v;
		this.data [ 2 * ( (w-1-j) + w*(w-i-1) ) + 1 ]  *= v;
	    }
	}
    }
    
    // get minimum and maximum (abs)
    getAbsMinMax() {

	var min = Number.MAX_VALUE;
	var max = Number.MIN_VALUE;
	var data = this.data;

	for ( var i=0; i<this.length; i++) {
	    var abs = Math.sqrt( data[2*i]*data[2*i] + data[2*i+1]*data[2*i+1]);
	    if (abs > max ) max = abs;
	    if (abs < min ) min = abs;
	}
	
	return [min, max];

    }
    
    // get minimum and maximum only from real-value'd componenet
    getRealMinMax() {

	var min = Number.MAX_VALUE;
	var max = Number.MIN_VALUE;
	var data = this.data;

	for ( var i=0; i<this.length; i++) {
	    var abs = data[2*i];
	    if (abs > max ) max = abs;
	    if (abs < min ) min = abs;
	}
	
	return [min, max];

    }


    // cross-correlate with given vector and offset
    crossCorrelate( vec, kx, ky ) {

	var sumRe=0, sumIm=0;

	for ( var y = 0; y< this.size; y++)
	for ( var x = 0; x< this.size; x++) {

		var pha = 2 * Math.PI * (kx*x + ky*y) / this.size;    
		var co  = Math.cos( pha );
		var si  = Math.sin( pha );
	
		var i = this.size * y + x;	
		var r1 = vec.data[ i*2 + 0 ];	    
		var i1 = vec.data[ i*2 + 1 ];	    
		var rr = r1 * co - i1 * si ;
		var ir = r1 * si + co * i1 ; 

		var r2 = this.data[ i*2 + 0 ];
		var i2 = this.data[ i*2 + 1 ];

		sumRe += rr * r2 - ir * i2;
		sumIm += rr * i2 + r2 * ir;
	}
	
	return [ sumRe, sumIm, Math.sqrt( sumRe*sumRe + sumIm*sumIm ) ];

    }


    // swap the quandrands of this vector
    swapQuadrands() {

	var size= this.size;
	var tmpRe=0.0, tmpIm=0.0;
	var data = this.data;

	for ( var y=0 ; y<size/2 ; y++ ) {
	    for ( var x=0 ; x<size/2 ; x++ ) {
		// 1 <-> 3
		tmpRe = this.data[ 2* (x+size*y) + 0 ];
		tmpIm = this.data[ 2* (x+size*y) + 1 ];
		this.data[ 2* (x+size*y) + 0 ] = this.data[ 2 * (x+size/2 + (y+size/2)*size) + 0 ];
		this.data[ 2* (x+size*y) + 1 ] = this.data[ 2 * (x+size/2 + (y+size/2)*size) + 1 ];
		this.data[ 2 * (x+size/2 + (y+size/2)*size) + 0 ] = tmpRe;
		this.data[ 2 * (x+size/2 + (y+size/2)*size) + 1 ] = tmpIm; 
		// 2 <-> 4
		tmpRe = this.data[ 2* (x+size/2+ size*y) + 0 ];
		tmpIm = this.data[ 2* (x+size/2+ size*y) + 1 ];
		this.data[ 2* (x+size/2+ size*y) + 0 ] = this.data[ 2 * (x + (y+size/2)*size) + 0 ];
		this.data[ 2* (x+size/2+ size*y) + 1 ] = this.data[ 2 * (x + (y+size/2)*size) + 1 ];
		this.data[ 2 * (x+ (y+size/2)*size) + 0 ] = tmpRe;
		this.data[ 2 * (x+ (y+size/2)*size) + 1 ] = tmpIm; 
	    }
	}
    }






    // TODO: this should be the same as copyMirrored(), but on an existing
    // vector. However, it isn't, and I am too lazy to debug this now

    /*
    // mirror this vector
    // this keeps val[0,0] at ret[0,0], so suited for fft*(-k) <-> fft(k) - style
    // operations
    mirror() {

	var size= this.size;
	var tmpRe=0.0, tmpIm=0.0;
	var data = this.data;
		
	for (var c=1; c<this.size/2; c++) {
	    // lines
	    tmpRe = this.data[2*c+0];
	    tmpIm = this.data[2*c+1];
	    this.data[2*c+0] = this.data[2*(this.size-c)+0];
	    this.data[2*c+1] = this.data[2*(this.size-c)+1];
	    this.data[2*(this.size-c)+0] = tmpRe;
	    this.data[2*(this.size-c)+1] = tmpIm;
    	    
	    tmpRe = this.data[2*(c*this.size)+0]; 
	    tmpIm = this.data[2*(c*this.size)+1]; 
    	    this.data[2*(c*this.size)+0] = this.data[2*this.size*(this.size-c)+0];
	    this.data[2*(c*this.size)+1] = this.data[2*this.size*(this.size-c)+1];
	    this.data[2*this.size*(this.size-c)+0] = tmpRe;
	    this.data[2*this.size*(this.size-c)+1] = tmpIm; 
	}   

	for ( var y=1 ; y<size/2 ; y++ ) {
	    for ( var x=1 ; x<size/2 ; x++ ) {
		// 1 <-> 3
		tmpRe = this.data[ 2* (x+size*y) + 0 ];
		tmpIm = this.data[ 2* (x+size*y) + 1 ];
		this.data[ 2* (x+size*y) + 0 ] = this.data[ 2 * (size-x + (size-y)*size) + 0 ];
		this.data[ 2* (x+size*y) + 1 ] = this.data[ 2 * (size-x + (size-y)*size) + 1 ];
		this.data[ 2 * (size-x + (size-y)*size) + 0 ] = tmpRe;
		this.data[ 2 * (size-x + (size-y)*size) + 1 ] = tmpIm; 
		// 2 <-> 4
		tmpRe = this.data[ 2* (size-x + size*y) + 0 ];
		tmpIm = this.data[ 2* (size-x + size*y) + 1 ];
		this.data[ 2* (size-x + size*y) + 0 ] = this.data[ 2 * (x + (size-y)*size) + 0 ];
		this.data[ 2* (size-x + size*y) + 1 ] = this.data[ 2 * (x + (size-y)*size) + 1 ];
		this.data[ 2 * (x+ (size-y)*size) + 0 ] = tmpRe;
		this.data[ 2 * (x+ (size-y)*size) + 1 ] = tmpIm; 
	    }
	}

    }
    */
    

    // returns an x-y-mirror'd copy of this vector
    // this keeps val[0,0] at ret[0,0], so suited for fft*(-k) <-> fft(k) - style
    // operations
    duplicateMirrored() {
	var ret = new Vec2dCplx( this.size );
	for (var c=1; c<this.size; c++) {
	    // lines
	    ret.data[2*c+0] = this.data[2*(this.size-c)+0];
	    ret.data[2*c+1] = this.data[2*(this.size-c)+1];
	    ret.data[2*(c*this.size)+0] = this.data[2*this.size*(this.size-c)+0];
	    ret.data[2*(c*this.size)+1] = this.data[2*this.size*(this.size-c)+1];
	}

	for (var y=1; y<this.size; y++) {
	    // quadrands
	    for (var x=1; x<this.size; x++) {
		var i = x + y * this.size;
		var j = (this.size - x  )  + (this.size - y )*this.size;	    
		ret.data[2*j+0] = this.data[2*i+0];
		ret.data[2*j+1] = this.data[2*i+1];
	    }
	}
	
	// DC
	ret.data[0] = this.data[0] ;
	ret.data[1] = this.data[1] ;

	return ret;
    }




    // compute an (if requested, quadrand-swapped and/or log'd)
    // power spectrum
    getImg( swapQuadrand = true, compLog = true) {

	// create scaling
	var mm, min, max;
	mm = this.getAbsMinMax();
	
	if (compLog) {
	    //logger("log min/max "+mm[0]+" "+mm[1]);
	    min = Math.log(mm[0]);
	    max = Math.log(mm[1]);

	    if ( isNaN( min ) || max-min > 30) 
		min = max-30;
	} else {
	    min = mm[0];
	    max = mm[1];
	    //logger("min/max "+min+" "+max);
	}
	
	//logger( min+" "+max);

	var ret = new Float32Array( this.size * this.size );
	var data = this.data;

	var size= this.size;
	for ( var y=0 ; y<size ; y++ ) {
	for ( var x=0 ; x<size ; x++ ) {

	    var i = y*size + x;
	    var r = Math.sqrt( data[i*2] * data[i*2] + data[i*2+1] * data[i*2+1]);
	
	    if ( compLog ) {
		r = (Math.log(r) - min)/(max-min);
		if (isNaN(r) || r<0) r=0;
	    } else {
		//r = (r-min)/(max-min);
	    }
	    
	    if (swapQuadrand) {
		var xo = (x<size/2)?(x+size/2):(x-size/2);
		var yo = (y<size/2)?(y+size/2):(y-size/2);
		var io = xo + yo*size;
		ret[io] = r;
	    } else {
		ret[i] = r;
	    }
	}}
    
	return ret;

    }

    // multiply vector w. fourier shift phases
    fourierShift( kx, ky ) {
	
	for (var y=0; y<this.size; y++) {
	    for (var x=0; x<this.size; x++) {
	    
		var pha = 2 * Math.PI * (kx*x + ky*y) / this.size;    
		var co  = Math.cos( pha );
		var si  = Math.sin( pha );
	
		var i = this.size * y + x;	
		var r1 = this.data[ i*2 + 0 ];	    
		var i1 = this.data[ i*2 + 1 ];	    
		var rr = r1 * co - i1 * si ;
		var ir = r1 * si + co * i1 ; 
		this.data[ i * 2 + 0 ]  = rr;
		this.data[ i * 2 + 1 ]  = ir;

	    }
	}
    }

    // multiply each entry with a phase
    multPhase( pha ) {
	var si = Math.sin(pha);
	var co = Math.cos(pha);

	//logger("multiplying phase: "+pha);

	for (var i=0; i<this.length; i++) {
	    var r1 = this.data[ i*2 + 0 ];	    
	    var i1 = this.data[ i*2 + 1 ];	    
	    var rr = r1 * co - i1 * si ;
	    var ir = r1 * si + co * i1 ; 
	    this.data[ i * 2 + 0 ]  = rr;
	    this.data[ i * 2 + 1 ]  = ir;
	}

    }



    // add the content of a vector half in size, at position x,y
    paste( vec, kx, ky , invert =false) {
	if ( vec.size > this.size ) {
	    throw "vector size mismatch";
	}

	const hsi = vec.size/2;

	for (var y=0; y<vec.size; y++)
	for (var x=0; x<vec.size; x++) {
	   
	    var xo = (x<hsi)?(x):(x+this.size-vec.size);	
	    var yo = (y<hsi)?(y):(y+this.size-vec.size);	

	    xo = (xo+kx)%this.size;
	    yo = (yo+ky)%this.size;


	    var io = xo + yo * this.size;
	    var ii ;
	    if (!invert) {
		ii = x + y*vec.size;
	    } else {
		ii = (vec.size-x-1) + vec.size*(vec.size-y-1);
	    }

	    this.data[ io*2 + 0 ]  += vec.data[ ii*2 + 0 ];
	    this.data[ io*2 + 1 ]  += vec.data[ ii*2 + 1 ];

	}
    }

    // create a simple, 2D OTF
    createOtf( otfData , kx=0, ky=0, att=-1, coShift=1 ) {

	const cyclPxl   =  1./(this.size*otfData.pxlSize);
	const cutoff    =  ((2*otfData.objNA)/(otfData.emLambda/1000.))*coShift;
	const cutoffPxl =  cutoff/cyclPxl;

	for (var y=0; y<this.size; y++) {
	    for (var x=0; x<this.size; x++) {

		var xi = x + kx;
		var yi = y + ky;

		var xh = ((xi<this.size/2)?( xi):(xi-this.size)) ;
		var yh = ((yi<this.size/2)?( yi):(yi-this.size)) ;

		var dist = Math.sqrt( xh*xh+yh*yh );
		var val  = valOtf( dist/cutoffPxl );

		if ( att>0 && (dist/cutoffPxl) < 1 ) {
		    val *= 1.0-0.99*Math.exp( -(dist/cutoffPxl/att));
		}

		this.data[(x+y*this.size)*2+0] = val;
		this.data[(x+y*this.size)*2+1] = 0.0;

	    }
	}

    }


    // cut out regions beyond OTF support
    maskOtf( otfData , coShift = 1.) {

	const cyclPxl   =  1./(imageSize*pxlSize);
	const cutoff    =  ((2*otfData.objNA)/(otfData.emLambda/1000.))*coShift;
	const cutoffPxl =  cutoff/cyclPxl;

	for (var y=0; y<vec.size; y++) {
	    for (var x=0; x<vec.size; x++) {

		var xh = ((x<this.size/2)?( x):(x-this.size)) ;
		var yh = ((y<this.size/2)?( y):(y-this.size)) ;
		var dist = Math.sqrt( xh*xh+yh*yh );

		if ( (dist/cutoffPxl) >1.1) {
		    this.data[(x+y*this.size)*2+0] = 0.0;
		    this.data[(x+y*this.size)*2+1] = 0.0;
		} else if ( (dist/cutoffPxl) >1 ) {
		    var d2 = ((dist/cutoffPxl) -1.)*10.;
		    var v  = .5*(1+Math.cos(Math.PI*d2));
		    this.data[(x+y*this.size)*2+0] *= v;
		    this.data[(x+y*this.size)*2+1] *= v;
		}
	    }
	}

    }



}




