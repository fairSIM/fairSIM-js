
class Vec2dCplx {

    constructor( size ) {
	this.size = size;
	this.length = size*size;
        this.data = new Float32Array(this.length*2);

        this.fftForward   = new FFT.complex(size, false);  
	this.fftBackwards = new FFT.complex(size, true);  
    }
   
    copy() {
	var ret = new Vec2dCplx( this.size );
	for ( var i=0; i<this.length*2; i++)
	    ret.data[i] = this.data[i];
	return ret;
    }

 
    fft2d( dir = true )  {
	
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
	    // vert.
	    for ( var i = 0 ; i< this.size; i++ ) {
		this.fftBackwards.process( this.data, i, this.size, this.data, i, this.size, 'complex');
	    }
	    for ( var i = 0 ; i< 2*this.length; i++) {
		this.data[i] /= this.length;
	    }


	}

    }

    // set real-valued entries of this vector from 'input'
    // sets imaginary-valued entries to 0
    set(input) {

	for ( var i=0 ; i<this.length ; i++ ) {
	    this.data [ i*2 ] = input[i];
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


    set( x, y, val ) {
	this.data[ (this.size * y + x)*2 + 0 ] = val;
	this.data[ (this.size * y + x)*2 + 1 ] = 0;
    }

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


    // find the maximum, excluding regions close to origin
    findMax( minDist = 100 ) {

	var max = 0;
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
		    posX = x;
		    posY = y;
		}
	    }
	} 
    
	return [ posX, posY, max ];

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
    getMinMax() {

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


    // compute a (scaled) power spectrum
    getImg( swapQuadrand = true, compLog = true) {

	// create scaling
	var mm, min, max;
	mm = this.getMinMax();
	
	if (compLog) {
	    logger("log min/max "+mm[0]+" "+mm[1]);
	    min = Math.log(mm[0]);
	    max = Math.log(mm[1]);

	    if ( isNaN( min ) || max-min > 30) 
		min = max-30;
	} else {
	    min = mm[0];
	    max = mm[1];
	    logger("min/max "+min+" "+max);
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

    // add the content of a vector half in size, at position x,y
    paste( vec, kx, ky ) {
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
	    var ii = x + y*vec.size;

	    this.data[ io*2 + 0 ]  += vec.data[ ii*2 + 0 ];
	    this.data[ io*2 + 1 ]  += vec.data[ ii*2 + 1 ];

	}
    }

}




