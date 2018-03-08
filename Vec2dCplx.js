
class Vec2dCplx {

    constructor( size ) {
	this.size = size;
	this.length = size*size;
        this.data = new Float32Array(this.length*2);

        this.fftForward   = new FFT.complex(size, false);  
	this.fftBackwards = new FFT.complex(size, true);  
    }
    
    fft2d( dir )  {
	
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

    // compute a (scaled) power spectrum
    pwSpec( swapQuadrand = true) {

	// create scaling
	var mm = this.getMinMax();
	var min = Math.log(mm[0]);
	var max = Math.log(mm[1]);

	if ( isNaN( min ) || max-min > 30) 
	    min = max-30;
	
	logger( min+" "+max);

	var ret = new Float32Array( this.size * this.size );
	var data = this.data;

	var size= this.size;
	for ( var y=0 ; y<size ; y++ ) {
	for ( var x=0 ; x<size ; x++ ) {

	    var i = y*size + x;
	    var r = Math.sqrt( data[i*2] * data[i*2] + data[i*2+1] * data[i*2+1]);
	
	    r = (Math.log(r) - min)/(max-min);
	    if (isNaN(r) || r<0) r=0;
	    
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


    }


}




