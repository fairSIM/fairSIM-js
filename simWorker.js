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

/* A fun project, trying to get a simple 2D SIM reconstruction
   to run in the browser, completely in JavaScript
   See the README, this is not meant for productive use. */

// This file runs the actual reconstruction, it is used by
// 'simController.js' as a WebWorker

importScripts('./Vec2dCplx.js','./fft-runner.js');


// SIM parameters and results

var maxAng = 3;		    // number of angles
var maxPha = 5;		    // number of phases 

var objNA=1.4;		    // effective NA, for estimating an OTF
var emLambda=525;	    // emission wavelength
var attFactor = 0.4;	    // factor for OTF attenuation (in fraction of OTF cutoff, useful range 0..0.4, maybe)

var imageSize = 512;	    // image size, in pxl (currently fixed at 512 !)
var pxlSize   = 0.08;	    // physical size of proj. pixel, in micron

var inputFFTimg = null;	    // holds fft'd images once imported (as Vec2dCplx)
var bandImg	= null;	    // holds fft'd, separated SIM bands once imported (as Vec2dCplx)

var corrImg	= null;	    // holds the correlated bands (after correlation is computed)
var corrSubPxl  = null;	    // holds the subpixel data of the correlation
var maxCorr	= null;

var zImported   = -1;	    // holds which z-plane has been imported
var fullResult  = null;	    // holds the full SIM reconstruction (after reconstruction is computed)
var fullResultFFT  = null;  // holds the fft'd version of the reconstruction
var fullResultWidefield	    = null;  // holds the widefield
var fullResultWidefieldFFT  = null;  // holds the fft'd widefield

var showOtfInInput = false; // toggled if the OTFs are shown
var showWidefield  = false;

// some code to easily write to the js console
function logger(text) {
    console.log("fairSIM-js: "+text);
}

// this allows to provide progress feedback and such easily
function showStatus( text=null, progress =null) {
    postMessage( [ 'status', text, progress ] );
}



// this get called for communication with the worker
onmessage = function(e) {

    var what = e.data[0];

    // import images
    if (what=='import') {

	var inputPWspec = [];

	inputFFTimg = [];
	bandImg     = [];
	const bands = (maxPha+1)/2;

	for (var ang = 0 ; ang < maxAng ; ang++ ) {
	    for (var pha = 0 ; pha < maxPha ; pha++ ) {
	
		var pos   = ang*maxPha + pha;
		var fftIn = new Vec2dCplx( imageSize );
		fftIn.setFromTiff( e.data[1][pos] );
		fftIn.applyWindow(20);
		fftIn.fft2d(false);
		inputFFTimg.push( fftIn );
	    
		inputPWspec.push( fftIn.getImg());
   
		showStatus("fft'in input images", pos/(maxAng*maxPha*1.)); 
	    }
	
	    // band-separate 
	    var tmp = bandSeparate( inputFFTimg, ang*maxPha );
	    for (var b=0; b<bands; b++) {
		bandImg.push(tmp[b]);
	    }
	} 
    
	showStatus("Images imported",1);
	postMessage(['importComplete', inputPWspec ] );

    }

    // compute the correlation
    if (what=='computeCorrelation') {
	computeCorrImages();
    }
    
    // compute the reconstruction
    if (what=='computeReconstruction') {
	computeReconstruction();
    }

}

// called by 'importImages()', calculates a (very basic)
// band separation (only equi-distant phases, only #bands=(#phases+1)/2)
function bandSeparate(inVec, offset) {

    const bands = (maxPha+1)/2;
    var ret = [];
    
    for ( var b = 0 ; b<bands; b++) {

	bVec = new Vec2dCplx( inVec[0].size );
	//bVecN = new Vec2dCplx( inVec[0].size );

	for ( var i = 0 ; i < maxPha ; i++ ) {
	
	    var pha = 2*b*Math.PI * i / maxPha;	
	   
	    var addVec = inVec[i+offset].duplicate();  
	    addVec.mult( Math.cos(pha), Math.sin(pha));
	    bVec.add( addVec );
	   
	    /* 
	    // generate the conj. bands -1, -2 
	    if (b!=0) {
		var addVecN = inVec[i+offset].duplicate();  
		addVecN.mult( Math.cos(-pha), Math.sin(-pha));
		bVecN.add( addVecN );
	    } */
	
	}
    
	ret.push( bVec );
   
	/*
	// cross-check if the conj. bands are generated correctly
	if (b!=0) { 
	    // cross-check: there is no need to generate the conj.
	    // bands explicitly, they should just be band[i].conj().
	    // However, that failed before...
	    var bVecN2 = bVec.duplicateMirrored();
	    bVecN2.conj();	
	    logger( "cross-check: "+b+" "+bVecN2.comp( bVecN ));
	}
	*/

    }
    return ret;
}

// calculates the cross-correlation between bands
function computeCorrImages( minDist = 100 ) {
    
    var bands = (maxPha+1)/2;

    if (bandImg == null || bandImg.length != bands*maxAng ) {
	showStatus("import some images first...",-1 );
	return;
    }

    corrImg = [];
    corrSubPxl = [];
    maxCorr = [];

    var otf = new Vec2dCplx(imageSize);
    createOtf( otf,0,0,.1); 

    showStatus("Fitting SIM parameters..",0);
    
    for (var ang = 0 ; ang < maxAng ; ang++ ) {

	var bZero = bandImg[ ang*bands ].duplicate();
	bZero.times(otf);
	bZero.fft2d(true);

	for ( var b=1; b<bands; b++) {

	    var cImg = bandImg[ (b +ang*bands ) ].duplicate();

	    //c.times(otf);
    	    cImg.times(otf);
	    cImg.fft2d(true);
	    cImg.times( bZero );
	    cImg.fft2d(false);
    
	    // find the brightest pixel
	    var max = cImg.findMax(60*b);
	    
	    corrImg.push(cImg.getImg());

	    // do a subpixel-precision fit
	    if ( b==2) {

		// TODO: This should actually compute the common regions
		// between bands first ...

		bHigh = bandImg[ (b +ang*bands ) ].duplicate();
		bHigh.times(otf);
		bHigh.fft2d(true);
	
		// compute the correlation
		for ( var iter = 0; iter<2; iter++) {
		    
		    var maxX=0, maxY=0;
		    var spTmp = [], maxValue=0, minValue = Number.MAX_VALUE;

		    for ( var x=0; x<8; x++ ) {
			for ( var y=0; y<8; y++ ) {
			    var kx = max[0] + ((x-3.5)/3.5)*((iter==0)?(2):(0.5));
			    var ky = max[1] + ((y-3.5)/3.5)*((iter==0)?(2):(0.5));
			    var corr = bZero.crossCorrelate( bHigh, kx,ky);
			    if ( corr[2] > maxValue ) {
				maxValue = corr[2];
				maxX = kx;
				maxY = ky;
			    }
			    if ( corr[2] < minValue ) {
				minValue = corr[2];
			    }
			    spTmp.push( corr );
			}	
			showStatus(null, ((ang*128 + iter*64 + x*8)/(128.*maxAng)));
		    }

		    for ( var i=0; i<64; i++) {
			spTmp[i][2] = (spTmp[i][2]-minValue)/(maxValue-minValue);
		    }   

		    corrSubPxl.push( spTmp ); 
		    max[0] = maxX;
		    max[1] = maxY;
		}
	    
		// store the coordinates
		max[0] =  (max[0]<imageSize/2)?(max[0]):(max[0]-imageSize);
		max[1] =  (max[1]<imageSize/2)?(max[1]):(max[1]-imageSize);
		max[0] /= 2;		
		max[1] /= 2;		

		// now, retrieve the phase from correlation band 0 to band 1
		// wrap the position to -size/2 .. size/2
    		var bLow  = bandImg[ang*bands+1].duplicate();
		bLow.times( otf );
		bLow.fft2d(true);

		var valMagPha = bZero.crossCorrelate( bLow, max[0], max[1]);


		max[3] = -Math.atan2( valMagPha[1], valMagPha[0]);
		showStatus(("Fit: ang "+ang+" band "+b+" -> "+
		    " kx "+max[0].toFixed(2)+
		    " ky "+max[1].toFixed(2)+" pha "+max[3].toFixed(3)), null);
 
		maxCorr.push( max );
	    }

	}
    }    
    
    postMessage(['correlationComplete', corrImg, corrSubPxl, maxCorr ]); 
    showStatus("SIM parameters computed",1);
    

}

function setFixedParameters() {

    maxCorr = [];
    maxCorr.push( new Array( 137.433/2,  140.90/2, 0.8, 0.886 ) );
    maxCorr.push( new Array( -52.856/2,  189.47/2, 0.8, 2.730 ) );
    maxCorr.push( new Array( 190.078/2, -49.967/2, 0.8, 1.872 ) );
}



function computeReconstruction() {

    
    if ( maxCorr == null || maxCorr.length == 0 ) {
	showStatus("Please run the correlation estimation first",-1);
	return;
    } 

    showStatus("Reconstruction: starting",0);
    
    const bands = (maxPha+1)/2;
    fullResult		= new Vec2dCplx( 2*imageSize );
    fullResultWidefield = new Vec2dCplx( 2*imageSize );

    // sum up std. OTF
    var accOTF = new Vec2dCplx( imageSize *2 );
    
    // sum up att. OTF
    var accATT = new Vec2dCplx( imageSize *2 );
    
    var tmpOTF = new Vec2dCplx( imageSize *2 );

    // for all angles
    for ( var ang = 0; ang<maxAng; ang++) {

	
	createOtf( tmpOTF, 0,0, -1. );
	accOTF.add( tmpOTF );
	createOtf( tmpOTF, 0,0, attFactor );
	accATT.add( tmpOTF );
	
	var tmpZ   = new Vec2dCplx( imageSize *2 );
	tmpZ.paste( bandImg[ang*bands],0,0 );
	// save a widefield before multiplying an OTF to it
	fullResultWidefield.paste( tmpZ, 0, 0 );
	
	tmpZ.mult(.5,0.);
	tmpZ.times( tmpOTF );
	

	fullResult.paste( tmpZ, 0, 0 );
	
	// for all bands
	for ( var b=1; b<bands; b++) {
	    
	    showStatus("Reconstruction: ang "+ang+" band "+b, (ang*bands+b)*.9/(maxAng*bands));
	    
	    // multiply input w. otf
	    var bIdx = ang*bands+b;
	    var cIdx = ang*(bands-1)+b-1;

	    //bandImg[ bIdx ].multPhase( maxCorr[cIdx][3] );

	    var tmpP   = new Vec2dCplx( imageSize *2 );
	
	    tmpP.paste( bandImg[bIdx], 0,0);
	    tmpN = tmpP.duplicateMirrored();    
	    tmpN.conj();

	
	    /*
	    tmpP.paste( bandImg[bIdx+b*2-1], 0,0);
	    tmpN.paste( bandImg[bIdx+b*2-1], 0,0, true);
	    tmpP.multPhase( b*maxCorr[cIdx][3] );
	    tmpN.multPhase(-b*maxCorr[cIdx][3] );
	    */
	    //tmpN.conj();	    

	    // go to real space
	    //bandImg[bIdx].fft2d(true);
	    tmpP.fft2d(true);
	    tmpN.fft2d(true);

	    var kx = maxCorr[ang][0]*b;
	    var ky = maxCorr[ang][1]*b;

	    // subpixel fourier-shift the band
	    tmpP.fourierShift( kx,ky );
	    tmpN.fourierShift(-kx,-ky );
	    
	    // move back to freq. space
	    tmpP.fft2d(false);
	    tmpN.fft2d(false);
	    
	    tmpP.multPhase( b*maxCorr[ang][3] );
	    tmpN.multPhase(-b*maxCorr[ang][3] );

	    // multiply w. OTF
	    createOtf( tmpOTF, kx,ky, -1. );
	    accOTF.add( tmpOTF );
	    createOtf( tmpOTF, kx,ky, attFactor );
	    accATT.add( tmpOTF );
	    
	    tmpP.times( tmpOTF );
	    
	    createOtf( tmpOTF, -kx,-ky, -1. );
	    accOTF.add( tmpOTF );
	    createOtf( tmpOTF, -kx,-ky, attFactor );
	    accATT.add( tmpOTF );

    	    tmpN.times( tmpOTF );

	    // paste into full result
	    fullResult.paste( tmpP ,0,0);
	    fullResult.paste( tmpN ,0,0);
	    //logger("pasted angle "+ang+" band "+b+" ("+bIdx+","+cIdx+") at "+kx+" "+ky+" phase "
	    //	+(b*maxCorr[ang][3]));
	}

    }
    
    showStatus("Reconstruction: computing filters", 0.95);

    // Wiener filtering
    for ( var i=0; i<fullResult.length; i++) {
	// accOTF holds the OTF w/o attenuation
	// accATT holds the OTF multiplied by attenuation
	// so the product it OTF^2 * ATT as required
	var d = 1/(  accOTF.data[2*i] * accATT.data[2*i] + 0.5 );
	fullResult.data[2*i+0] *= d;
	fullResult.data[2*i+1] *= d;
    }

    // APO
    createOtf( tmpOTF, 0,0,-1,1.9);
    fullResult.times( tmpOTF );

    fullResultFFT = fullResult.duplicate();
    fullResult.fft2d(true);


    // widefield
    maskOtf( fullResultWidefield );
    fullResultWidefieldFFT = fullResultWidefield.duplicate();
    fullResultWidefield.fft2d(true);

    postMessage( ['reconstructionComplete', 
	fullResult.getImg(false,false), 
	fullResultFFT.getImg(), 
	fullResultWidefield.getImg(false,false),
	fullResultWidefieldFFT.getImg() ] );
    showStatus("Reconstruction: done",1);


}


function valOtf( dist ) {
    if ((dist<0)||(dist>=1))
	return 0.0;
    if (dist == 1 )
	return 1.0;
    return (2/Math.PI)*(Math.acos(dist) - dist*Math.sqrt(1-dist*dist));
}


// create a simple, 2D OTF
function createOtf( vec, kx=0, ky=0, att=-1, coShift=1 ) {

    const cyclPxl   =  1./(imageSize*pxlSize);
    const cutoff    =  ((2*objNA)/(emLambda/1000.))*coShift;
    const cutoffPxl =  cutoff/cyclPxl;

    for (var y=0; y<vec.size; y++) {
	for (var x=0; x<vec.size; x++) {

	    var xi = x + kx;
	    var yi = y + ky;

	    var xh = ((xi<vec.size/2)?( xi):(xi-vec.size)) ;
	    var yh = ((yi<vec.size/2)?( yi):(yi-vec.size)) ;

	    var dist = Math.sqrt( xh*xh+yh*yh );
	    var val  = valOtf( dist/cutoffPxl );

	    if ( att>0 && (dist/cutoffPxl) < 1 ) {
		val *= 1.0-0.99*Math.exp( -(dist/cutoffPxl/att));
	    }

	    vec.data[(x+y*vec.size)*2+0] = val;
	    vec.data[(x+y*vec.size)*2+1] = 0.0;

	}
    }

}


// cut out regions beyond OTF support
function maskOtf( vec , coShift = 1.) {

    const cyclPxl   =  1./(imageSize*pxlSize);
    const cutoff    =  ((2*objNA)/(emLambda/1000.))*coShift;
    const cutoffPxl =  cutoff/cyclPxl;

    for (var y=0; y<vec.size; y++) {
	for (var x=0; x<vec.size; x++) {

	    var xh = ((x<vec.size/2)?( x):(x-vec.size)) ;
	    var yh = ((y<vec.size/2)?( y):(y-vec.size)) ;
	    var dist = Math.sqrt( xh*xh+yh*yh );

	    if ( (dist/cutoffPxl) >1.1) {
		vec.data[(x+y*vec.size)*2+0] = 0.0;
		vec.data[(x+y*vec.size)*2+1] = 0.0;
	    } else if ( (dist/cutoffPxl) >1 ) {
		var d2 = ((dist/cutoffPxl) -1.)*10.;
		var v  = .5*(1+Math.cos(Math.PI*d2));
		vec.data[(x+y*vec.size)*2+0] *= v;
		vec.data[(x+y*vec.size)*2+1] *= v;
	    }
	}
    }

}





function updateFFTimage( pos ) {

    if ( inputFFTimg == null || inputFFTimg.length == 0 ) {
	return;
    }

    var imgCnv = document.getElementById("fftCanvas");
    var ctx = imgCnv.getContext("2d");
    var fftData = ctx.getImageData(0,0,imgCnv.width, imgCnv.height);
    var pwSpec  = inputFFTimg[pos].getImg();


    var data = fftData.data;
    for ( var i = 0 ; i<data.length/4; i++) {
	data[i*4+0] = pwSpec[i]*255 ;
	data[i*4+1] = pwSpec[i]*255 ;
	data[i*4+2] = pwSpec[i]*255 ;
	data[i*4+3] = 0xFF;
    } 
    ctx.putImageData( fftData,0,0);

}


function toggleOtfImage( val  ) {
    showOtfInInput = val;

    if (showOtfInInput) {
       updateOtfImage(document.getElementById("corrSlider").value);
    } else {
	updateFFTimage(document.getElementById("sSlider").value);
	updateResultImage(
	    document.getElementById("resMinSlider").value, 
	    document.getElementById("resMaxSlider").value);
    }
}


function updateOtfImage( pos ) {

    const bands = (maxPha+1)/2;
    
    var imgCnv = document.getElementById("fftCanvas");
    var ctx = imgCnv.getContext("2d");
    var fftData = ctx.getImageData(0,0,imgCnv.width, imgCnv.height);

    var otfVec = new Vec2dCplx( imageSize );
    var data = fftData.data;

    createOtf( otfVec , 0, 0, attFactor );
    var otfImg = otfVec.getImg(true,false);
    
    for ( var i = 0 ; i<data.length/4; i++) {
	data[i*4+0] = otfImg[i]*255 ;
	data[i*4+1] = otfImg[i]*255 ;
	data[i*4+2] = otfImg[i]*255 ;
	data[i*4+3] = 0xFF;
    } 

    ctx.putImageData( fftData,0,0);
   

    // output the overlay 
    if (maxCorr != null && maxCorr.length != 0 ) {
    
	var imgCnv = document.getElementById("resultCanvas");
	var ctx = imgCnv.getContext("2d");
	var fftData = ctx.getImageData(0,0,imgCnv.width, imgCnv.height);

	var otfVec = new Vec2dCplx( imageSize*2 );
	var data = fftData.data;
    
	var xo =  maxCorr[Math.floor(pos/2)][0] * (1+pos%2);
	var yo =  maxCorr[Math.floor(pos/2)][1] * (1+pos%2);

    	createOtf( otfVec , xo, yo, attFactor );
	var otfImg1 = otfVec.getImg(true,false);
	createOtf( otfVec , 0, 0, attFactor );
	var otfImg0 = otfVec.getImg(true,false);
	
	for ( var i = 0 ; i<data.length/4; i++) {
	    
	    var mark = ( otfImg0[i]>0.05 && otfImg1[i] > 0.05 )?(1):(0);
    
	    data[i*4+0] = otfImg0[i]*255 ;
	    data[i*4+1] = otfImg1[i]*255 ;
	    data[i*4+2] = mark*100;
	    data[i*4+3] = 0xFF;
	} 
	
	ctx.putImageData( fftData,0,0);
    
    }

    // compute the full output
    /*
    if (maxCorr != null && maxCorr.length != 0 ) {
	var imgCnv = document.getElementById("resultCanvasFFT");
	var ctx = imgCnv.getContext("2d");
	var fftData = ctx.getImageData(0,0,imgCnv.width, imgCnv.height);

	var otf0 = new Vec2dCplx( imageSize*2 );
	var otf1 = new Vec2dCplx( imageSize*2 );
	var otf2 = new Vec2dCplx( imageSize*2 );
	var tmp  = new Vec2dCplx( imageSize*2 );
	var data = fftData.data;
   
	createOtf( otf0 , 0, 0, attFactor );

	for (var ang =0; ang<maxAng; ang++) {
 
	    var xo =  maxCorr[ang][0] ;
	    var yo =  maxCorr[ang][1] ;

	    createOtf( tmp , xo, yo, attFactor );
	    otf1.add( tmp );	
	    createOtf( tmp , -xo, -yo, attFactor );
	    otf1.add( tmp );	
	    createOtf( tmp , 2*xo, 2*yo, attFactor );
	    otf2.add( tmp );	
	    createOtf( tmp , -2*xo, -2*yo, attFactor );
	    otf2.add( tmp );	
	}	

	var dat0 = otf0.getImg(true,false);
	var dat1 = otf1.getImg(true,false);
	var dat2 = otf2.getImg(true,false);

	for ( var i = 0 ; i<data.length/4; i++) {
	    
	    data[i*4+0] = dat0[i]*100 ;
	    data[i*4+1] = dat1[i]*100 ;
	    data[i*4+2] = dat2[i]*100 ;
	    data[i*4+3] = 0xFF;
	} 
	
	ctx.putImageData( fftData,0,0);
    } 
    */

}


function updateCorrelationImage( pos ) {

    if ( corrImg == null || corrImg.length == 0 ) {
	return;
    }

    if ( showOtfInInput ) {
	updateOtfImage(pos);
    }	 

    var imgCnv = document.getElementById("corrCanvas");
    var ctx = imgCnv.getContext("2d");
    var fftData = ctx.getImageData(0,0,imgCnv.width, imgCnv.height);
    var pwSpec  = corrImg[pos].getImg(true);


    var data = fftData.data;
    for ( var i = 0 ; i<data.length/4; i++) {
	data[i*4+0] = pwSpec[i]*255 ;
	data[i*4+1] = pwSpec[i]*255 ;
	data[i*4+2] = pwSpec[i]*255 ;
	data[i*4+3] = 0xFF;
    } 

    const b = 1 + (pos%2);
    const ang  = Math.floor(pos/2);

    if ( b ==2 ) {
	
	for ( var iter=0; iter<2; iter++) {
	    var bpos = Math.floor(pos/2)*2 +iter;
	    //logger( "pos -> "+pos+" bpos "+bpos);
	    for ( var x=0; x<8; x++)
	    for ( var y=0; y<8; y++) {

		var val = corrSubPxl[ bpos ][x+y*8][2];
		//logger( x+" "+y+" "+val);
		for ( var xz=0; xz<4; xz++)
		for ( var yz=0; yz<4; yz++) {
		    var i = (x*4+xz) + (y*4+yz)*imageSize + imageSize-(40*(iter+1));
		    data[i*4+1] = val*255;
		    data[i*4+2] = val*255;
		    data[i*4+0] = 0;
		}
	    }
	}

    }

    ctx.putImageData( fftData,0,0);

    //logger("displ: ang "+ang+"  band "+b);

    ctx.beginPath();
    var xo =  maxCorr[ang][0]*b+imageSize/2;
    var yo =  maxCorr[ang][1]*b+imageSize/2;
    ctx.stokeStyle = '#ff4400';
    ctx.arc( xo, yo, 4, 0, Math.PI*2);
    if ( pos%2 ==1 ) {
	ctx.moveTo(xo,yo-4),
	ctx.lineTo( 512-44, 32);
    }
    ctx.stroke();

}

function updateResultImage(slMin=-1, slMax=100) {

    if ( fullResult == null ) {
	return;
    }

    //logger("Updating image: "+slMax);
    var whatToDisplay = ( showWidefield )?(fullResultWidefield):(fullResult);
    var whatToDisplayFFT = ( showWidefield )?(fullResultWidefieldFFT):(fullResultFFT);

    var minMax = whatToDisplay.getRealMinMax();
    //var scal = 255./(minMax[1]-minMax[0]);
    var scal,min,max;
    if (slMin==-1){
	scal = 255./(minMax[1]*(slMax/100.));
	min = 0;
    } else {
	min = minMax[0]	+ (minMax[1]-minMax[0])*.5*slMin/100.;
	max = minMax[0] + (minMax[1]-minMax[0])*slMax/100.;
	if (min>max-5) max=min+5;
	scal = 255./(max-min);
    }

    var imgCnv = document.getElementById("resultCanvas");
    var ctx = imgCnv.getContext("2d");
    var imgData = ctx.getImageData(0,0,imgCnv.width, imgCnv.height);
    var resData  = whatToDisplay.data;
    
    var fftCnv = document.getElementById("resultCanvasFFT");
    var ctf = fftCnv.getContext("2d");
    var fftData = ctf.getImageData(0,0,imgCnv.width, imgCnv.height);
    var resFFT   = whatToDisplayFFT.getImg(true,true);



    var data   = imgData.data;
    var dataFFT = fftData.data;
    for ( var y = 0 ; y<fullResult.size; y++) {
        for ( var x = 0 ; x<fullResult.size; x++) {
	    var io  = x + y * imgCnv.width;
	    var ii  = x + y * fullResult.size;
	    var val = (resData[2*ii]-min)*scal ;
	    if (val<0) val = 0;
	    data[io*4+0] = val;
	    data[io*4+1] = val;
	    data[io*4+2] = val;
	    data[io*4+3] = 0xFF;
	    dataFFT[io*4+0] = resFFT[ii]*255;
	    dataFFT[io*4+1] = resFFT[ii]*255;
	    dataFFT[io*4+2] = resFFT[ii]*255;
	    dataFFT[io*4+3] = 0xFF;
	}
    } 

    ctx.putImageData( imgData,0,0);
    ctf.putImageData( fftData,0,0);
    ctf.beginPath();
    ctf.stokeStyle = '#ff4400';
    ctf.arc( 512, 512, 10, 0, Math.PI*2);
    ctf.stroke();


}
