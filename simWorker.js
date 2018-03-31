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

var workerOtfVals = new OtfVals();

var imageSize = 512;	    // image size, in pxl (currently fixed at 512 !)

var inputFFTimg = null;	    // holds fft'd images once imported (as Vec2dCplx)
var bandImg	= null;	    // holds fft'd, separated SIM bands once imported (as Vec2dCplx)

var corrImg	= null;	    // holds the correlated bands (after correlation is computed)
var corrSubPxl  = null;	    // holds the subpixel data of the correlation
var maxCorr	= null;	    // holds position and phase of max. correlation -> k0-vectors

var fullResult  = null;	    // holds the full SIM reconstruction (after reconstruction is computed)
var fullResultFFT  = null;  // holds the fft'd version of the reconstruction
var fullResultWidefield	    = null;  // holds the widefield
var fullResultWidefieldFFT  = null;  // holds the fft'd widefield


// some code to easily write to the js console
function logger(text) {
    console.log("fairSIM-js: "+text);
}

// this allows to provide progress feedback and status messages
function showStatus( text=null, progress =null) {
    postMessage( [ 'status', text, progress ] );
}



// this get called for communication with the worker
onmessage = function(e) {

    var what = e.data[0];

    // just confirm (by posting a message back) the worker works
    if (what=='confirmWorkerInit') {
	showStatus("fairSIM-js: ready",1)
    }


    // import images: e.data[1] should contain all the images (typ. as Uint...)
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

    // set new otf values: e.data[1] should contain an OtfVals object 
    if (what=='newOtfValues' ) {
	workerOtfVals = e.data[1];

	var msg = "New OTF: NA "+workerOtfVals.objNA+" wl "+workerOtfVals.emLambda+" (att "+
		    ((workerOtfVals.attFactor>0)?(workerOtfVals.attFactor):("off"))+")";
	showStatus(msg);
	//logger(e.data[1].objNA);

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
    otf.createOtf(workerOtfVals, 0,0,.1); 

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
    
    var tmpOtf = new Vec2dCplx( imageSize *2 );

    // for all angles
    for ( var ang = 0; ang<maxAng; ang++) {

	showStatus("Reconstruction: ang "+ang+" band 0", (ang*bands)*.9/(maxAng*bands));
	
	tmpOtf.createOtf(workerOtfVals, 0,0,-1); 
	accOTF.add( tmpOtf );
	tmpOtf.createOtf(workerOtfVals, 0,0,workerOtfVals.attFactor); 
	accATT.add( tmpOtf );
	
	var tmpZ   = new Vec2dCplx( imageSize *2 );
	tmpZ.paste( bandImg[ang*bands],0,0 );
	// save a widefield before multiplying an OTF to it
	fullResultWidefield.paste( tmpZ, 0, 0 );
	
	tmpZ.mult(.5,0.);
	tmpZ.times( tmpOtf );
	

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
	    tmpOtf.createOtf(workerOtfVals, kx,ky,-1); 
	    accOTF.add( tmpOtf );
	    tmpOtf.createOtf(workerOtfVals, kx,ky,workerOtfVals.attFactor); 
	    accATT.add( tmpOtf );
	    
	    tmpP.times( tmpOtf );
	    
	    tmpOtf.createOtf(workerOtfVals, -kx,-ky,-1); 
	    accOTF.add( tmpOtf );
	    tmpOtf.createOtf(workerOtfVals, -kx,-ky,workerOtfVals.attFactor); 
	    accATT.add( tmpOtf );

    	    tmpN.times( tmpOtf );

	    // paste into full result
	    fullResult.paste( tmpP ,0,0);
	    fullResult.paste( tmpN ,0,0);
	    //logger("pasted angle "+ang+" band "+b+" ("+bIdx+","+cIdx+") at "+kx+" "+ky+" phase "
	    //	+(b*maxCorr[ang][3]));
	}

    }
    
    showStatus("Reconstruction: computing filters", 0.9);

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
    tmpOtf.createOtf(workerOtfVals, 0, 0,-1,1.9); 
    fullResult.times( tmpOtf );

    fullResultFFT = fullResult.duplicate();
    fullResult.fft2d(true);


    // widefield
    fullResultWidefield.maskOtf( workerOtfVals );
    fullResultWidefieldFFT = fullResultWidefield.duplicate();
    fullResultWidefield.fft2d(true);

    postMessage( ['reconstructionComplete', 
	fullResult.getImg(false,false), 
	fullResultFFT.getImg(), 
	fullResultWidefield.getImg(false,false),
	fullResultWidefieldFFT.getImg() ] );
    showStatus("Reconstruction: done",1);


}

