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

// This file controls the DOM elements, the actual reconstruction
// is done via a WebWorker in 'simWorker.js'



var tiffPages   = null;	    // holds images once file is read
var maxPha	= 5;	    // number of angles
var maxAng	= 3;	    // number of phases
var imageSize	= 512;	    // width and height of input images

var contrOtfVals	= new OtfVals();    // otf values

var inputPWspec = null;	    // holds power spectra of the input images

var corrImg	= null;	    // holds spectra of correlated bands (after correlation is computed)
var corrSubPxl  = null;	    // holds the subpixel data of the correlation
var maxCorr	= null;	    // holds the maxima positions

var zImported   = -1;	    // holds which z-plane has been imported
var fullResult  = null;	    // holds the full SIM reconstruction (after reconstruction is computed)
var fullResultFFT  = null;  // holds the fft'd version of the reconstruction
var fullResultWidefield	    = null;  // holds the widefield
var fullResultWidefieldFFT  = null;  // holds the fft'd widefield

var showOtfInInput = false; // toggled if the OTFs are shown
var showWidefield  = false; // if widefield is shown or reconstruction
var shwOtfOutline  = false; // if the OTFs are shown as outlines in the result spectrum


// The worker running the actual reconstruction
var simWorker	= new Worker('./simWorker.js');

// Taking care of worker results
simWorker.onmessage = function(e) {

    var what = e.data[0];

    // allows to show status to the user
    if (what=='status') {
	showStatus( e.data[1], e.data[2]);
    }

    // import complete, so show the image
    if (what=='importComplete' ) {
	inputPWspec = e.data[1];
	updateInputImageDisplay();
    }

    // sim parameter estimation complete
    if (what=='correlationComplete') {
	corrImg	    = e.data[1];
	corrSubPxl  = e.data[2];
	maxCorr	    = e.data[3];
	updateCorrelationImageDisplay();
    }
    
    // sim reconstruction complete
    if (what=='reconstructionComplete') {
	fullResult	    = e.data[1];
	fullResultFFT	    = e.data[2];
	fullResultWidefield = e.data[3];
	fullResultWidefieldFFT = e.data[4];
	updateResultImageDisplay();
    }

}

window.onload = function() { simWorker.postMessage(["confirmWorkerInit"]); };

// some code to easily write to the js console
function logger(text) {
    console.log("fairSIM-js: "+text);
}

// some code to display status messages
function showStatus( text = null, compl = null ) {

    if (text!=null) {
	document.getElementById("fairsimlogger").innerHTML = "<pre>"+text+"</pre>";
    }

    if (compl!=null) {
	if (compl>0.99) {
	    document.getElementById("fairsimprogress").innerHTML = "<pre>OK</pre>";
	    document.getElementById("fairsimlogger").style.background = "#00ff00aa";
	} else if (compl<0) {
	    document.getElementById("fairsimprogress").innerHTML = "<pre>...</pre>";
	    document.getElementById("fairsimlogger").style.background = "#aaaa00aa"; 
	} else {
	    var pc = Math.floor( compl*100 );
	    document.getElementById("fairsimprogress").innerHTML = "<pre>"+pc+"%</pre>";
	    document.getElementById("fairsimlogger").style.background = 
		"linear-gradient(to right, #00ff00aa 0%, #ee3300aa "+pc+"%)";
	}
    }
}



// called to update which raw input image to display
function setImage(pos) {

    if (tiffPages==null) {
	showStatus("import: open a TIFF first",-1);
	return;
    }

    var phaAng =  document.getElementById("sSlider").value; 

    var imgCnv = document.getElementById("rawCanvas");
    var ctx = imgCnv.getContext("2d");
    var imgData = ctx.getImageData(0,0,imgCnv.width, imgCnv.height);
    
    var data = imgData.data;

    var pha = phaAng%5;
    var ang = Math.floor(phaAng/5);

    showStatus("Input, showing ang: "+ang+" pha: "+pha+" z: "+pos);

    var max=0;
    for ( var i = 0; i<data.length; i++) {
	var val = tiffPages[ pos*maxPha + pha + ang*(tiffPages.length/maxAng)  ].data[i];
	if (val>max) max=val;
    }
    max = 255./max;    


    for ( var i = 0 ; i<data.length/4; i++) {
	data[i*4+0] = tiffPages[ pos*maxPha + pha + ang*(tiffPages.length/maxAng)  ].data[i] *max;
	data[i*4+1] = data[i*4];
	data[i*4+2] = data[i*4];
	data[i*4+3] = 0xFF; 
    }
    ctx.putImageData( imgData,0,0);
    
    if ( pos == zImported ) {
	updateInputImageDisplay();
    }
}

// when the 'import' button is clicked...
function importImages() {
    if (tiffPages==null) {
	showStatus("import: open a TIFF first",-1);
	return;
    }
    
    var z =  document.getElementById("zSlider").value; 
    zImported = z;   

    // collect the images to import for the selected plane
    var dataImport = []; 
    
    for (var ang = 0 ; ang < maxAng ; ang++ ) {
	for (var pha = 0 ; pha < maxPha ; pha++ ) {
	    var pos   = z*maxPha + pha + ang*(tiffPages.length/maxAng);
	    dataImport.push( tiffPages[ pos ].data );
	}
    }	    

    // tell the worker to grab them & fft them
    simWorker.postMessage( ['import', dataImport] );

}


// when the 'compute correlation' button is clicked
function computeCorrImages() {
    simWorker.postMessage( ['computeCorrelation'] );
}


// when then 'compute reconstruction' button is clicked
function computeReconstruction() {
    simWorker.postMessage(['computeReconstruction']);
}


// when any of the OTF values change
function updateOtfValues() {

    contrOtfVals.objNA    = document.getElementById('inp_otfna').value;
    contrOtfVals.emLambda = document.getElementById('inp_emwl').value;
    
    if ( document.getElementById('inp_attOnOff').checked ) {
	contrOtfVals.attFactor= document.getElementById('inp_attval').value;
    } else {
	contrOtfVals.attFactor=-1;
    } 

    simWorker.postMessage(['newOtfValues', contrOtfVals]);

}



// set the display of power spectra
function updateInputImageDisplay() {

    if ( inputPWspec == null || inputPWspec.length == 0 ) {
	return;
    }

    var pos = document.getElementById("sSlider").value;
    var imgCnv = document.getElementById("fftCanvas");
    var ctx = imgCnv.getContext("2d");
    var fftData = ctx.getImageData(0,0,imgCnv.width, imgCnv.height);
    var data = fftData.data;
    
    for ( var i = 0 ; i<data.length/4; i++) {
	data[i*4+0] = inputPWspec[pos][i]*255 ;
	data[i*4+1] = inputPWspec[pos][i]*255 ;
	data[i*4+2] = inputPWspec[pos][i]*255 ;
	data[i*4+3] = 0xFF;
    } 
    ctx.putImageData( fftData,0,0);

}

// set if the OTFs are to be displayed and call the update functions
function toggleOtfImage( val  ) {
    showOtfInInput = val;

    if (showOtfInInput) {
	updateOtfImage(document.getElementById("corrSlider").value);
    } else {
	updateInputImageDisplay();
	updateResultImageDisplay();
    }
}


// computes the OTF displays
function updateOtfImage( pos ) {

    if (!showOtfInInput) {
	return;
    }

    const bands = (maxPha+1)/2;
    
    var imgCnv = document.getElementById("fftCanvas");
    var ctx = imgCnv.getContext("2d");
    var fftData = ctx.getImageData(0,0,imgCnv.width, imgCnv.height);

    var otfVec = new Vec2dCplx( imageSize );
    var data = fftData.data;

    otfVec.createOtf( contrOtfVals,0,0, contrOtfVals.attFactor );
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

    	otfVec.createOtf( contrOtfVals , xo, yo, contrOtfVals.attFactor );
	var otfImg1 = otfVec.getImg(true,false);
	otfVec.createOtf( contrOtfVals ,  0,  0, contrOtfVals.attFactor );
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

}

// compute the cross-correlation display
function updateCorrelationImageDisplay( ) {

    if ( corrImg == null || corrImg.length == 0 ) {
	return;
    }


    var pos = document.getElementById("corrSlider").value;
    
    updateOtfImage(pos);

    var imgCnv = document.getElementById("corrCanvas");
    var ctx = imgCnv.getContext("2d");
    var fftData = ctx.getImageData(0,0,imgCnv.width, imgCnv.height);
    var pwSpec  = corrImg[pos];

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

function updateResultImageDisplay() {

    var slMin = document.getElementById("resMinSlider").value;
    var slMax = document.getElementById("resMaxSlider").value;
    
    if ( fullResult == null ) {
	return;
    }

    //logger("Updating image: "+slMax);
    var resData = ( showWidefield )?(fullResultWidefield):(fullResult);
    var resFFT = ( showWidefield )?(fullResultWidefieldFFT):(fullResultFFT);

    var minMax = [ Number.MAX_VALUE, Number.MIN_VALUE ];
    for (var i =0; i<resData.length;i++) {
	if ( resData[i] < minMax[0] ) minMax[0] =resData[i];
	if ( resData[i] > minMax[1] ) minMax[1] =resData[i];
    }
    
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
    
    var fftCnv = document.getElementById("resultCanvasFFT");
    var ctf = fftCnv.getContext("2d");
    var fftData = ctf.getImageData(0,0,imgCnv.width, imgCnv.height);

    var data   = imgData.data;
    var dataFFT = fftData.data;
    for ( var y = 0 ; y<imageSize*2; y++) {
        for ( var x = 0 ; x<imageSize*2; x++) {
	    var io  = x + y * imgCnv.width;
	    var ii  = x + y * imageSize*2;
	    var val = (resData[ii]-min)*scal ;
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


    const cutoffPxl    =  ((2*contrOtfVals.objNA)/(contrOtfVals.emLambda/1000.)/contrOtfVals.pxlSize);
    
    const bands = (maxPha+1)/2;
    
    ctx.putImageData( imgData,0,0);
    ctf.putImageData( fftData,0,0);
    
    ctf.beginPath();
    ctf.strokeStyle = '#ff4400';
    ctf.arc( 512, 512, 10, 0, Math.PI*2);
    ctf.stroke();

    if (showOtfOutline == true) {    
	// widefield
	ctf.beginPath();
	ctf.strokeStyle = '#aa4400';
	ctf.arc( 512, 512, cutoffPxl, 0, Math.PI*2);
	ctf.stroke();
	
	// sim otfs
	for (var ang=0; ang<maxAng; ang++) {
	    for (var b=1; b<bands; b++) {
		if (b==1) {
		    ctf.strokeStyle = '#0044aa';
		} else {
		    ctf.strokeStyle = '#00aa44';
		}

		ctf.beginPath();
		ctf.arc( 512-maxCorr[ang][0]*b, 512-maxCorr[ang][1]*b, cutoffPxl, 0, Math.PI*2);
		ctf.stroke();
		ctf.beginPath();
		ctf.arc( 512+maxCorr[ang][0]*b, 512+maxCorr[ang][1]*b, cutoffPxl, 0, Math.PI*2);
		ctf.stroke();
	    }	
	}
    }
}
