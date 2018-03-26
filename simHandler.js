function logger(text) {

    console.log("fairSIM-js: "+text);
    document.getElementById("fairsimlogger").style.visibility = 'hidden';
    document.getElementById("fairsimlogger").innerHTML = "<pre>"+text+"<pre>";
    document.getElementById("fairsimlogger").style.visibility = 'visible';
}


var maxAng = 3;
var maxPha = 5;

var objNA=1.4;
var emLambda=525;
var attFactor = 0.4;

var imageSize = 512;
var pxlSize   = 0.08;

var tiffPages   = null;
var inputFFTimg = null;
var bandImg	= null;

var corrImg	= null;
var corrSubPxl  = null;
var maxCorr	= null;

var zImported   = -1;
var fullResult  = null;
var fullResultFFT  = null;

var showOtfInInput = false;

function setImage(pos) {

    if (tiffPages==null) {
	logger("open a TIFF first");
	return;
    }

    var phaAng =  document.getElementById("sSlider").value; 

    var imgCnv = document.getElementById("rawCanvas");
    var ctx = imgCnv.getContext("2d");
    var imgData = ctx.getImageData(0,0,imgCnv.width, imgCnv.height);
    
    var data = imgData.data;

    var pha = phaAng%5;
    var ang = Math.floor(phaAng/5);

    logger("ang: "+ang+" pha: "+pha+" z: "+pos);

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
	updateFFTimage(phaAng);
    }
}

// imports the images
function importImages() {
    if (tiffPages==null) {
	logger("Open a tiff first");
	return;
    }
    
    var z =  document.getElementById("zSlider").value; 
    zImported = z;   
 
    inputFFTimg = [];
    bandImg     = [];
    
    const bands = (maxPha+1)/2;

    for (var ang = 0 ; ang < maxAng ; ang++ ) {
	
	// copy input data
	for (var pha = 0 ; pha < maxPha ; pha++ ) {
	    var pos   = z*maxPha + pha + ang*(tiffPages.length/maxAng);
	    var fftIn = new Vec2dCplx( tiffPages[ pos  ].width );
	    fftIn.setFromTiff( tiffPages[ pos ].data );
	    fftIn.applyWindow(20);
	    inputFFTimg.push( fftIn );
	}
	    
	// fft input image
	for (var pha = 0 ; pha < maxPha ; pha++ ) {
	    inputFFTimg[ pha+ang*maxPha].fft2d(false);

	    /*
    	    // used to cross-check fft(k) = fft*(-k)
	    var tmp = inputFFTimg[ pha+ang*maxPha].duplicateMirrored();
	    tmp.conj();
	    logger("diff: "+tmp.comp( inputFFTimg[ pha+ang*maxPha] ));
	    tmp.mult(-1.,0.);    
	    inputFFTimg[ pha+ang*maxPha].add(tmp);
	    */
	} 
	
	// band-separate 
	tmp = bandSeparate( inputFFTimg, ang*maxPha );
	for (var b=0; b<bands; b++) {
	    bandImg.push(tmp[b]);
	}

    }

    updateFFTimage(document.getElementById("sSlider").value);
    logger("slice "+z+" imported");
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
	logger("import some images first..." );
	return;
    }

    corrImg = [];
    corrSubPxl = [];
    maxCorr = [];

    var otf = new Vec2dCplx(imageSize);
    createOtf( otf,0,0,.1); 

    for (var ang = 0 ; ang < maxAng ; ang++ ) {

	var bZero = bandImg[ ang*bands ].duplicate();
	bZero.times(otf);
	bZero.fft2d(true);

	for ( var b=1; b<bands; b++) {

	    var c = bandImg[ (b +ang*bands ) ].duplicate();

	    //c.times(otf);
    	    c.times(otf);
	    c.fft2d(true);
	    c.times( bZero );
	    c.fft2d(false);
	    corrImg.push(c);
    
	    // find the brightest pixel
	    var max = c.findMax(60*b);

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
		    logger( "Fitting angle "+ang+", iteration "+iter);

		    for ( var x=0; x<8; x++ )
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
		logger("paramFit: "+ang+" "+b+" -> "+max[0]+" "+max[1]+" pha "+max[3]);
 
		maxCorr.push( max );
	    }

	}
    }    
    
    //logger(maxCorr);

    updateCorrelationImage(document.getElementById("corrSlider").value);
}

function setFixedParameters() {

    maxCorr = [];
    maxCorr.push( new Array( 137.433/2,  140.90/2, 0.8, 0.886 ) );
    maxCorr.push( new Array( -52.856/2,  189.47/2, 0.8, 2.730 ) );
    maxCorr.push( new Array( 190.078/2, -49.967/2, 0.8, 1.872 ) );
}



function computeReconstruction() {

    
    if ( maxCorr == null || maxCorr.length == 0 ) {
	logger("please run the correlation estimation first");
	return;
    } 

    const bands = (maxPha+1)/2;
    fullResult = new Vec2dCplx( 2*imageSize );

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
	tmpZ.mult(.5,0.);
	tmpZ.times( tmpOTF );
	

	fullResult.paste( tmpZ, 0, 0 );
	
	// for all bands
	for ( var b=1; b<bands; b++) {
	    
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
	    logger("pasted angle "+ang+" band "+b+" ("+bIdx+","+cIdx+") at "+kx+" "+ky+" phase "
		+(b*maxCorr[ang][3]));
	}

    }

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

    updateResultImage();


}


function valOtf( dist ) {
    if ((dist<0)||(dist>=1))
	return 0.0;
    if (dist == 1 )
	return 1.0;
    return (2/Math.PI)*(Math.acos(dist) - dist*Math.sqrt(1-dist*dist));
}


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



function updateFFTimage( pos ) {

    showOtfInInput = false;

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


function showOtfImage() {
   showOtfInInput = true;
   updateOtfImage(document.getElementById("corrSlider").value);
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

    	createOtf( otfVec , xo, yo, .1 );
	var otfImg1 = otfVec.getImg(true,false);
	createOtf( otfVec , 0, 0, .1 );
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

    logger("displ: ang "+ang+"  band "+b);

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

    var minMax = fullResult.getRealMinMax();
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
    var resData  = fullResult.data;
    
    var fftCnv = document.getElementById("resultCanvasFFT");
    var ctf = fftCnv.getContext("2d");
    var fftData = ctf.getImageData(0,0,imgCnv.width, imgCnv.height);
    var resFFT   = fullResultFFT.getImg(true,true);



    var data   = imgData.data;
    var dataFFT = fftData.data;
    for ( var y = 0 ; y<fullResult.size; y++) {
        for ( var x = 0 ; x<fullResult.size; x++) {
	    var io  = x + y * imgCnv.width;
	    var ii  = x + y * fullResult.size;
	    data[io*4+0] = (resData[2*ii]-min)*scal ;
	    data[io*4+1] = (resData[2*ii]-min)*scal ;
	    data[io*4+2] = (resData[2*ii]-min)*scal ;
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
