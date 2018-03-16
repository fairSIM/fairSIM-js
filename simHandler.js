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
var imageSize = 512;

var tiffPages   = null;
var inputFFTimg = null;
var bandImg	= null;

var corrImg	= null;
var corrSubPxl  = null;
var maxCorr	= null;

var zImported   = -1;
var fullResult  = null;


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

/*
    var max=0;
    for ( var i = 0; i<data.length/4; i++) {
	var val = tiffPages[ pos*maxPha + pha + ang*(tiffPages.length/maxAng)  ].data[i];
	if (val>max) max=val;
    }
*/
    for ( var i = 0 ; i<data.length; i++) {
	data[i] = tiffPages[ pos*maxPha + pha + ang*(tiffPages.length/maxAng)  ].data[i] ;
/*	data[i*4+1] = data[i*4];
	data[i*4+2] = data[i*4];
	data[i*4+3] = 0xFF; */
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

    for (var ang = 0 ; ang < maxAng ; ang++ ) {
	
	// copy input data
	for (var pha = 0 ; pha < maxPha ; pha++ ) {
	    var pos   = z*maxPha + pha + ang*(tiffPages.length/maxAng);
	    var fftIn = new Vec2dCplx( tiffPages[ pos  ].width );
	    fftIn.setFromRGB( tiffPages[ pos ].data );
	    fftIn.applyWindow(20);
	    inputFFTimg.push( fftIn );
	}

	// band-separate (in real space)
	tmp = bandSeparate( inputFFTimg, ang*maxPha );
	for (var b=0; b<(maxPha+1)/2; b++) {
	    bandImg.push(tmp[b]);
	}
	for (var pha = 0 ; pha < maxPha ; pha++ ) {
	    inputFFTimg[ pha+ang*maxPha].fft2d();
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

	for ( var i = 0 ; i < maxPha ; i++ ) {
	
	    var pha = 2*b*Math.PI * i / maxPha;	
	   
	    addVec = inVec[i+offset].copy();  
	    addVec.mult( Math.cos(pha), Math.sin(pha));
	    bVec.add( addVec );
	}
    
	ret.push( bVec );
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

    for (var ang = 0 ; ang < maxAng ; ang++ ) {
	for ( var b=1; b<bands; b++) {

	    var c = bandImg[ (b +ang*bands ) ].copy();
	    c.times( bandImg[  ang*bands ]);
	    c.fft2d(false);
	    corrImg.push(c);
    
	    // find the brightest pixel
	    var max = c.findMax();

	    // do a subpixel-precision fit
	    if ( b==2) {
		for ( var iter = 0; iter<2; iter++) {
		    
		    var maxX=0, maxY=0;
		    var spTmp = [], maxValue=0, minValue = Number.MAX_VALUE;
		    logger( "Fitting angle "+ang+", iteration "+iter);
	
		    setTimeout(  function(bImg) {	
			for ( var x=0; x<10; x++ )
			for ( var y=0; y<10; y++ ) {
			    var kx = max[0] + ((x-4.5)/4.5)*((iter==0)?(2):(0.5));
			    var ky = max[1] + ((y-4.5)/4.5)*((iter==0)?(2):(0.5));
			    var corr = bImg[ (ang*bands ) ].crossCorrelate( bImg[ (b +ang*bands ) ], kx,ky)[2];
			    if ( corr > maxValue ) {
				maxValue = corr;
				maxX = kx;
				maxY = ky;
			    }
			    if ( corr < minValue ) {
				minValue = corr;
			    }
			    spTmp.push( corr );
			}

			for ( var i=0; i<100; i++) {
			    spTmp[i] = (spTmp[i]-minValue)/(maxValue-minValue);
			}   
		    }(bandImg), 200);

		    corrSubPxl.push( spTmp ); 
		    max[0] = maxX;
		    max[1] = maxY;
		}
	    }
	    
	    maxCorr.push( max );

	}
    }    
    
    //logger(maxCorr);

    updateCorrelationImage(document.getElementById("corrSlider").value);
}


function computeReconstruction() {

    
    if ( maxCorr == null || maxCorr.length == 0 ) {
	logger("please run the correlation estimation first");
	return;
    } 

    const bands = (maxPha+1)/2;
    fullResult = new Vec2dCplx( 2*imageSize );
    
    for ( var ang = 0; ang<maxAng; ang++) {

	bandImg[ang*bands].fft2d(true);
	fullResult.paste( bandImg[ang*bands] ,0,0);
	for ( var b=0; b<bands; b++) {
	    bandImg[ang*bands+b].fft2d(true);
	    fullResult.paste( bandImg[ang*bands+b] ,Math.floor( maxCorr[0]),Math.floor(maxCorr[0]));
	}

    }

    fullResult.fft2d(false);

    updateResultImage();


}




function createOtf( vec, kx=0, ky=0 ) {

    



}

function updateFFTimage( pos ) {

    if ( inputFFTimg == null || inputFFTimg.length == 0 ) {
	return;
    }

    var imgCnv = document.getElementById("fftCanvas");
    var ctx = imgCnv.getContext("2d");
    var fftData = ctx.getImageData(0,0,imgCnv.width, imgCnv.height);
    var pwSpec  = inputFFTimg[pos].pwSpec();

    //logger( pwSpec.length);

    var data = fftData.data;
    for ( var i = 0 ; i<data.length/4; i++) {
	data[i*4+0] = pwSpec[i]*255 ;
	data[i*4+1] = pwSpec[i]*255 ;
	data[i*4+2] = pwSpec[i]*255 ;
	/*data[i*4+0] = bandImg[pos%9].data[i*2] ;
	data[i*4+1] = bandImg[pos%9].data[i*2] ;
	data[i*4+2] = bandImg[pos%9].data[i*2] ;*/
	data[i*4+3] = 0xFF;
    } 
    ctx.putImageData( fftData,0,0);

}


function updateCorrelationImage( pos ) {

    if ( corrImg == null || corrImg.length == 0 ) {
	return;
    }

    var imgCnv = document.getElementById("corrCanvas");
    var ctx = imgCnv.getContext("2d");
    var fftData = ctx.getImageData(0,0,imgCnv.width, imgCnv.height);
    var pwSpec  = corrImg[pos].pwSpec();

    //logger( pwSpec.length);

    var data = fftData.data;
    for ( var i = 0 ; i<data.length/4; i++) {
	data[i*4+0] = pwSpec[i]*255 ;
	data[i*4+1] = pwSpec[i]*255 ;
	data[i*4+2] = pwSpec[i]*255 ;
	data[i*4+3] = 0xFF;
    } 


    if ( pos%2 ==1 ) {
	
	for ( var iter=0; iter<2; iter++) {
	    var bpos = Math.floor(pos/2)*2 +iter;
	    //logger( "pos -> "+pos+" bpos "+bpos);
	    for ( var x=0; x<10; x++)
	    for ( var y=0; y<10; y++) {

		var val = corrSubPxl[ bpos ][x+y*10];
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


    ctx.beginPath();
    var x =  maxCorr[pos][0];
    var y =  maxCorr[pos][1];
    var xo = (x<imageSize/2)?(x+imageSize/2):(x-imageSize/2);
    var yo = (y<imageSize/2)?(y+imageSize/2):(y-imageSize/2);
    ctx.stokeStyle = '#ff4400';
    ctx.arc( xo, yo, 4, 0, Math.PI*2);
    ctx.stroke();

}

function updateResultImage() {

    if ( fullResult == null ) {
	return;
    }

    var minMax = fullResult.getMinMax();
    var scal = 255./(minMax[1]-minMax[0]);
    var min = minMax[0];

    var imgCnv = document.getElementById("resultCanvas");
    var ctx = imgCnv.getContext("2d");
    var fftData = ctx.getImageData(0,0,imgCnv.width, imgCnv.height);
    //var resData  = fullResult.pwSpec();
    var resData  = fullResult.data;

    //logger( pwSpec.length);


    var data = fftData.data;
    for ( var y = 0 ; y<fullResult.size; y++) {
        for ( var x = 0 ; x<fullResult.size; x++) {
	    var io  = x + y * imgCnv.width;
	    var ii  = x + y * fullResult.size;

	    data[io*4+0] = (resData[2*ii]-min)*scal ;
	    data[io*4+1] = (resData[2*ii]-min)*scal ;
	    data[io*4+2] = (resData[2*ii]-min)*scal ;
	    /*
	    data[io*4+0] = resData[ii]*255;
	    data[io*4+1] = resData[ii]*255;
	    data[io*4+2] = resData[ii]*255;
	    */
	    data[io*4+3] = 0xFF;
	}
    } 

    ctx.putImageData( fftData,0,0);

}
