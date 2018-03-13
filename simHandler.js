function logger(text) {

    console.log("fairSIM-js: "+text);
    document.getElementById("fairsimlogger").innerHTML = "<pre>"+text+"<pre>";
}


const maxAng = 3;
const maxPha = 5;

var objNA=1.4;
var emLambda=525;

var tiffPages   = null;
var inputFFTimg = null;
var bandImg	= null;
var corrImg	= null;
var zImported   = -1;

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
function computeCorrImages() {
    
    const bands = (maxPha+1)/2;

    if (bandImg == null || bandImg.length != bands*maxAng ) {
	logger("import some images first..." );
	return;
    }

    corrImg = [];

    for (var ang = 0 ; ang < maxAng ; ang++ ) {
	for ( var b=1; b<bands; b++) {

	    var c = bandImg[ (b +ang*bands ) ].copy();
	    c.times( bandImg[  ang*bands ]);
	    c.fft2d();
	    corrImg.push(c);
	}
    }    

    updateCorrelationImage(document.getElementById("corrSlider").value);
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
    ctx.putImageData( fftData,0,0);

}


