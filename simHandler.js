function logger(text) {

    console.log("fairSIM-js: "+text);
    document.getElementById("fairsimlogger").innerHTML = "<pre>"+text+"<pre>";
}


const maxAng = 3;
const maxPha = 5;

var tiffPages   = null;
var inputFFTimg = null;
var zImported   = -1;

function setImage(pos) {

    if (tiffPages==null) {
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


function importImages() {
    if (tiffPages==null) {
	return;
    }
    
    var z =  document.getElementById("zSlider").value; 
    zImported = z;   
 
    inputFFTimg = [];

    for (var ang = 0 ; ang < maxAng ; ang++ ) {
    for (var pha = 0 ; pha < maxPha ; pha++ ) {

	var pos   = z*maxPha + pha + ang*(tiffPages.length/maxAng);
	var fftIn = new Vec2dCplx( tiffPages[ pos  ].width );
	fftIn.setFromRGB( tiffPages[ pos ].data );
	fftIn.fft2d( true );
	inputFFTimg.push( fftIn); 
    }}

    updateFFTimage(document.getElementById("sSlider").value);
    logger("Input complete, z: "+z);
}

function updateFFTimage( pos ) {

    if ( inputFFTimg == null || inputFFTimg.length == 0 ) {
	return;
    }

    var imgCnv = document.getElementById("fftCanvas");
    var ctx = imgCnv.getContext("2d");
    var fftData = ctx.getImageData(0,0,imgCnv.width, imgCnv.height);
    var pwSpec  = inputFFTimg[pos].pwSpec();

    logger( pwSpec.length);

    var data = fftData.data;
    for ( var i = 0 ; i<data.length/4; i++) {
	data[i*4+0] = pwSpec[i]*255 ;
	data[i*4+1] = pwSpec[i]*255 ;
	data[i*4+2] = pwSpec[i]*255 ;
	data[i*4+3] = 0xFF;
    } 
    ctx.putImageData( fftData,0,0);

}


