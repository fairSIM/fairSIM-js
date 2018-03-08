function logger(text) {

    console.log("fairSIM-js: "+text);
    document.getElementById("fairsimlogger").innerHTML = "<pre>"+text+"<pre>";
}

var tiffPages = null;

var reader = new FileReader();

function startRead() {
  // obtain input element through DOM

  var file = document.getElementById('file').files[0];
  if(file){
    getAsBuffer(file);
  }
}

function getAsBuffer(readFile) {


  // Read file into memory as UTF-16
  reader.readAsArrayBuffer(readFile);

  // Handle progress, success, and errors
  reader.onprogress = updateProgress;
  reader.onload = loaded;
  reader.onerror = errorHandler;
}

function updateProgress(evt) {
  if (evt.lengthComputable) {
    // evt.loaded and evt.total are ProgressEvent properties
    var loaded = (evt.loaded / evt.total);
    if (loaded <= 1) {
      // Increase the prog bar length
    	document.getElementById("status").innerHTML=Number.parseFloat(loaded*100).toFixed(2);
    }
  }
}

// file has loaded
function loaded(evt) {

    document.getElementById("status").innerHTML="100% done, decoding TIFF";	
    
    // decode tiff
    tiffPages = decode( reader.result, false );

    logger("loaded tiff: "+tiffPages.length+" slices @"
	+tiffPages[0].width+"x"+tiffPages[0].height);

    setImage(0);   
    document.getElementById("zSlider").max = (tiffPages.length/15)-1;
    document.getElementById("sSlider").max = 14;
    document.getElementById("status").innerHTML="100% done, TIFF decoded";	

    
}

function errorHandler(evt) {
  if(evt.target.error.name == "NotReadableError") {
    // The file could not be read
  }
}


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

    for ( var i = 0 ; i<data.length; i++) {
	data[i] = tiffPages[ pos*5 + pha + ang*(tiffPages.length/3)  ].data[i];
    }
    ctx.putImageData( imgData,0,0);

}

