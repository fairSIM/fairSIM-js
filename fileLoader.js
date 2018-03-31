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


// monitors progress when loading files from dist
function updateProgress(evt) {
  if (evt.lengthComputable) {
    
    var loaded = (evt.loaded / evt.total);
    if (loaded <= 1) {
	showStatus( "Loading TIFF ", loaded);
    }
  }
}

// file has loaded, convert the tiff
function loaded(evt) {

    showStatus("Loaded, now decoding TIFF...",.99);    

    // decode tiff
    tiffPages = decode( reader.result, false, true );

    setImage(0);   
    document.getElementById("zSlider").max = (tiffPages.length/15)-1;
    document.getElementById("sSlider").max = 14;
    document.getElementById("zSlider").value = 0;
    document.getElementById("sSlider").value = 0;
    
    showStatus("Loaded & decoded TIFF: "+tiffPages.length+" slices @"
	+tiffPages[0].width+"x"+tiffPages[0].height,1);
}

// handle file load errors
function errorHandler(evt) {
  if(evt.target.error.name == "NotReadableError") {
	showStatus("File load error",-1);  
    }
}


// handle downloading the example file
function downloadExample() {

    var xhr =  new XMLHttpRequest();
    xhr.open('GET', './examples.tif', true);

    xhr.responseType = 'arraybuffer';
    xhr.onprogress = updateProgress;


    xhr.onreadystatechange = function() {
	if (xhr.readyState == 4 ) {
	    tiffPages = decode( new Uint8Array(this.response), false, true );

	    setImage(0);   
	    document.getElementById("zSlider").max = (tiffPages.length/15)-1;
	    document.getElementById("sSlider").max = 14;
	    document.getElementById("zSlider").value = 0;
	    document.getElementById("sSlider").value = 0;
	    
	    showStatus("Dowloaded & decoded tiff: "+tiffPages.length+" slices @"
		+tiffPages[0].width+"x"+tiffPages[0].height,1);

	}
    }
    xhr.send(null);
}


