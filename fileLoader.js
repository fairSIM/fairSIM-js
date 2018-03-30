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
    	//document.getElementById("status").innerHTML=Number.parseFloat(loaded*100).toFixed(2);
	showStatus( "Loading TIFF ",Number.parseFloat(loaded*90).toFixed(2));
    }
  }
}

// file has loaded
function loaded(evt) {

    //document.getElementById("status").innerHTML="100% done, decoding TIFF";	
    showStatus("done, decoding TIFF",.9);    

    // decode tiff
    //tiffPages = decode( reader.result, false, true );
    tiffPages = decode( reader.result, false, true );


    setImage(0);   
    document.getElementById("zSlider").max = (tiffPages.length/15)-1;
    document.getElementById("sSlider").max = 14;
    //document.getElementById("status").innerHTML="100% done, TIFF decoded";	
    
    showStatus("loaded tiff: "+tiffPages.length+" slices @"
	+tiffPages[0].width+"x"+tiffPages[0].height,1);
}

function errorHandler(evt) {
  if(evt.target.error.name == "NotReadableError") {
	showStatus("File load error",-1);  
    }
}


function downloadExample() {

    var xhr =  new XMLHttpRequest();
    xhr.open('GET', './examples.tif', true);
    xhr.responseType = 'arraybuffer';
    xhr.onreadystatechange = function() {
	if (xhr.readyState == 4 ) {
	    tiffPages = decode( new Uint8Array(this.response), false, true );


	    setImage(0);   
	    document.getElementById("zSlider").max = (tiffPages.length/15)-1;
	    document.getElementById("sSlider").max = 14;
	    //document.getElementById("status").innerHTML="100% done, TIFF decoded";	
	    
	    showStatus("Dowloaded tiff: "+tiffPages.length+" slices @"
		+tiffPages[0].width+"x"+tiffPages[0].height,1);

	}
    }
    xhr.send(null);
}


