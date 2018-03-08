#!/bin/bash
sed '/module/ {N;N;s/module.exports.*}\;//}' decode-tiff/src/index.js > tiff-decoder.js
cp fft.js/lib/complex.js  fft-runner.js


