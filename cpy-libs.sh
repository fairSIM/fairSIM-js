#!/bin/bash
sed '/module/ {N;N;s/module.exports.*}\;//}' decode-tiff/src/index.js > tiff-decoder.js
