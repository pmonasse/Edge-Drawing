#!/bin/bash
# Arguments:
# -gradient: min gradient 
# -gap: anchor gap
# -length: min length
# -sigma: Gaussian blur
# -valid: a conrario validation (true/false)
# -leps: log10(epsilon) validation threshold (optional)
# -sub: sublines (0,1,2)
# -tjunc: compute T-junction (true/false)

grad=$1
gap=$2
length=$3
sigma=$4
valid="${5/false/}"
if [ "$5" = "true" ]; then
    valid="-e ${6:-0}" # default log10(eps)=0
    valid="$valid -S ${7:-0}";
fi

if [ "$8" = "true" ]; then
    tjunctions="-t tjunc.png";
fi

echo $bin/build/edgeDrawing $valid -g $grad -a $gap -l $length -s $sigma $tjunctions $input_0 edges.png
