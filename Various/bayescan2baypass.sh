#!/bin/bash

## 17/02/22 Javier Sala
## To run the script: sh bayescan2baypass.sh *name_of_bayescan_file*


awk '/^\[pop\]=/{out=substr($0, 7) "_tmp_baypass.txt"} {if(out) print $4, $5 > out}' $1
paste $(ls -v *_tmp_baypass.txt) | sed 's/\t/ /g' > baypass_infile.txt
printf "\nThe order in which 'paste' pasted each file to baypass_infile is:\n"
printf "$(ls -v *_tmp_baypass.txt) \n\n"
rm -f *_tmp_baypass*
