#!/bin/sh
#
# test.sh

return_code=0

for file in tests/physiboss/physiboss-cell-lines/*; 
do 
    if [ $file != "tests/physiboss/physiboss-cell-lines/test.sh" ]; then
        diff --strip-trailing-cr output/$(basename $file) $file; 
        if [ $? != 0 ]; then 
            return_code=1
        fi;
    fi;
done

exit $return_code