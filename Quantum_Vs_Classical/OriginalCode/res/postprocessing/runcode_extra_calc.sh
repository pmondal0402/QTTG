#!/bin/bash
nohup matlab -nosplash -nodisplay -nodesktop -r 'try; main_extra_calc; catch; end; quit' > output.log &
