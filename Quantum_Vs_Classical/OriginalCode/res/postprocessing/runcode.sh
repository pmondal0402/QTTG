#!/bin/bash
nohup matlab -nosplash -nodisplay -nodesktop -r 'try; main; catch; end; quit' > output.log &
