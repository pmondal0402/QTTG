#!/bin/bash
nohup matlab -nosplash -nodisplay -nodesktop -r 'try; main_get_entropy_elec;catch; end; quit' > output2.log &
