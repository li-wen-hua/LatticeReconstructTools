#!/bin/bash

echo '------------ Filter defects ------------'

cd AtomsNearDisline_f/
cp ../datafiles-input/disline.input ./
cp ../datafiles-input/ws.dump ./
bash Select.sh
cp awaydislocations.dump.output ../datafiles-output/ws-new.dump

echo 'DONE!'









