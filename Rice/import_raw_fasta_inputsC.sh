#!/bin/bash

module load sratoolkit/2.11.2

#inputJ
fasterq-dump SRR8746746 -S
#inputK
fasterq-dump SRR8746747 -S
#inputL
fasterq-dump SRR8746751 -S
#inputM
fasterq-dump SRR8746752 -S