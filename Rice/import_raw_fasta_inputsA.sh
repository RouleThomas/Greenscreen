#!/bin/bash

module load sratoolkit/2.11.2


#inputA
fasterq-dump SRR10423462 -S
#inputB
fasterq-dump SRR10423460 -S
#inputC
fasterq-dump SRR10423463 -S
#inputD
fasterq-dump SRR10423461 -S