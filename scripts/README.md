## How to Generate Synthetic Data

Generate DNA dataset with total length scaling

```bash
while read LEN; do python3 random_eds.py --length $LEN --alphabet D --deg-prob 0.1 --segment-size-avg 5 --segment-size-max 10 --element-len-avg 5 --element-len-max 10 --decorate-output "synth-dna"; done < lengths.list
```

Generate protein dataset with total length scaling

```bash
while read LEN; do python3 random_eds.py --length $LEN --freq-file protein.freq --deg-prob 0.1 --segment-size-avg 5 --segment-size-max 10 --element-len-avg 5 --element-len-max 10 --decorate-output "synth-protein"; done < lengths.list 
```

Generate DNA dataset with degenerate symbol probability scaling 

```bash
while read PROB; do python3 random_eds.py --length 1600000 --alphabet D --deg-prob $PROB --segment-size-avg 5 --segment-size-max 10 --element-len-avg 5 --element-len-max 10 --decorate-output "synth-prob"; done < deg-prob.list
```

Generate DNA dataset with average segment size scaling (the overall length here is set to a smaller amount to compensate for larger file sizes)

```bash
while read SEG_SIZE; do python3 random_eds.py --length 100000 --alphabet D --deg-prob 0.1 --segment-size-avg $SEG_SIZE --element-len-avg 5 --element-len-max 10 --decorate-output "synth-seg-size-avg"; done < seg-size-avg.list
```

Generate DNA dataset with average element length scaling 

```bash
while read ELM_LEN; do python3 random_eds.py --length 1600000 --alphabet D --deg-prob 0.1 --segment-size-avg 5 --segment-size-max 10 --element-len-avg $ELM_LEN --decorate-output "synth-elm-len-avg"; done < elm-len-avg.list
```
