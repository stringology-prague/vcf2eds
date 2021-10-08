## How to Generate Synthetic Data

### Length scaling datasets (baseline)

Generate DNA dataset with total length scaling

```bash
while read LEN; do python3 random_eds.py --length $LEN --alphabet D --deg-prob 0.1 --segment-size-avg 5 --segment-size-max 10 --element-len-avg 5 --element-len-max 10 --decorate-output "synth-dna"; done < lengths.list
```

Generate protein dataset with total length scaling

```bash
while read LEN; do python3 random_eds.py --length $LEN --freq-file protein.freq --deg-prob 0.1 --segment-size-avg 5 --segment-size-max 10 --element-len-avg 5 --element-len-max 10 --decorate-output "synth-protein"; done < lengths.list 
```

### Degenerate probability scaling datasets

Generate DNA dataset with degenerate symbol probability scaling 

```bash
while read PROB; do python3 random_eds.py --length 1600000 --alphabet D --deg-prob $PROB --segment-size-avg 5 --segment-size-max 10 --element-len-avg 5 --element-len-max 10 --decorate-output "synth-prob"; done < deg-prob.list
```

### Segment size scaling datasets

Generate DNA dataset with average segment size scaling, the overall LENGTH is fixed (and set to a smaller amount to compensate for larger file sizes)

```bash
while read SEG_SIZE; do python3 random_eds.py --length 100000 --alphabet D --deg-prob 0.1 --segment-size-avg $SEG_SIZE --element-len-avg 5 --element-len-max 10 --decorate-output "synth-seg-size-avg"; done < seg-size-avg.list
```

Generate DNA dataset with average segment size scaling, and degenerate probability scaling, the overall SIZE is fixed

```bash
while read PROB; do while read SEG_SIZE; do python3 random_eds.py --size 3200000 --alphabet D --deg-prob $PROB --segment-size-avg $SEG_SIZE --element-len-avg 5 --element-len-max 10 --decorate-output "synth-seg-size-avg-fix-size"; done < seg-size-avg.list; done < seg-size-avg.deg-prob.list
```

### Element length scaling datasets

Generate DNA dataset with average element length scaling and for different degenerate symbol probabilities, the overall LENGTH is fixed

```bash
while read PROB; do while read ELM_LEN; do python3 random_eds.py --length 1600000 --alphabet D --deg-prob $PROB --segment-size-avg 5 --segment-size-max 10 --element-len-avg $ELM_LEN --decorate-output "synth-elm-len-avg"; done < elm-len-avg.list; done < elm-len-avg.deg-prob.list 
```

Generate DNA dataset with average element length scaling and for different degenerate symbol probabilities, the overall SIZE is fixed

```bash
while read PROB; do while read ELM_LEN; do python3 random_eds.py --size 3200000 --alphabet D --deg-prob $PROB --segment-size-avg 5 --segment-size-max 10 --element-len-avg $ELM_LEN --decorate-output "synth-elm-len-avg-fix-size"; done < elm-len-avg.list; done < elm-len-avg.deg-prob.list 
```

### Special cases

Generate simple DNA string (no degenerate segments) - for testing BNDM only

```bash
while read LEN; do python3 random_eds.py --length $LEN --deg-prob 0.0 --segment-size-avg 0 --segment-size-max 0 --element-len-avg 0 --element-len-max 0 --decorate-output "synth-nodeg"; done < lengths.list
```


## Conver existing data to a format with --size

```bash
echo "bname,newname,size" > /tmp/eds.sizes.csv; find . -name '*.eds' -print0 | while read -d $'\0' FNAME; do SIZE=$(cat $FNAME | tr -d "{}," | wc -c); BNAME=$(basename $FNAME); NEW_FNAME=$(echo $BNAME | sed -n "s/\(.*\).eds/\1_N=$(printf "%010d" $SIZE).eds/p"); mv $FNAME $NEW_FNAME; gsutil mv "gs://prg-str-genome-data/eds/$BNAME" "gs://prg-str-genome-data/eds/$NEW_FNAME" || gsutil cp $NEW_FNAME "gs://prg-str-genome-data/eds/$NEW_FNAME"; done
```