#!/usr/bin/python3

import json
import pandas as pd

## JSON created directly from the GAMs
## vg view -a x.gam >x.json

ref_json = 'ref_NA12878_only_chr21.json' 
graph_json = 'graph_NA12878_only_chr21.json'

def json_readr(file):
	for line in open(file, mode="r"):
		yield json.loads(line)


## https://stackoverflow.com/questions/27907633/multiple-json-objects-in-one-file-extract-by-python
## "...Python allows you to put a generator function inside of a list and populate the list automatically..."

lst = list(json_readr(ref_json))

scores_list = []
identity_list = []
mapping_quality = []

for i in range(len(lst)):
	if 'score' in lst[i]:
		scores_list.append(lst[i]['score'])
	if 'identity' in lst[i]:
		identity_list.append(lst[i]['identity'])
	if 'mapping_quality' in lst[i]:
		mapping_quality.append(lst[i]['mapping_quality'])

## these lists may have different lengths...



mapp = pd.Series(mapping_quality)
## mapp.to_csv("ref_NA12878_only_chr21_mapping_quality.csv", index=False)
mapp.to_csv("graph_NA12878_only_chr21_mapping_quality.csv", index=False)


df = pd.DataFrame({'scores': scores_list, 'indentity': identity_list})
## df.to_csv("ref_NA12878_only_chr21_scores_identity.csv", index=False)
df.to_csv("graph_NA12878_only_chr21_scores_identity.csv", index=False)



## REFERENCE, 'ref_NA12878_only_chr21.json' 
## 
## >>> len(lst)
## 13443046
## >>> len(scores_list)
## 13410334
## >>> len(identity_list)
## 13410334
## >>> len(mapping_quality)
## 12913935
## 


## GRAPH, 'graph_NA12878_only_chr21.json'
## >>> len(lst)
## 13403870
## >>> len(scores_list)
## 13403870
## >>> len(identity_list)
## 13403870
## >>> len(mapping_quality)
## 12659777
##  






