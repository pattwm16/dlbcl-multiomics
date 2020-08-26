# Python script to convert a csv of Ensembl IDs to Gene Symbols

import mygene
import pandas
import numpy as np

mg = mygene.MyGeneInfo()

# read in csv and extract list of Ensembl ids
target = input("File path to csv containing genes: ")
genes = pandas.read_csv(str(target), header=0)
ens_ids = genes.x.tolist()

# query Ensembl REST API and return gene symbols
g_symbols = mg.querymany(ens_ids, scopes='ensembl.gene')

# check if any are members in published genes for TCGA
pg = pandas.read_csv('publishedGenes.csv', header=0)
pubGenes = pg.GeneSymbol.tolist()

# write these symbols to a .csv
df = pandas.DataFrame(g_symbols)
df["memberTCGA"] = np.where(df["symbol"].isin(pubGenes), "True", "False")

df.to_csv("./ProblemGeneSymbols.csv", sep=',', index=False,
            columns=['query', 'symbol', 'name', 'memberTCGA'])
