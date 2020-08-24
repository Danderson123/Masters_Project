#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 13:24:09 2020

@author: danielanderson
"""


import networkx as nx
import matplotlib.pyplot as plt

g = nx.read_gml("final_graph.gml")

unique = []
desc = []
mem = []
for node in g._node:
    y = g._node[node]
    if isinstance(y["members"], int):
        unique.append(y["name"])
        desc.append(y["description"])
        mem.append(y["members"])

sp = []
sp_d = []
for x in range(len(mem)):
    if mem[x] > 27:
        sp.append(mem[x])
        sp_d.append(desc[x])
        
core_cogs = []
annot_core_cogs = []

named_cogs = []
tot_cogs = []
named_genes = []
tot_genes = []
functionally_annotated = []
annot_mem = []
for node in g._node:
    y = g._node[node]
    tot_cogs.append(y["name"])
    if isinstance(y["members"], int):
        tot_genes.append(1)
    elif isinstance(y["members"], list):
        tot_genes.append(len(y["members"]))
    if not "group_" in y["name"]:
        named_cogs.append(y["name"])
        if isinstance(y["members"], int):
            named_genes.append(1)
        else:
            named_genes.append(len(y["members"]))
            if len(y["members"]) > 637:
               core_cogs.append(y["members"])
    if not y["description"] == "hypothetical protein" and not y["description"] == "" and not y["description"] == " ":
        functionally_annotated.append(y["description"])
        if isinstance(y["members"], int):
            annot_mem.append(1)
        else:
            annot_mem.append(len(y["members"]))
            if len(y["members"]) > 637:
               annot_core_cogs.append(y["members"])
               
named_genes = sum(named_genes)
annot_mem = sum(annot_mem)
tot_genes = sum(tot_genes)

named_len = []
for x in core_cogs:
    named_len.append(len(x))
named_len=sum(named_len)

leng = []
for x in annot_core_cogs:
    leng.append(len(x))
leng=sum(leng)


c = nx.read_gml("crouch_graph.gml")

c_unique = []
for node in c._node:
    y = c._node[node]
    if isinstance(y["members"], int):
        c_unique.append(y["name"])
        
c_core_cogs =[]
c_annot_core_cogs = []

c_named_cogs = []
c_tot_cogs = []
c_named_genes = []
c_tot_genes = []
c_functionally_annotated = []
c_annot_mem = []

for node in c._node:
    y = c._node[node]
    c_tot_cogs.append(y["name"])
    if isinstance(y["members"], int):
        c_tot_genes.append(1)
    elif isinstance(y["members"], list):
        c_tot_genes.append(len(y["members"]))
    if not "group_" in y["name"]:
        c_named_cogs.append(y["name"])
        if isinstance(y["members"], int):
            c_named_genes.append(1)
        else:
            c_named_genes.append(len(y["members"]))
            if len(y["members"]) > 637:
               c_core_cogs.append(y["members"])
    if not y["description"] == "hypothetical protein" and not y["description"] == "" and not y["description"] == " ":
            c_functionally_annotated.append(y["description"])
            if isinstance(y["members"], int):
                c_annot_mem.append(1)
            else:
                c_annot_mem.append(len(y["members"]))
                if len(y["members"]) > 637:
                    c_annot_core_cogs.append(y["members"])
c_named_genes = sum(c_named_genes)
c_annot_mem = sum(c_annot_mem)
c_tot_genes = sum(c_tot_genes)

c_named_len = []
for x in c_core_cogs:
    c_named_len.append(len(x))
c_named_len=sum(c_named_len)

c_len = []
for x in c_annot_core_cogs:
    c_len.append(len(x))
c_len=sum(c_len)

x = ["Reference annotated", "Manually annotated"]
#named
height_cogs = [471/18985*100, 491/4833*100]
#height_cogs_core = [394/18985*100, 390/4833*100]
height_cogs_core = [471/18985*100, 491/4833*100]
height_cogs_accessory = [77/18985*100, 111/4833*100]
height_genes = [292645500/1184477453*100, 292099748/1116592720*100]
#height_genes_core = [253677/1357236*100, 251080/1299989*100]
height_genes_core = [272887313/1184477453*100, 263962370/1116592720*100]
height_genes_accessory = [19758187/1184477453*100, 28137378/1116592720*100]

fig, ax = plt.subplots()
ax.bar(x, height_cogs_core, width=0.8, bottom=None, align='center', color= "deepskyblue")
ax.bar(x, height_cogs_accessory, width=0.8, bottom=None, align='center', color= "C1")
plt.ylim((0,100))
fig.savefig('named_cogs.pdf')

fig, ax = plt.subplots()
ax.bar(x, height_genes_core, width=0.8, bottom=None, align='center', color= "deepskyblue")
ax.bar(x, height_genes_accessory, width=0.8, bottom=None, align='center', color= "C1")
plt.ylim((0,100))
fig.savefig('named_sequence.pdf')

x = ["Reference annotated", "Manually annotated"]
#annotated
height_cogs = [3203/18985*100, 3663/4833*100]
height_cogs_core = [3203/18985*100, 3663/4833*100]
#height_cogs_accessory = [1742/18985*100, 2318/4833*100]
height_cogs_accessory = [1742/18985*100, 2318/4833*100]
height_genes = [1122359362/1184477453*100, 1066334521/1116592720*100]
#height_genes_core = [940440/1357236*100, 865587/1299989*100]
height_genes_core = [1122359362/1184477453*100, 1066334521/1116592720*100]
height_genes_accessory = [229541059/1184477453*100, 251065441/1116592720*100]

fig, ax = plt.subplots()
plt.bar(x, height_cogs_core, width=0.8, bottom=None, align='center', color= "deepskyblue")
plt.bar(x, height_cogs_accessory, width=0.8, bottom=None, align='center', color= "C1")
plt.ylim((0,100))
fig.savefig('annotated_cogs.pdf')

fig, ax = plt.subplots()
plt.bar(x, height_genes_core, width=0.8, bottom=None, align='center', color= "deepskyblue")
plt.bar(x, height_genes_accessory, width=0.8, bottom=None, align='center', color= "C1")
plt.ylim((0,100))
fig.savefig('annotated_sequence.pdf')

false_prop = (1122359362-14778246)/1184477453*100






import pandas as pd
import matplotlib.pyplot as plt
# Create a dataframe
value1= [100,100,0.2,7,16.6,11.4,99.5,98.4,19.0]
value2= [100,100,0.2,7.5,16.4,11.4,35.6,62.5,19.0]
df = pd.DataFrame({'Gene': ["pbp2a", 'pbp1a', "aph3IIIa",'ermB','mefA','tetM','pspA','pspC','zmpC'], 'Reference annotated':value1 , 'Manually annotated':value2 })

fig, ax = plt.subplots()

ordered_df = df.sort_values(by='Manually annotated')
my_range=range(1,len(df.index)+1)

import seaborn as sns
sns.set_style("whitegrid")
plt.hlines(y=my_range, xmin=ordered_df['Reference annotated'], xmax=ordered_df['Manually annotated'], color='grey', alpha=0.7)
plt.scatter(ordered_df['Manually annotated'], my_range, color='green' , label='Manually annotated')
plt.scatter(ordered_df['Reference annotated'], my_range, color='darkviolet',  label='Reference annotated')
plt.legend(loc = "lower right")
plt.xlim((-5,105))
# Add title and axis names
plt.yticks(fontsize= 10)
plt.yticks(my_range, ordered_df['Gene'])
#plt.title("Comparison of the value 1 and the value 2", loc='left')
ax.set_xlabel('Estimated frequency in the SPARC dataset (%)')
ax.set_ylabel('Gene')

fig.savefig('refvscrouchplot.pdf', bbox_inches='tight')


import pandas as pd
import matplotlib.pyplot as plt
# Create a dataframe
value1= [99.7,100,30.2,90,20.3,16.6,11.7,0,0.2]
value2= [100,100,35.6,62.5,19,16.4,11.4,7.5,0.2]
df = pd.DataFrame({'Gene': ["pbp2a", 'pbp1a',"pspA",'pspC','zmpC','mefA','tetM','ermB','aph3IIIa'], 'Annotated by consensus':value1 , 'Manually annotated':value2 })

fig, ax = plt.subplots()

ordered_df = df.sort_values(by='Manually annotated')

#ordered_df = df.sort_values(by='Manually annotated')
my_range=range(1,len(df.index)+1)

import seaborn as sns
sns.set_style("whitegrid")
plt.hlines(y=my_range, xmin=ordered_df['Annotated by consensus'], xmax=ordered_df['Manually annotated'], color='grey', alpha=0.7)
plt.scatter(ordered_df['Manually annotated'], my_range, color='green' , label='Manually annotated')
plt.scatter(ordered_df['Annotated by consensus'], my_range, color='plum',  label='Annotated by COBS')
plt.legend(loc = "lower right")
plt.xlim((-5,105))

# Add title and axis names
plt.yticks(fontsize= 10)
plt.yticks(my_range, ordered_df['Gene'])
#plt.title("Comparison of the value 1 and the value 2", loc='left')
ax.set_xlabel('Estimated frequency in the SPARC dataset (%)')
ax.set_ylabel('Gene')

fig.savefig('crouchvscobs.pdf', bbox_inches='tight')



import pandas as pd
import matplotlib.pyplot as plt
# Create a dataframe
value1= [78.7,96.6,20,42.9,14.6,13.1,9.7,0,0.2]
value2= [100,100,35.6,62.5,19,16.4,11.4,7.5,0.2]
df = pd.DataFrame({'Gene': ["pbp2a", 'pbp1a',"pspA",'pspC','zmpC','mefA','tetM','ermB','aph3IIIa'], 'Annotated by consensus':value1 , 'Manually annotated':value2 })

fig, ax = plt.subplots()

ordered_df = df.sort_values(by='Manually annotated')
my_range=range(1,len(df.index)+1)

import seaborn as sns
sns.set_style("whitegrid")
plt.hlines(y=my_range, xmin=ordered_df['Annotated by consensus'], xmax=ordered_df['Manually annotated'], color='grey', alpha=0.7)
plt.scatter(ordered_df['Manually annotated'], my_range, color='green' , label='Manually annotated')
plt.scatter(ordered_df['Annotated by consensus'], my_range, color='cornflowerblue',  label='Annotated by Panaroo-COBS consensus')
plt.legend(loc = "lower right")
plt.xlim((-5,105))

# Add title and axis names
plt.yticks(fontsize= 10)
plt.yticks(my_range, ordered_df['Gene'])
#plt.title("Comparison of the value 1 and the value 2", loc='left')
ax.set_xlabel('Estimated frequency in the SPARC dataset (%)')
ax.set_ylabel('Gene')

fig.savefig('crouchvsconsensus.pdf', bbox_inches='tight')





