#!/bin/bash

###############################################################################################################
## Define paths to software and reference files

RWD_MAIN=/Users/gosia/Dropbox/UZH/carsten_cytof
RCODE=/Users/gosia/Dropbox/UZH/carsten_cytof_code
METADATA=$RWD_MAIN/CK_metadata
PANELS=$RWD_MAIN/CK_panels


## Define which analysis to re-run
data_normalization=false
pcascores=false
select_observables=false
flowsom=false
flowsom_validation=false
heatmaps=true
runtsne=false
plottsne=false
plottsne_expr=false
frequencies=false
expression=false
cluster_merging=false
cluster_extracting=false
fcs_saving=false
cytokines_bimatrix=false
cytokines_bimatrix_main=false
cytokines_bimatrix_cluster_merging=false
pd1_bimatrix=false
pd1_bimatrix_main=false
pd1_bimatrix_cluster_merging=false
cytokines_expression=false
pd1_expression=false
cd69_bimatrix=false
cd69_bimatrix_main=false

## global parameters
tsne_pmin=1500

###############################################################################################################
# Analysis of CK_2016-06-29_01 data
# Use Analysis block 1
###############################################################################################################

DATA=29
PANEL=1
data_dir="CK_2016-06-29_01"

file_panel="panel1.xlsx"
file_metadata="metadata_29_01.xlsx"

pca_score_cutoff=3.5
rand_seed_consensus=1234
nmetaclusts=20

prefix_data="29_"
prefix_panel="01_"
prefix_pca="pca1_"
prefix_clust="cl20_"



# --------------------------------------------------
# Analysis of CK_2016-06-29_01 merging3
# Use Analysis block 2
# --------------------------------------------------

prefix_clust="cl20_"

file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging2.xlsx"
prefix_merging="merging2_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}


file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging3.xlsx"
prefix_merging="merging3_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}


file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging4.xlsx"
prefix_merging="merging4_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}


# --------------------------------------------------
# CK_2016-06-29_01 - CD4 and CD8 cluster extracting from merging2
# Use Analysis block 3
# --------------------------------------------------

prefix_merging="merging2_"

extract_cluster=("'CD4'" "'CD8'")
extract_dir=('CK_2016-06-29_01_CD4_merging2' 'CK_2016-06-29_01_CD8_merging2')

for i in 0 1
do
  ./Analysis_block_3_cluster_extracting.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_extracting ${cluster_extracting} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging} --extract_cluster ${extract_cluster[$i]} --extract_dir ${extract_dir[$i]}
done

# --------------------------------------------------
# Analysis of CK_2016-06-29_01_CD4_merging2 using panel1CD4.xlsx
# and CK_2016-06-29_01_CD8_merging2 using panel1CD8.xlsx
# Use Analysis block 1
# --------------------------------------------------

DATA=29
PANEL=1

file_metadata="metadata_29_01.xlsx"

rand_seed_consensus=1234
nmetaclusts=20

prefix_pca="pca1_"
prefix_clust="cl20_"

data_dir=('CK_2016-06-29_01_CD4_merging2' 'CK_2016-06-29_01_CD8_merging2')
prefix_data=('29CD4_' '29CD8_')
file_panel=('panel1CD4.xlsx' 'panel1CD8.xlsx')
prefix_panel=('01CD4_' '01CD8_')
pca_score_cutoff=(1.9 1.9)


# --------------------------------------------------
# Analysis of CK_2016-06-29_01_CD4_merging2 using panel1CD4.xlsx of cluster_merging2
# and CK_2016-06-29_01_CD8_merging2 using panel1CD8.xlsx of cluster_merging2
# Use Analysis block 2
# --------------------------------------------------

prefix_clust="cl20_"

file_merging=("${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_clust}cluster_merging2.xlsx" "${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_clust}cluster_merging2.xlsx")
prefix_merging=('merging2_' 'merging2_')

for i in 0 1
do
  ./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging[$i]} --file_merging ${file_merging[$i]}
done


prefix_clust=('cl8_' 'cl5_')

file_merging=("${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_clust[0]}cluster_merging4.xlsx" "${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_clust[1]}cluster_merging3.xlsx")
prefix_merging=('merging4_' 'merging3_')

for i in 0 1
do
  ./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust[$i]} --prefix_merging ${prefix_merging[$i]} --file_merging ${file_merging[$i]}
done


prefix_clust=('cl5_' 'cl8_')

file_merging=("${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_clust[0]}cluster_merging5.xlsx" "${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_clust[1]}cluster_merging4.xlsx")
prefix_merging=('merging5_' 'merging4_')

for i in 0 1
do
  ./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust[$i]} --prefix_merging ${prefix_merging[$i]} --file_merging ${file_merging[$i]}
done

###############################################################################################################
# Analysis of CK_2016-06-23_01 data
# Use Analysis block 1
###############################################################################################################

DATA=23
PANEL=1
data_dir="CK_2016-06-23_01"

file_panel="panel1.xlsx"
file_metadata="metadata_23_01.xlsx"

pca_score_cutoff=3
rand_seed_consensus=123
nmetaclusts=20

prefix_data="23_"
prefix_panel="01_"
prefix_pca="pca1_"
prefix_clust="cl20_"


# --------------------------------------------------
# Analysis of CK_2016-06-23_01 cluster_merging
# Use Analysis block 2
# --------------------------------------------------

prefix_clust="cl20_"

file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_mergingNEW2.xlsx"
prefix_merging="mergingNEW2_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}


file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging5.xlsx"
prefix_merging="merging5_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}


file_merging="${prefix_data}${prefix_panel}${prefix_pca}${prefix_clust}cluster_merging6.xlsx"
prefix_merging="merging6_"

./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging} --file_merging ${file_merging}


# --------------------------------------------------
# CK_2016-06-23_01 - CD4 and CD8 cluster extracting from mergingNEW2
# Use Analysis block 3
# --------------------------------------------------

prefix_merging="mergingNEW2_"

extract_cluster=("'CD4'" "'CD8'")
extract_dir=('CK_2016-06-23_01_CD4_mergingNEW2' 'CK_2016-06-23_01_CD8_mergingNEW2')

for i in 0 1
do
  ./Analysis_block_3_cluster_extracting.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir} --cluster_extracting ${cluster_extracting} --METADATA ${METADATA} --file_metadata ${file_metadata} --prefix_data ${prefix_data} --prefix_panel ${prefix_panel} --prefix_pca ${prefix_pca} --prefix_merging ${prefix_merging} --extract_cluster ${extract_cluster[$i]} --extract_dir ${extract_dir[$i]}
done


# --------------------------------------------------
# Analysis of CK_2016-06-23_01_CD4_mergingNEW2 using panel1CD4.xlsx
# and CK_2016-06-23_01_CD8_mergingNEW2 using panel1CD8.xlsx
# Use Analysis block 1
# --------------------------------------------------

DATA=23
PANEL=1

file_metadata="metadata_23_01.xlsx"

rand_seed_consensus=1234
nmetaclusts=20

prefix_pca="pca1_"
prefix_clust="cl20_"

data_dir=('CK_2016-06-23_01_CD4_mergingNEW2' 'CK_2016-06-23_01_CD8_mergingNEW2')
prefix_data=('23CD4_' '23CD8_')
file_panel=('panel1CD4.xlsx' 'panel1CD8.xlsx')
prefix_panel=('01CD4_' '01CD8_')
pca_score_cutoff=(2 2)




# --------------------------------------------------
# Analysis of CK_2016-06-23_01_CD4_mergingNEW2 using panel1CD4.xlsx for cluster_merging
# and CK_2016-06-23_01_CD8_mergingNEW2 using panel1CD8.xlsx for cluster_merging
# Use Analysis block 2
# --------------------------------------------------

prefix_clust="cl20_"

file_merging=("${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_clust}cluster_merging2.xlsx" "${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_clust}cluster_merging.xlsx")
prefix_merging=('merging2_' 'merging_')

for i in 0 1
do
  ./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust} --prefix_merging ${prefix_merging[$i]} --file_merging ${file_merging[$i]}
done


prefix_clust=('cl8_' 'cl5_')

file_merging=("${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_clust[0]}cluster_merging4.xlsx" "${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_clust[1]}cluster_merging3.xlsx")
prefix_merging=('merging4_' 'merging3_')

for i in 0 1
do
  ./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust[$i]} --prefix_merging ${prefix_merging[$i]} --file_merging ${file_merging[$i]}
done


prefix_clust=('cl5_' 'cl8_')

file_merging=("${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_clust[0]}cluster_merging5.xlsx" "${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_clust[1]}cluster_merging4.xlsx")
prefix_merging=('merging5_' 'merging4_')

for i in 0 1
do
  ./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust[$i]} --prefix_merging ${prefix_merging[$i]} --file_merging ${file_merging[$i]}
done

prefix_clust=('cl5_' 'cl8_')

file_merging=("${prefix_data[0]}${prefix_panel[0]}${prefix_pca}${prefix_clust[0]}cluster_merging5.xlsx" "${prefix_data[1]}${prefix_panel[1]}${prefix_pca}${prefix_clust[1]}cluster_merging5.xlsx")
prefix_merging=('merging5_' 'merging5_')

for i in  1 # run only for CD8
do
  ./Analysis_block_2_cluster_merging.sh --RCODE ${RCODE} --RWD_MAIN ${RWD_MAIN} --data_dir ${data_dir[$i]} --cluster_merging ${cluster_merging} --heatmaps ${heatmaps} --plottsne ${plottsne} --frequencies ${frequencies} --expression ${expression} --METADATA ${METADATA} --PANELS ${PANELS} --file_metadata ${file_metadata} --file_panel ${file_panel[$i]} --prefix_data ${prefix_data[$i]} --prefix_panel ${prefix_panel[$i]} --prefix_pca ${prefix_pca} --prefix_clust ${prefix_clust[$i]} --prefix_merging ${prefix_merging[$i]} --file_merging ${file_merging[$i]}
done
