export PATH=STAR-2.7.3a/bin/Linux_x86_64:$PATH
export PATH=subread-2.0.0-source/bin:$PATH

cd STAR/GRCm38_100

STAR -- runMode genomeGenerate \
     -- genomeDir . \
     -- genomeFastaFiles Mus_musculus.GRCm38.dna.primary_assembly.fa \
     -- sjdbGTFfile Mus_musculus.GRCm38.100.gtf \
     -- runThreadN 14 

export PATH=/data/tools/FastQC:$PATH

for f in *fastq.gz
do 
  fastqc -t 38 $f -o qc
done 

multiqc qc

TRIM='Trimmomatic-0.38/trimmomatic-0.38.jar'

for f in *.fastq.gz
do 
  o=${f/'.fastq.gz'/'_trimmed.fastq'}
  java -jar $TRIM SE -threads 38 \
    $f $o \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:30 \
    MINLEN:50
done

idx='STAR/GRCm38_100'

for f in *._trimmed.fastq
do
  root=${f/'.fastq'/''}
  STAR -- genomeDir $idx \
       -- readFilesIn $f \
       -- outFileNamePrefix 'GRCm38_100_'$root \
       -- outFilterMultimapNmax 1 \
       -- outReadsUnmapped Fastx \
       -- outSAMtype BAM SortedByCoordinate \
       -- twopassMode Basic \
       -- runThreadN 14 
done

featureCounts -O \
  -a $idx/Mus_musculus.GRCm38.100.gtf \
  -o featureCounts_GRCm38_100_O.txt \
  GRCm38_100*bam 
 
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-100/tsv/mus_musculus/Mus_musculus.GRCm38.100.uniprot.tsv.gz .
rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-100/tsv/mus_musculus/Mus_musculus.GRCm38.100.entrez.tsv.gz .
wget ftp://ftp.pantherdb.org/sequence_classifications/current_release/PANTHER_Sequence_Classification_files/PTHR15.0_mouse
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/MOUSE_10090_idmapping_selected.tab.gz
wget http://current.geneontology.org/annotations/mgi.gaf.gz
wget https://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz
wget https://ftp.ncbi.nih.gov/gene/DATA/gene2accession.gz


# R

library(org.Mm.eg.db)
library(GO.db)
library(KEGGREST)
library(edgeR)
library(DESeq2)
library(openxlsx)
library(tidyverse)

org.db = org.Mm.eg.db

uniprot.db_sel = read_tsv('MOUSE_10090_idmapping_selected.tab', col_names=c('UniProtKB_AC', 'UniProtKB_ID', 'GeneID_EntrezGene', 'RefSeq', 'GI', 'PDB', 'GO', 'UniRef100', 'UniRef90', 'UniRef50', 'UniParc', 'PIR', 'NCBI_taxon', 'MIM', 'UniGene', 'PubMed', 'EMBL', 'EMBL_CDS', 'Ensembl', 'Ensembl_TRS', 'Ensembl_PRO', 'Additional_PubMed')) 

PTHR = read_tsv(
    'PTHR15.0_mouse', 
    col_names=c('GeneIdentifier', 'ProteinID', 'SFID', 'FamilyName', 'SubfamilyName', 'MolecularFunction', 'BiologicalProcess', 'CellularComponents', 'ProteinClass', 'Pathway'), 
    quote='') %>% 
  separate(GeneIdentifier, c('GeneIdentifier', 'UniProt'), 'UniProtKB=')
    
ens2uniprot = read_tsv('Mus_musculus.GRCm38.100.uniprot.tsv', col_names=T) %>%
  rename(UniProt=xref, ens=gene_stable_id) %>%
  select(ens, UniProt) %>% 
  unique()

go = read_tsv('mgi.gaf', skip=43, col_names=F, 
  col_types = cols(.default = "c")) %>%
  select(2,5,7,9,10,11) %>%
  mutate(X9=case_when(X9=='C' ~ 'CC', X9=='P' ~ 'BP', X9=='F' ~ 'MF')) %>%
  rename(UniProt=1, GOid=2, evidence=3, ontology=4, GeneName=5, symbol=6)  

read.counts = read_delim('featureCounts_GRCm38_100_O.txt', delim='\t', col_names=T, skip=1) %>% 
  rename_at(vars(contains('GRCm38_100_')), list( ~ gsub('GRCm38_100_', '', .))) %>%  
  rename_at(vars(contains('_trimmedAligned.sortedByCoord.out.bam')), list( ~ gsub('_trimmedAligned.sortedByCoord.out.bam', '', .))) %>%  
  rename_at(vars(contains('-')), list( ~ paste0('s', gsub('-', '_', .)))) 
  
readcounts = read.counts %>% 
  select(grep('arc', tolower(colnames(.)))) %>% 
  as.data.frame()
rownames(readcounts) = read.counts$Geneid

sample_info = data.frame(smpl=colnames(readcounts))
rownames(sample_info) = sample_info$smpl
sample_info$grp = 'E2'
sample_info$grp[grep('OIL', sample_info$smpl)] = 'OIL'
sample_info$grp = factor(sample_info$grp)

dds = DESeqDataSetFromMatrix(
    countData = readcounts, 
    colData = sample_info, 
    design = ~ grp 
)

dds = estimateSizeFactors(dds)

counts_raw = readcounts
colnames(counts_raw) = paste0('raw_', colnames(counts_raw))
counts_raw$ens = rownames(counts_raw)
counts_raw = as_tibble(counts_raw)

counts_normalized = as.data.frame(counts(dds, normalized=T))
colnames(counts_normalized) = paste0('norm_', colnames(counts_normalized))
counts_normalized$ens = rownames(counts_normalized)
counts_normalized = as_tibble(counts_normalized)

m = as.matrix(readcounts)
counts_cpm = as.data.frame(cpm(m))
colnames(counts_cpm) = paste0('cpm_', colnames(counts_cpm))
counts_cpm$ens = rownames(counts_cpm)
counts_cpm = as_tibble(counts_cpm)

ids = rownames(readcounts)
n = 1
tmp = tibble(.rows=0, ens='', UniProt='')
for(ens in ids){
  tmp = rbind(tmp, 
    tibble(ens, UniProt= uniprot.db_sel %>% 
    filter(str_detect(Ensembl, ens)) %>%
    pull(UniProtKB_AC)  
    )
  )
  n=n+1
  print(n)
}

ens2uniprot = rbind(ens2uniprot, tmp) %>% 
  unique()

ens_uniprot_pthr = inner_join(ens2uniprot, PTHR)

tib = left_join(tibble(ens=rownames(readcounts)), ens_uniprot_pthr) %>% 
  select(ens, UniProt)
  
tib = inner_join(tib, counts_raw)
tib = inner_join(tib, counts_normalized)
tib = inner_join(tib, counts_cpm)


dds.dif = DESeq(dds)

res = results(dds.dif, contrast=c('grp', 'E2', 'OIL'))
fix = res[tib$ens,]
tib$log2FC = fix$log2FoldChange
tib$pvalue = fix$pvalue
tib$padj = fix$padj

annot = left_join(
  tib,
  PTHR %>% 
    select(UniProt, FamilyName, SubfamilyName, ProteinClass, BiologicalProcess, CellularComponents, MolecularFunction)
)

evs = sort(unique(go$evidence))
onts = sort(unique(go$ontology))

for(ont in onts){
  for(ev in evs){
    tmp = go %>% 
      filter(ontology==ont, evidence==ev) %>% 
      select(UniProt, GeneName)
    if(dim(tmp)[1]>0){ 
      lst = split(tmp$GeneName, tmp$UniProt) %>% 
        lapply(unique) %>% 
        lapply(sort) %>%  
        lapply(paste, collapse='\n')   
      annot = left_join(
        annot, 
        tibble(UniProt=names(lst), tmp=as.character(lst)) %>%
          rename_at(vars(tmp), ~ paste(ont, ev, sep='_'))
        )
    }
  }
}

lst = keggList('pathway', 'mmu') 
paths = tibble(
  PATH=gsub('path:mmu', '', names(lst)), 
  pathway=gsub(' - Mus musculus \\(mouse\\)', '', as.character(lst))
) 

i = 1
query = keggGet(paste0('mmu', paths$PATH[i]))
kegg = as_tibble(matrix(query[[1]]$GENE, nc=2, byrow=T)) %>% 
  rename(GeneID=1, descr=2) %>%
  mutate(PATH=paths$PATH[i])
  
for(i in 2:dim(paths)[1]){
  query = keggGet(paste0('mmu', paths$PATH[i]))
  if(!is.null(query[[1]]$GENE)){
    kegg = rbind(kegg, 
      as_tibble(matrix(query[[1]]$GENE, nc=2, byrow=T)) %>% 
      rename(GeneID=1, descr=2) %>%
      mutate(PATH=paths$PATH[i])
    )
  }
}

ens2entrez = read_tsv('Mus_musculus.GRCm38.100.entrez.tsv', 
  col_types = cols(.default = "c")) 
  
tmp = inner_join(
  inner_join(kegg,  
    inner_join(
      tibble(gene_stable_id=rownames(readcounts)),
      ens2entrez
    ) %>% 
    select(gene_stable_id, xref) %>% 
    unique() %>% 
    rename(ens=1, GeneID=2)
  ) %>% 
  select(ens, PATH),
  paths
)

lst = split(tmp$pathway, tmp$ens) %>% 
  lapply(unique) %>% 
  lapply(sort) %>%  
  lapply(paste, collapse='\n') 
kegg_res = tibble(ens=names(lst), KEGG=as.character(lst))  

annot = left_join(annot, kegg_res) %>% 
  rename(Ensembl=1)
  
cs = createStyle(wrapText=T)
wb = createWorkbook() 
addWorksheet(wb, 'ARC_with_annotation')
writeData(wb, 1, annot)  
addStyle(wb, 1, style=cs, rows=-1, cols=-1)
saveWorkbook(wb, 'ARC_with_annotation_GRCm38_100.xlsx', overwrite = TRUE)

####################################################################################
# R

library(org.Mm.eg.db)
library(GO.db)
library(KEGGREST)
library(edgeR)
library(DESeq2)
library(openxlsx)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)

gene_mean = counts_cpm %>% select(-7) %>% rowMeans()
sel_ens = counts_cpm$ens[which(gene_mean>10)]
sel_dds = dds[sel_ens]
sel_dds.dif = DESeq(sel_dds)

matcol = rev(colorRampPalette(brewer.pal(11, 'RdBu'))(100))

base = inner_join(
counts_normalized,
res %>% as_tibble() %>% mutate(ens=rownames(res))
)

wd = base %>% 
filter(padj<=0.05) %>% 
arrange(desc(log2FoldChange)) %>% 
column_to_rownames('ens') %>%
select(1:6) 

colnames(wd) = gsub('norm_', '', colnames(wd))

phann = sample_info %>% 
  tibble() %>% 
  rename(Group=2) %>% 
  mutate(smpl=gsub('_ARC', '', smpl)) %>% 
  mutate(smpl=substr(smpl,1,2)) %>% 
  mutate(Group=relevel(factor(Group), 'OIL')) %>% 
  column_to_rownames('smpl')
  
cord = rownames(phann)[c(which(phann$Group=='OIL'), which(phann$Group!='OIL'))]

ann_colors = list(Group=c(OIL="yellow",E2="firebrick"))

wd = wd[,cord]


matcol=colorRampPalette(brewer.pal(11, 'RdBu'))(100)

wd = inner_join(
counts_normalized,
res %>% as_tibble() %>% mutate(ens=rownames(res))
) %>% 
filter(padj<=0.05) %>% 
column_to_rownames('ens') %>%
select(1:6) 

colnames(wd) = gsub('norm_', '', colnames(wd))
colnames(wd) = substr(colnames(wd),1,2)

phann = sample_info %>% 
  tibble() %>% 
  rename(Group=2) %>% 
  mutate(smpl=gsub('norm_', '', smpl)) %>% 
  mutate(smpl=substr(smpl,1,2)) %>% 
  mutate(Group=relevel(factor(Group), 'OIL')) %>% 
  column_to_rownames('smpl')
  
cord = rownames(phann)[c(which(phann$Group=='OIL'), which(phann$Group!='OIL'))]

ann_colors = list(Group=c(OIL="yellow",E2="firebrick"))

wd = wd[,cord]

fontsize = 10

wd %>% 
  pheatmap(
  legend=T, 
  annotation_col=phann, annotation_colors=ann_colors,
  scale='row', 
  border_color=NA, 
  color=matcol, 
  width=7, height=10,
  show_rownames=F, fontfamily='sans', fontsize=fontsize*1.2, cluster_rows=F, cluster_cols=F,
  filename='figs/Fig1d.pdf'
)  


library(ComplexHeatmap)
library(Cairo)

lst = mapIds(org.Mm.eg.db, ens2entrez %>% pull(xref), 'SYMBOL', 'ENTREZID')
entrez_tab = tibble(xref=names(unlist(lst)), symbol=as.character(unlist(lst))) %>% unique()

wd = base %>% filter(padj<=0.05) 

tmp = left_join(
left_join(
wd,
ens2entrez %>% select(gene_stable_id, xref) %>% unique() %>% rename(ens=1)
),
entrez_tab
) %>% 
mutate(symbol=case_when(is.na(symbol) ~ ens, TRUE~symbol)) %>%
mutate(symbol=make.unique(symbol)) %>% 
filter(str_detect(symbol, '\\.', negate=T)) 

top_up = tmp %>% 
  filter(log2FoldChange>0) %>%
  arrange(desc(log2FoldChange)) %>% 
  slice_head(n=25)

top_down = tmp %>% 
  filter(log2FoldChange<0) %>%
  arrange(log2FoldChange) %>% 
  slice_head(n=25)

pm = rbind(top_up, top_down) %>%
  arrange(desc(log2FoldChange))  
  
pd = pm %>%
  column_to_rownames('symbol') %>%
  select(1:6)

colnames(pd) = gsub('norm_', '', colnames(pd))
colnames(pd) = substr(colnames(pd),1,2)

mat_scaled = t(scale(t(pd)))

ha = HeatmapAnnotation(df=phann, 
 col=list(Group=c(OIL='yellow', E2='firebrick')), 
 annotation_name_gp = gpar(fontsize = 12*1.2),
 annotation_legend_param = list(
   title_gp=gpar(fontsize = 10*1.2, fontface="bold"),
    labels_gp = gpar(fontsize = 10*1.2)
 )
)

ats = sort(c(0,round(range(mat_scaled),1)))

library(circlize)

cols = circlize::colorRamp2(ats, rev(brewer.pal(11, 'RdBu')[c(1,6,11)]))

ht1 = Heatmap(mat_scaled, 
  col=cols, 
  top_annotation=ha, 
  row_names_side='left', 
  name='Z-score', 
  column_order=order(phann$Group), 
  cluster_rows=F,
  row_names_gp = gpar(fontsize = 12*1.2),
  column_names_gp = gpar(fontsize = 12*1.2),
  heatmap_legend_param=list(at=ats, 
    legend_height=unit(5, 'cm'), 
    title_gp=gpar(fontsize = 10*1.2, fontface="bold"),
    labels_gp = gpar(fontsize = 10*1.2)
  )
)

ht_list = ht1 + 
rowAnnotation(
 log2FC=anno_barplot(
   pm$log2FoldChange, 
   width=unit(3, "cm"),
   axis_param = list(gp=gpar(fontsize=8*1.2))
 ), annotation_name_gp=gpar(fontsize=12*1.2)
)

CairoPDF('figs/Fig2b.pdf', height=11, width=9) 
draw(ht_list)
graphics.off()


library(EnhancedVolcano)

ppd = left_join(
 left_join(res %>% 
  as_tibble() %>% 
  mutate(ens=rownames(res)),
  ens2entrez %>% select(gene_stable_id, xref) %>% rename(ens=1) %>% unique()
), 
entrez_tab
) %>% 
mutate(symbol=case_when(is.na(symbol) ~ ens, TRUE~symbol)) %>%
mutate(symbol=make.unique(symbol)) %>%
mutate(semmi='')

CairoPDF('figs/Fig2d.pdf', height=10, width=13) 
EnhancedVolcano(ppd, title='', subtitle='', xlim = c(-8, 8),
    lab = ppd$semmi,
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    FCcutoff = 1, col = c('grey30', 'grey30', 'royalblue', 'red2'),
    pointSize = 3.0,
    legendLabels = c('NS', expression(Log[2]~FC), 'p-value', expression(p-value~and~log[2]~FC)),
    .legend = c('NS','Log2 FC','P','P & Log2 FC'),
    caption = bquote(~Log[2]~"fold change cutoff: 1; adjusted p-value cutoff: 0.05")
)
graphics.off()

library(readxl)

CairoPDF('figs/Fig2d_connected.pdf', height=10, width=13) 
EnhancedVolcano(ppd, title='', subtitle='', xlim = c(-8, 8),
    lab = ppd$symbol, selectLab=labs, 
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    FCcutoff = 1, col = c('grey30', 'grey30', 'royalblue', 'red2'),
    pointSize = 3.0,
    legendLabels = c('NS', expression(Log[2]~FC), 'p-value', expression(p-value~and~log[2]~FC)),
    .legend = c('NS','Log2 FC','P','P & Log2 FC'),
    caption = bquote(~Log[2]~"fold change cutoff: 1; adjusted p-value cutoff: 0.05"), 
    drawConnectors=T
)
graphics.off()



library(DOSE)
library(magrittr)
library(clusterProfiler)
library(org.Mm.eg.db)

tres = inner_join(res %>% 
  as_tibble() %>% 
  mutate(ens=rownames(res)),
  ens2entrez %>% select(gene_stable_id, xref) %>% rename(ens=1, entrezid=2) %>% unique()
) %>% 
filter(!is.na(log2FoldChange)) %>% 
arrange(desc(log2FoldChange))


prb = tres %>% group_by(entrezid) %>% summarize(fc=max(log2FoldChange)) %>% arrange(desc(fc))
ged = prb %>% pull(fc)
names(ged) = prb %>% pull(entrezid)      
enrich = gseKEGG(geneList=ged, organism='mmu', minGSSize=1, maxGSSize=5000, pvalueCutoff=1, eps=1e-20, verbose=F)                 
renrich = setReadable(enrich, 'org.Mm.eg.db', 'ENTREZID')

WriteXLS(as.data.frame(renrich), ExcelFileName='figs/ARC_GSEA_KEGG.xls', SheetNames='GSEA')



geneList = tres %>% pull(log2FoldChange)
names(geneList) = tres %>% pull(entrezid)

de = tres %>% filter(padj<0.05) %>% pull(entrezid)

kk = enrichKEGG(gene=de, organism='mmu', pvalueCutoff=0.05)

pkk = 
kk[
kk@result$Description=='Ribosome' |
kk@result$Description=='Protein processing in endoplasmic reticulum' |
kk@result$Description=='Axon guidance' |
kk@result$Description=='Neuroactive ligand-receptor interaction' |
kk@result$Description=='GABAergic synapse' |
kk@result$Description=='Glutamatergic synapse' |
kk@result$Description=='Cholinergic synapse' |
kk@result$Description=='Dopaminergic synapse' |
kk@result$Description=='Estrogen signaling pathway' |
kk@result$Description=='Synaptic vesicle cycle' |
kk@result$Description=='Gap junction' |
kk@result$Description=='cAMP signaling pathway' |
kk@result$Description=='Calcium signaling pathway' |
kk@result$Description=='ErbB signaling pathway' |
kk@result$Description=='cGMP-PKG signaling pathway' |
kk@result$Description=='FoxO signaling pathway' |
kk@result$Description=='PI3K-Akt signaling pathway' |
kk@result$Description=='MAPK signaling pathway' |
kk@result$Description=='Ras signaling pathway' |
kk@result$Description=='AMPK signaling pathway' |
kk@result$Description=='mTOR signaling pathway' |
kk@result$Description=='Rap1 signaling pathway', 
asis=T
]

mpl = 1.2
CairoPDF('figs/ARC_dotplot.pdf', height=10, width=10) 
dotplot(pkk, showCategory=dim(pkk)[1]) +
theme(
  axis.text.x=element_text(size=12*mpl),
  axis.text.y=element_text(size=12*mpl),
  legend.title=element_text(size=11*mpl),
  axis.title=element_text(size=12*mpl),
  legend.text=element_text(size=9*mpl)
)
graphics.off()    
 

edox = setReadable(kk, 'org.Mm.eg.db', 'ENTREZID')

library(WriteXLS)

WriteXLS(as.data.frame(edox), ExcelFileName='figs/ARC_ORA_KEGG.xls', SheetNames='ORA')

kka = enrichKEGG(gene=de, organism='mmu', pvalueCutoff=0.99)
edoxa = setReadable(kka, 'org.Mm.eg.db', 'ENTREZID')

CairoPDF('figs/ARC_net_7x7.pdf', height=7, width=7) 
cnetplot(
edoxa[
edoxa@result$Description=='Neuroactive ligand-receptor interaction' | 
edoxa@result$Description=='GABAergic synapse' | 
edoxa@result$Description=='Glutamatergic synapse' | 
edoxa@result$Description=='Cholinergic synapse' | 
edoxa@result$Description=='Dopaminergic synapse' | 
edoxa@result$Description=='Serotonergic synapse', 
asis=T
], 
foldChange=geneList, showCategory=6
) + theme(legend.position = c(.05, .3), legend.direction = "vertical", legend.box = "vertical") + scale_colour_gradientn(colours = rev(brewer.pal(11, 'RdBu')[-c(6)]), name ="log2FC") + 
guides(fill=guide_legend(order=0), size=guide_legend(order=1))
graphics.off()

library(ToPASeq)
library(graphite)
library(WriteXLS)

cmat = readcounts[rowSums(readcounts)>0,]
group = ifelse(as.character(sample_info$grp)=='E2', 1,0)

pwys = pathways(species="mmusculus", database="kegg")
pwys = graphite::convertIdentifiers(pwys, "ENSEMBL")

spi = SPIA(cmat, group, pwys, type="RNASeq", logFC.th=-1, test.method="DESeq2")

WriteXLS(res(spi)$results, row.names=T, ExcelFileName='figs/ARC_SPIA_KEGG.xls', SheetNames='SPIA')


