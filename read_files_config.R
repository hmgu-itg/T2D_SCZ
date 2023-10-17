traits <- c("T2D", "SCZ")
tissues <- c("islets")
qtl.type <- c("eQTL")

# for (trait in traits) {
#   assign(pa35e0(trait, ".filepath"), paste0("/project_data/data/", trait, "_multiancestry.txt"))
#   assign(paste0(trait, ".indep.signals"), paste0("/project_data/data/", trait, "_indep_signals.txt"))
# }
# assign(paste0(traits[1], ".ncases"), 180834)
# assign(paste0(traits[2], ".ncases"), 71554)
# assign(paste0(traits[1], "samplesize"), 1339889)
# assign(paste0(traits[2], "samplesize"), 169417)

pp4.thres <- 0.8

data.path <- paste0(project_folder, "data/")
tmp.path <- paste0(project_folder, "tempdata/")
deg.list <- c(paste0(data.path, "DEG_T2D.csv"), paste0(data.path, "DEG_SCZ.csv"))
genes.data <- paste0(data.path, "UCSC_GRCh37_Genes_UniqueList.txt")
hc.list <- c(paste0(data.path, "T2D_hc.txt"), paste0(data.path, "SCZ_hc.txt"))

output.path <- paste0(project_folder, "GWASColoc/")
coloc.plots.path <- paste0(output.path, "plots/")
qtl.coloc.plots.path <- paste0(output.path, "plots/qtl/")
ld.path <- paste0(output.path, "LDvariants/")
deg.path <- paste0(output.path, "lookups/DEG.csv")
omim.path <- paste0(output.path, "lookups/OMIM.csv")
ko.path <- paste0(output.path, "lookups/KOmice.csv")
missense.path <- paste0(output.path, "lookups/missense.csv")
mgi.result <- paste0(output.path, "lookups/MGIBatchReport_20221012_082318.txt")
hc.path <- paste0(output.path, "lookups/HC.csv")
all.genes <- paste0(output.path, "all_genes.csv")

# QTL data paths
islets.eqtl.path <- paste0(data.path, "bulk/eqtl/eQTL_PancreaticIslets.csv")
liver.eqtl.path <- paste0(data.path, "bulk/eqtl/eQTL_Liver.csv")
brain.eqtl.path <- paste0(data.path, "bulk/eqtl/eQTL_brain.csv")
fetal.eqtl.path <- paste0(data.path, "bulk/eqtl/eQTL_fetal_brain.csv")
cortex_eur.eqtl.path <- paste0(data.path, "bulk/eqtl/eQTL_cortex-EUR-80.csv")
cortex_eas.eqtl.path <- paste0(data.path, "bulk/eqtl/eQTL_cortex-EAS-30.csv")
cortex_afr.eqtl.path <- paste0(data.path, "bulk/eqtl/eQTL_cortex-AFR-40.csv")
hippocampus.eqtl.path <- paste0(data.path, "bulk/eqtl/eQTL_hippocampus-EUR-30.csv")
cerebellum.eqtl.path <- paste0(data.path, "bulk/eqtl/eQTL_cerebellum-EUR-60.csv")
basalganglia.eqtl.path <- paste0(data.path, "bulk/eqtl/eQTL_basalganglia-EUR-30.csv")
adipose.eqtl.path <- paste0(data.path, "bulk/eqtl/eQTL_adipose.csv")
adipose_GTEx.eqtl.path <- paste0(data.path, "bulk/eqtl/eQTL_adipose_GTEx.csv")

adipose.sqtl.path <- paste0(data.path, "bulk/sqtl/sQTL_adipose.csv")

brain.pqtl.path <- paste0(data.path, "bulk/pqtl/pQTL_brain.csv")
# liver.pqtl.path <- paste0(data.path, "bulk/pqtl/pQTL_liver.csv")
liver.pqtl.path <- paste0(data.path, "bulk/pqtl/pQTL_liver_raw.csv")

brain.caqtl.path <- paste0(data.path, "bulk/caqtl/caQTL_brain.csv")

brain.mqtl.path <- paste0(data.path, "bulk/mqtl/mQTL_brain.csv")
fetal.mqtl.path <- paste0(data.path, "bulk/mqtl/mQTL_fetal.csv")

sc_eqtl <- c("Astrocytes", "Endothelial.cells", "Excitatory.neurons", "Inhibitory.neurons", "Microglia", "Oligodendrocytes", "OPCs...COPs", "Pericytes")
for (cell.type in sc_eqtl){
  assign(paste0(cell.type, ".sc_eqtl.path"), paste0(data.path, "single_cell/brain/sc_eqtl_", cell.type, ".csv"))
}

# T2D.filepath <- "/project_data/data/T2D_multiancestry.txt"
# T2D.indep.signals <- "/project_data/data/T2D_indep_signals.txt"
# T2D.ncases <- 180834
# T2D.samplesize <- 1339889
# SCZ.filepath <- "/project_data/data/SCZ_primary.txt"
# SCZ.indep.signals <- "/project_data/data/SCZ_indep_signals.txt"
# SCZ.ncases <- 71554
# SCZ.samplesize <- 169417

t1.filepath <- paste0(data.path, traits[1], "_multiancestry.txt")
t1.indep.signals <- paste0(data.path, traits[1], "_indep_signals.txt")
t1.ncases <- 180834
t1.samplesize <- 1339889
t2d.omim <- c("Aicardi-Goutieres syndrome", "Desbuquois dysplasia 1", "Mitchell-Riley syndrome")
t1.terms <- c("insulin", "glycemia", "glucose", "diabetes", "pancreas", "beta-cell", "beta cell", "beta cells",
              "glucosuria", "body weight", "obesity", "BMI", "body mass", "body fat", "hyperglycemia", "pancreatic",
              "triglycerides", "adipose tissue", "triglyceride")
t1.terms <- c(t1.terms, t2d.omim)

t2.filepath <- paste0(data.path, traits[2], "_multiancestry.txt")
t2.indep.signals <- paste0(data.path, traits[2], "_indep_signals.txt")
t2.ncases <- 71554
t2.samplesize <- 169417

motoric <- c("locomotor",  "tremors", "hyperlocomotion", "grip strength", "limb grasping", "reflex", 
             "no spontaneous movement", "abnormal tail movements", "aphagia ", "ataxia", "dystonia", 
             "hemiparesis", "lethargy", "paralysis", "tactile stimuli",  "abnormal gait", "coordination")
brain <- c("cerebellum", "cerebellar", "cranium", "dendrite", "dendritic", "cranial", "neuromuscular", "rhombomere", 
           "nervous system", "cerebral", "action potential", "nerve conduction", "temporal lobe", "neural",
           "neuron", "synaptic", "synapses", "synapse",  "brain", "spike wave", "hippocampus", "neurotransmitter",
           "astrocytosis",  "myelination", "retrotrapezoid nucleus morphology", "craniopharyngeal ducts",
           "anterior commissure", "corpus callosum", "retrosplenial region", "myelin sheath", "paired-pulse facilitation")
other <- c("eye",  "optic cup", "mandible", "nerve",  "eyelid", "iris", "cornea", "mouth",  "lens", "retina", "craniofacial",
           "vitreous body", "optic disk", "malocclusion")
behaviour <- c("abrainttention", "impulsivity", "learning", "failure to thrive")
psych <- c("schizophrenia", "intelligence", "bipolar", "cognitive", "cognition", "epilepsy", "psychosis", 
           "psychological", "anxiety",  "prepulse inhibition", 
           "hyperactivity",  "novel object", "behavior", "behaviour", "seizures",  "seizure", "kindling response", 
           "delusion", "hallucinations", "paranoia", "disorganized speech", 
           "new environment", "habituation", "response to novelty", "depression", "thigmotaxis", 
           "abnormal social investigation", "abnormal thermal nociception", "abnormal maternal nurturing", "sensory gating", 
           "sensorimotor gating")
scz.omim <- c("Keutel syndrome", "Mitochondrial complex I deficiency, nuclear type 21", "PEHO syndrome", "Peroxisome biogenesis disorder", 
              "Glycosylphosphatidylinositol biosynthesis defect 11", "Acetyl-CoA carboxylase deficiency", "Charcot-Marie-Tooth disease", 
              "Hypomagnesemia", "Spinocerebellar ataxia", "Developmental and epileptic encephalopathy", "Acromesomelic dysplasia 1, Maroteaux type", 
              "Intellectual developmental disorder",  "Lissencephaly 10", "O'Donnell-Luria-Rodan syndrome", "Desbuquois dysplasia 1", 
              "Vertebral, cardiac, tracheoesophageal, renal, and limb defects",  "Aicardi-Goutieres syndrome")
t2.terms <- c(brain, psych, behaviour, other, scz.omim)  # motoric, "concentration", "memory" 


# Supplemental table of related traits
t2d.dt <- list(KO_mice=t1.terms, OMIM=t2d.omim)
t2d.dt <- setDT(lapply(t2d.dt, `length<-`, max(lengths(t2d.dt))))[]
data.table::fwrite(t2d.dt, "t2d_ko_omim_terms.csv")
  
scz.dt <- list(KO_mice=c(brain, other, behaviour, psych), KO_mice_motoric=motoric, OMIM=scz.omim)
scz.dt <- setDT(lapply(scz.dt, `length<-`, max(lengths(scz.dt))))[]
data.table::fwrite(scz.dt, "scz_ko_omim_terms.csv")
