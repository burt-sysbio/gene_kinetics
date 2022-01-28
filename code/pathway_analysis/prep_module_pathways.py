# load tfs, cytos, cytoR, th2 th1 genes
import pandas as pd

sign_dir = "../../genesets_literature/"
tf_genes = pd.read_csv(sign_dir + "references/signaling/stubbington_tf.csv", names= ["gene"])
cyto_genes = pd.read_csv(sign_dir + "references/signaling/stubbington_cyto.csv", names= ["gene"])
cytor_genes = pd.read_csv(sign_dir + "references/signaling/stubbington_cytor.csv", names= ["gene"])

th1_genes = pd.read_csv(sign_dir + "references/effector_states/stubbington_th1.txt", names = ["gene"])
th2_genes = pd.read_csv(sign_dir + "references/effector_states/stubbington_th2.txt", names = ["gene"])
th17_genes = pd.read_csv(sign_dir + "references/effector_states/stubbington_th17.txt", names = ["gene"])

prolif_genes = pd.read_csv(sign_dir + "tcell_module/go_tcell_proliferation.txt", sep = "\t", usecols =[1], header=0, names = ["gene"])
diff_genes = pd.read_csv(sign_dir + "tcell_module/go_thelper_differentiation.txt", sep = "\t", usecols =[1], header=0, names = ["gene"])
il2_genes = pd.read_csv(sign_dir + "tcell_module/hallmark_il2_stat5_signaling.txt", skiprows = 2, header = 0, names = ["gene"])
il2_genes["gene"] = il2_genes["gene"].str.title()

tcell_activation = pd.read_csv(sign_dir + "tcell_module/go_tcell_activation.txt", sep = "\t", usecols =[1], header=0, names = ["gene"])


my_modules = [tf_genes, cyto_genes, cytor_genes, th1_genes, th17_genes, th2_genes,
              prolif_genes, diff_genes, il2_genes, tcell_activation]

module_names = ["TF", "CytoK", "CytoR", "Th1",
                "Th2", "Th17", "Proliferation", "Differentiation", "IL2_signal", "Activation"]

for df, name in zip(my_modules, module_names):
    df.drop_duplicates(inplace=True)
    df["module"] = name
# from each data frame, remove duplicates and collapse to list
my_modules = [df.drop_duplicates() for df in my_modules]
df_modules = pd.concat(my_modules)

df_modules.to_csv("../../genesets_literature/gene_module_summary.csv", index=False)