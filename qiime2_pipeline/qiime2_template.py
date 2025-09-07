import os
import subprocess
from pathlib import Path
import shutil

# Helper to run shell commands
def run_cmd(cmd):
    cmd = [str(c) for c in cmd]
    subprocess.run(cmd, check=True)

# Define project root 
project_root = Path("/home/patwuch/projects/microbiome")
batch_serial = "20250905"  

# Define commonly used subdirectories
data_base = project_root / "data" 
processed_base = project_root / "data" / "processed" 
reference_base = project_root / "reference"
raw_base = project_root / "data" / "raw"

# Set specific directories for this batch of analysis
processed_dir = processed_base / batch_serial
raw_dir = raw_base / batch_serial

# Define separate directories for qza and qzv in the processed folder
qza_dir = os.path.join(processed_dir, "qza")
qzv_dir = os.path.join(processed_dir, "qzv")
# Specify metadata files
metadata_file = raw_dir / "metadata.tsv"
# And classifier paths
classifier_path = reference_base / "silva-138-99-nb-classifier.qza"

def create_manifest_from_subdirectories(raw_dir_path: Path, output_manifest_path: Path):
    """
    Scans subdirectories for paired-end fastq files and generates a manifest file.

    Args:
        raw_dir_path: Path to the directory containing sample subdirectories.
        output_manifest_path: Path where the manifest file will be written.
    """
    with open(output_manifest_path, "w") as f:
        # if possible use only lower case for naming the id column, qiime2 is picky
        f.write("sampleid\tforward-absolute-filepath\treverse-absolute-filepath\n")
        
        # Iterate over all subdirectories in the raw data directory
        for sample_folder in os.listdir(raw_dir_path):
            sample_folder_path = raw_dir_path / sample_folder
            
            # Check if the path is a directory
            if os.path.isdir(sample_folder_path):
                forward_file = None
                reverse_file = None
                
                # Iterate over files within the sample subdirectory
                for file in os.listdir(sample_folder_path):
                    if file.endswith('_1.fq.gz'):
                        forward_file = file
                    elif file.endswith('_2.fq.gz'):
                        reverse_file = file

                # Ensure both files were found and write to manifest
                if forward_file and reverse_file:
                    sample_id = sample_folder  # Use the folder name as the sample ID
                    forward_path = sample_folder_path / forward_file
                    reverse_path = sample_folder_path / reverse_file
                    
                    f.write(f"{sample_id}\t{forward_path}\t{reverse_path}\n")

# Define and create manifest file
manifest_file = raw_dir / "manifest.tsv"
# create_manifest_from_subdirectories(raw_dir, manifest_file)

# Define other output directories
core_metrics_dir = os.path.join(qza_dir, "core-metrics-results")
# Export tree and taxa for phyloseq construction if needed
exported_tree_dir = os.path.join(qza_dir, "exported_tree")
exported_taxonomy_dir = os.path.join(qza_dir, "exported_taxonomy")

# Ensure all directories exist
os.makedirs(qza_dir, exist_ok=True)
os.makedirs(qzv_dir, exist_ok=True)
os.makedirs(exported_tree_dir, exist_ok=True)
os.makedirs(exported_taxonomy_dir, exist_ok=True)

# Helpers to put outputs in correct dirs, either qiime artifacts or visualizations
def out_qza(filename):
    return os.path.join(qza_dir, filename)
def out_qzv(filename):
    return os.path.join(qzv_dir, filename)

# ----------------------------
# QIIME2 pipeline starts here
# ----------------------------

# Import
run_cmd([
    "qiime", "tools", "import",
    "--type", "SampleData[PairedEndSequencesWithQuality]",
    "--input-format", "PairedEndFastqManifestPhred33V2",
    "--input-path", manifest_file,
    "--output-path", out_qza("demux.qza")
])

# Summarize demux
run_cmd([
    "qiime", "demux", "summarize",
    "--i-data", out_qza("demux.qza"),
    "--o-visualization", out_qzv("demux.qzv")
])

# DADA2 denoise, change the trim length depending on QC plot from demux.qzv
run_cmd([
    "qiime", "dada2", "denoise-paired",
    "--i-demultiplexed-seqs", out_qza("demux.qza"),
    "--p-trim-left-f", "10",
    "--p-trim-left-r", "10",
    "--p-trunc-len-f", "240",
    "--p-trunc-len-r", "240",
    "--o-table", out_qza("table-dada2.qza"),
    "--o-representative-sequences", out_qza("rep-seqs-dada2.qza"),
    "--o-denoising-stats", out_qza("stats-dada2.qza"),
    "--p-n-threads", "0"
])

# Feature table and metadata summaries
run_cmd([
    "qiime", "feature-table", "summarize",
    "--i-table", out_qza("table-dada2.qza"),
    "--o-visualization", out_qzv("table-dada2.qzv"),
    "--m-sample-metadata-file", metadata_file
])

# At this point, check the table-dada2.qzv to make sure trimming length was appropriate
run_cmd([
    "qiime", "feature-table", "tabulate-seqs",
    "--i-data", out_qza("rep-seqs-dada2.qza"),
    "--o-visualization", out_qzv("rep-seqs-dada2.qzv")
])

run_cmd([
    "qiime", "metadata", "tabulate",
    "--m-input-file", out_qza("stats-dada2.qza"),
    "--o-visualization", out_qzv("stats-dada2.qzv")
])

# Taxonomy classification
run_cmd([
    "qiime", "feature-classifier", "classify-sklearn",
    "--i-classifier", classifier_path,
    "--i-reads", out_qza("rep-seqs-dada2.qza"),
    "--o-classification", out_qza("taxonomy.qza")
])

run_cmd([
    "qiime", "metadata", "tabulate",
    "--m-input-file", out_qza("taxonomy.qza"),
    "--o-visualization", out_qzv("taxonomy.qzv")
])

# Taxa barplots
run_cmd([
    "qiime", "taxa", "barplot",
    "--i-table", out_qza("table-dada2.qza"),
    "--i-taxonomy", out_qza("taxonomy.qza"),
    "--m-metadata-file", metadata_file,
    "--o-visualization", out_qzv("taxa-bar-plots.qzv")
])

run_cmd([
    "qiime", "krona", "collapse-and-plot",
    "--i-table", out_qza("table-dada2.qza"),
    "--i-taxonomy", out_qza("taxonomy.qza"),
    "--o-krona-plot", out_qzv("krona.qzv")
])

# Phylogeny
run_cmd([
    "qiime", "phylogeny", "align-to-tree-mafft-fasttree",
    "--i-sequences", out_qza("rep-seqs-dada2.qza"),
    "--o-alignment", out_qza("aligned-rep-seqs.qza"),
    "--o-masked-alignment", out_qza("masked-aligned-rep-seqs.qza"),
    "--o-tree", out_qza("unrooted-tree.qza"),
    "--o-rooted-tree", out_qza("rooted-tree.qza")
])


# Alpha rarefaction, choose depth based on QC plot
# Because of the low amount of sample in this batch, set to highest value that will include all samples
run_cmd([
    "qiime", "diversity", "alpha-rarefaction",
    "--i-table", out_qza("table-dada2.qza"),
    "--i-phylogeny", out_qza("rooted-tree.qza"),
    "--p-max-depth", "30000",
    "--m-metadata-file", metadata_file,
    "--o-visualization", out_qzv("alpha-rarefaction.qzv")
])

# Delete core-metrics folder if it already exists - qiime would not overwrite it
if os.path.exists(core_metrics_dir):
    shutil.rmtree(core_metrics_dir)

# Cun core metrics, using the sample depth based on previous step
# Change depending on if comparing subgroups requires higher depth value
run_cmd([
    "qiime", "diversity", "core-metrics-phylogenetic",
    "--i-phylogeny", out_qza("rooted-tree.qza"),
    "--i-table", out_qza("table-dada2.qza"),
    "--p-sampling-depth", "30000",
    "--m-metadata-file", metadata_file,
    "--output-dir", core_metrics_dir
])

# Alpha and Beta diversity group significance
# For first version report just use 'Group'
alpha_files = [
    (os.path.join(core_metrics_dir, "faith_pd_vector.qza"), out_qzv("faith-pd-group-significance.qzv")),
    (os.path.join(core_metrics_dir, "evenness_vector.qza"), out_qzv("evenness-group-significance.qzv"))
]

for alpha_input, alpha_output in alpha_files:
    run_cmd([
        "qiime", "diversity", "alpha-group-significance",
        "--i-alpha-diversity", alpha_input,
        "--m-metadata-file", metadata_file,
        "--o-visualization", alpha_output
    ])

beta_tests = [
    (os.path.join(core_metrics_dir, "unweighted_unifrac_distance_matrix.qza"), "Group", out_qzv("unweighted-unifrac-Group-significance.qzv"))
]

for dist, col, out in beta_tests:
    run_cmd([
        "qiime", "diversity", "beta-group-significance",
        "--i-distance-matrix", dist,
        "--m-metadata-file", metadata_file,
        "--m-metadata-column", col,
        "--o-visualization", out,
        "--p-pairwise"
    ])

# Run ANCOMBC2 on all MainType with Modifier as random effect
run_cmd([
    "qiime", "composition", "ancombc2",
    "--i-table", out_qza("table-dada2.qza"),   
    "--m-metadata-file", metadata_file,
    "--p-fixed-effects-formula", "MainType",
    "--o-ancombc2-output", out_qza("ancombc2-MainType-Modifier-results.qza")
])

# Visualization with taxonomy
run_cmd([
    "qiime", "composition", "ancombc2-visualizer",
    "--i-data", out_qza("ancombc2-MainType-Modifier-results.qza"),
    "--i-taxonomy", out_qza("taxonomy.qza"),
    "--o-visualization", out_qzv("ancombc2-MainType-Modifier.qzv")
])

# Beta-group significance ADONIS(permanova) or BETADISPER(permdisp)
columns = ["Group","MainType","Modifier"]
distances = [f for f in os.listdir(core_metrics_dir) if f.endswith("_distance_matrix.qza")]

for dist in distances:
    dist_path = os.path.join(core_metrics_dir, dist)
    base = dist.replace("_distance_matrix.qza", "")
    for col in columns:
        # PERMANOVA
        run_cmd([
            "qiime", "diversity", "beta-group-significance",
            "--i-distance-matrix", dist_path,
            "--m-metadata-file", metadata_file,
            "--m-metadata-column", col,
            "--p-method", "permanova",
            "--o-visualization", out_qzv(f"permanova_{base}_{col}.qzv")
        ])
        # PERMDISP
        run_cmd([
            "qiime", "diversity", "beta-group-significance",
            "--i-distance-matrix", dist_path,
            "--m-metadata-file", metadata_file,
            "--m-metadata-column", col,
            "--p-method", "permdisp",
            "--o-visualization", out_qzv(f"permdisp_{base}_{col}.qzv")
        ])

# Export tree and taxonomy (directories -> keep under qza_dir)
run_cmd([
    "qiime", "tools", "export",
    "--input-path", out_qza("rooted-tree.qza"),
    "--output-path", os.path.join(qza_dir, "exported_tree")
])

run_cmd([
    "qiime", "tools", "export",
    "--input-path", out_qza("taxonomy.qza"),
    "--output-path", os.path.join(qza_dir, "exported_taxonomy")
])

print("All QIIME2 commands finished successfully!")
