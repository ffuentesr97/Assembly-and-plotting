import sys
import os
import numpy as np
from Bio import SeqIO
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from pycirclize import Circos
from pycirclize.parser import Genbank, Gff

# Ensure the script is called with the correct number of arguments
if len(sys.argv) != 4:
    print("Usage: python plot_gbk.py <gbk_file> <genome_name> <gff_file>")
    sys.exit(1)

# Load GenBank file and genome name from command line arguments
gbk_file = sys.argv[1]
genome_name = sys.argv[2]
gff_file = sys.argv[3]


# Load the GenBank file
gbk = Genbank(gbk_file)
gff = Gff(gff_file)

# List of gene names to highlight
highlight_genes_nb2 = {"nodZ", "noeI", "nolK"}  # Replace with your target gene names
highlight_genes_nb8 = {"nodA", "nodB", "nodC", "nodI", "nodJ", "nolO"} 
highlight_genes_reg = {"nodD1", "nodD2", "syrM", "ttsI"}
highlight_genes_fix = {"fixA", "fixB", "fixC", "fixX"}
 
# Calculate genome size in base pairs
genome_size_bp = sum(len(record) for record in gbk.records)

circos = Circos(sectors={gbk.name: gbk.range_size})
circos.text(f"Sinorhizobium fredii HH103\n({genome_name})\n{genome_size_bp} bp", size=12, r=21)
sector = circos.get_sector(gbk.name)

# Plot outer track with xticks
major_ticks_interval = 50000
minor_ticks_interval = 25000
outer_track = sector.add_track((98, 100))
outer_track.axis(fc="dimgrey")
outer_track.xticks_by_interval(
    major_ticks_interval, label_formatter=lambda v: f"{v/ 1000:.01f} kb"
)
outer_track.xticks_by_interval(minor_ticks_interval, tick_length=1, show_label=False)

# Plot Forward CDS, Reverse CDS, rRNA, tRNA
f_cds_track = sector.add_track((90, 97), r_pad_ratio=0.1)
f_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=1), fc="tomato")
               
r_cds_track = sector.add_track((83, 90), r_pad_ratio=0.1)
r_cds_track.genomic_features(gbk.extract_features("CDS", target_strand=-1), fc="deepskyblue")

#rrna_track = sector.add_track((69, 76), r_pad_ratio=0.1)
#rrna_track.genomic_features(gbk.extract_features("rRNA"), fc="green", lw=0.1)

trna_track = sector.add_track((76, 83), r_pad_ratio=0.1)
trna_track.genomic_features(gbk.extract_features("tRNA"), color="magenta", lw=0.8)

#Highlight_genes_positive
f_cds_track = sector.add_track((70, 76), r_pad_ratio=0.1)
for f_cds_feature in gff.extract_features("CDS", target_strand=1):
    gene_name = f_cds_feature.qualifiers.get("gene_name", [None])[0]
    if gene_name in highlight_genes_nb2:
        f_cds_track.genomic_features(f_cds_feature, fc="limegreen", ec="limegreen", lw=2)  
    if gene_name in highlight_genes_nb8:
        f_cds_track.genomic_features(f_cds_feature, fc="chocolate", ec="chocolate", lw=2)
    if gene_name in highlight_genes_reg:
        f_cds_track.genomic_features(f_cds_feature, fc="darkgreen", ec="darkgreen", lw=2)
    if gene_name in highlight_genes_fix:
        f_cds_track.genomic_features(f_cds_feature, fc="crimson", ec="crimson", lw=2)
       
#Highlight_genes_negative
r_cds_track = sector.add_track((70, 76), r_pad_ratio=0.1)
for r_cds_feature in gff.extract_features("CDS", target_strand=-1):
    gene_name = r_cds_feature.qualifiers.get("gene_name", [None])[0]
    if gene_name in highlight_genes_nb2:
        r_cds_track.genomic_features(r_cds_feature, fc="limegreen", ec="limegreen", lw=2)
    if gene_name in highlight_genes_nb8:
        r_cds_track.genomic_features(r_cds_feature, fc="chocolate", ec="chocolate", lw=2)        
    if gene_name in highlight_genes_reg:
        r_cds_track.genomic_features(r_cds_feature, fc="darkgreen", ec="darkgreen", lw=2)        
    if gene_name in highlight_genes_fix:
        r_cds_track.genomic_features(r_cds_feature, fc="crimson", ec="crimson", lw=2)        
           

# Plot GC content
gc_content_track = sector.add_track((55, 70))

pos_list, gc_contents = gbk.calc_gc_content()
gc_contents = gc_contents - gbk.calc_genome_gc_content()
positive_gc_contents = np.where(gc_contents > 0, gc_contents, 0)
negative_gc_contents = np.where(gc_contents < 0, gc_contents, 0)
abs_max_gc_content = np.max(np.abs(gc_contents))
vmin, vmax = -abs_max_gc_content, abs_max_gc_content
gc_content_track.fill_between(
    pos_list, positive_gc_contents, 0, vmin=vmin, vmax=vmax, color="black"
)
gc_content_track.fill_between(
    pos_list, negative_gc_contents, 0, vmin=vmin, vmax=vmax, color="grey"
)

# Plot GC skew
gc_skew_track = sector.add_track((40, 55))

pos_list, gc_skews = gbk.calc_gc_skew()
positive_gc_skews = np.where(gc_skews > 0, gc_skews, 0)
negative_gc_skews = np.where(gc_skews < 0, gc_skews, 0)
abs_max_gc_skew = np.max(np.abs(gc_skews))
vmin, vmax = -abs_max_gc_skew, abs_max_gc_skew
gc_skew_track.fill_between(
    pos_list, positive_gc_skews, 0, vmin=vmin, vmax=vmax, color="olive"
)
gc_skew_track.fill_between(
    pos_list, negative_gc_skews, 0, vmin=vmin, vmax=vmax, color="purple"
)

fig = circos.plotfig()

# Add legend
handles = [
    Patch(color="tomato", label="Forward CDS"),
    Patch(color="deepskyblue", label="Reverse CDS"),
    #Patch(color="green", label="rRNA"),
    Patch(color="magenta", label="tRNA"),
    Patch(color="darkgreen", label="Main symbiotic regulators"),
    Patch(color="limegreen", label="NB2 controlled genes"),
    Patch(color="chocolate", label="NB8 controlled genes"),
    Patch(color="crimson", label="fix genes"),
    Line2D([], [], color="black", label="Positive GC Content", marker="^", ms=6, ls="None"),
    Line2D([], [], color="grey", label="Negative GC Content", marker="v", ms=6, ls="None"),
    Line2D([], [], color="olive", label="Positive GC Skew", marker="^", ms=6, ls="None"),
    Line2D([], [], color="purple", label="Negative GC Skew", marker="v", ms=6, ls="None"),
]
_ = circos.ax.legend(handles=handles, bbox_to_anchor=(0.52, 0.44), loc="center", fontsize=8)


# To save the figure
fig.savefig(f"result_plot_{genome_name}.png", dpi=1200)
