#!/usr/bin/env python3
"""
GUI wrapper for uht-tooling using Gradio.
This is a NON-DESTRUCTIVE wrapper that calls existing scripts.
"""

import os
import sys
import logging
import shutil
from pathlib import Path

# Import existing script modules
import scripts.nextera_designer as nextera_designer
import scripts.umi_hunter as umi_hunter
import scripts.design_slim as design_slim
import scripts.mutation_caller as mutation_caller
import scripts.design_gibson as design_gibson
import scripts.profile_inserts as profile_inserts
import scripts.mut_rate as mut_rate

try:
    import gradio as gr
    import pandas as pd
    from Bio import SeqIO
    import matplotlib
    matplotlib.use('Agg')
except ImportError as e:
    print(f"Missing dependency: {e}")
    print("Please install: pip install gradio")
    sys.exit(1)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("gui")

def run_nextera_designer(forward_primer, reverse_primer):
    """Run Nextera designer with uploaded primers."""
    try:
        os.makedirs("data/nextera_designer", exist_ok=True)
        csv_path = "data/nextera_designer/nextera_designer.csv"
        
        with open(csv_path, 'w') as f:
            f.write("binding_region\n")
            f.write(f"{forward_primer}\n")
            f.write(f"{reverse_primer}\n")
        
        logger.info("Running nextera_designer...")
        nextera_designer.main()
        
        results_file = "results/nextera_designer/nextera_xt_primers.csv"
        if os.path.exists(results_file):
            df = pd.read_csv(results_file)
            result_text = f"Generated {len(df)} primers\n\nFirst 5 primers:\n"
            result_text += df.head().to_string(index=False)
            
            log_file = "results/nextera_designer/nextera_designer.log"
            if os.path.exists(log_file):
                with open(log_file, 'r') as f:
                    log_content = f.read()
                    result_text += f"\n\n--- Log ---\n{log_content[-1000:]}"
            
            return result_text
        else:
            return "Error: Results file not found"
    except Exception as e:
        logger.exception(f"Error: {e}")
        return f"Error: {str(e)}"

def run_design_slim(template_gene_content, context_content, mutations_text):
    """Run SLIM designer."""
    try:
        os.makedirs("data/design_slim", exist_ok=True)
        
        with open("data/design_slim/slim_template_gene.fasta", 'w') as f:
            f.write(">template_gene\n")
            f.write(template_gene_content + "\n")
        
        with open("data/design_slim/slim_context.fasta", 'w') as f:
            f.write(">context\n")
            f.write(context_content + "\n")
        
        mutation_lines = [m.strip() for m in mutations_text.split('\n') if m.strip()]
        with open("data/design_slim/slim_target_mutations.csv", 'w') as f:
            f.write("mutations\n")
            for mut in mutation_lines:
                f.write(f"{mut}\n")
        
        logger.info("Running design_slim...")
        design_slim.main()
        
        results_file = "results/design_slim/SLIM_primers.csv"
        if os.path.exists(results_file):
            df = pd.read_csv(results_file)
            result_text = f"Generated primers for {len(df)//4} mutations\n\n"
            result_text += df.to_string(index=False)
            
            log_file = "results/design_slim/design_slim.log"
            if os.path.exists(log_file):
                with open(log_file, 'r') as f:
                    log_content = f.read()
                    result_text += f"\n\n--- Log ---\n{log_content[-1000:]}"
            
            return result_text
        else:
            return "Error: Results file not found"
    except Exception as e:
        logger.exception(f"Error: {e}")
        return f"Error: {str(e)}"

def run_design_gibson(template_gene_content, context_content, mutations_text):
    """Run Gibson designer."""
    try:
        os.makedirs("data/design_gibson", exist_ok=True)
        
        with open("data/design_gibson/gibson_template_gene.fasta", 'w') as f:
            f.write(">template_gene\n")
            f.write(template_gene_content + "\n")
        
        with open("data/design_gibson/gibson_context.fasta", 'w') as f:
            f.write(">context\n")
            f.write(context_content + "\n")
        
        mutation_lines = [m.strip() for m in mutations_text.split('\n') if m.strip()]
        with open("data/design_gibson/gibson_target_mutations.csv", 'w') as f:
            f.write("mutations\n")
            for mut in mutation_lines:
                f.write(f"{mut}\n")
        
        logger.info("Running design_gibson...")
        design_gibson.main()
        
        results_dir = "results/design_gibson"
        if os.path.exists(results_dir):
            result_text = "Gibson assembly plan generated\n\n"
            
            csv_files = list(Path(results_dir).glob("*.csv"))
            if csv_files:
                df = pd.read_csv(csv_files[0])
                result_text += f"Generated {len(df)} primer pairs\n\n"
                result_text += df.head().to_string(index=False)
            
            log_file = "results/design_gibson/design_gibson.log"
            if os.path.exists(log_file):
                with open(log_file, 'r') as f:
                    log_content = f.read()
                    result_text += f"\n\n--- Log ---\n{log_content[-1000:]}"
            
            return result_text
        else:
            return "Error: Results directory not found"
    except Exception as e:
        logger.exception(f"Error: {e}")
        return f"Error: {str(e)}"

def run_mutation_caller(fastq_file, template_file, config_csv_file):
    """Run mutation caller with uploaded files."""
    try:
        os.makedirs("data/mutation_caller", exist_ok=True)
        
        # Copy uploaded files
        if fastq_file:
            shutil.copy(fastq_file, "data/mutation_caller/")
        if template_file:
            shutil.copy(template_file, "data/mutation_caller/mutation_caller_template.fasta")
        if config_csv_file:
            shutil.copy(config_csv_file, "data/mutation_caller/mutation_caller.csv")
        
        logger.info("Running mutation_caller...")
        mutation_caller.main()
        
        results_dir = "results/mutation_caller"
        if os.path.exists(results_dir):
            result_text = "Mutation caller analysis complete\n\n"
            result_text += "Results saved to: results/mutation_caller/\n"
            
            log_file = "results/mutation_caller/mutation_caller.log"
            if os.path.exists(log_file):
                with open(log_file, 'r') as f:
                    log_content = f.read()
                    result_text += f"\n--- Log ---\n{log_content[-1000:]}"
            
            return result_text
        else:
            return "Error: Results directory not found"
    except Exception as e:
        logger.exception(f"Error: {e}")
        return f"Error: {str(e)}"

def run_umi_hunter(fastq_file, template_file, config_csv_file):
    """Run UMI hunter with uploaded files."""
    try:
        os.makedirs("data/umi_hunter", exist_ok=True)
        
        # Copy uploaded files
        if fastq_file:
            shutil.copy(fastq_file, "data/umi_hunter/")
        if template_file:
            shutil.copy(template_file, "data/umi_hunter/template.fasta")
        if config_csv_file:
            shutil.copy(config_csv_file, "data/umi_hunter/umi_hunter.csv")
        
        logger.info("Running umi_hunter...")
        umi_hunter.main()
        
        results_dir = "results/umi_hunter"
        if os.path.exists(results_dir):
            result_text = "UMI hunter analysis complete\n\n"
            result_text += "Results saved to: results/umi_hunter/\n"
            return result_text
        else:
            return "Error: Results directory not found"
    except Exception as e:
        logger.exception(f"Error: {e}")
        return f"Error: {str(e)}"

def run_ep_library_profile(fastq_files, region_file, plasmid_file):
    """Run ep-library-profile with uploaded files."""
    try:
        os.makedirs("data/ep-library-profile", exist_ok=True)
        
        # Handle fastq file uploads
        if fastq_files:
            for f in fastq_files:
                if f:
                    shutil.copy(f, "data/ep-library-profile/")
        
        if region_file:
            shutil.copy(region_file, "data/ep-library-profile/region_of_interest.fasta")
        if plasmid_file:
            shutil.copy(plasmid_file, "data/ep-library-profile/plasmid.fasta")
        
        logger.info("Running ep-library-profile...")
        mut_rate.main()
        
        results_dir = "results/ep-library-profile"
        if os.path.exists(results_dir):
            # Find the latest result directory
            result_dirs = [d for d in os.listdir(results_dir) if os.path.isdir(os.path.join(results_dir, d))]
            result_dirs.sort(key=lambda x: os.path.getmtime(os.path.join(results_dir, x)), reverse=True)
            
            images = []
            if result_dirs:
                latest_dir = os.path.join(results_dir, result_dirs[0])
                
                # Look for summary_panels.png
                summary_panels = os.path.join(latest_dir, "summary_panels.png")
                if os.path.exists(summary_panels):
                    images.append(summary_panels)
                
                # Look for qc_plot
                qc_files = [f for f in os.listdir(latest_dir) if f.startswith("qc_plot") and f.endswith(".png")]
                if qc_files:
                    qc_plot = os.path.join(latest_dir, qc_files[0])
                    if os.path.exists(qc_plot):
                        images.append(qc_plot)
            
            result_text = "EP library profile analysis complete\n\n"
            result_text += f"Results saved to: {os.path.join(results_dir, result_dirs[0]) if result_dirs else 'results/ep-library-profile/'}\n"
            
            log_file = "results/ep-library-profile/master_summary.txt"
            if os.path.exists(log_file):
                with open(log_file, 'r') as f:
                    log_content = f.read()
                    result_text += f"\n--- Summary ---\n{log_content}"
            
            # Return text and images
            return result_text, images if images else None
        else:
            return "Error: Results directory not found", None
    except Exception as e:
        logger.exception(f"Error: {e}")
        return f"Error: {str(e)}", None

def create_gui():
    """Create the Gradio interface with modern design."""
    
    custom_css = """
    .gradio-container {
        max-width: 1400px;
        margin: 0 auto;
    }
    .header {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 2rem;
        border-radius: 10px;
        margin-bottom: 2rem;
        color: white;
        text-align: center;
    }
    .header h1 {
        margin: 0;
        font-size: 2.5rem;
        font-weight: 600;
    }
    .header p {
        margin: 0.5rem 0 0 0;
        opacity: 0.9;
    }
    """
    
    with gr.Blocks(
        title="uht-tooling",
        theme=gr.themes.Soft(primary_hue="purple"),
        css=custom_css
    ) as demo:
        
        with gr.Column(elem_classes="header"):
            gr.Markdown(
                """
                # uht-tooling
                ### Ultra-High Throughput Screening Tools
                Automated molecular biology pipeline for primer design and mutation analysis
                """,
                elem_classes="header"
            )
        
        with gr.Tabs() as tabs:
            
            # Tab 1: Nextera Designer
            with gr.TabItem("Nextera Primer Design"):
                with gr.Row():
                    gr.Markdown("""
                    ### Illumina-Compatible Primer Design
                    
                    Generate Illumina-compatible primers with pre-loaded i5/i7 indices.
                    Supports up to 144 unique amplicons with twelve i5 and twelve i7 indices.
                    
                    **Enter primer sequences in 5'→3' orientation**
                    """)
                
                with gr.Row():
                    with gr.Column(scale=1):
                        forward_primer = gr.Textbox(
                            label="Forward Primer Sequence",
                            placeholder="cgcgtcatgcatgcgactcact",
                            max_lines=1
                        )
                        reverse_primer = gr.Textbox(
                            label="Reverse Primer Sequence", 
                            placeholder="gactaccggaagcagtgtgacc",
                            max_lines=1
                        )
                        nextera_btn = gr.Button("Generate Primers", variant="primary", scale=1)
                
                nextera_output = gr.Textbox(
                    label="Generated Primers",
                    lines=25,
                    show_copy_button=True
                )
                
                nextera_btn.click(
                    fn=run_nextera_designer,
                    inputs=[forward_primer, reverse_primer],
                    outputs=nextera_output
                )
            
            # Tab 2: SLIM Designer
            with gr.TabItem("SLIM Designer"):
                with gr.Row():
                    gr.Markdown("""
                    ### SLIM Cloning Primer Design
                    
                    Design primers for site-directed mutagenesis using SLIM (Sequence-Ligation Independent Mutagenesis).
                    
                    **Mutation Format Examples:**
                    - Substitution: `A123G`
                    - Deletion: `T241Del`
                    - Insertion: `T241TS`
                    """)
                
                with gr.Row():
                    with gr.Column(scale=1):
                        template_gene = gr.Textbox(
                            label="Template Gene Sequence (FASTA)",
                            placeholder="ATGGGC...",
                            lines=4
                        )
                        context_seq = gr.Textbox(
                            label="Full Plasmid Context (FASTA)",
                            placeholder="ATGGCG...",
                            lines=4
                        )
                        mutations = gr.Textbox(
                            label="Mutations (one per line)",
                            placeholder="A123G\nL67V\nW165C",
                            lines=6
                        )
                        slim_btn = gr.Button("Design SLIM Primers", variant="primary")
                
                slim_output = gr.Textbox(
                    label="Designed Primers",
                    lines=25,
                    show_copy_button=True
                )
                
                slim_btn.click(
                    fn=run_design_slim,
                    inputs=[template_gene, context_seq, mutations],
                    outputs=slim_output
                )
            
            # Tab 3: Gibson Assembly
            with gr.TabItem("Gibson Assembly"):
                with gr.Row():
                    gr.Markdown("""
                    ### Gibson Assembly Primer Design
                    
                    Design primers for Gibson assembly mutagenesis with support for multi-mutations.
                    Use `+` to combine multiple mutations (e.g., `A123G+T241Del`).
                    """)
                
                with gr.Row():
                    with gr.Column(scale=1):
                        gibson_template = gr.Textbox(
                            label="Template Gene Sequence (FASTA)",
                            placeholder="ATGGGC...",
                            lines=4
                        )
                        gibson_context = gr.Textbox(
                            label="Full Plasmid Context (FASTA)",
                            placeholder="ATGGCG...",
                            lines=4
                        )
                        gibson_mutations = gr.Textbox(
                            label="Mutations (use + for multi-mutations)",
                            placeholder="A123G\nL67V\nA123G+T241Del",
                            lines=6
                        )
                        gibson_btn = gr.Button("Design Gibson Primers", variant="primary")
                
                gibson_output = gr.Textbox(
                    label="Designed Primers",
                    lines=25,
                    show_copy_button=True
                )
                
                gibson_btn.click(
                    fn=run_design_gibson,
                    inputs=[gibson_template, gibson_context, gibson_mutations],
                    outputs=gibson_output
                )
            
            # Tab 4: Mutation Caller
            with gr.TabItem("Mutation Caller"):
                with gr.Row():
                    gr.Markdown("""
                    ### Mutation Analysis from Long-Read Data
                    
                    Identify and count mutations from long-read sequencing data without UMIs.
                    Provides mutation counts and co-occurrence data.
                    """)
                
                with gr.Row():
                    with gr.Column(scale=1):
                        mc_fastq = gr.File(
                            label="FASTQ File (.fastq.gz)",
                            file_types=[".fastq", ".gz"],
                            type="filepath"
                        )
                        mc_template = gr.File(
                            label="Template Gene (FASTA)",
                            file_types=[".fasta", ".fa"],
                            type="filepath"
                        )
                        mc_config = gr.File(
                            label="Configuration CSV",
                            file_types=[".csv"],
                            type="filepath"
                        )
                        mc_btn = gr.Button("Run Mutation Caller", variant="primary")
                
                mc_output = gr.Textbox(
                    label="Analysis Results",
                    lines=20,
                    show_copy_button=True
                )
                
                mc_btn.click(
                    fn=run_mutation_caller,
                    inputs=[mc_fastq, mc_template, mc_config],
                    outputs=mc_output
                )
            
            # Tab 5: UMI Hunter
            with gr.TabItem("UMI Hunter"):
                with gr.Row():
                    gr.Markdown("""
                    ### UMI-Gene Pairing Analysis
                    
                    Identify, count, and generate consensus UMI-gene pairings from long-read sequencing libraries.
                    Outputs UMI-gene clusters, counts, and consensus sequences.
                    """)
                
                with gr.Row():
                    with gr.Column(scale=1):
                        umi_fastq = gr.File(
                            label="FASTQ File (.fastq.gz)",
                            file_types=[".fastq", ".gz"],
                            type="filepath"
                        )
                        umi_template = gr.File(
                            label="Template Gene (FASTA)",
                            file_types=[".fasta", ".fa"],
                            type="filepath"
                        )
                        umi_config = gr.File(
                            label="Configuration CSV",
                            file_types=[".csv"],
                            type="filepath"
                        )
                        umi_btn = gr.Button("Run UMI Hunter", variant="primary")
                
                umi_output = gr.Textbox(
                    label="Analysis Results",
                    lines=20,
                    show_copy_button=True
                )
                
                umi_btn.click(
                    fn=run_umi_hunter,
                    inputs=[umi_fastq, umi_template, umi_config],
                    outputs=umi_output
                )
            
            # Tab 6: EP Library Profile
            with gr.TabItem("EP Library Profile"):
                with gr.Row():
                    gr.Markdown("""
                    ### Library Profiling Without UMIs
                    
                    Profile DNA libraries to calculate coverage, mutation rates, mutation spectrum,
                    and associated error rates from sequencing data.
                    """)
                
                with gr.Row():
                    with gr.Column(scale=1):
                        ep_fastq = gr.File(
                            label="FASTQ Files (.fastq.gz)",
                            file_types=[".fastq", ".gz"],
                            file_count="multiple",
                            type="filepath"
                        )
                        ep_region = gr.File(
                            label="Region of Interest (FASTA)",
                            file_types=[".fasta", ".fa"],
                            type="filepath"
                        )
                        ep_plasmid = gr.File(
                            label="Plasmid Reference (FASTA)",
                            file_types=[".fasta", ".fa"],
                            type="filepath"
                        )
                        ep_btn = gr.Button("Run EP Library Profile", variant="primary")
                
                ep_output = gr.Textbox(
                    label="Analysis Results",
                    lines=10,
                    show_copy_button=True
                )
                
                ep_images = gr.Gallery(
                    label="Result Plots",
                    show_label=True,
                    elem_id="gallery",
                    columns=2,
                    rows=1,
                    height="auto"
                )
                
                ep_btn.click(
                    fn=run_ep_library_profile,
                    inputs=[ep_fastq, ep_region, ep_plasmid],
                    outputs=[ep_output, ep_images]
                )
            
        
        gr.Markdown(
            """
            ---
            
            **About this interface:**
            
            - This GUI provides a user-friendly wrapper for the command-line tools
            - All existing `make` commands continue to work
            - Results are saved to the `results/` directory
            - Input files are never modified or deleted
            - For batch processing and automation, use the command-line interface
            
            **Keyboard shortcuts:** Press `Ctrl+C` in terminal to stop the GUI
            """
        )
    
    return demo

def main():
    """Launch the GUI."""
    logger.info("Starting uht-tooling GUI...")
    
    demo = create_gui()
    
    logger.info("Launching GUI on http://127.0.0.1:7860")
    demo.launch(
        server_name="127.0.0.1",
        server_port=7860,
        share=False,
        show_error=True
    )

if __name__ == "__main__":
    main()
