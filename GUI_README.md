# uht-tooling GUI Quick Start

This is an **optional** prototype GUI wrapper for the uht-tooling command-line tools.

## Installation

Install the GUI dependency (already in requirements.txt):
```bash
pip install gradio
```

## Launching the GUI

From the project root directory:
```bash
make gui
```

This will launch a web interface at: **http://127.0.0.1:7860**

## Features

### What the GUI provides:

1. **Nextera Primer Design**
   - Enter forward and reverse primer sequences directly
   - Get Illumina-compatible primers with i5/i7 indices instantly
   - No need to edit CSV files

2. **SLIM Designer** 
   - Enter template gene and context sequences
   - Input mutations in simple format (e.g., `A123G`, `T241Del`)
   - Get primer sequences for SLIM cloning

3. **Gibson Assembly Designer**
   - Design primers for Gibson assembly mutagenesis
   - Supports multi-mutations with `+` syntax

### How It Works

The GUI is a **non-destructive wrapper** that:
- Calls Python scripts under the hood
- Saves temporary input files to `data/` folders
- Generates outputs in `results/` folders (same as command line)
- **Never modifies or deletes existing files**

### Workflow Compatibility

**You can still use everything via command line:**
```bash
# All these still work exactly as before:
make nextera_primers
make design_slim
make design_gibson
make mutation_caller
# etc.
```

**Or use the GUI for interactive work:**
```bash
make gui
```

Both approaches produce identical results!

## Important Notes

### File-based Tools

Tools that require FASTQ files (mutation_caller, umi_hunter, ep-library-profile, profile_inserts) are best used via command line due to large file sizes. The GUI includes documentation on how to use these tools.

### For Best Results

- **Quick tests and interactive design**: Use GUI
- **Batch processing and automation**: Use command line
- **Large files**: Use command line

## Troubleshooting

**Port 7860 already in use?**
The GUI will automatically find an available port.

**Module not found errors?**
Make sure Gradio is installed: `pip install gradio`

**Want to stop the GUI?**
Press `Ctrl+C` in the terminal where you launched it.

## Getting Help

- Check the main README.md for detailed tool documentation
- All output files are in `results/` folders
- Log files are also in `results/` folders for debugging

