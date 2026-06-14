# 📊 duoqc 

`duoqc` is a fast, lightweight, and interactive paired-end FASTQ Quality Control pipeline engineered entirely from scratch in Python. It features a live-rendering split-screen terminal UI, adaptive file decompression, and a MultiQC-style multi-sample batch aggregator.

## ✨ Features
* **Zero-Config Auto-Discovery:** Automatically scans your working directory and pairs matching `R1` and `R2` sequencing files.
* **Hybrid Compression Engine:** Streams raw `.fastq` or compressed `.fastq.gz` datasets natively without extra dependencies.
* **Staggered Block/Stride Sampling:** Samples data chunks ultra-efficiently to maintain a near-zero memory footprint.
* **Live Quality Maps:** Dynamically computes and displays Unicode positional quality profiles (`▆▇▇▆▆`) across read cycles.
* **MultiQC Reports:** Aggregates multi-sample batch runs into a single, beautifully formatted summary table.

## 🚀 Installation & Usage

1. Clone the repository:
   $ git clone https://github.com/SayaGarud/duoqc.git
   $ cd duoqc

2. Install the required interface library:
   $ pip install -r requirements.txt

3. Run the pipeline (Autodetects any sequencing files sitting in the directory):
   $ python duoqc.py

4. Pass custom performance parameters via the CLI:
   $ python duoqc.py -b 100 -s 500 -o giant_stride_report.txt

## 🛠️ Command-Line Interface Flags

| Flag | Full Argument | Type | Default | Description |
|------|---------------|------|---------|-------------|
| -1   | --forward     | STR  | None    | Explicit path to Forward Read file (R1) |
| -2   | --reverse     | STR  | None    | Explicit path to Reverse Read file (R2) |
| -b   | --block       | INT  | 50      | Number of reads to analyze per batch block |
| -s   | --stride      | INT  | 200     | Number of records to skip over during strides |
| -o   | --output      | STR  | duofq_multiqc_summary.txt | Target path for the final master batch report |

## 📜 License
Distributed under the MIT License. See the repository page for details.
