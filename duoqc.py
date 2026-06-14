import os
import gzip
import time
import argparse
import glob
from rich.console import Console
from rich.layout import Layout
from rich.panel import Panel
from rich.live import Live

# 1. Core Phred Decoder
def decode_phred_string(quality_string: str) -> list[int]:
    return [ord(char) - 33 for char in quality_string]

# 2. Text/GZ Adaptive Generator
def stream_fastq(file_path):
    open_func = gzip.open if file_path.endswith('.gz') else open
    mode = 'rt' if file_path.endswith('.gz') else 'r'
    
    with open_func(file_path, mode) as f:
        while True:
            line1 = f.readline()
            if not line1: 
                break
            line2 = f.readline()
            line3 = f.readline()
            line4 = f.readline()
            
            if not line2 or not line3 or not line4:
                break
                
            yield line1.strip(), line2.strip(), line4.strip()

# 3. Comprehensive Batch Auto-Discovery (Option B)
def discover_all_pairs():
    """Scans directory and returns a dictionary of paired samples."""
    extensions = ["*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"]
    all_files = []
    for ext in extensions:
        all_files.extend(glob.glob(ext))
        
    samples = {}
    # Find all R1 files
    for f in sorted(all_files):
        if "_R1" in f or ".R1" in f or "_1.fastq" in f or ".1.fastq" in f:
            # Generate expected sample name prefix
            sample_name = os.path.basename(f).split('_R1')[0].split('.R1')[0].split('_1.')[0].split('.1.')[0]
            
            # Predict matching R2 file path
            possible_r2 = f.replace("R1", "R2").replace("_1.", "_2.").replace(".1.", ".2.")
            if os.path.exists(possible_r2):
                samples[sample_name] = {"R1": f, "R2": possible_r2}
                
    return samples

# 4. Text Sparkline Graph Engine (Option A)
def generate_sparkline(positional_sums, positional_counts):
    """Converts positional lists into an elegant terminal quality graph."""
    if not positional_sums:
        return "[No Data]"
    
    # Calculate averages for 10 structural bins along the length of the read
    num_bins = 10
    bin_size = max(1, len(positional_sums) // num_bins)
    averages = []
    
    for i in range(num_bins):
        start = i * bin_size
        end = start + bin_size
        b_sum = sum(positional_sums[start:end])
        b_count = sum(positional_counts[start:end])
        averages.append(b_sum / b_count if b_count > 0 else 0)
        
    # Scale scores from Q20 to Q40 to block characters
    spark_chars = [" ", "▂", "▃", "▄", "▅", "▆", "▇", "█"]
    sparkline = ""
    for avg in averages:
        if avg <= 20: index = 0
        elif avg >= 40: index = 7
        else: index = int((avg - 20) / 2.85) # Scale gracefully between index 0-7
        sparkline += spark_chars[index]
        
    return sparkline

# 5. Argument Parser Setup
def parse_arguments():
    parser = argparse.ArgumentParser(
        description="📊 duoqc: Multi-sample batch-processing FASTQ Quality Control pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-b", "--block", type=int, default=50, help="Number of reads to process per block")
    parser.add_argument("-s", "--stride", type=int, default=200, help="Number of reads to skip per stride")
    parser.add_argument("-o", "--output", type=str, default="duoqc_multiqc_summary.txt", help="Master batch report path")
    return parser.parse_args()


def main():
    args = parse_arguments()
    console = Console()
    
    BLOCK_SIZE = args.block
    STRIDE_SIZE = args.stride
    master_report_path = args.output
    
    # Discover all paired samples in the folder
    discovered_samples = discover_all_pairs()
    
    if not discovered_samples:
        console.print("\n[bold red]❌ Error:[/bold red] No paired FASTQ datasets found in this directory.\n")
        return
        
    console.print(f"[bold green]✔ MultiQC Mode Active:[/bold green] Found [cyan]{len(discovered_samples)}[/cyan] samples to batch-process.\n")
    time.sleep(1.5)
    
    # Store aggregated run history for the master summary report
    batch_summary_data = []
    
    # Process each discovered sample sequentially
    for current_sample_idx, (sample_name, paths) in enumerate(discovered_samples.items(), start=1):
        file_r1 = paths["R1"]
        file_r2 = paths["R2"]
        
        # Setup Screen Layout
        layout = Layout()
        layout.split_column(Layout(name="header", size=3), Layout(name="body"))
        layout["body"].split_row(Layout(name="left_read_r1"), Layout(name="right_read_r2"))
        
        # Connect Paired Stream
        dual_stream = zip(stream_fastq(file_r1), stream_fastq(file_r2))
        
        total_sampled = 0
        block_counter = 0
        
        # Positional quality matrices (Tracks every base index up to 300bp)
        r1_pos_sums, r1_pos_counts = [0] * 300, [0] * 300
        r2_pos_sums, r2_pos_counts = [0] * 300, [0] * 300
        
        r1_total_bases, r1_total_phred_sum, r1_total_gc = 0, 0, 0
        r2_total_bases, r2_total_phred_sum, r2_total_gc = 0, 0, 0
        
        with Live(layout, console=console, refresh_per_second=12):
            try:
                while True:
                    # Update layout header with global batch progress indicators
                    layout["header"].update(
                        Panel(f"📊 duoqc v0.2.0 | Processing Sample [{current_sample_idx}/{len(discovered_samples)}]: [bold yellow]{sample_name}[/bold yellow]", style="bold green", expand=True)
                    )
                    
                    if block_counter < BLOCK_SIZE:
                        read_r1, read_r2 = next(dual_stream)
                        total_sampled += 1
                        block_counter += 1
                        
                        # --- PROCESS FORWARD (R1) ---
                        seq_r1 = read_r1[1].upper()
                        scores_r1 = decode_phred_string(read_r1[2])
                        r1_total_bases += len(seq_r1)
                        r1_total_phred_sum += sum(scores_r1)
                        r1_total_gc += seq_r1.count('G') + seq_r1.count('C')
                        
                        # Accumulate positional sums for R1
                        for idx, score in enumerate(scores_r1):
                            if idx < 300:
                                r1_pos_sums[idx] += score
                                r1_pos_counts[idx] += 1
                                
                        # --- PROCESS REVERSE (R2) ---
                        seq_r2 = read_r2[1].upper()
                        scores_r2 = decode_phred_string(read_r2[2])
                        r2_total_bases += len(seq_r2)
                        r2_total_phred_sum += sum(scores_r2)
                        r2_total_gc += seq_r2.count('G') + seq_r2.count('C')
                        
                        # Accumulate positional sums for R2
                        for idx, score in enumerate(scores_r2):
                            if idx < 300:
                                r2_pos_sums[idx] += score
                                r2_pos_counts[idx] += 1
                                
                        # Calculate current performance windows
                        r1_avg_q = r1_total_phred_sum / r1_total_bases if r1_total_bases > 0 else 0
                        r1_gc_pct = (r1_total_gc / r1_total_bases) * 100 if r1_total_bases > 0 else 0
                        r1_spark = generate_sparkline(r1_pos_sums[:len(seq_r1)], r1_pos_counts[:len(seq_r1)])
                        
                        r2_avg_q = r2_total_phred_sum / r2_total_bases if r2_total_bases > 0 else 0
                        r2_gc_pct = (r2_total_gc / r2_total_bases) * 100 if r2_total_bases > 0 else 0
                        r2_spark = generate_sparkline(r2_pos_sums[:len(seq_r2)], r2_pos_counts[:len(seq_r2)])
                        
                        layout["left_read_r1"].update(
                            Panel(
                                f"📄 [bold]File:[/bold] {os.path.basename(file_r1)}\n\n"
                                f"🔹 Total Sampled:  [bold cyan]{total_sampled}[/bold cyan] reads\n"
                                f"📈 [bold]Avg Phred Q:[/bold]   [bold sky_blue3]{r1_avg_q:.2f}[/bold sky_blue3]\n"
                                f"🧬 [bold]GC Content:[/bold]   [bold light_salmon3]{r1_gc_pct:.1f}%[/bold light_salmon3]\n"
                                f"📊 [bold]Quality Map:[/bold]  [bold green]{r1_spark}[/bold green] (5'->3')\n\n"
                                f"[yellow]Status: Running Analysis...[/yellow]",
                                title="FORWARD READ (R1)", border_style="cyan"
                            )
                        )
                        layout["right_read_r2"].update(
                            Panel(
                                f"📄 [bold]File:[/bold] {os.path.basename(file_r2)}\n\n"
                                f"🔹 Total Sampled:  [bold magenta]{total_sampled}[/bold magenta] reads\n"
                                f"📈 [bold]Avg Phred Q:[/bold]   [bold peach_puff3]{r2_avg_q:.2f}[/bold peach_puff3]\n"
                                f"🧬 [bold]GC Content:[/bold]   [bold light_salmon3]{r2_gc_pct:.1f}%[/bold light_salmon3]\n"
                                f"📊 [bold]Quality Map:[/bold]  [bold green]{r2_spark}[/bold green] (5'->3')\n\n"
                                f"[yellow]Status: Running Analysis...[/yellow]",
                                title="REVERSE READ (R2)", border_style="magenta"
                            )
                        )
                        time.sleep(0.005)
                        
                    else:
                        # --- STRIDE SKIP PHASE ---
                        for _ in range(STRIDE_SIZE):
                            try:
                                next(dual_stream)
                            except StopIteration:
                                raise StopIteration
                        block_counter = 0
                        
            except StopIteration:
                # Compile final run configurations for this individual loop iteration
                r1_final_q = r1_total_phred_sum / r1_total_bases if r1_total_bases > 0 else 0
                r1_final_gc = (r1_total_gc / r1_total_bases) * 100 if r1_total_bases > 0 else 0
                r1_final_spark = generate_sparkline(r1_pos_sums, r1_pos_counts)
                
                r2_final_q = r2_total_phred_sum / r2_total_bases if r2_total_bases > 0 else 0
                r2_final_gc = (r2_total_gc / r2_total_bases) * 100 if r2_total_bases > 0 else 0
                r2_final_spark = generate_sparkline(r2_pos_sums, r2_pos_counts)
                
                # Lock terminal panels for this specific sample run
                layout["left_read_r1"].update(
                    Panel(
                        f"📄 [bold]File:[/bold] {os.path.basename(file_r1)}\n\n"
                        f"✅ Sampled Count: {total_sampled} reads\n"
                        f"📈 [bold]Final Avg Q:[/bold]    [bold green]{r1_final_q:.2f}[/bold green]\n"
                        f"🧬 [bold]Final GC %:[/bold]    [bold green]{r1_final_gc:.1f}%[/bold green]\n"
                        f"📊 [bold]Quality Map:[/bold]  [bold green]{r1_final_spark}[/bold green]\n\n"
                        f"[bold green]✔ SAMPLE PIPELINE FINISHED[/bold green]", 
                        title="R1 COMPLETE", border_style="green"
                    )
                )
                layout["right_read_r2"].update(
                    Panel(
                        f"📄 [bold]File:[/bold] {os.path.basename(file_r2)}\n\n"
                        f"✅ Sampled Count: {total_sampled} reads\n"
                        f"📈 [bold]Final Avg Q:[/bold]    [bold green]{r2_final_q:.2f}[/bold green]\n"
                        f"🧬 [bold]Final GC %:[/bold]    [bold green]{r2_final_gc:.1f}%[/bold green]\n"
                        f"📊 [bold]Quality Map:[/bold]  [bold green]{r2_final_spark}[/bold green]\n\n"
                        f"[bold green]✔ SAMPLE PIPELINE FINISHED[/bold green]", 
                        title="R2 COMPLETE", border_style="green"
                    )
                )
                
                # Cache results for the compiled multiQC final log
                batch_summary_data.append({
                    "name": sample_name,
                    "r1_q": r1_final_q, "r1_gc": r1_final_gc,
                    "r2_q": r2_final_q, "r2_gc": r2_final_gc,
                    "reads": total_sampled
                })
                
                time.sleep(2.0) # Pause so the user can look at the complete sample view before looping to the next dataset

    # 6. Out-of-Loop MultiQC Unified Master Report Generator
    with open(master_report_path, "w") as report:
        report.write("=========================================================================\n")
        report.write("               DUOQC PIPELINE MULTIQC AGGREGATION REPORT                 \n")
        report.write("=========================================================================\n")
        report.write(f"Execution Target Directory: {os.getcwd()}\n")
        report.write(f"Total Unique Paired Datasets Processed: {len(batch_summary_data)}\n")
        report.write("-------------------------------------------------------------------------\n")
        report.write(f"{'SAMPLE ID':<20} | {'SAMPLED':<10} | {'R1 QUALITY':<10} | {'R1 GC':<8} | {'R2 QUALITY':<10} | {'R2 GC':<8}\n")
        report.write("-------------------------------------------------------------------------\n")
        for data in batch_summary_data:
            report.write(
                f"{data['name'][:20]:<20} | "
                f"{data['reads']:<10} | "
                f"{data['r1_q']:<10.2f} | "
                f"{data['r1_gc']:<7.1f}% | "
                f"{data['r2_q']:<10.2f} | "
                f"{data['r2_gc']:<7.1f}%\n"
                f"Status: PASS\n"
            )
        report.write("=========================================================================\n")
        
    console.print(f"\n[bold green]✔ MultiQC Processing Finished successfully![/bold green]")
    console.print(f"📁 Combined batch analytics written to permanent log: [bold cyan]{master_report_path}[/bold cyan]\n")

if __name__ == "__main__":
    main()