#!/usr/bin/env python
"""
pubmed_noise_index_threshold_sweep.py
Extended PubMed search analysis:
1) Calculate noise index for each keyword.
2) Sweep noise threshold (1–100%) to refine noisy keywords.
3) Compare baseline vs refined total results and plot efficiency curve.
4) Generate all keyword × filter combinations with review-only restriction.
"""

from Bio import Entrez
import pandas as pd
import argparse
import time
import matplotlib.pyplot as plt

# REQUIRED: Set your email for NCBI Entrez
Entrez.email = "your_email@example.com"
# Optional: Add your NCBI API key for higher rate limits
# Entrez.api_key = "your_ncbi_api_key"

# Default files
DEFAULT_KEYWORDS_FILE = "keywords.txt"
DEFAULT_FILTERS_FILE = "filters.txt"
DEFAULT_REVIEW_FILE = "review_filter.txt"

# -----------------------------------------
def load_terms(file_path):
    with open(file_path, "r") as f:
        return [line.strip() for line in f if line.strip()]

def get_pubmed_count(query, retries=3, delay=0.5):
    """Return the number of PubMed hits for a query with retries."""
    for attempt in range(retries):
        try:
            handle = Entrez.esearch(db="pubmed", term=query, retmode="xml")
            record = Entrez.read(handle)
            handle.close()
            return int(record["Count"])
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(delay * (attempt + 1))  # exponential backoff
            else:
                print(f"Failed query: {query} | Error: {e}")
                return 0

# -----------------------------------------
def analyze_noise(keywords, health_filter, noise_threshold=50):
    results = []
    refined_queries = []

    for kw in keywords:
        base = get_pubmed_count(kw)
        filt = get_pubmed_count(f"{kw} AND {health_filter}")

        noise_index = (base - filt) / base * 100 if base > 0 else 0

        if noise_index >= noise_threshold:
            refined = f"{kw} AND {health_filter}"
            refined_queries.append(refined)
            note = "Refined"
        else:
            refined_queries.append(kw)
            note = "Kept broad"

        results.append({
            "Keyword": kw,
            "Base Hits": base,
            "With Health Filter": filt,
            "Noise Index %": round(noise_index, 1),
            "Decision": note
        })

    return pd.DataFrame(results), refined_queries

def compare_efficiency(original_keywords, refined_queries):
    base_total = sum(get_pubmed_count(kw) for kw in original_keywords)
    refined_total = sum(get_pubmed_count(kw) for kw in refined_queries)

    efficiency_gain = (base_total - refined_total) / base_total * 100 if base_total > 0 else 0

    return {
        "Original Total Hits": base_total,
        "Refined Total Hits": refined_total,
        "Efficiency Gain %": round(efficiency_gain, 1)
    }

# -----------------------------------------
def run_combinations(keywords, filters, review_filter, out_file="combination_review_hits.csv"):
    results = []
    health_filter = " OR ".join([f"{flt}[tiab]" for flt in filters])

    for kw in keywords:
        query = f'({kw}[tiab]) AND ({health_filter}) AND {review_filter}'
        hits = get_pubmed_count(query)
        results.append({
            "Keyword": kw,
            "Query": query,
            "Review Hits": hits
        })
        print(f"{query} → {hits} hits")


    df = pd.DataFrame(results)
    df.to_csv(out_file, index=False)
    print(f"\nCombination results saved to {out_file}")
    return df

# -----------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PubMed Noise Index Analysis with Threshold Sweep and Combination Search")
    parser.add_argument("--keywords", default=DEFAULT_KEYWORDS_FILE, help="Path to keywords.txt")
    parser.add_argument("--filters", default=DEFAULT_FILTERS_FILE, help="Path to filters.txt")
    parser.add_argument("--review", default=DEFAULT_REVIEW_FILE, help="Path to review_filter.txt")
    parser.add_argument("--out", default="noise_report.csv", help="CSV output file for noise analysis")
    args = parser.parse_args()

    keywords = load_terms(args.keywords)
    filters = load_terms(args.filters)
    review_filter = " ".join(load_terms(args.review))

    # Build health filter for noise analysis
    health_terms = " OR ".join(filters)
    health_filter = f"({health_terms})"

    # -------- Stage 1: Noise Analysis + Threshold Sweep --------
    sweep_results = []
    for threshold in range(10, 101, 10):  # increments of 10
        df, refined_queries = analyze_noise(keywords, health_filter, threshold)
        report = compare_efficiency(keywords, refined_queries)
        sweep_results.append({
            "Threshold %": threshold,
            "Original Total Hits": report["Original Total Hits"],
            "Refined Total Hits": report["Refined Total Hits"],
            "Efficiency Gain %": report["Efficiency Gain %"]
        })
        print(f"Threshold {threshold}% → Efficiency Gain: {report['Efficiency Gain %']}%")

    sweep_df = pd.DataFrame(sweep_results)
    sweep_df.to_csv("threshold_sweep_results.csv", index=False)
    print("\nThreshold sweep results saved to threshold_sweep_results.csv")

    # Plot efficiency curve
    plt.figure(figsize=(8, 5))
    plt.plot(sweep_df["Threshold %"], sweep_df["Efficiency Gain %"], marker="o")
    plt.xlabel("Noise Threshold (%)")
    plt.ylabel("Efficiency Gain (%)")
    plt.title("Noise Threshold vs Efficiency Gain in PubMed Search")
    plt.grid(True)
    plt.savefig("threshold_efficiency_curve.png", dpi=300)
    print("\nEfficiency curve saved to threshold_efficiency_curve.png")

    # -------- Stage 2: Keyword × Filter Combination Search --------
    run_combinations(keywords, filters, review_filter, out_file="combination_review_hits.csv")
