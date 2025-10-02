#!/usr/bin/env python
"""
pubmed_noise_index.py
Quantifies noisiness of search keywords in PubMed by comparing raw hits vs hits with health filters.
"""

from Bio import Entrez
import pandas as pd

Entrez.email = "your_email@example.com"  # REQUIRED

# -----------------------------------------
# Keywords from digital patient spectrum
keywords = [
    '"digital twin"[tiab]', '"digital twins"[tiab]', '"digital twinning"[tiab]',
    '"virtual twin"[tiab]', '"virtual twins"[tiab]', '"virtual twinning"[tiab]',
    '"digital replica"[tiab]', '"digital replicas"[tiab]', '"digital replication"[tiab]',
    '"virtual human*"[tiab]', '"virtual physiological human*"[tiab]', 
    '"digital patient*"[tiab]', '"virtual patient*"[tiab]', '"in silico patient*"[tiab]',
    '"synthetic patient*"[tiab]', '"virtual sibling*"[tiab]', 
    '"patient-specific computational model*"[tiab]', '"personalized computational model*"[tiab]',
    '"patient avatar*"[tiab]', '"digital avatar*"[tiab]', '"virtual avatar*"[tiab]',
    '"medical twin*"[tiab]', '"health twin*"[tiab]', '"organ twin*"[tiab]', 
    '"surrogate model*"[tiab]', '"biomimetic twin*"[tiab]'
]

# Health anchor terms
health_terms = '(patient*[tiab] OR clinical[tiab] OR medical[tiab] OR health[tiab] OR medicine[tiab])'

# -----------------------------------------
def get_pubmed_count(query):
    """Return the number of PubMed hits for a query."""
    handle = Entrez.esearch(db="pubmed", term=query, retmode="xml")
    record = Entrez.read(handle)
    handle.close()
    return int(record["Count"])

# -----------------------------------------
results = []
for kw in keywords:
    base_count = get_pubmed_count(kw)
    filtered_count = get_pubmed_count(f"{kw} AND {health_terms}")
    
    # Noise index = % drop in hits
    if base_count > 0:
        noise_index = (base_count - filtered_count) / base_count * 100
    else:
        noise_index = 0
    
    results.append({
        "Keyword": kw,
        "Base Hits": base_count,
        "With Health Filter": filtered_count,
        "Noise Index %": round(noise_index, 1)
    })

# -----------------------------------------
# Display results
import pandas as pd
df = pd.DataFrame(results)
print(df.to_string(index=False))

