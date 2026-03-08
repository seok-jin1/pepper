import pandas as pd
from Bio import SeqIO

def identify_top_asvs():
    print("Identifying Top 5 Dominant ASVs in Spaceflight...")
    
    # 1. Load ASV Table
    table = pd.read_csv('version-2_integrated/exported_table_space/feature-table.tsv', sep='	', skiprows=1, index_col=0)
    
    # Calculate Mean Relative Abundance
    rel_abun = table.div(table.sum(axis=0), axis=1) * 100
    mean_abun = rel_abun.mean(axis=1).sort_values(ascending=False)
    
    top5_ids = mean_abun.head(5).index.tolist()
    
    # 2. Load Taxonomy
    tax = pd.read_csv('version-2_integrated/exported_taxonomy/taxonomy.tsv', sep='	', index_col=0)
    
    # 3. Load Sequences
    seq_dict = SeqIO.to_dict(SeqIO.parse("version-2_integrated/exported_sequences/dna-sequences.fasta", "fasta"))
    
    # 4. Report
    results = []
    print("\n--- TOP 5 DOMINANT ASVs IN SPACE ---")
    for i, asv_id in enumerate(top5_ids):
        abundance = mean_abun[asv_id]
        taxonomy = tax.loc[asv_id, 'Taxon'] if asv_id in tax.index else "Unknown"
        sequence = str(seq_dict[asv_id].seq) if asv_id in seq_dict else "Sequence Not Found"
        
        # Simple species hint based on taxonomy string
        species_hint = taxonomy.split(';')[-1].strip()
        
        print(f"\nRank {i+1}: {asv_id}")
        print(f"  Mean Abundance: {abundance:.2f}%")
        print(f"  Taxonomy: {taxonomy}")
        print(f"  Sequence (V4): {sequence}")
        
        results.append({
            'Rank': i+1,
            'ASV_ID': asv_id,
            'Abundance': f"{abundance:.2f}%",
            'Taxonomy': taxonomy,
            'Sequence': sequence
        })
        
    # Save to CSV for reference
    pd.DataFrame(results).to_csv("version-2_integrated/Top5_Space_ASVs.csv", index=False)
    print("\nIdentification completed. Results saved to Top5_Space_ASVs.csv")

if __name__ == "__main__":
    identify_top_asvs()
