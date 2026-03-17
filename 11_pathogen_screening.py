import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def screen_pathogens():
    print("Screening for potential plant pathogens...")
    # 1. Load Genus level table (reusing logic from correlation script)
    table = pd.read_csv('version-2_integrated/exported_table_clean/feature-table.tsv', sep='	', skiprows=1, index_col=0)
    tax = pd.read_csv('version-2_integrated/exported_taxonomy/taxonomy.tsv', sep='	', index_col=0)
    meta = pd.read_csv('version-2_integrated/integrated_metadata.tsv', sep='	', index_col=0)
    
    def get_genus(taxon_str):
        if pd.isna(taxon_str): return 'Unassigned'
        parts = taxon_str.split(';')
        for part in reversed(parts):
            if 'g__' in part and len(part.split('__')) > 1:
                name = part.split('__')[1].strip()
                return name if name != "" else "Unassigned"
        return 'Unassigned'
    
    tax['Genus'] = tax['Taxon'].apply(get_genus)
    table_tax = table.join(tax[['Genus']], how='inner')
    genus_table = table_tax.groupby('Genus').sum()
    genus_rel = genus_table.div(genus_table.sum(axis=0), axis=1) * 100
    
    # 2. Target Pathogen List (Common plant pathogens)
    target_pathogens = [
        'Ralstonia', 'Xanthomonas', 'Erwinia', 'Pectobacterium', 
        'Dickeya', 'Agrobacterium', 'Clavibacter', 'Pantoea', 'Pseudomonas'
    ]
    # Note: Pseudomonas can be beneficial or pathogenic, but we include it for context.
    
    pathogen_data = genus_rel.loc[genus_rel.index.intersection(target_pathogens)].T
    pathogen_data = pathogen_data.join(meta[['Study_Group']], how='inner')
    
    # 3. Visualization
    melted = pathogen_data.melt(id_vars=['Study_Group'], var_name='Pathogen_Genus', value_name='Relative Abundance (%)')
    # Filter for relevant study groups
    melted = melted[melted['Study_Group'].isin(['Space_Flight', 'Terrestrial_Soil'])]
    
    plt.figure(figsize=(14, 8))
    sns.boxplot(x='Pathogen_Genus', y='Relative Abundance (%)', hue='Study_Group', data=melted, palette=['#3498db', '#2ecc71'])
    plt.yscale('log') # 병원균은 보통 소수이므로 로그 스케일 권장
    plt.title('Potential Plant Pathogen Screening (Space vs Terrestrial)', fontsize=15, fontweight='bold')
    plt.ylabel('Relative Abundance (%) - Log Scale')
    plt.xlabel('Bacterial Genus')
    plt.xticks(rotation=45)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig('version-2_integrated/Pathogen_Screening.png', dpi=300)
    
    # Summary Statistics
    summary = melted.groupby(['Pathogen_Genus', 'Study_Group'])['Relative Abundance (%)'].mean().unstack()
    print("\n--- Mean Relative Abundance (%) of Potential Pathogens ---")
    print(summary)
    
    print("\nPathogen screening completed. Pathogen_Screening.png generated.")

if __name__ == "__main__":
    try:
        screen_pathogens()
    except Exception as e:
        print(f"Error during screening: {e}")
