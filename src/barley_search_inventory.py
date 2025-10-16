import pandas as pd
from Bio import SeqIO

# === CONFIGURATION ===
protein_file = r'E:\Guido\sibel\import\ENB Barley import\storage_protein_inventory_template_OG_24_07_deconvoluted.xlsx'
sheet_name = 'Barley Inventory'  # Excel sheet name
fasta_file = r'E:\Guido\sibel\Database_Search_Ouput\ENB Barley\GCF_904849725.1_Hordeum_vulgare_combined_FASTA.faa'
output_fasta = "matched_sequences.fasta"

# === STEP 1: Read protein identifiers from the specified sheet in Excel ===
df = pd.read_excel(protein_file, sheet_name=sheet_name)
protein_ids = set(df['Protein'].astype(str).str.strip())

# === STEP 2: Parse the FASTA file and extract matching sequences ===
matched_records = []

for record in SeqIO.parse(fasta_file, "fasta"):
    fasta_id = record.id.split()[0]  # Get only the main ID before any spaces
    if fasta_id in protein_ids:
        matched_records.append(record)

# === STEP 3: Write matched sequences to new FASTA file ===
SeqIO.write(matched_records, output_fasta, "fasta")

print(f"✅ Done! {len(matched_records)} matching sequences written to '{output_fasta}'.")
