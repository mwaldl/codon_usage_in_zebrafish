
rule get_codon_usage_per_isoform:
    input:
        fasta = 'raw/gene_sequences/'+ config['gene_isoforms_fasta'],
    output:
        'resources/codon_usage_per_isoform/cu.tsv'
    run:
        from Bio import SeqIO
        from collections import Counter
        import pandas as pd

        data = []

        for record in SeqIO.parse(input.fasta, "fasta"):
            seq = str(record.seq)

            # Test if sequence length is a multiple of 3
            mod3 = (len(seq) % 3 == 0)

            # Split sequence into codons
            codons = [seq[i:i+3] for i in range(0, len(seq), 3)]

            # Count codons
            codon_counts = Counter(codons)

            # Differentiate start codon (ATG) from internal ATG
            first_codon = codons[0]
            if first_codon == "ATG":
                codon_counts['ATG'] -= 1
                codon_counts['iATG'] = 1

            # Differentiate and count stop codons (or missanotated last codons) with a leading 's'
            last_codon = codons[-1]
            codon_counts[f's{last_codon}'] = 1
            codon_counts[last_codon] -= 1

            # Parse gene information
            isoform = record.id.split('.')[0]
            gene = next((part.split('gene:')[1].split('.')[0] for part in record.description.split() if part.startswith('gene:')), 'unknown')

            # Prepare the record dictionary
            record_dict = {
                'isoform': isoform,
                'gene': gene,
                'length is multiple of 3': mod3,
                'nt length': len(seq)
            }

            # Update record dictionary with codon counts
            record_dict.update(codon_counts)
            data.append(record_dict)

        # Create DataFrame and write to TSV
        df = pd.DataFrame(data)
        df.to_csv(str(output), sep='\t', index=False)


rule include_detailed_appris_annotation_per_isoform:
    input:
        codon_usage = 'resources/codon_usage_per_isoform/cu.tsv',
        appris = 'raw/appris_data/'+config['appris_annotation'],
    output:
        detailed_tsv = 'resources/codon_usage_per_isoform/appris-detailed_cu.tsv',
    run:
        import pandas as pd
        cu_df = pd.read_csv(input.codon_usage, sep='\t', index_col = 'isoform')
        appris_df = pd.read_csv(input.appris, sep='\t', index_col='Transcript ID')
        cu_df=cu_df.join(appris_df, on=None, how='left')
        cu_df.to_csv(str(output.detailed_tsv ), sep='\t')


rule include_appris_annotation_per_isoform:
    input:
        codon_usage = 'resources/codon_usage_per_isoform/cu.tsv',
        appris = 'raw/appris_data/'+config['appris_annotation'],
    output:
        tsv = 'resources/codon_usage_per_isoform/appris_cu.tsv',
    run:
        import pandas as pd

        # Load codon usage data
        codon_usage_df = pd.read_csv(input.codon_usage, sep='\t', index_col='isoform')

        # Load appris data with specified columns
        ##appris_columns = ['Transcript ID', 'Transcript type', 'Not found tag', 'Protein length', 'APPRIS Annotation']
        appris_df = pd.read_csv(input.appris, sep='\t', usecols=[2,5,6,9,19], index_col='Transcript ID')

        # Join the dataframes on the index with a left join
        merged_df = codon_usage_df.join(appris_df, how='left')

        # Fill NaN values in 'APPRIS Annotation' with 'undefined'
        merged_df['APPRIS Annotation'] = merged_df['APPRIS Annotation'].fillna('undefined')

        # Create 'appris principal' column based on 'APPRIS Annotation'
        merged_df['appris principal'] = merged_df['APPRIS Annotation'].str.contains('PRINCIPAL')

        # Save the merged dataframe to a TSV file
        merged_df.to_csv(output.tsv, sep='\t')
