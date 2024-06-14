

rule get_read_counts_after_filtering_for_coverage_coding_mt_in_RPM:
    input:
        tsv = 'raw/sequencing_gene_counts/' + config['raw_sequencing_counts']
    output:
        tsv = "resources/normalized_counts/RPM_{cutoff}.tsv",
        log = "resources/normalized_counts/RPM_{cutoff}.txt"
    run:
        import pandas as pd

        # Convert cutoff from str to int
        cutoff = int(wildcards.cutoff)

        # Read the input TSV file
        df = pd.read_csv(input.tsv, sep='\t')

        # Filter rows where GeneType is 'protein_coding'
        df = df[df['GeneType'] == 'protein_coding']

        # Filter out rows where GeneSymbol starts with 'mt-'
        df = df[~df['GeneSymbol'].str.startswith('mt-')]

        # Normalize count columns
        non_gene_columns = ['GeneType', 'GeneSymbol', 'ENSG']
        for c in df.columns:
            if c not in non_gene_columns:
                column_sum = df[c].sum()
                df[c] = df[c] * 1000000 / column_sum

        # Calculate the max value across the row, ignoring non-gene columns
        df['max'] = df.drop(columns=non_gene_columns).max(axis=1)

        # Filter rows/genes based on the cutoff
        df = df[df['max'] >= cutoff]

        # Drop the 'max' column as it's no longer needed
        df.drop(columns=['max'], inplace=True)

        # Write the output log and normalize counts to TPM
        with open(output.log, "w") as f:
            for c in df.columns:
                if c not in non_gene_columns:
                    column_sum = df[c].sum()
                    f.write(f"{c}: {column_sum}\n")
                    df[c] = df[c] * 1000000 / column_sum

            genes = len(df)
            f.write(f"number of genes: {genes}\n")

        # Save the processed DataFrame to a new TSV file
        df.to_csv(str(output.tsv), sep='\t', index=False)


rule include_average_abundance:
    input:
        abundance_tsv = 'resources/normalized_counts/RPM_{cutoff}.tsv'
    output:
        abundance_tsv = 'resources/normalized_counts/RPM_{cutoff}_mean.tsv',
    run:
        import pandas as pd
        abundance_df = pd.read_csv(input.abundance_tsv, sep='\t',  index_col = 'ENSG')
        abundance_df['mean'] = abundance_df.loc[:, SAMPLES].mean(axis=1)
        abundance_df.to_csv(str(output.abundance_tsv ), sep='\t')
