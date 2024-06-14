import pandas as pd

def filter_isoforms(input_file, output_file, conditions, sort_columns, sort_order):
    # Load codon usage data
    cu_df = pd.read_csv(input_file, sep='\t', index_col='isoform')

    # Apply filter conditions
    for condition in conditions:
        cu_df = cu_df.query(condition)

    # Sort values based on the specified columns
    cu_df.sort_values(by=sort_columns, ascending=sort_order, inplace=True, kind='stable', na_position='last', ignore_index=False, key=None)

    # Drop duplicates to keep only one isoform per gene
    cu_df.drop_duplicates(subset='gene', keep='first', inplace=True)

    # Save the filtered dataframe to the output file
    cu_df.to_csv(output_file, sep='\t')

rule cu_shortest_isoform:
    input:
        tsv='resources/codon_usage_per_isoform/appris_cu.tsv'
    output:
        tsv='resources/codon_usage_per_gene/shortest_isoform_cu.tsv'
    run:
        filter_isoforms(input.tsv, output.tsv, [], ['nt length'], [True])

rule cu_longest_isoform:
    input:
        tsv='resources/codon_usage_per_isoform/appris_cu.tsv'
    output:
        tsv='resources/codon_usage_per_gene/longest_isoform_cu.tsv'
    run:
        filter_isoforms(input.tsv, output.tsv, [], ['nt length'], [False])

rule cu_start_only_longest_isoform:
    input:
        tsv='resources/codon_usage_per_isoform/appris_cu.tsv'
    output:
        tsv='resources/codon_usage_per_gene/start-only_longest_isoform_cu.tsv'
    run:
        filter_isoforms(input.tsv, output.tsv, ['`start found`'], ['nt length'], [False])

rule cu_triple_only_longest_isoform:
    input:
        tsv='resources/codon_usage_per_isoform/appris_cu.tsv'
    output:
        tsv='resources/codon_usage_per_gene/triple-only_longest_isoform_cu.tsv'
    run:
        filter_isoforms(input.tsv, output.tsv, ['`length is multiple of 3`'], ['nt length'], [False])

rule cu_priority_principal_triple_length_isoform:
    input:
        tsv='resources/codon_usage_per_isoform/appris_cu.tsv'
    output:
        tsv='resources/codon_usage_per_gene/priority_principal_triple_longest_isoform_cu.tsv'
    run:
        filter_isoforms(input.tsv, output.tsv, [],
                        ['appris principal', 'length is multiple of 3', 'nt length'], [False, False, False])

rule cu_priority_principal_start_triple_length_isoform:
    input:
        tsv='resources/codon_usage_per_isoform/appris_cu.tsv'
    output:
        tsv='resources/codon_usage_per_gene/priority_principal_start_triple_longest_isoform_cu.tsv'
    run:
        filter_isoforms(input.tsv, output.tsv, [],
                        ['appris principal','start found' ,'length is multiple of 3', 'nt length'], [False, False, False, False])

rule cu_triple_only_principal_length_isoform:
    input:
        tsv='resources/codon_usage_per_isoform/appris_cu.tsv'
    output:
        tsv='resources/codon_usage_per_gene/triple-only_principal_longest_isoform_cu.tsv'
    run:
        filter_isoforms(input.tsv, output.tsv, ['`length is multiple of 3`'],
                        ['appris principal', 'nt length'], [False, False])

rule cu_triple_only_principal_only_length_isoform:
    input:
        tsv='resources/codon_usage_per_isoform/appris_cu.tsv'
    output:
        tsv='resources/codon_usage_per_gene/triple-only_principal-only_longest_isoform_cu.tsv'
    run:
        filter_isoforms(input.tsv, output.tsv, ['`length is multiple of 3`', '`appris principal`'],
                        ['nt length'], [False])
