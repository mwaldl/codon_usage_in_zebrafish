

rule average_abundance_per_method:
    input:
        abundance_tsv = 'resources/normalized_counts/RPM_'+str(config['RPM_cutoff'])+'_mean.tsv',
        methods_tsvs = expand('resources/codon_usage_per_gene/{method}_isoform_cu.tsv', method = ['shortest', 'longest', 'triple-only_longest', 'start-only_longest', 'priority_principal_triple_longest', 'triple-only_principal_longest','triple-only_principal-only_longest', 'priority_principal_start_triple_longest', 'triple-only_principal_longest'])
    output:
        cu_tsv = 'resources/methods_comparison/method_cu_count.tsv',
        pdf = 'results/method_comparison/method_cu_count.pdf',
        f_cu_tsv = 'resources/methods_comparison/method_cu_fraction.tsv',
        f_pdf = 'results/method_comparison/method_cu_fraction.pdf',
    run:
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sn

        abundance_df = pd.read_csv(input.abundance_tsv, sep='\t',  index_col = 'ENSG')
        abundance_df.drop(columns = SAMPLES, inplace=True)
        abundance_df.drop(columns= ['GeneType', 'GeneSymbol'], inplace = True)

        data = []
        for tsv in input.methods_tsvs:
            print(tsv)
            df =  pd.read_csv(tsv, sep='\t',  index_col = 'gene')
            for c in  ['isoform', 'length is multiple of 3', 'nt length','appris' ,'appris principal', 'Transcript type', 'Not found tag', 'Protein length', 'APPRIS Annotation', 'start found']:
                if c in df.columns:
                    df.drop(columns=[c], inplace = True)
            df=df.join(abundance_df, on=None, how='inner')
            print(df.columns)
            for c in df.columns:
                if c not in ['mean']:
                    df[c] = df[c]*df['mean']
            df.drop(columns = 'mean', inplace = True)
            df=df.fillna(0)
            df['method'] = tsv.split('/')[-1].split('.')[0]
            df.reset_index(inplace = True)
            data.append(df)
        df = pd.concat(data)
        df = df.groupby('method').sum()
        df = df.transpose()
        df.to_csv(str(output.cu_tsv ), sep='\t')

        fig, axs = plt.subplots(figsize=(7, 25.5), nrows=1, ncols=1)
        sn.heatmap(df, ax=axs, square = True, linewidths= 0.0, cmap = 'Blues')
        fig.savefig(output.pdf, bbox_inches="tight")


        for c in df.columns:
            column_sum = df[c].sum()
            df[c] = df[c] / column_sum
        df.to_csv(str(output.f_cu_tsv ), sep='\t')

        fig, axs = plt.subplots(figsize=(7, 25.5), nrows=1, ncols=1)
        sn.heatmap(df, ax=axs, square = True, linewidths= 0.0, cmap = 'Blues')
        fig.savefig(output.f_pdf, bbox_inches="tight")
