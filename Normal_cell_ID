def get_cell_types(GEM):
#     Stat_sig_DEG_cibersort_blood_cells = pd.read_csv('../../../../ALL/Data/FerrandoLab/SingleCell/Stat_sig_DEG_Cibersort_Blood_cells.txt', header = 0, index_col=0, sep = '\t')
    Stat_sig_DEG_cibersort_blood_cells = pd.read_csv('../../../../ALL/Data/FerrandoLab/SingleCell/Stat_sig_DEG_Cibersort_Blood_cells_plusCD34.txt', header = 0, index_col=0, sep = '\t')
    diff_gene = {}
    nb_clusters = len(set(GEM.obs['leiden']))
    
    for i in range(nb_clusters):
        print(i)
        sc.tl.rank_genes_groups(GEM, groupby='leiden', groups = [i], reference= 'rest', n_genes = GEM.X.shape[0])
        diff_gene_exp_cluster = pd.DataFrame(data = list(zip(GEM.uns['rank_genes_groups']['names'], GEM.uns['rank_genes_groups']['logfoldchanges'], GEM.uns['rank_genes_groups']['scores'], GEM.uns['rank_genes_groups']['pvals_adj'])), 
                                             columns= ('names', 'logfoldchanges', 'scores', 'pvals_adj'))
        diff_gene_exp_cluster['names'] = diff_gene_exp_cluster['names'].str[0]
        diff_gene_exp_cluster['logfoldchanges'] = diff_gene_exp_cluster['logfoldchanges'].str[0]
        diff_gene_exp_cluster['scores'] = diff_gene_exp_cluster['scores'].str[0]
        diff_gene_exp_cluster['pvals_adj'] = diff_gene_exp_cluster['pvals_adj'].str[0]
        diff_gene[str(i)] = diff_gene_exp_cluster 
    
    df_clust_marker_genes = pd.DataFrame(data = 0, index = Stat_sig_DEG_cibersort_blood_cells.index, columns=['clust_' + str(i) for i in range(nb_clusters)])

    for clust_nb in range(nb_clusters):
        DE_genes = diff_gene[str(clust_nb)].loc[(diff_gene[str(clust_nb)]['pvals_adj'] < 0.01) & (diff_gene[str(clust_nb)]['scores'] > 0)]['names']
        for top_genes in DE_genes:
            if top_genes in Stat_sig_DEG_cibersort_blood_cells.index:
                df_clust_marker_genes.loc[top_genes][clust_nb] = 1


    result_kendall = pd.DataFrame(data = 0, columns = Stat_sig_DEG_cibersort_blood_cells.columns, index= df_clust_marker_genes.columns, dtype=float)
    for clust in result_kendall.index:
        for cell_type in Stat_sig_DEG_cibersort_blood_cells.columns:
#             res, pval = scipy.stats.kendalltau(x = df_clust_marker_genes[clust], y = Stat_sig_DEG_cibersort_blood_cells[cell_type])
            res, pval = scipy.stats.pearsonr(x = df_clust_marker_genes[clust], y = Stat_sig_DEG_cibersort_blood_cells[cell_type])
            #print(res, pval)

            result_kendall.loc[clust][cell_type] = res

    cell_id_per_cluster_kendall = pd.DataFrame(data=0, columns = ('cell_id', 'correlation'), index = result_kendall.index)
    cell_id_per_cluster_kendall['cell_id'] = result_kendall.idxmax(axis=1)
    cell_id_per_cluster_kendall['correlation'] = result_kendall.max(axis=1)
    return cell_id_per_cluster_kendall, result_kendall, df_clust_marker_genes, Stat_sig_DEG_cibersort_blood_cells
    
    
    def normalization_pipeline(GEM, sample = 'sample_name', preprocessing = True, pca = True, umap = True):
    if preprocessing == True:
        sc.pp.filter_cells(GEM, min_genes=200)
        sc.pp.filter_genes(GEM, min_cells=3)
        sc.pp.normalize_total(GEM, target_sum = 1.0e6, exclude_highly_expressed=True)
        sc.pp.log1p(GEM)
    if pca == True:
        sc.tl.pca(GEM)
        sc.pl.pca(GEM, title = 'PCA - ' + sample)
        if umap == True:
            sc.pp.neighbors(GEM, n_pcs = 20, n_neighbors = 15 )
            sc.tl.umap(GEM, n_components = 2, min_dist = 1)
            sc.pl.umap(GEM, title = 'UMAP - ' + sample)
    return GEM
    
    def subset_anndata(anndata_object, observation = 'leiden', value = ['0']):
    import anndata
    indexing = []
    counter = 0
    for i in value:
        counter += 1
        clustering = pd.Index(list(anndata_object.obs[observation]))
        if counter == 1:
            indexing = clustering.get_loc(i)
        if counter > 1:
            indexing_1plus = clustering.get_loc(i)
            indexing = indexing | indexing_1plus

    new_obj = anndata.AnnData(X = anndata_object.raw.X[indexing], 
                         obs = anndata_object.obs[indexing], 
                         var = anndata_object.raw.var)
    
    return new_obj
    
    def normalization_pipeline(GEM, sample = 'sample_name', preprocessing = True, pca = True, umap = True, batch_correction = False, batches = 'orig.ident'):
    if preprocessing == True:
        print('Preprocessing')
        sc.pp.filter_cells(GEM, min_genes=200)
        sc.pp.filter_genes(GEM, min_cells=3)
        if batch_correction == True:
            sc.pp.combat(GEM, batches)
        sc.pp.normalize_total(GEM, target_sum = 1.0e6, exclude_highly_expressed=True)
        sc.pp.log1p(GEM)
    if pca == True:
        print('PCA')
        sc.tl.pca(GEM)
        sc.pl.pca(GEM, title = 'PCA - ' + sample)
        if umap == True:
            print('UMAP')
            sc.pp.neighbors(GEM, n_pcs = 20, n_neighbors = 15 )
            sc.tl.umap(GEM, n_components = 2, min_dist = 1)
            sc.pl.umap(GEM, title = 'UMAP - ' + sample)
    return GEM
    
    def get_cell_types(GEM, observation = 'leiden', ref = 'cibersort', ref_for_DGE = 'rest'):
    
    if type(ref) is str:
        if ref == 'cibersort':    
            refmat = pd.read_csv('../../../../ALL/Data/FerrandoLab/SingleCell/Stat_sig_DEG_Cibersort_Blood_cells.txt', header = 0, index_col=0, sep = '\t')
        elif ref == 'cibersort+CD34':
            refmat = pd.read_csv('../../../../ALL/Data/FerrandoLab/SingleCell/Stat_sig_DEG_Cibersort_Blood_cells_plusCD34.txt', header = 0, index_col=0, sep = '\t')
        elif ref == 'RNAseq':
            refmat = pd.read_csv('C:/Users/Alexandre/Desktop/Columbia/CalifanoLab/Scripts/R/PBMC_deconvolution/Results/Diff_exp_binary_immune_cells.txt', header = 0, index_col=0, sep = '\t')
        elif ref == 'RNAseq+leukemia':
            refmat = pd.read_csv('C:/Users/Alexandre/Desktop/Columbia/CalifanoLab/Scripts/R/PBMC_deconvolution/Results/Diff_exp_binary_immune_cells_with_leukemia.txt', header = 0, index_col=0, sep = '\t')
        elif ref == 'single_cell_immune':
            refmat = pd.read_csv('C:/Users/Alexandre/Desktop/Columbia/CalifanoLab/Scripts/Python/Jupyter_Notebook/Single_cell_GEM/Results/DiffGeneExp_single_cell_ref.txt', header = 0, index_col=0, sep = '\t')
        elif ref == 'single_cell_immune_noCD34':
            refmat = pd.read_csv('C:/Users/Alexandre/Desktop/Columbia/CalifanoLab/Scripts/Python/Jupyter_Notebook/Single_cell_GEM/Results/DiffGeneExp_single_cell_ref_noCD34.txt', header = 0, index_col=0, sep = '\t')
        elif ref == 'single_cell_PBMC':
            refmat = PBMC_DGE
        elif ref == 'single_cell_immune_prolif_markers':
            refmat = pd.read_csv('./Results/DiffGeneExp_single_cell_ref_DNTT_MKI67_MYC_CCND1_TFRC.txt', header = 0, index_col=0, sep = '\t')
        elif ref == 'single_cell_immune_prolif_markers_highLeuk_lowNormal':
            refmat = pd.read_csv('./Results/DiffGeneExp_single_cell_ref_DNTT_MKI67_MYC_CCND1_TFRC_vs_lowExpNormal.txt', header = 0, index_col=0, sep = '\t')
        elif ref == 'single_cell_immune_prolif_markers_highLeuk_vsNormalNoT_noMono_noB':
            refmat = pd.read_csv('./Results/DiffGeneExp_single_cell_ref_DNTT_MKI67_MYC_CCND1_TFRC_vs_lowExpNormal_notT_NK_mono_B.txt', header = 0, index_col=0, sep = '\t')
        else:
            raise ValueError(str(ref) + ' is not a valid reference.')
    else:
        refmat = ref 

    
    diff_gene = {}
    nb_clusters = len(set(GEM.obs[observation]))
    cell_types = list(set(GEM.obs[observation]))
    
    for i in cell_types:
        if str(i) == ref_for_DGE:
            cell_types.remove(str(i))
    
    for i in cell_types:
        print(i)
        sc.tl.rank_genes_groups(GEM, groupby = observation, groups = [i], reference= ref_for_DGE, n_genes = GEM.X.shape[0])
        diff_gene_exp_cluster = pd.DataFrame(data = list(zip(GEM.uns['rank_genes_groups']['names'], GEM.uns['rank_genes_groups']['logfoldchanges'], GEM.uns['rank_genes_groups']['scores'], GEM.uns['rank_genes_groups']['pvals_adj'])), 
                                             columns= ('names', 'logfoldchanges', 'scores', 'pvals_adj'))
        diff_gene_exp_cluster['names'] = diff_gene_exp_cluster['names'].str[0]
        diff_gene_exp_cluster['logfoldchanges'] = diff_gene_exp_cluster['logfoldchanges'].str[0]
        diff_gene_exp_cluster['scores'] = diff_gene_exp_cluster['scores'].str[0]
        diff_gene_exp_cluster['pvals_adj'] = diff_gene_exp_cluster['pvals_adj'].str[0]
        diff_gene[str(i)] = diff_gene_exp_cluster 
    
    df_clust_marker_genes = pd.DataFrame(data = 0, index = refmat.index, columns=[str(i) for i in cell_types])

    for clust in cell_types:
        DE_genes = diff_gene[str(clust)].loc[(diff_gene[str(clust)]['pvals_adj'] < 0.01) & (diff_gene[str(clust)]['scores'] > 0)]['names']
        for top_genes in DE_genes:
            if top_genes in refmat.index:
                df_clust_marker_genes.loc[top_genes][clust] = 1


    result_kendall = pd.DataFrame(data = 0, columns = refmat.columns, index= df_clust_marker_genes.columns, dtype=float)
    for clust in result_kendall.index:
        for cell_type in refmat.columns:
            res, pval = scipy.stats.kendalltau(x = df_clust_marker_genes[clust], y = refmat[cell_type])
            #print(res, pval)

            result_kendall.loc[clust][cell_type] = res

    cell_id_per_cluster_kendall = pd.DataFrame(data=0, columns = ('cell_id', 'correlation'), index = result_kendall.index)
    cell_id_per_cluster_kendall['cell_id'] = result_kendall.idxmax(axis=1)
    cell_id_per_cluster_kendall['correlation'] = result_kendall.max(axis=1)
    
    return cell_id_per_cluster_kendall, result_kendall, df_clust_marker_genes, refmat
    
    def correlations_vs_number_of_merged_clusters(GEM, cell_types, observation = 'Cell_id', DGE_ref = 'single_cell_immune'):

    GEM.obs[observation] = 'Other'
    GEM.obs[observation] = GEM.obs[observation].astype('str')
    res = []
    nb_cells = []
    list_correlations = []
    list_indices = []
    counter = 0
    GEM.obs
    for i,j in zip(cell_types[1]['Leukemia'].sort_values(ascending = False), cell_types[1]['Leukemia'].sort_values(ascending = False).index):
        list_correlations.append(i)
        list_indices.append(j)
    for corr, ind in zip(list_correlations, list_indices):
        print(corr)
        counter += 1
        if counter < len(list_correlations):
            GEM.obs.loc[GEM.obs['leiden'] == ind, observation] = 'Leukemia'
            corr_leuk = get_cell_types(GEM, observation, DGE_ref)
            res.append(corr_leuk[0].loc['Leukemia']['correlation'])
            nb_cells.append(GEM.obs[observation].value_counts()['Leukemia'])

    return res, nb_cells
    
    def optimize_leiden_by_leukemic_cell_correlation(GEM, resolution_range=(0.5, 1.5), nb_iterations = 10):
    result = []
    res_range = []
    counter = 0
    
    if nb_iterations == None:
        resolution_range_to_test = resolution_range
    else:
        resolution_range_to_test = np.arange(start = resolution_range[0], stop = resolution_range[1], step = (resolution_range[1]-resolution_range[0])/nb_iterations)
    
    for res in resolution_range_to_test:
        
    #compute leiden clustering at resolution i, then assign cell types
        sc.tl.leiden(GEM, resolution = res)    
        Cell_types = get_cell_types(GEM, 'leiden')

    #Merge identical cell types
        GEM.obs['Cell_id'] = 0
        GEM.obs['Cell_id'] = GEM.obs['Cell_id'].astype('str')
        for i in Cell_types[0].index:
            GEM.obs.loc[GEM.obs['leiden'] == str(i),'Cell_id'] = Cell_types[0].loc[str(i)]['cell_id']

        GEM.obs['Cell_id'] = GEM.obs['Cell_id'].astype('category')
        merged_cell_types = get_cell_types(GEM, 'Cell_id')

        res_range.append(res)
        result.append(merged_cell_types[0].loc['Leukemia']['correlation'])
    
    return res_range, result                   
    
    def merge_cell_types(cell_types, GEM):
    GEM.obs['Cell_id'] = 0
    GEM.obs['Cell_id'] = GEM.obs['Cell_id'].astype('str')
    for i in cell_types.index:
        GEM.obs.loc[GEM.obs['leiden'] == i,'Cell_id'] = cell_types.loc[i]['cell_id']
    sc.pl.umap(GEM, color = 'Cell_id')
    GEM.obs['Cell_id'].value_counts()
    return GEM
    
    def merge_binary_cell_types(cell_types, GEM, threshold = 0.3, category = 'Leukemia'):
    GEM.obs['Cell_id'] = 0
    GEM.obs['Cell_id'] = GEM.obs['Cell_id'].astype('str')
    for i in cell_types[1].loc[category].index:
        if cell_types[1].loc[str(i)][category] >= threshold:
            GEM.obs.loc[GEM.obs['leiden'] == i,'Cell_id'] = category
        else:
            GEM.obs.loc[GEM.obs['leiden'] == i,'Cell_id'] = 'Other'
    
    sc.pl.umap(GEM, color = 'Cell_id')
    GEM.obs['Cell_id'].value_counts()
    return GEM
    
    def create_differential_expression_reference_mat(GEM, 
                                                 category_for_cells_in_marker_sel = 'Disease',
                                                 subset_for_marker_selection = 'B-ALL',
                                                 category_for_cells_in_known_cell_types = 'Cell_type',
                                                 subset_for_known_cell_types = ['T_cell_cytotoxic_naive', 'T_cell_CD4_helper', 'T_cell_regulatory', 'T_cell_memory', 'T_cell_cytotoxic', 'B_cells', 'CD34', 'CD56_NK', 'CD14_monocytes'],
                                                 SD_above_average = 1, 
                                                 markers = ['DNTT', 'MKI67', 'MYC', 'CCND1', 'TFRC'],
                                                reference_against_opposite = False,
                                                equilibrate_classes = True):

    
    anndata_indexes = GEM.obs.loc[GEM.obs[category_for_cells_in_marker_sel] == subset_for_marker_selection].index
    boolean_index = [i in subset_for_known_cell_types for i in GEM.obs[category_for_cells_in_known_cell_types]]
    boolean_index = GEM.obs.loc[boolean_index].index
    anndata_indexes = anndata_indexes | boolean_index
    anndata_indexes = GEM.obs.loc[anndata_indexes.tolist()].index
    anndata_indexes = anndata_indexes.tolist()
    
    GEM.obs['cell_selection'] = 'Other_cell'
    GEM.obs['cell_selection'] = GEM.obs['cell_selection'].astype('str')
    GEM.obs.loc[anndata_indexes,'cell_selection'] = 'subset'
    GEM = subset_anndata(GEM, 'cell_selection', ['subset'])

    GEM.obs['High_exp_markers'] = 'remove'
    GEM.obs['High_exp_markers'] =  GEM.obs['High_exp_markers'].astype('str')

    for marker in markers:
        index_marker = GEM.var_names.get_loc(marker)
        marker_expression = GEM.X.T[index_marker].todense()
        indices_nonzero_expr_markers = [i > 0 for i in marker_expression]
        marker_expression_nonzero = marker_expression[indices_nonzero_expr_markers]
        marker_expression_nonzero = marker_expression_nonzero.tolist()[0]
        mean_marker_expr = scipy.stats.describe(marker_expression_nonzero)[2]
        sd_marker_expr = scipy.stats.describe(marker_expression_nonzero)[3]
        marker_expression = marker_expression.tolist()[0]
        
        high_marker_expr = marker_expression > (mean_marker_expr + SD_above_average*sd_marker_expr)
        GEM.obs.loc[high_marker_expr,'High_exp_markers'] = 'Leukemia'
    
    GEM.obs['random'] = 0
    GEM.obs['random'] = GEM.obs['random'].astype('int')
    
    for cell_types in subset_for_known_cell_types:
        GEM.obs.loc[GEM.obs[category_for_cells_in_known_cell_types] == cell_types, 'High_exp_markers'] = cell_types

        #subsetting the majority classes
        popsize = sum(GEM.obs[category_for_cells_in_known_cell_types] == cell_types)
        random_sampling = random.sample(range(popsize), popsize)
        GEM.obs.loc[GEM.obs[category_for_cells_in_known_cell_types] == cell_types, 'random'] = random_sampling

    #equal sizes for each class (threshold is the number of leukemic cells in the differential gene expression calculation)
    threshold_subsampling = sum(GEM.obs['High_exp_markers'] == 'Leukemia')
    GEM.obs['subsampling_thresold'] = 'above'
    GEM.obs['subsampling_thresold'] = GEM.obs['subsampling_thresold'].astype('str')
    GEM.obs.loc[(GEM.obs['random'] < threshold_subsampling),'subsampling_thresold'] = 'below'
    
    GEM.raw = GEM
    GEM = subset_anndata(GEM, 'High_exp_markers', subset_for_known_cell_types + ['Leukemia'])
    GEM.raw = GEM
    if equilibrate_classes == True:
        GEM = subset_anndata(GEM, 'subsampling_thresold', ['below'])

    GEM = normalization_pipeline(GEM, pca = False, umap = False)
    
    print(GEM.obs['High_exp_markers'].value_counts())
    
    #create differential expression profile 
    dict_DGE_ref = {}
    cell_types = GEM.obs['High_exp_markers'].value_counts().index.tolist()

    for i in cell_types:
        if reference_against_opposite == True:
            if i == 'Leukemia':
                ref_group = 'rest'
            else:
                ref_group = 'Leukemia'
        else:
            ref_group = 'rest'
            
        print(ref_group)        
        DGE = sc.tl.rank_genes_groups(GEM, groupby='High_exp_markers', groups = [i], reference= ref_group, n_genes = GEM.X.shape[0])
        dict_DGE_ref[i] = pd.DataFrame(data = list(zip(GEM.uns['rank_genes_groups']['names'], 
                                                         GEM.uns['rank_genes_groups']['logfoldchanges'], 
                                                         GEM.uns['rank_genes_groups']['scores'], 
                                                         GEM.uns['rank_genes_groups']['pvals_adj'])), 
                                             columns= ('names', 'logfoldchanges', 'scores', 'pvals_adj'))
        dict_DGE_ref[i]['names'] = dict_DGE_ref[i]['names'].str[0]
        dict_DGE_ref[i]['logfoldchanges'] = dict_DGE_ref[i]['logfoldchanges'].str[0]
        dict_DGE_ref[i]['scores'] = dict_DGE_ref[i]['scores'].str[0]
        dict_DGE_ref[i]['pvals_adj'] = dict_DGE_ref[i]['pvals_adj'].str[0]

    DGE_ref = pd.DataFrame(data = 0, index = GEM.var_names.tolist(), columns = [i for i in dict_DGE_ref.keys()])
    for i in dict_DGE_ref.keys():
        threshold = 5
        for j in dict_DGE_ref[i]['names'][dict_DGE_ref[i]['scores'] > threshold]:
            DGE_ref.loc[j][i] = 1

    return DGE_ref
