def _highly_variable_pearson_residuals(  # noqa: PLR0912, PLR0915
    adata: AnnData,
    *,
    theta: float = 100,
    clip: float | None = None,
    n_top_genes: int = 1000,
    batch_key: str | None = None,
    chunksize: int = 1000,
    check_values: bool = True,
    layer: str | None = None,
    subset: bool = False,
    inplace: bool = True,
) -> pd.DataFrame | None:
    view_to_actual(adata)
    x = _get_obs_rep(adata, layer=layer)
    computed_on = layer if layer else "adata.X"

    # Check for raw counts
    if check_values and not check_nonnegative_integers(x):
        warnings.warn(
            "`flavor='pearson_residuals'` expects raw count data, but non-integers were found.",
            UserWarning,
            stacklevel=3,
        )
    # check theta
    if theta <= 0:
        # TODO: would "underdispersion" with negative theta make sense?
        # then only theta=0 were undefined..
        msg = "Pearson residuals require theta > 0"
        raise ValueError(msg)
    # prepare clipping

    if batch_key is None:
        batch_info = np.zeros(adata.shape[0], dtype=int)
    else:
        batch_info = adata.obs[batch_key].values
    n_batches = len(np.unique(batch_info))

    # Get pearson residuals for each batch separately
    residual_gene_vars = []
    for batch in _batch_order(np.unique(batch_info)):
        adata_subset_prefilter = adata[batch_info == batch]
        x_batch_prefilter = _get_obs_rep(adata_subset_prefilter, layer=layer)

        # Filter out zero genes
        with settings.verbosity.override(Verbosity.error):
            nonzero_genes = np.ravel(x_batch_prefilter.sum(axis=0)) != 0
        adata_subset = adata_subset_prefilter[:, nonzero_genes]
        x_batch = _get_obs_rep(adata_subset, layer=layer)

        # Prepare clipping
        batch_clip = np.sqrt(x_batch.shape[0]) if clip is None else clip
        if batch_clip < 0:
            msg = "Pearson residuals require `clip>=0` or `clip=None`."
            raise ValueError(msg)

        if isinstance(x_batch, CSBase):
            x_batch = x_batch.tocsc()
            x_batch.eliminate_zeros()
            calculate_res = partial(
                _calculate_res_sparse,
                x_batch.indptr,
                x_batch.indices,
                x_batch.data.astype(np.float64),
            )
        else:
            x_batch = np.array(x_batch, dtype=np.float64, order="F")
            calculate_res = partial(_calculate_res_dense, x_batch)

        sums_genes = np.array(x_batch.sum(axis=0)).ravel()
        sums_cells = np.array(x_batch.sum(axis=1)).ravel()
        sum_total = np.sum(sums_genes)

        residual_gene_var = calculate_res(
            sums_genes=sums_genes,
            sums_cells=sums_cells,
            sum_total=np.float64(sum_total),
            clip=np.float64(batch_clip),
            theta=np.float64(theta),
            n_genes=x_batch.shape[1],
            n_cells=x_batch.shape[0],
        )

        # Add 0 values for genes that were filtered out
        unmasked_residual_gene_var = np.zeros(len(nonzero_genes))
        unmasked_residual_gene_var[nonzero_genes] = residual_gene_var
        _observe(batch, batch_clip, unmasked_residual_gene_var)
        residual_gene_vars.append(unmasked_residual_gene_var.reshape(1, -1))

    residual_gene_vars = np.concatenate(residual_gene_vars, axis=0)

    # Get rank per gene within each batch
    # argsort twice gives ranks, small rank means most variable
    ranks_residual_var = np.argsort(np.argsort(-residual_gene_vars, axis=1), axis=1)
    ranks_residual_var = ranks_residual_var.astype(np.float32)
    # count in how many batches a genes was among the n_top_genes
    highly_variable_nbatches = np.sum(
        (ranks_residual_var < n_top_genes).astype(int), axis=0
    )
    # set non-top genes within each batch to nan
    ranks_residual_var[ranks_residual_var >= n_top_genes] = np.nan
    ranks_masked_array = np.ma.masked_invalid(ranks_residual_var)
    # Median rank across batches, ignoring batches in which gene was not selected
    medianrank_residual_var = np.ma.median(ranks_masked_array, axis=0).filled(np.nan)

    means, variances = materialize_as_ndarray(_get_mean_var(x))
    df = pd.DataFrame.from_dict(
        dict(
            means=means,
            variances=variances,
            residual_variances=np.mean(residual_gene_vars, axis=0),
            highly_variable_rank=medianrank_residual_var,
            highly_variable_nbatches=highly_variable_nbatches.astype(np.int64),
            highly_variable_intersection=highly_variable_nbatches == n_batches,
        )
    )
    df = df.set_index(adata.var_names)

    # Sort genes by how often they selected as hvg within each batch and
    # break ties with median rank of residual variance across batches
    df.sort_values(
        ["highly_variable_nbatches", "highly_variable_rank"],
        ascending=[False, True],
        na_position="last",
        inplace=True,
    )

    high_var = np.zeros(df.shape[0], dtype=bool)
    high_var[:n_top_genes] = True
    df["highly_variable"] = high_var
    df = df.loc[adata.var_names, :]

    if inplace:
        adata.uns["hvg"] = {"flavor": "pearson_residuals", "computed_on": computed_on}
        logg.hint(
            "added\n"
            "    'highly_variable', boolean vector (adata.var)\n"
            "    'highly_variable_rank', float vector (adata.var)\n"
            "    'highly_variable_nbatches', int vector (adata.var)\n"
            "    'highly_variable_intersection', boolean vector (adata.var)\n"
            "    'means', float vector (adata.var)\n"
            "    'variances', float vector (adata.var)\n"
            "    'residual_variances', float vector (adata.var)"
        )
        adata.var["means"] = df["means"].values
        adata.var["variances"] = df["variances"].values
        adata.var["residual_variances"] = df["residual_variances"]
        adata.var["highly_variable_rank"] = df["highly_variable_rank"].values
        if batch_key is not None:
            adata.var["highly_variable_nbatches"] = df[
                "highly_variable_nbatches"
            ].values
            adata.var["highly_variable_intersection"] = df[
                "highly_variable_intersection"
            ].values
        adata.var["highly_variable"] = df["highly_variable"].values

        if subset:
            adata._inplace_subset_var(df["highly_variable"].values)

    else:
        if batch_key is None:
            df = df.drop(
                ["highly_variable_nbatches", "highly_variable_intersection"], axis=1
            )
        if subset:
            df = df.iloc[df.highly_variable.values, :]

        return df
