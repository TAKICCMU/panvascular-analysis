# Small, sourceable helpers used by the active nsNMF adapter.
# The inherited nsNMF fitting, rank-selection and many-to-one matching are retained.

descriptive_profile_spearman <- function(m, transform_name, donor_keys=NULL) {
  if (!is.matrix(m) || is.null(rownames(m)) || anyDuplicated(rownames(m)))
    stop("Named unique node matrix required")
  if (any(!is.finite(m))) stop("Nonfinite profile input")
  if (nrow(m)<2L) stop("At least two nodes required")
  if (!is.null(donor_keys) && (length(donor_keys)!=ncol(m) || anyNA(donor_keys)))
    stop("Resampling identity blocks must align with profiles")
  n_blocks <- if (is.null(donor_keys)) NA_integer_ else length(unique(donor_keys))
  pairs <- combn(seq_len(nrow(m)),2L,simplify=FALSE)
  do.call(rbind,lapply(pairs,function(idx) data.frame(
    node1=rownames(m)[idx[1L]],node2=rownames(m)[idx[2L]],
    rho=suppressWarnings(cor(m[idx[1L],],m[idx[2L],],method="spearman")),
    n_profiles=ncol(m),n_resampling_identity_blocks=n_blocks,
    block_unit="known_donor_or_unresolved_biological_sample_not_assumed_person",
    transform=transform_name,
    inference="DESCRIPTIVE_NO_P_OR_FDR",
    validation_family="same_profile_rank_association_not_independent_validation",
    stringsAsFactors=FALSE)))
}

# Rectangular Hungarian shortest-augmenting-path algorithm, O(n^2*m).
# Computes maximum total similarity for one distinct candidate per reference.
# This is a labelled sensitivity and does not replace inherited many-to-one best matches.
maximum_one_to_one <- function(similarity) {
  s <- as.matrix(similarity);n <- nrow(s);m <- ncol(s)
  if (!n || m<n || any(!is.finite(s)))
    stop("One-to-one sensitivity requires finite matrix and candidates >= references")
  cost <- max(s)-s
  u <- numeric(n+1L);v <- numeric(m+1L);p <- integer(m+1L);way <- integer(m+1L)
  for (i in seq_len(n)) {
    p[1L] <- i;j0 <- 1L;minv <- rep(Inf,m+1L);used <- rep(FALSE,m+1L)
    repeat {
      used[j0] <- TRUE;i0 <- p[j0];delta <- Inf;j1 <- 0L
      for (j in seq.int(2L,m+1L)) if (!used[j]) {
        cur <- cost[i0,j-1L]-u[i0+1L]-v[j]
        if (cur<minv[j]) {minv[j]<-cur;way[j]<-j0}
        if (minv[j]<delta) {delta<-minv[j];j1<-j}
      }
      for (j in seq_len(m+1L)) {
        if (used[j]) {u[p[j]+1L]<-u[p[j]+1L]+delta;v[j]<-v[j]-delta}
        else minv[j] <- minv[j]-delta
      }
      j0<-j1
      if (p[j0]==0L) break
    }
    repeat {
      j1<-way[j0];p[j0]<-p[j1];j0<-j1
      if (j0==1L) break
    }
  }
  assignment<-integer(n)
  for (j in seq.int(2L,m+1L)) if (p[j]>0L) assignment[p[j]]<-j-1L
  if (any(assignment==0L) || anyDuplicated(assignment)) stop("Invalid one-to-one solution")
  assignment
}

cosine_similarity_review <- function(A,B) {
  A<-sweep(A,2L,sqrt(colSums(A^2)),"/")
  B<-sweep(B,2L,sqrt(colSums(B^2)),"/")
  A[!is.finite(A)]<-0;B[!is.finite(B)]<-0
  t(A)%*%B
}

compare_basis_review <- function(reference,candidate,run_type,omitted=NA_character_,top_nodes=10L) {
  if (is.null(rownames(reference)) || is.null(rownames(candidate)) ||
      anyDuplicated(rownames(reference)) || anyDuplicated(rownames(candidate)) ||
      !setequal(rownames(reference),rownames(candidate)))
    stop("Exact matching state universe required for basis comparison")
  candidate<-candidate[rownames(reference),,drop=FALSE]
  if (any(!is.finite(reference)) || any(!is.finite(candidate))) stop("Nonfinite basis")
  sim<-cosine_similarity_review(reference,candidate)
  best<-apply(sim,1L,which.max)
  unique_match<-maximum_one_to_one(sim)
  build<-function(matches,mode) {
    counts<-tabulate(matches,nbins=ncol(candidate))
    do.call(rbind,lapply(seq_len(ncol(reference)),function(k) {
      ref_top<-rownames(reference)[order(reference[,k],decreasing=TRUE)[seq_len(min(top_nodes,nrow(reference)))]]
      cand_top<-rownames(candidate)[order(candidate[,matches[k]],decreasing=TRUE)[seq_len(min(top_nodes,nrow(candidate)))]]
      data.frame(run_type=run_type,omitted_or_bootstrap=omitted,
        matching_mode=mode,reference_module=colnames(reference)[k],matched_module=matches[k],
        cosine_similarity=sim[k,matches[k]],
        top_node_jaccard=length(intersect(ref_top,cand_top))/length(union(ref_top,cand_top)),
        n_reference_modules_matching_candidate=counts[matches[k]],
        duplicate_match=counts[matches[k]]>1L,
        stringsAsFactors=FALSE)
    }))
  }
  rbind(build(best,"inherited_many_to_one_primary"),
        build(unique_match,"one_to_one_hungarian_sensitivity"))
}
