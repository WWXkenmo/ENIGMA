###################################S_mode_correction##################################################
require("sva")
require("purrr")

remove_batch_effect <- function(bulk_eset, ref_eset, varname_main_ct, varnames_sub_ct=NULL, n_pseudo_bulk=1000) {
    cat(date(), "generating pseudo bulk... \n")
    pseudo_bulk_main_ct = generate_pseudo_bulk_from_scRNA( ref_eset, varname_main_ct, n = n_pseudo_bulk )
    colnames(pseudo_bulk_main_ct$mat) = paste(varname_main_ct, seq_len(ncol(pseudo_bulk_main_ct$mat)), sep = "_")
    
    if(is.null(varnames_sub_ct)) pseudo_bulk_subs = NULL
    else {
        pseudo_bulk_subs = lapply( varnames_sub_ct, function(x) {
            sample_ids = subset( pData(ref_eset), eval(parse(text = paste0("!is.na(", x, ")") )) ) %>% rownames
            tmp = generate_pseudo_bulk_from_scRNA( ref_eset[,sample_ids], x, n = n_pseudo_bulk )
            colnames(tmp$mat) = paste(x, seq_len(ncol(tmp$mat)), sep = "_")
            return(tmp)
        } )
    }
    
    cat(date(), "do ComBat... \n")
    M_list = c(
        list( exprs(bulk_eset), pseudo_bulk_main_ct$mat ),
        lapply( pseudo_bulk_subs, function(x) x$mat )
    )
    gene_id = M_list %>% lapply(rownames) %>% Reduce(intersect, .)
    M_mix = M_list %>% lapply(function(x) x[gene_id,] ) %>% do.call(cbind, .)
    M_mix_log2 = log2(M_mix+1)
    
    M_correct = ComBat(
        M_mix_log2[rowSums(M_mix_log2)>0,],
        rep( c("bulk", "pseudo_bulk"), c( ncol(bulk_eset), ncol(M_mix_log2)-ncol(bulk_eset) ) )
    )
    M_correct = 2^M_correct - 1
    
    cat(date(), "restore reference... \n")
    res = c(
        list(
            M_correct[,colnames(bulk_eset)],
            restore_ref( pseudo_bulk_main_ct$frac, M_correct[,colnames(pseudo_bulk_main_ct$mat)] )
        ),
        lapply( pseudo_bulk_subs, function(x) {
            restore_ref( x$frac, M_correct[,colnames(x$mat)] )
        } )
    )
    names(res) = c( "bulk", varname_main_ct, varnames_sub_ct )
    
    return(res)
    cat(date(), "done. \n")
}


generate_pseudo_bulk_from_scRNA <- function(ref_eset, ct_varname, n=1000) {
    frac_init = table( pData(ref_eset)[,ct_varname] )/nrow(pData(ref_eset))
    
    frac = map2_df(frac_init, frac_init*2, rnorm, n=n) %>% sapply(function(x) x)
    
    # set nagative equal to zero
    frac[frac<0] <- 0
    frac <- frac[rowSums(frac)>0,]
    # normalization
    frac <- t(t(frac) %*% diag(1/rowSums(frac)))
    
    M_star <- NULL
    for(i in 1:nrow(frac)){
        Ma = lapply( pData(ref_eset)[,ct_varname] %>% unique, function(ct) {
            sample_ids = subset( pData(ref_eset), eval(parse( text = paste0(ct_varname, "==\"", ct, "\"") )) ) %>%
                rownames %>%
                sample( 100*frac[i,ct], replace = TRUE )
            exprs(ref_eset)[,colnames(ref_eset) %in% sample_ids]
        } ) %>%
            do.call(cbind, .)
        
        M_star <- cbind(M_star,rowSums(Ma))
    }
    
    return(list( frac=frac, mat=M_star ))
}


restore_ref <- function(frac, pseudo_bulk) {
    ref_exp <- NULL
    for(i in 1:nrow(pseudo_bulk)){
        coef <- nnls(frac, pseudo_bulk[i,])
        
        # pval = sapply( names(coef$coefficients), function(x) f.robftest(coef, x)$p.value )
        coef.v = coef$x
        # coef.v[which(pval>0.05)] = 0
        
        ref_exp <- rbind(ref_exp,coef.v)
    }
    rownames(ref_exp) <- rownames(pseudo_bulk)
    ref_exp[ref_exp<0] = 0
    
    return(ref_exp)
}


###########################B_mode_correction#######################################
B_mode_batch_effect_remove <- function(X,ref,fra){
    cat("Run B-mode to correct batch effect...")
    X_hat <- ref %*% t(fra)
    X_log <- log2(X+1)
    X_hat_log <- log2(X_hat+1)
    
    #
    cat("\n do Combat...")
    correct <- sva::ComBat(cbind(X_log,X_hat_log),batch=c(rep("1",ncol(X)),rep("2",nrow(fra))))
    X_correct <- 2^correct[,1:ncol(X)]-1
    
    correct <- 2^correct[,(ncol(X)+1):ncol(correct)]-1
    
    #
    #cat("/nRestore the gene specific expression...")
    ##Using nnls
    #gene_ref_exp <- NULL
    #for(i in 1:nrow(X_correct)){
    #    gene_exp <- rlm(fra, correct[i,])$coefficients
    #    gene_ref_exp <- rbind(gene_ref_exp,gene_exp)
    #}
    
    #rownames(gene_ref_exp) <- rownames(ref)
    #colnames(gene_ref_exp) <- colnames(ref)
    #gene_ref_exp[gene_ref_exp<0] <-0
    
    cat("\n Done")
    #res <- list()
    res <- X_correct
    res
}

















