# Common functions for PRSet analyses

standardize <- function(x, n) {
  (x - mean(x)) / (sqrt(sum(x^2) / n))
}

scale_this <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

compute_nagelkerkeR2 <- function(rr){
  n <- nrow(rr$model)
  R2 <- (1 - exp((rr$dev - rr$null)/n))/(1 - exp(-rr$null/n))
  return(R2)
}

convert_pathway_id <- function(pathway, term) {
  #' available term options = go, mgi, kegg, pid, reactome
  if(!term %in% c("go", "mgi", "kegg", "pid", "reactome", "biocarta")) {
    message("Did you specify one of the following accepted terms: go, mgi, kegg, pid, reactome, biocarta")
    stop()
  }
  
  switch(term,
         go = Term(pathway),
         kegg = gsub("^KEGG_$", "", pathway) %>% tolower(),
         pid = gsub("^PID_$", "", pathway) %>% tolower(),
         reactome = gsub("^REACTOME_$", "", pathway) %>% tolower(),
         biocarta = gsub("^BIOCARTA_$", "", pathway) %>% tolower())
}

save_model_coef <- function(model) {
  lasso_selected_pathways <- model %>%
    coef(., s = "lambda.1se") %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    setNames(c("go_id", "coef")) %>% 
    filter(!is.na(coef)) %>%
    filter(coef != 0) %>%
    arrange(desc(coef))
  
  # write_csv(lasso_selected_pathways, paste0(getwd(), "/", filename))
  
  return(lasso_selected_pathways)
}

evaluate_single_classifier <- function(results_df) {
  
  # requires caret
  pred_funs <- list(sens = sensitivity,
                    spec = specificity,
                    ppv = posPredValue,
                    npv = negPredValue,
                    n_fn = function(pred, true) sum(true == 1 & pred == 0),
                    n_fp = function(pred, true) sum(true == 0 & pred == 1),
                    n = function(pred, true) length(true))
  out <- map(pred_funs, 
             ~.x(factor(as.character(results_df$predicted_label), levels = c("0", "1")), 
                 factor(as.character(results_df$outcome), levels = c("0", "1")))) %>% 
    bind_cols()
  
  out
}

plot_roc <- function(data=data,bin=0.01,roccol="green",sp=19,output="roc.pdf")
{
  pn <- nrow(subset(data,data[,2]== 1))
  fn <- nrow(data) - pn
  diag = data.frame(x = seq(0,1,by=bin),y = seq(0,1,by=bin))
  cutoffs <- seq(0,1,by=bin)
  x = 0
  y = 0
  for (i in cutoffs)
  {
    tpn <- nrow(subset(data, data[,1] >= i & data[,2] == 1))
    fpn <- nrow(subset(data, data[,1] >= i & data[,2] == 0))
    fnn <- pn - tpn
    tnn <- fn - fpn
    tpr <- tpn/(tpn + fnn)
    fpr <- fpn/(fpn + tnn)
    x <-c(x,fpr)
    y <- c(y,tpr)
  }
  FPR = ''
  TPR = ''
  rocdata <- data.frame("FPR" = x, "TPR" = y)
  p <- ggplot(data=rocdata, aes(x = FPR, y = TPR)) + geom_point(color=roccol) + geom_line(color=roccol) + geom_line(data=diag, aes(x=x,y=y),color="red") 
  f <- p + geom_point(data=diag, aes(x=x,y=y),color="red",shape=sp) + theme(axis.text = element_text(size=16), title=element_text(size=18)) + labs(x="False Positive Rate", y="True Positive Rate", title="ROC curve")
  ggsave(f, filename = output, width=8,height=6,units=c("in"))
}