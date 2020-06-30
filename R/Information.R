#' @title Correction Factor
#' @description If a numerical vector has a value of `0`, a correction factor is added to that value.
#' `0` `1` normalization is applied to the vector.
#' @param x a numerical vector
#' @param correction_factor correction factor that will be added to value of `0`. The default value is `0.00001`
#'
#' @return a numerical vector
#' @export
#'
#' @examples
#' vector.1 <- c(1,22,4,0)
#' correction(vector.1)
correction <- function(x, correction_factor = 0.00001){
  x[x == 0] = correction_factor
  x = x / sum(x)
}



#' @title Entropy of Variable
#' @description Calculate total entropy of each variable in a numerical vector. However, if a vector contains value of `0`.
#' \code{correction()} function is applied.
#'
#' @param x a numerical vector
#'
#' @return Total entropy of numerical vector.
#' @export
#' @references B. McMillan and D. Slepian, “Information Theory,” Proc. IRE, vol. 50, no. 5, pp. 1151–1157, 1962.
#' @examples
#' vector.1 <- c(1,4,5,0)
#' entropy(vector.1)
entropy <- function(x) {
  total_entropy <- data.frame()

  if(any(x == 0)){
    x = correction(x)
  }
  total_entropy = - sum(x * log(x))
}




#' @title Conditional Entropy of a Variable
#' @description Calculates conditional entropy of a target variable in Bayesian network given other variable(s).
#'
#' @param target variable whose conditional entropy to be calculated. Variable type is a character.
#' @param evidence.e evidence variable. Variable type is a character.
#' @param bn.model Bayesian network model including all variables
#'
#' @return Conditional entropy of target variable
#' @export
#' @references C. E. Shannon, “A Mathematical Theory of Communication,” Bell Syst. Tech. J., vol. 27, no. 3, pp. 379–423, Jul. 1948.
#'
#' @importFrom gRain querygrain
#' @importFrom gRain setEvidence
#'
#' @seealso \code{\link[gRain]{querygrain}} \code{\link[gRain]{setEvidence}}
#'
#' @examples
cond.ent <- function(target, evidence.e, bn.model) {
  target_t <- querygrain(bn.model, nodes = c(target), type = "marginal")
  evidence_e_t <- querygrain(bn.model, nodes = c(evidence.e), type = "marginal")
  t_states <- sapply(target_t,length)
  e_states <- length(evidence_e_t[[1]])
  e_state_names <- names(evidence_e_t[[1]])

  evidence_e_t <- as.data.frame(evidence_e_t)
  evidence_list <- list()
  evidence_prob <- list()
  evidence_entered_prob <- data.frame()
  cond_ent = rep(NA, length(target))
  names(cond_ent) = paste(names(target_t),"given",evidence.e,sep = "_")

  for (a in 1:e_states) {
    evidence_list[[a]] <- setEvidence(bn.model, nodes = c(evidence.e),
                                      states =c(e_state_names[a]))
    evidence_prob[[a]] <- querygrain(evidence_list[[a]],
                                     nodes = c(target),
                                     type = "marginal")
  }

  for(j in 1:length(target)) {
    evidence_entered_prob <- matrix(NA,nrow=t_states[j],ncol = e_states)
      for(i in 1:e_states) {
        # Correction for 0 probabilities
        if(0 %in% evidence_prob[[i]][[j]]){
        evidence_prob[[i]][[j]] = correction(evidence_prob[[i]][[j]])
        }
      evidence_prob[[i]][[j]]
      evidence_entered_prob[,i] = evidence_prob[[i]][[j]] * log(evidence_prob[[i]][[j]])
      }
    sum_evidence <- apply(evidence_entered_prob, 2, sum)
    cond_ent[j] = - sum(sum_evidence * evidence_e_t)
  }
  return(cond_ent)
}



#' @title Mutual Information of a Variable
#' @description Calculate mutual information value of variable(s) in Bayesian network given other variables
#'
#' @param target variable(s) whose mutual information to be calculated. Variable type is a character.
#' @param evidence evidence variable(s). Variable type is a character.
#' @param bn.model Bayesian network model including all variables
#'
#' @return a list containing probability distribution of target variable(s), mutual information between target(s) and evidence variable(s),
#'         and variable(s) which has maximum mutual information
#' @export
#' @references P. Baudot, M. Tapia, D. Bennequin, and J. M. Goaillard, “Topological information data analysis,” Entropy, vol. 21, no. 9, pp. 1–31, 2019.
#' @details Multiple target and evidence variables can be defined as vectors provided that the variable type is a character. Mutual information between
#'          target and evidence variables are calculated and returned.
#'
#' @importFrom gRain querygrain
#' @seealso \code{\link[gRain]{querygrain}} \code{\link{entropy}} \cr \code{\link{cond.ent}}
#'
#'
#' @examples
mut.inf <- function(target, evidence, bn.model) {
  num_ev = length(evidence)
  num_t = length(target)
  mut_i = matrix(NA, ncol=num_ev, nrow=num_t)
  dimnames(mut_i) = list(target, evidence)

  target_t <- querygrain(bn.model, nodes = c(target), type = "marginal")
  if(length(target) > 1) {
    eh <-   sapply(target_t, entropy)
  }else {
    target_t_t <- as.numeric(unlist(target_t))
    eh <- entropy(target_t_t)
   }

  for(i in 1:num_ev) {
    this_evidence = evidence[i]
    cond <- cond.ent(target, this_evidence, bn.model)
    temp_mut <- eh - cond
    # Match row names
    mut_i[, i] <- temp_mut[match(target,names(target_t))]
  }
  max_mut <- apply(mut_i,1,max)
  max_mut_ind <- apply(mut_i,1,which.max)
  max_mut_vars <- evidence[max_mut_ind]
  res_eh <- eh[match(target,names(target_t))]
  return(list(marg = target_t, ent = res_eh, max_vars = max_mut_vars, max_mut_vals = max_mut, max_mut_ind = max_mut_ind, mut_df = mut_i))
}


#' @title Algorithm for Posterior Distribution of variable
#' @description  The posterior distribution of the target variable(s) is calculated using the entropy and mutual information values.
#'
#' @param x a observation vector containing all variables in Bayesian network as characters. Each row in data frame defined as observation vector.
#' @param target variable(s) whose posterior distribution to be calculated. Variable type is a character.
#' @param query input variable(s) to be used when calculating the posterior distribution of target variables. Variable type is a character.
#' @param evidence a vector containing observed variables in the Bayesian network. Default value is `NULL`.
#' @param bn.model Bayesian network model including all variables.
#' @param num.iter an integer that shows number of query variables which will be used during calculation of
#'                 posterior distribution of target variable. Default value is `1`.
#'
#' @return a list containing posterior distribution of target variable.
#' @export
#'
#' @importFrom gRain querygrain
#' @importFrom gRain setEvidence
#'
#' @seealso \code{\link{mut.inf}} \code{\link[gRain]{setEvidence}}
#'
#' @examples
information.algorithm <- function(x, target, query, evidence = NULL, bn.model, num.iter = 1) {
  entered_ev <- c()

  # Enter Evidence
  if(!is.null(evidence)) {
    bn_model <- setEvidence(bn.model, evidence = as.list(evidence))
  }

  for(i in 1:num.iter){
    #Get nodes with maximum mutual information for each target
    mi <- mut.inf(target, query, bn.model)
    entered_ev <- c(entered_ev,unique(mi$max_vars))
    query <- query[!(query %in% mi$max_vars)]

    #Enter evidence for those nodes and remove them from query nodes
    query_ev <- x[entered_ev]
    bn_model <- setEvidence(bn.model, evidence = as.list(query_ev))
  }
  querygrain(bn_model, nodes = c(target), type = "marginal")
}




