#' Title
#' @title Learn Bayesian Network models
#' @description Learn the structure and parameters of Bayesian network either using algorithms
#'in R or using user-defined algorithms. If the structure of Bayesian network defined in advance,
#'parameters of it learned and model added to list.
#'
#' @param Train_Data a data frame containing discrete variables and factors used for learning Bayesian model
#' @param Algorithm.Names a list of the algorithms used in Bayesian network model learning
#' @param Structure.List If there is pre-defined structure of Bayesian network, parameters of it are learned
#'
#'
#'
#' @return List of models
#' @export
#'
#'
#'
#'
#'
Model_Learning <- function(Train_Data, Algorithm.Names = NULL, Structure.List = NULL) {

  Models <- list()
  Train_Data <- na.omit(Train_Data)
  Train_Data <- apply(Train_Data, 2, as.factor)
  Train_Data <- as.data.frame(Train_Data, stringsAsFactors = TRUE)

  if(!is.null(Algorithm_Names)) {
    Algorithms <- Algorithm_Names

    for (a in 1:length(Algorithms)) {

      model <- Algorithms[[a]](Train_Data)
      modelB <- bn.fit(model, Train_Data, method = "bayes")
      modelgr <- as.grain(modelB)
      Models[[a]] <- modelgr
      names(Models)[a] <- model$learning$algo
      rm(model, modelB, modelgr)

    }

  }

  Model_Expert <- list()
  if(!is.null(Arc_Set)) {

    if(names(Arc_Set[[1]])[1] == "learning"){

      for (b in 1:length(Arc_Set)) {

        Model_GAS <- bn.fit(Arc_Set[[b]], Train_Data, method = "bayes")
        Model_GAS_Gr  <- as.grain(Model_GAS)
        Model_Expert[[b]] <- Model_GAS_Gr

        if(is.null(Algorithm_Names)) {

          Models[[b]] <- Model_GAS_Gr
          names(Models)[b] <- paste("Expert_Model", b, sep = "_")

        }else {

          Models[[length(Algorithm_Names)+b]] <- Model_Expert[[b]]
          names(Models)[length(Algorithm_Names)+b] <-  paste("Expert_Model", b, sep = "_")

        }

      }

    }else {

      Model_GAS <- bn.fit(Arc_Set, Train_Data, method = "bayes")
      Model_GAS_Gr  <- as.grain(Model_GAS)
      if(is.null(Algorithm_Names)){
        Models[[1]] <- Model_GAS_Gr
        names(Models)[1] <- "Expert_Model"

      }else {

        Models[[length(Algorithm_Names)+1]] <- Model_GAS_Gr
        names(Models)[length(Algorithm_Names)+1] <- "Expert_Model"

      }

    }

  }

  return(Models)

}
