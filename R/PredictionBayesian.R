globalVariables(c("Algorithm.Names", "Structure.List"))
#' @title Learn Bayesian Network models
#' @description Learn the structure and parameters of Bayesian network either using algorithms
#' in R or using user-defined algorithms. If the structure of Bayesian network defined in advance,
#' parameters of it learned and model added to list.
#'
#' @param train.data a data frame containing discrete variables and factors used for learning Bayesian model
#' @param algorithm.names a vector containing the algorithms used in Bayesian network model learning. The default value is `NULL`
#' @param structure.list If there is pre-defined structure of Bayesian network, parameters of it are learned. The default value is `NULL`
#'
#' @details Available algorithms in R could be used to learn structure and parameters. See \code{\link[bnlearn]{structure-learning}}
#' for available algorithms. Also user-defined algorithms could be defined and used. If the structure(s) of Bayesian network defined in advance
#' parameters of it learned and model(s) added to list.
#'
#' @return List of Bayesian network models. Each element of the list is represented by the name of the algorithm.
#' If the list element is the predefined model, it is named \code{Predefined_Model_1}, \code{Predefined_Model_2}...
#'
#' @export
#'
#' @importFrom stats na.omit
#' @importFrom bnlearn bn.fit
#' @importFrom bnlearn as.grain
#'
#' @seealso \code{\link[bnlearn]{bn.fit}} \code{\link[bnlearn]{as.grain}}
#'
Model_Learning <- function(train.data, algorithm.names = NULL, structure.list = NULL) {

  Models <- list()
  Train_Data <- na.omit(train.data)
  Train_Data <- apply(Train_Data, 2, as.factor)
  Train_Data <- as.data.frame(Train_Data, stringsAsFactors = TRUE)

  if(!is.null(algorithm.names)) {
    Algorithms <- algorithm.names

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
  if(!is.null(structure.list)) {

    if(names(structure.list[[1]])[1] == "learning"){

      for (b in 1:length(structure.list)) {

        Model_GAS <- bn.fit(structure.list[[b]], Train_Data, method = "bayes")
        Model_GAS_Gr  <- as.grain(Model_GAS)
        Model_Expert[[b]] <- Model_GAS_Gr

        if(is.null(algorithm.names)) {

          Models[[b]] <- Model_GAS_Gr
          names(Models)[b] <- paste("Predefined_Model", b, sep = "_")

        }else {

          Models[[length(algorithm.names)+b]] <- Model_Expert[[b]]
          names(Models)[length(algorithm.names)+b] <-  paste("Predefined_Model", b, sep = "_")

        }

      }

    }else {

      Model_GAS <- bn.fit(structure.list, Train_Data, method = "bayes")
      Model_GAS_Gr  <- as.grain(Model_GAS)
      if(is.null(algorithm.names)){
        Models[[1]] <- Model_GAS_Gr
        names(Models)[1] <- "Predefined_Model"

      }else {

        Models[[length(algorithm.names)+1]] <- Model_GAS_Gr
        names(Models)[length(algorithm.names)+1] <- "Predefined_Model"

      }

    }

  }

  return(Models)

}





#' @title Calculate Posterior Distribution of Variable and Store It in List
#' @description posterior distribution of target variable(s) is calculated. \code{\link{information.algorithm}} is used
#'              in order to obtain distributions. Obtained distributions are stored in lists.
#'
#' @param test.data a data frame containing all variables. Columns represent variables. Rows represent observations.
#' @param target variable(s) whose posterior distribution to be calculated.
#' @param query input variable(s) to be used when calculating the posterior distribution of target variables.
#' @param bn.model Bayesian network model including all variables.
#' @param num.iter an integer that shows number of query variables which will be used during calculation of
#'                 posterior distribution of target variable. Default value is `1`.
#'
#' @return a list containing posterior distribution of target variable. List has two index. Posterior distributions
#'         obtained for each target variable and each Bayesian network model. First index is for number of variables
#'         second one is for number of Bayesian network models.
#' @export
#'
#' @seealso \code{\link{mut.inf}} \code{\link{information.algorithm}}
#' @importFrom stats na.omit
#'
#'
#' @examples
Predictions <- function(test.data, target, query, bn.model, num.iter = 1) {

  #Putting the target variable names in same order with the data frame.
  if(length(target) > 1) {

    Variable_Names <- target
    sorted_response <- sort(match(target, names(test.data)))
    Variable_Names <- colnames(test.data[,sorted_response])

  }else {

    Variable_Names <- target

  }


  data <- na.omit(data)

  State_Length <- list()
  for (a in 1:length(Variable_Names)) {

    State_Length[[a]] <- length(bn.model[[1]]$universe$levels[[Variable_Names[a]]])

  }

  Number_Of_Rows <- nrow(data)
  Probability_Informative <- list()
  Prediction <- list()

  for (b in 1:length(State_Length)) {

    Prediction[[b]] <- list()

  }
  #1st index is for each discrete state of variable. 2nd index is for each Bayes model.
  for (d in 1:length(State_Length)) {

    for (c in 1:length(bn.model)) {

      Prediction[[d]][[c]] <- matrix(NA, ncol = State_Length[[d]], nrow = Number_Of_Rows)

    }
  }

  for (f in 1:length(Variable_Names)) {

    for (g in 1:length(bn.model)) {

      for(h in 1:Number_Of_Rows) {

        test_data <- test.data[h,]
        test_data <- apply(test_data,2,as.character)
        Probability_Informative <- information.algorithm(x = test_data, target = Variable_Names[f], query = query, bn.model = bn.model[[g]], num.iter = num.iter)
        Prediction[[f]][[g]][h,] <- unlist(Probability_Informative)

      }

      colnames(Prediction[[f]][[g]]) <- names(unlist(Probability_Informative))

    }

  }

  return(Prediction)

}



#' @title Actual Values of Variables
#' @description Obtain actual values of target variables.Target variable(s) has different discrete states.
#'              In each step, a state becomes boolean value `TRUE`, while the rest of the states remain boolean value `FALSE`.
#'              They are assigned to indexes of list for use in the model validation.
#'
#'
#' @param test.data a data frame containing all variables. Columns represent variables. Rows represent observations.
#' @param target target variable(s)
#' @param bn.model Bayesian network model including all variables.
#'
#' @return a list containing actual values of target variable. List has two index. First index represents number of variables.
#'         Second one represents number of discrete states of variables. For example, if there is one variable that has two different
#'         discrete states, resulting list would be \code{example.list[1][1]} for first state and \code{example.list[1][2]} for second
#'         state respectively.
#' @export
#'
#' @importFrom stats na.omit
#'
#' @seealso \code{\link{ROC_Calculation}}
#'
#' @examples
Actual.Values <- function(test.data, target, bn.model) {
  #Putting the target variable names in same order with the data frame.
  if(length(target) > 1) {

    Variable_Names <- target
    sorted_response <- sort(match(target, names(test.data)))
    Variable_Names <- colnames(test.data[,sorted_response])

  }else {

    Variable_Names <- target

   }

  Test_Data <- as.data.frame(test.data, stringsAsFactors = TRUE)
  Test_Data <- na.omit(Test_Data)

  Output_Names <- list()
  for (a in 1:length(Variable_Names)) {

    Output_Names[[a]] <- Test_Data[,Variable_Names[a]]

  }

  Output_States <- list()
  for (b in 1:length(Variable_Names)) {

    Output_States[[b]] <- list()

  }

  for (c in 1:length(Variable_Names)) {

    for (d in 1:length(bn.model[[1]]$universe$levels[[Variable_Names[c]]])) {

      Output_States[[c]][[d]] = (Output_Names[[c]] == bn.model[[1]]$universe$levels[[Variable_Names[c]]][d])

    }

  }

  return(Output_States)

}



#' @title Obtain ROC curve and AUC
#' @description This function created in order to obtain ROC curve and AUC values of model(s). Posterior Distributions of variable(s)
#'              is compared to actual values of variable(s).
#'
#'
#' @param predictions posterior Distributions of target variable(s).
#' @param a.Values actual values of target variable(s).
#'
#' @return a list containing ROC curve and AUC. List has two index. First index contains AUC values of model(s). Second index contains
#'         contains ROC curve(s) of model(s).
#'
#' @details Posterior distributions of target variable(s) obtained through \code{\link{Predictions}} function is compared to
#'          actual values of target variable(s) obtained through \code{\link{Actual.Values}} function.
#'
#' @export
#'
#' @seealso \code{\link[pROC]{roc}} \code{\link[pROC]{auc}}
#'
#' @importFrom stats na.omit
#' @importFrom pROC roc
#' @importFrom pROC auc
#'
#' @examples
ROC_Calculation <- function(predictions, a.Values){


  AUC <- list()
  ROC <- list()
  AUC_ROC <- list()

  for (a in 1:length(predictions[[1]])) {

    AUC[[a]] <- list()
    ROC[[a]] <- list()

  }

  for (b in 1:length(predictions[[1]])) {

    for (c in 1:length(predictions)) {

      AUC[[b]][[c]] <- matrix(nrow = 1, ncol = length(predictions[[c]][[1]][1,]))
      ROC[[b]][[c]] <- list()

    }

  }


  for (d in 1:length(predictions[[1]])) {

    for (f in 1:length(predictions)) {

      for (g in 1:length(predictions[[f]][[1]][1,])) {

        if(all(!a.Values[[f]][[g]]) == TRUE) {

          AUC[[d]][[f]][,g] = 0
          ROC[[d]][[f]][[g]] = 0

        }else {

          AUC[[d]][[f]][,g] <- auc(roc(a.Values[[f]][[g]],
                                       predictions[[f]][[d]][,g],
                                       ci = TRUE))

          ROC[[d]][[f]][[g]] <- roc(a.Values[[f]][[g]],
                                    predictions[[f]][[d]][,g],
                                    ci = TRUE)
        }

      }

    }

  }

  AUC <- t(sapply(sapply(AUC, rbind), rbind))

  for (x in 1:2) {

    AUC_ROC[[x]] <- list()

  }

  AUC_ROC[[1]] <- AUC
  AUC_ROC[[2]] <- ROC
  return(AUC_ROC)

}








#' @title Cross Validation Method for Bayesian Network Models
#' @description This function is used to perform k-fold cross validation method for a Bayesian network models.
#'
#' @param fold an integer value shows the fold number.
#' @param data.tt a data frame containing the variable(s) in the Bayesian network model(s). Values of the variables
#'                in the data frame are assumed to be discrete.
#' @param target target variable(s) whose actual values and posterior distributions will be obtained.
#' @param input input variable(s) to be used when calculating the posterior distribution of target variables.
#' @param num.iter an integer that shows number of query variables which will be used during calculation of
#'                 posterior distribution of target variable. Default value is `1`.
#' @param str.algorithms a list containing different algorithms \code{\link[bnlearn]{structure-learning}} that are used to learn
#'                       different Bayesian networks. The default value is `NULL`.
#' @param structure.list If there is pre-defined structure of Bayesian network, parameters of it are learned. The default value is `NULL`.
#'
#' @details by using \code{\link{Predictions}} and \code{\link{Actual.Values}} functions posterior distribution and actual values of target
#'          variable(s) obtained. With the \code{\link{Model_Learning}} function Bayesian network model(s) is learned. In the last step
#'          ROC curves and AUC values are obtained through \code{\link{ROC_Calculation}} and model performances are evaluated. K-fold
#'          cross validation method is applied.
#'
#' @return a list containing AUC and ROC curve results. List has two index. First index consists of data frame containing AUC results of
#'         all folds. Second index consists of ROC curve results of all folds.
#' @export
#'
#' @references M. Stone, “Cross-Validatory Choice and Assessment of Statistical Predictions,” J. R. Stat. Soc. Ser. B, vol. 36, no. 2, pp. 111–133, 1974.
#'
#' @seealso \code{\link{information.algorithm}} \code{\link[bnlearn]{bn.cv}}
#' @importFrom stats na.omit
#'
#'
#' @examples
Cross_Validation <- function(fold, data.tt, target, input, num.iter, str.algorithms = NULL, structure.list = NULL){

  if(length(target) > 1) {
    Variable_Names <- target
    sorted_response <- sort(match(target, names(data.tt)))
    Variable_Names <- colnames(data.tt[,sorted_response])
  }else {
    Variable_Names <- target
   }

  Factorized_data <- apply(data.tt, 2, as.factor)
  Factorized_data <- as.data.frame(Factorized_data, stringsAsFactors = TRUE)

  CV_Index <- list()
  for (a in 1:fold) {

    CV_Index[[a]] <- rep(a, round(nrow(Factorized_data)/fold))

  }

  CV_Index <- unlist(CV_Index)
  CV_Results <- list()

  for (x in 1:fold) {

    CV_Results[[x]] <- list()

  }
  #Cross validatation loop
  for (k_fold in 1:fold) {

    TrainIndex = CV_Index != k_fold
    TestIndex  = CV_Index == k_fold

    train_Data  =  Factorized_data[TrainIndex,]
    test_Data   =  Factorized_data[TestIndex,]

    train_Data <- na.omit(train_Data)
    test_Data <- na.omit(test_Data)


    Models_K <- Model_Learning(train.data = train_Data, str.algorithms, structure.list)
    Predictions_K <- Predictions(test.data = test_Data, target = Variable_Names, query = input, bn.model = Models_K, num.iter = num.iter)
    Real_Values_K <- Actual.Values(test.data = test_Data, target = Variable_Names, bn.model = Models_K)
    AUROC_K <- ROC_Calculation(Predictions_K, Real_Values_K)

    for (b in 1:length(Predictions_K[[1]])) {

      for (c in 1:length(Predictions_K)) {

        for (d in 1:length(Predictions_K[[c]][[1]][1,])) {

          names(AUROC_K[[2]])[b] <- names(Models_K)[b]
          names(AUROC_K[[2]][[b]])[c] <- Variable_Names[c]
          names(AUROC_K[[2]][[b]][[c]])[d] <- Models_K[[1]]$universe$levels[[Variable_Names[c]]][d]

        }

      }

    }

    f <- 1
    R_Names <- c()

    for (g in 1:length(Models_K)) {

      for (h in 1:length(Variable_Names)) {


        R_Names[f] <- paste(Variable_Names[h], names(Models_K)[g], sep = "_")
        f = f + 1
        if(f == (length(Models_K)*length(Variable_Names)+1))
          break

      }

    }

    rownames(AUROC_K[[1]]) <- R_Names
    C_Names <- list()

    for (i in 1:length(Variable_Names)) {

      C_Names[[i]] <- Models_K[[1]]$universe$levels[[Variable_Names[i]]]

    }

    C_Names <- unique(unlist(C_Names))
    colnames(AUROC_K[[1]]) <- paste(C_Names, "CV", k_fold, sep = "_")
    CV_Results[[k_fold]] <- AUROC_K

  }

  AUC_Dataframe <- c()
  for (t in 1:fold) {

    AUC_Dataframe <- cbind(AUC_Dataframe, CV_Results[[t]][[1]])


  }

  ROC_List <- list()
  for (k in 1:fold) {

    ROC_List[[k]] <- CV_Results[[k]][[2]]

  }

  AUROC <- list()
  AUROC[[1]] <- AUC_Dataframe
  AUROC[[2]] <- ROC_List

  for (y in 1:fold) {

    names(AUROC[[2]])[y] <- paste("Fold", y, sep = "_")

  }

  return(AUROC)

}



#' @title Convert Continuous Variable into Factor
#' @description Variables which have continuous values converted into categorical variables which have different intervals.
#'
#' @param vector.d a vector containing values of continuous variable.
#' @param cutoff.d cutting points required to divide the vector into intervals
#' @param states.d a vector containing characters which are used to label intervals. Default value is `NULL`. If value is `NULL`,
#'                 intervals are labeled as \code{State_1} \code{State_2}...
#'
#' @details minimum and maximum values of \code{vector.d} equal to the lowest and highest point of result vector respectively.
#'          Interval ranging from minimum point in the \code{vector.d} to the leftmost value defined in the \code{cutoff.d} represents first interval.
#'          Interval ranging from the rightmost value defined in the \code{cutoff.d} to the maximum point in \code{vector.d} represents last interval.
#'
#' @return a vector containing state labels representing different intervals.
#' @export
#'
#' @seealso \code{\link[bnlearn]{discretize}} \code{\link{cut}}
#'
#' @examples
Discretization <- function(vector.d, cutoff.d, states.d = NULL){

  Min_D <- min(vector.d)
  Max_D <- max(vector.d)

  if(is.null(states.d)){

    states.d <- list()
    for (a in 1:(length(cutoff.d)+1)) {

      states.d[[a]] <- paste("State", a, sep = "_")

    }

    states.d <- unlist(states.d)

  }

  Discretized_Vector <- cut(vector.d, breaks = c(Min_D, cutoff.d, Max_D), include.lowest = TRUE, labels = states.d)

  return(Discretized_Vector)

}



