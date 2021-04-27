#' Train MMRDetect
#'
#' @param mutationVariable A list of input variables,"Del_rep_mean","RepIndel_num","MMR_sum","maxcossim"
#' @param classification Sample classification
#' @return trained model

#' @export
MMRDetect.train <- function(mutationVariable, classification) {
  
  
  ## match the data with classification
  trainset = mutationVariable[,c("Sample","Del_rep_mean","RepIndel_num","MMR_sum","maxcossim")]
  trainset = merge(trainset, classification[,c("Sample","MSI_status")], by="Sample")
  
  if(nrow(trainset)<50){
    warning('training set size < 50')  
  }
  
  # normalize RepIndel_num and MMR_sum
  trainset$RepIndel_num <- trainset$RepIndel_num/max(trainset$RepIndel_num)
  trainset$MMR_sum <- trainset$RepIndel_num/max(trainset$MMR_sum)
  
  ## build model with trainset
  trainset$MSI_status<-as.factor(trainset$MSI_status)
  glm_model_logit = stats::glm(MSI_status~., data = trainset, family = binomial(link="logit"))
  glm_model_logit
  
}

