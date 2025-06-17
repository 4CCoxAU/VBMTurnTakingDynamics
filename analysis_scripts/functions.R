# Functions

#Function to extract relevant hypothesis data from models:
extract_relevant_data <- function(model, hypothesisformula, TestDescription) {
  
  hypothesis <- hypothesis(model, hypothesisformula)
  
  relevant_data <- hypothesis$hypothesis[, c("Estimate", "CI.Lower", "CI.Upper", "Evid.Ratio")] %>%
    mutate(Estimate = Estimate * 1000) %>%
    mutate(CI.Lower = CI.Lower * 1000) %>%
    mutate(CI.Upper = CI.Upper * 1000)
  
  relevant_data$TestDescription <- TestDescription
  
  return(relevant_data)
}

#Function to extract relevant hypothesis data from models:
extract_relevant_nonlatency_data <- function(model, hypothesisformula, TestDescription) {
  
  hypothesis <- hypothesis(model, hypothesisformula)
  
  relevant_data <- hypothesis$hypothesis[, c("Estimate", "CI.Lower", "CI.Upper", "Evid.Ratio")]
  
  relevant_data$TestDescription <- TestDescription
  
  return(relevant_data)
}

#Function to extract relevant hypothesis data from models:
extract_relevant_sd_data <- function(model, hypothesisformula, TestDescription) {
  
  hypothesis <- hypothesis(model, hypothesisformula, class = "sd", group = "ID")
  
  relevant_data <- hypothesis$hypothesis[, c("Estimate", "CI.Lower", "CI.Upper", "Evid.Ratio")] %>%
    mutate(Estimate = Estimate * 1000) %>%
    mutate(CI.Lower = CI.Lower * 1000) %>%
    mutate(CI.Upper = CI.Upper * 1000)
  
  relevant_data$TestDescription <- TestDescription
  
  return(relevant_data)
}

#Function to extract relevant hypothesis data from models:
extract_relevant_sdvisit_data <- function(model, hypothesisformula, TestDescription) {
  
  hypothesis <- hypothesis(model, hypothesisformula, class = "sd", group = "Visit")
  
  relevant_data <- hypothesis$hypothesis[, c("Estimate", "CI.Lower", "CI.Upper", "Evid.Ratio")] %>%
    mutate(Estimate = Estimate * 1000) %>%
    mutate(CI.Lower = CI.Lower * 1000) %>%
    mutate(CI.Upper = CI.Upper * 1000)
  
  relevant_data$TestDescription <- TestDescription
  
  return(relevant_data)
}

#Function to extract relevant hypothesis data from models:
extract_relevant_sigma_data <- function(model, hypothesisformula, TestDescription) {
  
  hypothesis <- hypothesis(model, hypothesisformula)
  
  relevant_data <- hypothesis$hypothesis[, c("Estimate", "CI.Lower", "CI.Upper", "Evid.Ratio")] %>%
    mutate(Estimate = exp(Estimate) * 1000) %>%
    mutate(CI.Lower = exp(CI.Lower) * 1000) %>%
    mutate(CI.Upper = exp(CI.Upper) * 1000)
  
  relevant_data$TestDescription <- TestDescription
  
  return(relevant_data)
}

#Function to extract relevant posterior for visits from model:
Visitextract_and_rename <- function(post, diagnosis, task, familiarity) {
  col_names <- paste0("r_Visit[", 1:7, ",", diagnosis, ":", task, ":", familiarity, "]")
  new_names <- paste0("Visit", 1:7)
  
  post[, col_names] %>%
    dplyr::rename_at(vars(all_of(col_names)), ~ new_names)
}

# Function to calculate test-retest reliability with bootstrap CIs
calculate_reliability_ci <- function(data, task_name, familiarity_level, visits, n_bootstrap = 1000) {
  
  # Correlation function for bootstrap
  cor_function <- function(data, indices) {
    boot_data <- data[indices, ]
    visit_cols <- boot_data %>% select(starts_with("Visit"))
    cor_matrix <- cor(visit_cols, use = "pairwise.complete.obs")
    mean(cor_matrix[upper.tri(cor_matrix)], na.rm = TRUE)
  }
  
  # Prepare data for the specific condition
  condition_data <- data %>%
    group_by(Task, Familiarity, ID, Diagnosis, Visit) %>%
    summarise(Latency = mean(Latency, na.rm = T), .groups = "drop") %>%
    filter(Task == task_name, Familiarity == familiarity_level) %>%
    filter(Visit %in% visits) %>%
    select(ID, Diagnosis, Visit, Latency) %>%
    pivot_wider(names_from = Visit, values_from = Latency, names_prefix = "Visit") %>%
    filter(complete.cases(.))
  
  # Bootstrap for each diagnosis group
  results <- map_dfr(c("ASD", "TDC"), function(diag) {
    diag_data <- condition_data %>% 
      filter(Diagnosis == diag) %>% 
      select(-ID, -Diagnosis)
    
    if(nrow(diag_data) < 3) {
      return(tibble(Diagnosis = diag, Estimate = NA, CI.Lower = NA, CI.Upper = NA))
    }
    
    boot_result <- boot(diag_data, cor_function, R = n_bootstrap)
    ci_result <- boot.ci(boot_result, type = "perc")
    
    tibble(
      Diagnosis = diag,
      Estimate = boot_result$t0,
      CI.Lower = ci_result$percent[4],
      CI.Upper = ci_result$percent[5]
    )
  })
  
  return(results)
}

# Updated correlation function with error handling
cross_context_cor <- function(data, indices) {
  # Sample the data using bootstrap indices
  boot_data <- data[indices, ]
  
  # Handle potential duplicates by taking mean
  boot_data <- boot_data %>%
    group_by(ID, Condition) %>%
    summarise(Estimate = mean(Estimate), .groups = "drop")
  
  # Convert to wide format for correlation
  wide_data <- boot_data %>%
    pivot_wider(names_from = Condition, values_from = Estimate) %>%
    select(-ID) # Remove ID column
  
  # Check if we have enough complete cases
  if(nrow(wide_data) < 3 || any(colSums(!is.na(wide_data)) < 3)) {
    return(NA) # Return NA if insufficient data
  }
  
  # Calculate correlation matrix
  cor_matrix <- cor(wide_data, use = "pairwise.complete.obs")
  
  # Check if correlation matrix is valid
  if(any(is.na(cor_matrix))) {
    return(NA)
  }
  
  # Return mean of upper triangle
  mean(cor_matrix[upper.tri(cor_matrix)])
}

# Bootstrap function for one diagnosis group
bootstrap_cross_context <- function(diagnosis_group, n_bootstrap = 1000) {
  # Filter data for the specific diagnosis
  group_data <- clean_data %>% 
    filter(Diagnosis == diagnosis_group)
  
  # Run bootstrap
  boot_result <- boot(data = group_data, 
                      statistic = cross_context_cor, 
                      R = n_bootstrap)
  
  # Calculate confidence intervals
  ci_result <- boot.ci(boot_result, type = "perc")
  
  # Return results
  tibble(
    Diagnosis = diagnosis_group,
    Estimate = boot_result$t0,
    CI.Lower = ci_result$percent[4],
    CI.Upper = ci_result$percent[5]
  )
}

extract_all_hypothesis_tests_m1_child <- function(model) {
  # Overall Estimates in Autism Group
  OverallLatencyThreeConditionsASD <- extract_relevant_data(model,
                                                            "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar +DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
                                                            "OverallLatencyThreeConditionsASD")
  
  MGOverallLatencyASD <- extract_relevant_data(model,
                                               "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > 0",
                                               "MGOverallLatencyASD>0")
  
  Q_F_OverallLatencyASD <- extract_relevant_data(model,
                                                 "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > 0",
                                                 "Q_F_OverallLatencyASD>0")
  
  Q_UF_OverallLatencyASD <- extract_relevant_data(model,
                                                  "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                  "Q_UF_OverallLatencyASD>0")
  
  # Overall Estimates in TD Group
  OverallLatencyThreeConditionsTD <- extract_relevant_data(model,
                                                           "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar +DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3 > 0",
                                                           "OverallLatencyThreeConditionsTD")
  
  MGOverallLatencyTDC <- extract_relevant_data(model,
                                               "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > 0",
                                               "MGOverallLatencyTDC>0")
  
  Q_F_OverallLatencyTDC <- extract_relevant_data(model,
                                                 "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar > 0",
                                                 "Q_F_OverallLatencyTDC>0")
  
  Q_UF_OverallLatencyTDC <- extract_relevant_data(model,
                                                  "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                  "Q_UF_OverallLatencyTDC>0")
  
  # Overall Differences between Autism Group and Typical Development Group
  OverallLatencyThreeConditionsASDFasterTD <- extract_relevant_data(model,
                                                                    "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar +DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 
                                <
                                (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar +DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3",
                                                                    "OverallLatencyThreeConditionsASD<TD")
  
  
  # Overall values for Beta parameter
  BetaParameterASD_MG <- extract_relevant_nonlatency_data(model,
                                                          "beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                          "BetaParameterASD_MG")
  
  BetaParameterASD_Q_F <- extract_relevant_nonlatency_data(model,
                                                           "beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < 0",
                                                           "BetaParameterASD_Q_F")
  
  BetaParameterASD_Q_UF <- extract_relevant_nonlatency_data(model,
                                                            "beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < 0",
                                                            "BetaParameterASD_Q_UF")
  
  BetaParameterASD_Overall <- extract_relevant_nonlatency_data(model,
                                                               "(beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar +
                 beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar +
                 beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar) / 3 < 0",
                                                               "BetaParameterASD_Overall")
  
  # Differences, Typical Development Group with Long Pauses
  BetaParameterTDC_MG <- extract_relevant_nonlatency_data(model,
                                                          "beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                          "BetaParameterTDC_MG")
  
  BetaParameterTDC_Q_F <- extract_relevant_nonlatency_data(model,
                                                           "beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar < 0",
                                                           "BetaParameterTDC_Q_F")
  
  BetaParameterTDC_Q_UF <- extract_relevant_nonlatency_data(model,
                                                            "beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar < 0",
                                                            "BetaParameterTDC_Q_UF")
  
  BetaParameterTDC_Overall <- extract_relevant_nonlatency_data(model,
                                                               "(beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar +
                 beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar +
                 beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3 < 0",
                                                               "BetaParameterTDC_Overall")
  
  
  
  
  WithinGroupBetaParameterASD_MG_Slower_Q_F <- extract_relevant_nonlatency_data(model,
                                                                                "beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                                "WithinGroupBetaParameterASD_MG_Slower_Q_F")
  
  WithinGroupBetaParameterASD_Q_F_Slower_Q_UF <- extract_relevant_nonlatency_data(model,
                                                                                  "beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                                  "WithinGroupBetaParameterASD_Q_F_Slower_Q_UF")
  
  WithinGroupBetaParameterTDC_MG_Slower_Q_F <- extract_relevant_nonlatency_data(model,
                                                                                "beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                                "WithinGroupBetaParameterTDC_MG_Slower_Q_F")
  
  WithinGroupBetaParameterTDC_Q_F_Slower_Q_UF <- extract_relevant_nonlatency_data(model,
                                                                                  "beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                                  "WithinGroupBetaParameterTDC_Q_F_Slower_Q_UF")
  
  DifferenceInDifferenceBetaParameterBetweenGroupsMGBiggerQ_F <- extract_relevant_nonlatency_data(model,
                                                                                                  "(beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar - beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) > 
                 (beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar - beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar)",
                                                                                                  "DifferenceInDifferenceBetaParameterBetweenGroupsMGBiggerQ_F")
  
  DifferenceInDifferenceBetaParameterBetweenGroupsFamiliarityBiggerUnfamiliarity <- extract_relevant_nonlatency_data(model,
                                                                                                                     "(beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar - beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar) < 
                 (beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar - beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)",
                                                                                                                     "DifferenceInDifferenceBetaParameterBetweenGroupsFamiliarityBiggerUnfamiliarity")
  
  # Differences between Autism Group and Typical Development Group with Long Pauses
  OverallBetaParameterASDLowerTD <- extract_relevant_nonlatency_data(model,
                                                                     "(beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar + beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar) / 3 < ((beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar + beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3)",
                                                                     "OverallBetaParameterASD<TD")
  
  BetaParameterASDFasterTD_MG <- extract_relevant_nonlatency_data(model,
                                                                  "beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar < beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar",
                                                                  "BetaParameter_ASD_MG<TD_MG")
  
  BetaParameterASDFasterTD_Q_F <- extract_relevant_nonlatency_data(model,
                                                                   "beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                   "BetaParameter_ASD_MG<TD_Q_F")
  
  BetaParameterASDFasterTD_Q_UF <- extract_relevant_nonlatency_data(model,
                                                                    "beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                    "BetaParameter_ASD_MG<TD_Q_UF")
  
  # Differences within Autism Group across Conditions:
  WithinASDGroupMGSlowerQuestionsFamiliar <- extract_relevant_data(model,
                                                                   "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                   "WithinASDGroupMG>QuestionsFamiliar")
  
  # Differences within TD Group across Conditions:
  WithinTDGroupMGSlowerQuestionsFamiliar <- extract_relevant_data(model,
                                                                  "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                  "WithinTDGroupMG>QuestionsFamiliar")
  
  # Differences within Autism Group across Familiarity:
  WithinAutismGroupQuestionsUnfamiliarSlowerFamiliar <- extract_relevant_data(model,
                                                                              "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                              "WithinAutismGroupQuestionsUnfamiliar>Familiar")
  
  # Differences within TD Group across Familiarity:
  WithinTDGroupQuestionsUnfamiliarSlowerFamiliar <- extract_relevant_data(model,
                                                                          "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                          "WithinTDGroupQuestionsUnfamiliar>Familiar")
  
  # Differences between Autism Group and Typical Development Group Within MatchingGame
  WithinConditionsMatchingGameASDFasterTD <- extract_relevant_data(model,
                                                                   "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar < DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar",
                                                                   "WithinConditionsMatchingGameASD<TD")
  
  # Differences between Autism Group and Typical Development Group Within Familiar Questions:
  WithinConditionsQuestionsFamiliarASDFasterTD <- extract_relevant_data(model,
                                                                        "DiagnosisASD:TaskQuestions:FamiliarityFamiliar < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                        "WithinConditionsQuestionsFamiliarASD<TD")
  
  # Differences between Autism Group and Typical Development Group Within Unfamiliar Questions:
  WithinConditionsQuestionsUnfamiliarASDFasterTD <- extract_relevant_data(model,
                                                                          "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                          "WithinConditionsQuestionsUnfamiliarASD<TD")
  
  # Differences between Autism Group and Typical Development Group in Differences within Familiarity
  DifferenceInDifferenceBetweenGroupsUnfamiliarBiggerFamiliarity <- extract_relevant_data(model,
                                                                                          "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) < 
                      (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                          "DifferenceInDifferenceBetweenGroupsUnfamiliar>Familiarity")
  
  # Differences between Autism Group and Typical Development Group in Differences within Conditions
  DifferenceInDifferenceBetweenGroupsMGBiggerThanConversations <- extract_relevant_data(model,
                                                                                        "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar) > 
                      (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                        "DifferenceInDifferenceBetweenGroupsMG>Conversations")
  
  # Overall Subject Variability in Autismm Group
  SubjectVariabilityASD <- extract_relevant_sd_data(model, 
                                                    "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisASD + TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD)/3 > 0",
                                                    "SubjectVariabilityASD")
  
  SubjectVariabilityASD_MG <- extract_relevant_sd_data(model, 
                                                       "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD > 0",
                                                       "SubjectVariabilityASD_MG")
  
  SubjectVariabilityASD_Q_F <- extract_relevant_sd_data(model, 
                                                        "TaskQuestions:FamiliarityFamiliar:DiagnosisASD > 0",
                                                        "SubjectVariabilityASD_Q_F")
  
  SubjectVariabilityASD_Q_UF <- extract_relevant_sd_data(model, 
                                                         "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD > 0",
                                                         "SubjectVariabilityASD_Q_UF")
  
  # Overall Subject Variability in TD Group
  SubjectVariabilityTD <- extract_relevant_sd_data(model, 
                                                   "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisTDC + TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC)/3 > 0",
                                                   "SubjectVariabilityTD")
  
  SubjectVariabilityTD_MG <- extract_relevant_sd_data(model, 
                                                      "TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC > 0",
                                                      "SubjectVariabilityTD_MG")
  
  SubjectVariabilityTD_Q_F <- extract_relevant_sd_data(model, 
                                                       "TaskQuestions:FamiliarityFamiliar:DiagnosisTDC > 0",
                                                       "SubjectVariabilityTD_Q_F")
  
  SubjectVariabilityTD_Q_UF <- extract_relevant_sd_data(model, 
                                                        "TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC > 0",
                                                        "SubjectVariabilityTD_Q_UF")
  
  # Overall Subject Variability in TD Group
  SubjectVariabilityASDSmallerThanTD <- extract_relevant_sd_data(model,
                                                                 "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisASD + TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD) / 3 <
                                              (TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisTDC + TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC) / 3",
                                                                 "SubjectVariabilityASD<TD")
  
  
  SubjectVariabilityASDSmallerThanTD_MG <- extract_relevant_sd_data(model,
                                                                    "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD  <
                                              TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC",
                                                                    "SubjectVariabilityASD_MG<TD_MG")
  
  SubjectVariabilityASDSmallerThanTD_Q_F <- extract_relevant_sd_data(model,
                                                                     "TaskQuestions:FamiliarityFamiliar:DiagnosisASD  <
                                              TaskQuestions:FamiliarityFamiliar:DiagnosisTDC",
                                                                     "SubjectVariabilityASD_Q_F<TD_Q_F")
  
  SubjectVariabilityASDSmallerThanTD_Q_UF <- extract_relevant_sd_data(model,
                                                                      "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD  <
                                              TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC",
                                                                      "SubjectVariabilityASD_Q_UF<TD_Q_UF")
  
  SubjectVariabilityASD_MG_SmallerThan_Q_F <- extract_relevant_sd_data(model,
                                                                       "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD > TaskQuestions:FamiliarityFamiliar:DiagnosisASD",
                                                                       "SubjectVariabilityASD_MG_SmallerThan_Q_F")
  
  SubjectVariabilityASD_Q_F_SmallerThan_Q_UF <- extract_relevant_sd_data(model,
                                                                         "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD > TaskQuestions:FamiliarityFamiliar:DiagnosisASD",
                                                                         "SubjectVariabilityASD_Q_F_SmallerThan_Q_UF")
  
  SubjectVariabilityTDC_MG_SmallerThan_Q_F <- extract_relevant_sd_data(model,
                                                                       "TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC < TaskQuestions:FamiliarityFamiliar:DiagnosisTDC",
                                                                       "SubjectVariabilityTDC_MG_SmallerThan_Q_F")
  
  SubjectVariabilityTDC_Q_F_SmallerThan_Q_UF <- extract_relevant_sd_data(model,
                                                                         "TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC > TaskQuestions:FamiliarityFamiliar:DiagnosisTDC",
                                                                         "SubjectVariabilityTDC_Q_F_SmallerThan_Q_UF")
  
  DifferenceInDifferenceSubjectVariabilityUnfamiliarSmallerFamiliarity <- extract_relevant_sd_data(model,
                                                                                                   "(TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD - TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC) 
                                              < 
                                              (TaskQuestions:FamiliarityFamiliar:DiagnosisASD - TaskQuestions:FamiliarityFamiliar:DiagnosisTDC)",
                                                                                                   "DifferenceInDifferenceSubjectVariabilityUnfamiliarBiggerFamiliarity")
  
  DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions <- extract_relevant_sd_data(model,
                                                                                                  "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC - TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD) 
                                              < 
                                              (TaskQuestions:FamiliarityFamiliar:DiagnosisTDC - TaskQuestions:FamiliarityFamiliar:DiagnosisASD)",
                                                                                                  "DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions")
  
  # Overall Visit Variability across visits in Autism Group
  VisitVariabilityASDOverall <- extract_relevant_sdvisit_data(model, 
                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
                                                              "VisitVariabilityASDOverall")
  
  # Overall Visit Variability across visits in TD Group
  VisitVariabilityTDCOverall <- extract_relevant_sdvisit_data(model, 
                                                              "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
                                                              "VisitVariabilityTDCOverall")
  
  
  VisitVariabilityASDOverallHigherTDSOverall <- extract_relevant_sdvisit_data(model, 
                                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 >
                                               (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3",
                                                                              "VisitVariabilityASDOverallLowerTDSOverall")
  
  # Subject Variability across visits MG in Autism Group
  VisitVariabilityASDMG <- extract_relevant_sdvisit_data(model, 
                                                         "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                         "VisitVariabilityASDMG")
  
  # Overall Subject Variability across visits in TD Group
  VisitVariabilityTDMG <- extract_relevant_sdvisit_data(model, 
                                                        "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                        "VisitVariabilityTDMG")
  
  # Subject Variability across visits MG in Autism Group
  VisitVariabilityASD_Q_F <- extract_relevant_sdvisit_data(model, 
                                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > 0",
                                                           "VisitVariabilityASD_Q_F")
  
  # Overall Subject Variability across visits in TD Group
  VisitVariabilityTDC_Q_F <- extract_relevant_sdvisit_data(model, 
                                                           "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar > 0",
                                                           "VisitVariabilityTDC_Q_F")
  
  # Subject Variability across visits MG in Autism Group
  VisitVariabilityASD_Q_UF <- extract_relevant_sdvisit_data(model, 
                                                            "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                            "VisitVariabilityASD_Q_UF")
  
  # Overall Subject Variability across visits in TD Group
  VisitVariabilityTDC_Q_UF <- extract_relevant_sdvisit_data(model, 
                                                            "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                            "VisitVariabilityTDC_Q_UF")
  
  # Overall Subject Variability across visits, differences
  VisitVariabilityASDMGSlowerThanTDMG <- extract_relevant_sdvisit_data(model, 
                                                                       "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar < DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar",
                                                                       "VisitVariabilityASDMG>TDMG")
  
  # Overall Subject Variability across visits, differences
  VisitVariabilityASD_Q_F_SlowerThanTD_Q_F <- extract_relevant_sdvisit_data(model, 
                                                                            "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                            "VisitVariabilityASD_Q_F_SlowerThanTD_Q_F")
  
  # Overall Subject Variability across visits, differences
  VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF <- extract_relevant_sdvisit_data(model, 
                                                                              "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                              "VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF")
  
  #DifferenceInDifference:
  DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity <- extract_relevant_sdvisit_data(model, 
                                                                                                       "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)
                                               < (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                                       "DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity")
  
  DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG <- extract_relevant_sdvisit_data(model, 
                                                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar)
                                               > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                              "DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG")
  
  
  # Overall Sigma Variability across visits in Autism Group
  SigmaASDOverall <- extract_relevant_sigma_data(model, 
                                                 "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 < 0",
                                                 "SigmaASDOverall")
  
  # Overall Sigma Variability across visits in TD Group
  SigmaTDCOverall <- extract_relevant_sigma_data(model, 
                                                 "(sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) /3 < 0",
                                                 "SigmaTDCOverall")
  
  
  SigmaASDOverallHigherTDSOverall <- extract_relevant_nonlatency_data(model, 
                                                                      "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 >
                                               (sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3",
                                                                      "SigmaASDOverallHigherTDSOverall")
  
  # Subject Variability across visits MG in Autism Group
  SigmaASDMG <- extract_relevant_sigma_data(model, 
                                            "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar < 0",
                                            "SigmaASDMG")
  
  # Overall Subject Variability across visits in TD Group
  SigmaTDMG <- extract_relevant_sigma_data(model, 
                                           "sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar < 0",
                                           "SigmaTDMG")
  
  # Subject Variability across visits MG in Autism Group
  SigmaASD_Q_F <- extract_relevant_sigma_data(model, 
                                              "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < 0",
                                              "SigmaASD_Q_F")
  
  # Overall Subject Variability across visits in TD Group
  SigmaTDC_Q_F <- extract_relevant_sigma_data(model, 
                                              "sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar < 0",
                                              "SigmaTDC_Q_F")
  
  
  # Subject Variability across visits MG in Autism Group
  SigmaASD_Q_UF <- extract_relevant_sigma_data(model, 
                                               "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < 0 ",
                                               "SigmaASD_Q_UF")
  
  # Overall Subject Variability across visits in TD Group
  SigmaTDC_Q_UF <- extract_relevant_sigma_data(model, 
                                               "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar < 0",
                                               "SigmaTDC_Q_UF")
  
  # Overall Subject Variability across visits, differences
  SigmaASDMGSlowerThanTDMG <- extract_relevant_nonlatency_data(model, 
                                                               "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar",
                                                               "SigmaASDMGSlowerThanTDMG")
  
  # Overall Subject Variability across visits, differences
  SigmaASD_Q_F_SlowerThanTD_Q_F <- extract_relevant_nonlatency_data(model, 
                                                                    "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                    "SigmaASD_Q_F_SlowerThanTD_Q_F")
  
  # Overall Subject Variability across visits, differences
  SigmaASD_Q_UF_SlowerThanTD_Q_UF <- extract_relevant_nonlatency_data(model, 
                                                                      "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                      "SigmaASD_Q_UF_SlowerThanTD_Q_UF")
  
  # Sigma across visits, differences
  SigmaASD_MG_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
                                                                 "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit < 0",
                                                                 "SigmaASD_MG_VisitLowerZero")
  
  SigmaTDC_MG_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
                                                                 "sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit > 0",
                                                                 "SigmaTDC_MG_VisitLowerZero")
  
  SigmaASD_MG_VisitLowerTDC <- extract_relevant_nonlatency_data(model, 
                                                                "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit < sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
                                                                "SigmaASD_MG_VisitLowerTDC")
  
  SigmaASD_Overall_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
                                                                      "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit + 
                                               sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit +
                                               sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit) / 3 < 0",
                                                                      "SigmaASD_Overall_VisitLowerZero")
  
  SigmaTDC_Overall_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
                                                                      "(sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit + 
                                               sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit +
                                               sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit) / 3 < 0",
                                                                      "SigmaTDC_Overall_VisitLowerZero")
  
  # Sigma across visits, differences
  SigmaASD_Q_F_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
                                                                  "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit > 0",
                                                                  "SigmaASD_Q_F_VisitLowerZero")
  
  SigmaTDC_Q_F_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
                                                                  "sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit < 0",
                                                                  "SigmaTDC_Q_F_VisitLowerZero")
  
  SigmaASD_Q_F_VisitLowerTDC <- extract_relevant_nonlatency_data(model, 
                                                                 "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit > sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
                                                                 "SigmaASD_Q_F_VisitLowerTDC")
  
  # Sigma across visits, differences
  SigmaASD_Q_UF_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
                                                                   "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < 0",
                                                                   "SigmaASD_Q_UF_VisitLowerZero")
  
  SigmaTDC_Q_UF_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
                                                                   "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > 0",
                                                                   "SigmaTDC_Q_UF_VisitLowerZero")
  
  SigmaASD_Q_UF_VisitLowerTDC <- extract_relevant_nonlatency_data(model, 
                                                                  "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
                                                                  "SigmaASD_Q_UF_VisitLowerTDC")
  
  SigmaVisitOverall_ASDLowerTDC <- extract_relevant_nonlatency_data(model, 
                                                                    "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit + 
                                               sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit +
                                               sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit) / 3 > (sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit + 
                                               sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit +
                                               sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit) /3",
                                                                    "SigmaVisitOverall_ASDLowerTDC")
  
  #DifferenceInDifference:
  DifferenceInDifferenceSigmaASDFamiliarityBiggerTDUnfamiliarity <- extract_relevant_nonlatency_data(model, 
                                                                                                     "(sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)
                                               > (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                                     "DifferenceInDifferenceSigmaASDFamiliarityBiggerTDUnfamiliarity")
  
  DifferenceInDifferenceSigmaASDMGSmallerThanTDMG <- extract_relevant_nonlatency_data(model, 
                                                                                      "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar)
                                               > (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                      "DifferenceInDifferenceSigmaASDMGSmallerThanTDMG")
  
  #WithinGroups
  WithinGroupSigmaASD_Q_F_Faster_UF <- extract_relevant_nonlatency_data(model, 
                                                                        "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaASD_Q_FFaster_UF")
  
  WithinGroupSigmaTDC_Q_F_Faster_UF <- extract_relevant_nonlatency_data(model, 
                                                                        "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaTDC_Q_F_Faster_UF")
  
  WithinGroupSigmaASD_Q_F_Faster_MG <- extract_relevant_nonlatency_data(model, 
                                                                        "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaASD_Q_F_Faster_MG")
  
  WithinGroupSigmaTDC_Q_F_Faster_MG <- extract_relevant_nonlatency_data(model, 
                                                                        "sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar < sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaTDC_Q_F_Faster_MG")
  
  
  WithinGroupSigmaASD_Q_F_Lower_UF_Visit <- extract_relevant_nonlatency_data(model, 
                                                                             "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit",
                                                                             "WithinGroupSigmaASD_Q_F_Lower_UF_Visit")
  
  WithinGroupSigmaTDC_Q_F_Lower_UF_Visit <- extract_relevant_nonlatency_data(model, 
                                                                             "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit",
                                                                             "WithinGroupSigmaTDC_Q_F_Lower_UF_Visit")
  
  WithinGroupSigmaASD_Q_F_Lower_MG_Visit <- extract_relevant_nonlatency_data(model, 
                                                                             "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit > sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit",
                                                                             "WithinGroupSigmaASD_Q_F_Lower_MG_Visit")
  
  WithinGroupSigmaTDC_Q_F_Lower_MG_Visit <- extract_relevant_nonlatency_data(model, 
                                                                             "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
                                                                             "WithinGroupSigmaTDC_Q_F_Lower_MG_Visit")
  
  
  SigmaTDC_Q_UF_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
                                                                   "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > 0",
                                                                   "SigmaTDC_Q_UF_VisitLowerZero")
  
  SigmaASD_Q_UF_VisitLowerTDC <- extract_relevant_nonlatency_data(model, 
                                                                  "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
                                                                  "SigmaASD_Q_UF_VisitLowerTDC")
  
  
  DifferenceInDifferenceSigmaVisit_Q_F_Lower_Q_UF <- extract_relevant_nonlatency_data(model, 
                                                                                      "(sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit - sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit) <
                                               (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit)", "DifferenceInDifferenceSigmaVisit_Q_F_Lower_Q_UF")
  
  
  DifferenceInDifferenceSigmaVisit_MG_Lower_Q_F <- extract_relevant_nonlatency_data(model, 
                                                                                    "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit - sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit) <
                                               (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit)", "DifferenceInDifferenceSigmaVisit_MG_Lower_Q_F")
  
  # Calculate bootstrapped correlations between visits
  parent_reliability_ci <- calculate_reliability_ci(
    data = d, 
    task_name = "Questions", 
    familiarity_level = "Familiar", 
    visits = c(2, 4, 6)
  ) %>%
    mutate(Evid.Ratio = NA,
           TestDescription = case_when(
             Diagnosis == "ASD" ~ "CorrelationAcrossVisitsASDFamiliarQuestions",
             Diagnosis == "TDC" ~ "CorrelationAcrossVisitsTDFamiliarQuestions")) %>%
    select(Estimate, CI.Lower, CI.Upper, Evid.Ratio, TestDescription)
  
  experimenter_reliability_ci <- calculate_reliability_ci(
    data = d, 
    task_name = "Questions", 
    familiarity_level = "Unfamiliar", 
    visits = c(1, 3, 5, 7)
  ) %>%
    mutate(Evid.Ratio = NA,
           TestDescription = case_when(
             Diagnosis == "ASD" ~ "CorrelationAcrossVisitsASDUnfamiliarQuestions",
             Diagnosis == "TDC" ~ "CorrelationAcrossVisitsTDUnfamiliarQuestions")) %>%
    select(Estimate, CI.Lower, CI.Upper, Evid.Ratio, TestDescription)
  
  matching_reliability_ci <- calculate_reliability_ci(
    data = d, 
    task_name = "MatchingGame", 
    familiarity_level = "Familiar", 
    visits = c(1, 2, 3, 4, 5, 6, 7)
  ) %>%
    mutate(Evid.Ratio = NA,
           TestDescription = case_when(
             Diagnosis == "ASD" ~ "CorrelationAcrossVisitsASDMG",
             Diagnosis == "TDC" ~ "CorrelationAcrossVisitsTDMG")) %>%
    select(Estimate, CI.Lower, CI.Upper, Evid.Ratio, TestDescription)
  
  CorrelationSummaryStatistics <- rbind(parent_reliability_ci, experimenter_reliability_ci, matching_reliability_ci)
  
  HypothesisResultsChildLatenciesModel1 <- rbind(OverallLatencyThreeConditionsASD,
                                                 OverallLatencyThreeConditionsTD,
                                                 MGOverallLatencyASD,
                                                 Q_F_OverallLatencyASD,
                                                 Q_UF_OverallLatencyASD,
                                                 MGOverallLatencyTDC,
                                                 Q_F_OverallLatencyTDC,
                                                 Q_UF_OverallLatencyTDC,
                                                 BetaParameterASD_MG,
                                                 BetaParameterASD_Q_F,
                                                 BetaParameterASD_Q_UF,
                                                 BetaParameterASD_Overall,
                                                 BetaParameterTDC_MG,
                                                 BetaParameterTDC_Q_F,
                                                 BetaParameterTDC_Q_UF,
                                                 BetaParameterTDC_Overall,
                                                 OverallBetaParameterASDLowerTD,
                                                 BetaParameterASDFasterTD_MG,
                                                 BetaParameterASDFasterTD_Q_F,
                                                 BetaParameterASDFasterTD_Q_UF,
                                                 WithinGroupBetaParameterASD_MG_Slower_Q_F,
                                                 WithinGroupBetaParameterTDC_MG_Slower_Q_F,
                                                 DifferenceInDifferenceBetaParameterBetweenGroupsMGBiggerQ_F,
                                                 WithinGroupBetaParameterASD_Q_F_Slower_Q_UF,
                                                 WithinGroupBetaParameterTDC_Q_F_Slower_Q_UF,
                                                 DifferenceInDifferenceBetaParameterBetweenGroupsFamiliarityBiggerUnfamiliarity,
                                                 OverallLatencyThreeConditionsASDFasterTD,
                                                 WithinASDGroupMGSlowerQuestionsFamiliar,
                                                 WithinTDGroupMGSlowerQuestionsFamiliar,
                                                 WithinAutismGroupQuestionsUnfamiliarSlowerFamiliar,
                                                 WithinTDGroupQuestionsUnfamiliarSlowerFamiliar,
                                                 WithinConditionsMatchingGameASDFasterTD,
                                                 WithinConditionsQuestionsFamiliarASDFasterTD,
                                                 WithinConditionsQuestionsUnfamiliarASDFasterTD,
                                                 DifferenceInDifferenceBetweenGroupsUnfamiliarBiggerFamiliarity,
                                                 DifferenceInDifferenceBetweenGroupsMGBiggerThanConversations,
                                                 SubjectVariabilityASD,
                                                 SubjectVariabilityTD,
                                                 SubjectVariabilityASDSmallerThanTD,
                                                 SubjectVariabilityASD_MG,
                                                 SubjectVariabilityASD_Q_F,
                                                 SubjectVariabilityASD_Q_UF,
                                                 SubjectVariabilityTD_MG,
                                                 SubjectVariabilityTD_Q_F,
                                                 SubjectVariabilityTD_Q_UF,
                                                 SubjectVariabilityASDSmallerThanTD_MG,
                                                 SubjectVariabilityASDSmallerThanTD_Q_F,
                                                 SubjectVariabilityASDSmallerThanTD_Q_UF,
                                                 SubjectVariabilityASD_MG_SmallerThan_Q_F,
                                                 SubjectVariabilityASD_Q_F_SmallerThan_Q_UF,
                                                 SubjectVariabilityTDC_MG_SmallerThan_Q_F,
                                                 SubjectVariabilityTDC_Q_F_SmallerThan_Q_UF,
                                                 DifferenceInDifferenceSubjectVariabilityUnfamiliarSmallerFamiliarity,
                                                 DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions,
                                                 VisitVariabilityASDMG,
                                                 VisitVariabilityTDMG,
                                                 VisitVariabilityASD_Q_F,
                                                 VisitVariabilityTDC_Q_F,
                                                 VisitVariabilityASD_Q_UF,
                                                 VisitVariabilityTDC_Q_UF,
                                                 VisitVariabilityASDOverall,
                                                 VisitVariabilityTDCOverall,
                                                 VisitVariabilityASDMGSlowerThanTDMG,
                                                 VisitVariabilityASD_Q_F_SlowerThanTD_Q_F,
                                                 VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF,
                                                 VisitVariabilityASDOverallHigherTDSOverall,
                                                 DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG,
                                                 DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity,
                                                 SigmaASDOverall,
                                                 SigmaTDCOverall,
                                                 SigmaASDOverallHigherTDSOverall,
                                                 SigmaASDMG,
                                                 SigmaTDMG,
                                                 SigmaASD_Q_F,
                                                 SigmaTDC_Q_F,
                                                 SigmaASD_Q_UF,
                                                 SigmaTDC_Q_UF,
                                                 SigmaASDMGSlowerThanTDMG,
                                                 SigmaASD_Q_F_SlowerThanTD_Q_F,
                                                 SigmaASD_Q_UF_SlowerThanTD_Q_UF,
                                                 SigmaASD_MG_VisitLowerZero,
                                                 SigmaTDC_MG_VisitLowerZero,
                                                 SigmaASD_MG_VisitLowerTDC,
                                                 SigmaASD_Q_F_VisitLowerZero,
                                                 SigmaTDC_Q_F_VisitLowerZero,
                                                 SigmaASD_Q_F_VisitLowerTDC,
                                                 SigmaASD_Q_UF_VisitLowerZero,
                                                 SigmaTDC_Q_UF_VisitLowerZero,
                                                 SigmaASD_Q_UF_VisitLowerTDC,
                                                 SigmaASD_Overall_VisitLowerZero,
                                                 SigmaTDC_Overall_VisitLowerZero,
                                                 WithinGroupSigmaASD_Q_F_Lower_UF_Visit,
                                                 WithinGroupSigmaTDC_Q_F_Lower_UF_Visit,
                                                 WithinGroupSigmaASD_Q_F_Lower_MG_Visit,
                                                 WithinGroupSigmaTDC_Q_F_Lower_MG_Visit,
                                                 WithinGroupSigmaASD_Q_F_Faster_UF,
                                                 WithinGroupSigmaTDC_Q_F_Faster_UF,
                                                 WithinGroupSigmaASD_Q_F_Faster_MG,
                                                 WithinGroupSigmaTDC_Q_F_Faster_MG,
                                                 SigmaVisitOverall_ASDLowerTDC,
                                                 DifferenceInDifferenceSigmaVisit_Q_F_Lower_Q_UF,
                                                 DifferenceInDifferenceSigmaVisit_MG_Lower_Q_F,
                                                 DifferenceInDifferenceSigmaASDFamiliarityBiggerTDUnfamiliarity,
                                                 DifferenceInDifferenceSigmaASDMGSmallerThanTDMG,
                                                 CorrelationSummaryStatistics) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  
  return(HypothesisResultsChildLatenciesModel1)
}

extract_all_hypothesis_tests_m1_adult <- function(model) {
  
  # Overall Estimates in Autism Group
  OverallLatencyThreeConditionsASD <- extract_relevant_data(model,
                                                            "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar +DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
                                                            "OverallLatencyThreeConditionsASD")
  
  MGOverallLatencyASD <- extract_relevant_data(model,
                                               "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > 0",
                                               "MGOverallLatencyASD>0")
  
  Q_F_OverallLatencyASD <- extract_relevant_data(model,
                                                 "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > 0",
                                                 "Q_F_OverallLatencyASD>0")
  
  Q_UF_OverallLatencyASD <- extract_relevant_data(model,
                                                  "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                  "Q_UF_OverallLatencyASD>0")
  
  # Overall Estimates in TD Group
  OverallLatencyThreeConditionsTD <- extract_relevant_data(model,
                                                           "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar +DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3 > 0",
                                                           "OverallLatencyThreeConditionsTD")
  
  MGOverallLatencyTDC <- extract_relevant_data(model,
                                               "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > 0",
                                               "MGOverallLatencyTDC>0")
  
  Q_F_OverallLatencyTDC <- extract_relevant_data(model,
                                                 "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar > 0",
                                                 "Q_F_OverallLatencyTDC>0")
  
  Q_UF_OverallLatencyTDC <- extract_relevant_data(model,
                                                  "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                  "Q_UF_OverallLatencyTDC>0")
  
  # Overall Differences between Autism Group and Typical Development Group
  OverallLatencyThreeConditionsASDFasterTD <- extract_relevant_data(model,
                                                                    "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar +DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 
                                <
                                (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar +DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3",
                                                                    "OverallLatencyThreeConditionsASD<TD")
  
  
  # Overall values for Beta parameter
  BetaParameterASD_MG <- extract_relevant_nonlatency_data(model,
                                                          "beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                          "BetaParameterASD_MG")
  
  BetaParameterASD_Q_F <- extract_relevant_nonlatency_data(model,
                                                           "beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < 0",
                                                           "BetaParameterASD_Q_F")
  
  BetaParameterASD_Q_UF <- extract_relevant_nonlatency_data(model,
                                                            "beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < 0",
                                                            "BetaParameterASD_Q_UF")
  
  BetaParameterASD_Overall <- extract_relevant_nonlatency_data(model,
                                                               "(beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar +
                 beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar +
                 beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar) / 3 < 0",
                                                               "BetaParameterASD_Overall")
  
  # Differences, Typical Development Group with Long Pauses
  BetaParameterTDC_MG <- extract_relevant_nonlatency_data(model,
                                                          "beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                          "BetaParameterTDC_MG")
  
  BetaParameterTDC_Q_F <- extract_relevant_nonlatency_data(model,
                                                           "beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar < 0",
                                                           "BetaParameterTDC_Q_F")
  
  BetaParameterTDC_Q_UF <- extract_relevant_nonlatency_data(model,
                                                            "beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar < 0",
                                                            "BetaParameterTDC_Q_UF")
  
  BetaParameterTDC_Overall <- extract_relevant_nonlatency_data(model,
                                                               "(beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar +
                 beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar +
                 beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3 < 0",
                                                               "BetaParameterTDC_Overall")
  
  WithinGroupBetaParameterASD_MG_Slower_Q_F <- extract_relevant_nonlatency_data(model,
                                                                                "beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                                "WithinGroupBetaParameterASD_MG_Slower_Q_F")
  
  WithinGroupBetaParameterASD_Q_F_Slower_Q_UF <- extract_relevant_nonlatency_data(model,
                                                                                  "beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                                  "WithinGroupBetaParameterASD_Q_F_Slower_Q_UF")
  
  WithinGroupBetaParameterTDC_MG_Slower_Q_F <- extract_relevant_nonlatency_data(model,
                                                                                "beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                                "WithinGroupBetaParameterTDC_MG_Slower_Q_F")
  
  WithinGroupBetaParameterTDC_Q_F_Slower_Q_UF <- extract_relevant_nonlatency_data(model,
                                                                                  "beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                                  "WithinGroupBetaParameterTDC_Q_F_Slower_Q_UF")
  
  DifferenceInDifferenceBetaParameterBetweenGroupsMGBiggerQ_F <- extract_relevant_nonlatency_data(model,
                                                                                                  "(beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar - beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) > 
                 (beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar - beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar)",
                                                                                                  "DifferenceInDifferenceBetaParameterBetweenGroupsMGBiggerQ_F")
  
  DifferenceInDifferenceBetaParameterBetweenGroupsFamiliarityBiggerUnfamiliarity <- extract_relevant_nonlatency_data(model,
                                                                                                                     "(beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar - beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar) < 
                 (beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar - beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)",
                                                                                                                     "DifferenceInDifferenceBetaParameterBetweenGroupsFamiliarityBiggerUnfamiliarity")
  
  # Differences between Autism Group and Typical Development Group with Long Pauses
  OverallBetaParameterASDLowerTD <- extract_relevant_nonlatency_data(model,
                                                                     "(beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar + beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar) / 3 < ((beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar + beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3)",
                                                                     "OverallBetaParameterASD<TD")
  
  BetaParameterASDFasterTD_MG <- extract_relevant_nonlatency_data(model,
                                                                  "beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar < beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar",
                                                                  "BetaParameter_ASD_MG<TD_MG")
  
  BetaParameterASDFasterTD_Q_F <- extract_relevant_nonlatency_data(model,
                                                                   "beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                   "BetaParameter_ASD_MG<TD_Q_F")
  
  BetaParameterASDFasterTD_Q_UF <- extract_relevant_nonlatency_data(model,
                                                                    "beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                    "BetaParameter_ASD_MG<TD_Q_UF")
  
  # Differences within Autism Group across Conditions:
  WithinASDGroupMGSlowerQuestionsFamiliar <- extract_relevant_data(model,
                                                                   "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                   "WithinASDGroupMG>QuestionsFamiliar")
  
  # Differences within TD Group across Conditions:
  WithinTDGroupMGSlowerQuestionsFamiliar <- extract_relevant_data(model,
                                                                  "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                  "WithinTDGroupMG>QuestionsFamiliar")
  
  # Differences within Autism Group across Familiarity:
  WithinAutismGroupQuestionsUnfamiliarSlowerFamiliar <- extract_relevant_data(model,
                                                                              "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                              "WithinAutismGroupQuestionsUnfamiliar>Familiar")
  
  # Differences within TD Group across Familiarity:
  WithinTDGroupQuestionsUnfamiliarSlowerFamiliar <- extract_relevant_data(model,
                                                                          "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                          "WithinTDGroupQuestionsUnfamiliar>Familiar")
  
  # Differences between Autism Group and Typical Development Group Within MatchingGame
  WithinConditionsMatchingGameASDFasterTD <- extract_relevant_data(model,
                                                                   "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar < DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar",
                                                                   "WithinConditionsMatchingGameASD<TD")
  
  # Differences between Autism Group and Typical Development Group Within Familiar Questions:
  WithinConditionsQuestionsFamiliarASDFasterTD <- extract_relevant_data(model,
                                                                        "DiagnosisASD:TaskQuestions:FamiliarityFamiliar < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                        "WithinConditionsQuestionsFamiliarASD<TD")
  
  # Differences between Autism Group and Typical Development Group Within Unfamiliar Questions:
  WithinConditionsQuestionsUnfamiliarASDFasterTD <- extract_relevant_data(model,
                                                                          "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                          "WithinConditionsQuestionsUnfamiliarASD<TD")
  
  # Differences between Autism Group and Typical Development Group in Differences within Familiarity
  DifferenceInDifferenceBetweenGroupsUnfamiliarBiggerFamiliarity <- extract_relevant_data(model,
                                                                                          "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) < 
                      (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                          "DifferenceInDifferenceBetweenGroupsUnfamiliar>Familiarity")
  
  # Differences between Autism Group and Typical Development Group in Differences within Conditions
  DifferenceInDifferenceBetweenGroupsMGBiggerThanConversations <- extract_relevant_data(model,
                                                                                        "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar) > 
                      (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                        "DifferenceInDifferenceBetweenGroupsMG>Conversations")
  
  # Overall Subject Variability in Autismm Group
  SubjectVariabilityASD <- extract_relevant_sd_data(model, 
                                                    "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisASD + TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD)/3 > 0",
                                                    "SubjectVariabilityASD")
  
  SubjectVariabilityASD_MG <- extract_relevant_sd_data(model, 
                                                       "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD > 0",
                                                       "SubjectVariabilityASD_MG")
  
  SubjectVariabilityASD_Q_F <- extract_relevant_sd_data(model, 
                                                        "TaskQuestions:FamiliarityFamiliar:DiagnosisASD > 0",
                                                        "SubjectVariabilityASD_Q_F")
  
  SubjectVariabilityASD_Q_UF <- extract_relevant_sd_data(model, 
                                                         "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD > 0",
                                                         "SubjectVariabilityASD_Q_UF")
  
  # Overall Subject Variability in TD Group
  SubjectVariabilityTD <- extract_relevant_sd_data(model, 
                                                   "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisTDC + TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC)/3 > 0",
                                                   "SubjectVariabilityTD")
  
  SubjectVariabilityTD_MG <- extract_relevant_sd_data(model, 
                                                      "TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC > 0",
                                                      "SubjectVariabilityTD_MG")
  
  SubjectVariabilityTD_Q_F <- extract_relevant_sd_data(model, 
                                                       "TaskQuestions:FamiliarityFamiliar:DiagnosisTDC > 0",
                                                       "SubjectVariabilityTD_Q_F")
  
  SubjectVariabilityTD_Q_UF <- extract_relevant_sd_data(model, 
                                                        "TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC > 0",
                                                        "SubjectVariabilityTD_Q_UF")
  
  # Overall Subject Variability in TD Group
  SubjectVariabilityASDSmallerThanTD <- extract_relevant_sd_data(model,
                                                                 "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisASD + TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD) / 3 <
                                              (TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisTDC + TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC) / 3",
                                                                 "SubjectVariabilityASD<TD")
  
  
  SubjectVariabilityASDSmallerThanTD_MG <- extract_relevant_sd_data(model,
                                                                    "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD  <
                                              TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC",
                                                                    "SubjectVariabilityASD_MG<TD_MG")
  
  SubjectVariabilityASDSmallerThanTD_Q_F <- extract_relevant_sd_data(model,
                                                                     "TaskQuestions:FamiliarityFamiliar:DiagnosisASD  <
                                              TaskQuestions:FamiliarityFamiliar:DiagnosisTDC",
                                                                     "SubjectVariabilityASD_Q_F<TD_Q_F")
  
  SubjectVariabilityASDSmallerThanTD_Q_UF <- extract_relevant_sd_data(model,
                                                                      "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD  <
                                              TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC",
                                                                      "SubjectVariabilityASD_Q_UF<TD_Q_UF")
  
  SubjectVariabilityASD_MG_SmallerThan_Q_F <- extract_relevant_sd_data(model,
                                                                       "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD > TaskQuestions:FamiliarityFamiliar:DiagnosisASD",
                                                                       "SubjectVariabilityASD_MG_SmallerThan_Q_F")
  
  SubjectVariabilityASD_Q_F_SmallerThan_Q_UF <- extract_relevant_sd_data(model,
                                                                         "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD > TaskQuestions:FamiliarityFamiliar:DiagnosisASD",
                                                                         "SubjectVariabilityASD_Q_F_SmallerThan_Q_UF")
  
  SubjectVariabilityTDC_MG_SmallerThan_Q_F <- extract_relevant_sd_data(model,
                                                                       "TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC < TaskQuestions:FamiliarityFamiliar:DiagnosisTDC",
                                                                       "SubjectVariabilityTDC_MG_SmallerThan_Q_F")
  
  SubjectVariabilityTDC_Q_F_SmallerThan_Q_UF <- extract_relevant_sd_data(model,
                                                                         "TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC > TaskQuestions:FamiliarityFamiliar:DiagnosisTDC",
                                                                         "SubjectVariabilityTDC_Q_F_SmallerThan_Q_UF")
  
  DifferenceInDifferenceSubjectVariabilityUnfamiliarSmallerFamiliarity <- extract_relevant_sd_data(model,
                                                                                                   "(TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD - TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC) 
                                              < 
                                              (TaskQuestions:FamiliarityFamiliar:DiagnosisASD - TaskQuestions:FamiliarityFamiliar:DiagnosisTDC)",
                                                                                                   "DifferenceInDifferenceSubjectVariabilityUnfamiliarBiggerFamiliarity")
  
  DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions <- extract_relevant_sd_data(model,
                                                                                                  "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC - TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD) 
                                              < 
                                              (TaskQuestions:FamiliarityFamiliar:DiagnosisTDC - TaskQuestions:FamiliarityFamiliar:DiagnosisASD)",
                                                                                                  "DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions")
  
  # Overall Visit Variability across visits in Autism Group
  VisitVariabilityASDOverall <- extract_relevant_sdvisit_data(model, 
                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
                                                              "VisitVariabilityASDOverall")
  
  # Overall Visit Variability across visits in TD Group
  VisitVariabilityTDCOverall <- extract_relevant_sdvisit_data(model, 
                                                              "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
                                                              "VisitVariabilityTDCOverall")
  
  
  VisitVariabilityASDOverallHigherTDSOverall <- extract_relevant_sdvisit_data(model, 
                                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 >
                                               (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3",
                                                                              "VisitVariabilityASDOverallLowerTDSOverall")
  
  # Subject Variability across visits MG in Autism Group
  VisitVariabilityASDMG <- extract_relevant_sdvisit_data(model, 
                                                         "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                         "VisitVariabilityASDMG")
  
  # Overall Subject Variability across visits in TD Group
  VisitVariabilityTDMG <- extract_relevant_sdvisit_data(model, 
                                                        "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                        "VisitVariabilityTDMG")
  
  # Subject Variability across visits MG in Autism Group
  VisitVariabilityASD_Q_F <- extract_relevant_sdvisit_data(model, 
                                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > 0",
                                                           "VisitVariabilityASD_Q_F")
  
  # Overall Subject Variability across visits in TD Group
  VisitVariabilityTDC_Q_F <- extract_relevant_sdvisit_data(model, 
                                                           "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar > 0",
                                                           "VisitVariabilityTDC_Q_F")
  
  # Subject Variability across visits MG in Autism Group
  VisitVariabilityASD_Q_UF <- extract_relevant_sdvisit_data(model, 
                                                            "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                            "VisitVariabilityASD_Q_UF")
  
  # Overall Subject Variability across visits in TD Group
  VisitVariabilityTDC_Q_UF <- extract_relevant_sdvisit_data(model, 
                                                            "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                            "VisitVariabilityTDC_Q_UF")
  
  # Overall Subject Variability across visits, differences
  VisitVariabilityASDMGSlowerThanTDMG <- extract_relevant_sdvisit_data(model, 
                                                                       "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar < DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar",
                                                                       "VisitVariabilityASDMG>TDMG")
  
  # Overall Subject Variability across visits, differences
  VisitVariabilityASD_Q_F_SlowerThanTD_Q_F <- extract_relevant_sdvisit_data(model, 
                                                                            "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                            "VisitVariabilityASD_Q_F_SlowerThanTD_Q_F")
  
  # Overall Subject Variability across visits, differences
  VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF <- extract_relevant_sdvisit_data(model, 
                                                                              "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                              "VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF")
  
  #DifferenceInDifference:
  DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity <- extract_relevant_sdvisit_data(model, 
                                                                                                       "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)
                                               < (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                                       "DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity")
  
  DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG <- extract_relevant_sdvisit_data(model, 
                                                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar)
                                               > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                              "DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG")
  
  
  # Overall Sigma Variability across visits in Autism Group
  SigmaASDOverall <- extract_relevant_sigma_data(model, 
                                                 "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 < 0",
                                                 "SigmaASDOverall")
  
  # Overall Sigma Variability across visits in TD Group
  SigmaTDCOverall <- extract_relevant_sigma_data(model, 
                                                 "(sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) /3 < 0",
                                                 "SigmaTDCOverall")
  
  
  SigmaASDOverallHigherTDSOverall <- extract_relevant_nonlatency_data(model, 
                                                                      "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 >
                                               (sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3",
                                                                      "SigmaASDOverallHigherTDSOverall")
  
  # Subject Variability across visits MG in Autism Group
  SigmaASDMG <- extract_relevant_sigma_data(model, 
                                            "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar < 0",
                                            "SigmaASDMG")
  
  # Overall Subject Variability across visits in TD Group
  SigmaTDMG <- extract_relevant_sigma_data(model, 
                                           "sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar < 0",
                                           "SigmaTDMG")
  
  # Subject Variability across visits MG in Autism Group
  SigmaASD_Q_F <- extract_relevant_sigma_data(model, 
                                              "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < 0",
                                              "SigmaASD_Q_F")
  
  # Overall Subject Variability across visits in TD Group
  SigmaTDC_Q_F <- extract_relevant_sigma_data(model, 
                                              "sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar < 0",
                                              "SigmaTDC_Q_F")
  
  
  # Subject Variability across visits MG in Autism Group
  SigmaASD_Q_UF <- extract_relevant_sigma_data(model, 
                                               "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < 0 ",
                                               "SigmaASD_Q_UF")
  
  # Overall Subject Variability across visits in TD Group
  SigmaTDC_Q_UF <- extract_relevant_sigma_data(model, 
                                               "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar < 0",
                                               "SigmaTDC_Q_UF")
  
  # Overall Subject Variability across visits, differences
  SigmaASDMGSlowerThanTDMG <- extract_relevant_nonlatency_data(model, 
                                                               "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar",
                                                               "SigmaASDMGSlowerThanTDMG")
  
  # Overall Subject Variability across visits, differences
  SigmaASD_Q_F_SlowerThanTD_Q_F <- extract_relevant_nonlatency_data(model, 
                                                                    "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                    "SigmaASD_Q_F_SlowerThanTD_Q_F")
  
  # Overall Subject Variability across visits, differences
  SigmaASD_Q_UF_SlowerThanTD_Q_UF <- extract_relevant_nonlatency_data(model, 
                                                                      "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                      "SigmaASD_Q_UF_SlowerThanTD_Q_UF")
  
  # Sigma across visits, differences
  # SigmaASD_MG_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit < 0",
  #                                                "SigmaASD_MG_VisitLowerZero")
  # 
  # SigmaTDC_MG_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit > 0",
  #                                                "SigmaTDC_MG_VisitLowerZero")
  # 
  # SigmaASD_MG_VisitLowerTDC <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit < sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
  #                                                "SigmaASD_MG_VisitLowerTDC")
  # 
  # SigmaASD_Overall_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
  #                                                "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit + 
  #                                                sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit +
  #                                                sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit) / 3 < 0",
  #                                                "SigmaASD_Overall_VisitLowerZero")
  # 
  # SigmaTDC_Overall_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
  #                                                "(sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit + 
  #                                                sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit +
  #                                                sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit) / 3 < 0",
  #                                                "SigmaTDC_Overall_VisitLowerZero")
  # 
  # # Sigma across visits, differences
  # SigmaASD_Q_F_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit > 0",
  #                                                "SigmaASD_Q_F_VisitLowerZero")
  # 
  # SigmaTDC_Q_F_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit < 0",
  #                                                "SigmaTDC_Q_F_VisitLowerZero")
  # 
  # SigmaASD_Q_F_VisitLowerTDC <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit > sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
  #                                                "SigmaASD_Q_F_VisitLowerTDC")
  # 
  # # Sigma across visits, differences
  # SigmaASD_Q_UF_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < 0",
  #                                                "SigmaASD_Q_UF_VisitLowerZero")
  # 
  # SigmaTDC_Q_UF_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > 0",
  #                                                "SigmaTDC_Q_UF_VisitLowerZero")
  # 
  # SigmaASD_Q_UF_VisitLowerTDC <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
  #                                                "SigmaASD_Q_UF_VisitLowerTDC")
  # 
  # SigmaVisitOverall_ASDLowerTDC <- extract_relevant_nonlatency_data(model, 
  #                                                "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit + 
  #                                                sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit +
  #                                                sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit) / 3 > (sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit + 
  #                                                sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit +
  #                                                sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit) /3",
  #                                                "SigmaVisitOverall_ASDLowerTDC")
  
  #DifferenceInDifference:
  DifferenceInDifferenceSigmaASDFamiliarityBiggerTDUnfamiliarity <- extract_relevant_nonlatency_data(model, 
                                                                                                     "(sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)
                                               > (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                                     "DifferenceInDifferenceSigmaASDFamiliarityBiggerTDUnfamiliarity")
  
  DifferenceInDifferenceSigmaASDMGSmallerThanTDMG <- extract_relevant_nonlatency_data(model, 
                                                                                      "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar)
                                               > (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                      "DifferenceInDifferenceSigmaASDMGSmallerThanTDMG")
  
  #WithinGroups
  WithinGroupSigmaASD_Q_F_Faster_UF <- extract_relevant_nonlatency_data(model, 
                                                                        "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaASD_Q_FFaster_UF")
  
  WithinGroupSigmaTDC_Q_F_Faster_UF <- extract_relevant_nonlatency_data(model, 
                                                                        "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaTDC_Q_F_Faster_UF")
  
  WithinGroupSigmaASD_Q_F_Faster_MG <- extract_relevant_nonlatency_data(model, 
                                                                        "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaASD_Q_F_Faster_MG")
  
  WithinGroupSigmaTDC_Q_F_Faster_MG <- extract_relevant_nonlatency_data(model, 
                                                                        "sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar < sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaTDC_Q_F_Faster_MG")
  
  
  # WithinGroupSigmaASD_Q_F_Lower_UF_Visit <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit",
  #                                                "WithinGroupSigmaASD_Q_F_Lower_UF_Visit")
  # 
  # WithinGroupSigmaTDC_Q_F_Lower_UF_Visit <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit",
  #                                                "WithinGroupSigmaTDC_Q_F_Lower_UF_Visit")
  # 
  # WithinGroupSigmaASD_Q_F_Lower_MG_Visit <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit > sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit",
  #                                                "WithinGroupSigmaASD_Q_F_Lower_MG_Visit")
  # 
  # WithinGroupSigmaTDC_Q_F_Lower_MG_Visit <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
  #                                                "WithinGroupSigmaTDC_Q_F_Lower_MG_Visit")
  # 
  # 
  # SigmaTDC_Q_UF_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > 0",
  #                                                "SigmaTDC_Q_UF_VisitLowerZero")
  # 
  # SigmaASD_Q_UF_VisitLowerTDC <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
  #                                                "SigmaASD_Q_UF_VisitLowerTDC")
  # 
  # 
  # DifferenceInDifferenceSigmaVisit_Q_F_Lower_Q_UF <- extract_relevant_nonlatency_data(model, 
  #                                                "(sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit - sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit) <
  #                                                (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit)", "DifferenceInDifferenceSigmaVisit_Q_F_Lower_Q_UF")
  # 
  # 
  # DifferenceInDifferenceSigmaVisit_MG_Lower_Q_F <- extract_relevant_nonlatency_data(model, 
  #                                                "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit - sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit) <
  #                                                (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit)", "DifferenceInDifferenceSigmaVisit_MG_Lower_Q_F")
  
  HypothesisResultsAdultlatenciesModel1 <- rbind(OverallLatencyThreeConditionsASD,
                                                 OverallLatencyThreeConditionsTD,
                                                 MGOverallLatencyASD,
                                                 Q_F_OverallLatencyASD,
                                                 Q_UF_OverallLatencyASD,
                                                 MGOverallLatencyTDC,
                                                 Q_F_OverallLatencyTDC,
                                                 Q_UF_OverallLatencyTDC,
                                                 BetaParameterASD_MG,
                                                 BetaParameterASD_Q_F,
                                                 BetaParameterASD_Q_UF,
                                                 BetaParameterASD_Overall,
                                                 BetaParameterTDC_MG,
                                                 BetaParameterTDC_Q_F,
                                                 BetaParameterTDC_Q_UF,
                                                 BetaParameterTDC_Overall,
                                                 OverallBetaParameterASDLowerTD,
                                                 BetaParameterASDFasterTD_MG,
                                                 BetaParameterASDFasterTD_Q_F,
                                                 BetaParameterASDFasterTD_Q_UF,
                                                 WithinGroupBetaParameterASD_MG_Slower_Q_F,
                                                 WithinGroupBetaParameterTDC_MG_Slower_Q_F,
                                                 DifferenceInDifferenceBetaParameterBetweenGroupsMGBiggerQ_F,
                                                 WithinGroupBetaParameterASD_Q_F_Slower_Q_UF,
                                                 WithinGroupBetaParameterTDC_Q_F_Slower_Q_UF,
                                                 DifferenceInDifferenceBetaParameterBetweenGroupsFamiliarityBiggerUnfamiliarity,
                                                 OverallLatencyThreeConditionsASDFasterTD,
                                                 WithinASDGroupMGSlowerQuestionsFamiliar,
                                                 WithinTDGroupMGSlowerQuestionsFamiliar,
                                                 WithinAutismGroupQuestionsUnfamiliarSlowerFamiliar,
                                                 WithinTDGroupQuestionsUnfamiliarSlowerFamiliar,
                                                 WithinConditionsMatchingGameASDFasterTD,
                                                 WithinConditionsQuestionsFamiliarASDFasterTD,
                                                 WithinConditionsQuestionsUnfamiliarASDFasterTD,
                                                 DifferenceInDifferenceBetweenGroupsUnfamiliarBiggerFamiliarity,
                                                 DifferenceInDifferenceBetweenGroupsMGBiggerThanConversations,
                                                 SubjectVariabilityASD,
                                                 SubjectVariabilityTD,
                                                 SubjectVariabilityASDSmallerThanTD,
                                                 SubjectVariabilityASD_MG,
                                                 SubjectVariabilityASD_Q_F,
                                                 SubjectVariabilityASD_Q_UF,
                                                 SubjectVariabilityTD_MG,
                                                 SubjectVariabilityTD_Q_F,
                                                 SubjectVariabilityTD_Q_UF,
                                                 SubjectVariabilityASDSmallerThanTD_MG,
                                                 SubjectVariabilityASDSmallerThanTD_Q_F,
                                                 SubjectVariabilityASDSmallerThanTD_Q_UF,
                                                 SubjectVariabilityASD_MG_SmallerThan_Q_F,
                                                 SubjectVariabilityASD_Q_F_SmallerThan_Q_UF,
                                                 SubjectVariabilityTDC_MG_SmallerThan_Q_F,
                                                 SubjectVariabilityTDC_Q_F_SmallerThan_Q_UF,
                                                 DifferenceInDifferenceSubjectVariabilityUnfamiliarSmallerFamiliarity,
                                                 DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions,
                                                 VisitVariabilityASDMG,
                                                 VisitVariabilityTDMG,
                                                 VisitVariabilityASD_Q_F,
                                                 VisitVariabilityTDC_Q_F,
                                                 VisitVariabilityASD_Q_UF,
                                                 VisitVariabilityTDC_Q_UF,
                                                 VisitVariabilityASDOverall,
                                                 VisitVariabilityTDCOverall,
                                                 VisitVariabilityASDMGSlowerThanTDMG,
                                                 VisitVariabilityASD_Q_F_SlowerThanTD_Q_F,
                                                 VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF,
                                                 VisitVariabilityASDOverallHigherTDSOverall,
                                                 DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG,
                                                 DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity,
                                                 SigmaASDOverall,
                                                 SigmaTDCOverall,
                                                 SigmaASDOverallHigherTDSOverall,
                                                 SigmaASDMG,
                                                 SigmaTDMG,
                                                 SigmaASD_Q_F,
                                                 SigmaTDC_Q_F,
                                                 SigmaASD_Q_UF,
                                                 SigmaTDC_Q_UF,
                                                 SigmaASDMGSlowerThanTDMG,
                                                 SigmaASD_Q_F_SlowerThanTD_Q_F,
                                                 SigmaASD_Q_UF_SlowerThanTD_Q_UF,
                                                 # SigmaASD_MG_VisitLowerZero,
                                                 # SigmaTDC_MG_VisitLowerZero,
                                                 # SigmaASD_MG_VisitLowerTDC,
                                                 # SigmaASD_Q_F_VisitLowerZero,
                                                 # SigmaTDC_Q_F_VisitLowerZero,
                                                 # SigmaASD_Q_F_VisitLowerTDC,
                                                 # SigmaASD_Q_UF_VisitLowerZero,
                                                 # SigmaTDC_Q_UF_VisitLowerZero,
                                                 # SigmaASD_Q_UF_VisitLowerTDC,
                                                 # SigmaASD_Overall_VisitLowerZero,
                                                 # SigmaTDC_Overall_VisitLowerZero,
                                                 # WithinGroupSigmaASD_Q_F_Lower_UF_Visit,
                                                 # WithinGroupSigmaTDC_Q_F_Lower_UF_Visit,
                                                 # WithinGroupSigmaASD_Q_F_Lower_MG_Visit,
                                                 # WithinGroupSigmaTDC_Q_F_Lower_MG_Visit,
                                                 WithinGroupSigmaASD_Q_F_Faster_UF,
                                                 WithinGroupSigmaTDC_Q_F_Faster_UF,
                                                 WithinGroupSigmaASD_Q_F_Faster_MG,
                                                 WithinGroupSigmaTDC_Q_F_Faster_MG,
                                                 #SigmaVisitOverall_ASDLowerTDC,
                                                 # DifferenceInDifferenceSigmaVisit_Q_F_Lower_Q_UF,
                                                 # DifferenceInDifferenceSigmaVisit_MG_Lower_Q_F,
                                                 DifferenceInDifferenceSigmaASDFamiliarityBiggerTDUnfamiliarity,
                                                 DifferenceInDifferenceSigmaASDMGSmallerThanTDMG,
                                                 CorrelationSummaryStatistics) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  
  return(HypothesisResultsAdultlatenciesModel1)
}

extract_all_hypothesis_tests_m2_child <- function(model) {
  LanguageASD_TG_F <- extract_relevant_data(model,
                                            "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled > 0",
                                            "LanguageASD_TG_F")
  
  LanguageASD_Q_F <- extract_relevant_data(model,
                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled > 0",
                                           "LanguageASD_Q_F")
  
  LanguageASD_Q_UF <- extract_relevant_data(model,
                                            "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > 0",
                                            "LanguageASD_Q_UF")
  
  LanguageASD_overall <- extract_relevant_data(model,
                                               "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)/3 > 0",
                                               "LanguageASD_overall")
  
  LanguageTDC_TG_F <- extract_relevant_data(model,
                                            "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled > 0",
                                            "LanguageTDC_TG_F")
  
  LanguageTDC_Q_F <- extract_relevant_data(model,
                                           "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled > 0",
                                           "LanguageTDC_Q_F")
  
  LanguageTDC_Q_UF <- extract_relevant_data(model,
                                            "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > 0",
                                            "LanguageTDC_Q_UF")
  
  LanguageTDC_overall <- extract_relevant_data(model,
                                               "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)/3 > 0",
                                               "LanguageTDC_overall")
  
  MotorASD_TG_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled < 0", "MotorASD_TG_F")
  MotorASD_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled < 0", "MotorASD_Q_F")
  MotorASD_Q_UF <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled < 0", "MotorASD_Q_UF")
  
  #1
  WithinGroupMotorASD_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled", "WithinGroupMotorASD_MGSlower_Q_F")
  
  WithinGroupMotorASD_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled", "WithinGroupMotorASD_Q_UFSlower_Q_F")
  #2
  WithinGroupCognitionASD_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled < DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled", "WithinGroupCognitionASD_MGSlower_Q_F")
  
  WithinGroupCognitionASD_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled", "WithinGroupCognitionASD_Q_UFSlower_Q_F")
  #3
  WithinGroupAwarenessASD_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled < DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled", "WithinGroupAwarenessASD_MGSlower_Q_F")
  
  WithinGroupAwarenessASD_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled < DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled", "WithinGroupAwarenessASD_Q_UFSlower_Q_F")
  #4
  WithinGroupMotivationASD_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled", "WithinGroupMotivationASD_MGSlower_Q_F")
  
  WithinGroupMotivationASD_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled", "WithinGroupMotivationASD_Q_UFSlower_Q_F")
  #5
  WithinGroupLanguageASD_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled", "WithinGroupLanguageASD_MGSlower_Q_F")
  
  WithinGroupLanguageASD_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled", "WithinGroupLanguageASD_Q_UFSlower_Q_F")
  
  #1
  WithinGroupMotorTDC_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled", "WithinGroupMotorTDC_MGSlower_Q_F")
  
  WithinGroupMotorTDC_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled", "WithinGroupMotorTDC_Q_UFSlower_Q_F")
  #2
  WithinGroupCognitionTDC_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled", "WithinGroupCognitionTDC_MGSlower_Q_F")
  
  WithinGroupCognitionTDC_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled", "WithinGroupCognitionTDC_Q_UFSlower_Q_F")
  #3
  WithinGroupAwarenessTDC_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled", "WithinGroupAwarenessTDC_MGSlower_Q_F")
  
  WithinGroupAwarenessTDC_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled", "WithinGroupAwarenessTDC_Q_UFSlower_Q_F")
  #4
  WithinGroupMotivationTDC_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled", "WithinGroupMotivationTDC_MGSlower_Q_F")
  
  WithinGroupMotivationTDC_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled", "WithinGroupMotivationTDC_Q_UFSlower_Q_F")
  #5
  WithinGroupLanguageTDC_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled", "WithinGroupLanguageTDC_MGSlower_Q_F")
  
  WithinGroupLanguageTDC_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled", "WithinGroupLanguageTDC_Q_UFSlower_Q_F")
  
  MotorASD_overall <- extract_relevant_data(model, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled)/3 < 0",  "MotorASD_overall")
  
  MotorTDC_TG_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled < 0", "MotorTDC_TG_F")
  MotorTDC_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled < 0", "MotorTDC_Q_F")
  MotorTDC_Q_UF <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled < 0", "MotorTDC_Q_UF")
  MotorTDC_overall <- extract_relevant_data(model, "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled)/3 < 0", "MotorTDC_overall")
  
  CognitionASD_TG_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled < 0", "CognitionASD_TG_F")
  CognitionASD_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled < 0", "CognitionASD_Q_F")
  CognitionASD_Q_UF <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled < 0","CognitionASD_Q_UF")
  CognitionASD_overall <- extract_relevant_data(model, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled)/3 < 0","CognitionASD_overall")
  
  CognitionTDC_TG_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled < 0","CognitionTDC_TG_F")
  CognitionTDC_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled < 0","CognitionTDC_Q_F")
  CognitionTDC_Q_UF <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled < 0","CognitionTDC_Q_UF")
  CognitionTDC_overall <- extract_relevant_data(model, "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled)/3 < 0","CognitionTDC_overall")
  
  AwarenessASD_TG_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled > 0","AwarenessASD_TG_F")
  AwarenessASD_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled > 0","AwarenessASD_Q_F")
  AwarenessASD_Q_UF <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled < 0", "AwarenessASD_Q_UF")
  AwarenessASD_overall <- extract_relevant_data(model, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled)/3 > 0", "AwarenessASD_overall")
  
  AwarenessTDC_TG_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled > 0", "AwarenessTDC_TG_F")
  AwarenessTDC_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled > 0","AwarenessTDC_Q_F")
  AwarenessTDC_Q_UF <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled > 0","AwarenessTDC_Q_UF")
  AwarenessTDC_overall <- extract_relevant_data(model, "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled)/3 > 0","AwarenessTDC_overall")
  
  MotivationASD_TG_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled > 0","MotivationASD_TG_F")
  MotivationASD_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled > 0","MotivationASD_Q_F")
  MotivationASD_Q_UF <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled > 0","MotivationASD_Q_UF")
  MotivationASD_overall <-extract_relevant_data(model, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled)/3 > 0","MotivationASD_overall")
  
  MotivationTDC_TG_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled < 0","MotivationTDC_TG_F")
  MotivationTDC_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled > 0","MotivationTDC_Q_F")
  MotivationTDC_Q_UF <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled > 0","MotivationTDC_Q_UF")
  MotivationTDC_overall <- extract_relevant_data(model, "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled)/3 > 0","MotivationTDC_overall")
  
  
  LanguageASD_TG_F <- extract_relevant_data(model,
                                            "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled > 0",
                                            "LanguageASD_TG_F")
  
  LanguageASD_Q_F <- extract_relevant_data(model,
                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled > 0",
                                           "LanguageASD_Q_F")
  
  LanguageASD_Q_UF <- extract_relevant_data(model,
                                            "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > 0",
                                            "LanguageASD_Q_UF")
  
  LanguageASD_overall <- extract_relevant_data(model,
                                               "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)/3 > 0",
                                               "LanguageASD_overall")
  
  LanguageTDC_TG_F <- extract_relevant_data(model,
                                            "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled < 0",
                                            "LanguageTDC_TG_F")
  
  LanguageTDC_Q_F <- extract_relevant_data(model,
                                           "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled > 0",
                                           "LanguageTDC_Q_F")
  
  LanguageTDC_Q_UF <- extract_relevant_data(model,
                                            "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > 0",
                                            "LanguageTDC_Q_UF")
  
  LanguageTDC_overall <- extract_relevant_data(model,
                                               "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)/3 > 0",
                                               "LanguageTDC_overall")
  
  
  LanguageOverallASDSlowerTDC <- extract_relevant_data(model,
                                                       "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)/3 > (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)/3",
                                                       "LanguageOverallASDSlowerTDC")
  
  MotorOverallASDSlowerTDC <- extract_relevant_data(model, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled)/3 > 
                                          (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled)/3",  
                                                    "MotorOverallASDFasterTDC")
  
  CognitionOverallASDSlowerTDC <- extract_relevant_data(model, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled)/3 >
                                              (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled)/3",
                                                        "CognitionOverallASDSlowerTDC")
  
  AwarenessOverallASDFasterTDC <- extract_relevant_data(model, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled)/3 < (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled)/3", "AwarenessOverallASDSlowerTDC")
  
  MotivationOverallASDFasterTDC <-extract_relevant_data(model, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled)/3 > (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled)/3","MotivationOverallASDFasterTDC")
  
  LanguageOverallASDSlowerTDC_MG <- extract_relevant_data(model,
                                                          "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled > DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled",
                                                          "LanguageOverallASDSlowerTDC_MG")
  
  LanguageOverallASDSlowerTDC_Q_F <- extract_relevant_data(model,
                                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled",
                                                           "LanguageOverallASDSlowerTDC_Q_F")
  
  LanguageOverallASDSlowerTDC_Q_UF <- extract_relevant_data(model,
                                                            "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled",
                                                            "LanguageOverallASDSlowerTDC_Q_UF")
  
  
  MotorOverallASDSlowerTDC_MG <- extract_relevant_data(model,
                                                       "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled > DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled",
                                                       "MotorOverallASDSlowerTDC_MG")
  
  MotorOverallASDSlowerTDC_Q_F <- extract_relevant_data(model,
                                                        "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled",
                                                        "MotorOverallASDSlowerTDC_Q_F")
  
  MotorOverallASDSlowerTDC_Q_UF <- extract_relevant_data(model,
                                                         "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled",
                                                         "MotorOverallASDSlowerTDC_Q_UF")
  
  CognitionOverallASDSlowerTDC_MG <- extract_relevant_data(model,
                                                           "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled > DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled",
                                                           "CognitionOverallASDSlowerTDC_MG")
  
  CognitionOverallASDSlowerTDC_Q_F <- extract_relevant_data(model,
                                                            "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled",
                                                            "CognitionOverallASDSlowerTDC_Q_F")
  
  CognitionOverallASDSlowerTDC_Q_UF <- extract_relevant_data(model,
                                                             "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled",
                                                             "CognitionOverallASDSlowerTDC_Q_UF")
  
  AwarenessOverallASDSlowerTDC_MG <- extract_relevant_data(model,
                                                           "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled < DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled",
                                                           "AwarenessOverallASDSlowerTDC_MG")
  
  AwarenessOverallASDSlowerTDC_Q_F <- extract_relevant_data(model,
                                                            "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled",
                                                            "AwarenessOverallASDSlowerTDC_Q_F")
  
  AwarenessOverallASDSlowerTDC_Q_UF <- extract_relevant_data(model,
                                                             "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled",
                                                             "AwarenessOverallASDSlowerTDC_Q_UF")
  
  MotivationOverallASDSlowerTDC_MG <- extract_relevant_data(model,
                                                            "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled > DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled",
                                                            "MotivationOverallASDSlowerTDC_MG")
  
  MotivationOverallASDSlowerTDC_Q_F <- extract_relevant_data(model,
                                                             "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled",
                                                             "MotivationOverallASDSlowerTDC_Q_F")
  
  MotivationOverallASDSlowerTDC_Q_UF <- extract_relevant_data(model,
                                                              "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled",
                                                              "MotivationOverallASDSlowerTDC_Q_UF")
  
  
  DifferenceInDifferenceLanguageMGSlowerQ_F <- extract_relevant_data(model,
                                                                     "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled)",
                                                                     "DifferenceInDifferenceLanguageMGSlowerQ_F")
  
  DifferenceInDifferenceLanguageUFSlowerQ_F <- extract_relevant_data(model,
                                                                     "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled)",
                                                                     "DifferenceInDifferenceLanguageUFSlowerQ_F")
  
  DifferenceInDifferenceMotorMGSlowerQ_F <- extract_relevant_data(model,
                                                                  "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled)",
                                                                  "DifferenceInDifferenceMotorMGSlowerQ_F")
  
  DifferenceInDifferenceMotorUFSlowerQ_F <- extract_relevant_data(model,
                                                                  "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled)",
                                                                  "DifferenceInDifferenceMotorUFSlowerQ_F")
  
  DifferenceInDifferenceCognitionMGSlowerQ_F <- extract_relevant_data(model,
                                                                      "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled)",
                                                                      "DifferenceInDifferenceCognitionMGSlowerQ_F")
  
  DifferenceInDifferenceCognitionUFSlowerQ_F <- extract_relevant_data(model,
                                                                      "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled)",
                                                                      "DifferenceInDifferenceCognitionUFSlowerQ_F")
  
  DifferenceInDifferenceAwarenessMGSlowerQ_F <- extract_relevant_data(model,
                                                                      "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled)
                      < (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled)",
                                                                      "DifferenceInDifferenceAwarenessMGSlowerQ_F")
  
  DifferenceInDifferenceAwarenessUFSlowerQ_F <- extract_relevant_data(model,
                                                                      "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled)
                      < (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled)",
                                                                      "DifferenceInDifferenceAwarenessUFSlowerQ_F")
  
  DifferenceInDifferenceMotivationMGSlowerQ_F <- extract_relevant_data(model,
                                                                       "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled)",
                                                                       "DifferenceInDifferenceMotivationMGSlowerQ_F")
  
  DifferenceInDifferenceMotivationUFSlowerQ_F <- extract_relevant_data(model,
                                                                       "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled)
                      < (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled)",
                                                                       "DifferenceInDifferenceMotivationUFSlowerQ_F")
  
  LanguageOverallASDSlowerTDC_Q_F <- extract_relevant_data(model,
                                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled",
                                                           "LanguageOverallASDSlowerTDC_Q_F")
  
  LanguageOverallASDSlowerTDC_Q_UF <- extract_relevant_data(model,
                                                            "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled",
                                                            "LanguageOverallASDSlowerTDC_Q_UF")
  
  hypothesis_data_full <- rbind(
    MotorASD_TG_F, MotorASD_Q_F, MotorASD_Q_UF, MotorASD_overall,
    MotorTDC_TG_F, MotorTDC_Q_F, MotorTDC_Q_UF, MotorTDC_overall,
    
    CognitionASD_TG_F, CognitionASD_Q_F, CognitionASD_Q_UF, CognitionASD_overall,
    CognitionTDC_TG_F, CognitionTDC_Q_F, CognitionTDC_Q_UF, CognitionTDC_overall,
    
    LanguageASD_TG_F, LanguageASD_Q_F, LanguageASD_Q_UF, LanguageASD_overall,
    LanguageTDC_TG_F, LanguageTDC_Q_F, LanguageTDC_Q_UF, LanguageTDC_overall,
    
    AwarenessASD_TG_F, AwarenessASD_Q_F, AwarenessASD_Q_UF, AwarenessASD_overall,
    AwarenessTDC_TG_F, AwarenessTDC_Q_F, AwarenessTDC_Q_UF, AwarenessTDC_overall,
    
    MotivationASD_TG_F, MotivationASD_Q_F, MotivationASD_Q_UF, MotivationASD_overall,
    MotivationTDC_TG_F, MotivationTDC_Q_F, MotivationTDC_Q_UF, MotivationTDC_overall,
    
    LanguageOverallASDSlowerTDC, MotorOverallASDSlowerTDC, CognitionOverallASDSlowerTDC, AwarenessOverallASDFasterTDC, MotivationOverallASDFasterTDC,
    
    LanguageOverallASDSlowerTDC_MG, LanguageOverallASDSlowerTDC_Q_F, LanguageOverallASDSlowerTDC_Q_UF,
    CognitionOverallASDSlowerTDC_MG, CognitionOverallASDSlowerTDC_Q_F, CognitionOverallASDSlowerTDC_Q_UF,
    MotorOverallASDSlowerTDC_MG, MotorOverallASDSlowerTDC_Q_F, MotorOverallASDSlowerTDC_Q_UF,
    AwarenessOverallASDSlowerTDC_MG, AwarenessOverallASDSlowerTDC_Q_F, AwarenessOverallASDSlowerTDC_Q_UF,
    MotivationOverallASDSlowerTDC_MG, MotivationOverallASDSlowerTDC_Q_F, MotivationOverallASDSlowerTDC_Q_UF,
    
    DifferenceInDifferenceLanguageMGSlowerQ_F, DifferenceInDifferenceLanguageUFSlowerQ_F,
    DifferenceInDifferenceMotorMGSlowerQ_F, DifferenceInDifferenceMotorUFSlowerQ_F,
    DifferenceInDifferenceMotivationMGSlowerQ_F, DifferenceInDifferenceMotivationUFSlowerQ_F,
    DifferenceInDifferenceCognitionMGSlowerQ_F, DifferenceInDifferenceCognitionUFSlowerQ_F,
    DifferenceInDifferenceAwarenessMGSlowerQ_F, DifferenceInDifferenceAwarenessUFSlowerQ_F,
    
    WithinGroupLanguageASD_MGSlower_Q_F, WithinGroupLanguageASD_Q_UFSlower_Q_F,
    WithinGroupCognitionASD_MGSlower_Q_F, WithinGroupCognitionASD_Q_UFSlower_Q_F,
    WithinGroupMotorASD_MGSlower_Q_F, WithinGroupMotorASD_Q_UFSlower_Q_F,
    WithinGroupAwarenessASD_MGSlower_Q_F, WithinGroupAwarenessASD_Q_UFSlower_Q_F,
    WithinGroupMotivationASD_MGSlower_Q_F, WithinGroupMotivationASD_Q_UFSlower_Q_F,
    
    WithinGroupLanguageTDC_MGSlower_Q_F, WithinGroupLanguageTDC_Q_UFSlower_Q_F,
    WithinGroupCognitionTDC_MGSlower_Q_F, WithinGroupCognitionTDC_Q_UFSlower_Q_F,
    WithinGroupMotorTDC_MGSlower_Q_F, WithinGroupMotorTDC_Q_UFSlower_Q_F,
    WithinGroupAwarenessTDC_MGSlower_Q_F, WithinGroupAwarenessTDC_Q_UFSlower_Q_F,
    WithinGroupMotivationTDC_MGSlower_Q_F, WithinGroupMotivationTDC_Q_UFSlower_Q_F
    
    
  ) %>%
    rename("value" = TestDescription) %>%
    #mutate(Estimate = round(Estimate, 1)) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 1 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -1 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 1 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 1 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -1 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -1 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  
  return(hypothesis_data_full)
}

extract_all_hypothesis_tests_m2_adult <- function(model) {
  
  LanguageASD_TG_F <- extract_relevant_data(model,
                                            "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled > 0",
                                            "LanguageASD_TG_F")
  
  LanguageASD_Q_F <- extract_relevant_data(model,
                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled > 0",
                                           "LanguageASD_Q_F")
  
  LanguageASD_Q_UF <- extract_relevant_data(model,
                                            "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > 0",
                                            "LanguageASD_Q_UF")
  
  LanguageASD_overall <- extract_relevant_data(model,
                                               "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)/3 > 0",
                                               "LanguageASD_overall")
  
  LanguageTDC_TG_F <- extract_relevant_data(model,
                                            "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled > 0",
                                            "LanguageTDC_TG_F")
  
  LanguageTDC_Q_F <- extract_relevant_data(model,
                                           "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled > 0",
                                           "LanguageTDC_Q_F")
  
  LanguageTDC_Q_UF <- extract_relevant_data(model,
                                            "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > 0",
                                            "LanguageTDC_Q_UF")
  
  LanguageTDC_overall <- extract_relevant_data(model,
                                               "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)/3 > 0",
                                               "LanguageTDC_overall")
  
  MotorASD_TG_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled < 0", "MotorASD_TG_F")
  MotorASD_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled < 0", "MotorASD_Q_F")
  MotorASD_Q_UF <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled < 0", "MotorASD_Q_UF")
  
  #1
  WithinGroupMotorASD_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled", "WithinGroupMotorASD_MGSlower_Q_F")
  
  WithinGroupMotorASD_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled", "WithinGroupMotorASD_Q_UFSlower_Q_F")
  #2
  WithinGroupCognitionASD_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled < DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled", "WithinGroupCognitionASD_MGSlower_Q_F")
  
  WithinGroupCognitionASD_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled", "WithinGroupCognitionASD_Q_UFSlower_Q_F")
  #3
  WithinGroupAwarenessASD_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled < DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled", "WithinGroupAwarenessASD_MGSlower_Q_F")
  
  WithinGroupAwarenessASD_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled < DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled", "WithinGroupAwarenessASD_Q_UFSlower_Q_F")
  #4
  WithinGroupMotivationASD_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled", "WithinGroupMotivationASD_MGSlower_Q_F")
  
  WithinGroupMotivationASD_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled", "WithinGroupMotivationASD_Q_UFSlower_Q_F")
  #5
  WithinGroupLanguageASD_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled", "WithinGroupLanguageASD_MGSlower_Q_F")
  
  WithinGroupLanguageASD_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled", "WithinGroupLanguageASD_Q_UFSlower_Q_F")
  
  #1
  WithinGroupMotorTDC_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled", "WithinGroupMotorTDC_MGSlower_Q_F")
  
  WithinGroupMotorTDC_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled", "WithinGroupMotorTDC_Q_UFSlower_Q_F")
  #2
  WithinGroupCognitionTDC_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled", "WithinGroupCognitionTDC_MGSlower_Q_F")
  
  WithinGroupCognitionTDC_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled", "WithinGroupCognitionTDC_Q_UFSlower_Q_F")
  #3
  WithinGroupAwarenessTDC_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled", "WithinGroupAwarenessTDC_MGSlower_Q_F")
  
  WithinGroupAwarenessTDC_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled", "WithinGroupAwarenessTDC_Q_UFSlower_Q_F")
  #4
  WithinGroupMotivationTDC_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled", "WithinGroupMotivationTDC_MGSlower_Q_F")
  
  WithinGroupMotivationTDC_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled", "WithinGroupMotivationTDC_Q_UFSlower_Q_F")
  #5
  WithinGroupLanguageTDC_MGSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled", "WithinGroupLanguageTDC_MGSlower_Q_F")
  
  WithinGroupLanguageTDC_Q_UFSlower_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled", "WithinGroupLanguageTDC_Q_UFSlower_Q_F")
  
  MotorASD_overall <- extract_relevant_data(model, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled)/3 < 0",  "MotorASD_overall")
  
  MotorTDC_TG_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled < 0", "MotorTDC_TG_F")
  MotorTDC_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled < 0", "MotorTDC_Q_F")
  MotorTDC_Q_UF <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled < 0", "MotorTDC_Q_UF")
  MotorTDC_overall <- extract_relevant_data(model, "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled)/3 < 0", "MotorTDC_overall")
  
  CognitionASD_TG_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled < 0", "CognitionASD_TG_F")
  CognitionASD_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled < 0", "CognitionASD_Q_F")
  CognitionASD_Q_UF <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled < 0","CognitionASD_Q_UF")
  CognitionASD_overall <- extract_relevant_data(model, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled)/3 < 0","CognitionASD_overall")
  
  CognitionTDC_TG_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled < 0","CognitionTDC_TG_F")
  CognitionTDC_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled < 0","CognitionTDC_Q_F")
  CognitionTDC_Q_UF <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled < 0","CognitionTDC_Q_UF")
  CognitionTDC_overall <- extract_relevant_data(model, "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled)/3 < 0","CognitionTDC_overall")
  
  AwarenessASD_TG_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled > 0","AwarenessASD_TG_F")
  AwarenessASD_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled > 0","AwarenessASD_Q_F")
  AwarenessASD_Q_UF <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled > 0", "AwarenessASD_Q_UF")
  AwarenessASD_overall <- extract_relevant_data(model, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled)/3 > 0", "AwarenessASD_overall")
  
  AwarenessTDC_TG_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled > 0", "AwarenessTDC_TG_F")
  AwarenessTDC_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled > 0","AwarenessTDC_Q_F")
  AwarenessTDC_Q_UF <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled > 0","AwarenessTDC_Q_UF")
  AwarenessTDC_overall <- extract_relevant_data(model, "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled)/3 > 0","AwarenessTDC_overall")
  
  MotivationASD_TG_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled > 0","MotivationASD_TG_F")
  MotivationASD_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled > 0","MotivationASD_Q_F")
  MotivationASD_Q_UF <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled > 0","MotivationASD_Q_UF")
  MotivationASD_overall <-extract_relevant_data(model, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled)/3 > 0","MotivationASD_overall")
  
  MotivationTDC_TG_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled < 0","MotivationTDC_TG_F")
  MotivationTDC_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled < 0","MotivationTDC_Q_F")
  MotivationTDC_Q_UF <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled < 0","MotivationTDC_Q_UF")
  MotivationTDC_overall <- extract_relevant_data(model, "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled)/3 > 0","MotivationTDC_overall")
  
  
  LanguageASD_TG_F <- extract_relevant_data(model,
                                            "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled > 0",
                                            "LanguageASD_TG_F")
  
  LanguageASD_Q_F <- extract_relevant_data(model,
                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled > 0",
                                           "LanguageASD_Q_F")
  
  LanguageASD_Q_UF <- extract_relevant_data(model,
                                            "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > 0",
                                            "LanguageASD_Q_UF")
  
  LanguageASD_overall <- extract_relevant_data(model,
                                               "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)/3 > 0",
                                               "LanguageASD_overall")
  
  LanguageTDC_TG_F <- extract_relevant_data(model,
                                            "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled > 0",
                                            "LanguageTDC_TG_F")
  
  LanguageTDC_Q_F <- extract_relevant_data(model,
                                           "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled > 0",
                                           "LanguageTDC_Q_F")
  
  LanguageTDC_Q_UF <- extract_relevant_data(model,
                                            "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > 0",
                                            "LanguageTDC_Q_UF")
  
  LanguageTDC_overall <- extract_relevant_data(model,
                                               "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)/3 > 0",
                                               "LanguageTDC_overall")
  
  
  LanguageOverallASDSlowerTDC <- extract_relevant_data(model,
                                                       "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)/3 > (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)/3",
                                                       "LanguageOverallASDSlowerTDC")
  
  MotorOverallASDSlowerTDC <- extract_relevant_data(model, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled)/3 > 
                                          (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled)/3",  
                                                    "MotorOverallASDFasterTDC")
  
  CognitionOverallASDSlowerTDC <- extract_relevant_data(model, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled)/3 >
                                              (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled)/3",
                                                        "CognitionOverallASDSlowerTDC")
  
  AwarenessOverallASDFasterTDC <- extract_relevant_data(model, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled)/3 < (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled)/3", "AwarenessOverallASDSlowerTDC")
  
  MotivationOverallASDFasterTDC <-extract_relevant_data(model, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled)/3 > (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled)/3","MotivationOverallASDFasterTDC")
  
  LanguageOverallASDSlowerTDC_MG <- extract_relevant_data(model,
                                                          "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled > DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled",
                                                          "LanguageOverallASDSlowerTDC_MG")
  
  LanguageOverallASDSlowerTDC_Q_F <- extract_relevant_data(model,
                                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled",
                                                           "LanguageOverallASDSlowerTDC_Q_F")
  
  LanguageOverallASDSlowerTDC_Q_UF <- extract_relevant_data(model,
                                                            "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled",
                                                            "LanguageOverallASDSlowerTDC_Q_UF")
  
  
  MotorOverallASDSlowerTDC_MG <- extract_relevant_data(model,
                                                       "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled > DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled",
                                                       "MotorOverallASDSlowerTDC_MG")
  
  MotorOverallASDSlowerTDC_Q_F <- extract_relevant_data(model,
                                                        "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled",
                                                        "MotorOverallASDSlowerTDC_Q_F")
  
  MotorOverallASDSlowerTDC_Q_UF <- extract_relevant_data(model,
                                                         "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled",
                                                         "MotorOverallASDSlowerTDC_Q_UF")
  
  CognitionOverallASDSlowerTDC_MG <- extract_relevant_data(model,
                                                           "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled > DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled",
                                                           "CognitionOverallASDSlowerTDC_MG")
  
  CognitionOverallASDSlowerTDC_Q_F <- extract_relevant_data(model,
                                                            "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled",
                                                            "CognitionOverallASDSlowerTDC_Q_F")
  
  CognitionOverallASDSlowerTDC_Q_UF <- extract_relevant_data(model,
                                                             "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled",
                                                             "CognitionOverallASDSlowerTDC_Q_UF")
  
  AwarenessOverallASDSlowerTDC_MG <- extract_relevant_data(model,
                                                           "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled < DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled",
                                                           "AwarenessOverallASDSlowerTDC_MG")
  
  AwarenessOverallASDSlowerTDC_Q_F <- extract_relevant_data(model,
                                                            "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled",
                                                            "AwarenessOverallASDSlowerTDC_Q_F")
  
  AwarenessOverallASDSlowerTDC_Q_UF <- extract_relevant_data(model,
                                                             "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled",
                                                             "AwarenessOverallASDSlowerTDC_Q_UF")
  
  MotivationOverallASDSlowerTDC_MG <- extract_relevant_data(model,
                                                            "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled > DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled",
                                                            "MotivationOverallASDSlowerTDC_MG")
  
  MotivationOverallASDSlowerTDC_Q_F <- extract_relevant_data(model,
                                                             "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled",
                                                             "MotivationOverallASDSlowerTDC_Q_F")
  
  MotivationOverallASDSlowerTDC_Q_UF <- extract_relevant_data(model,
                                                              "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled",
                                                              "MotivationOverallASDSlowerTDC_Q_UF")
  
  
  DifferenceInDifferenceLanguageMGSlowerQ_F <- extract_relevant_data(model,
                                                                     "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled)",
                                                                     "DifferenceInDifferenceLanguageMGSlowerQ_F")
  
  DifferenceInDifferenceLanguageUFSlowerQ_F <- extract_relevant_data(model,
                                                                     "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled)",
                                                                     "DifferenceInDifferenceLanguageUFSlowerQ_F")
  
  DifferenceInDifferenceMotorMGSlowerQ_F <- extract_relevant_data(model,
                                                                  "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled)",
                                                                  "DifferenceInDifferenceMotorMGSlowerQ_F")
  
  DifferenceInDifferenceMotorUFSlowerQ_F <- extract_relevant_data(model,
                                                                  "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled)",
                                                                  "DifferenceInDifferenceMotorUFSlowerQ_F")
  
  DifferenceInDifferenceCognitionMGSlowerQ_F <- extract_relevant_data(model,
                                                                      "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled)",
                                                                      "DifferenceInDifferenceCognitionMGSlowerQ_F")
  
  DifferenceInDifferenceCognitionUFSlowerQ_F <- extract_relevant_data(model,
                                                                      "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled)",
                                                                      "DifferenceInDifferenceCognitionUFSlowerQ_F")
  
  DifferenceInDifferenceAwarenessMGSlowerQ_F <- extract_relevant_data(model,
                                                                      "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled)
                      < (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled)",
                                                                      "DifferenceInDifferenceAwarenessMGSlowerQ_F")
  
  DifferenceInDifferenceAwarenessUFSlowerQ_F <- extract_relevant_data(model,
                                                                      "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled)
                      < (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled)",
                                                                      "DifferenceInDifferenceAwarenessUFSlowerQ_F")
  
  DifferenceInDifferenceMotivationMGSlowerQ_F <- extract_relevant_data(model,
                                                                       "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled)",
                                                                       "DifferenceInDifferenceMotivationMGSlowerQ_F")
  
  DifferenceInDifferenceMotivationUFSlowerQ_F <- extract_relevant_data(model,
                                                                       "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled)
                      < (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled)",
                                                                       "DifferenceInDifferenceMotivationUFSlowerQ_F")
  
  LanguageOverallASDSlowerTDC_Q_F <- extract_relevant_data(model,
                                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled",
                                                           "LanguageOverallASDSlowerTDC_Q_F")
  
  LanguageOverallASDSlowerTDC_Q_UF <- extract_relevant_data(model,
                                                            "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled",
                                                            "LanguageOverallASDSlowerTDC_Q_UF")
  
  hypothesis_data_full <- rbind(
    MotorASD_TG_F, MotorASD_Q_F, MotorASD_Q_UF, MotorASD_overall,
    MotorTDC_TG_F, MotorTDC_Q_F, MotorTDC_Q_UF, MotorTDC_overall,
    
    CognitionASD_TG_F, CognitionASD_Q_F, CognitionASD_Q_UF, CognitionASD_overall,
    CognitionTDC_TG_F, CognitionTDC_Q_F, CognitionTDC_Q_UF, CognitionTDC_overall,
    
    LanguageASD_TG_F, LanguageASD_Q_F, LanguageASD_Q_UF, LanguageASD_overall,
    LanguageTDC_TG_F, LanguageTDC_Q_F, LanguageTDC_Q_UF, LanguageTDC_overall,
    
    AwarenessASD_TG_F, AwarenessASD_Q_F, AwarenessASD_Q_UF, AwarenessASD_overall,
    AwarenessTDC_TG_F, AwarenessTDC_Q_F, AwarenessTDC_Q_UF, AwarenessTDC_overall,
    
    MotivationASD_TG_F, MotivationASD_Q_F, MotivationASD_Q_UF, MotivationASD_overall,
    MotivationTDC_TG_F, MotivationTDC_Q_F, MotivationTDC_Q_UF, MotivationTDC_overall,
    
    LanguageOverallASDSlowerTDC, MotorOverallASDSlowerTDC, CognitionOverallASDSlowerTDC, AwarenessOverallASDFasterTDC, MotivationOverallASDFasterTDC,
    
    LanguageOverallASDSlowerTDC_MG, LanguageOverallASDSlowerTDC_Q_F, LanguageOverallASDSlowerTDC_Q_UF,
    CognitionOverallASDSlowerTDC_MG, CognitionOverallASDSlowerTDC_Q_F, CognitionOverallASDSlowerTDC_Q_UF,
    MotorOverallASDSlowerTDC_MG, MotorOverallASDSlowerTDC_Q_F, MotorOverallASDSlowerTDC_Q_UF,
    AwarenessOverallASDSlowerTDC_MG, AwarenessOverallASDSlowerTDC_Q_F, AwarenessOverallASDSlowerTDC_Q_UF,
    MotivationOverallASDSlowerTDC_MG, MotivationOverallASDSlowerTDC_Q_F, MotivationOverallASDSlowerTDC_Q_UF,
    
    DifferenceInDifferenceLanguageMGSlowerQ_F, DifferenceInDifferenceLanguageUFSlowerQ_F,
    DifferenceInDifferenceMotorMGSlowerQ_F, DifferenceInDifferenceMotorUFSlowerQ_F,
    DifferenceInDifferenceMotivationMGSlowerQ_F, DifferenceInDifferenceMotivationUFSlowerQ_F,
    DifferenceInDifferenceCognitionMGSlowerQ_F, DifferenceInDifferenceCognitionUFSlowerQ_F,
    DifferenceInDifferenceAwarenessMGSlowerQ_F, DifferenceInDifferenceAwarenessUFSlowerQ_F,
    
    WithinGroupLanguageASD_MGSlower_Q_F, WithinGroupLanguageASD_Q_UFSlower_Q_F,
    WithinGroupCognitionASD_MGSlower_Q_F, WithinGroupCognitionASD_Q_UFSlower_Q_F,
    WithinGroupMotorASD_MGSlower_Q_F, WithinGroupMotorASD_Q_UFSlower_Q_F,
    WithinGroupAwarenessASD_MGSlower_Q_F, WithinGroupAwarenessASD_Q_UFSlower_Q_F,
    WithinGroupMotivationASD_MGSlower_Q_F, WithinGroupMotivationASD_Q_UFSlower_Q_F,
    
    WithinGroupLanguageTDC_MGSlower_Q_F, WithinGroupLanguageTDC_Q_UFSlower_Q_F,
    WithinGroupCognitionTDC_MGSlower_Q_F, WithinGroupCognitionTDC_Q_UFSlower_Q_F,
    WithinGroupMotorTDC_MGSlower_Q_F, WithinGroupMotorTDC_Q_UFSlower_Q_F,
    WithinGroupAwarenessTDC_MGSlower_Q_F, WithinGroupAwarenessTDC_Q_UFSlower_Q_F,
    WithinGroupMotivationTDC_MGSlower_Q_F, WithinGroupMotivationTDC_Q_UFSlower_Q_F
    
    
  ) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  
  return(hypothesis_data_full)
}

extract_all_hypothesis_tests_m3_child <- function(model) {
  
  ASD_OtherPredictability_MG_F <- extract_relevant_data(model,
                                                        "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability < 0",
                                                        "ASD_OtherPredictability_MG_F")
  
  TDC_OtherPredictability_MG_F <- extract_relevant_data(model, 
                                                        "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability < 0",
                                                        "TDC_OtherPredictability_MG_F")
  
  ASD_OtherPredictability_Q_F <- extract_relevant_data(model, 
                                                       "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability < 0",
                                                       "ASD_OtherPredictability_Q_F")
  TDC_OtherPredictability_Q_F <- extract_relevant_data(model, 
                                                       "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability < 0",
                                                       "TDC_OtherPredictability_Q_F")
  
  ASD_OtherPredictability_Q_UF <- extract_relevant_data(model, 
                                                        "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability < 0",
                                                        "ASD_OtherPredictability_Q_UF")
  TDC_OtherPredictability_Q_UF <- extract_relevant_data(model,
                                                        "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability < 0",
                                                        "TDC_OtherPredictability_Q_UF")
  
  ASD_OtherPredictability_overall <- extract_relevant_data(model, 
                                                           "((DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability) / 3) < 0", "ASD_OtherPredictability_overall")
  
  TDC_OtherPredictability_overall <- extract_relevant_data(model,
                                                           "((DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability) / 3) < 0",
                                                           "TDC_OtherPredictability_overall")
  
  WithinGroupASD_OtherPredictability_MG_Higher_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability", 
                                                                            "WithinGroupASD_OtherPredictability_MG_Higher_Q_F")
  
  WithinGroupASD_OtherPredictability_Q_F_Higher_Q_UF <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability < DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability", 
                                                                              "WithinGroupASD_OtherPredictability_Q_F_Higher_Q_UF")
  
  WithinGroupTDC_OtherPredictability_MG_Higher_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability", 
                                                                            "WithinGroupTDC_OtherPredictability_MG_Higher_Q_F")
  
  WithinGroupTDC_OtherPredictability_Q_F_Higher_Q_UF <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability", 
                                                                              "WithinGroupTDC_OtherPredictability_Q_F_Higher_Q_UF")
  
  OtherPredictabilityHypothesesFull <- rbind(ASD_OtherPredictability_MG_F, ASD_OtherPredictability_Q_F, ASD_OtherPredictability_Q_UF,
                                             TDC_OtherPredictability_MG_F, TDC_OtherPredictability_Q_F, TDC_OtherPredictability_Q_UF,
                                             ASD_OtherPredictability_overall, TDC_OtherPredictability_overall, WithinGroupASD_OtherPredictability_MG_Higher_Q_F, WithinGroupASD_OtherPredictability_Q_F_Higher_Q_UF, WithinGroupTDC_OtherPredictability_MG_Higher_Q_F, WithinGroupTDC_OtherPredictability_Q_F_Higher_Q_UF)
  
  
  OtherPredictability_MG_F_ASDHigherTDC <- extract_relevant_data(model,
                                                                 "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability < 
           DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability",
                                                                 "OtherPredictability_MG_F_ASDHigherTDC")
  
  OtherPredictability_Q_F_ASDHigherTDC <- extract_relevant_data(model,
                                                                "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability > 
           DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability",
                                                                "OtherPredictability_Q_F_ASDHigherTDC")
  
  OtherPredictability_Q_UF_ASDHigherTDC <- extract_relevant_data(model,
                                                                 "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability > 
           DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability",
                                                                 "OtherPredictability_Q_UF_ASDHigherTDC")
  
  OtherPredictability_overall_ASDHigherTDC <- extract_relevant_data(model, 
                                                                    "((DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability +
                                                       DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability +
                                                       DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability) / 3) >
                                                       ((DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability +
                                                       DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability +
                                                       DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability) / 3)",
                                                                    "OtherPredictability_overall_ASDHigherTDC")
  
  DifferenceOtherPredictability_ASDLowerMatchingGameTDCVersusQuestions <- extract_relevant_data(model,
                                                                                                "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability) <
            (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability)",
                                                                                                "DifferenceOtherPredictability_ASDLowerMatchingGameTDCVersusQuestions")
  
  DifferenceOtherPredictability_ASDLowerFamiliarityTDCVersusUnfamiliarity <- extract_relevant_data(model,
                                                                                                   "(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability) >
            (DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability)",
                                                                                                   "DifferenceOtherPredictability_ASDHigherFamiliarityTDCVersusUnfamiliarity")
  
  OtherPredictabilityComparisonsHypothesesFull <- rbind(
    OtherPredictabilityHypothesesFull,
    OtherPredictability_MG_F_ASDHigherTDC,
    OtherPredictability_Q_F_ASDHigherTDC,
    OtherPredictability_Q_UF_ASDHigherTDC,
    OtherPredictability_overall_ASDHigherTDC,
    DifferenceOtherPredictability_ASDLowerMatchingGameTDCVersusQuestions,
    DifferenceOtherPredictability_ASDLowerFamiliarityTDCVersusUnfamiliarity) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 1 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -1 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 1 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 1 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -1 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -1 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  
  #Let's add the same hypotheses for OwnPredictability:
  ASD_OwnPredictability_MG_F <- extract_relevant_data(model,
                                                      "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability < 0",
                                                      "ASD_OwnPredictability_MG_F")
  
  TDC_OwnPredictability_MG_F <- extract_relevant_data(model, 
                                                      "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability > 0",
                                                      "TDC_OwnPredictability_MG_F")
  
  ASD_OwnPredictability_Q_F <- extract_relevant_data(model, 
                                                     "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability > 0",
                                                     "ASD_OwnPredictability_Q_F")
  TDC_OwnPredictability_Q_F <- extract_relevant_data(model, 
                                                     "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability > 0",
                                                     "TDC_OwnPredictability_Q_F")
  
  ASD_OwnPredictability_Q_UF <- extract_relevant_data(model, 
                                                      "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability > 0",
                                                      "ASD_OwnPredictability_Q_UF")
  TDC_OwnPredictability_Q_UF <- extract_relevant_data(model,
                                                      "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability > 0",
                                                      "TDC_OwnPredictability_Q_UF")
  
  
  ASD_OwnPredictability_overall <- extract_relevant_data(model, 
                                                         "((DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability) / 3) > 0", "ASD_OwnPredictability_overall")
  
  TDC_OwnPredictability_overall <- extract_relevant_data(model,
                                                         "((DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability) / 3) > 0",
                                                         "TDC_OwnPredictability_overall")
  
  WithinGroupASD_OwnPredictability_MG_Higher_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability < DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability", 
                                                                          "WithinGroupASD_OwnPredictability_MG_Higher_Q_F")
  
  WithinGroupASD_OwnPredictability_Q_F_Higher_Q_UF <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability < DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability", 
                                                                            "WithinGroupASD_OwnPredictability_Q_F_Higher_Q_UF")
  
  WithinGroupTDC_OwnPredictability_MG_Higher_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability", 
                                                                          "WithinGroupTDC_OwnPredictability_MG_Higher_Q_F")
  
  WithinGroupTDC_OwnPredictability_Q_F_Higher_Q_UF <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability", 
                                                                            "WithinGroupTDC_OwnPredictability_Q_F_Higher_Q_UF")
  
  
  OwnPredictabilityHypothesesFull <- rbind(ASD_OwnPredictability_MG_F, ASD_OwnPredictability_Q_F, ASD_OwnPredictability_Q_UF,
                                           TDC_OwnPredictability_MG_F, TDC_OwnPredictability_Q_F, TDC_OwnPredictability_Q_UF,
                                           ASD_OwnPredictability_overall, TDC_OwnPredictability_overall, WithinGroupASD_OwnPredictability_MG_Higher_Q_F, WithinGroupASD_OwnPredictability_Q_F_Higher_Q_UF, WithinGroupTDC_OwnPredictability_MG_Higher_Q_F, WithinGroupTDC_OwnPredictability_Q_F_Higher_Q_UF)
  
  
  OwnPredictability_MG_F_ASDHigherTDC <- extract_relevant_data(model,
                                                               "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability < 
           DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability",
                                                               "OwnPredictability_MG_F_ASDHigherTDC")
  
  OwnPredictability_Q_F_ASDHigherTDC <- extract_relevant_data(model,
                                                              "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability > 
           DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability",
                                                              "OwnPredictability_Q_F_ASDHigherTDC")
  
  OwnPredictability_Q_UF_ASDHigherTDC <- extract_relevant_data(model,
                                                               "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability > 
           DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability",
                                                               "OwnPredictability_Q_UF_ASDHigherTDC")
  
  OwnPredictability_overall_ASDHigherTDC <- extract_relevant_data(model, 
                                                                  "((DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability +
                                                       DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability +
                                                       DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability) / 3) >
                                                       ((DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability +
                                                       DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability +
                                                       DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability) / 3)",
                                                                  "OwnPredictability_overall_ASDHigherTDC")
  
  DifferenceOwnPredictability_ASDLowerMatchingGameTDCVersusQuestions <- extract_relevant_data(model,
                                                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability) <
            (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability)",
                                                                                              "DifferenceOwnPredictability_ASDLowerMatchingGameTDCVersusQuestions")
  
  DifferenceOwnPredictability_ASDLowerFamiliarityTDCVersusUnfamiliarity <- extract_relevant_data(model,
                                                                                                 "(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability) <
            (DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability)",
                                                                                                 "DifferenceOwnPredictability_ASDHigherFamiliarityTDCVersusUnfamiliarity")
  
  OwnPredictabilityComparisonsHypothesesFull <- rbind(
    OwnPredictabilityHypothesesFull,
    OwnPredictability_MG_F_ASDHigherTDC,
    OwnPredictability_Q_F_ASDHigherTDC,
    OwnPredictability_Q_UF_ASDHigherTDC,
    OwnPredictability_overall_ASDHigherTDC,
    DifferenceOwnPredictability_ASDLowerMatchingGameTDCVersusQuestions,
    DifferenceOwnPredictability_ASDLowerFamiliarityTDCVersusUnfamiliarity) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 1 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -1 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 1 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 1 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -1 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -1 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  
  PredictabilityComparisonsHypothesesFull <- rbind(OtherPredictabilityComparisonsHypothesesFull, OwnPredictabilityComparisonsHypothesesFull)
  
  return(PredictabilityComparisonsHypothesesFull)
}

extract_all_hypothesis_tests_m3_child <- function(model) {
  
  ASD_OtherPredictability_MG_F <- extract_relevant_data(model,
                                                        "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability < 0",
                                                        "ASD_OtherPredictability_MG_F")
  
  TDC_OtherPredictability_MG_F <- extract_relevant_data(model, 
                                                        "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability < 0",
                                                        "TDC_OtherPredictability_MG_F")
  
  ASD_OtherPredictability_Q_F <- extract_relevant_data(model, 
                                                       "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability < 0",
                                                       "ASD_OtherPredictability_Q_F")
  TDC_OtherPredictability_Q_F <- extract_relevant_data(model, 
                                                       "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability < 0",
                                                       "TDC_OtherPredictability_Q_F")
  
  ASD_OtherPredictability_Q_UF <- extract_relevant_data(model, 
                                                        "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability < 0",
                                                        "ASD_OtherPredictability_Q_UF")
  TDC_OtherPredictability_Q_UF <- extract_relevant_data(model,
                                                        "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability < 0",
                                                        "TDC_OtherPredictability_Q_UF")
  
  ASD_OtherPredictability_overall <- extract_relevant_data(model, 
                                                           "((DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability) / 3) < 0", "ASD_OtherPredictability_overall")
  
  TDC_OtherPredictability_overall <- extract_relevant_data(model,
                                                           "((DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability) / 3) < 0",
                                                           "TDC_OtherPredictability_overall")
  
  WithinGroupASD_OtherPredictability_MG_Higher_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability", 
                                                                            "WithinGroupASD_OtherPredictability_MG_Higher_Q_F")
  
  WithinGroupASD_OtherPredictability_Q_F_Higher_Q_UF <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability < DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability", 
                                                                              "WithinGroupASD_OtherPredictability_Q_F_Higher_Q_UF")
  
  WithinGroupTDC_OtherPredictability_MG_Higher_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability", 
                                                                            "WithinGroupTDC_OtherPredictability_MG_Higher_Q_F")
  
  WithinGroupTDC_OtherPredictability_Q_F_Higher_Q_UF <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability", 
                                                                              "WithinGroupTDC_OtherPredictability_Q_F_Higher_Q_UF")
  
  OtherPredictabilityHypothesesFull <- rbind(ASD_OtherPredictability_MG_F, ASD_OtherPredictability_Q_F, ASD_OtherPredictability_Q_UF,
                                             TDC_OtherPredictability_MG_F, TDC_OtherPredictability_Q_F, TDC_OtherPredictability_Q_UF,
                                             ASD_OtherPredictability_overall, TDC_OtherPredictability_overall, WithinGroupASD_OtherPredictability_MG_Higher_Q_F, WithinGroupASD_OtherPredictability_Q_F_Higher_Q_UF, WithinGroupTDC_OtherPredictability_MG_Higher_Q_F, WithinGroupTDC_OtherPredictability_Q_F_Higher_Q_UF)
  
  
  OtherPredictability_MG_F_ASDHigherTDC <- extract_relevant_data(model,
                                                                 "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability < 
           DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability",
                                                                 "OtherPredictability_MG_F_ASDHigherTDC")
  
  OtherPredictability_Q_F_ASDHigherTDC <- extract_relevant_data(model,
                                                                "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability > 
           DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability",
                                                                "OtherPredictability_Q_F_ASDHigherTDC")
  
  OtherPredictability_Q_UF_ASDHigherTDC <- extract_relevant_data(model,
                                                                 "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability > 
           DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability",
                                                                 "OtherPredictability_Q_UF_ASDHigherTDC")
  
  OtherPredictability_overall_ASDHigherTDC <- extract_relevant_data(model, 
                                                                    "((DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability +
                                                       DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability +
                                                       DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability) / 3) >
                                                       ((DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability +
                                                       DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability +
                                                       DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability) / 3)",
                                                                    "OtherPredictability_overall_ASDHigherTDC")
  
  DifferenceOtherPredictability_ASDLowerMatchingGameTDCVersusQuestions <- extract_relevant_data(model,
                                                                                                "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability) <
            (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability)",
                                                                                                "DifferenceOtherPredictability_ASDLowerMatchingGameTDCVersusQuestions")
  
  DifferenceOtherPredictability_ASDLowerFamiliarityTDCVersusUnfamiliarity <- extract_relevant_data(model,
                                                                                                   "(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability) >
            (DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability)",
                                                                                                   "DifferenceOtherPredictability_ASDHigherFamiliarityTDCVersusUnfamiliarity")
  
  OtherPredictabilityComparisonsHypothesesFull <- rbind(
    OtherPredictabilityHypothesesFull,
    OtherPredictability_MG_F_ASDHigherTDC,
    OtherPredictability_Q_F_ASDHigherTDC,
    OtherPredictability_Q_UF_ASDHigherTDC,
    OtherPredictability_overall_ASDHigherTDC,
    DifferenceOtherPredictability_ASDLowerMatchingGameTDCVersusQuestions,
    DifferenceOtherPredictability_ASDLowerFamiliarityTDCVersusUnfamiliarity) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 1 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -1 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 1 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 1 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -1 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -1 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  
  #Let's add the same hypotheses for OwnPredictability:
  ASD_OwnPredictability_MG_F <- extract_relevant_data(model,
                                                      "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability < 0",
                                                      "ASD_OwnPredictability_MG_F")
  
  TDC_OwnPredictability_MG_F <- extract_relevant_data(model, 
                                                      "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability > 0",
                                                      "TDC_OwnPredictability_MG_F")
  
  ASD_OwnPredictability_Q_F <- extract_relevant_data(model, 
                                                     "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability > 0",
                                                     "ASD_OwnPredictability_Q_F")
  TDC_OwnPredictability_Q_F <- extract_relevant_data(model, 
                                                     "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability > 0",
                                                     "TDC_OwnPredictability_Q_F")
  
  ASD_OwnPredictability_Q_UF <- extract_relevant_data(model, 
                                                      "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability > 0",
                                                      "ASD_OwnPredictability_Q_UF")
  TDC_OwnPredictability_Q_UF <- extract_relevant_data(model,
                                                      "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability > 0",
                                                      "TDC_OwnPredictability_Q_UF")
  
  
  ASD_OwnPredictability_overall <- extract_relevant_data(model, 
                                                         "((DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability) / 3) > 0", "ASD_OwnPredictability_overall")
  
  TDC_OwnPredictability_overall <- extract_relevant_data(model,
                                                         "((DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability) / 3) > 0",
                                                         "TDC_OwnPredictability_overall")
  
  WithinGroupASD_OwnPredictability_MG_Higher_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability < DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability", 
                                                                          "WithinGroupASD_OwnPredictability_MG_Higher_Q_F")
  
  WithinGroupASD_OwnPredictability_Q_F_Higher_Q_UF <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability < DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability", 
                                                                            "WithinGroupASD_OwnPredictability_Q_F_Higher_Q_UF")
  
  WithinGroupTDC_OwnPredictability_MG_Higher_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability", 
                                                                          "WithinGroupTDC_OwnPredictability_MG_Higher_Q_F")
  
  WithinGroupTDC_OwnPredictability_Q_F_Higher_Q_UF <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability", 
                                                                            "WithinGroupTDC_OwnPredictability_Q_F_Higher_Q_UF")
  
  
  OwnPredictabilityHypothesesFull <- rbind(ASD_OwnPredictability_MG_F, ASD_OwnPredictability_Q_F, ASD_OwnPredictability_Q_UF,
                                           TDC_OwnPredictability_MG_F, TDC_OwnPredictability_Q_F, TDC_OwnPredictability_Q_UF,
                                           ASD_OwnPredictability_overall, TDC_OwnPredictability_overall, WithinGroupASD_OwnPredictability_MG_Higher_Q_F, WithinGroupASD_OwnPredictability_Q_F_Higher_Q_UF, WithinGroupTDC_OwnPredictability_MG_Higher_Q_F, WithinGroupTDC_OwnPredictability_Q_F_Higher_Q_UF)
  
  
  OwnPredictability_MG_F_ASDHigherTDC <- extract_relevant_data(model,
                                                               "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability < 
           DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability",
                                                               "OwnPredictability_MG_F_ASDHigherTDC")
  
  OwnPredictability_Q_F_ASDHigherTDC <- extract_relevant_data(model,
                                                              "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability > 
           DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability",
                                                              "OwnPredictability_Q_F_ASDHigherTDC")
  
  OwnPredictability_Q_UF_ASDHigherTDC <- extract_relevant_data(model,
                                                               "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability > 
           DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability",
                                                               "OwnPredictability_Q_UF_ASDHigherTDC")
  
  OwnPredictability_overall_ASDHigherTDC <- extract_relevant_data(model, 
                                                                  "((DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability +
                                                       DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability +
                                                       DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability) / 3) >
                                                       ((DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability +
                                                       DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability +
                                                       DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability) / 3)",
                                                                  "OwnPredictability_overall_ASDHigherTDC")
  
  DifferenceOwnPredictability_ASDLowerMatchingGameTDCVersusQuestions <- extract_relevant_data(model,
                                                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability) <
            (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability)",
                                                                                              "DifferenceOwnPredictability_ASDLowerMatchingGameTDCVersusQuestions")
  
  DifferenceOwnPredictability_ASDLowerFamiliarityTDCVersusUnfamiliarity <- extract_relevant_data(model,
                                                                                                 "(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability) <
            (DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability)",
                                                                                                 "DifferenceOwnPredictability_ASDHigherFamiliarityTDCVersusUnfamiliarity")
  
  OwnPredictabilityComparisonsHypothesesFull <- rbind(
    OwnPredictabilityHypothesesFull,
    OwnPredictability_MG_F_ASDHigherTDC,
    OwnPredictability_Q_F_ASDHigherTDC,
    OwnPredictability_Q_UF_ASDHigherTDC,
    OwnPredictability_overall_ASDHigherTDC,
    DifferenceOwnPredictability_ASDLowerMatchingGameTDCVersusQuestions,
    DifferenceOwnPredictability_ASDLowerFamiliarityTDCVersusUnfamiliarity) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 1 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -1 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 1 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 1 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -1 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -1 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  
  PredictabilityComparisonsHypothesesFull <- rbind(OtherPredictabilityComparisonsHypothesesFull, OwnPredictabilityComparisonsHypothesesFull)
  
  return(PredictabilityComparisonsHypothesesFull)
}

extract_all_hypothesis_tests_m3_adult <- function(model) {
  
  ASD_OtherPredictability_MG_F <- extract_relevant_data(model,
                                                        "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability < 0",
                                                        "ASD_OtherPredictability_MG_F")
  
  TDC_OtherPredictability_MG_F <- extract_relevant_data(model, 
                                                        "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability < 0",
                                                        "TDC_OtherPredictability_MG_F")
  
  ASD_OtherPredictability_Q_F <- extract_relevant_data(model, 
                                                       "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability < 0",
                                                       "ASD_OtherPredictability_Q_F")
  TDC_OtherPredictability_Q_F <- extract_relevant_data(model, 
                                                       "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability < 0",
                                                       "TDC_OtherPredictability_Q_F")
  
  ASD_OtherPredictability_Q_UF <- extract_relevant_data(model, 
                                                        "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability < 0",
                                                        "ASD_OtherPredictability_Q_UF")
  TDC_OtherPredictability_Q_UF <- extract_relevant_data(model,
                                                        "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability < 0",
                                                        "TDC_OtherPredictability_Q_UF")
  
  ASD_OtherPredictability_overall <- extract_relevant_data(model, 
                                                           "((DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability) / 3) < 0", "ASD_OtherPredictability_overall")
  
  TDC_OtherPredictability_overall <- extract_relevant_data(model,
                                                           "((DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability) / 3) < 0",
                                                           "TDC_OtherPredictability_overall")
  
  WithinGroupASD_OtherPredictability_MG_Higher_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability", 
                                                                            "WithinGroupASD_OtherPredictability_MG_Higher_Q_F")
  
  WithinGroupASD_OtherPredictability_Q_F_Higher_Q_UF <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability < DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability", 
                                                                              "WithinGroupASD_OtherPredictability_Q_F_Higher_Q_UF")
  
  WithinGroupTDC_OtherPredictability_MG_Higher_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability", 
                                                                            "WithinGroupTDC_OtherPredictability_MG_Higher_Q_F")
  
  WithinGroupTDC_OtherPredictability_Q_F_Higher_Q_UF <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability", 
                                                                              "WithinGroupTDC_OtherPredictability_Q_F_Higher_Q_UF")
  
  OtherPredictabilityHypothesesFull <- rbind(ASD_OtherPredictability_MG_F, ASD_OtherPredictability_Q_F, ASD_OtherPredictability_Q_UF,
                                             TDC_OtherPredictability_MG_F, TDC_OtherPredictability_Q_F, TDC_OtherPredictability_Q_UF,
                                             ASD_OtherPredictability_overall, TDC_OtherPredictability_overall, WithinGroupASD_OtherPredictability_MG_Higher_Q_F, WithinGroupASD_OtherPredictability_Q_F_Higher_Q_UF, WithinGroupTDC_OtherPredictability_MG_Higher_Q_F, WithinGroupTDC_OtherPredictability_Q_F_Higher_Q_UF)
  
  
  OtherPredictability_MG_F_ASDHigherTDC <- extract_relevant_data(model,
                                                                 "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability < 
           DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability",
                                                                 "OtherPredictability_MG_F_ASDHigherTDC")
  
  OtherPredictability_Q_F_ASDHigherTDC <- extract_relevant_data(model,
                                                                "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability > 
           DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability",
                                                                "OtherPredictability_Q_F_ASDHigherTDC")
  
  OtherPredictability_Q_UF_ASDHigherTDC <- extract_relevant_data(model,
                                                                 "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability > 
           DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability",
                                                                 "OtherPredictability_Q_UF_ASDHigherTDC")
  
  OtherPredictability_overall_ASDHigherTDC <- extract_relevant_data(model, 
                                                                    "((DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability +
                                                       DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability +
                                                       DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability) / 3) >
                                                       ((DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability +
                                                       DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability +
                                                       DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability) / 3)",
                                                                    "OtherPredictability_overall_ASDHigherTDC")
  
  DifferenceOtherPredictability_ASDLowerMatchingGameTDCVersusQuestions <- extract_relevant_data(model,
                                                                                                "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability) <
            (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability)",
                                                                                                "DifferenceOtherPredictability_ASDLowerMatchingGameTDCVersusQuestions")
  
  DifferenceOtherPredictability_ASDLowerFamiliarityTDCVersusUnfamiliarity <- extract_relevant_data(model,
                                                                                                   "(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability) >
            (DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability)",
                                                                                                   "DifferenceOtherPredictability_ASDHigherFamiliarityTDCVersusUnfamiliarity")
  
  OtherPredictabilityComparisonsHypothesesFull <- rbind(
    OtherPredictabilityHypothesesFull,
    OtherPredictability_MG_F_ASDHigherTDC,
    OtherPredictability_Q_F_ASDHigherTDC,
    OtherPredictability_Q_UF_ASDHigherTDC,
    OtherPredictability_overall_ASDHigherTDC,
    DifferenceOtherPredictability_ASDLowerMatchingGameTDCVersusQuestions,
    DifferenceOtherPredictability_ASDLowerFamiliarityTDCVersusUnfamiliarity) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  
  #Let's add the same hypotheses for OwnPredictability:
  ASD_OwnPredictability_MG_F <- extract_relevant_data(model,
                                                      "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability < 0",
                                                      "ASD_OwnPredictability_MG_F")
  
  TDC_OwnPredictability_MG_F <- extract_relevant_data(model, 
                                                      "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability > 0",
                                                      "TDC_OwnPredictability_MG_F")
  
  ASD_OwnPredictability_Q_F <- extract_relevant_data(model, 
                                                     "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability > 0",
                                                     "ASD_OwnPredictability_Q_F")
  TDC_OwnPredictability_Q_F <- extract_relevant_data(model, 
                                                     "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability > 0",
                                                     "TDC_OwnPredictability_Q_F")
  
  ASD_OwnPredictability_Q_UF <- extract_relevant_data(model, 
                                                      "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability > 0",
                                                      "ASD_OwnPredictability_Q_UF")
  TDC_OwnPredictability_Q_UF <- extract_relevant_data(model,
                                                      "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability > 0",
                                                      "TDC_OwnPredictability_Q_UF")
  
  
  ASD_OwnPredictability_overall <- extract_relevant_data(model, 
                                                         "((DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability) / 3) > 0", "ASD_OwnPredictability_overall")
  
  TDC_OwnPredictability_overall <- extract_relevant_data(model,
                                                         "((DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability) / 3) > 0",
                                                         "TDC_OwnPredictability_overall")
  
  WithinGroupASD_OwnPredictability_MG_Higher_Q_F <- extract_relevant_data(model, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability < DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability", 
                                                                          "WithinGroupASD_OwnPredictability_MG_Higher_Q_F")
  
  WithinGroupASD_OwnPredictability_Q_F_Higher_Q_UF <- extract_relevant_data(model, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability < DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability", 
                                                                            "WithinGroupASD_OwnPredictability_Q_F_Higher_Q_UF")
  
  WithinGroupTDC_OwnPredictability_MG_Higher_Q_F <- extract_relevant_data(model, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability", 
                                                                          "WithinGroupTDC_OwnPredictability_MG_Higher_Q_F")
  
  WithinGroupTDC_OwnPredictability_Q_F_Higher_Q_UF <- extract_relevant_data(model, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability", 
                                                                            "WithinGroupTDC_OwnPredictability_Q_F_Higher_Q_UF")
  
  
  OwnPredictabilityHypothesesFull <- rbind(ASD_OwnPredictability_MG_F, ASD_OwnPredictability_Q_F, ASD_OwnPredictability_Q_UF,
                                           TDC_OwnPredictability_MG_F, TDC_OwnPredictability_Q_F, TDC_OwnPredictability_Q_UF,
                                           ASD_OwnPredictability_overall, TDC_OwnPredictability_overall, WithinGroupASD_OwnPredictability_MG_Higher_Q_F, WithinGroupASD_OwnPredictability_Q_F_Higher_Q_UF, WithinGroupTDC_OwnPredictability_MG_Higher_Q_F, WithinGroupTDC_OwnPredictability_Q_F_Higher_Q_UF)
  
  
  OwnPredictability_MG_F_ASDHigherTDC <- extract_relevant_data(model,
                                                               "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability < 
           DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability",
                                                               "OwnPredictability_MG_F_ASDHigherTDC")
  
  OwnPredictability_Q_F_ASDHigherTDC <- extract_relevant_data(model,
                                                              "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability > 
           DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability",
                                                              "OwnPredictability_Q_F_ASDHigherTDC")
  
  OwnPredictability_Q_UF_ASDHigherTDC <- extract_relevant_data(model,
                                                               "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability > 
           DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability",
                                                               "OwnPredictability_Q_UF_ASDHigherTDC")
  
  OwnPredictability_overall_ASDHigherTDC <- extract_relevant_data(model, 
                                                                  "((DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability +
                                                       DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability +
                                                       DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability) / 3) >
                                                       ((DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability +
                                                       DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability +
                                                       DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability) / 3)",
                                                                  "OwnPredictability_overall_ASDHigherTDC")
  
  DifferenceOwnPredictability_ASDLowerMatchingGameTDCVersusQuestions <- extract_relevant_data(model,
                                                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability) <
            (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability)",
                                                                                              "DifferenceOwnPredictability_ASDLowerMatchingGameTDCVersusQuestions")
  
  DifferenceOwnPredictability_ASDLowerFamiliarityTDCVersusUnfamiliarity <- extract_relevant_data(model,
                                                                                                 "(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability) <
            (DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability)",
                                                                                                 "DifferenceOwnPredictability_ASDHigherFamiliarityTDCVersusUnfamiliarity")
  
  OwnPredictabilityComparisonsHypothesesFull <- rbind(
    OwnPredictabilityHypothesesFull,
    OwnPredictability_MG_F_ASDHigherTDC,
    OwnPredictability_Q_F_ASDHigherTDC,
    OwnPredictability_Q_UF_ASDHigherTDC,
    OwnPredictability_overall_ASDHigherTDC,
    DifferenceOwnPredictability_ASDLowerMatchingGameTDCVersusQuestions,
    DifferenceOwnPredictability_ASDLowerFamiliarityTDCVersusUnfamiliarity) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  
  PredictabilityComparisonsHypothesesFull <- rbind(OtherPredictabilityComparisonsHypothesesFull, OwnPredictabilityComparisonsHypothesesFull)
  return(PredictabilityComparisonsHypothesesFull)
}

extract_all_hypothesis_tests_m4_interpersonal_child <- function(model) {
  
  ASD_MG_F <- extract_relevant_data(model, 
                                    "ChildLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:PrevOtherLatency > 0",
                                    "ASD_MG_F")
  
  TDC_MG_F <- extract_relevant_data(model, 
                                    "ChildLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:PrevOtherLatency > 0",
                                    "TDC_MG_F")
  
  ASD_Q_F <- extract_relevant_data(model, 
                                   "ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:PrevOtherLatency > 0",
                                   "ASD_Q_F")
  
  TDC_Q_F <- extract_relevant_data(model, 
                                   "ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:PrevOtherLatency > 0",
                                   "TDC_Q_F")
  
  ASD_Q_UF <- extract_relevant_data(model, 
                                    "ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:PrevOtherLatency < 0",
                                    "ASD_Q_UF")
  
  TDC_Q_UF <- extract_relevant_data(model, 
                                    "ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:PrevOtherLatency > 0",
                                    "TDC_Q_UF")
  
  ASD_overall <- extract_relevant_data(model, 
                                       "((ChildLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:PrevOtherLatency +
           ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:PrevOtherLatency +
           ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:PrevOtherLatency)/3) > 0",
                                       "ASD_overall")
  
  TDC_overall <- extract_relevant_data(model, 
                                       "((ChildLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:PrevOtherLatency +
           ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:PrevOtherLatency +
           ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:PrevOtherLatency)/3) > 0",
                                       "TDC_overall")
  
  InterpersonalAdjustmentOverallHypotheses <- rbind(ASD_MG_F, ASD_Q_F, ASD_Q_UF, ASD_overall,
                                                    TDC_MG_F, TDC_Q_F, TDC_Q_UF, TDC_overall)
  
  InterpersonalAdjustment_MG_F_ASDHigherTDC <- extract_relevant_data(model,
                                                                     "ChildLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:PrevOtherLatency > 
           ChildLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:PrevOtherLatency",
                                                                     "InterpersonalAdjustment_MG_F_ASDHigherTDC")
  
  InterpersonalAdjustment_Q_F_ASDHigherTDC <- extract_relevant_data(model,
                                                                    "ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:PrevOtherLatency > 
           ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:PrevOtherLatency",
                                                                    "InterpersonalAdjustment_Q_F_ASDHigherTDC")
  
  InterpersonalAdjustment_Q_UF_ASDHigherTDC <- extract_relevant_data(model,
                                                                     "ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:PrevOtherLatency > 
           ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:PrevOtherLatency",
                                                                     "InterpersonalAdjustment_Q_UF_ASDHigherTDC")
  
  InterpersonalAdjustment_overall_ASDHigherTDC <- extract_relevant_data(model, 
                                                                        "((ChildLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:PrevOtherLatency +
                                                       ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:PrevOtherLatency +
                                                       ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:PrevOtherLatency) / 3) >
                                                       ((ChildLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:PrevOtherLatency +
                                                       ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:PrevOtherLatency +
                                                       ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:PrevOtherLatency) / 3)",
                                                                        "InterpersonalAdjustment_overall_ASDHigherTDC")
  
  DifferenceInterpersonalAdjustment_ASDLowerMatchingGameTDCVersusQuestions <- extract_relevant_data(model,
                                                                                                    "(ChildLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:PrevOtherLatency - ChildLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:PrevOtherLatency) <
            (ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:PrevOtherLatency - ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:PrevOtherLatency)",
                                                                                                    "DifferenceInterpersonalAdjustment_ASDLowerMatchingGameTDCVersusQuestions")
  
  DifferenceInterpersonalAdjustment_ASDHigherFamiliarityTDCVersusUnfamiliarity <- extract_relevant_data(model,
                                                                                                        "(ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:PrevOtherLatency - ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:PrevOtherLatency) >
            (ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:PrevOtherLatency - ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:PrevOtherLatency)",
                                                                                                        "DifferenceInterpersonalAdjustment_ASDHigherFamiliarityTDCVersusUnfamiliarity")
  
  TDCWithinGroupInterpersonalAdjustment_Q_F_HigherQ_UF <- extract_relevant_data(model,
                                                                                "ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:PrevOtherLatency > 
           ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:PrevOtherLatency",
                                                                                "TDCWithinGroupInterpersonalAdjustment_Q_F_HigherQ_UF")
  
  ASDWithinGroupInterpersonalAdjustment_Q_F_HigherQ_UF <- extract_relevant_data(model,
                                                                                "ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:PrevOtherLatency < ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:PrevOtherLatency",
                                                                                "ASDWithinGroupInterpersonalAdjustment_Q_F_HigherQ_UF")
  
  TDCWithinGroupInterpersonalAdjustment_MG_HigherQ_F <- extract_relevant_data(model,
                                                                              "ChildLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:PrevOtherLatency > 
           ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:PrevOtherLatency",
                                                                              "TDCWithinGroupInterpersonalAdjustment_MG_HigherQ_F")
  
  ASDWithinGroupInterpersonalAdjustment_MG_HigherQ_F <- extract_relevant_data(model,
                                                                              "ChildLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:PrevOtherLatency < 
           ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:PrevOtherLatency",
                                                                              "ASDWithinGroupInterpersonalAdjustment_MG_HigherQ_F")
  
  InterpersonalComparisonsHypothesesFull <- rbind(
    InterpersonalAdjustmentOverallHypotheses,
    InterpersonalAdjustment_MG_F_ASDHigherTDC,
    InterpersonalAdjustment_Q_F_ASDHigherTDC,
    InterpersonalAdjustment_Q_UF_ASDHigherTDC,
    InterpersonalAdjustment_overall_ASDHigherTDC,
    DifferenceInterpersonalAdjustment_ASDLowerMatchingGameTDCVersusQuestions,
    DifferenceInterpersonalAdjustment_ASDHigherFamiliarityTDCVersusUnfamiliarity,
    TDCWithinGroupInterpersonalAdjustment_Q_F_HigherQ_UF,
    ASDWithinGroupInterpersonalAdjustment_Q_F_HigherQ_UF,
    TDCWithinGroupInterpersonalAdjustment_MG_HigherQ_F,
    ASDWithinGroupInterpersonalAdjustment_MG_HigherQ_F) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 1 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -1 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 1 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 1 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -1 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -1 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  return(InterpersonalComparisonsHypothesesFull)
}

extract_all_hypothesis_tests_m4_self_child <- function(model) {
  
  ASD_MG_F <- extract_relevant_data(model, 
                                    "ChildLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:PrevSelfLatency > 0",
                                    "ASD_MG_F")
  
  TDC_MG_F <- extract_relevant_data(model, 
                                    "ChildLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:PrevSelfLatency > 0",
                                    "TDC_MG_F")
  
  ASD_Q_F <- extract_relevant_data(model, 
                                   "ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:PrevSelfLatency > 0",
                                   "ASD_Q_F")
  
  TDC_Q_F <- extract_relevant_data(model, 
                                   "ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:PrevSelfLatency > 0",
                                   "TDC_Q_F")
  
  ASD_Q_UF <- extract_relevant_data(model, 
                                    "ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:PrevSelfLatency > 0",
                                    "ASD_Q_UF")
  
  TDC_Q_UF <- extract_relevant_data(model, 
                                    "ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:PrevSelfLatency < 0",
                                    "TDC_Q_UF")
  
  ASD_overall <- extract_relevant_data(model, 
                                       "((ChildLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:PrevSelfLatency +
           ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:PrevSelfLatency +
           ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:PrevSelfLatency)/3) > 0",
                                       "ASD_overall")
  
  TDC_overall <- extract_relevant_data(model, 
                                       "((ChildLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:PrevSelfLatency +
           ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:PrevSelfLatency +
           ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:PrevSelfLatency)/3) > 0",
                                       "TDC_overall")
  
  SelfAdjustmentOverallHypotheses <- rbind(ASD_MG_F, ASD_Q_F, ASD_Q_UF, ASD_overall,
                                           TDC_MG_F, TDC_Q_F, TDC_Q_UF, TDC_overall)
  
  SelfAdjustment_MG_F_ASDHigherTDC <- extract_relevant_data(model,
                                                            "ChildLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:PrevSelfLatency < 
           ChildLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:PrevSelfLatency",
                                                            "SelfAdjustment_MG_F_ASDHigherTDC")
  
  SelfAdjustment_Q_F_ASDHigherTDC <- extract_relevant_data(model,
                                                           "ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:PrevSelfLatency > 
           ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:PrevSelfLatency",
                                                           "SelfAdjustment_Q_F_ASDHigherTDC")
  
  SelfAdjustment_Q_UF_ASDHigherTDC <- extract_relevant_data(model,
                                                            "ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:PrevSelfLatency > 
           ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:PrevSelfLatency",
                                                            "SelfAdjustment_Q_UF_ASDHigherTDC")
  
  SelfAdjustment_overall_ASDHigherTDC <- extract_relevant_data(model, 
                                                               "((ChildLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:PrevSelfLatency +
                                                       ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:PrevSelfLatency +
                                                       ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:PrevSelfLatency) / 3) >
                                                       ((ChildLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:PrevSelfLatency +
                                                       ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:PrevSelfLatency +
                                                       ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:PrevSelfLatency) / 3)",
                                                               "SelfAdjustment_overall_ASDHigherTDC")
  
  DifferenceSelfAdjustment_ASDLowerMatchingGameTDCVersusQuestions <- extract_relevant_data(model,
                                                                                           "(ChildLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:PrevSelfLatency - ChildLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:PrevSelfLatency) <
            (ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:PrevSelfLatency - ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:PrevSelfLatency)",
                                                                                           "DifferenceSelfAdjustment_ASDLowerMatchingGameTDCVersusQuestions")
  
  DifferenceSelfAdjustment_ASDHigherFamiliarityTDCVersusUnfamiliarity <- extract_relevant_data(model,
                                                                                               "(ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:PrevSelfLatency - ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:PrevSelfLatency) <
            (ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:PrevSelfLatency - ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:PrevSelfLatency)",
                                                                                               "DifferenceSelfAdjustment_ASDHigherFamiliarityTDCVersusUnfamiliarity")
  
  
  TDCWithinGroupSelfAdjustment_Q_F_HigherQ_UF <- extract_relevant_data(model,
                                                                       "ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:PrevSelfLatency > 
           ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:PrevSelfLatency",
                                                                       "TDCWithinGroupSelfAdjustment_Q_F_HigherQ_UF")
  
  ASDWithinGroupSelfAdjustment_Q_F_HigherQ_UF <- extract_relevant_data(model,
                                                                       "ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:PrevSelfLatency > 
           ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:PrevSelfLatency",
                                                                       "ASDWithinGroupSelfAdjustment_Q_F_HigherQ_UF")
  
  TDCWithinGroupSelfAdjustment_MG_HigherQ_F <- extract_relevant_data(model,
                                                                     "ChildLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:PrevSelfLatency > 
           ChildLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:PrevSelfLatency",
                                                                     "TDCWithinGroupSelfAdjustment_MG_HigherQ_F")
  
  ASDWithinGroupSelfAdjustment_MG_HigherQ_F <- extract_relevant_data(model,
                                                                     "ChildLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:PrevSelfLatency < 
           ChildLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:PrevSelfLatency",
                                                                     "ASDWithinGroupSelfAdjustment_MG_HigherQ_F")
  
  SelfAdjustmentComparisonsHypothesesFull <- rbind(
    SelfAdjustmentOverallHypotheses,
    SelfAdjustment_MG_F_ASDHigherTDC,
    SelfAdjustment_Q_F_ASDHigherTDC,
    SelfAdjustment_Q_UF_ASDHigherTDC,
    SelfAdjustment_overall_ASDHigherTDC,
    DifferenceSelfAdjustment_ASDLowerMatchingGameTDCVersusQuestions,
    DifferenceSelfAdjustment_ASDHigherFamiliarityTDCVersusUnfamiliarity,
    TDCWithinGroupSelfAdjustment_Q_F_HigherQ_UF,
    ASDWithinGroupSelfAdjustment_Q_F_HigherQ_UF,
    TDCWithinGroupSelfAdjustment_MG_HigherQ_F,
    ASDWithinGroupSelfAdjustment_MG_HigherQ_F) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 1 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -1 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 1 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 1 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -1 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -1 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  
  return(SelfAdjustmentComparisonsHypothesesFull)
}

extract_all_hypothesis_tests_m4_interpersonal_adult <- function(model) {
  
  ASD_MG_F <- extract_relevant_data(model, 
                                    "AdultLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AdultPrevOtherLatency > 0",
                                    "ASD_MG_F")
  
  TDC_MG_F <- extract_relevant_data(model, 
                                    "AdultLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AdultPrevOtherLatency > 0",
                                    "TDC_MG_F")
  
  ASD_Q_F <- extract_relevant_data(model, 
                                   "AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AdultPrevOtherLatency > 0",
                                   "ASD_Q_F")
  
  TDC_Q_F <- extract_relevant_data(model, 
                                   "AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AdultPrevOtherLatency > 0",
                                   "TDC_Q_F")
  
  ASD_Q_UF <- extract_relevant_data(model, 
                                    "AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AdultPrevOtherLatency < 0",
                                    "ASD_Q_UF")
  
  TDC_Q_UF <- extract_relevant_data(model, 
                                    "AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AdultPrevOtherLatency > 0",
                                    "TDC_Q_UF")
  
  ASD_overall <- extract_relevant_data(model, 
                                       "((AdultLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AdultPrevOtherLatency +
           AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AdultPrevOtherLatency +
           AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AdultPrevOtherLatency)/3) > 0",
                                       "ASD_overall")
  
  TDC_overall <- extract_relevant_data(model, 
                                       "((AdultLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AdultPrevOtherLatency +
           AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AdultPrevOtherLatency +
           AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AdultPrevOtherLatency)/3) > 0",
                                       "TDC_overall")
  
  InterpersonalAdjustmentOverallHypotheses <- rbind(ASD_MG_F, ASD_Q_F, ASD_Q_UF, ASD_overall,
                                                    TDC_MG_F, TDC_Q_F, TDC_Q_UF, TDC_overall)
  
  InterpersonalAdjustment_MG_F_ASDHigherTDC <- extract_relevant_data(model,
                                                                     "AdultLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AdultPrevOtherLatency > 
           AdultLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AdultPrevOtherLatency",
                                                                     "InterpersonalAdjustment_MG_F_ASDHigherTDC")
  
  InterpersonalAdjustment_Q_F_ASDHigherTDC <- extract_relevant_data(model,
                                                                    "AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AdultPrevOtherLatency > 
           AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AdultPrevOtherLatency",
                                                                    "InterpersonalAdjustment_Q_F_ASDHigherTDC")
  
  InterpersonalAdjustment_Q_UF_ASDHigherTDC <- extract_relevant_data(model,
                                                                     "AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AdultPrevOtherLatency > 
           AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AdultPrevOtherLatency",
                                                                     "InterpersonalAdjustment_Q_UF_ASDHigherTDC")
  
  InterpersonalAdjustment_overall_ASDHigherTDC <- extract_relevant_data(model, 
                                                                        "((AdultLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AdultPrevOtherLatency +
                                                       AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AdultPrevOtherLatency +
                                                       AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AdultPrevOtherLatency) / 3) >
                                                       ((AdultLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AdultPrevOtherLatency +
                                                       AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AdultPrevOtherLatency +
                                                       AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AdultPrevOtherLatency) / 3)",
                                                                        "InterpersonalAdjustment_overall_ASDHigherTDC")
  
  DifferenceInterpersonalAdjustment_ASDLowerMatchingGameTDCVersusQuestions <- extract_relevant_data(model,
                                                                                                    "(AdultLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AdultPrevOtherLatency - AdultLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AdultPrevOtherLatency) <
            (AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AdultPrevOtherLatency - AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AdultPrevOtherLatency)",
                                                                                                    "DifferenceInterpersonalAdjustment_ASDLowerMatchingGameTDCVersusQuestions")
  
  DifferenceInterpersonalAdjustment_ASDHigherFamiliarityTDCVersusUnfamiliarity <- extract_relevant_data(model,
                                                                                                        "(AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AdultPrevOtherLatency - AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AdultPrevOtherLatency) >
            (AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AdultPrevOtherLatency - AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AdultPrevOtherLatency)",
                                                                                                        "DifferenceInterpersonalAdjustment_ASDHigherFamiliarityTDCVersusUnfamiliarity")
  
  TDCWithinGroupInterpersonalAdjustment_Q_F_HigherQ_UF <- extract_relevant_data(model,
                                                                                "AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AdultPrevOtherLatency > 
           AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AdultPrevOtherLatency",
                                                                                "TDCWithinGroupInterpersonalAdjustment_Q_F_HigherQ_UF")
  
  ASDWithinGroupInterpersonalAdjustment_Q_F_HigherQ_UF <- extract_relevant_data(model,
                                                                                "AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AdultPrevOtherLatency > 
           AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AdultPrevOtherLatency",
                                                                                "ASDWithinGroupInterpersonalAdjustment_Q_F_HigherQ_UF")
  
  TDCWithinGroupInterpersonalAdjustment_MG_HigherQ_F <- extract_relevant_data(model,
                                                                              "AdultLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AdultPrevOtherLatency > 
           AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AdultPrevOtherLatency",
                                                                              "TDCWithinGroupInterpersonalAdjustment_MG_HigherQ_F")
  
  ASDWithinGroupInterpersonalAdjustment_MG_HigherQ_F <- extract_relevant_data(model,
                                                                              "AdultLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AdultPrevOtherLatency < 
           AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AdultPrevOtherLatency",
                                                                              "ASDWithinGroupInterpersonalAdjustment_MG_HigherQ_F")
  
  
  InterpersonalComparisonsHypothesesFull <- rbind(
    InterpersonalAdjustmentOverallHypotheses,
    InterpersonalAdjustment_MG_F_ASDHigherTDC,
    InterpersonalAdjustment_Q_F_ASDHigherTDC,
    InterpersonalAdjustment_Q_UF_ASDHigherTDC,
    InterpersonalAdjustment_overall_ASDHigherTDC,
    DifferenceInterpersonalAdjustment_ASDLowerMatchingGameTDCVersusQuestions,
    DifferenceInterpersonalAdjustment_ASDHigherFamiliarityTDCVersusUnfamiliarity,
    TDCWithinGroupInterpersonalAdjustment_Q_F_HigherQ_UF,
    ASDWithinGroupInterpersonalAdjustment_Q_F_HigherQ_UF,
    TDCWithinGroupInterpersonalAdjustment_MG_HigherQ_F,
    ASDWithinGroupInterpersonalAdjustment_MG_HigherQ_F) %>%
    rename("value" = TestDescription) %>%
    #mutate(Estimate = round(Estimate, 1)) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  return(InterpersonalComparisonsHypothesesFull)
}

extract_all_hypothesis_tests_m4_self_adult <- function(model) {
  # Let's extract the same for SelfLatency:
  ASD_MG_F <- extract_relevant_data(model, 
                                    "AdultLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AdultPrevSelfLatency > 0",
                                    "ASD_MG_F")
  
  TDC_MG_F <- extract_relevant_data(model, 
                                    "AdultLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AdultPrevSelfLatency > 0",
                                    "TDC_MG_F")
  
  ASD_Q_F <- extract_relevant_data(model, 
                                   "AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AdultPrevSelfLatency > 0",
                                   "ASD_Q_F")
  
  TDC_Q_F <- extract_relevant_data(model, 
                                   "AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AdultPrevSelfLatency > 0",
                                   "TDC_Q_F")
  
  ASD_Q_UF <- extract_relevant_data(model, 
                                    "AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AdultPrevSelfLatency > 0",
                                    "ASD_Q_UF")
  
  TDC_Q_UF <- extract_relevant_data(model, 
                                    "AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AdultPrevSelfLatency < 0",
                                    "TDC_Q_UF")
  
  ASD_overall <- extract_relevant_data(model, 
                                       "((AdultLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AdultPrevSelfLatency +
           AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AdultPrevSelfLatency +
           AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AdultPrevSelfLatency)/3) > 0",
                                       "ASD_overall")
  
  TDC_overall <- extract_relevant_data(model, 
                                       "((AdultLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AdultPrevSelfLatency +
           AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AdultPrevSelfLatency +
           AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AdultPrevSelfLatency)/3) > 0",
                                       "TDC_overall")
  
  SelfAdjustmentOverallHypotheses <- rbind(ASD_MG_F, ASD_Q_F, ASD_Q_UF, ASD_overall,
                                           TDC_MG_F, TDC_Q_F, TDC_Q_UF, TDC_overall)
  
  SelfAdjustment_MG_F_ASDHigherTDC <- extract_relevant_data(model,
                                                            "AdultLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AdultPrevSelfLatency < 
           AdultLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AdultPrevSelfLatency",
                                                            "SelfAdjustment_MG_F_ASDHigherTDC")
  
  SelfAdjustment_Q_F_ASDHigherTDC <- extract_relevant_data(model,
                                                           "AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AdultPrevSelfLatency > 
           AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AdultPrevSelfLatency",
                                                           "SelfAdjustment_Q_F_ASDHigherTDC")
  
  SelfAdjustment_Q_UF_ASDHigherTDC <- extract_relevant_data(model,
                                                            "AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AdultPrevSelfLatency > 
           AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AdultPrevSelfLatency",
                                                            "SelfAdjustment_Q_UF_ASDHigherTDC")
  
  SelfAdjustment_overall_ASDHigherTDC <- extract_relevant_data(model, 
                                                               "((AdultLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AdultPrevSelfLatency +
                                                       AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AdultPrevSelfLatency +
                                                       AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AdultPrevSelfLatency) / 3) >
                                                       ((AdultLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AdultPrevSelfLatency +
                                                       AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AdultPrevSelfLatency +
                                                       AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AdultPrevSelfLatency) / 3)",
                                                               "SelfAdjustment_overall_ASDHigherTDC")
  
  DifferenceSelfAdjustment_ASDLowerMatchingGameTDCVersusQuestions <- extract_relevant_data(model,
                                                                                           "(AdultLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AdultPrevSelfLatency - AdultLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AdultPrevSelfLatency) <
            (AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AdultPrevSelfLatency - AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AdultPrevSelfLatency)",
                                                                                           "DifferenceSelfAdjustment_ASDLowerMatchingGameTDCVersusQuestions")
  
  DifferenceSelfAdjustment_ASDHigherFamiliarityTDCVersusUnfamiliarity <- extract_relevant_data(model,
                                                                                               "(AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AdultPrevSelfLatency - AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AdultPrevSelfLatency) <
            (AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AdultPrevSelfLatency - AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AdultPrevSelfLatency)",
                                                                                               "DifferenceSelfAdjustment_ASDHigherFamiliarityTDCVersusUnfamiliarity")
  
  
  TDCWithinGroupSelfAdjustment_Q_F_HigherQ_UF <- extract_relevant_data(model,
                                                                       "AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AdultPrevSelfLatency > 
           AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AdultPrevSelfLatency",
                                                                       "TDCWithinGroupSelfAdjustment_Q_F_HigherQ_UF")
  
  ASDWithinGroupSelfAdjustment_Q_F_HigherQ_UF <- extract_relevant_data(model,
                                                                       "AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AdultPrevSelfLatency > 
           AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AdultPrevSelfLatency",
                                                                       "ASDWithinGroupSelfAdjustment_Q_F_HigherQ_UF")
  
  TDCWithinGroupSelfAdjustment_MG_HigherQ_F <- extract_relevant_data(model,
                                                                     "AdultLatency_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AdultPrevSelfLatency > 
           AdultLatency_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AdultPrevSelfLatency",
                                                                     "TDCWithinGroupSelfAdjustment_MG_HigherQ_F")
  
  ASDWithinGroupSelfAdjustment_MG_HigherQ_F <- extract_relevant_data(model,
                                                                     "AdultLatency_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AdultPrevSelfLatency < 
           AdultLatency_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AdultPrevSelfLatency",
                                                                     "ASDWithinGroupSelfAdjustment_MG_HigherQ_F")
  
  SelfAdjustmentComparisonsHypothesesFull <- rbind(
    SelfAdjustmentOverallHypotheses,
    SelfAdjustment_MG_F_ASDHigherTDC,
    SelfAdjustment_Q_F_ASDHigherTDC,
    SelfAdjustment_Q_UF_ASDHigherTDC,
    SelfAdjustment_overall_ASDHigherTDC,
    DifferenceSelfAdjustment_ASDLowerMatchingGameTDCVersusQuestions,
    DifferenceSelfAdjustment_ASDHigherFamiliarityTDCVersusUnfamiliarity,
    TDCWithinGroupSelfAdjustment_Q_F_HigherQ_UF,
    ASDWithinGroupSelfAdjustment_Q_F_HigherQ_UF,
    TDCWithinGroupSelfAdjustment_MG_HigherQ_F,
    ASDWithinGroupSelfAdjustment_MG_HigherQ_F) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  return(SelfAdjustmentComparisonsHypothesesFull)
}

extract_all_hypothesis_tests_m1_no_child <- function(model) {
  
  # Overall Estimates in Autism Group
  OverallLatencyThreeConditionsASD <- extract_relevant_data(model,
                                                            "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar +DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
                                                            "OverallLatencyThreeConditionsASD")
  
  MGOverallLatencyASD <- extract_relevant_data(model,
                                               "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > 0",
                                               "MGOverallLatencyASD>0")
  
  Q_F_OverallLatencyASD <- extract_relevant_data(model,
                                                 "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > 0",
                                                 "Q_F_OverallLatencyASD>0")
  
  Q_UF_OverallLatencyASD <- extract_relevant_data(model,
                                                  "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                  "Q_UF_OverallLatencyASD>0")
  
  # Overall Estimates in TD Group
  OverallLatencyThreeConditionsTD <- extract_relevant_data(model,
                                                           "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar +DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3 > 0",
                                                           "OverallLatencyThreeConditionsTD")
  
  MGOverallLatencyTDC <- extract_relevant_data(model,
                                               "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > 0",
                                               "MGOverallLatencyTDC>0")
  
  Q_F_OverallLatencyTDC <- extract_relevant_data(model,
                                                 "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar > 0",
                                                 "Q_F_OverallLatencyTDC>0")
  
  Q_UF_OverallLatencyTDC <- extract_relevant_data(model,
                                                  "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                  "Q_UF_OverallLatencyTDC>0")
  
  # Overall Differences between Autism Group and Typical Development Group
  OverallLatencyThreeConditionsASDFasterTD <- extract_relevant_data(model,
                                                                    "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar +DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 
                                <
                                (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar +DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3",
                                                                    "OverallLatencyThreeConditionsASD<TD")
  
  
  # Overall values for Beta parameter
  BetaParameterASD_MG <- extract_relevant_nonlatency_data(model,
                                                          "beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                          "BetaParameterASD_MG")
  
  BetaParameterASD_Q_F <- extract_relevant_nonlatency_data(model,
                                                           "beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < 0",
                                                           "BetaParameterASD_Q_F")
  
  BetaParameterASD_Q_UF <- extract_relevant_nonlatency_data(model,
                                                            "beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < 0",
                                                            "BetaParameterASD_Q_UF")
  
  BetaParameterASD_Overall <- extract_relevant_nonlatency_data(model,
                                                               "(beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar +
                 beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar +
                 beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar) / 3 < 0",
                                                               "BetaParameterASD_Overall")
  
  # Differences, Typical Development Group with Long Pauses
  BetaParameterTDC_MG <- extract_relevant_nonlatency_data(model,
                                                          "beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                          "BetaParameterTDC_MG")
  
  BetaParameterTDC_Q_F <- extract_relevant_nonlatency_data(model,
                                                           "beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar < 0",
                                                           "BetaParameterTDC_Q_F")
  
  BetaParameterTDC_Q_UF <- extract_relevant_nonlatency_data(model,
                                                            "beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar < 0",
                                                            "BetaParameterTDC_Q_UF")
  
  BetaParameterTDC_Overall <- extract_relevant_nonlatency_data(model,
                                                               "(beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar +
                 beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar +
                 beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3 < 0",
                                                               "BetaParameterTDC_Overall")
  
  
  
  
  WithinGroupBetaParameterASD_MG_Slower_Q_F <- extract_relevant_nonlatency_data(model,
                                                                                "beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                                "WithinGroupBetaParameterASD_MG_Slower_Q_F")
  
  WithinGroupBetaParameterASD_Q_F_Slower_Q_UF <- extract_relevant_nonlatency_data(model,
                                                                                  "beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                                  "WithinGroupBetaParameterASD_Q_F_Slower_Q_UF")
  
  WithinGroupBetaParameterTDC_MG_Slower_Q_F <- extract_relevant_nonlatency_data(model,
                                                                                "beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                                "WithinGroupBetaParameterTDC_MG_Slower_Q_F")
  
  WithinGroupBetaParameterTDC_Q_F_Slower_Q_UF <- extract_relevant_nonlatency_data(model,
                                                                                  "beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                                  "WithinGroupBetaParameterTDC_Q_F_Slower_Q_UF")
  
  DifferenceInDifferenceBetaParameterBetweenGroupsMGBiggerQ_F <- extract_relevant_nonlatency_data(model,
                                                                                                  "(beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar - beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) > 
                 (beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar - beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar)",
                                                                                                  "DifferenceInDifferenceBetaParameterBetweenGroupsMGBiggerQ_F")
  
  DifferenceInDifferenceBetaParameterBetweenGroupsFamiliarityBiggerUnfamiliarity <- extract_relevant_nonlatency_data(model,
                                                                                                                     "(beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar - beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar) < 
                 (beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar - beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)",
                                                                                                                     "DifferenceInDifferenceBetaParameterBetweenGroupsFamiliarityBiggerUnfamiliarity")
  
  # Differences between Autism Group and Typical Development Group with Long Pauses
  OverallBetaParameterASDLowerTD <- extract_relevant_nonlatency_data(model,
                                                                     "(beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar + beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar) / 3 < ((beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar + beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3)",
                                                                     "OverallBetaParameterASD<TD")
  
  BetaParameterASDFasterTD_MG <- extract_relevant_nonlatency_data(model,
                                                                  "beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar < beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar",
                                                                  "BetaParameter_ASD_MG<TD_MG")
  
  BetaParameterASDFasterTD_Q_F <- extract_relevant_nonlatency_data(model,
                                                                   "beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                   "BetaParameter_ASD_MG<TD_Q_F")
  
  BetaParameterASDFasterTD_Q_UF <- extract_relevant_nonlatency_data(model,
                                                                    "beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                    "BetaParameter_ASD_MG<TD_Q_UF")
  
  # Differences within Autism Group across Conditions:
  WithinASDGroupMGSlowerQuestionsFamiliar <- extract_relevant_data(model,
                                                                   "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                   "WithinASDGroupMG>QuestionsFamiliar")
  
  # Differences within TD Group across Conditions:
  WithinTDGroupMGSlowerQuestionsFamiliar <- extract_relevant_data(model,
                                                                  "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                  "WithinTDGroupMG>QuestionsFamiliar")
  
  # Differences within Autism Group across Familiarity:
  WithinAutismGroupQuestionsUnfamiliarSlowerFamiliar <- extract_relevant_data(model,
                                                                              "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                              "WithinAutismGroupQuestionsUnfamiliar>Familiar")
  
  # Differences within TD Group across Familiarity:
  WithinTDGroupQuestionsUnfamiliarSlowerFamiliar <- extract_relevant_data(model,
                                                                          "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                          "WithinTDGroupQuestionsUnfamiliar>Familiar")
  
  # Differences between Autism Group and Typical Development Group Within MatchingGame
  WithinConditionsMatchingGameASDFasterTD <- extract_relevant_data(model,
                                                                   "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar < DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar",
                                                                   "WithinConditionsMatchingGameASD<TD")
  
  # Differences between Autism Group and Typical Development Group Within Familiar Questions:
  WithinConditionsQuestionsFamiliarASDFasterTD <- extract_relevant_data(model,
                                                                        "DiagnosisASD:TaskQuestions:FamiliarityFamiliar < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                        "WithinConditionsQuestionsFamiliarASD<TD")
  
  # Differences between Autism Group and Typical Development Group Within Unfamiliar Questions:
  WithinConditionsQuestionsUnfamiliarASDFasterTD <- extract_relevant_data(model,
                                                                          "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                          "WithinConditionsQuestionsUnfamiliarASD<TD")
  
  # Differences between Autism Group and Typical Development Group in Differences within Familiarity
  DifferenceInDifferenceBetweenGroupsUnfamiliarBiggerFamiliarity <- extract_relevant_data(model,
                                                                                          "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) < 
                      (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                          "DifferenceInDifferenceBetweenGroupsUnfamiliar>Familiarity")
  
  # Differences between Autism Group and Typical Development Group in Differences within Conditions
  DifferenceInDifferenceBetweenGroupsMGBiggerThanConversations <- extract_relevant_data(model,
                                                                                        "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar) > 
                      (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                        "DifferenceInDifferenceBetweenGroupsMG>Conversations")
  
  # Overall Subject Variability in Autismm Group
  SubjectVariabilityASD <- extract_relevant_sd_data(model, 
                                                    "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisASD + TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD)/3 > 0",
                                                    "SubjectVariabilityASD")
  
  SubjectVariabilityASD_MG <- extract_relevant_sd_data(model, 
                                                       "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD > 0",
                                                       "SubjectVariabilityASD_MG")
  
  SubjectVariabilityASD_Q_F <- extract_relevant_sd_data(model, 
                                                        "TaskQuestions:FamiliarityFamiliar:DiagnosisASD > 0",
                                                        "SubjectVariabilityASD_Q_F")
  
  SubjectVariabilityASD_Q_UF <- extract_relevant_sd_data(model, 
                                                         "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD > 0",
                                                         "SubjectVariabilityASD_Q_UF")
  
  # Overall Subject Variability in TD Group
  SubjectVariabilityTD <- extract_relevant_sd_data(model, 
                                                   "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisTDC + TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC)/3 > 0",
                                                   "SubjectVariabilityTD")
  
  SubjectVariabilityTD_MG <- extract_relevant_sd_data(model, 
                                                      "TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC > 0",
                                                      "SubjectVariabilityTD_MG")
  
  SubjectVariabilityTD_Q_F <- extract_relevant_sd_data(model, 
                                                       "TaskQuestions:FamiliarityFamiliar:DiagnosisTDC > 0",
                                                       "SubjectVariabilityTD_Q_F")
  
  SubjectVariabilityTD_Q_UF <- extract_relevant_sd_data(model, 
                                                        "TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC > 0",
                                                        "SubjectVariabilityTD_Q_UF")
  
  # Overall Subject Variability in TD Group
  SubjectVariabilityASDSmallerThanTD <- extract_relevant_sd_data(model,
                                                                 "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisASD + TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD) / 3 <
                                              (TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisTDC + TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC) / 3",
                                                                 "SubjectVariabilityASD<TD")
  
  
  SubjectVariabilityASDSmallerThanTD_MG <- extract_relevant_sd_data(model,
                                                                    "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD  <
                                              TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC",
                                                                    "SubjectVariabilityASD_MG<TD_MG")
  
  SubjectVariabilityASDSmallerThanTD_Q_F <- extract_relevant_sd_data(model,
                                                                     "TaskQuestions:FamiliarityFamiliar:DiagnosisASD  <
                                              TaskQuestions:FamiliarityFamiliar:DiagnosisTDC",
                                                                     "SubjectVariabilityASD_Q_F<TD_Q_F")
  
  SubjectVariabilityASDSmallerThanTD_Q_UF <- extract_relevant_sd_data(model,
                                                                      "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD  <
                                              TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC",
                                                                      "SubjectVariabilityASD_Q_UF<TD_Q_UF")
  
  SubjectVariabilityASD_MG_SmallerThan_Q_F <- extract_relevant_sd_data(model,
                                                                       "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD > TaskQuestions:FamiliarityFamiliar:DiagnosisASD",
                                                                       "SubjectVariabilityASD_MG_SmallerThan_Q_F")
  
  SubjectVariabilityASD_Q_F_SmallerThan_Q_UF <- extract_relevant_sd_data(model,
                                                                         "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD > TaskQuestions:FamiliarityFamiliar:DiagnosisASD",
                                                                         "SubjectVariabilityASD_Q_F_SmallerThan_Q_UF")
  
  SubjectVariabilityTDC_MG_SmallerThan_Q_F <- extract_relevant_sd_data(model,
                                                                       "TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC < TaskQuestions:FamiliarityFamiliar:DiagnosisTDC",
                                                                       "SubjectVariabilityTDC_MG_SmallerThan_Q_F")
  
  SubjectVariabilityTDC_Q_F_SmallerThan_Q_UF <- extract_relevant_sd_data(model,
                                                                         "TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC > TaskQuestions:FamiliarityFamiliar:DiagnosisTDC",
                                                                         "SubjectVariabilityTDC_Q_F_SmallerThan_Q_UF")
  
  DifferenceInDifferenceSubjectVariabilityUnfamiliarSmallerFamiliarity <- extract_relevant_sd_data(model,
                                                                                                   "(TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD - TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC) 
                                              < 
                                              (TaskQuestions:FamiliarityFamiliar:DiagnosisASD - TaskQuestions:FamiliarityFamiliar:DiagnosisTDC)",
                                                                                                   "DifferenceInDifferenceSubjectVariabilityUnfamiliarBiggerFamiliarity")
  
  DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions <- extract_relevant_sd_data(model,
                                                                                                  "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC - TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD) 
                                              < 
                                              (TaskQuestions:FamiliarityFamiliar:DiagnosisTDC - TaskQuestions:FamiliarityFamiliar:DiagnosisASD)",
                                                                                                  "DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions")
  
  # Overall Visit Variability across visits in Autism Group
  VisitVariabilityASDOverall <- extract_relevant_sdvisit_data(model, 
                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
                                                              "VisitVariabilityASDOverall")
  
  # Overall Visit Variability across visits in TD Group
  VisitVariabilityTDCOverall <- extract_relevant_sdvisit_data(model, 
                                                              "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
                                                              "VisitVariabilityTDCOverall")
  
  
  VisitVariabilityASDOverallHigherTDSOverall <- extract_relevant_sdvisit_data(model, 
                                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 >
                                               (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3",
                                                                              "VisitVariabilityASDOverallLowerTDSOverall")
  
  # Subject Variability across visits MG in Autism Group
  VisitVariabilityASDMG <- extract_relevant_sdvisit_data(model, 
                                                         "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                         "VisitVariabilityASDMG")
  
  # Overall Subject Variability across visits in TD Group
  VisitVariabilityTDMG <- extract_relevant_sdvisit_data(model, 
                                                        "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                        "VisitVariabilityTDMG")
  
  # Subject Variability across visits MG in Autism Group
  VisitVariabilityASD_Q_F <- extract_relevant_sdvisit_data(model, 
                                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > 0",
                                                           "VisitVariabilityASD_Q_F")
  
  # Overall Subject Variability across visits in TD Group
  VisitVariabilityTDC_Q_F <- extract_relevant_sdvisit_data(model, 
                                                           "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar > 0",
                                                           "VisitVariabilityTDC_Q_F")
  
  # Subject Variability across visits MG in Autism Group
  VisitVariabilityASD_Q_UF <- extract_relevant_sdvisit_data(model, 
                                                            "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                            "VisitVariabilityASD_Q_UF")
  
  # Overall Subject Variability across visits in TD Group
  VisitVariabilityTDC_Q_UF <- extract_relevant_sdvisit_data(model, 
                                                            "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                            "VisitVariabilityTDC_Q_UF")
  
  # Overall Subject Variability across visits, differences
  VisitVariabilityASDMGSlowerThanTDMG <- extract_relevant_sdvisit_data(model, 
                                                                       "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar < DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar",
                                                                       "VisitVariabilityASDMG>TDMG")
  
  # Overall Subject Variability across visits, differences
  VisitVariabilityASD_Q_F_SlowerThanTD_Q_F <- extract_relevant_sdvisit_data(model, 
                                                                            "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                            "VisitVariabilityASD_Q_F_SlowerThanTD_Q_F")
  
  # Overall Subject Variability across visits, differences
  VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF <- extract_relevant_sdvisit_data(model, 
                                                                              "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                              "VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF")
  
  #DifferenceInDifference:
  DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity <- extract_relevant_sdvisit_data(model, 
                                                                                                       "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)
                                               < (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                                       "DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity")
  
  DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG <- extract_relevant_sdvisit_data(model, 
                                                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar)
                                               > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                              "DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG")
  
  
  # Overall Sigma Variability across visits in Autism Group
  SigmaASDOverall <- extract_relevant_sigma_data(model, 
                                                 "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 < 0",
                                                 "SigmaASDOverall")
  
  # Overall Sigma Variability across visits in TD Group
  SigmaTDCOverall <- extract_relevant_sigma_data(model, 
                                                 "(sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) /3 < 0",
                                                 "SigmaTDCOverall")
  
  
  SigmaASDOverallHigherTDSOverall <- extract_relevant_nonlatency_data(model, 
                                                                      "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 >
                                               (sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3",
                                                                      "SigmaASDOverallHigherTDSOverall")
  
  # Subject Variability across visits MG in Autism Group
  SigmaASDMG <- extract_relevant_sigma_data(model, 
                                            "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar < 0",
                                            "SigmaASDMG")
  
  # Overall Subject Variability across visits in TD Group
  SigmaTDMG <- extract_relevant_sigma_data(model, 
                                           "sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar < 0",
                                           "SigmaTDMG")
  
  # Subject Variability across visits MG in Autism Group
  SigmaASD_Q_F <- extract_relevant_sigma_data(model, 
                                              "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < 0",
                                              "SigmaASD_Q_F")
  
  # Overall Subject Variability across visits in TD Group
  SigmaTDC_Q_F <- extract_relevant_sigma_data(model, 
                                              "sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar < 0",
                                              "SigmaTDC_Q_F")
  
  
  # Subject Variability across visits MG in Autism Group
  SigmaASD_Q_UF <- extract_relevant_sigma_data(model, 
                                               "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < 0 ",
                                               "SigmaASD_Q_UF")
  
  # Overall Subject Variability across visits in TD Group
  SigmaTDC_Q_UF <- extract_relevant_sigma_data(model, 
                                               "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar < 0",
                                               "SigmaTDC_Q_UF")
  
  # Overall Subject Variability across visits, differences
  SigmaASDMGSlowerThanTDMG <- extract_relevant_nonlatency_data(model, 
                                                               "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar",
                                                               "SigmaASDMGSlowerThanTDMG")
  
  # Overall Subject Variability across visits, differences
  SigmaASD_Q_F_SlowerThanTD_Q_F <- extract_relevant_nonlatency_data(model, 
                                                                    "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                    "SigmaASD_Q_F_SlowerThanTD_Q_F")
  
  # Overall Subject Variability across visits, differences
  SigmaASD_Q_UF_SlowerThanTD_Q_UF <- extract_relevant_nonlatency_data(model, 
                                                                      "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                      "SigmaASD_Q_UF_SlowerThanTD_Q_UF")
  
  # Sigma across visits, differences
  # SigmaASD_MG_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit < 0",
  #                                                "SigmaASD_MG_VisitLowerZero")
  
  # SigmaTDC_MG_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit > 0",
  #                                                "SigmaTDC_MG_VisitLowerZero")
  # 
  # SigmaASD_MG_VisitLowerTDC <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit < sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
  #                                                "SigmaASD_MG_VisitLowerTDC")
  # 
  # SigmaASD_Overall_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
  #                                                "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit + 
  #                                                sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit +
  #                                                sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit) / 3 < 0",
  #                                                "SigmaASD_Overall_VisitLowerZero")
  # 
  # SigmaTDC_Overall_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
  #                                                "(sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit + 
  #                                                sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit +
  #                                                sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit) / 3 < 0",
  #                                                "SigmaTDC_Overall_VisitLowerZero")
  # 
  # # Sigma across visits, differences
  # SigmaASD_Q_F_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit > 0",
  #                                                "SigmaASD_Q_F_VisitLowerZero")
  # 
  # SigmaTDC_Q_F_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit < 0",
  #                                                "SigmaTDC_Q_F_VisitLowerZero")
  # 
  # SigmaASD_Q_F_VisitLowerTDC <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit > sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
  #                                                "SigmaASD_Q_F_VisitLowerTDC")
  # 
  # # Sigma across visits, differences
  # SigmaASD_Q_UF_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < 0",
  #                                                "SigmaASD_Q_UF_VisitLowerZero")
  # 
  # SigmaTDC_Q_UF_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > 0",
  #                                                "SigmaTDC_Q_UF_VisitLowerZero")
  # 
  # SigmaASD_Q_UF_VisitLowerTDC <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
  #                                                "SigmaASD_Q_UF_VisitLowerTDC")
  # 
  # SigmaVisitOverall_ASDLowerTDC <- extract_relevant_nonlatency_data(model, 
  #                                                "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit + 
  #                                                sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit +
  #                                                sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit) / 3 > (sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit + 
  #                                                sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit +
  #                                                sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit) /3",
  #                                                "SigmaVisitOverall_ASDLowerTDC")
  
  #DifferenceInDifference:
  DifferenceInDifferenceSigmaASDFamiliarityBiggerTDUnfamiliarity <- extract_relevant_nonlatency_data(model, 
                                                                                                     "(sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)
                                               > (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                                     "DifferenceInDifferenceSigmaASDFamiliarityBiggerTDUnfamiliarity")
  
  DifferenceInDifferenceSigmaASDMGSmallerThanTDMG <- extract_relevant_nonlatency_data(model, 
                                                                                      "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar)
                                               > (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                      "DifferenceInDifferenceSigmaASDMGSmallerThanTDMG")
  
  #WithinGroups
  WithinGroupSigmaASD_Q_F_Faster_UF <- extract_relevant_nonlatency_data(model, 
                                                                        "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaASD_Q_FFaster_UF")
  
  WithinGroupSigmaTDC_Q_F_Faster_UF <- extract_relevant_nonlatency_data(model, 
                                                                        "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaTDC_Q_F_Faster_UF")
  
  WithinGroupSigmaASD_Q_F_Faster_MG <- extract_relevant_nonlatency_data(model, 
                                                                        "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaASD_Q_F_Faster_MG")
  
  WithinGroupSigmaTDC_Q_F_Faster_MG <- extract_relevant_nonlatency_data(model, 
                                                                        "sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar < sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaTDC_Q_F_Faster_MG")
  
  
  # WithinGroupSigmaASD_Q_F_Lower_UF_Visit <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit",
  #                                                "WithinGroupSigmaASD_Q_F_Lower_UF_Visit")
  # 
  # WithinGroupSigmaTDC_Q_F_Lower_UF_Visit <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit",
  #                                                "WithinGroupSigmaTDC_Q_F_Lower_UF_Visit")
  # 
  # WithinGroupSigmaASD_Q_F_Lower_MG_Visit <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit > sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit",
  #                                                "WithinGroupSigmaASD_Q_F_Lower_MG_Visit")
  # 
  # WithinGroupSigmaTDC_Q_F_Lower_MG_Visit <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
  #                                                "WithinGroupSigmaTDC_Q_F_Lower_MG_Visit")
  # 
  # 
  # SigmaTDC_Q_UF_VisitLowerZero <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > 0",
  #                                                "SigmaTDC_Q_UF_VisitLowerZero")
  # 
  # SigmaASD_Q_UF_VisitLowerTDC <- extract_relevant_nonlatency_data(model, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
  #                                                "SigmaASD_Q_UF_VisitLowerTDC")
  # 
  # 
  # DifferenceInDifferenceSigmaVisit_Q_F_Lower_Q_UF <- extract_relevant_nonlatency_data(model, 
  #                                                "(sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit - sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit) <
  #                                                (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit)", "DifferenceInDifferenceSigmaVisit_Q_F_Lower_Q_UF")
  # 
  # 
  # DifferenceInDifferenceSigmaVisit_MG_Lower_Q_F <- extract_relevant_nonlatency_data(model, 
  #                                                "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit - sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit) <
  #                                                (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit)", "DifferenceInDifferenceSigmaVisit_MG_Lower_Q_F")
  
  
  
  # Correlation Estimates across visits in Autism Group
  post <- as_draws_df(model)
  
  
  # Apply the function for each case
  cor_post_mu_MG_ASD <- Visitextract_and_rename(post, "DiagnosisASD", "TaskMatchingGame", "FamiliarityFamiliar") %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_MG_TDC <- Visitextract_and_rename(post, "DiagnosisTDC", "TaskMatchingGame", "FamiliarityFamiliar") %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_Q_ASD_Familiar <- Visitextract_and_rename(post, "DiagnosisASD", "TaskQuestions", "FamiliarityFamiliar") %>%
    select(Visit2, Visit4, Visit6) %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_Q_TDC_Familiar <- Visitextract_and_rename(post, "DiagnosisTDC", "TaskQuestions", "FamiliarityFamiliar") %>%
    select(Visit2, Visit4, Visit6) %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_Q_ASD_Unfamiliar <- Visitextract_and_rename(post, "DiagnosisASD", "TaskQuestions", "FamiliarityUnfamiliar") %>%
    select(Visit1, Visit3, Visit5, Visit7) %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_Q_TDC_Unfamiliar <- Visitextract_and_rename(post, "DiagnosisTDC", "TaskQuestions", "FamiliarityUnfamiliar") %>%
    select(Visit1, Visit3, Visit5, Visit7) %>%
    cor() %>%
    as.matrix()
  
  diag(cor_post_mu_MG_ASD) <- NA
  diag(cor_post_mu_MG_TDC) <- NA
  diag(cor_post_mu_Q_ASD_Familiar) <- NA
  diag(cor_post_mu_Q_TDC_Familiar) <- NA
  diag(cor_post_mu_Q_ASD_Unfamiliar) <- NA
  diag(cor_post_mu_Q_TDC_Unfamiliar) <- NA
  
  CorrelationSummaryStatistics <- tibble(
    Estimate = c(mean(cor_post_mu_MG_ASD, na.rm = TRUE), mean(cor_post_mu_Q_ASD_Familiar, na.rm = TRUE), mean(cor_post_mu_Q_ASD_Unfamiliar, na.rm = TRUE), 
                 mean(cor_post_mu_MG_TDC, na.rm = TRUE), mean(cor_post_mu_Q_TDC_Familiar, na.rm = TRUE), mean(cor_post_mu_Q_TDC_Unfamiliar, na.rm = TRUE)),
    CI.Lower = c(quantile(cor_post_mu_MG_ASD, probs = 0.025, na.rm = TRUE), quantile(cor_post_mu_Q_ASD_Familiar, probs = 0.025, na.rm = TRUE), quantile(cor_post_mu_Q_ASD_Unfamiliar, probs = 0.025, na.rm = TRUE), 
                 quantile(cor_post_mu_MG_TDC, probs = 0.025, na.rm = TRUE), quantile(cor_post_mu_Q_TDC_Familiar, probs = 0.025, na.rm = TRUE), quantile(cor_post_mu_Q_TDC_Unfamiliar, probs = 0.025, na.rm = TRUE)),
    CI.Upper = c(quantile(cor_post_mu_MG_ASD, probs = 0.975, na.rm = TRUE), quantile(cor_post_mu_Q_ASD_Familiar, probs = 0.975, na.rm = TRUE), quantile(cor_post_mu_Q_ASD_Unfamiliar, probs = 0.975, na.rm = TRUE), 
                 quantile(cor_post_mu_MG_TDC, probs = 0.975, na.rm = TRUE), quantile(cor_post_mu_Q_TDC_Familiar, probs = 0.975, na.rm = TRUE), quantile(cor_post_mu_Q_TDC_Unfamiliar, probs = 0.975, na.rm = TRUE)),
    Evid.Ratio = c(NA, NA, NA, NA, NA, NA),
    TestDescription = c(
      "CorrelationAcrossVisitsASDMG", "CorrelationAcrossVisitsASDFamiliarQuestions", "CorrelationAcrossVisitsASDUnfamiliarQuestions",
      "CorrelationAcrossVisitsTDMG", "CorrelationAcrossVisitsTDFamiliarQuestions", "CorrelationAcrossVisitsTDUnfamiliarQuestions"),
  )
  
  HypothesisResultsChildLatenciesModel1 <- rbind(OverallLatencyThreeConditionsASD,
                                                 OverallLatencyThreeConditionsTD,
                                                 MGOverallLatencyASD,
                                                 Q_F_OverallLatencyASD,
                                                 Q_UF_OverallLatencyASD,
                                                 MGOverallLatencyTDC,
                                                 Q_F_OverallLatencyTDC,
                                                 Q_UF_OverallLatencyTDC,
                                                 BetaParameterASD_MG,
                                                 BetaParameterASD_Q_F,
                                                 BetaParameterASD_Q_UF,
                                                 BetaParameterASD_Overall,
                                                 BetaParameterTDC_MG,
                                                 BetaParameterTDC_Q_F,
                                                 BetaParameterTDC_Q_UF,
                                                 BetaParameterTDC_Overall,
                                                 OverallBetaParameterASDLowerTD,
                                                 BetaParameterASDFasterTD_MG,
                                                 BetaParameterASDFasterTD_Q_F,
                                                 BetaParameterASDFasterTD_Q_UF,
                                                 WithinGroupBetaParameterASD_MG_Slower_Q_F,
                                                 WithinGroupBetaParameterTDC_MG_Slower_Q_F,
                                                 DifferenceInDifferenceBetaParameterBetweenGroupsMGBiggerQ_F,
                                                 WithinGroupBetaParameterASD_Q_F_Slower_Q_UF,
                                                 WithinGroupBetaParameterTDC_Q_F_Slower_Q_UF,
                                                 DifferenceInDifferenceBetaParameterBetweenGroupsFamiliarityBiggerUnfamiliarity,
                                                 OverallLatencyThreeConditionsASDFasterTD,
                                                 WithinASDGroupMGSlowerQuestionsFamiliar,
                                                 WithinTDGroupMGSlowerQuestionsFamiliar,
                                                 WithinAutismGroupQuestionsUnfamiliarSlowerFamiliar,
                                                 WithinTDGroupQuestionsUnfamiliarSlowerFamiliar,
                                                 WithinConditionsMatchingGameASDFasterTD,
                                                 WithinConditionsQuestionsFamiliarASDFasterTD,
                                                 WithinConditionsQuestionsUnfamiliarASDFasterTD,
                                                 DifferenceInDifferenceBetweenGroupsUnfamiliarBiggerFamiliarity,
                                                 DifferenceInDifferenceBetweenGroupsMGBiggerThanConversations,
                                                 SubjectVariabilityASD,
                                                 SubjectVariabilityTD,
                                                 SubjectVariabilityASDSmallerThanTD,
                                                 SubjectVariabilityASD_MG,
                                                 SubjectVariabilityASD_Q_F,
                                                 SubjectVariabilityASD_Q_UF,
                                                 SubjectVariabilityTD_MG,
                                                 SubjectVariabilityTD_Q_F,
                                                 SubjectVariabilityTD_Q_UF,
                                                 SubjectVariabilityASDSmallerThanTD_MG,
                                                 SubjectVariabilityASDSmallerThanTD_Q_F,
                                                 SubjectVariabilityASDSmallerThanTD_Q_UF,
                                                 SubjectVariabilityASD_MG_SmallerThan_Q_F,
                                                 SubjectVariabilityASD_Q_F_SmallerThan_Q_UF,
                                                 SubjectVariabilityTDC_MG_SmallerThan_Q_F,
                                                 SubjectVariabilityTDC_Q_F_SmallerThan_Q_UF,
                                                 DifferenceInDifferenceSubjectVariabilityUnfamiliarSmallerFamiliarity,
                                                 DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions,
                                                 VisitVariabilityASDMG,
                                                 VisitVariabilityTDMG,
                                                 VisitVariabilityASD_Q_F,
                                                 VisitVariabilityTDC_Q_F,
                                                 VisitVariabilityASD_Q_UF,
                                                 VisitVariabilityTDC_Q_UF,
                                                 VisitVariabilityASDOverall,
                                                 VisitVariabilityTDCOverall,
                                                 VisitVariabilityASDMGSlowerThanTDMG,
                                                 VisitVariabilityASD_Q_F_SlowerThanTD_Q_F,
                                                 VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF,
                                                 VisitVariabilityASDOverallHigherTDSOverall,
                                                 DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG,
                                                 DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity,
                                                 SigmaASDOverall,
                                                 SigmaTDCOverall,
                                                 SigmaASDOverallHigherTDSOverall,
                                                 SigmaASDMG,
                                                 SigmaTDMG,
                                                 SigmaASD_Q_F,
                                                 SigmaTDC_Q_F,
                                                 SigmaASD_Q_UF,
                                                 SigmaTDC_Q_UF,
                                                 SigmaASDMGSlowerThanTDMG,
                                                 SigmaASD_Q_F_SlowerThanTD_Q_F,
                                                 SigmaASD_Q_UF_SlowerThanTD_Q_UF,
                                                 SigmaASD_MG_VisitLowerZero,
                                                 SigmaTDC_MG_VisitLowerZero,
                                                 SigmaASD_MG_VisitLowerTDC,
                                                 SigmaASD_Q_F_VisitLowerZero,
                                                 SigmaTDC_Q_F_VisitLowerZero,
                                                 SigmaASD_Q_F_VisitLowerTDC,
                                                 SigmaASD_Q_UF_VisitLowerZero,
                                                 SigmaTDC_Q_UF_VisitLowerZero,
                                                 SigmaASD_Q_UF_VisitLowerTDC,
                                                 SigmaASD_Overall_VisitLowerZero,
                                                 SigmaTDC_Overall_VisitLowerZero,
                                                 WithinGroupSigmaASD_Q_F_Lower_UF_Visit,
                                                 WithinGroupSigmaTDC_Q_F_Lower_UF_Visit,
                                                 WithinGroupSigmaASD_Q_F_Lower_MG_Visit,
                                                 WithinGroupSigmaTDC_Q_F_Lower_MG_Visit,
                                                 WithinGroupSigmaASD_Q_F_Faster_UF,
                                                 WithinGroupSigmaTDC_Q_F_Faster_UF,
                                                 WithinGroupSigmaASD_Q_F_Faster_MG,
                                                 WithinGroupSigmaTDC_Q_F_Faster_MG,
                                                 SigmaVisitOverall_ASDLowerTDC,
                                                 DifferenceInDifferenceSigmaVisit_Q_F_Lower_Q_UF,
                                                 DifferenceInDifferenceSigmaVisit_MG_Lower_Q_F,
                                                 DifferenceInDifferenceSigmaASDFamiliarityBiggerTDUnfamiliarity,
                                                 DifferenceInDifferenceSigmaASDMGSmallerThanTDMG,
                                                 CorrelationSummaryStatistics) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  return(HypothesisResultsChildLatenciesModel1)
}

extract_all_hypothesis_tests_m1_no_adult <- function(model) {
  
  # Overall Estimates in Autism Group
  OverallLatencyThreeConditionsASD <- extract_relevant_data(model,
                                                            "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar +DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
                                                            "OverallLatencyThreeConditionsASD")
  
  MGOverallLatencyASD <- extract_relevant_data(model,
                                               "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > 0",
                                               "MGOverallLatencyASD>0")
  
  Q_F_OverallLatencyASD <- extract_relevant_data(model,
                                                 "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > 0",
                                                 "Q_F_OverallLatencyASD>0")
  
  Q_UF_OverallLatencyASD <- extract_relevant_data(model,
                                                  "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                  "Q_UF_OverallLatencyASD>0")
  
  # Overall Estimates in TD Group
  OverallLatencyThreeConditionsTD <- extract_relevant_data(model,
                                                           "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar +DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3 > 0",
                                                           "OverallLatencyThreeConditionsTD")
  
  MGOverallLatencyTDC <- extract_relevant_data(model,
                                               "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > 0",
                                               "MGOverallLatencyTDC>0")
  
  Q_F_OverallLatencyTDC <- extract_relevant_data(model,
                                                 "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar > 0",
                                                 "Q_F_OverallLatencyTDC>0")
  
  Q_UF_OverallLatencyTDC <- extract_relevant_data(model,
                                                  "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                  "Q_UF_OverallLatencyTDC>0")
  
  # Overall Differences between Autism Group and Typical Development Group
  OverallLatencyThreeConditionsASDFasterTD <- extract_relevant_data(model,
                                                                    "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar +DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 
                                <
                                (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar +DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3",
                                                                    "OverallLatencyThreeConditionsASD<TD")
  
  
  # Overall values for Beta parameter
  BetaParameterASD_MG <- extract_relevant_nonlatency_data(model,
                                                          "beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                          "BetaParameterASD_MG")
  
  BetaParameterASD_Q_F <- extract_relevant_nonlatency_data(model,
                                                           "beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < 0",
                                                           "BetaParameterASD_Q_F")
  
  BetaParameterASD_Q_UF <- extract_relevant_nonlatency_data(model,
                                                            "beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < 0",
                                                            "BetaParameterASD_Q_UF")
  
  BetaParameterASD_Overall <- extract_relevant_nonlatency_data(model,
                                                               "(beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar +
                 beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar +
                 beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar) / 3 < 0",
                                                               "BetaParameterASD_Overall")
  
  # Differences, Typical Development Group with Long Pauses
  BetaParameterTDC_MG <- extract_relevant_nonlatency_data(model,
                                                          "beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                          "BetaParameterTDC_MG")
  
  BetaParameterTDC_Q_F <- extract_relevant_nonlatency_data(model,
                                                           "beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar < 0",
                                                           "BetaParameterTDC_Q_F")
  
  BetaParameterTDC_Q_UF <- extract_relevant_nonlatency_data(model,
                                                            "beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar < 0",
                                                            "BetaParameterTDC_Q_UF")
  
  BetaParameterTDC_Overall <- extract_relevant_nonlatency_data(model,
                                                               "(beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar +
                 beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar +
                 beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3 < 0",
                                                               "BetaParameterTDC_Overall")
  
  # Differences between Autism Group and Typical Development Group with Long Pauses
  OverallBetaParameterASDLowerTD <- extract_relevant_nonlatency_data(model,
                                                                     "(beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar + beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar) / 3 < ((beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar + beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3)",
                                                                     "OverallBetaParameterASD<TD")
  
  BetaParameterASDFasterTD_MG <- extract_relevant_nonlatency_data(model,
                                                                  "beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar",
                                                                  "BetaParameter_ASD_MG<TD_MG")
  
  BetaParameterASDFasterTD_Q_F <- extract_relevant_nonlatency_data(model,
                                                                   "beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                   "BetaParameter_ASD_MG<TD_Q_F")
  
  BetaParameterASDFasterTD_Q_UF <- extract_relevant_nonlatency_data(model,
                                                                    "beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                    "BetaParameter_ASD_MG<TD_Q_UF")
  
  # Differences within Autism Group across Conditions:
  WithinASDGroupMGSlowerQuestionsFamiliar <- extract_relevant_data(model,
                                                                   "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                   "WithinASDGroupMG>QuestionsFamiliar")
  
  # Differences within TD Group across Conditions:
  WithinTDGroupMGSlowerQuestionsFamiliar <- extract_relevant_data(model,
                                                                  "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                  "WithinTDGroupMG>QuestionsFamiliar")
  
  # Differences within Autism Group across Familiarity:
  WithinAutismGroupQuestionsUnfamiliarSlowerFamiliar <- extract_relevant_data(model,
                                                                              "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                              "WithinAutismGroupQuestionsUnfamiliar>Familiar")
  
  # Differences within TD Group across Familiarity:
  WithinTDGroupQuestionsUnfamiliarSlowerFamiliar <- extract_relevant_data(model,
                                                                          "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                          "WithinTDGroupQuestionsUnfamiliar>Familiar")
  
  # Differences between Autism Group and Typical Development Group Within MatchingGame
  WithinConditionsMatchingGameASDFasterTD <- extract_relevant_data(model,
                                                                   "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar",
                                                                   "WithinConditionsMatchingGameASD<TD")
  
  # Differences between Autism Group and Typical Development Group Within Familiar Questions:
  WithinConditionsQuestionsFamiliarASDFasterTD <- extract_relevant_data(model,
                                                                        "DiagnosisASD:TaskQuestions:FamiliarityFamiliar < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                        "WithinConditionsQuestionsFamiliarASD<TD")
  
  # Differences between Autism Group and Typical Development Group Within Unfamiliar Questions:
  WithinConditionsQuestionsUnfamiliarASDFasterTD <- extract_relevant_data(model,
                                                                          "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                          "WithinConditionsQuestionsUnfamiliarASD<TD")
  
  # Differences between Autism Group and Typical Development Group in Differences within Familiarity
  DifferenceInDifferenceBetweenGroupsUnfamiliarBiggerFamiliarity <- extract_relevant_data(model,
                                                                                          "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) > 
                      (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                          "DifferenceInDifferenceBetweenGroupsUnfamiliar>Familiarity")
  
  # Differences between Autism Group and Typical Development Group in Differences within Conditions
  DifferenceInDifferenceBetweenGroupsMGBiggerThanConversations <- extract_relevant_data(model,
                                                                                        "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar) > 
                      (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                        "DifferenceInDifferenceBetweenGroupsMG>Conversations")
  
  # Overall Subject Variability in Autismm Group
  SubjectVariabilityASD <- extract_relevant_sd_data(model, 
                                                    "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisASD + TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD)/3 > 0",
                                                    "SubjectVariabilityASD")
  
  SubjectVariabilityASD_MG <- extract_relevant_sd_data(model, 
                                                       "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD > 0",
                                                       "SubjectVariabilityASD_MG")
  
  SubjectVariabilityASD_Q_F <- extract_relevant_sd_data(model, 
                                                        "TaskQuestions:FamiliarityFamiliar:DiagnosisASD > 0",
                                                        "SubjectVariabilityASD_Q_F")
  
  SubjectVariabilityASD_Q_UF <- extract_relevant_sd_data(model, 
                                                         "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD > 0",
                                                         "SubjectVariabilityASD_Q_UF")
  
  # Overall Subject Variability in TD Group
  SubjectVariabilityTD <- extract_relevant_sd_data(model, 
                                                   "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisTDC + TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC)/3 > 0",
                                                   "SubjectVariabilityTD")
  
  SubjectVariabilityTD_MG <- extract_relevant_sd_data(model, 
                                                      "TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC > 0",
                                                      "SubjectVariabilityTD_MG")
  
  SubjectVariabilityTD_Q_F <- extract_relevant_sd_data(model, 
                                                       "TaskQuestions:FamiliarityFamiliar:DiagnosisTDC > 0",
                                                       "SubjectVariabilityTD_Q_F")
  
  SubjectVariabilityTD_Q_UF <- extract_relevant_sd_data(model, 
                                                        "TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC > 0",
                                                        "SubjectVariabilityTD_Q_UF")
  
  # Overall Subject Variability in TD Group
  SubjectVariabilityASDSmallerThanTD <- extract_relevant_sd_data(model,
                                                                 "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisASD + TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD) / 3 <
                                              (TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisTDC + TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC) / 3",
                                                                 "SubjectVariabilityASD<TD")
  
  
  SubjectVariabilityASDSmallerThanTD_MG <- extract_relevant_sd_data(model,
                                                                    "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD  <
                                              TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC",
                                                                    "SubjectVariabilityASD_MG<TD_MG")
  
  SubjectVariabilityASDSmallerThanTD_Q_F <- extract_relevant_sd_data(model,
                                                                     "TaskQuestions:FamiliarityFamiliar:DiagnosisASD  <
                                              TaskQuestions:FamiliarityFamiliar:DiagnosisTDC",
                                                                     "SubjectVariabilityASD_Q_F<TD_Q_F")
  
  SubjectVariabilityASDSmallerThanTD_Q_UF <- extract_relevant_sd_data(model,
                                                                      "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD  <
                                              TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC",
                                                                      "SubjectVariabilityASD_Q_UF<TD_Q_UF")
  
  DifferenceInDifferenceSubjectVariabilityUnfamiliarSmallerFamiliarity <- extract_relevant_sd_data(model,
                                                                                                   "(TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD - TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC) 
                                              < 
                                              (TaskQuestions:FamiliarityFamiliar:DiagnosisASD - TaskQuestions:FamiliarityFamiliar:DiagnosisTDC)",
                                                                                                   "DifferenceInDifferenceSubjectVariabilityUnfamiliarBiggerFamiliarity")
  
  DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions <- extract_relevant_sd_data(model,
                                                                                                  "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC - TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD) 
                                              < 
                                              (TaskQuestions:FamiliarityFamiliar:DiagnosisTDC - TaskQuestions:FamiliarityFamiliar:DiagnosisASD)",
                                                                                                  "DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions")
  
  # Overall Visit Variability across visits in Autism Group
  VisitVariabilityASDOverall <- extract_relevant_sdvisit_data(model, 
                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
                                                              "VisitVariabilityASDOverall")
  
  # Overall Visit Variability across visits in TD Group
  VisitVariabilityTDCOverall <- extract_relevant_sdvisit_data(model, 
                                                              "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
                                                              "VisitVariabilityTDCOverall")
  
  
  VisitVariabilityASDOverallHigherTDSOverall <- extract_relevant_sdvisit_data(model, 
                                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 >
                                               (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3",
                                                                              "VisitVariabilityASDOverallLowerTDSOverall")
  
  # Subject Variability across visits MG in Autism Group
  VisitVariabilityASDMG <- extract_relevant_sdvisit_data(model, 
                                                         "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                         "VisitVariabilityASDMG")
  
  # Overall Subject Variability across visits in TD Group
  VisitVariabilityTDMG <- extract_relevant_sdvisit_data(model, 
                                                        "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                        "VisitVariabilityTDMG")
  
  # Subject Variability across visits MG in Autism Group
  VisitVariabilityASD_Q_F <- extract_relevant_sdvisit_data(model, 
                                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > 0",
                                                           "VisitVariabilityASD_Q_F")
  
  # Overall Subject Variability across visits in TD Group
  VisitVariabilityTDC_Q_F <- extract_relevant_sdvisit_data(model, 
                                                           "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar > 0",
                                                           "VisitVariabilityTDC_Q_F")
  
  # Subject Variability across visits MG in Autism Group
  VisitVariabilityASD_Q_UF <- extract_relevant_sdvisit_data(model, 
                                                            "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                            "VisitVariabilityASD_Q_UF")
  
  # Overall Subject Variability across visits in TD Group
  VisitVariabilityTDC_Q_UF <- extract_relevant_sdvisit_data(model, 
                                                            "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                            "VisitVariabilityTDC_Q_UF")
  
  # Overall Subject Variability across visits, differences
  VisitVariabilityASDMGSlowerThanTDMG <- extract_relevant_sdvisit_data(model, 
                                                                       "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar",
                                                                       "VisitVariabilityASDMG>TDMG")
  
  # Overall Subject Variability across visits, differences
  VisitVariabilityASD_Q_F_SlowerThanTD_Q_F <- extract_relevant_sdvisit_data(model, 
                                                                            "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                            "VisitVariabilityASD_Q_F_SlowerThanTD_Q_F")
  
  # Overall Subject Variability across visits, differences
  VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF <- extract_relevant_sdvisit_data(model, 
                                                                              "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                              "VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF")
  
  
  #DifferenceInDifference:
  DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity <- extract_relevant_sdvisit_data(model, 
                                                                                                       "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)
                                               < (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                                       "DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity")
  
  DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG <- extract_relevant_sdvisit_data(model, 
                                                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar)
                                               > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                              "DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG")
  
  
  # Overall Sigma Variability across visits in Autism Group
  SigmaASDOverall <- extract_relevant_sigma_data(model, 
                                                 "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 < 0",
                                                 "SigmaASDOverall")
  
  # Overall Sigma Variability across visits in TD Group
  SigmaTDCOverall <- extract_relevant_sigma_data(model, 
                                                 "(sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) /3 < 0",
                                                 "SigmaTDCOverall")
  
  
  SigmaASDOverallHigherTDSOverall <- extract_relevant_nonlatency_data(model, 
                                                                      "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 <
                                               (sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3",
                                                                      "SigmaASDOverallHigherTDSOverall")
  
  # Subject Variability across visits MG in Autism Group
  SigmaASDMG <- extract_relevant_sigma_data(model, 
                                            "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar < 0",
                                            "SigmaASDMG")
  
  # Overall Subject Variability across visits in TD Group
  SigmaTDMG <- extract_relevant_sigma_data(model, 
                                           "sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar < 0",
                                           "SigmaTDMG")
  
  # Subject Variability across visits MG in Autism Group
  SigmaASD_Q_F <- extract_relevant_sigma_data(model, 
                                              "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < 0",
                                              "SigmaASD_Q_F")
  
  # Overall Subject Variability across visits in TD Group
  SigmaTDC_Q_F <- extract_relevant_sigma_data(model, 
                                              "sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar < 0",
                                              "SigmaTDC_Q_F")
  
  
  # Subject Variability across visits MG in Autism Group
  SigmaASD_Q_UF <- extract_relevant_sigma_data(model, 
                                               "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < 0 ",
                                               "SigmaASD_Q_UF")
  
  # Overall Subject Variability across visits in TD Group
  SigmaTDC_Q_UF <- extract_relevant_sigma_data(model, 
                                               "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar < 0",
                                               "SigmaTDC_Q_UF")
  
  # Overall Subject Variability across visits, differences
  SigmaASDMGSlowerThanTDMG <- extract_relevant_nonlatency_data(model, 
                                                               "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar",
                                                               "SigmaASDMGSlowerThanTDMG")
  
  # Overall Subject Variability across visits, differences
  SigmaASD_Q_F_SlowerThanTD_Q_F <- extract_relevant_nonlatency_data(model, 
                                                                    "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                    "SigmaASD_Q_F_SlowerThanTD_Q_F")
  
  # Overall Subject Variability across visits, differences
  SigmaASD_Q_UF_SlowerThanTD_Q_UF <- extract_relevant_nonlatency_data(model, 
                                                                      "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                      "SigmaASD_Q_UF_SlowerThanTD_Q_UF")
  
  #DifferenceInDifference:
  DifferenceInDifferenceSigmaASDFamiliarityBiggerTDUnfamiliarity <- extract_relevant_nonlatency_data(model, 
                                                                                                     "(sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)
                                               < (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                                     "DifferenceInDifferenceSigmaASDFamiliarityBiggerTDUnfamiliarity")
  
  DifferenceInDifferenceSigmaASDMGSmallerThanTDMG <- extract_relevant_nonlatency_data(model, 
                                                                                      "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar)
                                               > (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                      "DifferenceInDifferenceSigmaASDMGSmallerThanTDMG")
  
  
  # Correlation Estimates across visits in Autism Group
  post <- as_draws_df(model)
  
  
  # Apply the function for each case
  cor_post_mu_MG_ASD <- Visitextract_and_rename(post, "DiagnosisASD", "TaskMatchingGame", "FamiliarityFamiliar")
  cor_post_mu_MG_TDC <- Visitextract_and_rename(post, "DiagnosisTDC", "TaskMatchingGame", "FamiliarityFamiliar")
  
  cor_post_mu_Q_ASD_Familiar <- Visitextract_and_rename(post, "DiagnosisASD", "TaskQuestions", "FamiliarityFamiliar") %>%
    select(Visit2, Visit4, Visit6)
  
  cor_post_mu_Q_TDC_Familiar <- Visitextract_and_rename(post, "DiagnosisTDC", "TaskQuestions", "FamiliarityFamiliar") %>%
    select(Visit2, Visit4, Visit6)
  
  cor_post_mu_Q_ASD_Unfamiliar <- Visitextract_and_rename(post, "DiagnosisASD", "TaskQuestions", "FamiliarityUnfamiliar") %>%
    select(Visit1, Visit3, Visit5, Visit7)
  
  cor_post_mu_Q_TDC_Unfamiliar <- Visitextract_and_rename(post, "DiagnosisTDC", "TaskQuestions", "FamiliarityUnfamiliar") %>%
    select(Visit1, Visit3, Visit5, Visit7)
  
  
  CorrelationSummaryStatistics <- tibble(
    Estimate = c(mean(cor(cor_post_mu_MG_ASD)), mean(cor(cor_post_mu_Q_ASD_Familiar)), mean(cor(cor_post_mu_Q_ASD_Unfamiliar)), 
                 mean(cor(cor_post_mu_MG_TDC)), mean(cor(cor_post_mu_Q_TDC_Familiar)), mean(cor(cor_post_mu_Q_TDC_Unfamiliar))),
    CI.Lower = c(quantile(cor(cor_post_mu_MG_ASD), probs = 0.025), quantile(cor(cor_post_mu_Q_ASD_Familiar), probs = 0.025), quantile(cor(cor_post_mu_Q_ASD_Unfamiliar), probs = 0.025), 
                 quantile(cor(cor_post_mu_MG_TDC), probs = 0.025), quantile(cor(cor_post_mu_Q_TDC_Familiar), probs = 0.025), quantile(cor(cor_post_mu_Q_TDC_Unfamiliar), probs = 0.025)),
    CI.Upper = c(quantile(cor(cor_post_mu_MG_ASD), probs = 0.975), quantile(cor(cor_post_mu_Q_ASD_Familiar), probs = 0.975), quantile(cor(cor_post_mu_Q_ASD_Unfamiliar), probs = 0.975), 
                 quantile(cor(cor_post_mu_MG_TDC), probs = 0.975), quantile(cor(cor_post_mu_Q_TDC_Familiar), probs = 0.975), quantile(cor(cor_post_mu_Q_TDC_Unfamiliar), probs = 0.975)),
    Evid.Ratio = c(NA, NA, NA, NA, NA, NA),
    TestDescription = c(
      "CorrelationAcrossVisitsASDMG", "CorrelationAcrossVisitsASDFamiliarQuestions", "CorrelationAcrossVisitsASDUnfamiliarQuestions",
      "CorrelationAcrossVisitsTDMG", "CorrelationAcrossVisitsTDFamiliarQuestions", "CorrelationAcrossVisitsTDUnfamiliarQuestions"),
  )
  
  HypothesisResultsAdultLatenciesModel1 <- rbind(OverallLatencyThreeConditionsASD,
                                                 OverallLatencyThreeConditionsTD,
                                                 MGOverallLatencyASD,
                                                 Q_F_OverallLatencyASD,
                                                 Q_UF_OverallLatencyASD,
                                                 MGOverallLatencyTDC,
                                                 Q_F_OverallLatencyTDC,
                                                 Q_UF_OverallLatencyTDC,
                                                 BetaParameterASD_MG,
                                                 BetaParameterASD_Q_F,
                                                 BetaParameterASD_Q_UF,
                                                 BetaParameterASD_Overall,
                                                 BetaParameterTDC_MG,
                                                 BetaParameterTDC_Q_F,
                                                 BetaParameterTDC_Q_UF,
                                                 BetaParameterTDC_Overall,
                                                 OverallBetaParameterASDLowerTD,
                                                 BetaParameterASDFasterTD_MG,
                                                 BetaParameterASDFasterTD_Q_F,
                                                 BetaParameterASDFasterTD_Q_UF,
                                                 OverallLatencyThreeConditionsASDFasterTD,
                                                 WithinASDGroupMGSlowerQuestionsFamiliar,
                                                 WithinTDGroupMGSlowerQuestionsFamiliar,
                                                 WithinAutismGroupQuestionsUnfamiliarSlowerFamiliar,
                                                 WithinTDGroupQuestionsUnfamiliarSlowerFamiliar,
                                                 WithinConditionsMatchingGameASDFasterTD,
                                                 WithinConditionsQuestionsFamiliarASDFasterTD,
                                                 WithinConditionsQuestionsUnfamiliarASDFasterTD,
                                                 DifferenceInDifferenceBetweenGroupsUnfamiliarBiggerFamiliarity,
                                                 DifferenceInDifferenceBetweenGroupsMGBiggerThanConversations,
                                                 SubjectVariabilityASD,
                                                 SubjectVariabilityTD,
                                                 SubjectVariabilityASDSmallerThanTD,
                                                 SubjectVariabilityASD_MG,
                                                 SubjectVariabilityASD_Q_F,
                                                 SubjectVariabilityASD_Q_UF,
                                                 SubjectVariabilityTD_MG,
                                                 SubjectVariabilityTD_Q_F,
                                                 SubjectVariabilityTD_Q_UF,
                                                 SubjectVariabilityASDSmallerThanTD_MG,
                                                 SubjectVariabilityASDSmallerThanTD_Q_F,
                                                 SubjectVariabilityASDSmallerThanTD_Q_UF,
                                                 DifferenceInDifferenceSubjectVariabilityUnfamiliarSmallerFamiliarity,
                                                 DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions,
                                                 VisitVariabilityASDMG,
                                                 VisitVariabilityTDMG,
                                                 VisitVariabilityASD_Q_F,
                                                 VisitVariabilityTDC_Q_F,
                                                 VisitVariabilityASD_Q_UF,
                                                 VisitVariabilityTDC_Q_UF,
                                                 VisitVariabilityASDOverall,
                                                 VisitVariabilityTDCOverall,
                                                 VisitVariabilityASDMGSlowerThanTDMG,
                                                 VisitVariabilityASD_Q_F_SlowerThanTD_Q_F,
                                                 VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF,
                                                 VisitVariabilityASDOverallHigherTDSOverall,
                                                 DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG,
                                                 DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity,
                                                 SigmaASDOverall,
                                                 SigmaTDCOverall,
                                                 SigmaASDOverallHigherTDSOverall,
                                                 SigmaASDMG,
                                                 SigmaTDMG,
                                                 SigmaASD_Q_F,
                                                 SigmaTDC_Q_F,
                                                 SigmaASD_Q_UF,
                                                 SigmaTDC_Q_UF,
                                                 SigmaASDMGSlowerThanTDMG,
                                                 SigmaASD_Q_F_SlowerThanTD_Q_F,
                                                 SigmaASD_Q_UF_SlowerThanTD_Q_UF,
                                                 DifferenceInDifferenceSigmaASDFamiliarityBiggerTDUnfamiliarity,
                                                 DifferenceInDifferenceSigmaASDMGSmallerThanTDMG,
                                                 CorrelationSummaryStatistics) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  return(HypothesisResultsAdultLatenciesModel1)
}

extract_all_hypothesis_tests_m1_nosigma_child <- function(model) {
  
  
  # Overall Estimates in Autism Group
  OverallLatencyThreeConditionsASD <- extract_relevant_data(ChildLatency_m1_no_sigma_no_beta,
                                                            "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar +DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
                                                            "OverallLatencyThreeConditionsASD")
  
  MGOverallLatencyASD <- extract_relevant_data(ChildLatency_m1_no_sigma_no_beta,
                                               "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > 0",
                                               "MGOverallLatencyASD>0")
  
  Q_F_OverallLatencyASD <- extract_relevant_data(ChildLatency_m1_no_sigma_no_beta,
                                                 "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > 0",
                                                 "Q_F_OverallLatencyASD>0")
  
  Q_UF_OverallLatencyASD <- extract_relevant_data(ChildLatency_m1_no_sigma_no_beta,
                                                  "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                  "Q_UF_OverallLatencyASD>0")
  
  # Overall Estimates in TD Group
  OverallLatencyThreeConditionsTD <- extract_relevant_data(ChildLatency_m1_no_sigma_no_beta,
                                                           "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar +DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3 > 0",
                                                           "OverallLatencyThreeConditionsTD")
  
  MGOverallLatencyTDC <- extract_relevant_data(ChildLatency_m1_no_sigma_no_beta,
                                               "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > 0",
                                               "MGOverallLatencyTDC>0")
  
  Q_F_OverallLatencyTDC <- extract_relevant_data(ChildLatency_m1_no_sigma_no_beta,
                                                 "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar > 0",
                                                 "Q_F_OverallLatencyTDC>0")
  
  Q_UF_OverallLatencyTDC <- extract_relevant_data(ChildLatency_m1_no_sigma_no_beta,
                                                  "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                  "Q_UF_OverallLatencyTDC>0")
  
  # Overall Differences between Autism Group and Typical Development Group
  OverallLatencyThreeConditionsASDFasterTD <- extract_relevant_data(ChildLatency_m1_no_sigma_no_beta,
                                                                    "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar +DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 
                                <
                                (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar +DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3",
                                                                    "OverallLatencyThreeConditionsASD<TD")
  
  # Differences within Autism Group across Conditions:
  WithinASDGroupMGSlowerQuestionsFamiliar <- extract_relevant_data(ChildLatency_m1_no_sigma_no_beta,
                                                                   "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                   "WithinASDGroupMG>QuestionsFamiliar")
  
  # Differences within TD Group across Conditions:
  WithinTDGroupMGSlowerQuestionsFamiliar <- extract_relevant_data(ChildLatency_m1_no_sigma_no_beta,
                                                                  "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                  "WithinTDGroupMG>QuestionsFamiliar")
  
  # Differences within Autism Group across Familiarity:
  WithinAutismGroupQuestionsUnfamiliarSlowerFamiliar <- extract_relevant_data(ChildLatency_m1_no_sigma_no_beta,
                                                                              "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                              "WithinAutismGroupQuestionsUnfamiliar>Familiar")
  
  # Differences within TD Group across Familiarity:
  WithinTDGroupQuestionsUnfamiliarSlowerFamiliar <- extract_relevant_data(ChildLatency_m1_no_sigma_no_beta,
                                                                          "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                          "WithinTDGroupQuestionsUnfamiliar>Familiar")
  
  # Differences between Autism Group and Typical Development Group Within MatchingGame
  WithinConditionsMatchingGameASDFasterTD <- extract_relevant_data(ChildLatency_m1_no_sigma_no_beta,
                                                                   "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar < DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar",
                                                                   "WithinConditionsMatchingGameASD<TD")
  
  # Differences between Autism Group and Typical Development Group Within Familiar Questions:
  WithinConditionsQuestionsFamiliarASDFasterTD <- extract_relevant_data(ChildLatency_m1_no_sigma_no_beta,
                                                                        "DiagnosisASD:TaskQuestions:FamiliarityFamiliar < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                        "WithinConditionsQuestionsFamiliarASD<TD")
  
  # Differences between Autism Group and Typical Development Group Within Unfamiliar Questions:
  WithinConditionsQuestionsUnfamiliarASDFasterTD <- extract_relevant_data(ChildLatency_m1_no_sigma_no_beta,
                                                                          "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                          "WithinConditionsQuestionsUnfamiliarASD<TD")
  
  # Differences between Autism Group and Typical Development Group in Differences within Familiarity
  DifferenceInDifferenceBetweenGroupsUnfamiliarBiggerFamiliarity <- extract_relevant_data(ChildLatency_m1_no_sigma_no_beta,
                                                                                          "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) < 
                      (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                          "DifferenceInDifferenceBetweenGroupsUnfamiliar>Familiarity")
  
  # Differences between Autism Group and Typical Development Group in Differences within Conditions
  DifferenceInDifferenceBetweenGroupsMGBiggerThanConversations <- extract_relevant_data(ChildLatency_m1_no_sigma_no_beta,
                                                                                        "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar) > 
                      (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                        "DifferenceInDifferenceBetweenGroupsMG>Conversations")
  
  # Overall Subject Variability in Autismm Group
  SubjectVariabilityASD <- extract_relevant_sd_data(ChildLatency_m1_no_sigma_no_beta, 
                                                    "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisASD + TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD)/3 > 0",
                                                    "SubjectVariabilityASD")
  
  SubjectVariabilityASD_MG <- extract_relevant_sd_data(ChildLatency_m1_no_sigma_no_beta, 
                                                       "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD > 0",
                                                       "SubjectVariabilityASD_MG")
  
  SubjectVariabilityASD_Q_F <- extract_relevant_sd_data(ChildLatency_m1_no_sigma_no_beta, 
                                                        "TaskQuestions:FamiliarityFamiliar:DiagnosisASD > 0",
                                                        "SubjectVariabilityASD_Q_F")
  
  SubjectVariabilityASD_Q_UF <- extract_relevant_sd_data(ChildLatency_m1_no_sigma_no_beta, 
                                                         "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD > 0",
                                                         "SubjectVariabilityASD_Q_UF")
  
  # Overall Subject Variability in TD Group
  SubjectVariabilityTD <- extract_relevant_sd_data(ChildLatency_m1_no_sigma_no_beta, 
                                                   "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisTDC + TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC)/3 > 0",
                                                   "SubjectVariabilityTD")
  
  SubjectVariabilityTD_MG <- extract_relevant_sd_data(ChildLatency_m1_no_sigma_no_beta, 
                                                      "TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC > 0",
                                                      "SubjectVariabilityTD_MG")
  
  SubjectVariabilityTD_Q_F <- extract_relevant_sd_data(ChildLatency_m1_no_sigma_no_beta, 
                                                       "TaskQuestions:FamiliarityFamiliar:DiagnosisTDC > 0",
                                                       "SubjectVariabilityTD_Q_F")
  
  SubjectVariabilityTD_Q_UF <- extract_relevant_sd_data(ChildLatency_m1_no_sigma_no_beta, 
                                                        "TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC > 0",
                                                        "SubjectVariabilityTD_Q_UF")
  
  # Overall Subject Variability in TD Group
  SubjectVariabilityASDSmallerThanTD <- extract_relevant_sd_data(ChildLatency_m1_no_sigma_no_beta,
                                                                 "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisASD + TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD) / 3 <
                                              (TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisTDC + TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC) / 3",
                                                                 "SubjectVariabilityASD<TD")
  
  
  SubjectVariabilityASDSmallerThanTD_MG <- extract_relevant_sd_data(ChildLatency_m1_no_sigma_no_beta,
                                                                    "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD  <
                                              TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC",
                                                                    "SubjectVariabilityASD_MG<TD_MG")
  
  SubjectVariabilityASDSmallerThanTD_Q_F <- extract_relevant_sd_data(ChildLatency_m1_no_sigma_no_beta,
                                                                     "TaskQuestions:FamiliarityFamiliar:DiagnosisASD  <
                                              TaskQuestions:FamiliarityFamiliar:DiagnosisTDC",
                                                                     "SubjectVariabilityASD_Q_F<TD_Q_F")
  
  SubjectVariabilityASDSmallerThanTD_Q_UF <- extract_relevant_sd_data(ChildLatency_m1_no_sigma_no_beta,
                                                                      "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD  <
                                              TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC",
                                                                      "SubjectVariabilityASD_Q_UF<TD_Q_UF")
  
  SubjectVariabilityASD_MG_SmallerThan_Q_F <- extract_relevant_sd_data(ChildLatency_m1_no_sigma_no_beta,
                                                                       "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD > TaskQuestions:FamiliarityFamiliar:DiagnosisASD",
                                                                       "SubjectVariabilityASD_MG_SmallerThan_Q_F")
  
  SubjectVariabilityASD_Q_F_SmallerThan_Q_UF <- extract_relevant_sd_data(ChildLatency_m1_no_sigma_no_beta,
                                                                         "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD > TaskQuestions:FamiliarityFamiliar:DiagnosisASD",
                                                                         "SubjectVariabilityASD_Q_F_SmallerThan_Q_UF")
  
  SubjectVariabilityTDC_MG_SmallerThan_Q_F <- extract_relevant_sd_data(ChildLatency_m1_no_sigma_no_beta,
                                                                       "TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC < TaskQuestions:FamiliarityFamiliar:DiagnosisTDC",
                                                                       "SubjectVariabilityTDC_MG_SmallerThan_Q_F")
  
  SubjectVariabilityTDC_Q_F_SmallerThan_Q_UF <- extract_relevant_sd_data(ChildLatency_m1_no_sigma_no_beta,
                                                                         "TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC > TaskQuestions:FamiliarityFamiliar:DiagnosisTDC",
                                                                         "SubjectVariabilityTDC_Q_F_SmallerThan_Q_UF")
  
  DifferenceInDifferenceSubjectVariabilityUnfamiliarSmallerFamiliarity <- extract_relevant_sd_data(ChildLatency_m1_no_sigma_no_beta,
                                                                                                   "(TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD - TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC) 
                                              < 
                                              (TaskQuestions:FamiliarityFamiliar:DiagnosisASD - TaskQuestions:FamiliarityFamiliar:DiagnosisTDC)",
                                                                                                   "DifferenceInDifferenceSubjectVariabilityUnfamiliarBiggerFamiliarity")
  
  DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions <- extract_relevant_sd_data(ChildLatency_m1_no_sigma_no_beta,
                                                                                                  "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC - TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD) 
                                              < 
                                              (TaskQuestions:FamiliarityFamiliar:DiagnosisTDC - TaskQuestions:FamiliarityFamiliar:DiagnosisASD)",
                                                                                                  "DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions")
  
  # Overall Visit Variability across visits in Autism Group
  VisitVariabilityASDOverall <- extract_relevant_sdvisit_data(ChildLatency_m1_no_sigma_no_beta, 
                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
                                                              "VisitVariabilityASDOverall")
  
  # Overall Visit Variability across visits in TD Group
  VisitVariabilityTDCOverall <- extract_relevant_sdvisit_data(ChildLatency_m1_no_sigma_no_beta, 
                                                              "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
                                                              "VisitVariabilityTDCOverall")
  
  
  VisitVariabilityASDOverallHigherTDSOverall <- extract_relevant_sdvisit_data(ChildLatency_m1_no_sigma_no_beta, 
                                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 >
                                               (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3",
                                                                              "VisitVariabilityASDOverallLowerTDSOverall")
  
  # Subject Variability across visits MG in Autism Group
  VisitVariabilityASDMG <- extract_relevant_sdvisit_data(ChildLatency_m1_no_sigma_no_beta, 
                                                         "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                         "VisitVariabilityASDMG")
  
  # Overall Subject Variability across visits in TD Group
  VisitVariabilityTDMG <- extract_relevant_sdvisit_data(ChildLatency_m1_no_sigma_no_beta, 
                                                        "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                        "VisitVariabilityTDMG")
  
  # Subject Variability across visits MG in Autism Group
  VisitVariabilityASD_Q_F <- extract_relevant_sdvisit_data(ChildLatency_m1_no_sigma_no_beta, 
                                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > 0",
                                                           "VisitVariabilityASD_Q_F")
  
  # Overall Subject Variability across visits in TD Group
  VisitVariabilityTDC_Q_F <- extract_relevant_sdvisit_data(ChildLatency_m1_no_sigma_no_beta, 
                                                           "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar > 0",
                                                           "VisitVariabilityTDC_Q_F")
  
  # Subject Variability across visits MG in Autism Group
  VisitVariabilityASD_Q_UF <- extract_relevant_sdvisit_data(ChildLatency_m1_no_sigma_no_beta, 
                                                            "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                            "VisitVariabilityASD_Q_UF")
  
  # Overall Subject Variability across visits in TD Group
  VisitVariabilityTDC_Q_UF <- extract_relevant_sdvisit_data(ChildLatency_m1_no_sigma_no_beta, 
                                                            "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                            "VisitVariabilityTDC_Q_UF")
  
  # Overall Subject Variability across visits, differences
  VisitVariabilityASDMGSlowerThanTDMG <- extract_relevant_sdvisit_data(ChildLatency_m1_no_sigma_no_beta, 
                                                                       "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar < DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar",
                                                                       "VisitVariabilityASDMG>TDMG")
  
  # Overall Subject Variability across visits, differences
  VisitVariabilityASD_Q_F_SlowerThanTD_Q_F <- extract_relevant_sdvisit_data(ChildLatency_m1_no_sigma_no_beta, 
                                                                            "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                            "VisitVariabilityASD_Q_F_SlowerThanTD_Q_F")
  
  # Overall Subject Variability across visits, differences
  VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF <- extract_relevant_sdvisit_data(ChildLatency_m1_no_sigma_no_beta, 
                                                                              "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                              "VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF")
  
  #DifferenceInDifference:
  DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity <- extract_relevant_sdvisit_data(ChildLatency_m1_no_sigma_no_beta, 
                                                                                                       "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)
                                               < (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                                       "DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity")
  
  DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG <- extract_relevant_sdvisit_data(ChildLatency_m1_no_sigma_no_beta, 
                                                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar)
                                               > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                              "DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG")
  
  
  # Correlation Estimates across visits in Autism Group
  post <- as_draws_df(ChildLatency_m1_no_sigma_no_beta)
  
  
  # Apply the function for each case
  cor_post_mu_MG_ASD <- Visitextract_and_rename(post, "DiagnosisASD", "TaskMatchingGame", "FamiliarityFamiliar") %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_MG_TDC <- Visitextract_and_rename(post, "DiagnosisTDC", "TaskMatchingGame", "FamiliarityFamiliar") %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_Q_ASD_Familiar <- Visitextract_and_rename(post, "DiagnosisASD", "TaskQuestions", "FamiliarityFamiliar") %>%
    select(Visit2, Visit4, Visit6) %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_Q_TDC_Familiar <- Visitextract_and_rename(post, "DiagnosisTDC", "TaskQuestions", "FamiliarityFamiliar") %>%
    select(Visit2, Visit4, Visit6) %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_Q_ASD_Unfamiliar <- Visitextract_and_rename(post, "DiagnosisASD", "TaskQuestions", "FamiliarityUnfamiliar") %>%
    select(Visit1, Visit3, Visit5, Visit7) %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_Q_TDC_Unfamiliar <- Visitextract_and_rename(post, "DiagnosisTDC", "TaskQuestions", "FamiliarityUnfamiliar") %>%
    select(Visit1, Visit3, Visit5, Visit7) %>%
    cor() %>%
    as.matrix()
  
  diag(cor_post_mu_MG_ASD) <- NA
  diag(cor_post_mu_MG_TDC) <- NA
  diag(cor_post_mu_Q_ASD_Familiar) <- NA
  diag(cor_post_mu_Q_TDC_Familiar) <- NA
  diag(cor_post_mu_Q_ASD_Unfamiliar) <- NA
  diag(cor_post_mu_Q_TDC_Unfamiliar) <- NA
  
  CorrelationSummaryStatistics <- tibble(
    Estimate = c(mean(cor_post_mu_MG_ASD, na.rm = TRUE), mean(cor_post_mu_Q_ASD_Familiar, na.rm = TRUE), mean(cor_post_mu_Q_ASD_Unfamiliar, na.rm = TRUE), 
                 mean(cor_post_mu_MG_TDC, na.rm = TRUE), mean(cor_post_mu_Q_TDC_Familiar, na.rm = TRUE), mean(cor_post_mu_Q_TDC_Unfamiliar, na.rm = TRUE)),
    CI.Lower = c(quantile(cor_post_mu_MG_ASD, probs = 0.025, na.rm = TRUE), quantile(cor_post_mu_Q_ASD_Familiar, probs = 0.025, na.rm = TRUE), quantile(cor_post_mu_Q_ASD_Unfamiliar, probs = 0.025, na.rm = TRUE), 
                 quantile(cor_post_mu_MG_TDC, probs = 0.025, na.rm = TRUE), quantile(cor_post_mu_Q_TDC_Familiar, probs = 0.025, na.rm = TRUE), quantile(cor_post_mu_Q_TDC_Unfamiliar, probs = 0.025, na.rm = TRUE)),
    CI.Upper = c(quantile(cor_post_mu_MG_ASD, probs = 0.975, na.rm = TRUE), quantile(cor_post_mu_Q_ASD_Familiar, probs = 0.975, na.rm = TRUE), quantile(cor_post_mu_Q_ASD_Unfamiliar, probs = 0.975, na.rm = TRUE), 
                 quantile(cor_post_mu_MG_TDC, probs = 0.975, na.rm = TRUE), quantile(cor_post_mu_Q_TDC_Familiar, probs = 0.975, na.rm = TRUE), quantile(cor_post_mu_Q_TDC_Unfamiliar, probs = 0.975, na.rm = TRUE)),
    Evid.Ratio = c(NA, NA, NA, NA, NA, NA),
    TestDescription = c(
      "CorrelationAcrossVisitsASDMG", "CorrelationAcrossVisitsASDFamiliarQuestions", "CorrelationAcrossVisitsASDUnfamiliarQuestions",
      "CorrelationAcrossVisitsTDMG", "CorrelationAcrossVisitsTDFamiliarQuestions", "CorrelationAcrossVisitsTDUnfamiliarQuestions"),
  )
  
  HypothesisResultsChildLatenciesModel1 <- rbind(OverallLatencyThreeConditionsASD,
                                                 OverallLatencyThreeConditionsTD,
                                                 MGOverallLatencyASD,
                                                 Q_F_OverallLatencyASD,
                                                 Q_UF_OverallLatencyASD,
                                                 MGOverallLatencyTDC,
                                                 Q_F_OverallLatencyTDC,
                                                 Q_UF_OverallLatencyTDC,
                                                 OverallLatencyThreeConditionsASDFasterTD,
                                                 WithinASDGroupMGSlowerQuestionsFamiliar,
                                                 WithinTDGroupMGSlowerQuestionsFamiliar,
                                                 WithinAutismGroupQuestionsUnfamiliarSlowerFamiliar,
                                                 WithinTDGroupQuestionsUnfamiliarSlowerFamiliar,
                                                 WithinConditionsMatchingGameASDFasterTD,
                                                 WithinConditionsQuestionsFamiliarASDFasterTD,
                                                 WithinConditionsQuestionsUnfamiliarASDFasterTD,
                                                 DifferenceInDifferenceBetweenGroupsUnfamiliarBiggerFamiliarity,
                                                 DifferenceInDifferenceBetweenGroupsMGBiggerThanConversations,
                                                 SubjectVariabilityASD,
                                                 SubjectVariabilityTD,
                                                 SubjectVariabilityASDSmallerThanTD,
                                                 SubjectVariabilityASD_MG,
                                                 SubjectVariabilityASD_Q_F,
                                                 SubjectVariabilityASD_Q_UF,
                                                 SubjectVariabilityTD_MG,
                                                 SubjectVariabilityTD_Q_F,
                                                 SubjectVariabilityTD_Q_UF,
                                                 SubjectVariabilityASDSmallerThanTD_MG,
                                                 SubjectVariabilityASDSmallerThanTD_Q_F,
                                                 SubjectVariabilityASDSmallerThanTD_Q_UF,
                                                 SubjectVariabilityASD_MG_SmallerThan_Q_F,
                                                 SubjectVariabilityASD_Q_F_SmallerThan_Q_UF,
                                                 SubjectVariabilityTDC_MG_SmallerThan_Q_F,
                                                 SubjectVariabilityTDC_Q_F_SmallerThan_Q_UF,
                                                 DifferenceInDifferenceSubjectVariabilityUnfamiliarSmallerFamiliarity,
                                                 DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions,
                                                 VisitVariabilityASDMG,
                                                 VisitVariabilityTDMG,
                                                 VisitVariabilityASD_Q_F,
                                                 VisitVariabilityTDC_Q_F,
                                                 VisitVariabilityASD_Q_UF,
                                                 VisitVariabilityTDC_Q_UF,
                                                 VisitVariabilityASDOverall,
                                                 VisitVariabilityTDCOverall,
                                                 VisitVariabilityASDMGSlowerThanTDMG,
                                                 VisitVariabilityASD_Q_F_SlowerThanTD_Q_F,
                                                 VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF,
                                                 VisitVariabilityASDOverallHigherTDSOverall,
                                                 DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG,
                                                 DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity,
                                                 CorrelationSummaryStatistics) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  return(HypothesisResultsChildLatenciesModel1)
}

extract_all_hypothesis_tests_m1_nobetasigma_child <- function(model) {
  # Overall Estimates in Autism Group
  OverallLatencyThreeConditionsASD <- extract_relevant_data(ChildLatency_m1_no_beta_with_sigma,
                                                            "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar +DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
                                                            "OverallLatencyThreeConditionsASD")
  
  MGOverallLatencyASD <- extract_relevant_data(ChildLatency_m1_no_beta_with_sigma,
                                               "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > 0",
                                               "MGOverallLatencyASD>0")
  
  Q_F_OverallLatencyASD <- extract_relevant_data(ChildLatency_m1_no_beta_with_sigma,
                                                 "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > 0",
                                                 "Q_F_OverallLatencyASD>0")
  
  Q_UF_OverallLatencyASD <- extract_relevant_data(ChildLatency_m1_no_beta_with_sigma,
                                                  "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                  "Q_UF_OverallLatencyASD>0")
  
  # Overall Estimates in TD Group
  OverallLatencyThreeConditionsTD <- extract_relevant_data(ChildLatency_m1_no_beta_with_sigma,
                                                           "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar +DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3 > 0",
                                                           "OverallLatencyThreeConditionsTD")
  
  MGOverallLatencyTDC <- extract_relevant_data(ChildLatency_m1_no_beta_with_sigma,
                                               "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > 0",
                                               "MGOverallLatencyTDC>0")
  
  Q_F_OverallLatencyTDC <- extract_relevant_data(ChildLatency_m1_no_beta_with_sigma,
                                                 "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar > 0",
                                                 "Q_F_OverallLatencyTDC>0")
  
  Q_UF_OverallLatencyTDC <- extract_relevant_data(ChildLatency_m1_no_beta_with_sigma,
                                                  "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                  "Q_UF_OverallLatencyTDC>0")
  
  # Overall Differences between Autism Group and Typical Development Group
  OverallLatencyThreeConditionsASDFasterTD <- extract_relevant_data(ChildLatency_m1_no_beta_with_sigma,
                                                                    "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar +DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 
                                <
                                (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar +DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3",
                                                                    "OverallLatencyThreeConditionsASD<TD")
  
  # Differences within Autism Group across Conditions:
  WithinASDGroupMGSlowerQuestionsFamiliar <- extract_relevant_data(ChildLatency_m1_no_beta_with_sigma,
                                                                   "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                   "WithinASDGroupMG>QuestionsFamiliar")
  
  # Differences within TD Group across Conditions:
  WithinTDGroupMGSlowerQuestionsFamiliar <- extract_relevant_data(ChildLatency_m1_no_beta_with_sigma,
                                                                  "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                  "WithinTDGroupMG>QuestionsFamiliar")
  
  # Differences within Autism Group across Familiarity:
  WithinAutismGroupQuestionsUnfamiliarSlowerFamiliar <- extract_relevant_data(ChildLatency_m1_no_beta_with_sigma,
                                                                              "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                              "WithinAutismGroupQuestionsUnfamiliar>Familiar")
  
  # Differences within TD Group across Familiarity:
  WithinTDGroupQuestionsUnfamiliarSlowerFamiliar <- extract_relevant_data(ChildLatency_m1_no_beta_with_sigma,
                                                                          "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                          "WithinTDGroupQuestionsUnfamiliar>Familiar")
  
  # Differences between Autism Group and Typical Development Group Within MatchingGame
  WithinConditionsMatchingGameASDFasterTD <- extract_relevant_data(ChildLatency_m1_no_beta_with_sigma,
                                                                   "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar < DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar",
                                                                   "WithinConditionsMatchingGameASD<TD")
  
  # Differences between Autism Group and Typical Development Group Within Familiar Questions:
  WithinConditionsQuestionsFamiliarASDFasterTD <- extract_relevant_data(ChildLatency_m1_no_beta_with_sigma,
                                                                        "DiagnosisASD:TaskQuestions:FamiliarityFamiliar < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                        "WithinConditionsQuestionsFamiliarASD<TD")
  
  # Differences between Autism Group and Typical Development Group Within Unfamiliar Questions:
  WithinConditionsQuestionsUnfamiliarASDFasterTD <- extract_relevant_data(ChildLatency_m1_no_beta_with_sigma,
                                                                          "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                          "WithinConditionsQuestionsUnfamiliarASD<TD")
  
  # Differences between Autism Group and Typical Development Group in Differences within Familiarity
  DifferenceInDifferenceBetweenGroupsUnfamiliarBiggerFamiliarity <- extract_relevant_data(ChildLatency_m1_no_beta_with_sigma,
                                                                                          "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) < 
                      (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                          "DifferenceInDifferenceBetweenGroupsUnfamiliar>Familiarity")
  
  # Differences between Autism Group and Typical Development Group in Differences within Conditions
  DifferenceInDifferenceBetweenGroupsMGBiggerThanConversations <- extract_relevant_data(ChildLatency_m1_no_beta_with_sigma,
                                                                                        "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar) > 
                      (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                        "DifferenceInDifferenceBetweenGroupsMG>Conversations")
  
  # Overall Subject Variability in Autismm Group
  SubjectVariabilityASD <- extract_relevant_sd_data(ChildLatency_m1_no_beta_with_sigma, 
                                                    "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisASD + TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD)/3 > 0",
                                                    "SubjectVariabilityASD")
  
  SubjectVariabilityASD_MG <- extract_relevant_sd_data(ChildLatency_m1_no_beta_with_sigma, 
                                                       "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD > 0",
                                                       "SubjectVariabilityASD_MG")
  
  SubjectVariabilityASD_Q_F <- extract_relevant_sd_data(ChildLatency_m1_no_beta_with_sigma, 
                                                        "TaskQuestions:FamiliarityFamiliar:DiagnosisASD > 0",
                                                        "SubjectVariabilityASD_Q_F")
  
  SubjectVariabilityASD_Q_UF <- extract_relevant_sd_data(ChildLatency_m1_no_beta_with_sigma, 
                                                         "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD > 0",
                                                         "SubjectVariabilityASD_Q_UF")
  
  # Overall Subject Variability in TD Group
  SubjectVariabilityTD <- extract_relevant_sd_data(ChildLatency_m1_no_beta_with_sigma, 
                                                   "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisTDC + TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC)/3 > 0",
                                                   "SubjectVariabilityTD")
  
  SubjectVariabilityTD_MG <- extract_relevant_sd_data(ChildLatency_m1_no_beta_with_sigma, 
                                                      "TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC > 0",
                                                      "SubjectVariabilityTD_MG")
  
  SubjectVariabilityTD_Q_F <- extract_relevant_sd_data(ChildLatency_m1_no_beta_with_sigma, 
                                                       "TaskQuestions:FamiliarityFamiliar:DiagnosisTDC > 0",
                                                       "SubjectVariabilityTD_Q_F")
  
  SubjectVariabilityTD_Q_UF <- extract_relevant_sd_data(ChildLatency_m1_no_beta_with_sigma, 
                                                        "TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC > 0",
                                                        "SubjectVariabilityTD_Q_UF")
  
  # Overall Subject Variability in TD Group
  SubjectVariabilityASDSmallerThanTD <- extract_relevant_sd_data(ChildLatency_m1_no_beta_with_sigma,
                                                                 "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisASD + TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD) / 3 <
                                              (TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisTDC + TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC) / 3",
                                                                 "SubjectVariabilityASD<TD")
  
  
  SubjectVariabilityASDSmallerThanTD_MG <- extract_relevant_sd_data(ChildLatency_m1_no_beta_with_sigma,
                                                                    "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD  <
                                              TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC",
                                                                    "SubjectVariabilityASD_MG<TD_MG")
  
  SubjectVariabilityASDSmallerThanTD_Q_F <- extract_relevant_sd_data(ChildLatency_m1_no_beta_with_sigma,
                                                                     "TaskQuestions:FamiliarityFamiliar:DiagnosisASD  <
                                              TaskQuestions:FamiliarityFamiliar:DiagnosisTDC",
                                                                     "SubjectVariabilityASD_Q_F<TD_Q_F")
  
  SubjectVariabilityASDSmallerThanTD_Q_UF <- extract_relevant_sd_data(ChildLatency_m1_no_beta_with_sigma,
                                                                      "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD  <
                                              TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC",
                                                                      "SubjectVariabilityASD_Q_UF<TD_Q_UF")
  
  SubjectVariabilityASD_MG_SmallerThan_Q_F <- extract_relevant_sd_data(ChildLatency_m1_no_beta_with_sigma,
                                                                       "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD > TaskQuestions:FamiliarityFamiliar:DiagnosisASD",
                                                                       "SubjectVariabilityASD_MG_SmallerThan_Q_F")
  
  SubjectVariabilityASD_Q_F_SmallerThan_Q_UF <- extract_relevant_sd_data(ChildLatency_m1_no_beta_with_sigma,
                                                                         "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD > TaskQuestions:FamiliarityFamiliar:DiagnosisASD",
                                                                         "SubjectVariabilityASD_Q_F_SmallerThan_Q_UF")
  
  SubjectVariabilityTDC_MG_SmallerThan_Q_F <- extract_relevant_sd_data(ChildLatency_m1_no_beta_with_sigma,
                                                                       "TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC < TaskQuestions:FamiliarityFamiliar:DiagnosisTDC",
                                                                       "SubjectVariabilityTDC_MG_SmallerThan_Q_F")
  
  SubjectVariabilityTDC_Q_F_SmallerThan_Q_UF <- extract_relevant_sd_data(ChildLatency_m1_no_beta_with_sigma,
                                                                         "TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC > TaskQuestions:FamiliarityFamiliar:DiagnosisTDC",
                                                                         "SubjectVariabilityTDC_Q_F_SmallerThan_Q_UF")
  
  DifferenceInDifferenceSubjectVariabilityUnfamiliarSmallerFamiliarity <- extract_relevant_sd_data(ChildLatency_m1_no_beta_with_sigma,
                                                                                                   "(TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD - TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC) 
                                              < 
                                              (TaskQuestions:FamiliarityFamiliar:DiagnosisASD - TaskQuestions:FamiliarityFamiliar:DiagnosisTDC)",
                                                                                                   "DifferenceInDifferenceSubjectVariabilityUnfamiliarBiggerFamiliarity")
  
  DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions <- extract_relevant_sd_data(ChildLatency_m1_no_beta_with_sigma,
                                                                                                  "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC - TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD) 
                                              < 
                                              (TaskQuestions:FamiliarityFamiliar:DiagnosisTDC - TaskQuestions:FamiliarityFamiliar:DiagnosisASD)",
                                                                                                  "DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions")
  
  # Overall Visit Variability across visits in Autism Group
  VisitVariabilityASDOverall <- extract_relevant_sdvisit_data(ChildLatency_m1_no_beta_with_sigma, 
                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
                                                              "VisitVariabilityASDOverall")
  
  # Overall Visit Variability across visits in TD Group
  VisitVariabilityTDCOverall <- extract_relevant_sdvisit_data(ChildLatency_m1_no_beta_with_sigma, 
                                                              "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
                                                              "VisitVariabilityTDCOverall")
  
  
  VisitVariabilityASDOverallHigherTDSOverall <- extract_relevant_sdvisit_data(ChildLatency_m1_no_beta_with_sigma, 
                                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 >
                                               (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3",
                                                                              "VisitVariabilityASDOverallLowerTDSOverall")
  
  # Subject Variability across visits MG in Autism Group
  VisitVariabilityASDMG <- extract_relevant_sdvisit_data(ChildLatency_m1_no_beta_with_sigma, 
                                                         "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                         "VisitVariabilityASDMG")
  
  # Overall Subject Variability across visits in TD Group
  VisitVariabilityTDMG <- extract_relevant_sdvisit_data(ChildLatency_m1_no_beta_with_sigma, 
                                                        "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > 0",
                                                        "VisitVariabilityTDMG")
  
  # Subject Variability across visits MG in Autism Group
  VisitVariabilityASD_Q_F <- extract_relevant_sdvisit_data(ChildLatency_m1_no_beta_with_sigma, 
                                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > 0",
                                                           "VisitVariabilityASD_Q_F")
  
  # Overall Subject Variability across visits in TD Group
  VisitVariabilityTDC_Q_F <- extract_relevant_sdvisit_data(ChildLatency_m1_no_beta_with_sigma, 
                                                           "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar > 0",
                                                           "VisitVariabilityTDC_Q_F")
  
  # Subject Variability across visits MG in Autism Group
  VisitVariabilityASD_Q_UF <- extract_relevant_sdvisit_data(ChildLatency_m1_no_beta_with_sigma, 
                                                            "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                            "VisitVariabilityASD_Q_UF")
  
  # Overall Subject Variability across visits in TD Group
  VisitVariabilityTDC_Q_UF <- extract_relevant_sdvisit_data(ChildLatency_m1_no_beta_with_sigma, 
                                                            "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                            "VisitVariabilityTDC_Q_UF")
  
  # Overall Subject Variability across visits, differences
  VisitVariabilityASDMGSlowerThanTDMG <- extract_relevant_sdvisit_data(ChildLatency_m1_no_beta_with_sigma, 
                                                                       "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar < DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar",
                                                                       "VisitVariabilityASDMG>TDMG")
  
  # Overall Subject Variability across visits, differences
  VisitVariabilityASD_Q_F_SlowerThanTD_Q_F <- extract_relevant_sdvisit_data(ChildLatency_m1_no_beta_with_sigma, 
                                                                            "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                            "VisitVariabilityASD_Q_F_SlowerThanTD_Q_F")
  
  # Overall Subject Variability across visits, differences
  VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF <- extract_relevant_sdvisit_data(ChildLatency_m1_no_beta_with_sigma, 
                                                                              "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                              "VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF")
  
  #DifferenceInDifference:
  DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity <- extract_relevant_sdvisit_data(ChildLatency_m1_no_beta_with_sigma, 
                                                                                                       "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)
                                               < (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                                       "DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity")
  
  DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG <- extract_relevant_sdvisit_data(ChildLatency_m1_no_beta_with_sigma, 
                                                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar)
                                               > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                              "DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG")
  
  
  # Overall Sigma Variability across visits in Autism Group
  SigmaASDOverall <- extract_relevant_sigma_data(ChildLatency_m1_no_beta_with_sigma, 
                                                 "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 < 0",
                                                 "SigmaASDOverall")
  
  # Overall Sigma Variability across visits in TD Group
  SigmaTDCOverall <- extract_relevant_sigma_data(ChildLatency_m1_no_beta_with_sigma, 
                                                 "(sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) /3 < 0",
                                                 "SigmaTDCOverall")
  
  
  SigmaASDOverallHigherTDSOverall <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma, 
                                                                      "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 >
                                               (sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3",
                                                                      "SigmaASDOverallHigherTDSOverall")
  
  # Subject Variability across visits MG in Autism Group
  SigmaASDMG <- extract_relevant_sigma_data(ChildLatency_m1_no_beta_with_sigma, 
                                            "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar < 0",
                                            "SigmaASDMG")
  
  # Overall Subject Variability across visits in TD Group
  SigmaTDMG <- extract_relevant_sigma_data(ChildLatency_m1_no_beta_with_sigma, 
                                           "sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar < 0",
                                           "SigmaTDMG")
  
  # Subject Variability across visits MG in Autism Group
  SigmaASD_Q_F <- extract_relevant_sigma_data(ChildLatency_m1_no_beta_with_sigma, 
                                              "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < 0",
                                              "SigmaASD_Q_F")
  
  # Overall Subject Variability across visits in TD Group
  SigmaTDC_Q_F <- extract_relevant_sigma_data(ChildLatency_m1_no_beta_with_sigma, 
                                              "sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar < 0",
                                              "SigmaTDC_Q_F")
  
  
  # Subject Variability across visits MG in Autism Group
  SigmaASD_Q_UF <- extract_relevant_sigma_data(ChildLatency_m1_no_beta_with_sigma, 
                                               "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < 0 ",
                                               "SigmaASD_Q_UF")
  
  # Overall Subject Variability across visits in TD Group
  SigmaTDC_Q_UF <- extract_relevant_sigma_data(ChildLatency_m1_no_beta_with_sigma, 
                                               "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar < 0",
                                               "SigmaTDC_Q_UF")
  
  # Overall Subject Variability across visits, differences
  SigmaASDMGSlowerThanTDMG <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma, 
                                                               "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar",
                                                               "SigmaASDMGSlowerThanTDMG")
  
  # Overall Subject Variability across visits, differences
  SigmaASD_Q_F_SlowerThanTD_Q_F <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma, 
                                                                    "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                    "SigmaASD_Q_F_SlowerThanTD_Q_F")
  
  # Overall Subject Variability across visits, differences
  SigmaASD_Q_UF_SlowerThanTD_Q_UF <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma, 
                                                                      "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                      "SigmaASD_Q_UF_SlowerThanTD_Q_UF")
  
  #Sigma across visits, differences
  SigmaASD_MG_VisitLowerZero <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                 "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit < 0",
                                                                 "SigmaASD_MG_VisitLowerZero")
  
  SigmaTDC_MG_VisitLowerZero <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                 "sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit > 0",
                                                                 "SigmaTDC_MG_VisitLowerZero")
  
  SigmaASD_MG_VisitLowerTDC <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit < sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
                                                                "SigmaASD_MG_VisitLowerTDC")
  
  SigmaASD_Overall_VisitLowerZero <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                      "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit +
                                               sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit +
                                               sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit) / 3 < 0",
                                                                      "SigmaASD_Overall_VisitLowerZero")
  
  SigmaTDC_Overall_VisitLowerZero <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                      "(sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit +
                                               sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit +
                                               sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit) / 3 < 0",
                                                                      "SigmaTDC_Overall_VisitLowerZero")
  
  # Sigma across visits, differences
  SigmaASD_Q_F_VisitLowerZero <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                  "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit > 0",
                                                                  "SigmaASD_Q_F_VisitLowerZero")
  
  SigmaTDC_Q_F_VisitLowerZero <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                  "sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit < 0",
                                                                  "SigmaTDC_Q_F_VisitLowerZero")
  
  SigmaASD_Q_F_VisitLowerTDC <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                 "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit > sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
                                                                 "SigmaASD_Q_F_VisitLowerTDC")
  
  # Sigma across visits, differences
  SigmaASD_Q_UF_VisitLowerZero <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                   "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < 0",
                                                                   "SigmaASD_Q_UF_VisitLowerZero")
  
  SigmaTDC_Q_UF_VisitLowerZero <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                   "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > 0",
                                                                   "SigmaTDC_Q_UF_VisitLowerZero")
  
  SigmaASD_Q_UF_VisitLowerTDC <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                  "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
                                                                  "SigmaASD_Q_UF_VisitLowerTDC")
  
  SigmaVisitOverall_ASDLowerTDC <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                    "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit +
                                               sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit +
                                               sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit) / 3 > (sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit +
                                               sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit +
                                               sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit) /3",
                                                                    "SigmaVisitOverall_ASDLowerTDC")
  
  #DifferenceInDifference:
  DifferenceInDifferenceSigmaASDFamiliarityBiggerTDUnfamiliarity <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma, 
                                                                                                     "(sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)
                                               > (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                                     "DifferenceInDifferenceSigmaASDFamiliarityBiggerTDUnfamiliarity")
  
  DifferenceInDifferenceSigmaASDMGSmallerThanTDMG <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma, 
                                                                                      "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar)
                                               > (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                      "DifferenceInDifferenceSigmaASDMGSmallerThanTDMG")
  
  #WithinGroups
  WithinGroupSigmaASD_Q_F_Faster_UF <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma, 
                                                                        "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaASD_Q_FFaster_UF")
  
  WithinGroupSigmaTDC_Q_F_Faster_UF <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma, 
                                                                        "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaTDC_Q_F_Faster_UF")
  
  WithinGroupSigmaASD_Q_F_Faster_MG <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma, 
                                                                        "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaASD_Q_F_Faster_MG")
  
  WithinGroupSigmaTDC_Q_F_Faster_MG <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma, 
                                                                        "sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar < sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaTDC_Q_F_Faster_MG")
  
  
  WithinGroupSigmaASD_Q_F_Lower_UF_Visit <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                             "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit",
                                                                             "WithinGroupSigmaASD_Q_F_Lower_UF_Visit")
  
  WithinGroupSigmaTDC_Q_F_Lower_UF_Visit <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                             "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit",
                                                                             "WithinGroupSigmaTDC_Q_F_Lower_UF_Visit")
  
  WithinGroupSigmaASD_Q_F_Lower_MG_Visit <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                             "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit > sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit",
                                                                             "WithinGroupSigmaASD_Q_F_Lower_MG_Visit")
  
  WithinGroupSigmaTDC_Q_F_Lower_MG_Visit <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                             "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
                                                                             "WithinGroupSigmaTDC_Q_F_Lower_MG_Visit")
  
  
  SigmaTDC_Q_UF_VisitLowerZero <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                   "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > 0",
                                                                   "SigmaTDC_Q_UF_VisitLowerZero")
  
  SigmaASD_Q_UF_VisitLowerTDC <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                  "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
                                                                  "SigmaASD_Q_UF_VisitLowerTDC")
  
  
  DifferenceInDifferenceSigmaVisit_Q_F_Lower_Q_UF <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                                      "(sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit - sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit) <
                                               (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit)", "DifferenceInDifferenceSigmaVisit_Q_F_Lower_Q_UF")
  
  
  DifferenceInDifferenceSigmaVisit_MG_Lower_Q_F <- extract_relevant_nonlatency_data(ChildLatency_m1_no_beta_with_sigma,
                                                                                    "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit - sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit) <
                                               (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit)", "DifferenceInDifferenceSigmaVisit_MG_Lower_Q_F")
  
  
  
  # Correlation Estimates across visits in Autism Group
  post <- as_draws_df(ChildLatency_m1_no_beta_with_sigma)
  
  
  # Apply the function for each case
  cor_post_mu_MG_ASD <- Visitextract_and_rename(post, "DiagnosisASD", "TaskMatchingGame", "FamiliarityFamiliar") %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_MG_TDC <- Visitextract_and_rename(post, "DiagnosisTDC", "TaskMatchingGame", "FamiliarityFamiliar") %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_Q_ASD_Familiar <- Visitextract_and_rename(post, "DiagnosisASD", "TaskQuestions", "FamiliarityFamiliar") %>%
    select(Visit2, Visit4, Visit6) %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_Q_TDC_Familiar <- Visitextract_and_rename(post, "DiagnosisTDC", "TaskQuestions", "FamiliarityFamiliar") %>%
    select(Visit2, Visit4, Visit6) %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_Q_ASD_Unfamiliar <- Visitextract_and_rename(post, "DiagnosisASD", "TaskQuestions", "FamiliarityUnfamiliar") %>%
    select(Visit1, Visit3, Visit5, Visit7) %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_Q_TDC_Unfamiliar <- Visitextract_and_rename(post, "DiagnosisTDC", "TaskQuestions", "FamiliarityUnfamiliar") %>%
    select(Visit1, Visit3, Visit5, Visit7) %>%
    cor() %>%
    as.matrix()
  
  diag(cor_post_mu_MG_ASD) <- NA
  diag(cor_post_mu_MG_TDC) <- NA
  diag(cor_post_mu_Q_ASD_Familiar) <- NA
  diag(cor_post_mu_Q_TDC_Familiar) <- NA
  diag(cor_post_mu_Q_ASD_Unfamiliar) <- NA
  diag(cor_post_mu_Q_TDC_Unfamiliar) <- NA
  
  CorrelationSummaryStatistics <- tibble(
    Estimate = c(mean(cor_post_mu_MG_ASD, na.rm = TRUE), mean(cor_post_mu_Q_ASD_Familiar, na.rm = TRUE), mean(cor_post_mu_Q_ASD_Unfamiliar, na.rm = TRUE), 
                 mean(cor_post_mu_MG_TDC, na.rm = TRUE), mean(cor_post_mu_Q_TDC_Familiar, na.rm = TRUE), mean(cor_post_mu_Q_TDC_Unfamiliar, na.rm = TRUE)),
    CI.Lower = c(quantile(cor_post_mu_MG_ASD, probs = 0.025, na.rm = TRUE), quantile(cor_post_mu_Q_ASD_Familiar, probs = 0.025, na.rm = TRUE), quantile(cor_post_mu_Q_ASD_Unfamiliar, probs = 0.025, na.rm = TRUE), 
                 quantile(cor_post_mu_MG_TDC, probs = 0.025, na.rm = TRUE), quantile(cor_post_mu_Q_TDC_Familiar, probs = 0.025, na.rm = TRUE), quantile(cor_post_mu_Q_TDC_Unfamiliar, probs = 0.025, na.rm = TRUE)),
    CI.Upper = c(quantile(cor_post_mu_MG_ASD, probs = 0.975, na.rm = TRUE), quantile(cor_post_mu_Q_ASD_Familiar, probs = 0.975, na.rm = TRUE), quantile(cor_post_mu_Q_ASD_Unfamiliar, probs = 0.975, na.rm = TRUE), 
                 quantile(cor_post_mu_MG_TDC, probs = 0.975, na.rm = TRUE), quantile(cor_post_mu_Q_TDC_Familiar, probs = 0.975, na.rm = TRUE), quantile(cor_post_mu_Q_TDC_Unfamiliar, probs = 0.975, na.rm = TRUE)),
    Evid.Ratio = c(NA, NA, NA, NA, NA, NA),
    TestDescription = c(
      "CorrelationAcrossVisitsASDMG", "CorrelationAcrossVisitsASDFamiliarQuestions", "CorrelationAcrossVisitsASDUnfamiliarQuestions",
      "CorrelationAcrossVisitsTDMG", "CorrelationAcrossVisitsTDFamiliarQuestions", "CorrelationAcrossVisitsTDUnfamiliarQuestions"),
  )
  
  HypothesisResultsChildLatenciesModel1 <- rbind(OverallLatencyThreeConditionsASD,
                                                 OverallLatencyThreeConditionsTD,
                                                 MGOverallLatencyASD,
                                                 Q_F_OverallLatencyASD,
                                                 Q_UF_OverallLatencyASD,
                                                 MGOverallLatencyTDC,
                                                 Q_F_OverallLatencyTDC,
                                                 Q_UF_OverallLatencyTDC,
                                                 OverallLatencyThreeConditionsASDFasterTD,
                                                 WithinASDGroupMGSlowerQuestionsFamiliar,
                                                 WithinTDGroupMGSlowerQuestionsFamiliar,
                                                 WithinAutismGroupQuestionsUnfamiliarSlowerFamiliar,
                                                 WithinTDGroupQuestionsUnfamiliarSlowerFamiliar,
                                                 WithinConditionsMatchingGameASDFasterTD,
                                                 WithinConditionsQuestionsFamiliarASDFasterTD,
                                                 WithinConditionsQuestionsUnfamiliarASDFasterTD,
                                                 DifferenceInDifferenceBetweenGroupsUnfamiliarBiggerFamiliarity,
                                                 DifferenceInDifferenceBetweenGroupsMGBiggerThanConversations,
                                                 SubjectVariabilityASD,
                                                 SubjectVariabilityTD,
                                                 SubjectVariabilityASDSmallerThanTD,
                                                 SubjectVariabilityASD_MG,
                                                 SubjectVariabilityASD_Q_F,
                                                 SubjectVariabilityASD_Q_UF,
                                                 SubjectVariabilityTD_MG,
                                                 SubjectVariabilityTD_Q_F,
                                                 SubjectVariabilityTD_Q_UF,
                                                 SubjectVariabilityASDSmallerThanTD_MG,
                                                 SubjectVariabilityASDSmallerThanTD_Q_F,
                                                 SubjectVariabilityASDSmallerThanTD_Q_UF,
                                                 SubjectVariabilityASD_MG_SmallerThan_Q_F,
                                                 SubjectVariabilityASD_Q_F_SmallerThan_Q_UF,
                                                 SubjectVariabilityTDC_MG_SmallerThan_Q_F,
                                                 SubjectVariabilityTDC_Q_F_SmallerThan_Q_UF,
                                                 DifferenceInDifferenceSubjectVariabilityUnfamiliarSmallerFamiliarity,
                                                 DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions,
                                                 VisitVariabilityASDMG,
                                                 VisitVariabilityTDMG,
                                                 VisitVariabilityASD_Q_F,
                                                 VisitVariabilityTDC_Q_F,
                                                 VisitVariabilityASD_Q_UF,
                                                 VisitVariabilityTDC_Q_UF,
                                                 VisitVariabilityASDOverall,
                                                 VisitVariabilityTDCOverall,
                                                 VisitVariabilityASDMGSlowerThanTDMG,
                                                 VisitVariabilityASD_Q_F_SlowerThanTD_Q_F,
                                                 VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF,
                                                 VisitVariabilityASDOverallHigherTDSOverall,
                                                 DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG,
                                                 DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity,
                                                 SigmaASDOverall,
                                                 SigmaTDCOverall,
                                                 SigmaASDOverallHigherTDSOverall,
                                                 SigmaASDMG,
                                                 SigmaTDMG,
                                                 SigmaASD_Q_F,
                                                 SigmaTDC_Q_F,
                                                 SigmaASD_Q_UF,
                                                 SigmaTDC_Q_UF,
                                                 SigmaASDMGSlowerThanTDMG,
                                                 SigmaASD_Q_F_SlowerThanTD_Q_F,
                                                 SigmaASD_Q_UF_SlowerThanTD_Q_UF,
                                                 SigmaASD_MG_VisitLowerZero,
                                                 SigmaTDC_MG_VisitLowerZero,
                                                 SigmaASD_MG_VisitLowerTDC,
                                                 SigmaASD_Q_F_VisitLowerZero,
                                                 SigmaTDC_Q_F_VisitLowerZero,
                                                 SigmaASD_Q_F_VisitLowerTDC,
                                                 SigmaASD_Q_UF_VisitLowerZero,
                                                 SigmaTDC_Q_UF_VisitLowerZero,
                                                 SigmaASD_Q_UF_VisitLowerTDC,
                                                 SigmaASD_Overall_VisitLowerZero,
                                                 SigmaTDC_Overall_VisitLowerZero,
                                                 WithinGroupSigmaASD_Q_F_Lower_UF_Visit,
                                                 WithinGroupSigmaTDC_Q_F_Lower_UF_Visit,
                                                 WithinGroupSigmaASD_Q_F_Lower_MG_Visit,
                                                 WithinGroupSigmaTDC_Q_F_Lower_MG_Visit,
                                                 WithinGroupSigmaASD_Q_F_Faster_UF,
                                                 WithinGroupSigmaTDC_Q_F_Faster_UF,
                                                 WithinGroupSigmaASD_Q_F_Faster_MG,
                                                 WithinGroupSigmaTDC_Q_F_Faster_MG,
                                                 SigmaVisitOverall_ASDLowerTDC,
                                                 DifferenceInDifferenceSigmaVisit_Q_F_Lower_Q_UF,
                                                 DifferenceInDifferenceSigmaVisit_MG_Lower_Q_F,
                                                 DifferenceInDifferenceSigmaASDFamiliarityBiggerTDUnfamiliarity,
                                                 DifferenceInDifferenceSigmaASDMGSmallerThanTDMG,
                                                 CorrelationSummaryStatistics) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  return(HypothesisResultsChildLatenciesModel1)
}

extract_all_hypothesis_tests_m2_child <- function(model) {
  LanguageASD_TG_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                            "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled > 0",
                                            "LanguageASD_TG_F")
  
  LanguageASD_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled > 0",
                                           "LanguageASD_Q_F")
  
  LanguageASD_Q_UF <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                            "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > 0",
                                            "LanguageASD_Q_UF")
  
  LanguageASD_overall <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                               "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)/3 > 0",
                                               "LanguageASD_overall")
  
  LanguageTDC_TG_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                            "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled > 0",
                                            "LanguageTDC_TG_F")
  
  LanguageTDC_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                           "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled > 0",
                                           "LanguageTDC_Q_F")
  
  LanguageTDC_Q_UF <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                            "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > 0",
                                            "LanguageTDC_Q_UF")
  
  LanguageTDC_overall <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                               "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)/3 > 0",
                                               "LanguageTDC_overall")
  
  MotorASD_TG_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled < 0", "MotorASD_TG_F")
  MotorASD_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled < 0", "MotorASD_Q_F")
  MotorASD_Q_UF <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled < 0", "MotorASD_Q_UF")
  
  #1
  WithinGroupMotorASD_MGSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled", "WithinGroupMotorASD_MGSlower_Q_F")
  
  WithinGroupMotorASD_Q_UFSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled", "WithinGroupMotorASD_Q_UFSlower_Q_F")
  #2
  WithinGroupCognitionASD_MGSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled < DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled", "WithinGroupCognitionASD_MGSlower_Q_F")
  
  WithinGroupCognitionASD_Q_UFSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled", "WithinGroupCognitionASD_Q_UFSlower_Q_F")
  #3
  WithinGroupAwarenessASD_MGSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled < DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled", "WithinGroupAwarenessASD_MGSlower_Q_F")
  
  WithinGroupAwarenessASD_Q_UFSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled < DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled", "WithinGroupAwarenessASD_Q_UFSlower_Q_F")
  #4
  WithinGroupMotivationASD_MGSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled", "WithinGroupMotivationASD_MGSlower_Q_F")
  
  WithinGroupMotivationASD_Q_UFSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled", "WithinGroupMotivationASD_Q_UFSlower_Q_F")
  #5
  WithinGroupLanguageASD_MGSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled", "WithinGroupLanguageASD_MGSlower_Q_F")
  
  WithinGroupLanguageASD_Q_UFSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled", "WithinGroupLanguageASD_Q_UFSlower_Q_F")
  
  #1
  WithinGroupMotorTDC_MGSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled", "WithinGroupMotorTDC_MGSlower_Q_F")
  
  WithinGroupMotorTDC_Q_UFSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled", "WithinGroupMotorTDC_Q_UFSlower_Q_F")
  #2
  WithinGroupCognitionTDC_MGSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled", "WithinGroupCognitionTDC_MGSlower_Q_F")
  
  WithinGroupCognitionTDC_Q_UFSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled", "WithinGroupCognitionTDC_Q_UFSlower_Q_F")
  #3
  WithinGroupAwarenessTDC_MGSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled", "WithinGroupAwarenessTDC_MGSlower_Q_F")
  
  WithinGroupAwarenessTDC_Q_UFSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled", "WithinGroupAwarenessTDC_Q_UFSlower_Q_F")
  #4
  WithinGroupMotivationTDC_MGSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled", "WithinGroupMotivationTDC_MGSlower_Q_F")
  
  WithinGroupMotivationTDC_Q_UFSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled", "WithinGroupMotivationTDC_Q_UFSlower_Q_F")
  #5
  WithinGroupLanguageTDC_MGSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled", "WithinGroupLanguageTDC_MGSlower_Q_F")
  
  WithinGroupLanguageTDC_Q_UFSlower_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled", "WithinGroupLanguageTDC_Q_UFSlower_Q_F")
  
  MotorASD_overall <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled)/3 < 0",  "MotorASD_overall")
  
  MotorTDC_TG_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled < 0", "MotorTDC_TG_F")
  MotorTDC_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled < 0", "MotorTDC_Q_F")
  MotorTDC_Q_UF <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled < 0", "MotorTDC_Q_UF")
  MotorTDC_overall <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled)/3 < 0", "MotorTDC_overall")
  
  CognitionASD_TG_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled < 0", "CognitionASD_TG_F")
  CognitionASD_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled < 0", "CognitionASD_Q_F")
  CognitionASD_Q_UF <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled < 0","CognitionASD_Q_UF")
  CognitionASD_overall <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled)/3 < 0","CognitionASD_overall")
  
  CognitionTDC_TG_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled < 0","CognitionTDC_TG_F")
  CognitionTDC_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled < 0","CognitionTDC_Q_F")
  CognitionTDC_Q_UF <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled < 0","CognitionTDC_Q_UF")
  CognitionTDC_overall <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled)/3 < 0","CognitionTDC_overall")
  
  AwarenessASD_TG_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled > 0","AwarenessASD_TG_F")
  AwarenessASD_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled > 0","AwarenessASD_Q_F")
  AwarenessASD_Q_UF <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled > 0", "AwarenessASD_Q_UF")
  AwarenessASD_overall <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled)/3 > 0", "AwarenessASD_overall")
  
  AwarenessTDC_TG_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled > 0", "AwarenessTDC_TG_F")
  AwarenessTDC_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled > 0","AwarenessTDC_Q_F")
  AwarenessTDC_Q_UF <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled > 0","AwarenessTDC_Q_UF")
  AwarenessTDC_overall <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled)/3 > 0","AwarenessTDC_overall")
  
  MotivationASD_TG_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled > 0","MotivationASD_TG_F")
  MotivationASD_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled > 0","MotivationASD_Q_F")
  MotivationASD_Q_UF <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled > 0","MotivationASD_Q_UF")
  MotivationASD_overall <-extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled)/3 > 0","MotivationASD_overall")
  
  MotivationTDC_TG_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled < 0","MotivationTDC_TG_F")
  MotivationTDC_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled < 0","MotivationTDC_Q_F")
  MotivationTDC_Q_UF <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled < 0","MotivationTDC_Q_UF")
  MotivationTDC_overall <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled)/3 > 0","MotivationTDC_overall")
  
  
  LanguageASD_TG_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                            "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled > 0",
                                            "LanguageASD_TG_F")
  
  LanguageASD_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled > 0",
                                           "LanguageASD_Q_F")
  
  LanguageASD_Q_UF <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                            "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > 0",
                                            "LanguageASD_Q_UF")
  
  LanguageASD_overall <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                               "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)/3 > 0",
                                               "LanguageASD_overall")
  
  LanguageTDC_TG_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                            "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled > 0",
                                            "LanguageTDC_TG_F")
  
  LanguageTDC_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                           "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled > 0",
                                           "LanguageTDC_Q_F")
  
  LanguageTDC_Q_UF <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                            "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > 0",
                                            "LanguageTDC_Q_UF")
  
  LanguageTDC_overall <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                               "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)/3 > 0",
                                               "LanguageTDC_overall")
  
  
  LanguageOverallASDSlowerTDC <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                       "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)/3 > (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)/3",
                                                       "LanguageOverallASDSlowerTDC")
  
  MotorOverallASDSlowerTDC <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled)/3 > 
                                          (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled)/3",  
                                                    "MotorOverallASDFasterTDC")
  
  CognitionOverallASDSlowerTDC <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled)/3 >
                                              (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled)/3",
                                                        "CognitionOverallASDSlowerTDC")
  
  AwarenessOverallASDFasterTDC <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled)/3 < (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled)/3", "AwarenessOverallASDSlowerTDC")
  
  MotivationOverallASDFasterTDC <-extract_relevant_data(ChildLatency_m2_ID_NoOverlaps, "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled)/3 > (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled)/3","MotivationOverallASDFasterTDC")
  
  LanguageOverallASDSlowerTDC_MG <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                          "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled > DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled",
                                                          "LanguageOverallASDSlowerTDC_MG")
  
  LanguageOverallASDSlowerTDC_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled",
                                                           "LanguageOverallASDSlowerTDC_Q_F")
  
  LanguageOverallASDSlowerTDC_Q_UF <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                            "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled",
                                                            "LanguageOverallASDSlowerTDC_Q_UF")
  
  
  MotorOverallASDSlowerTDC_MG <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                       "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled > DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled",
                                                       "MotorOverallASDSlowerTDC_MG")
  
  MotorOverallASDSlowerTDC_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                        "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled",
                                                        "MotorOverallASDSlowerTDC_Q_F")
  
  MotorOverallASDSlowerTDC_Q_UF <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                         "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled",
                                                         "MotorOverallASDSlowerTDC_Q_UF")
  
  CognitionOverallASDSlowerTDC_MG <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                           "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled > DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled",
                                                           "CognitionOverallASDSlowerTDC_MG")
  
  CognitionOverallASDSlowerTDC_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                            "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled",
                                                            "CognitionOverallASDSlowerTDC_Q_F")
  
  CognitionOverallASDSlowerTDC_Q_UF <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                             "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled",
                                                             "CognitionOverallASDSlowerTDC_Q_UF")
  
  AwarenessOverallASDSlowerTDC_MG <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                           "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled < DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled",
                                                           "AwarenessOverallASDSlowerTDC_MG")
  
  AwarenessOverallASDSlowerTDC_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                            "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled",
                                                            "AwarenessOverallASDSlowerTDC_Q_F")
  
  AwarenessOverallASDSlowerTDC_Q_UF <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                             "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled",
                                                             "AwarenessOverallASDSlowerTDC_Q_UF")
  
  MotivationOverallASDSlowerTDC_MG <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                            "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled > DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled",
                                                            "MotivationOverallASDSlowerTDC_MG")
  
  MotivationOverallASDSlowerTDC_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                             "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled",
                                                             "MotivationOverallASDSlowerTDC_Q_F")
  
  MotivationOverallASDSlowerTDC_Q_UF <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                              "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled",
                                                              "MotivationOverallASDSlowerTDC_Q_UF")
  
  
  DifferenceInDifferenceLanguageMGSlowerQ_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                                     "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled)",
                                                                     "DifferenceInDifferenceLanguageMGSlowerQ_F")
  
  DifferenceInDifferenceLanguageUFSlowerQ_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                                     "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled)",
                                                                     "DifferenceInDifferenceLanguageUFSlowerQ_F")
  
  DifferenceInDifferenceMotorMGSlowerQ_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                                  "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotorS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled)",
                                                                  "DifferenceInDifferenceMotorMGSlowerQ_F")
  
  DifferenceInDifferenceMotorUFSlowerQ_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                                  "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotorS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotorS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotorS_scaled)",
                                                                  "DifferenceInDifferenceMotorUFSlowerQ_F")
  
  DifferenceInDifferenceCognitionMGSlowerQ_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                                      "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled)",
                                                                      "DifferenceInDifferenceCognitionMGSlowerQ_F")
  
  DifferenceInDifferenceCognitionUFSlowerQ_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                                      "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled)",
                                                                      "DifferenceInDifferenceCognitionUFSlowerQ_F")
  
  DifferenceInDifferenceAwarenessMGSlowerQ_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                                      "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:AwarenessS_scaled)
                      < (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled)",
                                                                      "DifferenceInDifferenceAwarenessMGSlowerQ_F")
  
  DifferenceInDifferenceAwarenessUFSlowerQ_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                                      "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:AwarenessS_scaled)
                      < (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:AwarenessS_scaled)",
                                                                      "DifferenceInDifferenceAwarenessUFSlowerQ_F")
  
  DifferenceInDifferenceMotivationMGSlowerQ_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                                       "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:MotivationS_scaled)
                      > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled)",
                                                                       "DifferenceInDifferenceMotivationMGSlowerQ_F")
  
  DifferenceInDifferenceMotivationUFSlowerQ_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                                       "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:MotivationS_scaled)
                      < (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:MotivationS_scaled)",
                                                                       "DifferenceInDifferenceMotivationUFSlowerQ_F")
  
  LanguageOverallASDSlowerTDC_Q_F <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled",
                                                           "LanguageOverallASDSlowerTDC_Q_F")
  
  LanguageOverallASDSlowerTDC_Q_UF <- extract_relevant_data(ChildLatency_m2_ID_NoOverlaps,
                                                            "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled > DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled",
                                                            "LanguageOverallASDSlowerTDC_Q_UF")
  
  hypothesis_data_full <- rbind(
    MotorASD_TG_F, MotorASD_Q_F, MotorASD_Q_UF, MotorASD_overall,
    MotorTDC_TG_F, MotorTDC_Q_F, MotorTDC_Q_UF, MotorTDC_overall,
    
    CognitionASD_TG_F, CognitionASD_Q_F, CognitionASD_Q_UF, CognitionASD_overall,
    CognitionTDC_TG_F, CognitionTDC_Q_F, CognitionTDC_Q_UF, CognitionTDC_overall,
    
    LanguageASD_TG_F, LanguageASD_Q_F, LanguageASD_Q_UF, LanguageASD_overall,
    LanguageTDC_TG_F, LanguageTDC_Q_F, LanguageTDC_Q_UF, LanguageTDC_overall,
    
    AwarenessASD_TG_F, AwarenessASD_Q_F, AwarenessASD_Q_UF, AwarenessASD_overall,
    AwarenessTDC_TG_F, AwarenessTDC_Q_F, AwarenessTDC_Q_UF, AwarenessTDC_overall,
    
    MotivationASD_TG_F, MotivationASD_Q_F, MotivationASD_Q_UF, MotivationASD_overall,
    MotivationTDC_TG_F, MotivationTDC_Q_F, MotivationTDC_Q_UF, MotivationTDC_overall,
    
    LanguageOverallASDSlowerTDC, MotorOverallASDSlowerTDC, CognitionOverallASDSlowerTDC, AwarenessOverallASDFasterTDC, MotivationOverallASDFasterTDC,
    
    LanguageOverallASDSlowerTDC_MG, LanguageOverallASDSlowerTDC_Q_F, LanguageOverallASDSlowerTDC_Q_UF,
    CognitionOverallASDSlowerTDC_MG, CognitionOverallASDSlowerTDC_Q_F, CognitionOverallASDSlowerTDC_Q_UF,
    MotorOverallASDSlowerTDC_MG, MotorOverallASDSlowerTDC_Q_F, MotorOverallASDSlowerTDC_Q_UF,
    AwarenessOverallASDSlowerTDC_MG, AwarenessOverallASDSlowerTDC_Q_F, AwarenessOverallASDSlowerTDC_Q_UF,
    MotivationOverallASDSlowerTDC_MG, MotivationOverallASDSlowerTDC_Q_F, MotivationOverallASDSlowerTDC_Q_UF,
    
    DifferenceInDifferenceLanguageMGSlowerQ_F, DifferenceInDifferenceLanguageUFSlowerQ_F,
    DifferenceInDifferenceMotorMGSlowerQ_F, DifferenceInDifferenceMotorUFSlowerQ_F,
    DifferenceInDifferenceMotivationMGSlowerQ_F, DifferenceInDifferenceMotivationUFSlowerQ_F,
    DifferenceInDifferenceCognitionMGSlowerQ_F, DifferenceInDifferenceCognitionUFSlowerQ_F,
    DifferenceInDifferenceAwarenessMGSlowerQ_F, DifferenceInDifferenceAwarenessUFSlowerQ_F,
    
    WithinGroupLanguageASD_MGSlower_Q_F, WithinGroupLanguageASD_Q_UFSlower_Q_F,
    WithinGroupCognitionASD_MGSlower_Q_F, WithinGroupCognitionASD_Q_UFSlower_Q_F,
    WithinGroupMotorASD_MGSlower_Q_F, WithinGroupMotorASD_Q_UFSlower_Q_F,
    WithinGroupAwarenessASD_MGSlower_Q_F, WithinGroupAwarenessASD_Q_UFSlower_Q_F,
    WithinGroupMotivationASD_MGSlower_Q_F, WithinGroupMotivationASD_Q_UFSlower_Q_F,
    
    WithinGroupLanguageTDC_MGSlower_Q_F, WithinGroupLanguageTDC_Q_UFSlower_Q_F,
    WithinGroupCognitionTDC_MGSlower_Q_F, WithinGroupCognitionTDC_Q_UFSlower_Q_F,
    WithinGroupMotorTDC_MGSlower_Q_F, WithinGroupMotorTDC_Q_UFSlower_Q_F,
    WithinGroupAwarenessTDC_MGSlower_Q_F, WithinGroupAwarenessTDC_Q_UFSlower_Q_F,
    WithinGroupMotivationTDC_MGSlower_Q_F, WithinGroupMotivationTDC_Q_UFSlower_Q_F
    
    
  ) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  return(hypothesis_data_full)
}

extract_all_hypothesis_tests_m1_child_gender <- function(model) {
  # Overall Estimates in Autism Group
  OverallLatencyThreeConditionsASD_Male <- extract_relevant_data(ChildLatency_m1_Gender,
                                                                 "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderMale + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderMale +DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderMale)/3 > 0",
                                                                 "OverallLatencyThreeConditionsASD_Male")
  
  OverallLatencyThreeConditionsASD_Female <- extract_relevant_data(ChildLatency_m1_Gender,
                                                                   "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderFemale + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderFemale + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderFemale)/3 > 0",
                                                                   "OverallLatencyThreeConditionsASD_Female")
  
  OverallLatencyThreeConditionsASD_NonBinary <- extract_relevant_data(ChildLatency_m1_Gender,
                                                                      "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderNonBinary + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderNonBinary + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderNonBinary)/3 > 0",
                                                                      "OverallLatencyThreeConditionsASD_NonBinary")
  
  MGOverallLatencyASD_Male <- extract_relevant_data(ChildLatency_m1_Gender,
                                                    "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderMale > 0",
                                                    "MGOverallLatencyASD_Male>0")
  
  MGOverallLatencyASD_Female <- extract_relevant_data(ChildLatency_m1_Gender,
                                                      "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderFemale > 0",
                                                      "MGOverallLatencyASD_Female>0")
  
  MGOverallLatencyASD_NonBinary <- extract_relevant_data(ChildLatency_m1_Gender,
                                                         "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderNonBinary > 0",
                                                         "MGOverallLatencyASD_NonBinary>0")
  
  Q_F_OverallLatencyASD_Male <- extract_relevant_data(ChildLatency_m1_Gender,
                                                      "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderMale > 0",
                                                      "Q_F_OverallLatencyASD_Male>0")
  Q_F_OverallLatencyASD_Female <- extract_relevant_data(ChildLatency_m1_Gender,
                                                        "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderFemale > 0",
                                                        "Q_F_OverallLatencyASD_Female>0")
  Q_F_OverallLatencyASD_NonBinary <- extract_relevant_data(ChildLatency_m1_Gender,
                                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderNonBinary > 0",
                                                           "Q_F_OverallLatencyASD_NonBinary>0")
  
  Q_UF_OverallLatencyASD_Male <- extract_relevant_data(ChildLatency_m1_Gender,
                                                       "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderMale > 0",
                                                       "Q_UF_OverallLatencyASD_Male>0")
  Q_UF_OverallLatencyASD_Female <- extract_relevant_data(ChildLatency_m1_Gender,
                                                         "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderFemale > 0",
                                                         "Q_UF_OverallLatencyASD_Female>0")
  Q_UF_OverallLatencyASD_NonBinary <- extract_relevant_data(ChildLatency_m1_Gender,
                                                            "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderNonBinary > 0",
                                                            "Q_UF_OverallLatencyASD_NonBinary>0")
  
  # Overall Estimates in TD Group
  OverallLatencyThreeConditionsTD_Male <- extract_relevant_data(ChildLatency_m1_Gender,
                                                                "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderMale + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderMale + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderMale) / 3 > 0",
                                                                "OverallLatencyThreeConditionsTD_Male")
  OverallLatencyThreeConditionsTD_Female <- extract_relevant_data(ChildLatency_m1_Gender,
                                                                  "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderFemale + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderFemale + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderFemale) / 3 > 0",
                                                                  "OverallLatencyThreeConditionsTD_Female")
  OverallLatencyThreeConditionsTD_NonBinary <- extract_relevant_data(ChildLatency_m1_Gender,
                                                                     "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderNonBinary + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderNonBinary + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderNonBinary) / 3 > 0",
                                                                     "OverallLatencyThreeConditionsTD_NonBinary")
  
  MGOverallLatencyTDC_Male <- extract_relevant_data(ChildLatency_m1_Gender,
                                                    "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderMale > 0",
                                                    "MGOverallLatencyTDC_Male>0")
  MGOverallLatencyTDC_Female <- extract_relevant_data(ChildLatency_m1_Gender,
                                                      "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderFemale > 0",
                                                      "MGOverallLatencyTDC_Female>0")
  MGOverallLatencyTDC_NonBinary <- extract_relevant_data(ChildLatency_m1_Gender,
                                                         "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderNonBinary > 0",
                                                         "MGOverallLatencyTDC_NonBinary>0")
  
  Q_F_OverallLatencyTDC_Male <- extract_relevant_data(ChildLatency_m1_Gender,
                                                      "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderMale > 0",
                                                      "Q_F_OverallLatencyTDC_Male>0")
  Q_F_OverallLatencyTDC_Female <- extract_relevant_data(ChildLatency_m1_Gender,
                                                        "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderFemale > 0",
                                                        "Q_F_OverallLatencyTDC_Female>0")
  Q_F_OverallLatencyTDC_NonBinary <- extract_relevant_data(ChildLatency_m1_Gender,
                                                           "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderNonBinary > 0",
                                                           "Q_F_OverallLatencyTDC_NonBinary>0")
  
  Q_UF_OverallLatencyTDC_Male <- extract_relevant_data(ChildLatency_m1_Gender,
                                                       "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderMale > 0",
                                                       "Q_UF_OverallLatencyTDC_Male>0")
  Q_UF_OverallLatencyTDC_Female <- extract_relevant_data(ChildLatency_m1_Gender,
                                                         "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderFemale > 0",
                                                         "Q_UF_OverallLatencyTDC_Female>0")
  Q_UF_OverallLatencyTDC_NonBinary <- extract_relevant_data(ChildLatency_m1_Gender,
                                                            "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderNonBinary > 0",
                                                            "Q_UF_OverallLatencyTDC_NonBinary>0")
  
  # Overall Differences between Autism Group and Typical Development Group
  OverallLatencyThreeConditionsASDMaleSlowerASDFemale <- extract_relevant_data(ChildLatency_m1_Gender,
                                                                               "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderMale + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderMale + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderMale)/3 
                                >
                                (DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderFemale + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderFemale +DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderFemale) / 3",
                                                                               "OverallLatencyThreeConditionsASDMaleSlowerASDFemale")
  
  # Overall Differences between Autism Group and Typical Development Group
  ASDMaleSlowerASDFemale_MG <- extract_relevant_data(ChildLatency_m1_Gender,
                                                     "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderMale > DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderFemale",
                                                     "ASDMaleSlowerASDFemale_MG")
  
  ASDMaleSlowerASDFemale_Q_F <- extract_relevant_data(ChildLatency_m1_Gender,
                                                      "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderMale > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderFemale",
                                                      "ASDMaleSlowerASDFemale_Q_F")
  
  ASDMaleSlowerASDFemale_Q_UF <- extract_relevant_data(ChildLatency_m1_Gender,
                                                       "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderMale > DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderFemale",
                                                       "ASDMaleSlowerASDFemale_Q_UF")
  
  OverallLatencyThreeConditionsTDCMaleSlowerTDCFemale <- extract_relevant_data(ChildLatency_m1_Gender,
                                                                               "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderMale + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderMale + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderMale)/3 
                                >
                                (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderFemale + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderFemale + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderFemale) / 3",
                                                                               "OverallLatencyThreeConditionsTDCMaleSlowerTDCFemale")
  
  TDCMaleSlowerTDCFemale_MG <- extract_relevant_data(ChildLatency_m1_Gender,
                                                     "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderMale > DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderFemale",
                                                     "TDCMaleSlowerTDCFemale_MG")
  
  TDCMaleSlowerTDCFemale_Q_F <- extract_relevant_data(ChildLatency_m1_Gender,
                                                      "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderMale > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderFemale",
                                                      "TDCMaleSlowerTDCFemale_Q_F")
  
  TDCMaleSlowerTDCFemale_Q_UF <- extract_relevant_data(ChildLatency_m1_Gender,
                                                       "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderMale < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderFemale",
                                                       "TDCMaleSlowerTDCFemale_Q_UF")
  
  OverallLatencyThreeConditionsASDMaleSlowerASDNonBinary <- extract_relevant_data(ChildLatency_m1_Gender,
                                                                                  "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderMale + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderMale + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderMale)/3 
                                >
                                (DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderNonBinary + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderNonBinary +DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderNonBinary) / 3",
                                                                                  "OverallLatencyThreeConditionsASDMaleSlowerASDNonBinary")
  
  ASDMaleSlowerASDNonBinary_MG <- extract_relevant_data(ChildLatency_m1_Gender,
                                                        "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderMale > DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderNonBinary",
                                                        "ASDMaleSlowerASDNonBinary_MG")
  
  ASDMaleSlowerASDNonBinary_Q_F <- extract_relevant_data(ChildLatency_m1_Gender,
                                                         "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderMale > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderNonBinary",
                                                         "ASDMaleSlowerASDNonBinary_Q_F")
  
  ASDMaleSlowerASDNonBinary_Q_UF <- extract_relevant_data(ChildLatency_m1_Gender,
                                                          "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderMale > DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderNonBinary",
                                                          "ASDMaleSlowerASDNonBinary_Q_UF")
  
  OverallLatencyThreeConditionsASDFemaleSlowerASDNonBinary <- extract_relevant_data(ChildLatency_m1_Gender,
                                                                                    "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderFemale + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderFemale + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderFemale)/3 
                                >
                                (DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderNonBinary + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderNonBinary +DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderNonBinary) / 3",
                                                                                    "OverallLatencyThreeConditionsASDFemaleSlowerASDNonBinary")
  
  ASDFemaleSlowerASDNonBinary_MG <- extract_relevant_data(ChildLatency_m1_Gender,
                                                          "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderFemale > DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderNonBinary",
                                                          "ASDFemaleSlowerASDNonBinary_MG")
  
  ASDFemaleSlowerASDNonBinary_Q_F <- extract_relevant_data(ChildLatency_m1_Gender,
                                                           "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderFemale > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderNonBinary",
                                                           "ASDFemaleSlowerASDNonBinary_Q_F")
  
  ASDFemaleSlowerASDNonBinary_MG_Q_UF <- extract_relevant_data(ChildLatency_m1_Gender,
                                                               "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderFemale > DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderNonBinary",
                                                               "ASDFemaleSlowerASDNonBinary_MG_Q_UF")
  
  OverallLatencyThreeConditionsASDMaleSlowerTDCMale <- extract_relevant_data(ChildLatency_m1_Gender,
                                                                             "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderMale + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderMale + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderMale)/3 
                                <
                                (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderMale + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderMale + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderMale) / 3",
                                                                             "OverallLatencyThreeConditionsASDMaleSlowerTDCMale")
  
  ASDMaleSlowerTDCMale_MG <- extract_relevant_data(ChildLatency_m1_Gender,
                                                   "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderMale < DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderMale",
                                                   "ASDMaleSlowerTDCMale_MG")
  
  ASDMaleSlowerTDCMale_Q_F <- extract_relevant_data(ChildLatency_m1_Gender,
                                                    "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderMale < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderMale",
                                                    "ASDMaleSlowerTDCMale_Q_F")
  
  ASDMaleSlowerTDCMale_Q_UF <- extract_relevant_data(ChildLatency_m1_Gender,
                                                     "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderMale < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderMale",
                                                     "ASDMaleSlowerTDCMale_Q_UF")
  
  
  OverallLatencyThreeConditionsASDFemaleSlowerTDCFemale <- extract_relevant_data(ChildLatency_m1_Gender,
                                                                                 "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderFemale + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderFemale + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderFemale)/3 
                                <
                                (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderFemale + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderFemale + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderFemale) / 3",
                                                                                 "OverallLatencyThreeConditionsASDFemaleSlowerTDCFemale")
  
  ASDFemaleSlowerTDCFemale_MG <- extract_relevant_data(ChildLatency_m1_Gender,
                                                       "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderFemale < DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderFemale",
                                                       "ASDFemaleSlowerTDCFemale_MG")
  
  ASDFemaleSlowerTDCFemale_Q_F <- extract_relevant_data(ChildLatency_m1_Gender,
                                                        "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderFemale < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderFemale",
                                                        "ASDFemaleSlowerTDCFemale_Q_F")
  
  ASDFemaleSlowerTDCFemale_Q_UF <- extract_relevant_data(ChildLatency_m1_Gender,
                                                         "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderFemale < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderFemale",
                                                         "ASDFemaleSlowerTDCFemale_Q_UF")
  
  
  ASDMaleSlowerTDCFemale_MG <- extract_relevant_data(ChildLatency_m1_Gender,
                                                     "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderMale < DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderFemale",
                                                     "ASDMaleSlowerTDCFemale_MG")
  
  ASDMaleSlowerTDCFemale_Q_F <- extract_relevant_data(ChildLatency_m1_Gender,
                                                      "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderMale < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderFemale",
                                                      "ASDMaleSlowerTDCFemale_Q_F")
  
  ASDMaleSlowerTDCFemale_Q_UF <- extract_relevant_data(ChildLatency_m1_Gender,
                                                       "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderMale < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderFemale",
                                                       "ASDMaleSlowerTDCFemale_Q_UF")
  
  ASDFemaleSlowerTDCMale_MG <- extract_relevant_data(ChildLatency_m1_Gender,
                                                     "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderFemale < DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderMale",
                                                     "ASDFemaleSlowerTDCMale_MG")
  
  ASDFemaleSlowerTDCMale_Q_F <- extract_relevant_data(ChildLatency_m1_Gender,
                                                      "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderFemale < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderMale",
                                                      "ASDFemaleSlowerTDCMale_Q_F")
  
  ASDFemaleSlowerTDCMale_Q_UF <- extract_relevant_data(ChildLatency_m1_Gender,
                                                       "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderFemale < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderMale",
                                                       "ASDFemaleSlowerTDCMale_Q_UF")
  
  ASDMaleSlowerTDCMale_AcrossConditions <- extract_relevant_data(ChildLatency_m1_Gender,
                                                                 "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderMale + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderMale + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderMale) / 3
                      < (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderMale + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderMale + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderMale) / 3",
                                                                 "ASDMaleSlowerTDCMale_AcrossConditions")
  
  ASDFemaleSlowerTDCFemale_AcrossConditions <- extract_relevant_data(ChildLatency_m1_Gender,
                                                                     "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderFemale + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderFemale + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderFemale) / 3
                      < (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderFemale + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderFemale + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderFemale) / 3",
                                                                     "ASDFemaleSlowerTDCFemale_AcrossConditions")
  
  ASDMaleSlowerTDCFemale_AcrossConditions <- extract_relevant_data(ChildLatency_m1_Gender,
                                                                   "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderMale + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderMale + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderMale) / 3
                      < (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderFemale + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderFemale + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderFemale) / 3",
                                                                   "ASDMaleSlowerTDCFemale_AcrossConditions")
  
  ASDFemaleSlowerTDCMale_AcrossConditions <- extract_relevant_data(ChildLatency_m1_Gender,
                                                                   "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderFemale + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderFemale + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderFemale) / 3
                      < (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderMale + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderMale + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderMale) / 3",
                                                                   "ASDFemaleSlowerTDCMale_AcrossConditions")
  
  ASDNonBinSlowerTDCMale_MG <- extract_relevant_data(ChildLatency_m1_Gender,
                                                     "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderNonBinary < DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderMale",
                                                     "ASDNonBinSlowerTDCMale_MG")
  
  ASDNonBinSlowerTDCFemale_MG <- extract_relevant_data(ChildLatency_m1_Gender,
                                                       "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderNonBinary < DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderFemale",
                                                       "ASDNonBinSlowerTDCFemale_MG")
  
  ASDNonBinSlowerTDCMale_Q_F <- extract_relevant_data(ChildLatency_m1_Gender,
                                                      "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderNonBinary < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderMale",
                                                      "ASDNonBinSlowerTDCMale_Q_F")
  
  ASDNonBinSlowerTDCFemale_Q_F <- extract_relevant_data(ChildLatency_m1_Gender,
                                                        "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderNonBinary < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderFemale",
                                                        "ASDNonBinSlowerTDCFemale_Q_F")
  
  ASDNonBinSlowerTDCMale_Q_UF <- extract_relevant_data(ChildLatency_m1_Gender,
                                                       "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderNonBinary < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderMale",
                                                       "ASDNonBinSlowerTDCMale_Q_UF")
  
  ASDNonBinSlowerTDCFemale_Q_UF <- extract_relevant_data(ChildLatency_m1_Gender,
                                                         "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderNonBinary < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderFemale",
                                                         "ASDNonBinSlowerTDCFemale_Q_UF")
  
  
  
  ASDNonBinSlowerTDCMale_Overall <- extract_relevant_data(ChildLatency_m1_Gender,
                                                          "(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderNonBinary + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderNonBinary + DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderNonBinary) / 3 < 
                      (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderMale + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderMale + DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderMale) / 3",
                                                          "ASDNonBinSlowerTDCMale_Overall")
  
  ASDNonBinSlowerTDCFemale_Overall <- extract_relevant_data(ChildLatency_m1_Gender,
                                                            "(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderNonBinary + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderNonBinary + DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderNonBinary) / 3 < 
                      (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderFemale + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderFemale + DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderFemale) / 3",
                                                            "ASDNonBinSlowerTDCFemale_Overall")
  
  # Overall values for Beta parameter
  BetaParameterASD_MG_Male <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                               "beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderMale > 0",
                                                               "BetaParameterASD_MG_Male")
  
  BetaParameterASD_MG_Female <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                 "beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderFemale > 0",
                                                                 "BetaParameterASD_MG_Female")
  
  BetaParameterASD_MG_NonBinary <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                    "beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderNonBinary > 0",
                                                                    "BetaParameterASD_MG_NonBinary")
  
  BetaParameterASD_Q_F_Male <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                "beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderMale < 0",
                                                                "BetaParameterASD_Q_F_Male")
  
  BetaParameterASD_Q_F_Female <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                  "beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderFemale < 0",
                                                                  "BetaParameterASD_Q_F_Female")
  
  BetaParameterASD_Q_F_NonBinary <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                     "beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderNonBinary < 0",
                                                                     "BetaParameterASD_Q_F_Nonbinary")
  
  BetaParameterASD_Q_UF_Male <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                 "beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderMale < 0",
                                                                 "BetaParameterASD_Q_UF_Male")
  
  BetaParameterASD_Q_UF_Female <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                   "beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderFemale < 0",
                                                                   "BetaParameterASD_Q_UF_Female")
  
  BetaParameterASD_Q_UF_NonBinary <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                      "beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderNonBinary < 0",
                                                                      "BetaParameterASD_Q_UF_NonBinary")
  
  BetaParameterASD_Overall_Male <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                    "(beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderMale +
                 beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderMale +
                 beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderMale) / 3 < 0",
                                                                    "BetaParameterASD_Overall_Male")
  
  BetaParameterASD_Overall_Female <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                      "(beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderFemale +
                 beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderFemale +
                 beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderFemale) / 3 < 0",
                                                                      "BetaParameterASD_Overall_Female")
  
  BetaParameterASD_Overall_NonBinary <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                         "(beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:GenderNonBinary +
                 beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:GenderNonBinary +
                 beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:GenderNonBinary) / 3 < 0",
                                                                         "BetaParameterASD_Overall_NonBinary")
  
  BetaParameterTDC_MG_Male <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                               "beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderMale > 0",
                                                               "BetaParameterTDC_MG_Male")
  
  BetaParameterTDC_MG_Female <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                 "beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderFemale > 0",
                                                                 "BetaParameterTDC_MG_Female")
  
  BetaParameterTDC_MG_NonBinary <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                    "beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderNonBinary > 0",
                                                                    "BetaParameterTDC_MG_NonBinary")
  
  BetaParameterTDC_Q_F_Male <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                "beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderMale < 0",
                                                                "BetaParameterTDC_Q_F_Male")
  
  BetaParameterTDC_Q_F_Female <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                  "beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderFemale < 0",
                                                                  "BetaParameterTDC_Q_F_Female")
  
  BetaParameterTDC_Q_F_NonBinary <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                     "beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderNonBinary < 0",
                                                                     "BetaParameterTDC_Q_F_Nonbinary")
  
  BetaParameterTDC_Q_UF_Male <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                 "beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderMale < 0",
                                                                 "BetaParameterTDC_Q_UF_Male")
  
  BetaParameterTDC_Q_UF_Female <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                   "beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderFemale < 0",
                                                                   "BetaParameterTDC_Q_UF_Female")
  
  BetaParameterTDC_Q_UF_NonBinary <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                      "beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderNonBinary < 0",
                                                                      "BetaParameterTDC_Q_UF_NonBinary")
  
  BetaParameterTDC_Overall_Male <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                    "(beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderMale +
                 beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderMale +
                 beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderMale) / 3 < 0",
                                                                    "BetaParameterTDC_Overall")
  
  BetaParameterTDC_Overall_Female <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                      "(beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderFemale +
                 beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderFemale +
                 beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderFemale) / 3 < 0",
                                                                      "BetaParameterTDC_Overall")
  
  BetaParameterTDC_Overall_NonBinary <- extract_relevant_nonlatency_data(ChildLatency_m1_Gender,
                                                                         "(beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:GenderNonBinary +
                 beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:GenderNonBinary +
                 beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:GenderNonBinary) / 3 < 0",
                                                                         "BetaParameterTDC_Overall")
  
  
  HypothesisResultsChildLatenciesModel1_Gender <- rbind(OverallLatencyThreeConditionsASD_Male,
                                                        OverallLatencyThreeConditionsTD_Male,
                                                        MGOverallLatencyASD_Male,
                                                        Q_F_OverallLatencyASD_Male,
                                                        Q_UF_OverallLatencyASD_Male,
                                                        MGOverallLatencyTDC_Male,
                                                        Q_F_OverallLatencyTDC_Male,
                                                        Q_UF_OverallLatencyTDC_Male,
                                                        
                                                        OverallLatencyThreeConditionsASD_Female,
                                                        OverallLatencyThreeConditionsTD_Female,
                                                        MGOverallLatencyASD_Female,
                                                        Q_F_OverallLatencyASD_Female,
                                                        Q_UF_OverallLatencyASD_Female,
                                                        MGOverallLatencyTDC_Female,
                                                        Q_F_OverallLatencyTDC_Female,
                                                        Q_UF_OverallLatencyTDC_Female,
                                                        
                                                        OverallLatencyThreeConditionsASD_NonBinary,
                                                        OverallLatencyThreeConditionsTD_NonBinary,
                                                        MGOverallLatencyASD_NonBinary,
                                                        Q_F_OverallLatencyASD_NonBinary,
                                                        Q_UF_OverallLatencyASD_NonBinary,
                                                        MGOverallLatencyTDC_NonBinary,
                                                        Q_F_OverallLatencyTDC_NonBinary,
                                                        Q_UF_OverallLatencyTDC_NonBinary,
                                                        
                                                        OverallLatencyThreeConditionsASDMaleSlowerASDFemale,
                                                        OverallLatencyThreeConditionsTDCMaleSlowerTDCFemale,
                                                        OverallLatencyThreeConditionsASDMaleSlowerASDNonBinary,
                                                        OverallLatencyThreeConditionsASDFemaleSlowerASDNonBinary,
                                                        OverallLatencyThreeConditionsASDMaleSlowerTDCMale,
                                                        OverallLatencyThreeConditionsASDFemaleSlowerTDCFemale,
                                                        
                                                        ASDMaleSlowerASDFemale_MG, 
                                                        ASDMaleSlowerASDFemale_Q_F,
                                                        ASDMaleSlowerASDFemale_Q_UF,
                                                        TDCMaleSlowerTDCFemale_MG,
                                                        TDCMaleSlowerTDCFemale_Q_F,
                                                        TDCMaleSlowerTDCFemale_Q_UF,
                                                        ASDMaleSlowerASDNonBinary_MG,
                                                        ASDMaleSlowerASDNonBinary_Q_F,
                                                        ASDMaleSlowerASDNonBinary_Q_UF,
                                                        ASDFemaleSlowerASDNonBinary_MG,
                                                        ASDFemaleSlowerASDNonBinary_Q_F,
                                                        ASDFemaleSlowerASDNonBinary_MG_Q_UF,
                                                        ASDMaleSlowerTDCMale_MG,
                                                        ASDMaleSlowerTDCMale_Q_F,
                                                        ASDMaleSlowerTDCMale_Q_UF,
                                                        ASDFemaleSlowerTDCFemale_MG,
                                                        ASDFemaleSlowerTDCFemale_Q_F,
                                                        ASDFemaleSlowerTDCFemale_Q_UF,
                                                        ASDMaleSlowerTDCMale_AcrossConditions,
                                                        ASDFemaleSlowerTDCFemale_AcrossConditions,
                                                        ASDMaleSlowerTDCFemale_AcrossConditions,
                                                        ASDFemaleSlowerTDCMale_AcrossConditions,
                                                        ASDMaleSlowerTDCFemale_MG,
                                                        ASDMaleSlowerTDCFemale_Q_F,
                                                        ASDMaleSlowerTDCFemale_Q_UF,
                                                        ASDFemaleSlowerTDCMale_MG,
                                                        ASDFemaleSlowerTDCMale_Q_F,
                                                        ASDFemaleSlowerTDCMale_Q_UF,
                                                        ASDNonBinSlowerTDCMale_MG,
                                                        ASDNonBinSlowerTDCFemale_MG,
                                                        ASDNonBinSlowerTDCMale_Q_F,
                                                        ASDNonBinSlowerTDCFemale_Q_F,
                                                        ASDNonBinSlowerTDCMale_Q_UF,
                                                        ASDNonBinSlowerTDCFemale_Q_UF,
                                                        ASDNonBinSlowerTDCMale_Overall,
                                                        ASDNonBinSlowerTDCFemale_Overall,
                                                        
                                                        BetaParameterASD_MG_Male,
                                                        BetaParameterASD_Q_F_Male,
                                                        BetaParameterASD_Q_UF_Male,
                                                        BetaParameterASD_Overall_Male,
                                                        BetaParameterTDC_MG_Male,
                                                        BetaParameterTDC_Q_F_Male,
                                                        BetaParameterTDC_Q_UF_Male,
                                                        BetaParameterTDC_Overall_Male,
                                                        
                                                        BetaParameterASD_MG_Female,
                                                        BetaParameterASD_Q_F_Female,
                                                        BetaParameterASD_Q_UF_Female,
                                                        BetaParameterASD_Overall_Female,
                                                        BetaParameterTDC_MG_Female,
                                                        BetaParameterTDC_Q_F_Female,
                                                        BetaParameterTDC_Q_UF_Female,
                                                        BetaParameterTDC_Overall_Female,
                                                        
                                                        BetaParameterASD_MG_NonBinary,
                                                        BetaParameterASD_Q_F_NonBinary,
                                                        BetaParameterASD_Q_UF_NonBinary,
                                                        BetaParameterASD_Overall_NonBinary,
                                                        BetaParameterTDC_MG_NonBinary,
                                                        BetaParameterTDC_Q_F_NonBinary,
                                                        BetaParameterTDC_Q_UF_NonBinary,
                                                        BetaParameterTDC_Overall_NonBinary) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  return(HypothesisResultsChildLatenciesModel1_Gender)
}

extract_all_hypothesis_tests_m1_surrogate <- function(model) {
  # Overall Estimates in Autism Group
  OverallLatencyThreeConditionsASD <- extract_relevant_data(ChildLatency_m1_surrogate,
                                                            "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar +DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
                                                            "OverallLatencyThreeConditionsASD")
  
  MGOverallLatencyASD <- extract_relevant_data(ChildLatency_m1_surrogate,
                                               "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > 0",
                                               "MGOverallLatencyASD>0")
  
  Q_F_OverallLatencyASD <- extract_relevant_data(ChildLatency_m1_surrogate,
                                                 "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > 0",
                                                 "Q_F_OverallLatencyASD>0")
  
  Q_UF_OverallLatencyASD <- extract_relevant_data(ChildLatency_m1_surrogate,
                                                  "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                  "Q_UF_OverallLatencyASD>0")
  
  # Overall Estimates in TD Group
  OverallLatencyThreeConditionsTD <- extract_relevant_data(ChildLatency_m1_surrogate,
                                                           "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar +DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3 > 0",
                                                           "OverallLatencyThreeConditionsTD")
  
  MGOverallLatencyTDC <- extract_relevant_data(ChildLatency_m1_surrogate,
                                               "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > 0",
                                               "MGOverallLatencyTDC>0")
  
  Q_F_OverallLatencyTDC <- extract_relevant_data(ChildLatency_m1_surrogate,
                                                 "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar > 0",
                                                 "Q_F_OverallLatencyTDC>0")
  
  Q_UF_OverallLatencyTDC <- extract_relevant_data(ChildLatency_m1_surrogate,
                                                  "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > 0",
                                                  "Q_UF_OverallLatencyTDC>0")
  
  # Overall Differences between Autism Group and Typical Development Group
  OverallLatencyThreeConditionsASDFasterTD <- extract_relevant_data(ChildLatency_m1_surrogate,
                                                                    "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar +DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 
                                <
                                (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar +DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3",
                                                                    "OverallLatencyThreeConditionsASD<TD")
  
  
  # Overall values for Beta parameter
  # BetaParameterASD_MG <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate,
  #                  "beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > 0",
  #                  "BetaParameterASD_MG")
  # 
  # BetaParameterASD_Q_F <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate,
  #                  "beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < 0",
  #                  "BetaParameterASD_Q_F")
  # 
  # BetaParameterASD_Q_UF <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate,
  #                  "beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < 0",
  #                  "BetaParameterASD_Q_UF")
  # 
  # BetaParameterASD_Overall <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate,
  #                  "(beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar +
  #                  beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar +
  #                  beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar) / 3 < 0",
  #                  "BetaParameterASD_Overall")
  # 
  # # Differences, Typical Development Group with Long Pauses
  # BetaParameterTDC_MG <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate,
  #                  "beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > 0",
  #                  "BetaParameterTDC_MG")
  # 
  # BetaParameterTDC_Q_F <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate,
  #                  "beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar < 0",
  #                  "BetaParameterTDC_Q_F")
  # 
  # BetaParameterTDC_Q_UF <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate,
  #                  "beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar < 0",
  #                  "BetaParameterTDC_Q_UF")
  # 
  # BetaParameterTDC_Overall <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate,
  #                  "(beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar +
  #                  beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar +
  #                  beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3 < 0",
  #                  "BetaParameterTDC_Overall")
  # 
  # 
  # 
  # 
  # WithinGroupBetaParameterASD_MG_Slower_Q_F <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate,
  #                  "beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
  #                  "WithinGroupBetaParameterASD_MG_Slower_Q_F")
  # 
  # WithinGroupBetaParameterASD_Q_F_Slower_Q_UF <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate,
  #                  "beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
  #                  "WithinGroupBetaParameterASD_Q_F_Slower_Q_UF")
  # 
  # WithinGroupBetaParameterTDC_MG_Slower_Q_F <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate,
  #                  "beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
  #                  "WithinGroupBetaParameterTDC_MG_Slower_Q_F")
  # 
  # WithinGroupBetaParameterTDC_Q_F_Slower_Q_UF <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate,
  #                  "beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
  #                  "WithinGroupBetaParameterTDC_Q_F_Slower_Q_UF")
  # 
  # DifferenceInDifferenceBetaParameterBetweenGroupsMGBiggerQ_F <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate,
  #                  "(beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar - beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) > 
  #                  (beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar - beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar)",
  #                  "DifferenceInDifferenceBetaParameterBetweenGroupsMGBiggerQ_F")
  # 
  # DifferenceInDifferenceBetaParameterBetweenGroupsFamiliarityBiggerUnfamiliarity <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate,
  #                  "(beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar - beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar) < 
  #                  (beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar - beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)",
  #                  "DifferenceInDifferenceBetaParameterBetweenGroupsFamiliarityBiggerUnfamiliarity")
  # 
  # # Differences between Autism Group and Typical Development Group with Long Pauses
  # OverallBetaParameterASDLowerTD <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate,
  #                  "(beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar + beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar) / 3 < ((beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar + beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) / 3)",
  #                  "OverallBetaParameterASD<TD")
  # 
  # BetaParameterASDFasterTD_MG <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate,
  #                       "beta_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar < beta_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar",
  #                       "BetaParameter_ASD_MG<TD_MG")
  # 
  # BetaParameterASDFasterTD_Q_F <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate,
  #                       "beta_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < beta_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
  #                       "BetaParameter_ASD_MG<TD_Q_F")
  # 
  # BetaParameterASDFasterTD_Q_UF <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate,
  #                       "beta_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < beta_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
  #                       "BetaParameter_ASD_MG<TD_Q_UF")
  
  # Differences within Autism Group across Conditions:
  WithinASDGroupMGSlowerQuestionsFamiliar <- extract_relevant_data(ChildLatency_m1_surrogate,
                                                                   "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                   "WithinASDGroupMG>QuestionsFamiliar")
  
  # Differences within TD Group across Conditions:
  WithinTDGroupMGSlowerQuestionsFamiliar <- extract_relevant_data(ChildLatency_m1_surrogate,
                                                                  "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                  "WithinTDGroupMG>QuestionsFamiliar")
  
  # Differences within Autism Group across Familiarity:
  WithinAutismGroupQuestionsUnfamiliarSlowerFamiliar <- extract_relevant_data(ChildLatency_m1_surrogate,
                                                                              "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > DiagnosisASD:TaskQuestions:FamiliarityFamiliar",
                                                                              "WithinAutismGroupQuestionsUnfamiliar>Familiar")
  
  # Differences within TD Group across Familiarity:
  WithinTDGroupQuestionsUnfamiliarSlowerFamiliar <- extract_relevant_data(ChildLatency_m1_surrogate,
                                                                          "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                          "WithinTDGroupQuestionsUnfamiliar>Familiar")
  
  # Differences between Autism Group and Typical Development Group Within MatchingGame
  WithinConditionsMatchingGameASDFasterTD <- extract_relevant_data(ChildLatency_m1_surrogate,
                                                                   "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar < DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar",
                                                                   "WithinConditionsMatchingGameASD<TD")
  
  # Differences between Autism Group and Typical Development Group Within Familiar Questions:
  WithinConditionsQuestionsFamiliarASDFasterTD <- extract_relevant_data(ChildLatency_m1_surrogate,
                                                                        "DiagnosisASD:TaskQuestions:FamiliarityFamiliar < DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                        "WithinConditionsQuestionsFamiliarASD<TD")
  
  # Differences between Autism Group and Typical Development Group Within Unfamiliar Questions:
  WithinConditionsQuestionsUnfamiliarASDFasterTD <- extract_relevant_data(ChildLatency_m1_surrogate,
                                                                          "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                          "WithinConditionsQuestionsUnfamiliarASD<TD")
  
  # Differences between Autism Group and Typical Development Group in Differences within Familiarity
  DifferenceInDifferenceBetweenGroupsUnfamiliarBiggerFamiliarity <- extract_relevant_data(ChildLatency_m1_surrogate,
                                                                                          "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) < 
                      (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                          "DifferenceInDifferenceBetweenGroupsUnfamiliar>Familiarity")
  
  # Differences between Autism Group and Typical Development Group in Differences within Conditions
  DifferenceInDifferenceBetweenGroupsMGBiggerThanConversations <- extract_relevant_data(ChildLatency_m1_surrogate,
                                                                                        "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar) > 
                      (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                        "DifferenceInDifferenceBetweenGroupsMG>Conversations")
  
  # Overall Subject Variability in Autismm Group
  SubjectVariabilityASD <- extract_relevant_sd_data(ChildLatency_m1_surrogate, 
                                                    "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisASD + TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD)/3 > 0",
                                                    "SubjectVariabilityASD")
  
  SubjectVariabilityASD_MG <- extract_relevant_sd_data(ChildLatency_m1_surrogate, 
                                                       "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD > 0",
                                                       "SubjectVariabilityASD_MG")
  
  SubjectVariabilityASD_Q_F <- extract_relevant_sd_data(ChildLatency_m1_surrogate, 
                                                        "TaskQuestions:FamiliarityFamiliar:DiagnosisASD > 0",
                                                        "SubjectVariabilityASD_Q_F")
  
  SubjectVariabilityASD_Q_UF <- extract_relevant_sd_data(ChildLatency_m1_surrogate, 
                                                         "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD > 0",
                                                         "SubjectVariabilityASD_Q_UF")
  
  # Overall Subject Variability in TD Group
  SubjectVariabilityTD <- extract_relevant_sd_data(ChildLatency_m1_surrogate, 
                                                   "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisTDC + TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC)/3 > 0",
                                                   "SubjectVariabilityTD")
  
  SubjectVariabilityTD_MG <- extract_relevant_sd_data(ChildLatency_m1_surrogate, 
                                                      "TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC > 0",
                                                      "SubjectVariabilityTD_MG")
  
  SubjectVariabilityTD_Q_F <- extract_relevant_sd_data(ChildLatency_m1_surrogate, 
                                                       "TaskQuestions:FamiliarityFamiliar:DiagnosisTDC > 0",
                                                       "SubjectVariabilityTD_Q_F")
  
  SubjectVariabilityTD_Q_UF <- extract_relevant_sd_data(ChildLatency_m1_surrogate, 
                                                        "TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC > 0",
                                                        "SubjectVariabilityTD_Q_UF")
  
  # Overall Subject Variability in TD Group
  SubjectVariabilityASDSmallerThanTD <- extract_relevant_sd_data(ChildLatency_m1_surrogate,
                                                                 "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisASD + TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD) / 3 <
                                              (TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC + 
                                               TaskQuestions:FamiliarityFamiliar:DiagnosisTDC + TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC) / 3",
                                                                 "SubjectVariabilityASD<TD")
  
  
  SubjectVariabilityASDSmallerThanTD_MG <- extract_relevant_sd_data(ChildLatency_m1_surrogate,
                                                                    "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD  <
                                              TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC",
                                                                    "SubjectVariabilityASD_MG<TD_MG")
  
  SubjectVariabilityASDSmallerThanTD_Q_F <- extract_relevant_sd_data(ChildLatency_m1_surrogate,
                                                                     "TaskQuestions:FamiliarityFamiliar:DiagnosisASD  <
                                              TaskQuestions:FamiliarityFamiliar:DiagnosisTDC",
                                                                     "SubjectVariabilityASD_Q_F<TD_Q_F")
  
  SubjectVariabilityASDSmallerThanTD_Q_UF <- extract_relevant_sd_data(ChildLatency_m1_surrogate,
                                                                      "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD  <
                                              TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC",
                                                                      "SubjectVariabilityASD_Q_UF<TD_Q_UF")
  
  SubjectVariabilityASD_MG_SmallerThan_Q_F <- extract_relevant_sd_data(ChildLatency_m1_surrogate,
                                                                       "TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD > TaskQuestions:FamiliarityFamiliar:DiagnosisASD",
                                                                       "SubjectVariabilityASD_MG_SmallerThan_Q_F")
  
  SubjectVariabilityASD_Q_F_SmallerThan_Q_UF <- extract_relevant_sd_data(ChildLatency_m1_surrogate,
                                                                         "TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD > TaskQuestions:FamiliarityFamiliar:DiagnosisASD",
                                                                         "SubjectVariabilityASD_Q_F_SmallerThan_Q_UF")
  
  SubjectVariabilityTDC_MG_SmallerThan_Q_F <- extract_relevant_sd_data(ChildLatency_m1_surrogate,
                                                                       "TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC < TaskQuestions:FamiliarityFamiliar:DiagnosisTDC",
                                                                       "SubjectVariabilityTDC_MG_SmallerThan_Q_F")
  
  SubjectVariabilityTDC_Q_F_SmallerThan_Q_UF <- extract_relevant_sd_data(ChildLatency_m1_surrogate,
                                                                         "TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC > TaskQuestions:FamiliarityFamiliar:DiagnosisTDC",
                                                                         "SubjectVariabilityTDC_Q_F_SmallerThan_Q_UF")
  
  DifferenceInDifferenceSubjectVariabilityUnfamiliarSmallerFamiliarity <- extract_relevant_sd_data(ChildLatency_m1_surrogate,
                                                                                                   "(TaskQuestions:FamiliarityUnfamiliar:DiagnosisASD - TaskQuestions:FamiliarityUnfamiliar:DiagnosisTDC) 
                                              < 
                                              (TaskQuestions:FamiliarityFamiliar:DiagnosisASD - TaskQuestions:FamiliarityFamiliar:DiagnosisTDC)",
                                                                                                   "DifferenceInDifferenceSubjectVariabilityUnfamiliarBiggerFamiliarity")
  
  DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions <- extract_relevant_sd_data(ChildLatency_m1_surrogate,
                                                                                                  "(TaskMatchingGame:FamiliarityFamiliar:DiagnosisTDC - TaskMatchingGame:FamiliarityFamiliar:DiagnosisASD) 
                                              < 
                                              (TaskQuestions:FamiliarityFamiliar:DiagnosisTDC - TaskQuestions:FamiliarityFamiliar:DiagnosisASD)",
                                                                                                  "DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions")
  
  # # Overall Visit Variability across visits in Autism Group
  # VisitVariabilityASDOverall <- extract_relevant_sdvisit_data(ChildLatency_m1_surrogate, 
  #                                                "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
  #                                                "VisitVariabilityASDOverall")
  # 
  # # Overall Visit Variability across visits in TD Group
  # VisitVariabilityTDCOverall <- extract_relevant_sdvisit_data(ChildLatency_m1_surrogate, 
  #                                                "(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3 > 0",
  #                                                "VisitVariabilityTDCOverall")
  # 
  # 
  # VisitVariabilityASDOverallHigherTDSOverall <- extract_relevant_sdvisit_data(ChildLatency_m1_surrogate, 
  #                                                "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityFamiliar + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 >
  #                                                (DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3",
  #                                                "VisitVariabilityASDOverallLowerTDSOverall")
  # 
  # # Subject Variability across visits MG in Autism Group
  # VisitVariabilityASDMG <- extract_relevant_sdvisit_data(ChildLatency_m1_surrogate, 
  #                                                "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > 0",
  #                                                "VisitVariabilityASDMG")
  # 
  # # Overall Subject Variability across visits in TD Group
  # VisitVariabilityTDMG <- extract_relevant_sdvisit_data(ChildLatency_m1_surrogate, 
  #                                                "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar > 0",
  #                                                "VisitVariabilityTDMG")
  # 
  # # Subject Variability across visits MG in Autism Group
  # VisitVariabilityASD_Q_F <- extract_relevant_sdvisit_data(ChildLatency_m1_surrogate, 
  #                                                "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > 0",
  #                                                "VisitVariabilityASD_Q_F")
  # 
  # # Overall Subject Variability across visits in TD Group
  # VisitVariabilityTDC_Q_F <- extract_relevant_sdvisit_data(ChildLatency_m1_surrogate, 
  #                                                "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar > 0",
  #                                                "VisitVariabilityTDC_Q_F")
  # 
  # # Subject Variability across visits MG in Autism Group
  # VisitVariabilityASD_Q_UF <- extract_relevant_sdvisit_data(ChildLatency_m1_surrogate, 
  #                                                "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > 0",
  #                                                "VisitVariabilityASD_Q_UF")
  # 
  # # Overall Subject Variability across visits in TD Group
  # VisitVariabilityTDC_Q_UF <- extract_relevant_sdvisit_data(ChildLatency_m1_surrogate, 
  #                                                "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > 0",
  #                                                "VisitVariabilityTDC_Q_UF")
  # 
  # # Overall Subject Variability across visits, differences
  # VisitVariabilityASDMGSlowerThanTDMG <- extract_relevant_sdvisit_data(ChildLatency_m1_surrogate, 
  #                                                "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar < DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar",
  #                                                "VisitVariabilityASDMG>TDMG")
  # 
  # # Overall Subject Variability across visits, differences
  # VisitVariabilityASD_Q_F_SlowerThanTD_Q_F <- extract_relevant_sdvisit_data(ChildLatency_m1_surrogate, 
  #                                                "DiagnosisASD:TaskQuestions:FamiliarityFamiliar > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
  #                                                "VisitVariabilityASD_Q_F_SlowerThanTD_Q_F")
  # 
  # # Overall Subject Variability across visits, differences
  # VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF <- extract_relevant_sdvisit_data(ChildLatency_m1_surrogate, 
  #                                                "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
  #                                                "VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF")
  # 
  # #DifferenceInDifference:
  # DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity <- extract_relevant_sdvisit_data(ChildLatency_m1_surrogate, 
  #                                                "(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)
  #                                                < (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
  #                                                "DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity")
  # 
  # DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG <- extract_relevant_sdvisit_data(ChildLatency_m1_surrogate, 
  #                                                "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar)
  #                                                > (DiagnosisASD:TaskQuestions:FamiliarityFamiliar - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
  #                                                "DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG")
  
  
  # Overall Sigma Variability across visits in Autism Group
  SigmaASDOverall <- extract_relevant_sigma_data(ChildLatency_m1_surrogate, 
                                                 "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 < 0",
                                                 "SigmaASDOverall")
  
  # Overall Sigma Variability across visits in TD Group
  SigmaTDCOverall <- extract_relevant_sigma_data(ChildLatency_m1_surrogate, 
                                                 "(sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) /3 < 0",
                                                 "SigmaTDCOverall")
  
  
  SigmaASDOverallHigherTDSOverall <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
                                                                      "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)/3 >
                                               (sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar + sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar  + sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)/3",
                                                                      "SigmaASDOverallHigherTDSOverall")
  
  # Subject Variability across visits MG in Autism Group
  SigmaASDMG <- extract_relevant_sigma_data(ChildLatency_m1_surrogate, 
                                            "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar < 0",
                                            "SigmaASDMG")
  
  # Overall Subject Variability across visits in TD Group
  SigmaTDMG <- extract_relevant_sigma_data(ChildLatency_m1_surrogate, 
                                           "sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar < 0",
                                           "SigmaTDMG")
  
  # Subject Variability across visits MG in Autism Group
  SigmaASD_Q_F <- extract_relevant_sigma_data(ChildLatency_m1_surrogate, 
                                              "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < 0",
                                              "SigmaASD_Q_F")
  
  # Overall Subject Variability across visits in TD Group
  SigmaTDC_Q_F <- extract_relevant_sigma_data(ChildLatency_m1_surrogate, 
                                              "sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar < 0",
                                              "SigmaTDC_Q_F")
  
  
  # Subject Variability across visits MG in Autism Group
  SigmaASD_Q_UF <- extract_relevant_sigma_data(ChildLatency_m1_surrogate, 
                                               "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar < 0 ",
                                               "SigmaASD_Q_UF")
  
  # Overall Subject Variability across visits in TD Group
  SigmaTDC_Q_UF <- extract_relevant_sigma_data(ChildLatency_m1_surrogate, 
                                               "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar < 0",
                                               "SigmaTDC_Q_UF")
  
  # Overall Subject Variability across visits, differences
  SigmaASDMGSlowerThanTDMG <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
                                                               "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar",
                                                               "SigmaASDMGSlowerThanTDMG")
  
  # Overall Subject Variability across visits, differences
  SigmaASD_Q_F_SlowerThanTD_Q_F <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
                                                                    "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar < sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar",
                                                                    "SigmaASD_Q_F_SlowerThanTD_Q_F")
  
  # Overall Subject Variability across visits, differences
  SigmaASD_Q_UF_SlowerThanTD_Q_UF <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
                                                                      "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar",
                                                                      "SigmaASD_Q_UF_SlowerThanTD_Q_UF")
  
  # Sigma across visits, differences
  # SigmaASD_MG_VisitLowerZero <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit < 0",
  #                                                "SigmaASD_MG_VisitLowerZero")
  
  # SigmaTDC_MG_VisitLowerZero <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit > 0",
  #                                                "SigmaTDC_MG_VisitLowerZero")
  # 
  # SigmaASD_MG_VisitLowerTDC <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit < sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
  #                                                "SigmaASD_MG_VisitLowerTDC")
  # 
  # SigmaASD_Overall_VisitLowerZero <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit + 
  #                                                sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit +
  #                                                sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit) / 3 < 0",
  #                                                "SigmaASD_Overall_VisitLowerZero")
  # 
  # SigmaTDC_Overall_VisitLowerZero <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "(sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit + 
  #                                                sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit +
  #                                                sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit) / 3 < 0",
  #                                                "SigmaTDC_Overall_VisitLowerZero")
  # 
  # # Sigma across visits, differences
  # SigmaASD_Q_F_VisitLowerZero <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit > 0",
  #                                                "SigmaASD_Q_F_VisitLowerZero")
  # 
  # SigmaTDC_Q_F_VisitLowerZero <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit < 0",
  #                                                "SigmaTDC_Q_F_VisitLowerZero")
  # 
  # SigmaASD_Q_F_VisitLowerTDC <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit > sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
  #                                                "SigmaASD_Q_F_VisitLowerTDC")
  # 
  # # Sigma across visits, differences
  # SigmaASD_Q_UF_VisitLowerZero <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < 0",
  #                                                "SigmaASD_Q_UF_VisitLowerZero")
  # 
  # SigmaTDC_Q_UF_VisitLowerZero <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > 0",
  #                                                "SigmaTDC_Q_UF_VisitLowerZero")
  # 
  # SigmaASD_Q_UF_VisitLowerTDC <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
  #                                                "SigmaASD_Q_UF_VisitLowerTDC")
  # 
  # SigmaVisitOverall_ASDLowerTDC <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit + 
  #                                                sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit +
  #                                                sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit) / 3 > (sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit + 
  #                                                sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit +
  #                                                sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit) /3",
  #                                                "SigmaVisitOverall_ASDLowerTDC")
  
  #DifferenceInDifference:
  DifferenceInDifferenceSigmaASDFamiliarityBiggerTDUnfamiliarity <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
                                                                                                     "(sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar - sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)
                                               > (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                                     "DifferenceInDifferenceSigmaASDFamiliarityBiggerTDUnfamiliarity")
  
  DifferenceInDifferenceSigmaASDMGSmallerThanTDMG <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
                                                                                      "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar - sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar)
                                               > (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                                      "DifferenceInDifferenceSigmaASDMGSmallerThanTDMG")
  
  #WithinGroups
  WithinGroupSigmaASD_Q_F_Faster_UF <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
                                                                        "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar > sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaASD_Q_FFaster_UF")
  
  WithinGroupSigmaTDC_Q_F_Faster_UF <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
                                                                        "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar > sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaTDC_Q_F_Faster_UF")
  
  WithinGroupSigmaASD_Q_F_Faster_MG <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
                                                                        "sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar > sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaASD_Q_F_Faster_MG")
  
  WithinGroupSigmaTDC_Q_F_Faster_MG <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
                                                                        "sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar < sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar ",
                                                                        "WithinGroupSigmaTDC_Q_F_Faster_MG")
  
  
  # WithinGroupSigmaASD_Q_F_Lower_UF_Visit <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit",
  #                                                "WithinGroupSigmaASD_Q_F_Lower_UF_Visit")
  # 
  # WithinGroupSigmaTDC_Q_F_Lower_UF_Visit <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit",
  #                                                "WithinGroupSigmaTDC_Q_F_Lower_UF_Visit")
  # 
  # WithinGroupSigmaASD_Q_F_Lower_MG_Visit <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit > sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit",
  #                                                "WithinGroupSigmaASD_Q_F_Lower_MG_Visit")
  # 
  # WithinGroupSigmaTDC_Q_F_Lower_MG_Visit <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
  #                                                "WithinGroupSigmaTDC_Q_F_Lower_MG_Visit")
  # 
  # 
  # SigmaTDC_Q_UF_VisitLowerZero <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit > 0",
  #                                                "SigmaTDC_Q_UF_VisitLowerZero")
  # 
  # SigmaASD_Q_UF_VisitLowerTDC <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit < sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit",
  #                                                "SigmaASD_Q_UF_VisitLowerTDC")
  # 
  # 
  # DifferenceInDifferenceSigmaVisit_Q_F_Lower_Q_UF <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "(sigma_DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:Visit - sigma_DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:Visit) <
  #                                                (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit)", "DifferenceInDifferenceSigmaVisit_Q_F_Lower_Q_UF")
  # 
  # 
  # DifferenceInDifferenceSigmaVisit_MG_Lower_Q_F <- extract_relevant_nonlatency_data(ChildLatency_m1_surrogate, 
  #                                                "(sigma_DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:Visit - sigma_DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:Visit) <
  #                                                (sigma_DiagnosisASD:TaskQuestions:FamiliarityFamiliar:Visit - sigma_DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:Visit)", "DifferenceInDifferenceSigmaVisit_MG_Lower_Q_F")
  
  
  
  # Correlation Estimates across visits in Autism Group
  post <- as_draws_df(ChildLatency_m1_surrogate)
  
  
  # Apply the function for each case
  cor_post_mu_MG_ASD <- Visitextract_and_rename(post, "DiagnosisASD", "TaskMatchingGame", "FamiliarityFamiliar") %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_MG_TDC <- Visitextract_and_rename(post, "DiagnosisTDC", "TaskMatchingGame", "FamiliarityFamiliar") %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_Q_ASD_Familiar <- Visitextract_and_rename(post, "DiagnosisASD", "TaskQuestions", "FamiliarityFamiliar") %>%
    select(Visit2, Visit4, Visit6) %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_Q_TDC_Familiar <- Visitextract_and_rename(post, "DiagnosisTDC", "TaskQuestions", "FamiliarityFamiliar") %>%
    select(Visit2, Visit4, Visit6) %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_Q_ASD_Unfamiliar <- Visitextract_and_rename(post, "DiagnosisASD", "TaskQuestions", "FamiliarityUnfamiliar") %>%
    select(Visit1, Visit3, Visit5, Visit7) %>%
    cor() %>%
    as.matrix()
  
  cor_post_mu_Q_TDC_Unfamiliar <- Visitextract_and_rename(post, "DiagnosisTDC", "TaskQuestions", "FamiliarityUnfamiliar") %>%
    select(Visit1, Visit3, Visit5, Visit7) %>%
    cor() %>%
    as.matrix()
  
  diag(cor_post_mu_MG_ASD) <- NA
  diag(cor_post_mu_MG_TDC) <- NA
  diag(cor_post_mu_Q_ASD_Familiar) <- NA
  diag(cor_post_mu_Q_TDC_Familiar) <- NA
  diag(cor_post_mu_Q_ASD_Unfamiliar) <- NA
  diag(cor_post_mu_Q_TDC_Unfamiliar) <- NA
  
  CorrelationSummaryStatistics <- tibble(
    Estimate = c(mean(cor_post_mu_MG_ASD, na.rm = TRUE), mean(cor_post_mu_Q_ASD_Familiar, na.rm = TRUE), mean(cor_post_mu_Q_ASD_Unfamiliar, na.rm = TRUE), 
                 mean(cor_post_mu_MG_TDC, na.rm = TRUE), mean(cor_post_mu_Q_TDC_Familiar, na.rm = TRUE), mean(cor_post_mu_Q_TDC_Unfamiliar, na.rm = TRUE)),
    CI.Lower = c(quantile(cor_post_mu_MG_ASD, probs = 0.025, na.rm = TRUE), quantile(cor_post_mu_Q_ASD_Familiar, probs = 0.025, na.rm = TRUE), quantile(cor_post_mu_Q_ASD_Unfamiliar, probs = 0.025, na.rm = TRUE), 
                 quantile(cor_post_mu_MG_TDC, probs = 0.025, na.rm = TRUE), quantile(cor_post_mu_Q_TDC_Familiar, probs = 0.025, na.rm = TRUE), quantile(cor_post_mu_Q_TDC_Unfamiliar, probs = 0.025, na.rm = TRUE)),
    CI.Upper = c(quantile(cor_post_mu_MG_ASD, probs = 0.975, na.rm = TRUE), quantile(cor_post_mu_Q_ASD_Familiar, probs = 0.975, na.rm = TRUE), quantile(cor_post_mu_Q_ASD_Unfamiliar, probs = 0.975, na.rm = TRUE), 
                 quantile(cor_post_mu_MG_TDC, probs = 0.975, na.rm = TRUE), quantile(cor_post_mu_Q_TDC_Familiar, probs = 0.975, na.rm = TRUE), quantile(cor_post_mu_Q_TDC_Unfamiliar, probs = 0.975, na.rm = TRUE)),
    Evid.Ratio = c(NA, NA, NA, NA, NA, NA),
    TestDescription = c(
      "CorrelationAcrossVisitsASDMG", "CorrelationAcrossVisitsASDFamiliarQuestions", "CorrelationAcrossVisitsASDUnfamiliarQuestions",
      "CorrelationAcrossVisitsTDMG", "CorrelationAcrossVisitsTDFamiliarQuestions", "CorrelationAcrossVisitsTDUnfamiliarQuestions"),
  )
  
  HypothesisResultsChildLatenciesModel1 <- rbind(OverallLatencyThreeConditionsASD,
                                                 OverallLatencyThreeConditionsTD,
                                                 MGOverallLatencyASD,
                                                 Q_F_OverallLatencyASD,
                                                 Q_UF_OverallLatencyASD,
                                                 MGOverallLatencyTDC,
                                                 Q_F_OverallLatencyTDC,
                                                 Q_UF_OverallLatencyTDC,
                                                 # BetaParameterASD_MG,
                                                 # BetaParameterASD_Q_F,
                                                 # BetaParameterASD_Q_UF,
                                                 # BetaParameterASD_Overall,
                                                 # BetaParameterTDC_MG,
                                                 # BetaParameterTDC_Q_F,
                                                 # BetaParameterTDC_Q_UF,
                                                 # BetaParameterTDC_Overall,
                                                 # OverallBetaParameterASDLowerTD,
                                                 # BetaParameterASDFasterTD_MG,
                                                 # BetaParameterASDFasterTD_Q_F,
                                                 # BetaParameterASDFasterTD_Q_UF,
                                                 # WithinGroupBetaParameterASD_MG_Slower_Q_F,
                                                 # WithinGroupBetaParameterTDC_MG_Slower_Q_F,
                                                 # DifferenceInDifferenceBetaParameterBetweenGroupsMGBiggerQ_F,
                                                 # WithinGroupBetaParameterASD_Q_F_Slower_Q_UF,
                                                 # WithinGroupBetaParameterTDC_Q_F_Slower_Q_UF,
                                                 # DifferenceInDifferenceBetaParameterBetweenGroupsFamiliarityBiggerUnfamiliarity,
                                                 OverallLatencyThreeConditionsASDFasterTD,
                                                 WithinASDGroupMGSlowerQuestionsFamiliar,
                                                 WithinTDGroupMGSlowerQuestionsFamiliar,
                                                 WithinAutismGroupQuestionsUnfamiliarSlowerFamiliar,
                                                 WithinTDGroupQuestionsUnfamiliarSlowerFamiliar,
                                                 WithinConditionsMatchingGameASDFasterTD,
                                                 WithinConditionsQuestionsFamiliarASDFasterTD,
                                                 WithinConditionsQuestionsUnfamiliarASDFasterTD,
                                                 DifferenceInDifferenceBetweenGroupsUnfamiliarBiggerFamiliarity,
                                                 DifferenceInDifferenceBetweenGroupsMGBiggerThanConversations,
                                                 SubjectVariabilityASD,
                                                 SubjectVariabilityTD,
                                                 SubjectVariabilityASDSmallerThanTD,
                                                 SubjectVariabilityASD_MG,
                                                 SubjectVariabilityASD_Q_F,
                                                 SubjectVariabilityASD_Q_UF,
                                                 SubjectVariabilityTD_MG,
                                                 SubjectVariabilityTD_Q_F,
                                                 SubjectVariabilityTD_Q_UF,
                                                 SubjectVariabilityASDSmallerThanTD_MG,
                                                 SubjectVariabilityASDSmallerThanTD_Q_F,
                                                 SubjectVariabilityASDSmallerThanTD_Q_UF,
                                                 SubjectVariabilityASD_MG_SmallerThan_Q_F,
                                                 SubjectVariabilityASD_Q_F_SmallerThan_Q_UF,
                                                 SubjectVariabilityTDC_MG_SmallerThan_Q_F,
                                                 SubjectVariabilityTDC_Q_F_SmallerThan_Q_UF,
                                                 DifferenceInDifferenceSubjectVariabilityUnfamiliarSmallerFamiliarity,
                                                 DifferenceInDifferenceSubjectVariabilityMatchingGameBiggerQuestions,
                                                 VisitVariabilityASDMG,
                                                 VisitVariabilityTDMG,
                                                 VisitVariabilityASD_Q_F,
                                                 VisitVariabilityTDC_Q_F,
                                                 VisitVariabilityASD_Q_UF,
                                                 VisitVariabilityTDC_Q_UF,
                                                 VisitVariabilityASDOverall,
                                                 VisitVariabilityTDCOverall,
                                                 VisitVariabilityASDMGSlowerThanTDMG,
                                                 VisitVariabilityASD_Q_F_SlowerThanTD_Q_F,
                                                 VisitVariabilityASD_Q_UF_SlowerThanTD_Q_UF,
                                                 VisitVariabilityASDOverallHigherTDSOverall,
                                                 DifferenceInDifferenceVisitVariabilityASDMGSmallerThanTDMG,
                                                 DifferenceInDifferenceVisitVariabilityASDFamiliarityTDUnfamiliarity,
                                                 SigmaASDOverall,
                                                 SigmaTDCOverall,
                                                 SigmaASDOverallHigherTDSOverall,
                                                 SigmaASDMG,
                                                 SigmaTDMG,
                                                 SigmaASD_Q_F,
                                                 SigmaTDC_Q_F,
                                                 SigmaASD_Q_UF,
                                                 SigmaTDC_Q_UF,
                                                 SigmaASDMGSlowerThanTDMG,
                                                 SigmaASD_Q_F_SlowerThanTD_Q_F,
                                                 SigmaASD_Q_UF_SlowerThanTD_Q_UF,
                                                 # SigmaASD_MG_VisitLowerZero,
                                                 # SigmaTDC_MG_VisitLowerZero,
                                                 # SigmaASD_MG_VisitLowerTDC,
                                                 # SigmaASD_Q_F_VisitLowerZero,
                                                 # SigmaTDC_Q_F_VisitLowerZero,
                                                 # SigmaASD_Q_F_VisitLowerTDC,
                                                 # SigmaASD_Q_UF_VisitLowerZero,
                                                 # SigmaTDC_Q_UF_VisitLowerZero,
                                                 # SigmaASD_Q_UF_VisitLowerTDC,
                                                 # SigmaASD_Overall_VisitLowerZero,
                                                 # SigmaTDC_Overall_VisitLowerZero,
                                                 # WithinGroupSigmaASD_Q_F_Lower_UF_Visit,
                                                 # WithinGroupSigmaTDC_Q_F_Lower_UF_Visit,
                                                 # WithinGroupSigmaASD_Q_F_Lower_MG_Visit,
                                                 # WithinGroupSigmaTDC_Q_F_Lower_MG_Visit,
                                                 WithinGroupSigmaASD_Q_F_Faster_UF,
                                                 WithinGroupSigmaTDC_Q_F_Faster_UF,
                                                 WithinGroupSigmaASD_Q_F_Faster_MG,
                                                 WithinGroupSigmaTDC_Q_F_Faster_MG,
                                                 # SigmaVisitOverall_ASDLowerTDC,
                                                 # DifferenceInDifferenceSigmaVisit_Q_F_Lower_Q_UF,
                                                 # DifferenceInDifferenceSigmaVisit_MG_Lower_Q_F,
                                                 DifferenceInDifferenceSigmaASDFamiliarityBiggerTDUnfamiliarity,
                                                 DifferenceInDifferenceSigmaASDMGSmallerThanTDMG,
                                                 CorrelationSummaryStatistics) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  return(HypothesisResultsChildLatenciesModel1)
}

extract_all_hypothesis_tests_m3_shortutterances <- function(model) {
  ASD_OtherPredictability_MG_F <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords,
                                                        "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability < 0",
                                                        "ASD_OtherPredictability_MG_F")
  
  TDC_OtherPredictability_MG_F <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, 
                                                        "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability < 0",
                                                        "TDC_OtherPredictability_MG_F")
  
  ASD_OtherPredictability_Q_F <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, 
                                                       "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability < 0",
                                                       "ASD_OtherPredictability_Q_F")
  TDC_OtherPredictability_Q_F <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, 
                                                       "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability < 0",
                                                       "TDC_OtherPredictability_Q_F")
  
  ASD_OtherPredictability_Q_UF <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, 
                                                        "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability < 0",
                                                        "ASD_OtherPredictability_Q_UF")
  TDC_OtherPredictability_Q_UF <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords,
                                                        "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability < 0",
                                                        "TDC_OtherPredictability_Q_UF")
  
  ASD_OtherPredictability_overall <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, 
                                                           "((DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability) / 3) < 0", "ASD_OtherPredictability_overall")
  
  TDC_OtherPredictability_overall <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords,
                                                           "((DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability) / 3) < 0",
                                                           "TDC_OtherPredictability_overall")
  
  WithinGroupASD_OtherPredictability_MG_Higher_Q_F <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability > DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability", 
                                                                            "WithinGroupASD_OtherPredictability_MG_Higher_Q_F")
  
  WithinGroupASD_OtherPredictability_Q_F_Higher_Q_UF <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability < DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability", 
                                                                              "WithinGroupASD_OtherPredictability_Q_F_Higher_Q_UF")
  
  WithinGroupTDC_OtherPredictability_MG_Higher_Q_F <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability", 
                                                                            "WithinGroupTDC_OtherPredictability_MG_Higher_Q_F")
  
  WithinGroupTDC_OtherPredictability_Q_F_Higher_Q_UF <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability", 
                                                                              "WithinGroupTDC_OtherPredictability_Q_F_Higher_Q_UF")
  
  OtherPredictabilityHypothesesFull <- rbind(ASD_OtherPredictability_MG_F, ASD_OtherPredictability_Q_F, ASD_OtherPredictability_Q_UF,
                                             TDC_OtherPredictability_MG_F, TDC_OtherPredictability_Q_F, TDC_OtherPredictability_Q_UF,
                                             ASD_OtherPredictability_overall, TDC_OtherPredictability_overall, WithinGroupASD_OtherPredictability_MG_Higher_Q_F, WithinGroupASD_OtherPredictability_Q_F_Higher_Q_UF, WithinGroupTDC_OtherPredictability_MG_Higher_Q_F, WithinGroupTDC_OtherPredictability_Q_F_Higher_Q_UF)
  
  
  OtherPredictability_MG_F_ASDHigherTDC <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords,
                                                                 "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability < 
           DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability",
                                                                 "OtherPredictability_MG_F_ASDHigherTDC")
  
  OtherPredictability_Q_F_ASDHigherTDC <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords,
                                                                "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability > 
           DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability",
                                                                "OtherPredictability_Q_F_ASDHigherTDC")
  
  OtherPredictability_Q_UF_ASDHigherTDC <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords,
                                                                 "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability > 
           DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability",
                                                                 "OtherPredictability_Q_UF_ASDHigherTDC")
  
  OtherPredictability_overall_ASDHigherTDC <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, 
                                                                    "((DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability +
                                                       DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability +
                                                       DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability) / 3) >
                                                       ((DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability +
                                                       DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability +
                                                       DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability) / 3)",
                                                                    "OtherPredictability_overall_ASDHigherTDC")
  
  DifferenceOtherPredictability_ASDLowerMatchingGameTDCVersusQuestions <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords,
                                                                                                "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OtherPredictability) <
            (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability)",
                                                                                                "DifferenceOtherPredictability_ASDLowerMatchingGameTDCVersusQuestions")
  
  DifferenceOtherPredictability_ASDLowerFamiliarityTDCVersusUnfamiliarity <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords,
                                                                                                   "(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OtherPredictability - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OtherPredictability) >
            (DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OtherPredictability)",
                                                                                                   "DifferenceOtherPredictability_ASDHigherFamiliarityTDCVersusUnfamiliarity")
  
  OtherPredictabilityComparisonsHypothesesFull <- rbind(
    OtherPredictabilityHypothesesFull,
    OtherPredictability_MG_F_ASDHigherTDC,
    OtherPredictability_Q_F_ASDHigherTDC,
    OtherPredictability_Q_UF_ASDHigherTDC,
    OtherPredictability_overall_ASDHigherTDC,
    DifferenceOtherPredictability_ASDLowerMatchingGameTDCVersusQuestions,
    DifferenceOtherPredictability_ASDLowerFamiliarityTDCVersusUnfamiliarity) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  
  #Let's add the same hypotheses for OwnPredictability:
  ASD_OwnPredictability_MG_F <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords,
                                                      "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability < 0",
                                                      "ASD_OwnPredictability_MG_F")
  
  TDC_OwnPredictability_MG_F <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, 
                                                      "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability > 0",
                                                      "TDC_OwnPredictability_MG_F")
  
  ASD_OwnPredictability_Q_F <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, 
                                                     "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability > 0",
                                                     "ASD_OwnPredictability_Q_F")
  TDC_OwnPredictability_Q_F <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, 
                                                     "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability > 0",
                                                     "TDC_OwnPredictability_Q_F")
  
  ASD_OwnPredictability_Q_UF <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, 
                                                      "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability > 0",
                                                      "ASD_OwnPredictability_Q_UF")
  TDC_OwnPredictability_Q_UF <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords,
                                                      "DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability > 0",
                                                      "TDC_OwnPredictability_Q_UF")
  
  
  ASD_OwnPredictability_overall <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, 
                                                         "((DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability + DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability + DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability) / 3) > 0", "ASD_OwnPredictability_overall")
  
  TDC_OwnPredictability_overall <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords,
                                                         "((DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability + DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability + DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability) / 3) > 0",
                                                         "TDC_OwnPredictability_overall")
  
  WithinGroupASD_OwnPredictability_MG_Higher_Q_F <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability < DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability", 
                                                                          "WithinGroupASD_OwnPredictability_MG_Higher_Q_F")
  
  WithinGroupASD_OwnPredictability_Q_F_Higher_Q_UF <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability < DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability", 
                                                                            "WithinGroupASD_OwnPredictability_Q_F_Higher_Q_UF")
  
  WithinGroupTDC_OwnPredictability_MG_Higher_Q_F <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, "DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability > DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability", 
                                                                          "WithinGroupTDC_OwnPredictability_MG_Higher_Q_F")
  
  WithinGroupTDC_OwnPredictability_Q_F_Higher_Q_UF <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, "DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability < DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability", 
                                                                            "WithinGroupTDC_OwnPredictability_Q_F_Higher_Q_UF")
  
  
  OwnPredictabilityHypothesesFull <- rbind(ASD_OwnPredictability_MG_F, ASD_OwnPredictability_Q_F, ASD_OwnPredictability_Q_UF,
                                           TDC_OwnPredictability_MG_F, TDC_OwnPredictability_Q_F, TDC_OwnPredictability_Q_UF,
                                           ASD_OwnPredictability_overall, TDC_OwnPredictability_overall, WithinGroupASD_OwnPredictability_MG_Higher_Q_F, WithinGroupASD_OwnPredictability_Q_F_Higher_Q_UF, WithinGroupTDC_OwnPredictability_MG_Higher_Q_F, WithinGroupTDC_OwnPredictability_Q_F_Higher_Q_UF)
  
  
  OwnPredictability_MG_F_ASDHigherTDC <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords,
                                                               "DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability < 
           DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability",
                                                               "OwnPredictability_MG_F_ASDHigherTDC")
  
  OwnPredictability_Q_F_ASDHigherTDC <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords,
                                                              "DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability > 
           DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability",
                                                              "OwnPredictability_Q_F_ASDHigherTDC")
  
  OwnPredictability_Q_UF_ASDHigherTDC <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords,
                                                               "DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability > 
           DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability",
                                                               "OwnPredictability_Q_UF_ASDHigherTDC")
  
  OwnPredictability_overall_ASDHigherTDC <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords, 
                                                                  "((DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability +
                                                       DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability +
                                                       DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability) / 3) >
                                                       ((DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability +
                                                       DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability +
                                                       DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability) / 3)",
                                                                  "OwnPredictability_overall_ASDHigherTDC")
  
  DifferenceOwnPredictability_ASDLowerMatchingGameTDCVersusQuestions <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords,
                                                                                              "(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability - DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:OwnPredictability) <
            (DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability)",
                                                                                              "DifferenceOwnPredictability_ASDLowerMatchingGameTDCVersusQuestions")
  
  DifferenceOwnPredictability_ASDLowerFamiliarityTDCVersusUnfamiliarity <- extract_relevant_data(ChildLatency_pred_m3_morethanthreewords,
                                                                                                 "(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:OwnPredictability - DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:OwnPredictability) <
            (DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability - DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:OwnPredictability)",
                                                                                                 "DifferenceOwnPredictability_ASDHigherFamiliarityTDCVersusUnfamiliarity")
  
  OwnPredictabilityComparisonsHypothesesFull <- rbind(
    OwnPredictabilityHypothesesFull,
    OwnPredictability_MG_F_ASDHigherTDC,
    OwnPredictability_Q_F_ASDHigherTDC,
    OwnPredictability_Q_UF_ASDHigherTDC,
    OwnPredictability_overall_ASDHigherTDC,
    DifferenceOwnPredictability_ASDLowerMatchingGameTDCVersusQuestions,
    DifferenceOwnPredictability_ASDLowerFamiliarityTDCVersusUnfamiliarity) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  
  PredictabilityComparisonsHypothesesFull <- rbind(OtherPredictabilityComparisonsHypothesesFull, OwnPredictabilityComparisonsHypothesesFull)
  return(PredictabilityComparisonsHypothesesFull)
}

run_single_crqa <- function(id_val, task_val) {
  cat("\nProcessing ID:", id_val, "Task:", task_val, "\n")
  
  # Filter data
  dyad_data <- d %>% 
    filter(ID == id_val, Task == task_val)
  
  # Extract latencies
  adult_latencies <- dyad_data %>% 
    filter(Speaker == "Adult") %>%
    arrange(StartTime) %>%
    pull(Latency) %>%
    na.omit()
  
  child_latencies <- dyad_data %>% 
    filter(Speaker == "Child") %>%
    arrange(StartTime) %>%
    pull(Latency) %>%
    na.omit()
  
  # Check if we have enough data
  if(length(adult_latencies) < 5 || length(child_latencies) < 5) {
    cat("Not enough data points\n")
    return(NULL)
  }
  
  # Ensure equal length
  min_length <- min(length(adult_latencies), length(child_latencies))
  adult_latencies <- adult_latencies[1:min_length]
  child_latencies <- child_latencies[1:min_length]
  
  cat("Number of data points:", min_length, "\n")
  
  # Check for data issues
  if(any(is.infinite(adult_latencies)) || any(is.infinite(child_latencies))) {
    cat("Found infinite values\n")
    return(NULL)
  }
  
  # Try CRQA
  tryCatch({
    results <- crqa(adult_latencies, child_latencies,
                    delay = 1, embed = 1,
                    radius = 0.1,
                    normalize = TRUE,
                    rescale = 0,
                    mindiagline = 2, minvertline = 2)
    
    # Return results as a data frame row
    data.frame(
      ID = id_val,
      Task = task_val,
      DataPoints = min_length,
      RR = results$RR,
      DET = results$DET,
      LAM = results$LAM,
      L = results$L,
      TT = results$TT,
      ENTR = results$ENTR
    )
  }, error = function(e) {
    cat("CRQA error:", e$message, "\n")
    return(NULL)
  })
}

extract_all_hypothesis_tests_m1_overlaps <- function(model) {
  
  ASD_OverlapsHigherThanZero <- extract_relevant_nonlatency_data(OverlappingChild_m, 
                                                                 "(inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) +
           inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) +
           inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar))/3 >
           0",
                                                                 "ASD_OverlapsHigherThanZero")
  
  ASD_MG_F <- extract_relevant_nonlatency_data(OverlappingChild_m, 
                                               "inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) > 0",
                                               "ASD_MG_F")
  
  ASD_Q_F <- extract_relevant_nonlatency_data(OverlappingChild_m, 
                                              "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) > 0",
                                              "ASD_Q_F")
  
  ASD_Q_UF <- extract_relevant_nonlatency_data(OverlappingChild_m, 
                                               "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar) > 0",
                                               "ASD_Q_UF")
  
  TDC_OverlapsHigherThanZero <- extract_relevant_nonlatency_data(OverlappingChild_m, 
                                                                 "(inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar) +
           inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar) +
           inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar))/3 >
           0",
                                                                 "TDC_OverlapsHigherThanZero")
  
  TD_MG_F <- extract_relevant_nonlatency_data(OverlappingChild_m, 
                                              "inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar) > 0",
                                              "TD_MG_F")
  
  TD_Q_F <- extract_relevant_nonlatency_data(OverlappingChild_m, 
                                             "inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar) > 0",
                                             "TD_Q_F")
  
  TD_Q_UF <- extract_relevant_nonlatency_data(OverlappingChild_m, 
                                              "inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) > 0",
                                              "TD_Q_UF")
  
  
  
  WithinGroup_ASD_MG_F_Higher_Q_F <- extract_relevant_nonlatency_data(OverlappingChild_m, 
                                                                      "inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) < inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar)",
                                                                      "WithinGroup_ASD_MG_F_Higher_Q_F")
  
  WithinGroup_ASD_Q_F_Higher_Q_UF <- extract_relevant_nonlatency_data(OverlappingChild_m, 
                                                                      "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) < inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)",
                                                                      "WithinGroup_ASD_Q_F_Higher_Q_UF")
  
  WithinGroup_TDC_MG_F_Higher_Q_F <- extract_relevant_nonlatency_data(OverlappingChild_m, 
                                                                      "inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar) > inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                      "WithinGroup_TDC_MG_F_Higher_Q_F")
  
  WithinGroup_TDC_Q_F_Higher_Q_UF <- extract_relevant_nonlatency_data(OverlappingChild_m, 
                                                                      "inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar) < inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)",
                                                                      "WithinGroup_TDC_Q_F_Higher_Q_UF")
  
  
  ASDOverlapsHigherThanTDCArossConditions <- extract_relevant_nonlatency_data(OverlappingChild_m, 
                                                                              "(inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) +
            inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) +
            inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar))/3 >
            (inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar) +
            inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar) +
            inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar))/3",
                                                                              "ASDOverlapsHigherThanTDCArossConditions")
  
  ASDOverlapsHigherThanTDCWithinMatchingGame <- extract_relevant_nonlatency_data(OverlappingChild_m, 
                                                                                 "inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) > inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar)",
                                                                                 "ASDOverlapsHigherThanTDCWithinMatchingGame")
  
  ASDOverlapsHigherThanTDCWithinQFamiliar <- extract_relevant_nonlatency_data(OverlappingChild_m, 
                                                                              "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) > inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                              "ASDOverlapsHigherThanTDCWithinQFamiliar")
  
  ASDOverlapsHigherThanTDCWithinQUnfamiliar <- extract_relevant_nonlatency_data(OverlappingChild_m, 
                                                                                "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar) > inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)",
                                                                                "ASDOverlapsHigherThanTDCWithinQUnfamiliar")
  
  #Differences between groups according to Condition
  DifferencesBetweenGroupsQuestionsHigherMatchingGame <- extract_relevant_nonlatency_data(OverlappingChild_m, 
                                                                                          "(inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) - inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar)) <
           (inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) - inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar))",
                                                                                          "DifferenceInDifferenceBetweenGroupsASDHigherTDCMatchingGame")
  
  #Differences between groups according to Familiarty
  DifferencesBetweenGroupsUnfamiliarityHigherFamiliarity <- extract_relevant_nonlatency_data(OverlappingChild_m, 
                                                                                             "(inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar) - inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)) <
           (inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) - inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar))",
                                                                                             "DifferenceInDifferenceBetweenGroupsUnfamiliarityHigherFamiliarity")
  
  
  HypothesisResultsChildOverlaps <- rbind(ASD_OverlapsHigherThanZero,
                                          TDC_OverlapsHigherThanZero,
                                          ASD_MG_F,
                                          ASD_Q_F,
                                          ASD_Q_UF,
                                          TD_MG_F,
                                          TD_Q_F,
                                          TD_Q_UF,
                                          ASDOverlapsHigherThanTDCArossConditions,
                                          ASDOverlapsHigherThanTDCWithinMatchingGame,
                                          ASDOverlapsHigherThanTDCWithinQFamiliar,
                                          ASDOverlapsHigherThanTDCWithinQUnfamiliar,
                                          DifferencesBetweenGroupsQuestionsHigherMatchingGame, 
                                          DifferencesBetweenGroupsUnfamiliarityHigherFamiliarity,
                                          WithinGroup_ASD_MG_F_Higher_Q_F,
                                          WithinGroup_ASD_Q_F_Higher_Q_UF,
                                          WithinGroup_TDC_MG_F_Higher_Q_F,
                                          WithinGroup_TDC_Q_F_Higher_Q_UF) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  return(HypothesisResultsChildOverlaps)
}

extract_all_hypothesis_tests_m1_overlaps_adults <- function(model) {
  ASD_OverlapsHigherThanZero <- extract_relevant_nonlatency_data(OverlappingAdult_m, 
                                                                 "(inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) +
           inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) +
           inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar))/3 >
           0",
                                                                 "ASD_OverlapsHigherThanZero")
  
  ASD_MG_F <- extract_relevant_nonlatency_data(OverlappingAdult_m, 
                                               "inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) > 0",
                                               "ASD_MG_F")
  
  ASD_Q_F <- extract_relevant_nonlatency_data(OverlappingAdult_m, 
                                              "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) > 0",
                                              "ASD_Q_F")
  
  ASD_Q_UF <- extract_relevant_nonlatency_data(OverlappingAdult_m, 
                                               "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar) > 0",
                                               "ASD_Q_UF")
  
  TDC_OverlapsHigherThanZero <- extract_relevant_nonlatency_data(OverlappingAdult_m, 
                                                                 "(inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar) +
           inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar) +
           inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar))/3 >
           0",
                                                                 "TDC_OverlapsHigherThanZero")
  
  TD_MG_F <- extract_relevant_nonlatency_data(OverlappingAdult_m, 
                                              "inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar) > 0",
                                              "TD_MG_F")
  
  TD_Q_F <- extract_relevant_nonlatency_data(OverlappingAdult_m, 
                                             "inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar) > 0",
                                             "TD_Q_F")
  
  TD_Q_UF <- extract_relevant_nonlatency_data(OverlappingAdult_m, 
                                              "inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) > 0",
                                              "TD_Q_UF")
  
  
  
  WithinGroup_ASD_MG_F_Higher_Q_F <- extract_relevant_nonlatency_data(OverlappingAdult_m, 
                                                                      "inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) < inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar)",
                                                                      "WithinGroup_ASD_MG_F_Higher_Q_F")
  
  WithinGroup_ASD_Q_F_Higher_Q_UF <- extract_relevant_nonlatency_data(OverlappingAdult_m, 
                                                                      "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) < inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar)",
                                                                      "WithinGroup_ASD_Q_F_Higher_Q_UF")
  
  WithinGroup_TDC_MG_F_Higher_Q_F <- extract_relevant_nonlatency_data(OverlappingAdult_m, 
                                                                      "inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar) > inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                      "WithinGroup_TDC_MG_F_Higher_Q_F")
  
  WithinGroup_TDC_Q_F_Higher_Q_UF <- extract_relevant_nonlatency_data(OverlappingAdult_m, 
                                                                      "inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar) < inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)",
                                                                      "WithinGroup_TDC_Q_F_Higher_Q_UF")
  
  
  ASDOverlapsHigherThanTDCArossConditions <- extract_relevant_nonlatency_data(OverlappingAdult_m, 
                                                                              "(inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) +
            inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) +
            inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar))/3 >
            (inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar) +
            inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar) +
            inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar))/3",
                                                                              "ASDOverlapsHigherThanTDCArossConditions")
  
  ASDOverlapsHigherThanTDCWithinMatchingGame <- extract_relevant_nonlatency_data(OverlappingAdult_m, 
                                                                                 "inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) > inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar)",
                                                                                 "ASDOverlapsHigherThanTDCWithinMatchingGame")
  
  ASDOverlapsHigherThanTDCWithinQFamiliar <- extract_relevant_nonlatency_data(OverlappingAdult_m, 
                                                                              "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) > inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                              "ASDOverlapsHigherThanTDCWithinQFamiliar")
  
  ASDOverlapsHigherThanTDCWithinQUnfamiliar <- extract_relevant_nonlatency_data(OverlappingAdult_m, 
                                                                                "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar) > inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)",
                                                                                "ASDOverlapsHigherThanTDCWithinQUnfamiliar")
  
  #Differences between groups according to Condition
  DifferencesBetweenGroupsQuestionsHigherMatchingGame <- extract_relevant_nonlatency_data(OverlappingAdult_m, 
                                                                                          "(inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) - inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar)) <
           (inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) - inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar))",
                                                                                          "DifferenceInDifferenceBetweenGroupsASDHigherTDCMatchingGame")
  
  #Differences between groups according to Familiarty
  DifferencesBetweenGroupsUnfamiliarityHigherFamiliarity <- extract_relevant_nonlatency_data(OverlappingAdult_m, 
                                                                                             "(inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar) - inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)) <
           (inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) - inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar))",
                                                                                             "DifferenceInDifferenceBetweenGroupsUnfamiliarityHigherFamiliarity")
  
  
  HypothesisResultsAdultOverlaps <- rbind(ASD_OverlapsHigherThanZero,
                                          TDC_OverlapsHigherThanZero,
                                          ASD_MG_F,
                                          ASD_Q_F,
                                          ASD_Q_UF,
                                          TD_MG_F,
                                          TD_Q_F,
                                          TD_Q_UF,
                                          ASDOverlapsHigherThanTDCArossConditions,
                                          ASDOverlapsHigherThanTDCWithinMatchingGame,
                                          ASDOverlapsHigherThanTDCWithinQFamiliar,
                                          ASDOverlapsHigherThanTDCWithinQUnfamiliar,
                                          DifferencesBetweenGroupsQuestionsHigherMatchingGame, 
                                          DifferencesBetweenGroupsUnfamiliarityHigherFamiliarity,
                                          WithinGroup_ASD_MG_F_Higher_Q_F,
                                          WithinGroup_ASD_Q_F_Higher_Q_UF,
                                          WithinGroup_TDC_MG_F_Higher_Q_F,
                                          WithinGroup_TDC_Q_F_Higher_Q_UF) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  return(HypothesisResultsAdultOverlaps)
}

extract_all_hypothesis_tests_m2_overlaps_child <- function(model) {
  ASD_OverlapsHigherThanZero <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                 "(inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) +
           inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) +
           inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar))/3 >
           0",
                                                                 "ASD_OverlapsHigherThanZero")
  
  ASD_MG_F <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                               "inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) > 0",
                                               "ASD_MG_F")
  
  ASD_Q_F <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                              "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) > 0",
                                              "ASD_Q_F")
  
  ASD_Q_UF <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                               "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar) > 0",
                                               "ASD_Q_UF")
  
  TDC_OverlapsHigherThanZero <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                 "(inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar) +
           inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar) +
           inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar))/3 >
           0",
                                                                 "TDC_OverlapsHigherThanZero")
  
  TD_MG_F <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                              "inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar) > 0",
                                              "TD_MG_F")
  
  TD_Q_F <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                             "inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar) > 0",
                                             "TD_Q_F")
  
  TD_Q_UF <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                              "inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar) > 0",
                                              "TD_Q_UF")
  
  ASDOverlapsHigherThanTDCArossConditions <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                              "(inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) +
            inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) +
            inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar))/3 >
            (inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar) +
            inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar) +
            inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar))/3",
                                                                              "ASDOverlapsHigherThanTDCArossConditions")
  
  ASDOverlapsHigherThanTDCWithinMatchingGame <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                                 "inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) > inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar)",
                                                                                 "ASDOverlapsHigherThanTDCWithinMatchingGame")
  
  ASDOverlapsHigherThanTDCWithinQFamiliar <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                              "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) > inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar)",
                                                                              "ASDOverlapsHigherThanTDCWithinQFamiliar")
  
  ASDOverlapsHigherThanTDCWithinQUnfamiliar <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                                "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar) > inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)",
                                                                                "ASDOverlapsHigherThanTDCWithinQUnfamiliar")
  
  #Differences between groups according to Condition
  DifferencesBetweenGroupsQuestionsHigherMatchingGame <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                                          "(inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar) - inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar)) <
           (inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) - inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar))",
                                                                                          "DifferenceInDifferenceBetweenGroupsASDHigherTDCMatchingGame")
  
  #Differences between groups according to Familiarty
  DifferencesBetweenGroupsUnfamiliarityHigherFamiliarity <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                                             "(inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar) - inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar)) <
           (inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar) - inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar))",
                                                                                             "DifferenceInDifferenceBetweenGroupsUnfamiliarityHigherFamiliarity")
  
  #Cognitive Skills
  ASD_OverlapsHigherThanZero_Cognition <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                           "(inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled) +
           inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled) +
           inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled))/3 >
           0",
                                                                           "ASD_OverlapsHigherThanZero_Cognition")
  
  ASD_MG_F_Cognition <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                         "inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled) > 0",
                                                         "ASD_MG_F_Cognition")
  
  ASD_Q_F_Cognition <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                        "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled) > 0",
                                                        "ASD_Q_F_Cognition")
  
  ASD_Q_UF_Cognition <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                         "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled) > 0",
                                                         "ASD_Q_UF_Cognition")
  
  TDC_OverlapsHigherThanZero_Cognition <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                           "(inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled) +
           inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled) +
           inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled))/3 >
           0",
                                                                           "TDC_OverlapsHigherThanZero_Cognition")
  
  TD_MG_F_Cognition <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                        "inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled) > 0",
                                                        "TD_MG_F_Cognition")
  
  TD_Q_F_Cognition <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                       "inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled) > 0",
                                                       "TD_Q_F_Cognition")
  
  TD_Q_UF_Cognition <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                        "inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled) > 0",
                                                        "TD_Q_UF_Cognition")
  
  ASDOverlapsHigherThanTDCArossConditions_Cognition <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                                        "(inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled) +
            inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled) +
            inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled))/3 >
            (inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled) +
            inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled) +
            inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled))/3",
                                                                                        "ASDOverlapsHigherThanTDCArossConditions_Cognition")
  
  ASDOverlapsHigherThanTDCWithinMatchingGame_Cognition <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                                           "inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled) > inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled)",
                                                                                           "ASDOverlapsHigherThanTDCWithinMatchingGame_Cognition")
  
  ASDOverlapsHigherThanTDCWithinQFamiliar_Cognition <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                                        "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled) > inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled)",
                                                                                        "ASDOverlapsHigherThanTDCWithinQFamiliar_Cognition")
  
  ASDOverlapsHigherThanTDCWithinQUnfamiliar_Cognition <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                                          "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled) > inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled)",
                                                                                          "ASDOverlapsHigherThanTDCWithinQUnfamiliar_Cognition")
  
  #Differences between groups according to Condition
  DifferencesBetweenGroupsQuestionsHigherMatchingGame_Cognition <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                                                    "(inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled) - inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:CognitionS_scaled)) <
           (inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled) - inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled))",
                                                                                                    "DifferencesBetweenGroupsQuestionsHigherMatchingGame_Cognition")
  
  #Differences between groups according to Familiarty
  DifferencesBetweenGroupsUnfamiliarityHigherFamiliarity_Cognition <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                                                       "(inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled) - inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:CognitionS_scaled)) <
           (inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled) - inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:CognitionS_scaled))",
                                                                                                       "DifferencesBetweenGroupsUnfamiliarityHigherFamiliarity_Cognition")
  
  #Language Skills
  ASD_OverlapsHigherThanZero_Language <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                          "(inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled) +
           inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled) +
           inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled))/3 >
           0",
                                                                          "ASD_OverlapsHigherThanZero_Language")
  
  ASD_MG_F_Language <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                        "inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled) > 0",
                                                        "ASD_MG_F_Language")
  
  ASD_Q_F_Language <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                       "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled) > 0",
                                                       "ASD_Q_F_Language")
  
  ASD_Q_UF_Language <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                        "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled) > 0",
                                                        "ASD_Q_UF_Language")
  
  TDC_OverlapsHigherThanZero_Language <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                          "(inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled) +
           inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled) +
           inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled))/3 >
           0",
                                                                          "TDC_OverlapsHigherThanZero_Language")
  
  TD_MG_F_Language <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                       "inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled) > 0",
                                                       "TD_MG_F_Language")
  
  TD_Q_F_Language <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                      "inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled) > 0",
                                                      "TD_Q_F_Language")
  
  TD_Q_UF_Language <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                       "inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled) > 0",
                                                       "TD_Q_UF_Language")
  
  ASDOverlapsHigherThanTDCArossConditions_Language <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                                       "(inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled) +
            inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled) +
            inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled))/3 >
            (inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled) +
            inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled) +
            inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled))/3",
                                                                                       "ASDOverlapsHigherThanTDCArossConditions_Language")
  
  ASDOverlapsHigherThanTDCWithinMatchingGame_Language <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                                          "inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled) > inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled)",
                                                                                          "ASDOverlapsHigherThanTDCWithinMatchingGame_Language")
  
  ASDOverlapsHigherThanTDCWithinQFamiliar_Language <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                                       "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled) > inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled)",
                                                                                       "ASDOverlapsHigherThanTDCWithinQFamiliar_Language")
  
  ASDOverlapsHigherThanTDCWithinQUnfamiliar_Language <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                                         "inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled) > inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)",
                                                                                         "ASDOverlapsHigherThanTDCWithinQUnfamiliar_Language")
  
  #Differences between groups according to Condition
  DifferencesBetweenGroupsQuestionsHigherMatchingGame_Language <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                                                   "(inv_logit_scaled(DiagnosisASD:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled) - inv_logit_scaled(DiagnosisTDC:TaskMatchingGame:FamiliarityFamiliar:LanguageS_scaled)) <
           (inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled) - inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled))",
                                                                                                   "DifferencesBetweenGroupsQuestionsHigherMatchingGame_Language")
  
  #Differences between groups according to Familiarty
  DifferencesBetweenGroupsUnfamiliarityHigherFamiliarity_Language <- extract_relevant_nonlatency_data(OverlappingChild_m2_IndividualDifferences, 
                                                                                                      "(inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled) - inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityUnfamiliar:LanguageS_scaled)) <
           (inv_logit_scaled(DiagnosisASD:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled) - inv_logit_scaled(DiagnosisTDC:TaskQuestions:FamiliarityFamiliar:LanguageS_scaled))",
                                                                                                      "DifferencesBetweenGroupsUnfamiliarityHigherFamiliarity_Language")
  
  HypothesisResultsChildOverlaps <- rbind(ASD_OverlapsHigherThanZero,
                                          TDC_OverlapsHigherThanZero,
                                          ASD_MG_F,
                                          ASD_Q_F,
                                          ASD_Q_UF,
                                          TD_MG_F,
                                          TD_Q_F,
                                          TD_Q_UF,
                                          ASDOverlapsHigherThanTDCArossConditions,
                                          ASDOverlapsHigherThanTDCWithinMatchingGame,
                                          ASDOverlapsHigherThanTDCWithinQFamiliar,
                                          ASDOverlapsHigherThanTDCWithinQUnfamiliar,
                                          DifferencesBetweenGroupsQuestionsHigherMatchingGame, 
                                          DifferencesBetweenGroupsUnfamiliarityHigherFamiliarity,
                                          ASD_OverlapsHigherThanZero_Cognition,
                                          TDC_OverlapsHigherThanZero_Cognition,
                                          ASD_MG_F_Cognition,
                                          ASD_Q_F_Cognition,
                                          ASD_Q_UF_Cognition,
                                          TD_MG_F_Cognition,
                                          TD_Q_F_Cognition,
                                          TD_Q_UF_Cognition,
                                          ASDOverlapsHigherThanTDCArossConditions_Cognition,
                                          ASDOverlapsHigherThanTDCWithinMatchingGame_Cognition,
                                          ASDOverlapsHigherThanTDCWithinQFamiliar_Cognition,
                                          ASDOverlapsHigherThanTDCWithinQUnfamiliar_Cognition,
                                          DifferencesBetweenGroupsQuestionsHigherMatchingGame_Cognition, 
                                          DifferencesBetweenGroupsUnfamiliarityHigherFamiliarity_Cognition,
                                          
                                          ASD_OverlapsHigherThanZero_Language,
                                          TDC_OverlapsHigherThanZero_Language,
                                          ASD_MG_F_Language,
                                          ASD_Q_F_Language,
                                          ASD_Q_UF_Language,
                                          TD_MG_F_Language,
                                          TD_Q_F_Language,
                                          TD_Q_UF_Language,
                                          ASDOverlapsHigherThanTDCArossConditions_Language,
                                          ASDOverlapsHigherThanTDCWithinMatchingGame_Language,
                                          ASDOverlapsHigherThanTDCWithinQFamiliar_Language,
                                          ASDOverlapsHigherThanTDCWithinQUnfamiliar_Language,
                                          DifferencesBetweenGroupsQuestionsHigherMatchingGame_Language, 
                                          DifferencesBetweenGroupsUnfamiliarityHigherFamiliarity_Language) %>%
    rename("value" = TestDescription) %>%
    mutate(Evid.Ratio = round(Evid.Ratio, 1)) %>%
    mutate(
      Estimate = case_when(Estimate > 5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      Estimate = case_when(Estimate < -5 ~ round(Estimate, 0), TRUE ~ round(Estimate, 2)),
      CI.Lower = case_when(CI.Lower > 5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower > 5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)),
      CI.Lower = case_when(CI.Lower < -5 ~ round(CI.Lower, 0), TRUE ~ round(CI.Lower, 2)),
      CI.Upper = case_when(CI.Lower < -5 ~ round(CI.Upper, 0), TRUE ~ round(CI.Upper, 2)))
  return(HypothesisResultsChildOverlaps)
}
