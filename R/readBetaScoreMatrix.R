readBetaScoreMatrix <- function(betaScore_file_path)
{
  if (endsWith(betaScore_file_path, suffix = "NA")!=TRUE)
  {
    con <- tryCatch({
      file(betaScore_file_path, "r")
    },
    warning = function(w)
    {
      message(paste0("Warning: Please check your id or internect access!"))
      stop()
    },
    error = function(err)
    {
      message(paste0("Error: Please check your id or internect access"))
    })

    if (!is.null(con))
    {
      # beta score vetor to store all beta scores
      beta_score <- c()
      nposition <- 0
      while(TRUE){
        line <- readLines(con, n=1)
        if (length(line)==0){
          break
        }
        line_vector <- unlist(strsplit(as.character(line), split="\t"))
        nposition <- length(line_vector[-1])
        beta_score <- c(beta_score, line_vector[-1])
      }
      close(con)
      beta_score_matrix <- matrix(as.integer(beta_score), ncol  = nposition, byrow = TRUE)
      beta_score_matrix <- beta_score_matrix[-1,]
      return(beta_score_matrix)
    }
    else
    {
      message(paste0("Beta score matrix for is NULL!"))
      return(NULL)
    }
  }
  else
  {
    return(as.matrix(NA))
  }
}
