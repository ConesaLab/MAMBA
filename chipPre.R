# ChIP-seq data transformation to reaction-level

# Read the GPRs of your metabolic model
grRules <- read.table("gpr_complete_model.csv", header = TRUE, sep = ";") 

# Load your ChIP-seq data matrix
ch_data <- load("ChIP.RData")

# Compute GPRs to get reaction-level values
ch_vals <- NULL
for ( j in 1:dim(ch_data)[2]) {
  
  chip_data <- ch_data[,j,drop = FALSE]
  #exp_data <- rna_data[,j,drop = FALSE]
  ch_vals_temp <- c()
  #rna_vals_temp <- c()
  
  for (x in 1:dim(grRules)[1]) {
    print(x)
    r <- paste0("R_",as.character(grRules[x,1]))
    gpr <- as.character(grRules[x,2])
    
    # GPR
    if (!(is.na(gpr) || is.null(gpr) || gpr == "")) {
      chip <- getValues(gpr, chip_data)
      resch <- as.numeric(eval(chip))
    } else {
      resch <- NA
    }
    ch_vals_temp <- c(ch_vals_temp, resch)
  }
  ch_vals <- cbind(ch_vals, ch_vals_temp)
}

# Compute logFC between conditions and save the table

################################################################################
## Aux function

# Get numerical value of a gpr
getValues <- function(gpr, table) {
  str2 <- strsplit(gpr, " ")[[1]]
  #browser()
  for ( i in 1:length(str2)) {
    if (!is.element(str2[i], c("(", ")","and", "or", "UNKNOWN"))) {
      if (str2[max(1,i-1)]=="(") {
        #str3 <- strsplit(str2[i], "\\(")[[1]]
        name <- str2[i]
        if(length(grep("@", name))>0) {
          str4 <- strsplit(name, "@")[[1]]
          name <- checkName (str4[2])
          #eid <- table[which(table[,1] == name), 5]
          eid <- table[name, 1]*1/as.numeric(str4[1])
        } else {
          name <- checkName (str2[i])
          eid <- table[name, 1]
        }
        str2[i] <- eid
      } else if (str2[min(length(str2),i+1)]==")") {
        name <- str2[i]
        if(length(grep("@", name))>0) {
          str4 <- strsplit(name, "@")[[1]]
          name <- checkName (str4[2])
          #eid <- table[which(table[,1] == name), 5]
          eid <- table[name, 1]*1/as.numeric(str4[1])
        } else {
          #str3 <- strsplit(str2[i], "\\)")[[1]]
          name <- checkName (str2[i])
          eid <- table[name, 1]
        }
        str2[i] <- eid
        #str2[i] <- paste0(paste(str3, collapse = ")"), ")")
      } else {
        name <- str2[i]
        if(length(grep("@", name))>0) {
          str4 <- strsplit(name, "@")[[1]]
          name <- checkName (str4[2])
          #eid <- table[which(table[,1] == name), 5]
          eid <- table[name, 1]*1/as.numeric(str4[1])
        } else {
          name <- checkName (str2[i])
          eid <- table[name, 1]
        }
        str2[i] <- eid
      }
    }
  }
  return(paste(str2, collapse = " "))
}

# Solve boolean relationships between genes
eval <- function(exp) {
  
  if (str_count(exp, "\\(")>0) {
    open <- str_locate_all(exp, "\\(")[[1]]
    close <- str_locate_all(exp, "\\)")[[1]]
    table <- NULL
    i = 1
    while ( i <= dim(close)[1]) {
      sel <- which(open[i:dim(open)[1],1] < close[i,1])
      table <- rbind(table, c(min(open[i:dim(open)[1],1][sel]), max(close[i:dim(open)[1],1][sel])))
      i = i + length(sel)
    }
    parts <- apply(table, 1, function(x) {
      eval(substr(exp, x[1]+1, x[2]-1))
    })
    gval <- parts[1]
    if ( dim(table)[1] > 1) {
      for ( i in 1:(dim(table)[1]-1)) {
        subs <- substr(exp, table[i,2]+1, table[i+1,1] - 1)
        if (str_count(subs, "and")>0) {
          gval <- min(gval, parts[i+1])
          if ( gval == 0) {
            gval <- max(gval, parts[i+1])
          }
        } else if (str_count(exp, "or")>0) {
          gval <- max(gval, parts[i+1])
        }
      } 
    }
    if (table[dim(table)[1],2]+1 < nchar(exp)) {
      subs <- substr(exp, table[dim(table)[1],2]+1, nchar(exp))
      if (str_count(subs, "and")>0) {
        gval <- min(gval, as.numeric(strsplit(subs, "and")[[1]]), na.rm = TRUE)
        if ( gval == 0) {
          gval <- max(gval, as.numeric(strsplit(subs, "and")[[1]]), na.rm = TRUE)
        }
      } else if (str_count(exp, "or")>0) {
        gval <- max(gval, as.numeric(strsplit(subs, "and")[[1]]), na.rm = TRUE)
      }
    }
    gval
    
  } else {
    if (str_count(exp, "and")>0) {
      res <- min(as.numeric(strsplit(exp, "and")[[1]]), na.rm = TRUE)
      if ( res == 0) {
        res <- max(as.numeric(strsplit(exp, "and")[[1]]), na.rm = TRUE)
      }
      res
    } else if (str_count(exp, "or") >0 ) {
      max(as.numeric(strsplit(exp, "or")[[1]]), na.rm = TRUE)
    } else {
      as.numeric(exp)
    }
  }
}