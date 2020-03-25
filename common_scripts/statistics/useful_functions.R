aggregate_by_mean <- function(data, xs) {
  if(dim(data)[2]!=length(xs))stop("length between data and variable 'xs' do not match!")
  library(reshape)
  library(data.table)
  # Convert to data.table.
  dat <- data.table(t(data))
  # Append the vector of group names as an extra column.
  dat$agg_var <- as.character(xs)
  # Melt the data.table so all values are in one column called "value".
  dat <- melt(dat, id.vars = "agg_var")
  # Cast the data.table back into the original shape, and take the mean.
  dat <- dcast.data.table(
    dat, agg_var ~ variable, value.var = "value",
    fun.aggregate = mean, na.rm = TRUE
  )
  
  namecol=dat$agg_var
  
  # Delete the extra column.
  dat[ , agg_var := NULL]
  
  dat=t(dat)
  
  rownames(dat) <- rownames(data)
  colnames(dat)=namecol
  return(dat)
}

aggregate_by_median <- function(data, xs) {
  if(dim(data)[2]!=length(xs))stop("length between data and variable 'xs' do not match!")
  library(reshape)
  library(data.table)
  # Convert to data.table.
  dat <- data.table(t(data))
  # Append the vector of group names as an extra column.
  dat$agg_var <- as.character(xs)
  # Melt the data.table so all values are in one column called "value".
  dat <- melt(dat, id.vars = "agg_var")
  # Cast the data.table back into the original shape, and take the mean.
  dat <- dcast.data.table(
    dat, agg_var ~ variable, value.var = "value",
    fun.aggregate = median, na.rm = TRUE
  )
  
  namecol=dat$agg_var
  
  # Delete the extra column.
  dat[ , agg_var := NULL]
  
  dat=t(dat)
  
  rownames(dat) <- rownames(data)
  colnames(dat)=namecol
  return(dat)
}