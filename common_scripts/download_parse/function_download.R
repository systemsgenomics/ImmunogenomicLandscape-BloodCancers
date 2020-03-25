getFTPfiles=function(FTP, PATTERN="*"){
  files_found=scan(pipe(paste0("lftp -c 'connect ", FTP, "; find | grep ", PATTERN, "'")), "a", quiet = T)
  files_found_abs=paste0(FTP, gsub("^\\.", "", files_found))
}

download.file.ftp=function(filep, out=getwd()){
  download.file(filep, destfile = file.path(out, sapply(strsplit(filep, "/"), tail, 1)))
}
