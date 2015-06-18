.onLoad <- function(lib, pkg)
{
  library.dynam("DMarker", pkg, lib);
}

dmarker = function(data = c("ZZZ3","AAA","ABC", "ZPBP2", "ZPBPL", "ADAMTS4"), pre = "Gene", file = "FALSE", type = "both", Seq = "F")
{
  ktype = 0;
  kgorp = 0;
  kf = 0;
  if(type == "both")
    ktype = 1;
  if( (type == "B") || (type == "blood") || (type == "Blood") || (type == "BLOOD") || (type == "b") )
    ktype = 2;
  if( (type == "U") || (type == "u") || (type == "Urine") || (type == "urine") || (type == "URINE") )
    ktype = 3;
  if( (pre == "G") || (pre == "g") || (pre == "gene") || (pre == "Gene") || (pre == "GENE") )
    kgorp = 1;
  if( (pre == "P") || (pre == "p") || (pre == "protein") || (pre == "Protein") || (pre == "PROTEIN") )
    kgorp = 2;
  if(file != "FALSE")
  {
    kf = 1;
    cat("The results of prediction will write to the file : ");
    cat(file);
    cat(" .\n\n");
  }
  if(ktype == 0)
  {
    cat("The data type can not be recognised.\n");
    ktype = 1;
  }
  if(kgorp == 0)
  {
    cat("The prediction type can not be recognised.\n");
    kgorp = 1;
  }
  cat("DMarker will predict the");
  if(kgorp == 1)
    cat(" GENE ");
  if(kgorp == 2)
    cat(" PROTEIN ");
  cat("whether can be detected");
  if(ktype == 1)
    cat(" both in blood and urine .\n");
  if(ktype == 2)
    cat(" in blood .\n");
  if(ktype == 3)
    cat(" in urine .\n");
  data = as.character(data);
  if( (Seq == "t") || (Seq == "T") || (Seq == "True") || (Seq == "TRUE") || (Seq == "true") )
  {
    if(!require(DMarkerData))
    {
      cat("DMarker needs DMarkerData package to predict the sequence.\nNow trying to install DMarkerData ... \n");
      if(!require(devtools))
      {
        cat("DMarkerData package require devtools package.\nNow installing ... \n");
        install.packages("devtools");
        if(!require(devtools))
        {
          cat("devtools can not be installed.\nDMarker can not predict the sequence.\nPlease use protein or gene names instead.\nBlast can be used for sequence to names.\n");
          return (FALSE);
        }
      }
      require(devtools);
      install_github('yu-shang/DMarkerData');
      if(!require(DMarkerData))
      {
        cat("DMarkerData can not be installed.\nDMarker can not predict the sequence.\nPlease use protein or gene names instead.\nBlast can be used for sequence to names.\n");
        return (FALSE);
      }
    }
    if(!require(seqinr))
    {
      cat("Require seqinr package.\nNow installing ... \n");
      source("http://bioconductor.org/biocLite.R");
      biocLite("seqinr");
      if(!require(seqinr))
      {
        cat("Can not install seqinr from bioconductor.\nDMarker can not predict the sequence.\nPlease use protein or gene names instead.\nBlast can be used for sequence to names.\n");
        return (FALSE);
      }
    }
    require(seqinr);
    require(DMarkerData);
    temp = c();
    data(proseq);
    for(i in 1:length(data))
    {
      tt = s2c(data[i]);
      for(j in 1:length(proseq))
        if(length(tt) == length(proseq[[j]]))
          if(sum(tt != as.character(proseq[[i]])) == 0)
            temp = c(temp, attr(proseq[[j]], "name"));
    }
    data0 = data;
    data = temp;
    kgorp = 2;
  }
  data(uniprot)
  n = length(data);
  uniprot = as.matrix(uniprot);
  rows = nrow(uniprot);
  if(kgorp == 1)
    temp = as.character(uniprot[,2]);
  if(kgorp == 2)
    temp = as.character(uniprot[,1]);
  tempb = temp[which(uniprot[,3] == "1")];
  tempu = temp[which(uniprot[,4] == " 1")];
  result = matrix(0, max(length(tempb), length(tempu)), 2);
  file = c(file, file);
  for(i in 1:n)
  {
    if(length(which(tempb == data[i]) > 0))
      result[which(tempb == data[i]), 1] = 1;
    if(length(which(tempu == data[i]) > 0))
      result[which(tempu == data[i]), 2] = 1;
  }
  tb = tempb[which(result[,1] == 1)];
  tu = tempu[which(result[,2] == 1)];
  results = list();
  if( (ktype == 1) || (ktype == 2) )
  {
    if(length(tb) == 0)
    {cat("No ");cat(pre);cat(" can be detected in BLOOD.\n");}
    cat("\n");
    if(length(tb) > 0)
    {
      cat(length(tb));cat(" ");cat(pre);cat(" can be detected in BLOOD.\n");
      if(kgorp == 1)
      {cat("[GENE]: ");results$BLOOD$GENE = tb;cat(tb);cat(" \n[PROTEIN]: ");results$BLOOD$PROTEIN = as.character(unique(uniprot[which(uniprot[,2]==tb),1]));cat(as.character(unique(uniprot[which(uniprot[,2]==tb),1])));cat(" \n");}
      if(kgorp == 2)
      {cat("[GENE]: ");results$BLOOD$GENE = as.character(unique(uniprot[which(uniprot[,1]==tb),2]));cat(as.character(unique(uniprot[which(uniprot[,1]==tb),2])));cat(" \n[PROTEIN]: ");results$BLOOD$PROTEIN = tb;cat(tb);cat(" \n");}
    }
  }
  if( (ktype == 1) || (ktype == 3) )
  {
    if(length(tu) == 0)
      {cat("No ");cat(pre);cat(" can be detected in URINE\n");}
    if(length(tu) > 0)
    {
      cat(length(tu));cat(" ");cat(pre);cat(" can be detected in URINE.\n");
      if(kgorp == 1)
        {cat("[GENE]: ");results$URINE$GENE = tu;cat(tu);cat(" \n[PROTEIN]: ");results$URINE$PROTEIN = as.character(unique(uniprot[which(uniprot[,2]==tu),1]));cat(as.character(unique(uniprot[which(uniprot[,2]==tu),1])));cat(" \n");}
      if(kgorp == 2)
        {cat("[GENE]: ");results$URINE$GENE = as.character(unique(uniprot[which(uniprot[,1]==tu),2]));cat(as.character(unique(uniprot[which(uniprot[,1]==tu),2])));cat(" \n[PROTEIN]: ");results$URINE$PROTEIN = tu;cat(tu);cat(" \n");}
    }
  }
  rm(uniprot);
  cat("\n\n");
  if(kf == 1)  .C( "r_dmarker", file = as.character(file), BG = as.character(c(results$BLOOD$GENE, "*", "*")), BP = as.character(c(results$BLOOD$PROTEIN, "*", "*")), UG = as.character(c(results$URINE$GENE, "*", "*")), UP = as.character(c(results$URINE$PROTEIN, "*", "*")), iBG = as.integer(length(results$BLOOD$GENE)), iBP = as.integer(length(results$BLOOD$PROTEIN)), iUG = as.integer(length(results$URINE$GENE)), iUP = as.integer(length(results$URINE$PROTEIN)), ktype = as.integer(ktype), kgorp = as.integer(kgorp));
  return(results);
}

demo = function()
{
  results = dmarker(c("ZZZ3","AAA","ABC", "ZPBP2", "ZPBPL", "ADAMTS4"));
  return(results);
}
