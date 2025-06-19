### Function to read in new format metabolon data ###

## read in parameters - choose accordingly
source("parameter_files/04.parameters.R")

## capture output
sink(file=paste0(log_output_dir, "04.read.in.metabolon.new.log"),split=TRUE)

## record todays date
today = Sys.Date()
today = gsub("-","_",today)

## create new directory to save flat text files in
dd = data_output_dir
dd = gsub(" ","\\\\ ", dd)
cmd = paste0("mkdir -p ", dd, "metaboprep_release_", today)
system(cmd)
dir.create(path = paste0(data_output_dir, "metaboprep_release_", today, "/extracted_data"), showWarnings = TRUE, recursive = TRUE, mode = "0777")

read.in.metabolon.new <- function (file2process, data_dir, data_output_dir, projectname) {
  
  ### Is the data file an excel sheet or a flat text file ??
  ftype = c(".xls",".xlsx", ".XLS", ".XLSX")
  isexcel = sum(unlist( sapply(ftype, function(x){ grep(x, file2process) } ) ) )
  
  ## check that you passed an excel sheet
  if( isexcel < 1 ){
    stop( 
      paste0("You must provide the name of the commercial Metabolon excel sheet to process your data using this function.\n"),
      call.=FALSE)
  } else {
    
    ### it is an excel sheet so read in with readxl
    cat(paste0("\t- Your File Being Processed is: ", file2process, "\n"))
    n = paste0(data_dir, file2process)
    
    ## read xls sheet names
    sheetnames = readxl::excel_sheets( n )
    
    
    cat(paste0("\t- There is/are ", length(sheetnames), " sheet(s) in your excel file\n"))
    cat( paste0("\t\t- They are:\n") )
    cat( paste0("\t\t\t- ", sheetnames, "\n") )
    cat( paste0("\t- You have declared that the data being processed is from Metabolon\n") )
    cat( paste0("\t- Reading in your metabolite data\n") )
    
    #########################
    ## TASK 1
    ## read in the data
    #########################
    ## Progress Statement
    cat( paste0("\t- Reading in sheet 4\n") )
    d = readxl::read_excel( n , sheet = 4, na = c("") )
    d = as.data.frame(d)
    
    #########################
    ## TASK 4
    ##  Process Metabolite Data
    ##  on tab 4  
    ##  Peak Area Data
    #########################
    ## Progress Statement
    cat( paste0("\t- Processing Metabolite data on sheet 4: ", sheetnames[4],"\n") )
    
    ## Output the Original Scale metabolite data
    metabolitedata = d
    metabolitedata[metabolitedata == ""] = NA ## turn "" into NAs - think this has already been done in the excel file

    ## redefine column (metabolite) names to COMPIDs
    featuresheet = readxl::read_excel( n , sheet = 2, na = c("") )
    featuresheet = as.data.frame(featuresheet)
    for(i in 1:ncol(featuresheet)){
      featuresheet[,i] = sapply(featuresheet[,i], function(x){ gsub("\r","",x) })
      featuresheet[,i] = sapply(featuresheet[,i], function(x){ gsub("\n","",x) })
    }
    colnames(metabolitedata)[2:ncol(metabolitedata)] = paste0("compid_",featuresheet$COMP_ID)
    
    rownames(metabolitedata) = metabolitedata[,1]
    metabolitedata = metabolitedata[,c(2:ncol(metabolitedata))]
    
    ## write table to file
    metabo_out_name = paste0(data_output_dir, "metaboprep_release_", today, "/extracted_data/", project,"_", today, "_", sheetnames[4], "_Metabolon_metabolitedata.txt")
    metabo_out_name = gsub(" ", "_", metabo_out_name)
    write.table(metabolitedata, file = metabo_out_name, row.names = TRUE, col.names = TRUE, sep = "\t", quote = TRUE)
    
    ## same but for the batch normalised data:
    cat( paste0("\t- Reading in sheet 5\n") )
    d = readxl::read_excel( n , sheet = 5, na = c("") )
    d = as.data.frame(d)
    
    ## Progress Statement
    cat( paste0("\t- Processing Metabolite data on sheet 5: ", sheetnames[5],"\n") )
    
    ## Output the Original Scale metabolite data
    metabolitedata = d
    metabolitedata[metabolitedata == ""] = NA ## turn "" into NAs - think this has already been done in the excel file
    
    ## redefine column (metabolite) names to COMPIDs
    colnames(metabolitedata)[2:ncol(metabolitedata)] = paste0("compid_",featuresheet$COMP_ID)
    
    rownames(metabolitedata) = metabolitedata[,1]
    metabolitedata = metabolitedata[,c(2:ncol(metabolitedata))]
    
    # save this data to an object for later use
    unimputed_data <- metabolitedata
    
    ## write table to file
    metabo_out_name = paste0(data_output_dir, "metaboprep_release_", today, "/extracted_data/", project,"_", today, "_", sheetnames[5], "_Metabolon_metabolitedata.txt")
    metabo_out_name = gsub(" ", "_", metabo_out_name)
    write.table(metabolitedata, file = metabo_out_name, row.names = TRUE, col.names = TRUE, sep = "\t", quote = TRUE)
    
    
    ## same but for the vol-extracted-normalised data:
    cat( paste0("\t- Reading in sheet 8\n") )
    d = readxl::read_excel( n , sheet = 8, na = c("") )
    d = as.data.frame(d)
    
    ## Progress Statement
    cat( paste0("\t- Processing Metabolite data on sheet 8: ", sheetnames[8],"\n") )
    
    ## Output the Original Scale metabolite data
    metabolitedata = d
    metabolitedata[metabolitedata == ""] = NA ## turn "" into NAs - think this has already been done in the excel file
    
    ## redefine column (metabolite) names to COMPIDs
    colnames(metabolitedata)[2:ncol(metabolitedata)] = paste0("compid_",featuresheet$COMP_ID)
    
    rownames(metabolitedata) = metabolitedata[,1]
    metabolitedata = metabolitedata[,c(2:ncol(metabolitedata))]
    
    ### unimpute the data
    ndata <- metabolitedata
    ndata[is.na(unimputed_data)] = NA
    
    ## write table to file
    metabo_out_name = paste0(data_output_dir, "metaboprep_release_", today, "/extracted_data/", project,"_", today, "_", sheetnames[8], "_Metabolon_metabolitedata.txt")
    metabo_out_name = gsub(" ", "_", metabo_out_name)
    write.table(ndata, file = metabo_out_name, row.names = TRUE, col.names = TRUE, sep = "\t", quote = TRUE)
    
    
    
    ## write feature sheet to file
    
    feature_out_name = paste0(data_output_dir, "metaboprep_release_", today, "/extracted_data/", project,"_", today, "_", sheetnames[2], "_Metabolon_featuredata.txt")
    feature_out_name = gsub(" ", "_", feature_out_name)
    write.table(featuresheet, file = feature_out_name, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    
    ## save primary raw metabolite data to a new element
    primary_metabolitedata = ndata
    
    
    #########################
    ## TASK 5
    ##  Process Metabolite Data
    ##  on all other tabs
    #########################
    for(i in 6:7){
      ## Progress Statement
      cat( paste0("\t- Processing Metabolite data on sheet ",i, ": ", sheetnames[i],"\n") )
      
      d = readxl::read_excel( n , sheet = i, na = c("") )
      d = as.data.frame(d)
      
      
      ## Transpose and redefine column (metabolite) names to COMPIDs
      colnames(d) = paste0("compid_",featuresheet$COMP_ID) ## comp_IDs are in sheet 2 column 3
      
      ## write table to file
      outname = paste0(data_output_dir, "metaboprep_release_", today, "/extracted_data/", project,"_", today, "_", sheetnames[i], "_Metabolon_metabolitedata.txt")
      outname = gsub(" ", "_", outname)
      write.table(d, file = outname, row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)
    }
    
    samplesheet = readxl::read_excel( n , sheet = 3, na = c("") )
    samplesheet = as.data.frame(samplesheet)
    
    ## write samplesheet to file
    sample_out_name = paste0(data_output_dir, "metaboprep_release_", today, "/extracted_data/", project,"_", today, "_", sheetnames[3], "_Metabolon_sampledata.txt")
    sample_out_name = gsub(" ", "_", sample_out_name)
    write.table(samplesheet, file = sample_out_name, row.names = TRUE, col.names = TRUE, sep = "\t", quote = TRUE)
    
    ## create featuredata and sampledata flat text files
    sampledata = data.frame(matrix(NA, nrow = nrow(samplesheet), ncol = 3))
    colnames(sampledata) = c("SAMPLE_NAME", "RUN_DAY", "BOX_ID")
    sampledata[,1] = samplesheet$PARENT_SAMPLE_NAME
    sampledata[,3] = samplesheet$BOX_NUMBER
    sampledata_out_name = paste0(data_output_dir, "metaboprep_release_", today, "/extracted_data/", project,"_", today, "_Metabolon_sampledata.txt")
    write.table(sampledata, file = sampledata_out_name, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    
    featuredata = data.frame(matrix(NA, nrow = 1331, ncol = 18))
    colnames(featuredata) = c("feature_names", "SUPER_PATHWAY", "SUB_PATHWAY", "PUBCHEM", "PLATFORM", "KEGG", "HMDB",
                              "CHEM_ID", "LIB_ID", "CHRO_LIB_ENTRY_ID", "PATHWAY_SORTORDER", "TYPE", "INCHIKEY", "SMILES",
                              "CHEMICAL_NAME", "PLOT_NAME", "CAS", "CHEMSPIDER")
    featuredata[,1] = paste0("compid_", featuresheet$COMP_ID)
    featuredata[,2] = featuresheet$SUPER_PATHWAY
    featuredata[,3] = featuresheet$SUB_PATHWAY
    featuredata[,4] = featuresheet$PUBCHEM
    featuredata[,5] = featuresheet$PLATFORM
    featuredata[,6] = featuresheet$KEGG
    featuredata[,7] = featuresheet$HMDB
    featuredata[,8] = featuresheet$CHEM_ID
    featuredata[,9] = featuresheet$LIB_ID
    featuredata[,10] = featuresheet$CHRO_LIB_ENTRY_ID
    featuredata[,11] = featuresheet$PATHWAY_SORTORDER
    featuredata[,12] = featuresheet$TYPE
    featuredata[,13] = featuresheet$INCHIKEY
    featuredata[,14] = featuresheet$SMILES
    featuredata[,15] = featuresheet$CHEMICAL_NAME
    featuredata[,16] = featuresheet$PLOT_NAME
    featuredata[,17] = featuresheet$CAS
    featuredata[,18] = featuresheet$CHEMSPIDER
   ## remove \r (extra lines in some cells due to carriage return in excel)
     for(i in 1:ncol(featuredata)){
      featuredata[,i] = sapply(featuredata[,i], function(x){ gsub("\r","",x) })
      featuredata[,i] = sapply(featuredata[,i], function(x){ gsub("\n","",x) })
    }
    
    featuredata_out_name = paste0(data_output_dir, "metaboprep_release_", today, "/extracted_data/", project,"_", today, "_Metabolon_featuredata.txt")
    write.table(featuredata, file = featuredata_out_name, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
    
    ### return data to user
    mydata = list(metabolitedata = primary_metabolitedata, sampledata = sampledata, featuredata = featuredata )
    return(mydata)
    
  }  ### END OF ELSE STATAMENT, yes you passed the function an excel sheet
  
} ### END OF READ IN RAW METABOLON DATA FUNCTION

## run on data - to make flat text files to feed into metabolon
read.in.metabolon.new(file2process = file2process, data_dir = data_dir, data_output_dir = data_output_dir, projectname = project)

# capture session info
print("Session information:")
sessionInfo()

# save output
sink()
