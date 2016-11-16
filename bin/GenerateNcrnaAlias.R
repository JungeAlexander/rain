# Script for generation of non-miR non-coding RNA alias file
# Author: Christian Garde, garde@sund.ku.dk, 
# Version 1, Nov 1 2016
#------------------------------------------
options(warn=1)

dedup_alias = function(x){
    #----------------------
    # Remove duplicates if any
    #----------------------
    x$dedup_id = paste(x$RAIN,x$Alias,sep="_")
    x = data.table(x,key="dedup_id")
    x = unique(x,by = "dedup_id")
    x$dedup_id=NULL
    return(data.frame(x,stringsAsFactors = F))
}

rm_spurious = function(x){
    #----------------------
    # Remove entries that point to multiple RAIN.IDs
    #----------------------
    spurious.ids = table(x$Alias)
    spurious.ids = names(spurious.ids[spurious.ids>1])
    x = x[!(x$Alias %in% spurious.ids),]
    return(x)
    
}

sanity_correction = function(x){
    return(rm_spurious(dedup_alias(x)))
}

return_spurious = function(x){
    x = dedup_alias(x)
    x[x$RAIN=="FLJ35934",]
    spurious.ids = table(x$Alias)
    spurious.ids = names(spurious.ids[spurious.ids>1])
    x = x[(x$Alias %in% spurious.ids),]
    return(x[order(x$Alias),])
}

addmissing2alias = function(x){
    tmp=unique(x[!(x$RAIN %in% x$Alias),"RAIN"])
    if(length(tmp)>0){
        x = rbind(x,
                  data.frame(Texonomy=taxonomy,
                             RAIN=tmp,Alias=tmp,Database="RAIN", 
                             stringsAsFactors = F)
        )
    }
    return(x)
}


rm_paralog_trail = function(x){
    #-------------------------
    # Has checked that that only paralogs end with dash or underscore followed by digit, hence remove that trailing pattern
    #------------------------- 
    #unique(alias.list$RAIN[grepl("-[0-9]*$",alias.list$RAIN)])
    #unique(alias.list$RAIN[grepl("_[0-9]*$",alias.list$RAIN)])
    x$RAIN = gsub("-[0-9]*$","",x$RAIN)
    x$RAIN = gsub("_[0-9]*$","",x$RAIN)
    x = addmissing2alias(x)
    return(x)
}

paralog_correction = function(alias.list,ensembl,
                              paralog_db_attr="hsapiens_paralog_ensembl_gene",
                              paralog_conf_attr="hsapiens_paralog_paralogy_confidence"
                              ){
    #----------------------
    # Correct for paralogs
    # THIS FUNCTION DOES NOT SOLVE THE ISSUE DUE TO FLAWED ENSEMBLE PARALOG ANNOTRATON.
    # DO NOT USE
    #----------------------

    # Retrieve paralog list
    biomart_table3 = getBM(attributes = c("external_gene_name",
                                          'ensembl_gene_id',
                                          paralog_db_attr,
                                          paralog_conf_attr
                                          ),
                           mart = ensembl)
    
    biomart_table3 = biomart_table3[biomart_table3$ensembl_gene_id %in% alias.list$Alias &
                                    biomart_table3$hsapiens_paralog_ensembl_gene %in% alias.list$Alias &
                                    biomart_table3$hsapiens_paralog_paralogy_confidence == 1 &
                                    !is.na(biomart_table3$hsapiens_paralog_paralogy_confidence),]
    
    # Rename RAIN.ids of secondary paralogs, primnary paralog defined by sort
    # but prioritizing curated ids, then database specific id, and finallaly others.
    # Better ideas for prioritizing paralog ids for RAIN ids are welcome!
    
    biomart_table3 = data.table(biomart_table3)
    biomart_table3$keep = T
    setkey(biomart_table3,external_gene_name)
    
    biomart_table4 = data.table(biomart_table3)
    setkey(biomart_table4,ensembl_gene_id)
    
    curated_ids  = c("5S_rRNA","5_8S_rRNA","18S_rRNA","28S_rRNA","45S_rRNA")
    
    RAIN.ids = c(curated_ids,
                 sort(unique(alias.list$RAIN[alias.list$Database==taxon.db.name])),
                 sort(unique(alias.list$RAIN[alias.list$Database!=taxon.db.name])))
    
    alias.list=data.table(alias.list,key=c("RAIN","Alias"))
    for(RAIN.id in RAIN.ids){
    
        tmp = biomart_table3[ RAIN.id ]
        tmp = tmp[tmp$keep]
        tmp = tmp[!rowSums(is.na(tmp))==ncol(tmp)]
    
        if(nrow(tmp) == 0 ){
            next
        }
    
        Ensembl2set = unique(c(tmp$hsapiens_paralog_ensembl_gene,tmp$ensembl_gene_id))
        GeneName2set = unique(biomart_table4[ Ensembl2set ]$external_gene_name)
    
        setkeyv(biomart_table3[ GeneName2set, "keep":=F ],key(biomart_table3))
        setkeyv(alias.list[ GeneName2set,RAIN:=RAIN.id],c("RAIN","Alias"))
    }
    return(alias.list)
}

return_manualdb = function(alias.list,ManualDB,taxonomy){
    
    manual.aliases = read.table(ManualDB,sep="\t", header=F, stringsAsFactors = F)
    colnames(manual.aliases) = c("Texonomy","RAIN","Alias","Database")
    manual.aliases = manual.aliases[manual.aliases$Texonomy==taxonomy,]
    manual.aliases = dedup_alias(manual.aliases)
    
    # Quality check
    # return_spurious(manual.aliases)
    
    # Update current alias list based on manual annotation
    ids2correct = unique(alias.list[alias.list$RAIN %in% manual.aliases$Alias,"RAIN"])
    for(id in ids2correct){
        alias.list[alias.list$RAIN == id,"RAIN"] =  manual.aliases[manual.aliases$Alias == id,"RAIN"]
    }
    
    # Add manual db to alias list and remove potential duplicated entries
    alias.list = rbind(alias.list,manual.aliases)
    alias.list = dedup_alias(alias.list)

    return( alias.list )
}

retrieve_biomart = function(taxonomy="9606",taxon.db="hsapiens_gene_ensembl", 
                            taxon.attr="hgnc_id",taxon.db.name="HGNC",
                            ManualDB = "ManualAliasDict.tsv",
                            HGNC_noncoding="ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_groups/non-coding_RNA.txt"){

    #----------------------
    # Retrieve biomart lists
    #----------------------
    ensembl <- useEnsembl('ensembl',dataset=taxon.db, version = 79)
    
    taxon_special_treatment = c("4932","6239")
    if(taxonomy=="4932"){
        table1_attr = c("external_gene_name",'ensembl_gene_id', taxon.attr,'description')
        db.list =c("Ensembl","Ensembl",taxon.db.name,"Ensembl_Description")
    }else if(taxonomy %in% c("6239","7227")){
        table1_attr = c("external_gene_name",'ensembl_gene_id', taxon.attr,
                        'refseq_ncrna','description')
        db.list =c("Ensembl","Ensembl",taxon.db.name,"RefSeq","Ensembl_Description")
    }else{
        table1_attr = c("external_gene_name",'ensembl_gene_id', taxon.attr,
                        'refseq_ncrna_predicted','refseq_ncrna','description')
        db.list =c("Ensembl","Ensembl",taxon.db.name,"RefSeq","RefSeq","Ensembl_Description")
    }
    
    if(taxonomy %in% taxon_special_treatment){
        table2_attr = c("external_gene_name",'gene_biotype','ensembl_gene_id')
    }else{
        table2_attr = c("external_gene_name",'gene_biotype','ensembl_gene_id',"mirbase_accession")
    }
    
    biomart_table1 = getBM(attributes = table1_attr,mart = ensembl)
    biomart_table1$description = trimws(gsub(" $","",gsub("\\[Source:.*","",biomart_table1$description)),"both")
    
    biomart_table2 = getBM(attributes = table2_attr,mart = ensembl)

    biotypeclass2keep =  c("Mt_tRNA","Mt_rRNA","lincRNA","antisense",
                           "TEC","snoRNA","misc_RNA","rRNA","snRNA",                                                          
                           "sRNA","ribozyme","sense_overlapping",
                           "3prime_overlapping_ncrna","non_coding",                         
                           "vaultRNA", "macro_lncRNA","scaRNA")
    
    GenesToKeep = (if (taxonomy %in% taxon_special_treatment) T else biomart_table2$mirbase_accession=="") &
        biomart_table2$gene_biotype %in% biotypeclass2keep &  
        biomart_table2$ensembl_gene_id != "" 
    
    biomart_table1 = biomart_table1[biomart_table1$external_gene_name %in% biomart_table2[GenesToKeep,"external_gene_name"],]
    biomart_table1 = biomart_table1[!grepl("pseudogene",biomart_table1$description),] 
    
    #----------------------
    # Assemble alias list
    #----------------------
    alias.list = cbind(rep(taxonomy,ncol(biomart_table1)*nrow(biomart_table1)),
                       rep(biomart_table1[,1],each=ncol(biomart_table1)) ,
                       unlist(lapply(t(biomart_table1) , function(x) as.character(x) )),
                       rep(db.list,nrow(biomart_table1)))
    rownames(alias.list) = NULL
    alias.list = data.frame(alias.list[alias.list[,2]!="" & alias.list[,3]!="",],stringsAsFactors=F)
    colnames(alias.list) = c("Texonomy","RAIN","Alias","Database")
    if(!(taxonomy %in% c("6239","7227"))){ # for worm, I cannot tell whether these are separate gnees or paralogs
        alias.list = rm_paralog_trail(alias.list)
    }
    alias.list = sanity_correction(alias.list)
    alias.list = alias.list[!is.na(alias.list$RAIN),]
    
    #-----------------------
    # Complement with hgnc table and manual annotation since 
    # ensembl seem to miss some curated targets
    #-----------------------
    
    if(taxonomy=="9606"){
        
        # manual correction -- sub-optimal solution but required in order to avoid mis-matches (spurious links)
        alias.list[alias.list$RAIN %in% c("RNU1","RNU1-1","RNU1-2","RNU1-3","RNU1-4"),"RAIN"] = "U1"
        alias.list[alias.list$RAIN  %in% c("RNU2","RNU2-1","snoU2"),"RAIN"] = "U2"
        alias.list[alias.list$RAIN  %in% c("RNU4","RNU4-2","RNU4-1"),"RAIN"] = "U4"
        alias.list[alias.list$RAIN  %in% c("RNU5A-1","RNU5B-1","RNU5D-1","RNU5E-1","RNU5F-1"),"RAIN"] = "U5"
        alias.list[alias.list$RAIN  %in% c("RNU5A","RNU5B","RNU5D","RNU5E","RNU5F"),"RAIN"] = "U5"
        alias.list[alias.list$RAIN  %in% c("RNU6","RNU6-1","RNU6-2","RNU6-7","RNU6-8","RNU6-9"),"RAIN"] = "U6"
        alias.list[alias.list$RAIN  %in% c("RNU7","RNU7-1"),"RAIN"] = "U7"    
        alias.list[alias.list$RAIN  == "RNU11","RAIN"] = "U11"    
        alias.list[alias.list$RAIN  == "RNU12","RAIN"] = "U12"    
        alias.list[alias.list$RAIN  ==  "RP11-145M4.3","RAIN"] = "LINC01358"
        alias.list[alias.list$RAIN  ==  "RP11-166D19.1","RAIN"] = "MIR100HG"
        alias.list[alias.list$RAIN  ==  "RNase_MRP","RAIN"] = "RMRP"
        alias.list[alias.list$RAIN  ==  "RNaseP_nuc","RAIN"] = "RPPH1"
        alias.list[alias.list$RAIN  ==  "RN7SK","RAIN"] = "7SK"
        alias.list[alias.list$RAIN  ==  "U3","RAIN"] = "SNORD3A"
        alias.list = return_manualdb(alias.list, ManualDB,taxonomy)
        
        
        # expand the known aliases based on hgnc symbol
        if (grepl('\\.gz$', HGNC_noncoding)) {
          db = suppressWarnings(fread(paste("zcat", HGNC_noncoding)))
        } else {
          db = suppressWarnings(fread(HGNC_noncoding))
        }
        db = as.data.frame(db[db$gene_family!="MicroRNAs" & db$mirbase=="" & db$status=="Approved"])
        
        tmp = unlist(lapply(db$alias_symbol, function(x) unlist(strsplit(x,"|",fixed=T))))
        tmp = table(tmp)
        tmp = names(tmp[tmp>1])

        tmp2 = unlist(lapply(db$prev_symbol, function(x) unlist(strsplit(x,"|",fixed=T))))
        tmp2 = table(tmp2)
        tmp2 = names(tmp2[tmp2>1])
        
        disqualified_ids = unique(c(tmp,tmp2))
      
        for(i in 1:nrow(db)){
            record = as.character(db[i,c("hgnc_id","symbol","alias_symbol","prev_symbol","ensembl_gene_id")])
            record = unlist(lapply(record, function(x) unlist(strsplit(x,"|",fixed=T))))
        
            if(length(record)>3){
                record = c(record[1:2],record[3:length(record)][ !(record[3:length(record)] %in% disqualified_ids) ])
            }
            
            
            RAIN.id = unique(alias.list[alias.list$Alias %in% record,"RAIN"])
            
            if( record[2] %in% paste("RNA5-8S",1:5,sep="") ){
                    RAIN.id = "5_8S_rRNA"  
                    alias.list[alias.list$RAIN %in% record[2],"RAIN" ] = RAIN.id
            }else if( record[2] %in% paste("RNA5S",1:17,sep="") ){
                    RAIN.id = "5S_rRNA" 
                    alias.list[alias.list$RAIN %in% record[2],"RAIN" ] = RAIN.id
            }else if( record[2] %in% paste("RNA18S",1:5,sep="") ){
                    RAIN.id = "18S_rRNA"
                    alias.list[alias.list$RAIN %in% record[2],"RAIN" ] = RAIN.id
            }else if( record[2] %in% paste("RNA28S",1:5,sep="") ){
                    RAIN.id = "28S_rRNA"
                    alias.list[alias.list$RAIN %in% record[2],"RAIN" ] = RAIN.id
            }else if( record[2] %in% paste("RNA45S",1:5,sep="") ){
                    RAIN.id = "45S_rRNA"
                    alias.list[alias.list$RAIN %in% record[2],"RAIN" ] = RAIN.id
            }else if(length(RAIN.id)>1){
                tmp = gsub("-[0-9]*$","",record[2])
                tmp = gsub("_[0-9]*$","",tmp)
                alias.list[alias.list$RAIN %in% RAIN.id,"RAIN" ] = tmp
                RAIN.id = tmp
            }else{
                    RAIN.id = RAIN.id 
            }
            
            # add new entry
            if(length(RAIN.id)==0){
                RAIN.id = record[2]
                RAIN.id = gsub("-[0-9]*$","",RAIN.id)
                RAIN.id = gsub("_[0-9]*$","",RAIN.id)
            }
            
            record = c(RAIN.id,record)
                
            new_alias = data.frame(cbind("9606",
                              rep(RAIN.id, each=length(record)),
                              rep(record,length(RAIN.id)),
                              "HGNC"),stringsAsFactors = F)
            colnames(new_alias) = c("Texonomy","RAIN","Alias","Database")
            
            alias.list = rbind(new_alias,alias.list) # this step could be optimized through preallocation potentially
            
        }
        
      
    }else if(taxonomy=="10090"){
        alias.list[alias.list$RAIN %in% paste("n-R5s",1:500,sep=""),"RAIN"] = "Rn5s"
    }else if(taxonomy=="10116"){
        alias.list$RAIN[alias.list$RAIN=="5S_rRNA"] = "Rn5s"
    }else if(taxonomy=="6239"){
        alias.list = alias.list[alias.list$RAIN!="mir-1831",]   
    }else{
        Q=1
        #do nothing?!?
    }
    # Update with manual annotations
    alias.list = return_manualdb(alias.list, ManualDB,taxonomy)
    alias.list = sanity_correction(alias.list)
    
    # Ensure that main ids always have a pointer to themselves
    alias.list = addmissing2alias(alias.list)
    
    if( sum(!(unique(alias.list$RAIN) %in% alias.list$Alias))>0 ){
        print(paste(taxonomy,": We have a problem"))
    }
    
    return(alias.list)
    
}

print_helppage = function(){
cat(
"------------------------------------
Generation of ncRNA (non-miR) alias file
Version 1, November 3, 2016
------------------------------------

Positional Arugments:

ManualDB         Path to the manually curated aliases
CircDB           Path to the CircDB aliases
HGNC_noncoding   Path to the HGNC non-coding aliases
OutFile          Path to output file     
")
}

message('Started generating non-miR non-coding RNA alias file.')
design = cbind(
    c(9606,10090,10116,4932,6239,7227,7955),
    c('hsapiens_gene_ensembl', 'mmusculus_gene_ensembl','rnorvegicus_gene_ensembl',
      "scerevisiae_gene_ensembl","celegans_gene_ensembl","dmelanogaster_gene_ensembl",
      "drerio_gene_ensembl"),
    c("hgnc_id",'mgi_id',"rgd","sgd_gene","wormbase_gene","flybase_gene_id","zfin_id"),
    c("HGNC","MGI","RGD","SGD","WormBase","FlyBase","ZFIN")
    )
colnames(design) = c("Taxon","biomartdb","orgndb","orgdbname")

args = commandArgs(trailingOnly=TRUE)
ManualDB = args[1]
CircDB = args[2]
HGNC_noncoding = args[3]
OutFile = args[4]

if(!(ManualDB %in% c("-h","--help")) & length(args) < 4){
    cat(
"------------------------------------------------
>>>> Please set the positional arugments <<<<<
------------------------------------------------\n"
)
    print_helppage()
    quit(save = "no", status = 1, runLast = T)
}else if(ManualDB %in% c("-h","--help")){
    print_helppage()
    quit(save = "no", status = 0, runLast = T)
}else{
    
    library(biomaRt)
    library(data.table)

    alias.lists = list()
    for(i in 1:nrow(design)){
        taxonomy=as.character(design[i,1])
        taxon.db=as.character(design[i,2])
        taxon.attr=as.character(design[i,3])
        taxon.db.name=as.character(design[i,4])
        alias.lists[[taxonomy]] = retrieve_biomart(taxonomy,taxon.db, 
                                        taxon.attr,taxon.db.name,
                                        ManualDB,HGNC_noncoding)
    }

    alias.lists = do.call("rbind", alias.lists)
    rownames(alias.lists) = NULL
    
    CircDB.list = fread(paste("zcat",CircDB),sep="\t", header = F, 
                        stringsAsFactors = F, data.table=F)
    colnames(CircDB.list) = c("Texonomy","RAIN","Alias","Database")
    alias.lists = rbind(alias.lists, CircDB.list)
    
    if (grepl('\\.gz$', OutFile)) {
      gz1 <- gzfile(OutFile, "w")
      write.table(alias.lists,file=gz1,col.names = F, row.names = F,quote=F, sep="\t")
      close(gz1)
    } else {
      write.table(alias.lists,file=OutFile,col.names = F, row.names = F,quote=F, sep="\t")
    }
}
message('Done.')
