#make table of which ctenophore sequences have mismatches to coral primers and where

library("tidyverse")
setwd("/Users/quattrinia/Desktop/AW RStudio/data/primer-matches-sponges-ctenos") #set working directory to folder with sponge data

#import files
ogcteno <- read_csv('28S_Ctenos_with_Outgroups.csv') #import table of original ctenophore sequences
cteno.barcode <- read_csv('28S_Ctenos_with_Outgroups_barcode.csv', show_col_types = FALSE) # import table of ctenophore sequences containing barcode
cteno.trimprim <- read_csv('28S_Cteno_trimprimers.csv', show_col_types = FALSE) #import table of ctenophore sequences with primers trimmed out
cteno.retainprim <- read_csv('28S_Cteno_retainprimers.csv', show_col_types = FALSE) #import table of ctenophore sequences with primer mismatch locations

#create barcode column
table.barcode <- left_join(ogcteno, cteno.barcode, by = "Name") #combined table of original ctenophore sequences with ctenophore sequences containing barcode
table.barcode$barcodematch <- ifelse(is.na(table.barcode$Sequence.y), "no", "yes") #create column saying no if sequence doesn't have barcode and yes if sequence has it

#create primer match column
table.primers <- left_join(table.barcode, cteno.trimprim, by = "Name") #combined working table with ctenophore sequences that match primers with <= 2 mismatches
table.primers$primermatch <- ifelse(is.na(table.primers$Sequence), "no", "yes") #create column saying no if sequence doesn't match primers and yes if it does

#create mismatches column for F and R primer
table.mismatch <- filter(cteno.retainprim, Type == "mismatch") #filter out the rows that have mismatch in the type column from the primer table
table.mismatchinfo <- table.mismatch %>%
  group_by(`Sequence Name`) %>% #grouping the sequences by primer name since some have more than one mismatch
  summarise(Number=sum(Length), Positions=toString(unique(Minimum))) #summarizing by amount and location of mismatches
table.mismatchinfo <- table.mismatchinfo %>% rename("Name" = "Sequence Name") #rename name column to match working table
table.ctenosummary <- left_join(table.primers, table.mismatchinfo, by = "Name") #combine column of mismatch information with working table
view(table.ctenosummary)

#polishing touches
table.ctenosummary$NumberMismatch <- ifelse(table.ctenosummary$barcodematch == "yes", 
                                     ifelse(table.ctenosummary$primermatch == "yes", ifelse(is.na(table.ctenosummary$Number), "0", table.ctenosummary$Number), table.ctenosummary$Number), 
                                     table.ctenosummary$Number) #if the barcode is present and the primer matches but there are no mismatches, change the NA values to 0
view(table.ctenosummary)

#export final table as csv
write.csv(table.ctenosummary, "/Users/quattrinia/Desktop/AW RStudio/28Sctenoprimermatches.csv") 


