library(RSQLite)
library(data.table)
sqlite    <- dbDriver("SQLite")
db <- dbConnect(sqlite,'/Volumes/Seagate Backup Plus Drive/Data Tables/MM Stromal Clustering/Whole Dataset copy 2/DefaultDB.db')

# Fix NA values in Table
dbListTables(db)
duh <- dbReadTable(db, 'sqlite_sequence')
duh[duh$field == 'db_sqlite_file','value'] <- 'DefaultDB.db'
dbWriteTable(db, 'Experiment_Properties', duh, overwrite=T)

myTable <- data.table(dbReadTable( db, 'MyExpt_Per_Object' ))
badCols <- names(myTable)[as.logical(as.vector(myTable[,lapply(.SD, function(x){length(which(is.na(x))) > 0}),]))]
badRows <- paste0(myTable$ImageNumber, '_', myTable$ObjectNumber)[is.na(myTable$Nuc_AreaShape_Area)]
myTable <- myTable[!is.na(Nuc_AreaShape_Area)]
badCols <- names(myTable)[as.logical(as.vector(myTable[,lapply(.SD, function(x){length(which(is.na(x))) > 0}),]))]
dbWriteTable(db, 'MyExpt_Per_Object', myTable, overwrite=T, row.names=F)

# Fix Image Paths
myTable <- data.table(dbReadTable(db, 'MyExpt_Per_Image'))
pathCols <- names(myTable)[grepl('Path',names(myTable))]
fixPaths <- function(x)
{
     y <- gsub('\\','/',x, fixed=T)
     z <- gsub('H:','/Volumes/Seagate Backup Plus Drive',y, fixed=T)
     return(z)
}
fixedCols <- myTable[,lapply(.SD, fixPaths), .SDcols=pathCols]
myTable[,c(pathCols)] <- fixedCols
dbWriteTable(db, 'MyExpt_Per_Image', myTable, overwrite=T, row.names=F)

dbDisconnect(db)

# Do classification in CellProfiler Analyst
# Score All

# Read assigned CellTypes
dbListTables(db)
dbListFields(db, "MyExpt_Per_Object")
myTable <- data.table(dbReadTable( db, 'MyExpt_Per_Object' ))
badCols <- names(myTable)[as.logical(as.vector(myTable[,lapply(.SD, function(x){length(which(is.na(x))) > 0}),]))]
badRows <- paste0(myTable$ImageNumber, '_', myTable$ObjectNumber)[is.na(myTable$Nuc_AreaShape_Area)]
myTable <- myTable[!is.na(Nuc_AreaShape_Area)]
badCols <- names(myTable)[as.logical(as.vector(myTable[,lapply(.SD, function(x){length(which(is.na(x))) > 0}),]))]
dbWriteTable(db, 'MyExpt_Per_Object', myTable, overwrite=T, row.names=F)

