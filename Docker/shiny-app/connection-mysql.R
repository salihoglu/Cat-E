library(RMySQL)

con <-reactive({
  dbConnect(RMySQL::MySQL(),
                 user = 'Tarcin',
                 host = 'db',
                 port = 3306,
                 dbname='CATE',
                 password = "DumanRS")



}) 
#on.exit(dbDisconnect(con), add = TRUE)

