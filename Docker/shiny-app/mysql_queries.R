e<-reactiveValues(data=NULL)


observe({
  multipleInput()
})




multipleInput<- eventReactive(input$searchProtein,{
  if(input$searchProtein > 0){
    withProgress(message = "Running...",{
    
  isolate({
  dbGetQuery(con(), "USE CATE;")
  dbGetQuery(con(), "SET SESSION max_statement_time = 900;")
  dbGetQuery(con(), "DROP VIEW IF EXISTS temptable;")


  query1<-paste0(
    "CREATE VIEW temptable AS          
SELECT *
FROM Gene2Gene
WHERE GeneInfoIDfrom IN (
    SELECT UniprotKB
    FROM GeneInfo
    WHERE Gene_Symbol = '",input$searchProtein,"' AND Interaction_ID IN ('activation', 'inhibition')
) OR GeneInfoIDto IN (
    SELECT UniprotKB
    FROM GeneInfo
    WHERE Gene_Symbol = '",input$searchProtein,"' AND Interaction_ID IN ('activation', 'inhibition')
);")
  
  qmain1 <-dbGetQuery(con(), query1)
  
  dbGetQuery(con(), "DROP TABLE IF exists set1;")
  query2<-paste0("               
CREATE TEMPORARY TABLE set1 AS
SELECT 
    t.GeneInfoIDfrom,
    t.Interaction_ID,
    t.GeneInfoIDto,
    g.Gene_Symbol AS Gene_1
FROM 
    temptable AS t
    INNER JOIN GeneInfo AS g ON t.GeneInfoIDfrom = g.UniProtKB;")
  qmain2 <-dbGetQuery(con(), query2)
    
  query3<-paste0("SELECT 
    set1.Gene_1,
    set1.Interaction_ID,
    g.Gene_Symbol AS Gene_2
FROM 
    set1 
    INNER JOIN GeneInfo AS g ON set1.GeneInfoIDto = g.UniProtKB;")
  qmain3 <-dbGetQuery(con(), query3)
  e$data <- as.data.frame(qmain3)
  
  })
    })
  }
  
})

observeEvent(input$runVisPPI,{
  if (input$virusID >0 && input$virusID != "10" && input$virusID != "0"){
    isolate({
      dbGetQuery(con(), "USE CATE;")
      dbGetQuery(con(), "SET SESSION max_statement_time = 900;")
      
      query1 <- paste0("
      SELECT 
    g1.GeneInfoIDfrom AS GeneInfoIDfrom,
    gi1.Gene_Symbol AS Gene_1,
    gi1.Protein_Name AS Protein_Des_1,
    gi1.Taxonomy_ID AS Tax_1,
    t1.Organism_Name AS Organism_Name_1,
    g1.Interaction_ID AS Interaction_ID, 
    g1.GeneInfoIDto AS GeneInfoIDto,
    gi2.Gene_Symbol AS Gene_2,
    gi2.Protein_Name AS Protein_Des_2,
    gi2.Taxonomy_ID AS Tax_2,
    t2.Organism_Name AS Organism_Name_2
    FROM Gene2Gene g1
    INNER JOIN GeneInfo gi1 ON g1.GeneInfoIDfrom = gi1.UniProtKB
    INNER JOIN GeneInfo gi2 ON g1.GeneInfoIDto = gi2.UniProtKB
    INNER JOIN Taxonomy t1 ON gi1.Taxonomy_ID = t1.ID
    INNER JOIN Taxonomy t2 ON gi2.Taxonomy_ID = t2.ID
    WHERE (gi1.Taxonomy_ID IN (",input$virusID,")
           AND gi2.Taxonomy_ID = 9606)
    OR (gi1.Taxonomy_ID = 9606
        AND gi2.Taxonomy_ID IN (",input$virusID,"));
                       ")
      
      qmain1 <-dbGetQuery(con(), query1)
      
      e$data <- as.data.frame(qmain1)
      
    })
  }
  
  ################################################################################################################  
  if (input$virusID == "1"){
    isolate({
      edges <- read.csv("data/GLV-1h68.csv")
      e$data <- as.data.frame(edges)
      
      
    })
    
  }
  
  ################################################################################################################  
  ################################################################################################################  
  if (input$virusID == "10"){
    isolate({
      edges <- read.csv("data/VVTKN1L.csv")
      e$data <- as.data.frame(edges)
      
      
    })
    
  }
  
  ################################################################################################################  
})
##################################################################################################################

observeEvent(input$runVisPPI2,{
  if (input$virusID2 >0 && input$virusID2 != "10" && input$virusID2 != "0"){
    isolate({
      dbGetQuery(con(), "USE CATE;")
      dbGetQuery(con(), "SET SESSION max_statement_time = 900;")
      
      query1 <- paste0("
      SELECT 
    g1.GeneInfoIDfrom AS GeneInfoIDfrom,
    gi1.Gene_Symbol AS Gene_1,
    gi1.Protein_Name AS Protein_Des_1,
    gi1.Taxonomy_ID AS Tax_1,
    t1.Organism_Name AS Organism_Name_1,
    g1.Interaction_ID AS Interaction_ID, 
    g1.GeneInfoIDto AS GeneInfoIDto,
    gi2.Gene_Symbol AS Gene_2,
    gi2.Protein_Name AS Protein_Des_2,
    gi2.Taxonomy_ID AS Tax_2,
    t2.Organism_Name AS Organism_Name_2
    FROM Gene2Gene g1
    INNER JOIN GeneInfo gi1 ON g1.GeneInfoIDfrom = gi1.UniProtKB
    INNER JOIN GeneInfo gi2 ON g1.GeneInfoIDto = gi2.UniProtKB
    INNER JOIN Taxonomy t1 ON gi1.Taxonomy_ID = t1.ID
    INNER JOIN Taxonomy t2 ON gi2.Taxonomy_ID = t2.ID
    WHERE (gi1.Taxonomy_ID IN (",input$virusID2,")
           AND gi2.Taxonomy_ID = 9606)
    OR (gi1.Taxonomy_ID = 9606
        AND gi2.Taxonomy_ID IN (",input$virusID2,"));
                       ")
      
      qmain1 <-dbGetQuery(con(), query1)
      
      e$data <- as.data.frame(qmain1)
      
    })
  }
  
  ################################################################################################################  
  if (input$virusID2 == "1"){
    isolate({
      edges <- read.csv("data/GLV-1h68.csv")
      e$data <- as.data.frame(edges)
      
      
    })
    
  }
  
  ################################################################################################################  
  ################################################################################################################  
  if (input$virusID2 == "10"){
    isolate({
      edges <- read.csv("data/VVTKN1L.csv")
      e$data <- as.data.frame(edges)
      
      
    })
    
  }
  
  ################################################################################################################  
})

