library(shiny)
library(qtl)
library(qtlcharts)
library(intermediate)
source("miniDOQTL.R")

load("shiny_annotation.RData")

# for given MGI symbol, find Ensembl ids
other.ids <- function(gene.name, level) {
  if (level == "mRNA") {
    sel <- which(mRNA.list$symbol == gene.name)[1]
    if (!is.na(sel)) return(mRNA.list[sel,]) else return(c(NA,NA,NA))
  } 
  if (level == "protein") {
    sel <- which(protein.list$symbol == gene.name)[1]
    if (!is.na(sel)) return(protein.list[sel,]) else return(c(NA,NA,NA))
  }
}

# no data plot
nodata <- function() {
  plot(0, type="n", xaxt='n', yaxt='n', xlab="", ylab="")
  text(1, 0, "No data", cex=5)
}


# for given gene, level and plotType, return file with full path
file.path <- function(gene.name, level, plotType) {

  gene.triplet <- other.ids(gene.name, level)

  if (level == "mRNA") {
    output <- switch(as.character(plotType),
      "0" = paste0("scanone_mrna/",gene.triplet[1],"_",gene.triplet[2]),
      "1" = paste0("scanone_mrna_ciscov/",gene.triplet[1],"_",gene.triplet[2]),
      "2" = paste0("scanone_mrna_prot/",gene.triplet[3],"_",gene.triplet[1],"_",gene.triplet[2]),
      "3" = paste0("mediation_mrna_mrna/",gene.triplet[1],"_",gene.triplet[2]),
      "4" = paste0("mediation_mrna_prot/",gene.triplet[1],"_",gene.triplet[2]),
      "9" = paste0("scanint_mrna/",gene.triplet[1],"_",gene.triplet[2]),
      paste0("data_dump_mrna/",gene.triplet[1],"_",gene.triplet[2])
    )
  }
  if (level == "protein") {
    output <- switch(as.character(plotType),
      "0" = paste0("scanone_protein/",gene.triplet[3],"_",gene.triplet[1],"_",gene.triplet[2]),
      "1" = paste0("scanone_protein_ciscov/",gene.triplet[3],"_",gene.triplet[1],"_",gene.triplet[2]),
      "2" = paste0("scanone_prot_mrna/",gene.triplet[3],"_",gene.triplet[1],"_",gene.triplet[2]),
      "3" = paste0("mediation_prot_mrna/",gene.triplet[3],"_",gene.triplet[2]),
      "4" = paste0("mediation_prot_prot/",gene.triplet[3],"_",gene.triplet[2]),
      "9" = paste0("scanint_protein/",gene.triplet[3],"_",gene.triplet[2]),
      paste0("data_dump_prot/",gene.triplet[3],"_",gene.triplet[2])
   )
  }
 
  return(paste0("/QTL_mapping/kidney_combined/",output,".rds"))
}


# DOQTL type LOD plot and effect plot
create.plot <- function(gene.name, level, plotType, chr) {
  # if not ready, plot nothing
  if (is.null(gene.name) | is.null(level) | is.null(plotType) | is.null(chr)) return(NULL)

  file.location <- file.path(gene.name, level, plotType)
  if (!file.exists(file.location)) {nodata(); return(NULL)} 
  fit <- readRDS(file.location)

  title <- gene.name
  if (level=="protein") title <- toupper(title)

  # lod plots
  if (plotType == 0 | plotType == 1 | plotType == 2) {

    # sign. thresholds based on 10000 permutations (100 permutations x 100 randomly selected gnes)
    if (level=="protein") {
      sig.thr <- structure(c(6.33, 8.16, 9.22,
                           9.22, 11.45, 12.58), .Dim = c(3L,
                           2L), .Dimnames = list(c("0.63", "0.05", "0.01"), c("A", "X")))
    } else {
      sig.thr <- structure(c(6.28, 8.13, 9.14,
                           9.08, 11.32, 12.57), .Dim = c(3L,
                           2L), .Dimnames = list(c("0.63", "0.05", "0.01"), c("A", "X")))
    }

    if (chr == "all") {
      plot.doqtl(fit, sig.thr = sig.thr, main=title, sig.col = c("yellow", "blue", "red"))
    } else {
      if (chr!="X") {
        coefplot(fit, chr, main=title)
      } else {
        coefplot(fit, chr, main=paste(title, "females"), sex="F")
      }
    }
  }

  # mediation plots
  if (plotType == 3 | plotType == 4) {
    if (chr == "all") {
      med <- fit[fit$chr %in% c(as.character(1:19),"X"),] # to avoid problems with Y, M and NA
    } else {
      med <- fit[fit$chr == chr,]
    }

    plot(med, main=title)
  }
  
  # profile plots
  if (plotType == 5 | plotType == 6 | plotType == 7 | plotType == 8) {
    library(ggplot2)
    fit$Age <- factor(fit$Age, levels=c("6", "12", "18"))
    fit$Generation <- factor(fit$Generation, levels=c("G8","G9","G10","G11","G12"))

    # profile - Age
    if (plotType == 5) {
      p1 <- ggplot(aes(y=y, x=Age, colour=Sex, group=Sex), data=fit) +
          geom_jitter(position = position_jitter(width = .15)) +
          stat_summary(fun.y=mean, geom="line", size=1) +
          stat_summary(fun.y=mean, geom="point", size=2) +
          ylab("rankZ - expression")
    } 
    # profile - Sex
    if (plotType == 6) {
      p1 <- ggplot(aes(y=y, x=Sex, colour=Age, group=Age), data=fit) +
          geom_jitter(position = position_jitter(width = .15)) +
          stat_summary(fun.y=mean, geom="line", size=1) +
          stat_summary(fun.y=mean, geom="point", size=2) +
          ylab("rankZ - expression")
    }
    # profile - Generation
    if (plotType == 7) {
      p1 <- ggplot(aes(y=y, x=Generation,colour=Sex,shape=Age,group=Generation),data=fit)+
          geom_boxplot(outlier.shape = NA,show_guide=FALSE) + # remove outliers  
          geom_jitter(position = position_jitter(width = .15)) +
          ylab("rankZ - expression")
    }
    if (plotType == 8) {
      file.location2 <- file.path(gene.name, ifelse(level=="protein","mRNA","protein"), plotType)
      if (!file.exists(file.location2)) {nodata(); return(NULL)}
      fit2 <- readRDS(file.location2)
      fit$x <- fit2$y
 
      p1 <- ggplot(aes(y=y, x=x, colour=Sex, shape=Age, group=Sex),data=fit) + 
              geom_point() +
              geom_smooth(method="lm", fill=NA) +
              xlab(ifelse(level=="protein","mRNA","protein")) +
              ylab(ifelse(level=="protein","protein","mRNA")) 
    }
    print(p1 + ggtitle(title) + theme_bw())
  }
  if (plotType == 9) {
    if (chr == "all") {
      plot(fit, main=title, incl.markers=FALSE, bandcol="grey75", lodcolumn=c(1,3,2),
         col=c("green","yellow","black")) 
    } else {
      library(ggplot2)
      fit.chr <- fit[fit$chr==chr,]
      max.marker <- rownames(fit.chr)[which.max(fit.chr$lod2)[1]]
      if (level=="mRNA") ens <- other.ids(gene.name, level)[[1]] else ens <- other.ids(gene.name, level)[[3]]
      
      # if not yet loaded, load the data
      if (!exists("expr.mrna")) load("/QTL_mapping/kidney_combined/DO192_kidney.RData")
      
      if (level=="mRNA") y <- expr.mrna[,ens] else y <- expr.protein[,ens]
      y[pheno$Sex=="M"] <- y[pheno$Sex=="M"] - mean(y[pheno$Sex=="M"], na.rm=TRUE)
      y[pheno$Sex=="F"] <- y[pheno$Sex=="F"] - mean(y[pheno$Sex=="F"], na.rm=TRUE)
      coef <- se <- tstat <-  NULL
      for (age in c("6", "12", "18")) {
        sel <- pheno$Age == age
        lm.fit <- summary(lm(y[sel] ~ 0 + probs[sel,,max.marker]))
        coef <- cbind(coef, lm.fit$coef[1:8,1])
        se   <- cbind(se, lm.fit$coef[1:8,2])
        tstat   <- cbind(tstat, lm.fit$coef[1:8,3])
      }
      dt <- data.frame(Allele = LETTERS[rep(1:8,3)], 
                 Age = factor(rep(c(6,12,18), each=8)),
                 beta = as.vector(coef), 
                 beta.se = as.vector(se), 
                 Tstat=as.vector(tstat))
      mytheme <- theme(panel.border = element_rect(linetype = "solid", , fill=NA, colour = "grey"))                            
      p <- ggplot(aes(x=as.numeric(Allele)+as.numeric(Age)/15-2/15,
                      y=beta,colour=Age), data=dt) + 
                  geom_errorbar(aes(ymin=beta-beta.se, ymax=beta+beta.se), width=.1) + 
                  geom_point(aes(size=abs(Tstat))) + xlab("Allele") +
                  ylab("beta +/- SE") +
                  scale_x_continuous(breaks=c(1:8), labels=LETTERS[1:8]) + 
                  geom_abline(intercept = 0, slope = 0, colour=I("grey")) +
                  ggtitle(max.marker) + mytheme
       print(p)  
    }
  } 
}  

# interactive plots
create.iplot <- function(gene.name, level, plotType, chr) {
  # if not ready, plot nothing
  if (is.null(gene.name) | is.null(level) | is.null(plotType) | is.null(chr)) return(NULL)

  file.location <- file.path(gene.name, level, plotType)
  # if (!file.exists(file.location)) {nodata(); return(NULL)} # not yet implemented
  fit <- readRDS(file.location)

  title <- gene.name
  if (level=="protein") title <- toupper(title)

  # lod plots
  if (plotType == 0 | plotType == 1 | plotType == 2) {
    chr = c(fit$lod$A$chr, fit$lod$X$chr)
    pos = c(fit$lod$A$bp, fit$lod$X$bp) 
    x = intermediate::gmb.coordinates(chr=chr,pos=pos)
    y = c(fit$lod$A$lod, fit$lod$X$lod)
    snp = c(fit$lod$A$marker, fit$lod$X$marker)    

    p0 <- iplot(x=x, y=y, group=factor(chr), indID=snp, chartOpts=list(title=title))
    return(p0)
  }

  # mediation plots
  if (plotType == 3 | plotType == 4) {
    med <- fit[fit$chr %in% c(as.character(1:19),"X"),] # to avoid problems with Y, M and NA
    return(kplot(med, "symbol"))
  }

  # profile plots
  if (plotType == 5 | plotType == 6 | plotType == 7 | plotType == 8) {
    fit$Age <- as.numeric(as.character(fit$Age))
    fit$Generation <- as.numeric(factor(fit$Generation))
    fit$sex <- as.numeric(factor(fit$Sex))
    n <- nrow(fit)

    # profile - Age
    if (plotType == 5) {
      p1 <- with(fit, iplot(x=Age+runif(n, min=-0.5, max=0.5), y=y, group=Sex, 
                            indID=paste0(Sex,Sample.Number), chartOpts=list(title=title)))
      return(p1)
    }
    # profile - Sex
    if (plotType == 6) {
      p2 <- with(fit, iplot(x=sex+runif(n, min=-0.1, max=0.1), y=y, group=factor(Age),
                            indID=paste0(Sex,Sample.Number), chartOpts=list(title=title)))
      return(p2)
    }
    # profile - Generation
    if (plotType == 7) {
      p3 <- with(fit, iplot(x=Generation+runif(n, min=-0.1, max=0.1), y=y, group=Sex,
                            indID=paste0(Sex,Sample.Number), chartOpts=list(title=title)))
      return(p3)
    }
    if (plotType == 8) {
      file.location2 <- file.path(gene.name, ifelse(level=="protein","mRNA","protein"), plotType)
      if (!file.exists(file.location2)) return(NULL)
      fit2 <- readRDS(file.location2)
      fit$x <- fit2$y

      p4 <- with(fit, iplot(x=x, y=y, group=Sex,
                            indID=paste0(Sex,Sample.Number), chartOpts=list(title=title)))
      return(p4)
    }
  }
}



# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output, session) {

  session$allowReconnect(TRUE)
 
  # MGI symbol input, default from URL address if given
  output$selectTextInput <- renderUI({
    default <- parseQueryString(session$clientData$url_search)$gene
    if (is.null(default)) default <- "Ace"
    textInput("geneId", "MGI Symbol", value=default)
  })

  output$selectLevel <- renderUI({
    default <- parseQueryString(session$clientData$url_search)$level
    if (is.null(default)) default <- "mRNA"
    radioButtons("level", "Level", c("mRNA", "protein"), inline=TRUE, selected=default)
  })

  output$selectType <- renderUI({
    default <- parseQueryString(session$clientData$url_search)$type
    if (is.null(default)) default <- 0
    radioButtons("plotType", "Scan",
                c("lod - default" = 0, "lod - conditioned on local genotype" = 1,
                  "lod - conditioned on mRNA/protein" = 2,
                  "mediation - all mRNAs" = 3, "mediation - all proteins" = 4,
                  "profile - Age" = 5, "profile - Sex" = 6, "profile - Generation" = 7,
                  "comparison to protein/mRNA" = 8, "Experimental: Age*QTL interaction" = 9), selected=default)
  })

  output$selectChr <- renderUI({
    default <- parseQueryString(session$clientData$url_search)$chr
    if (is.null(default)) default <- "all"
    selectInput("chr","Chromosome",c("all",as.character(1:19),"X"), default)
  })
 
  # plot
  output$lodPlot <- renderPlot({create.plot(input$geneId, input$level, input$plotType, input$chr)})

  # iplot
  output$interactivePlot <- iplot_render({create.iplot(input$geneId, input$level, input$plotType, input$chr)})

  # download RDS file
  output$downloadRDS <- downloadHandler(
    filename = function() {basename(file.path(input$geneId, input$level, input$plotType))},
    content = function(file) {
      file.location <- file.path(input$geneId, input$level, input$plotType)
      if (!file.exists(file.location)) return(NULL)
      tmp <- readRDS(file.location)
      if (input$plotType==8) {
        file.location2 <- file.path(input$geneId, ifelse(input$level=="protein","mRNA","protein"), input$plotType)
        if (!file.exists(file.location2)) return(NULL)
        tmp2 <- readRDS(file.location2)
        tmp$x <- tmp2$y
      }
      
      saveRDS(tmp, file)
    }
  )
  
# download RDS file
  output$downloadCSV <- downloadHandler(
    filename = function() {sub("rds$", "csv", basename(file.path(input$geneId, input$level, input$plotType)))},
    content = function(file) {
      file.location <- file.path(input$geneId, input$level, input$plotType)
      if (!file.exists(file.location)) return(NULL)
      tmp <- readRDS(file.location)
      if (input$plotType==8) {
        file.location2 <- file.path(input$geneId, ifelse(input$level=="protein","mRNA","protein"), input$plotType)
        if (!file.exists(file.location2)) return(NULL)
        tmp2 <- readRDS(file.location2)
        tmp$x <- tmp2$y
      }
      
      if (input$plotType==0 | input$plotType==1 | input$plotType==2){
        chr = c(tmp$lod$A$chr, tmp$lod$X$chr)
        pos = c(tmp$lod$A$bp, tmp$lod$X$bp)
        lod = c(tmp$lod$A$lod, tmp$lod$X$lod)
        marker = c(tmp$lod$A$marker, tmp$lod$X$marker)
        tmp <- data.frame(marker=marker, chr=chr, pos=pos, lod=lod)
      }

      write.csv(tmp, file, row.names=FALSE)
    }
  )


  # gene description
  output$description <- renderUI({ 
      # if not ready, no description
      if (is.null(input$geneId) | is.null(input$level) | is.null(input$plotType) | is.null(input$chr)) return(NULL)

      tmp <- other.ids(input$geneId, input$level)
      ensembl_link <- paste0("http://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=", tmp[1])
      mgi_link <- paste0("http://www.informatics.jax.org/searchtool/Search.do?query=", tmp[2])
      credible_link <- function() { #calculate CI and returns the link to Ensembl
        if (input$chr=="all" | as.numeric(input$plotType)>2) return(NULL)
        rds <- readRDS(file.path(input$geneId, input$level, input$plotType))
        tmp <- bayesint(rds, input$chr, expandtomarkers = TRUE)
        CI <- round(tmp[c(1,3),3]*10^6)
        link <- paste0("http://useast.ensembl.org/Mus_musculus/Location/View?r=",input$chr,
                       "%3A",CI[1],"-",CI[2])
        p("Bayesian 95% credible interval for QTL on this chromosome is ", 
          a(paste0("chr",input$chr,": ", CI[1],"-", CI[2]), href=link, target="_blank"))
      }
      HTML(ifelse(is.na(tmp[1]), "",paste(
                  br(), br(),
                  p(strong(tmp[2])," gene is located on chromosome ", tmp[4],
                  "(", round(tmp[5]/10^6,2), " - ", round(tmp[6]/10^6,2), "Mb,", 
                  ifelse(tmp[7]==1,"positive","negative")," strand), ",
                  a("[Ensembl]", href=ensembl_link, target="_blank"),
                  a("[MGI]", href=mgi_link, target="_blank")),
                  credible_link()
          )))})
})
