library(dplyr)
library(ggplot2)
library(biomaRt)
setwd("/projects/korstanje-lab/ytakemon/JAC_DO_Kidney/")
load("./shiny_annotation.RData")
load("./RNAseq_data/DO1045_kidney.Rdata")

# coarse into one column per phenotype in new object: pheno
# DO-0021 is a bug in the cross-sectional that needs to be removed
Upheno <- Upheno %>%
          filter(study == "Cross-sectional" & Mouse.ID != "DO-0021") %>%
          mutate(
            cr_all_1plog = log1p(coalesce(cr.u.6, cr.u.12, cr.u.18)),
            alb_all_1plog = log1p(coalesce(ma.u.6, ma.u.12, ma.u.18)),
            cr_all_log = log(coalesce(cr.u.6, cr.u.12, cr.u.18)),
            phs_all_log = log(coalesce(phs.u.6, phs.u.12, phs.u.18))) %>%
          select(Mouse.ID, cr_all_1plog, alb_all_1plog, cr_all_log, phs_all_log)

# Define the following variables
pheno <- "alb"
gene_select <- "Akt1"
level_select <- "mRNA"
chr_select <- "12"

if (pheno == "alb"){
  file <- paste0("./QTLscan/addscan_urine/Intscan_alb_all.rds")
} else if (pheno == "phs"){
  file <- paste0("./QTLscan/addscan_urine/Intscan_phs_all.rds")
}

fit <- readRDS(file)
fit <- as.data.frame(fit)
fit$chr <- MM_snps[MM_snps$marker %in% rownames(fit),]$chr
fit.chr <- fit[fit$chr==chr_select,]
max.marker <- rownames(fit.chr)[which.max(fit.chr[,1])]
max.pos <- MM_snps[MM_snps$marker == max.marker,]$pos

Upheno$Sex <- samples %>% filter(Mouse.ID %in% Upheno$Mouse.ID) %>% select(Sex, Cohort.Age.mo)



fit <- readRDS("")
    # fit
    fit <- as.data.frame(fit)
    fit$chr <- substr(rownames(fit),1,regexpr("_",rownames(fit)) - 1)
    fit.chr <- fit[fit$chr==chr_select,]
    max.marker <- rownames(fit.chr)[which.max(fit.chr$pheno1)]
    if (level_select=="mRNA") ens <- other.ids(symbol, level_select)[[1]] else ens <- other.ids(symbol, level_select)[[3]]
    # Validate existance before moving forward
    validate(
      need(!is.na(ens), "Query is not found in RNAseq or Proteomics dataset.")
    )
    if (level_select=="mRNA") y <- expr.mrna[,ens] else y <- expr.protein[,ens]
    y[annot.samples$Sex=="M"] <- y[annot.samples$Sex=="M"] - mean(y[annot.samples$Sex=="M"], na.rm=TRUE)
    y[annot.samples$Sex=="F"] <- y[annot.samples$Sex=="F"] - mean(y[annot.samples$Sex=="F"], na.rm=TRUE)
    coef <- se <- tstat <-  NULL

    # Get effect by age
    for (age in c("6", "12", "18")) {
      sel <- annot.samples$Age == age
      lm.fit <- summary(lm(y[sel] ~ 0 + genoprobs[sel,,max.marker]))
      coef <- cbind(coef, lm.fit$coef[1:8,1])
      se   <- cbind(se, lm.fit$coef[1:8,2])
      tstat   <- cbind(tstat, lm.fit$coef[1:8,3])
    }

    # Compile into df
    dt <- data.frame(Allele = LETTERS[rep(1:8,3)],
               Age = factor(rep(c(6,12,18), each=8)),
               beta = as.vector(coef),
               beta.se = as.vector(se),
               Tstat=as.vector(tstat))

    ggplot(aes(x=as.numeric(Allele)+as.numeric(Age)/15-2/15, y=beta,colour=Age), data=dt) +
          geom_errorbar(aes(ymin=beta-beta.se, ymax=beta+beta.se), width=.1) +
          geom_point(aes(size=abs(Tstat))) + xlab("Allele") +
          ylab("beta +/- SE") +
          scale_x_continuous(breaks=c(1:8), labels=LETTERS[1:8]) +
          geom_abline(intercept = 0, slope = 0, colour=I("grey")) +
          ggtitle(paste(max.marker, gene_select)) +
          theme(panel.border = element_rect(linetype = "solid", fill=NA, colour = "grey"))
  }

  # Render plot ---------------------------------------
  output$plot <- renderPlot({
      if (is.null(AgeQTLIntPlot())){
        return(NULL)
      }
      AgeQTLIntPlot()
    })

  # Download plot --------------------------
  output$download_plot <- downloadHandler(
    filename <- function(){
      paste0(input$gene_input, "_chr", input$chr,"_ageQTLInteraction.pdf")
      },
    # content must be a function with arguemtn files to write plot
    content <- function(file) {
      pdf(file, width = 11, height = 7) #open device
        print(AgeQTLIntPlot()) #print plot
      dev.off() # close device
    }
  )

  # Output links --------------------------------
  output$gene_links <- renderText({

    # if not ready
    if(is.null(AgeQTLIntPlot())){
      return(NULL)
    }

    # find gene links
    gene <- input$gene_input
    if ( substr(gene, 1, 7) == "ENSMUSG"){
      mart_extract <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "ensembl_transcript_id",
      																"chromosome_name", "start_position", "end_position"),
                                      filters = "ensembl_gene_id",
                                      values = gene,
      																mart = ensembl)
    } else {
      mart_extract <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "ensembl_transcript_id",
      																"chromosome_name", "start_position", "end_position"),
                                      filters = "mgi_symbol",
                                      values = gene,
      																mart = ensembl)
    }
    symbol <- mart_extract$mgi_symbol[1]
    ens_id <- mart_extract$ensembl_gene_id[1]
    chr <- mart_extract$chromosome_name[1]
    start <- mart_extract$start_position[1]
    end <- mart_extract$end_position[1]
    ensembl_link <- paste0("http://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=", ens_id)
    mgi_link <- paste0("http://www.informatics.jax.org/searchtool/Search.do?query=", ens_id)

    # output links
    paste(p(symbol, "is located on chromosome", chr, ":", start, "-", end),
          a("[Ensembl]", href = ensembl_link, target="_blank"),
          a("[MGI]", href = mgi_link, target = "_blank"),br(),
          a("Diversity Outbred founder strain guide", href = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4602074/figure/Fig1/", target = "_blank"))
  })
}

# Run the app -----------------------------------------------------------------
shinyApp(ui = ui, server = server)
