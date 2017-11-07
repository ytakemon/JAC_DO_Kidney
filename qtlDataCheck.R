# From: https://github.com/churchill-lab/qtl-viewer/blob/master/scripts/qtlDataCheck.R
# VERY rough
#

IsVariableOK <- function(varName, varClass, varRequired) {
    # Check if varable is good
    #
    #  Example: exists('genome_build')
    #           class(genome_build) == 'character'
    #
    if ((varRequired) && (!exists(varName))) {
        print(paste0('ERROR: ', varName, ' does not exist, but should'))
        return (FALSE)
    } else {
        classesFound <- class(get(varName))
        if (!any(varClass == classesFound)) {
            print(paste0('ERROR: ', varName, ' is type ', classesFound, ', not ', varClass))
            return (FALSE)
        }
    }
    return (TRUE)
}

CheckVariables <- function() {
    # Check to see if the variables exist and the of the correct type

    # grab the datasets in the environment
    datasets <- grep('^dataset*', apropos('dataset\\.'), value=TRUE)

    if (length(datasets) == 0) {
        print('ERROR: No datasets found')
    }

    # expected
    allNames <- c('genome.build', 'genoprobs', 'K', 'map', 'markers', datasets)
    allClasses <- c('character', 'calc_genoprob', 'list', 'list', 'data.frame', rep('list', length(datasets)))
    allRequired <- c(rep(TRUE, length(allNames)))

    # construct the data.frame
    dataCheck <- data.frame(name=allNames, class=allClasses, required=allRequired, stringsAsFactors=FALSE)

    # check the variables
    errors <- mapply(IsVariableOK, dataCheck$name, dataCheck$class, dataCheck$required)
}

CheckDatasets <- function(all_vars) {
    # Check to see if the names in each dataset exist and the of the correct class

    # grab the datasets in the environment
    datasets <- grep('^dataset*', apropos('dataset\\.'), value=TRUE)

    # expected elements
    phenoNames   <- c('annots',     'covar',  'covar.factors', 'datatype',  'display.name', 'lod.peaks',  'pheno',  'samples')
    phenoClasses <- c('data.frame', 'matrix', 'data.frame',    'character', 'character',    'data.frame', 'matrix', 'data.frame')
    phenoRequired <- c(TRUE,         TRUE,     TRUE,            TRUE,        FALSE,          TRUE,         TRUE,     TRUE)

    mrnaNames   <- c('annots',     'covar',      'covar.factors', 'datatype',  'display.name', 'ensembl.version', 'expr',   'lod.peaks',  'raw',    'samples')
    mrnaClasses <- c('data.frame', 'data.frame', 'data.frame',    'character', 'character',    'numeric',         'matrix', 'data.frame', 'matrix', 'data.frame')
    mrnaRequired <- c(TRUE,         TRUE,         TRUE,            TRUE,        FALSE,          FALSE,             TRUE,     TRUE,         FALSE,    TRUE)

    # construct the data.frame
    dataCheck <- data.frame(name=mrnaNames, class=mrnaClasses, required=mrnaRequired, stringsAsFactors=FALSE)

    # be explicit to the end user, we could have used mapply and some trickery, but this is simple
    for (d in datasets) {
        dataset <- get(d)

        if (!('datatype' %in% names(dataset))) {
            print(paste0('ERROR: data_type missing for dataset:', d))
            print(paste0('ERROR: Skipping rest of check for dataset:', d))
            next
        }

        dataType = dataset[['datatype']]

        if (!(dataType %in% c('protein', 'mRNA', 'phenotype'))) {
            print(paste0('ERROR: data_type of ', dataType, ' is invalid in dataset: ', d))
            print('ERROR: data_type should be mRNA, protein, or pheno')
            print(paste0('ERROR: Skipping rest of check for dataset: ', d))
            next
        }

        nameList <- phenoNames
        classList <- phenoClasses
        requiredList <- phenoRequired

        if (dataType %in% c('protein', 'mRNA')) {
            nameList <- mrnaNames
            classList <- mrnaClasses
            requiredList <- mrnaRequired
        }

        # look at the names and classes in the list
        for (n in c(1:length(nameList))) {
            varName <- nameList[n]
            className <- classList[n]
            required <- requiredList[n]

            if ((required) && (!(varName %in% names(dataset)))) {
                print(paste0('ERROR: ', varName, ' does not exist in dataset: ', d))
            } else {
                classesFound <- class(dataset[[varName]])
                if (!any(className == classesFound)) {
                    print(paste0('ERROR: ', varName, ' is type: ', classesFound, ', should be type: ', className, ', in dataset: ', d))
                }
            }
        }

        CheckDataNames(ds = dataset)
    }
}


# Verify that the sample IDs and marker IDs match between all of the objects in a dataset.
# This assumes that the dataset has already been checked to contain 'datatype'.
# Argument: ds: list that is a dataset.
CheckDataNames <- function(ds) {

    if(ds$datatype == "phenotype") {

        # Phenotype names in phenotype annotation.
        pheno.union = union(ds$annots$R_name, colnames(ds$pheno))
        pheno.inter = intersect(ds$annots$R_name, colnames(ds$pheno))
        if(length(pheno.union) != length(pheno.inter)) {
            wh <- setdiff(pheno.union, pheno.inter)
            print("ERROR: The following phenotypes do not match between pheno and annots...")
            print(paste(wh, collapse = ", "))
        }

        # Sample IDs match everywhere.
        if(!all(rownames(ds$pheno) == rownames(ds$covar))) {
            print("ERROR: Samples IDs do not match between pheno and covar.")
        }

        if(!all(rownames(ds$pheno) %in% rownames(genoprobs[[1]]))) {
            print("ERROR: Samples IDs do not match between pheno and genoprobs.")
        }

    } else if(ds[["datatype"]] == "mRNA" | ds[["datatype"]] == "protein"){

        # TODO: Fill this in.

    } else {
        stop("ERROR: Bad datatype for dataset.")
    }
}


CheckExtraVars <- function(allNames) {

    # Check for extra variables
    expectedNames <- c('genome.build', 'genoprobs', 'K', 'map', 'markers',
                       grep('^dataset*', apropos('dataset\\.'), value=TRUE))

    extraNames <- setdiff(allNames, expectedNames)
    if (length(extraNames) > 0) {
        print('Warning: the following extra variables were found...')
        print(extraNames)
    }

    if(any(!expectedNames[1:5] %in% allNames)) {
        wh <- which(!(expectedNames[1:5] %in% allNames))
        print(paste0("ERROR: The followign required variables not found..."))
        print(paste(expectedNames[wh], collapse = ", "))
    }

    # At this point, genoprobs, K, map and markers are in the environment.

    # Check Sample IDs for genoprobs and K.
    if(length(genoprobs) != length(K)) {
        print(paste0("ERROR: genoprobs (", length(genoprobs), ") and K (", length(K), ") do not have the same length."))
    } else {
        if(any(names(genoprobs) != names(K))) {
            print("ERROR: names of genoprobs and K do not match.")
            print(paste("names(genoprobs) =", paste(names(genoprobs), collpase = ", ")))
            print(paste("names(K) =", paste(names(K), collpase = ", ")))
        }

        rownames.eq <- mapply(function(x, y) { all(rownames(x) == rownames(y)) }, genoprobs, K)
        if(any(rownames.eq == FALSE)) {
            print("ERROR: sample IDs do not match between genoprobs and K.")
        }
    }

    # Check Marker IDs for genoprobs and map.
    if(length(genoprobs) != length(map)) {
        print(paste0("ERROR: genoprobs (", length(genoprobs), ") and map (", length(map), ") do not have the same length."))
    } else {
        rownames.eq <- mapply(function(x, y) { all(dimnames(x)[[3]] == names(y)) }, genoprobs, map)
        if(any(rownames.eq == FALSE)) {
            print("ERROR: marker names do not match between genoprobs and map.")
        }
    }

    # Check dimensions of markers and map.
    map.length = sum(sapply(map, length))
    if(map.length != nrow(markers)) {
        print(paste("ERROR: Number of rows in markers (", nrow(markers), ") does not equal the number of markers in map (", map.length, ")",))
    }
}
