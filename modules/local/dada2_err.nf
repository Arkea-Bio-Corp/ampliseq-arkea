process DADA2_ERR {
    tag "$meta.run"
    label 'process_medium'

    conda "bioconda::bioconductor-dada2=1.34.0 conda-forge::r-base=4.4.3 conda-forge::tbb=2020.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-dada2:1.34.0--r44he5774e6_2' :
        'biocontainers/bioconductor-dada2:1.34.0--r44he5774e6_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.err.rds"), emit: errormodel
    tuple val(meta), path("*.err.pdf"), emit: pdf
    tuple val(meta), path("*.err.svg"), emit: svg
    tuple val(meta), path("*.err.log"), emit: log
    tuple val(meta), path("*.err.convergence.txt"), emit: convergence
    path "versions.yml"               , emit: versions
    path "*.args.txt"                 , emit: args


    script:
    def prefix = task.ext.prefix ?: "prefix"
    def args = task.ext.args ?: ''
    def seed = task.ext.seed ?: '100'
    def binned = params.binned_quality ? 'TRUE': 'FALSE'
    if (!meta.single_end) {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))
        set.seed($seed) # Initialize random number generator for reproducibility

        fnFs <- sort(list.files(".", pattern = "_1.filt.fastq.gz", full.names = TRUE), method = "radix")
        fnRs <- sort(list.files(".", pattern = "_2.filt.fastq.gz", full.names = TRUE), method = "radix")

        if ($binned){
            # Binned quality score error model
            binnedQs <- c(2, 11, 25, 37)
            binnedQualErrfun <- makeBinnedQualErrfun(binnedQs)
            errF <- learnErrors(fnFs, errorEstimationFunction = binnedQualErrfun, multithread = $task.cpus, verbose = TRUE)
            errR <- learnErrors(fnRs, errorEstimationFunction = binnedQualErrfun, multithread = $task.cpus, verbose = TRUE)
        } else {
            # Standard DADA2 error model
            errF <- learnErrors(fnFs, $args, multithread = $task.cpus, verbose = TRUE)
            errR <- learnErrors(fnRs, $args, multithread = $task.cpus, verbose = TRUE)
        }

        saveRDS(errF, "${prefix}_1.err.rds")
        saveRDS(errR, "${prefix}_2.err.rds")
        sink(file = NULL)

        pdf("${prefix}_1.err.pdf")
        plotErrors(errF, nominalQ = TRUE)
        dev.off()
        svg("${prefix}_1.err.svg")
        plotErrors(errF, nominalQ = TRUE)
        dev.off()

        pdf("${prefix}_2.err.pdf")
        plotErrors(errR, nominalQ = TRUE)
        dev.off()
        svg("${prefix}_2.err.svg")
        plotErrors(errR, nominalQ = TRUE)
        dev.off()

        sink(file = "${prefix}_1.err.convergence.txt")
        dada2:::checkConvergence(errF)
        sink(file = NULL)

        sink(file = "${prefix}_2.err.convergence.txt")
        dada2:::checkConvergence(errR)
        sink(file = NULL)

        write.table('learnErrors\t$args', file = "learnErrors.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
        """
    } else {
        """
        #!/usr/bin/env Rscript
        suppressPackageStartupMessages(library(dada2))
        set.seed($seed) # Initialize random number generator for reproducibility

        fnFs <- sort(list.files(".", pattern = ".filt.fastq.gz", full.names = TRUE))

        sink(file = "${prefix}.err.log")

        if ($binned){
            binnedQs <- c(2,11,25,37)
            binnedQualErrfun <- makeBinnedQualErrfun(binnedQs)
            errF <- learnErrors(fnFs, errorEstimationFunction = binnedQualErrfun, multithread = $task.cpus, verbose = TRUE)
        } else {
            errF <- learnErrors(fnFs, $args, multithread = $task.cpus, verbose = TRUE)
        }
        saveRDS(errF, "${prefix}.err.rds")
        sink(file = NULL)

        pdf("${prefix}.err.pdf")
        plotErrors(errF, nominalQ = TRUE)
        dev.off()
        svg("${prefix}.err.svg")
        plotErrors(errF, nominalQ = TRUE)
        dev.off()

        sink(file = "${prefix}.err.convergence.txt")
        dada2:::checkConvergence(errF)
        sink(file = NULL)

        write.table('learnErrors\t$args', file = "learnErrors.args.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, na = '')
        writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dada2: ", packageVersion("dada2")) ), "versions.yml")
        """
    }
}
