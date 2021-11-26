=============================================================
Tutorial: Using OCMS_16S to process amplicon sequencing data
=============================================================

Here we show you how to run OCMS_16S dada2 to process 16S rRNA amplicon sequencing data. We are going to use publically available data that were published `here`_. The data were generated from stool samples of either wild type or MMP9 -/- mice that were either treated with dextran sodium sulphate to induce colitis or vehicle. Sequencing was performed on the Illumina MiSeq platform which produced 2 x 250bp reads targetting the V4 region of the 16S rRNA gene. This step-by-step tutorial should help you to get a sense of the pipeline and how you might use it to process your own data.


Step 1: Download tutorial data
-------------------------------

We have provided links to raw data in the OCMS_16S repository (OCMS_16S/tutorial). To download the data, create a working directory where you will run the pipeline and download the files::

    mkdir ocms_16s_tutorial
    cd ocms_16s_tutorial
    bash <path-to-OCMS_16S/tutorial/download_files.sh>

This will place the raw fastq files in your working directory.

Step 2: Rename files
-----------------------

To allow the pipeline to run, the filenames must satisfy a standard format. This is <something>.fastq.1.gz (with <something>.fastq.2.gz when data are paired). Therefore in this example we will rename the files::

    rename .R1.fastq.gz .fastq.1.gz *.R1.fastq.gz
    rename .R2.fastq.gz .fastq.2.gz *.R2.fastq.gz

Step 3: Create and edit configuration file
-------------------------------------------

We need a pipeline.yml file to pass parameters to pipeline tasks. To create this we type::


    ocms_16s dada2 config

This produces the pipeline.yml file that we can now edit::


    #############################################################################
    #############################################################################
    #############################################################################
    # Configuration file for OCMS_16S dada2.
    #
    # The parameters and options that are passed to the wrappers around dada2
    # functions are the same as those in the dada2 package. For a full list of
    # the parameters and verbatim explanation of them (from dada2 package)
    # you can use the help functions
    #
    # e.g. ocms_16s dada2 help --filterAndTrim
    #      ocms_16s dada2 help --sampleInference
    #      ocms_16s dada2 help --assignTaxonomy
    #
    #############################################################################
    #############################################################################
    #############################################################################

    # Edit parameters

    # specify whether data are paired or single end. The
    # pipeline will pick up whether this is true but being
    # explicit here is helpful
    paired: 1

    # dada2 parameters
    trim:

        # parameters used for trimming reads. If the data are
        # paired-end then you need to specify 2 values for
        # maxee, truncLen and trimLeft. These parameters must be specified.
        maxn: 0
        maxee: 2,2
        truncq: 2
        trunclen: 250,250
        trimleft: 0,0

    sample_inference:

        memory: 10G

        # parameters for sample inference. This includes
        # error learning, de-replication, merging (if paired) and
        # sample inference.

        # number of reads to use (per sample) to estimate error
        # model
        nbases: 10000000

        # additional options
        options: ''

    taxonomy:

        memory: 10G

        # assigning taxonomy
        taxonomy_file: RefSeq-RDP16S_v2_May2018.fa.gz

        # This is the file that is used for the addSpecies function in
        # dada2 for exact matching and species assignment. It must therefore
        # be derived from the same database used as taxonomy_file above
        species_file: RefSeq-RDP_dada2_assignment_species.fa.gz

    report:
        #   whether to run diagnostics report. This is only necessary if after the
        # main report is built you want to get into more regarding the specifics of
        # how dada2 processed sequences. Specify as 1 if you wish to run it
        diagnostics:

    database:
        # name of the output database. This is a database that is built to
        # be compatible with the OCMSlooksy.
        name: output_db


For this example, we specify that the data are paired and we want the final length of both the forward and reverse reads to be 250bp. This will result in fully overlapping reads when pairs are merged during dada2 processing.

The default settings are taken from the defaults used by dada2. If you want an explanation of the parameters for the dada2 steps then you can type for example::

    ocms_16s help --filterAndTrim

This will show you the options that are passed from the pipeline.yml to the R scripts along with a description::


    -packages/ocms_16S-0.0.1-py3.8.egg/ocms16S/R/dada2_filter_and_trim.R [options]


    Options:
            -i INFILE, --infile=INFILE
                    input fastq file [default NA]

            -f FILTERED-DIRECTORY, --filtered-directory=FILTERED-DIRECTORY
                    directory for filtered fastq files [default filtered]

            -p, --paired
                    is it paired-end data [default FALSE]

            -n MAXN, --maxN=MAXN
                    maxN parameter [default 0]

            -e MAXEE, --maxEE=MAXEE
                    maximum number of expected errors [default 2,2]

            --truncQ=TRUNCQ
                    truncate reads at the first instance of a quality score less
                    than or equal to truncQ [default 2]

            --truncLen=TRUNCLEN
                   truncate  reads  after truncLen bases [default 250,250]

            --trimLeft=TRIMLEFT
                   trim left sequence (primers) [default 0,0]

            -h, --help
                   Show this help message and exit


We will leave the majority of settings as they are for this example. However, we need to specify annotation files that will be used to assign taxonomic information to the amplicon sequence variants (ASVs) that are produced by dada2. In this example we will download them into our working directory, however you may want to have them somewhere else for using with future data. E.g. to use RefSeq databases do::

        wget https://zenodo.org/record/2541239/files/RefSeq-RDP16S_v2_May2018.fa.gz
        wget https://zenodo.org/record/2658728/files/RefSeq-RDP_dada2_assignment_species.fa.gz


These are then specified in the pipeline.yml as above.

Step 4: Run the pipeline
--------------------------

Once you are happy with the parameterisation, you can run the pipeline::

    ocms_16s dada2 make full -v5 -p40

Here we are running the pipeline using 80 processors as this is the number of samples we have - they will be processed in parallel on the cluster. If you are running this on a laptop make sure to specify the --local flag::

    ocms_16s dada2 make full -v5 -p1 --local

The -v5 specifies the verbosity level of the logging information. At 5 this will be very verbose and useful for debugging. you can check how the pipeline is progressing by viewing the pipeline.log file that is created in the working directory::

    cat pipeline.log


When the pipeline has finished running, the log file will look like this::

    tail pipeline.log

    ...

    2020-02-06 12:08:28,597 INFO main task - Task enters queue = 'full' 
    2020-02-06 12:08:28,598 INFO main task -     Job no function to check if up-to-date 
    2020-02-06 12:08:28,799 INFO main task -     Job completed
    2020-02-06 12:08:28,799 INFO main control - {"task": "'full'", "task_status": "completed", "task_total": 0, "task_completed": 0, "task_completed_percent": 0}
    2020-02-06 12:08:28,799 INFO main task - Completed Task = 'full' 
    2020-02-06 12:08:28,801 INFO main experiment - job finished in 3002 seconds at Thu Feb  6 12:08:28 2020 -- 66.94 45.21  5.67  9.90 -- 8e55fc3a-59c9-4f94-b014-1062ee84c9bc


Step 5: Build the report
-------------------------

You can then build the report to inspect the performance of the processing::


    ocms_16s dada2 make build_report


This will produce the file report.dir/report.html that you can view in your browser.
 

Step 6: Downstream analysis with OCMSlooksy
--------------------------------------------

OCMSlooksy is an R/Shiny application that will take the output from OCMS_16S dada2 and enable the user to inspect the parameters of the dada2 run, the processing results as well as perform downstream viosualisation and statistical analyses. To format the outputs from OCMS_16S dada2 for OCMSlooksy you can run::


    ocms_16s dada2 make build_db

In the example above this will create the SQLite database, 'output_db', that serves as the main input to OCMSlooksy.


.. _here: https://www.nature.com/articles/s41522-018-0059-0
