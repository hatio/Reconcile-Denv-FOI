
    #   Theoretical analysis of seroincidence rates
    #   ...........................................
    
Rscript --vanilla -e "
    rmarkdown::render(
        'Scripts/theory/seroincidence.R'
        , output_format = 'pdf_document'
        , output_dir = 'Report/theory'
        , output_file = 'seroincidence.pdf'
        , knit_root_dir = getwd()
    )"

