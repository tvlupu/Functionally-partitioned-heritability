version 1.0

workflow run_scat_PCs_aou {
  input {
    Array[String] ancestry_array
  }
  scatter (ancestry in ancestry_array) {
    call run_PCs_aou { 
        input:
            ancestry = ancestry
        }
    }
}


task run_PCs_aou {

    input {
        String ancestry
    }

    command <<<
    set -euox pipefail

    # get R files
    gsutil -m cp aou_bucket/PCs/PCs.r  \
                aou_bucket/PCs/backingfiles/impute.r .


    ## PART 1 IMPUTATION ##
        # get QCed files
        gsutil -m cp aou_bucket/PCs/~{ancestry}QC_umap1* . 

        #plink2 --bfile ~{ancestry}QC --keep keep_umap1_~{ancestry}_fidiid.txt --make-bed --out ~{ancestry}QC_umap1
        #^ was done locally

        # convert to bigsnpR and impute
        Rscript impute.R ~{ancestry}

        rm *.bed *.bim *.fam # to save memory
         
        # upload backingfiles
        gsutil -m cp ~{ancestry}*bk ~{ancestry}*rds "aou_bucket/backingfiles/"
                      
    
    ## PART 2 PCS ##
        # compute 120 PCs
        Rscript PCs.R ~{ancestry}

        # upload PCs
        gsutil -m cp *.csv "aou_bucket/PCs/"

    >>>


    runtime {
        docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/tvlupu/bigsnpr_pca_computation:413" 
        memory: "128 GiB"
        cpu: 32
        disks: "local-disk 5000 HDD"
        zones: "us-central1-a"
        cpuPlatform: "Intel Ice Lake"
    }
}
