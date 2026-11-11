
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
        String umap_no
    }

    command <<<
    set -euox pipefail

    # get R file + QCed files
    gsutil -m cp aou_bucket/PCs/impute_runPCs.R \
                 aou_bucket/PCs/~{ancestry}_umap~{umap_no}* . 

# QC done locally with:
#plink2 --bfile ~{ancestry}QC --keep keep_umap~{umap_no}_~{ancestry}_fidiid.txt --autosome --make-bed --out ~{ancestry}_umap~{umap_no}

        # convert to bigsnpR + impute + calculate PCs
        Rscript impute_runPCs.R ~{ancestry} ~{umap_no}
        rm *.bed *.bim *.fam # to save memory
         
        # upload backingfiles
        gsutil -m cp  ~{ancestry}*bk ~{ancestry}*rds "aou_bucket/backingfiles/"
        
        # upload PCs
        gsutil -m cp *.csv "aou_bucket/PCs/"

    >>>


    runtime {
        docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/tvlupu/bigsnpr_pca_computation:413" 
        memory: "128 GiB"
        cpu: 32
        disks: "local-disk 1000 HDD"
        zones: "us-central1-a"
        cpuPlatform: "Intel Ice Lake"
    }
}
