
version 1.0

workflow run_scatd_LD_aou {
    input {
        Array[String] ancestry_array = ["afr", "amr", "eur"]
        Array[Array[String]] blocks 
    }
    scatter (ancestry in ancestry_array) {
        scatter (block in blocks) {
            call run_ld_aou { 
                input:
                    ancestry = ancestry,
                    chr = block[0],
                    start = block[1],
                    end = block[2]
            }
        }
    }
}


task run_ld_aou {

    input {
        String ancestry
        String chr
        String start
        String end
    }

    command <<<
        set -euox pipefail
        
        ### download block files, ancestry IID lists, and srWGS data
        gcloud config set storage/parallel_composite_upload_enabled False
        gcloud storage cp aou_bucket/LD/*~{ancestry}*fidiid.txt \
                          aou_bucket/wgs_v8_pgens/chr~{chr}_trimmed_merged.p* .


        ### separate chr bed files into 20Mb blocks with 1Mb overlaps using block files
        plink2 \
            --pfile chr~{chr}_trimmed_merged \
            --keep keep_unrelated_~{ancestry}_fidiid.txt \
            --chr ~{chr} --from-bp ~{start} --to-bp ~{end} \
            --mac 5 \
            --make-bed \
            --out ~{ancestry}_forgcta
        
        
        ### Calculating LD in 1MB windows with 0.5MB overlap for each ancestry group
        # in parallel per chromosome
        # compress LD output files while going

        gcta64 \
            --threads 32 \
            --bfile ~{ancestry}_forgcta \
            --ld-score \
            --ld-wind 1000 \
            --ld-rsq-cutoff 0.01 \
            --out ~{ancestry}_~{chr}_~{start}_~{end}_ld 
    

        ### Upload to cloud
        gcloud storage cp *score.ld *.log "aou_bucket/LD/"
    >>>


    runtime {
        docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/tvlupu/plink_plink2_gcta:19.20.64" 
        memory: "32 GiB"
        cpu: 8
        disks: "local-disk 1000 HDD"
        zones: "us-central1-a"
        cpuPlatform: "Intel Ice Lake"
        preemptible: 3
    }
}
