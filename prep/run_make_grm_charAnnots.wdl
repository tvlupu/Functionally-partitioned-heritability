version 1.0

workflow run_make_grm_charAnnots_aou {
    input {
        String ancestry
        String filtflag
        Array[Array[String]] ldsc 
        Array[String] parts
    }

    String max_parts = length(parts)

    scatter (part in parts) {
        scatter (bin in ldsc) {
            call make_grm_bins { 
                input:
                    ancestry = ancestry,
                    filtflag = filtflag, 
                    bin_num = bin[0],
                    annot = bin[1],
                    ld_bin = bin[2],
                    maf_low = bin[3],
                    maf_high = bin[4],
                    part = part,
                    max_parts = max_parts
            }
        }
    }
    

    scatter (bin in ldsc) {
        call cat_grm_parts { 
            input:
                parts_done = make_grm_bins.parts_done,
                ancestry = ancestry,
                filtflag = filtflag,
                bin_num = bin[0],
                annot = bin[1],
                ld_bin = bin[2],
                maf_low = bin[3],
                maf_high = bin[4]
        }
    }
}


task make_grm_bins {

    input {
        String ancestry
        String filtflag
        String bin_num
        String annot
        String ld_bin
        String maf_low
        String maf_high
        String part
        String max_parts
    }

    command <<<
        set -euox pipefail

        # download relevant files
        gsutil -m cp  aou_bucket/heritability/~{ancestry}_covs~{filtflag}.txt \
                        aou_bucket/heritability/grm/unmerged/~{ancestry}~{filtflag}_newvep~{bin_num}*varlist \
                        aou_bucket/heritability/merged_~{ancestry}.p* .

          # extract list of individuals from covariates file
        awk '{print $1,$2}' *covs*.txt | tail -n +2 > keep_indivs.txt

        # make bed file for each maf-LD bin
        plink2 \
            --pfile merged_~{ancestry} \
            --keep keep_indivs.txt \
            --extract ~{ancestry}~{filtflag}_newvep~{bin_num}_~{annot}_~{ld_bin}_~{maf_low}_~{maf_high}.varlist \
            --make-bed \
            --out ~{ancestry}~{filtflag}_~{bin_num}

        # make GRM for each maf-LD bin and upload
        gcta64 \
            --bfile ~{ancestry}~{filtflag}_~{bin_num} \
            --make-grm-part ~{max_parts} ~{part} \
            --thread-num 16 \
            --out newvep~{bin_num}_~{ancestry}~{filtflag}_~{annot}_~{ld_bin}_~{maf_low}_~{maf_high} \
            --make-grm-alg 0 \
            --autosome

      
        gsutil -m cp newvep~{bin_num}_~{ancestry}_* aou_bucket/heritability/grm/unmerged/ 
    >>>

    output{
        String parts_done = "done"
        }

    runtime {
        docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/tvlupu/r_plink2_gcta:4.20.64"
        memory: "64 GiB"
        cpu: 16
        disks: "local-disk 1000 HDD"
        zones: "us-central1-a"
        cpuPlatform: "Intel Ice Lake"
        preemptible: 3
    }
}


## Concatenate GRM parts per maf-LD bin

task cat_grm_parts {
    input {
        Array[Array[String]] parts_done
        String ancestry
        String filtflag
        String bin_num
        String annot
        String ld_bin
        String maf_low
        String maf_high
    }

    command <<<
        set -euox pipefail
        gsutil -m cp aou_bucket/heritability/grm/unmerged/newvep~{bin_num}*_~{ancestry}*~{annot}* .

        # concatenate all GRM parts from make_snp_bins
        cat newvep~{bin_num}*_~{ancestry}_*.grm.id > merged_~{ancestry}~{filtflag}_newvep~{bin_num}_~{annot}_~{ld_bin}_~{maf_low}_~{maf_high}.grm.id
        cat newvep~{bin_num}*_~{ancestry}_*.grm.bin > merged_~{ancestry}~{filtflag}_newvep~{bin_num}_~{annot}_~{ld_bin}_~{maf_low}_~{maf_high}.grm.bin
        cat newvep~{bin_num}*_~{ancestry}_*.grm.N.bin > merged_~{ancestry}~{filtflag}_newvep~{bin_num}_~{annot}_~{ld_bin}_~{maf_low}_~{maf_high}.grm.N.bin 
    
        gsutil -m cp merged_* aou_bucket/heritability/grm/
    >>>

   runtime {
        docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/tvlupu/r_plink2_gcta:4.20.64"
        memory: "64 GiB"
        cpu: 16
        disks: "local-disk 1000 HDD"
        zones: "us-central1-a"
        cpuPlatform: "Intel Ice Lake"
    }
}
