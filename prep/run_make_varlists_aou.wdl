# reg-other 

version 1.0

workflow run_make_varlists_aou {
    input {
        String ancestry
        String filtflag
        Array[Array[String]] ldsc 
}

    scatter (bin in ldsc) {
        call make_varlist { 
            input:
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


task make_varlist {

    input {
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

          # download relevant files
        gsutil -m cp aou_bucket/heritability/~{ancestry}_score.ld \
                        aou_bucket/heritability/grm/ldbins_charAnnots_vep_regother.py \
                        aou_bucket/heritability/grm/ldbins_charAnnots_vep.py \
                        aou_bucket/heritability/annots/chr*.zst .
        
                
        if [ ~{annot} == "regulatory" ]
        then
            $other_bin=$((~{bin_num} + 4))
            for j in ~{annot}; do
                for i in {1..22} ; do
                    zstd --decompress --rm chr${i}_${j}.zst
                    awk '{print $0, "'$j'"}' chr${i}_${j} > tmp && mv tmp chr${i}_${j}   # add ANNOT col
                done
            done

            # file to exclude non-reg annotations 
            for i in {1..22} ; do
                zstd --decompress --rm chr${i}_*.zst
                cat chr${i}_* > chr${i}_vep
            done

            # make score.ld per chr to not crash script
            awk '{print > ("chr" $2 "_score.ld")}' ~{ancestry}~{filtflag}_score.ld
            
            # make lists of vars in annotation-maf-LD bins using ldbins_charAnnots.py
            # outputs one var per row > ~{bin_num}_allchrs_vep.varlist
            for i in {1..22} ; do python3 ldbins_charAnnots_vep_regother.py $i ~{bin_num} ~{annot} ~{ld_bin} ~{maf_low} ~{maf_high} ${other_bin} ; done 
            cat *~{bin_num}_vep.varlist > ~{bin_num}_allchrs_vep.varlist
            cat *${other_bin}_vep.varlist > ${other_bin}_allchrs_vep.varlist

            # upload but give informative name
            gsutil -m cp ~{bin_num}_allchrs_vep.varlist aou_bucket/heritability/grm/unmerged/~{ancestry}_newvep~{bin_num}_~{annot}_~{ld_bin}_~{maf_low}_~{maf_high}.varlist
            gsutil -m cp ${other_bin}_allchrs_vep.varlist aou_bucket/heritability/grm/unmerged/~{ancestry}_newvep${other_bin}_other_~{ld_bin}_~{maf_low}_~{maf_high}.varlist
        
        else
            for j in ~{annot}; do
                for i in {1..22} ; do
                    zstd --decompress --rm chr${i}_${j}.zst
                    awk '{print $0, "'$j'"}' chr${i}_${j} > tmp && mv tmp chr${i}_${j}   # add ANNOT col
                done
            done

            # make score.ld per chr to not crash script
            awk '{print > ("chr" $2 "_score.ld")}' ~{ancestry}~{filtflag}_score.ld
            
            # make lists of vars in annotation-maf-LD bins using ldbins_charAnnots.py
            # outputs one var per row > ~{bin_num}_allchrs_vep.varlist
            for i in {1..22} ; do python3 ldbins_charAnnots_vep.py $i ~{bin_num} ~{annot} ~{ld_bin} ~{maf_low} ~{maf_high} ; done 
            cat *~{bin_num}_vep.varlist > ~{bin_num}_allchrs_vep.varlist

            # upload but give informative name
            gsutil -m cp ~{bin_num}_allchrs_vep.varlist aou_bucket/heritability/grm/unmerged/~{ancestry}_newvep~{bin_num}_~{annot}_~{ld_bin}_~{maf_low}_~{maf_high}.varlist
        fi
    >>>

    runtime {
        docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/tvlupu/python_numpy_pandas:3.10"
        memory: "64 GiB"
        cpu: 16
        disks: "local-disk 500 HDD"
        zones: "us-central1-a"
        cpuPlatform: "Intel Ice Lake"
        preemptible: 3
    }
}
