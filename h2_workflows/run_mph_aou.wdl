version 1.0

workflow run_mph_aou {
    input {
        String ancestry
        String pheno1
        String pheno2
        String outflag
    }

    call run_mph_her {
        input:
            ancestry = ancestry,
            pheno1 = pheno1,
            pheno2 = pheno2,
            outflag = outflag
    }
}   


task run_mph_her {
    input {
        String ancestry
        String pheno1
        String pheno2
        String outflag
    }

    command <<<
    
    # height raw: pheno ---> height
    # height rintres: pheno ---> height_rintres
    
        set -euox pipefail 
        
        fclist="aou_bucket/heritability"
        gsutil -m cp ${fclist}/mph/mph-v0.55.0-x86_64-linux.zip .
        unzip mph-v0.55.0-x86_64-linux.zip
        chmod +x mph
        
        # phenotypes            
        gsutil -m cp ${fclist}/phenos/~{pheno1}_~{ancestry}_mph.csv \
                    ${fclist}/covs_~{ancestry}_mph.csv \
                    ${fclist}/mph_convert_gcta.py . 
                    
        # get GRMs
            # coding
        gsutil -m cp ${fclist}/grm/mph/merged_~{ancestry}_newvep_alpha0{1..4}*grm* .
            # rna_exons / sepRNA (*RNA_exons)
        gsutil -m cp ${fclist}/grm/newvep_2/merged_~{ancestry}*newvep*rna_exons*mph*grm* .
        
            # sepcod
        #gsutil -m cp ${fclist}/grm/mph/merged_~{ancestry}_*sepcod2_lof*grm* \
                     #${fclist}/grm/mph/merged_~{ancestry}_*sepcod2*added*grm* \
                     #${fclist}/grm/mph/merged_~{ancestry}_*sepcod2*removed*grm* \
                     #${fclist}/grm/mph/merged_~{ancestry}_*sepcod2_miss*grm* \
                     #${fclist}/grm/mph/merged_~{ancestry}_*sepcod2_oth*grm* .   
            
            # prox, intronic, reg, non-annotated
        gsutil -m cp ${fclist}/grm/newvep_2/merged_~{ancestry}_newvep*prox*mph*grm* \
                     ${fclist}/grm/newvep_2/merged_~{ancestry}_newvep*intronic*mph*grm* \
                     ${fclist}/grm/mph/merged_~{ancestry}_newvep_alpha{21..28}*grm* .

            # chrX
        #gsutil -m cp ${fclist}/grm/mph/*~{ancestry}*chrX*mph*grm* .
                    
        
        
# modify GRMs from GCTA format to MPH format and upload
        #gsutil -m cp ${fclist}/grm/newvep_2/*~{ancestry}*sepcod2*grm* .
        #for file in *.grm.N.bin ; do
            #python3 mph_convert_gcta.py "${file%.grm.N.bin}" "${file%.grm.N.bin}_mph" 
            #rm $file "${file%.N.bin}.bin" "${file%.N.bin}.id" ;
            #gsutil -m cp "${file%.grm.N.bin}_mph.grm.iid" ${fclist}/grm/mph/ ;
            #gsutil -m cp "${file%.grm.N.bin}_mph.grm.bin" ${fclist}/grm/mph/ ;
        #done
        
        
  
# run MPH heritability for pheno 1        
        ls *.grm.iid | awk -F '.grm.iid' '{print $1}' > multi_GRMs.txt
        cat multi_GRMs.txt

        ./mph \
            --grm_list multi_GRMs.txt \
            --phenotype_file ~{pheno1}_~{ancestry}_mph.csv \
            --trait ~{pheno1} \
            --covariate_file covs_~{ancestry}_mph.csv \
            --minque \
            --save_memory \
            --num_threads 32 \
            --output_file ~{ancestry}_~{pheno1}_mph_~{outflag}
           
        gsutil -m cp *~{pheno1}*mq* ${fclist}/newvep_2/mph/
        
        
# run MPH heritability for pheno 2
        gsutil -m cp ${fclist}/phenos/~{pheno2}_~{ancestry}_mph.csv .
        ./mph \
            --grm_list multi_GRMs.txt \
            --phenotype_file ~{pheno2}_~{ancestry}_mph.csv \
            --trait ~{pheno2} \
            --covariate_file covs_~{ancestry}_mph.csv \
            --minque \
            --save_memory \
            --num_threads 32 \
            --output_file ~{ancestry}_~{pheno2}_mph_~{outflag}

        gsutil -m cp *~{pheno2}*mq* ${fclist}/newvep_2/mph/
       
        
    >>>

    runtime {
        docker: "us-central1-docker.pkg.dev/all-of-us-rw-prod/aou-rw-gar-remote-repo-docker-prod/tvlupu/python_numpy_pandas:3.10"
        memory: "256 GiB"
        cpu: 32
        disks: "local-disk 2000 HDD"
        zones: "us-central1-a"
        cpuPlatform: "Intel Ice Lake"
    }
}
