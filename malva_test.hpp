/**
 * MALVA - genotyping by Mapping-free ALternate-allele detection of known VAriants
 * Copyright (C) 2019  Giulia Bernardini, Luca Denti, Marco Previtali
 *
 * This file is part of MALVA.
 *
 * MALVA is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MALVA is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MALVA; see the file LICENSE. If not, see
 * <https://www.gnu.org/licenses/>.
 **/

#ifndef _MALVA_TEST_HPP_
#define _MALVA_TEST_HPP_

#include <iostream>
#include <stdexcept>
#include <string.h>
#include <sys/resource.h>
#include "vcf.h"

/*
 Returns the peak (maximum so far) resident set size (physical memory use) measured in Megabytes.
 */
long get_mem_usage(){
    struct rusage myusage;
    
    getrusage(RUSAGE_SELF, &myusage);
    //Return the maximum resident set size used (in kilobytes).
    return myusage.ru_maxrss;
}

void usage(const string name) {
    std::cerr << "Usage: " << name << std::endl;
}


/* INPUT -> TWO FILEs VCF
   OUTPUT -> % Precision VCFs
 */
int compare_genotypes(const char* sample_vcf, const char* geno_vcf){
    
    //PRINT USED FILEs
    usage(sample_vcf);
    usage(geno_vcf);
    
    //SAMPLE INIT
    htsFile *sample_bcf = NULL;
    bcf_hdr_t *sample_header = NULL;
    bcf1_t *sample_record = bcf_init();
    
        //GENO INIT
        htsFile *geno_bcf = NULL;
        bcf_hdr_t *geno_header = NULL;
        bcf1_t *geno_record = bcf_init();
    
    //SAMPLE OPEN FILE
    sample_bcf = bcf_open(sample_vcf, "r");
    if(sample_bcf == NULL) {
        throw std::runtime_error("Unable to open sample file.");
    }
        //GENO OPEN FILE
        geno_bcf = bcf_open(geno_vcf, "r");
        if(geno_bcf == NULL) {
            throw std::runtime_error("Unable to open genotype file.");
        }
    
    //SAMPLE HEALDER
    sample_header = bcf_hdr_read(sample_bcf);
    if(sample_header == NULL) {
        throw std::runtime_error("Unable to read sample header.");
    }
        //GENO HEALDER
        geno_header = bcf_hdr_read(geno_bcf);
        if(geno_header == NULL) {
            throw std::runtime_error("Unable to read genotype header.");
        }
    
    int rec_sa = 0; //number of sample records 
    int match = 0; //same record match
    
    /***
        * 1. read sample_RECORD
        * 2. check sample_RECORD are covered in             the geno_RECORD
        * 3. if covered match++
        * 4. output the precision of covered
    ***/
    while(bcf_read(sample_bcf, sample_header, sample_record) == 0) {
        /***
        * 1. extract sample_RECORD (SR)
        * 2. search SR in GENO
        * 3. every compared record: rec_sa++
        * 4. every succes match: match++
        ***/
        while(bcf_read(geno_bcf, geno_header, geno_record) == 0) {
            rec_sa++;
            //COMPARE CHROM (int64_t)
            if(sample_record->rid == geno_record->rid){
                //COMPARE POS (int64_t)
                if(sample_record->pos == geno_record->pos){
                    //COMPARE length of REF (int64_t)
                    if(sample_record->rlen == geno_record->rlen){
                        match++;
                        break;
                    }
                }
            }
        }
    }
    
    //SAMPLE DESTROY
    bcf_hdr_destroy(sample_header);
    bcf_destroy(sample_record); 
    bcf_close(sample_bcf);
        //GENO DESTROY
        bcf_hdr_destroy(geno_header);
        bcf_destroy(geno_record); 
        bcf_close(geno_bcf);
    
    if(rec_sa == 0){//sample record empty
        return 0;
    }
    return (match/rec_sa)*100; //% precision match
}



#endif
