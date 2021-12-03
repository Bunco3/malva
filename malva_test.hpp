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
#include <sys/time.h>
#include "vcf.h"

const int DEBUG = 0; //=1 FOR DEBUG PRINTs

typedef struct {
    htsFile *bcf;
    bcf_hdr_t *header;
    bcf1_t *record;
} VCFt;

long get_mem_usage();
double get_cpu_time();

void compare_genotypes(const char* sample_vcf, const char* geno_vcf);
VCFt read_vcf(const char* vcf);

/*
 Returns the peak (maximum so far) resident set size (physical memory use) measured in Megabytes.
 */
long get_mem_usage(){
    struct rusage myusage;
    
    getrusage(RUSAGE_SELF, &myusage);
    //Return the maximum resident set size used (in kilobytes).
    return myusage.ru_maxrss;
}

/*
 This is the total amount of time spent executing in user mode, expressed in a timeval structure. 
 struct timeval {
    time_t      tv_sec;     // seconds
    suseconds_t tv_usec;    // microseconds }; 
 */
double get_cpu_time(){
    struct rusage myusage;
    getrusage(RUSAGE_SELF, &myusage);

    long seconds = myusage.ru_stime.tv_sec;
    long microseconds = myusage.ru_stime.tv_usec;
    
    double time = seconds + (microseconds*1e-6);
    return time;
}

/* INPUT -> TWO FILEs VCF
   OUTPUT -> Print of % Precision match
 */
void compare_genotypes(const char* sample_vcf, const char* geno_vcf){
    
    //PRINT USED FILEs
    std::cerr << "Compare \"" << sample_vcf << "\" with \"" << geno_vcf << "\"" << std::endl;
    
    //LOAD SAMPLE
    VCFt sample = read_vcf(sample_vcf);
    
    double rec_sa = 0; //number of sample records 
    double match = 0; //records matched
    
    /***
        * 1. read sample_RECORD
        * 2. check sample_RECORD are covered in the geno_RECORD
        * 3. if covered match++
        * 4. output the precision of covered
    ***/
    while(bcf_read(sample.bcf, sample.header, sample.record) == 0) {
        /***
        * 1. extract sample_RECORD (SR)
        * 2. search SR in GENO
        * 3. every compared record: rec_sa++
        * 4. every succes match: match++
        ***/
        rec_sa++;
       
        //LOAD GENO
        VCFt geno = read_vcf(geno_vcf);
        
        while(bcf_read(geno.bcf, geno.header, geno.record) == 0) {
            
            //DEBUG, print all records comparation
            if(DEBUG){
                std::cerr << "Record sample #CHROM " << sample.record->rid+1 << ", #POS " << sample.record->pos+1 << ", LENGHT #REF " << sample.record->rlen << std::endl;
                
                std::cerr << "Record Geno #CHROM " << geno.record->rid+1 << ", #POS " << geno.record->pos+1 << ", LENGHT #REF " << geno.record->rlen << std::endl << std::endl;
            }
            
            //COMPARE CHROM (int64_t)
            if(sample.record->rid+1 == geno.record->rid+1){
                //COMPARE POS (int64_t)
                if(sample.record->pos+1 == geno.record->pos+1){
                    //COMPARE length of REF (int64_t)
                    if(sample.record->rlen == geno.record->rlen){
                        //DEBUG
                        if(DEBUG){
                            std::cerr << "<< RECORD MATCH FOUND >>" << std::endl << std::endl;
                        }
                        match++;
                        break;
                    }
                }
            }
        }
        //GENO DESTROY
        bcf_hdr_destroy(geno.header);
        bcf_destroy(geno.record); 
        bcf_close(geno.bcf);
    }
    
    //SAMPLE DESTROY
    bcf_hdr_destroy(sample.header);
    bcf_destroy(sample.record); 
    bcf_close(sample.bcf);
        
    //PRINT RESULTS
    if(rec_sa == 0){//sample records empty
        std::cerr << "SAMPLE RECORDs EMPTY!" << std::endl;
    }else{
        std::cerr << "Records Matched " << match << std::endl << "Records Processed " << rec_sa << std::endl;
        std::cerr << "Value of Precision: " << 100*(match/rec_sa) << "%" << endl;
    }
}

/* Initialize the vcf file and create a VCFf struct for easier reading.
 * INPUT -> VCF file path
   OUTPUT -> struct {bcf, header, record}
 */
VCFt read_vcf(const char* vcf){
    
    VCFt out = {NULL, NULL, NULL};
        
    //VCF INIT
    htsFile *vcf_bcf = NULL;
    //SAMPLE OPEN FILE
    vcf_bcf = bcf_open(vcf, "r");
    if(vcf_bcf == NULL) {
        throw std::runtime_error("Unable to open vcf file.");
    }
    out.bcf = vcf_bcf;
    
    bcf_hdr_t *vcf_header = NULL;
    //VCF HEALDER
    vcf_header = bcf_hdr_read(vcf_bcf);
    if(vcf_header == NULL) {
        throw std::runtime_error("Unable to read vcf header.");
    }
    out.header = vcf_header;
    
    bcf1_t *vcf_record = bcf_init();
    out.record = vcf_record;
    
    return out;
}

#endif
