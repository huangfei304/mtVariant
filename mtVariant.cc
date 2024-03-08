#include <algorithm>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>

#include "mtVariant.h"
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "Bam.h"
#include "CoverageCounter.h"
#include "Primer.h"
#include "EventScanner.h"
#include "ReferenceSequence.h"
#include "SnpEvent.h"
#include "IndelEvent.h"
#include "VcfWriter.h"
#include "MapStat.h"

namespace mtVariant {

static const char *help = "\
Required arguments:\n\
    -r, --ref REF           mitochondrion genome in FASTA format\n\
    -i, --in IN             Sorted and indexed input BAM or CRAM file\n\
    -s, --sample SN         Sample name to use in the output VCF file\n\
    -o, --outdir DIR        Output directory\n\
    -t, --type STR          Sequencing Type, fwd, rev or unknown(as both), [rev]\n\
\n\
Options:\n\
    -k, --rlen INT           The reads length [100]\n\
    -x, --extend INT         The extend n bp of chrM end with the first n bp bases for chrM genome re-build. [n=65]\n\
    -e, --read INT           Minimum reads length, after trim-primer. [40]\n\
    -p, --primer STR         Primer file for filter primer sequence and non-specific PCR sequence\n\
                             Format: chrM<tab>F_start<tab>F_end<tab>R_start<tab>R_end\n\
    -P, --pread  STR         Output file for read count of each primer\n\
    -O, --out STR            Output filtered read mapping result, bam format.\n\
    -n, --dist   INT         Flank distance of primer start allowed [5]\n\
    -q, --qual   INT         Minimum reads mapping quality [1]\n\
    -b, --base   INT         Minimum base quality for SNP [20]\n\
    -d, --depth  INT         Minimum variant depth [2]\n\
    -S, --site   INT         Minimum site depth [10]\n\
    -m, --mis    FLT         Maximum reads mismatch rate to filter read for non-HVR, HVR is 1.5-fold. [0.05]\n\
    -g, --gap    FLT         Maximum reads gap rate to filter read for non-HVR, HVR is 1.5-fold. [0.02]\n\
    -c, --cov    INT         Maximum read depth [65536]\n\
    -a, --rate   FLT         Minimum variant rate for indel [0.1]\n\
    //-l, --all                Output  all variant result\n\
    -h, --help               Show this help\n\
";

void usage()
{
    std::cerr << MTVARIANT_NAME << " v" << MTVARIANT_VERSION << std::endl
              << std::endl
              << help << std::endl;
}

int mtVariant_main(int argc, char **argv)
{
    char *ref = nullptr;
    char *bam_fn = nullptr;
    char *outdir = nullptr;
    char *sample_name = nullptr;
    char *primer_fn = nullptr;
    char *out_bam = nullptr;
    char *seq_type = nullptr;
    char *out_primer_file = nullptr;
    int tmp_min_mapq = 1;

    opts_s opts;
    int c, optidx;
    static struct option long_options[] = {
        {"ref", 1, nullptr, 0},                     // r
        {"in", 1, nullptr, 0},                      // i
        {"sample", 1, nullptr, 0},                  // s
        {"outdir", 1, nullptr, 0},                  // o
        {"type",1, nullptr, 0},                     // t
        {"rlen",1, nullptr, 0},                     // k
        {"extend",1, nullptr, 0},                   // x
       // {"seq",1, nullptr, 0},                    // j
        {"read", 1, nullptr, 0},                    // e
       //{"filter",0, nullptr,0},                    // f
        {"primer",1, nullptr, 0},                   // p
        {"pread",1, nullptr, 0},                    // P
        {"out",1, nullptr, 0},                      // O
        {"dist",1, nullptr, 0},                     // n
        {"qual", 1, nullptr, 0},                    // q
        {"base",1, nullptr, 0},                     // b
        {"depth",1, nullptr, 0},                    // d
        {"site", 1, nullptr, 0},                    // S
        {"mis", 1, nullptr, 0},                     // m
        {"gap", 1, nullptr, 0},                     // g
        {"cov", 1, nullptr, 0},                     // c
        {"rate", 1, nullptr, 0},                    // a
        //{"all", 0, nullptr, 0},                    // l
        {"help", 0, nullptr, 0},                    // h
        {nullptr, 0, nullptr, 0}};
    const char *short_options = "r:i:s:o:k:x:e:p:P:O:n:q:b:d:t:S:m:g:c:a:lh";
    const char *shorter_options = "risokxepPOnqbdtSmgcalh";

    if (argc == 1) {
        usage();
        return EXIT_SUCCESS;
    }

    while ((c = getopt_long(argc, argv, short_options, long_options, &optidx)) != -1) {
        if (c == 0) {
            c = shorter_options[optidx];
        }
        switch (c) {
        case 'r': // ref
            ref = optarg;
            break;
        case 'i': // in
            bam_fn = optarg;
            break;
        case 's': // sample-name
            sample_name = optarg;
            break;
        case 'b': // base quality
            opts.min_base_score = (uint8_t)std::atoi(optarg);
            break;
        case 'o': // prefix
            outdir = optarg;
            break;
        case 't': // sequencing type
            seq_type = optarg;
            break;
        case 'k': // read length
            opts.read_len_max = (int)std::atoi(optarg);
            break;
        case 'x': //extend
            opts.extend = (int)std::atoi(optarg);
            break;
        case 'e': // min read length
            opts.read_len_min = (int)std::atoi(optarg);
            break;
        case 'p': // primer file
            primer_fn = optarg;
            opts.filter_primer = true;
            break;
        case 'P': // output primer read count file
            out_primer_file = optarg;
            break;
        case 'O': // output bam file
            out_bam = optarg;
            break;
        case 'n': // primer distance
            opts.distance = (int)std::atoi(optarg);
            break;
        case 'q': // min-map-qual
            tmp_min_mapq = (uint8_t)std::atoi(optarg);
            break;
        case 'd': // min-variant-depth
            opts.min_var_dp = (uint8_t)std::atoi(optarg);
            break;
        case 'S': // min-site-depth
            opts.min_site_dp = (int)std::atoi(optarg);
            break;
        case 'm': // max read mismatch rate
            opts.read_max_mis_ratio = std::atof(optarg);
            break;
        case 'g': // max read gap rate
            opts.read_max_gap_ratio = std::atof(optarg);
            break;
        case 'c': // max-coverage
            opts.max_cov = std::atoi(optarg);
            break;
        case 'a': // min variand read rate for indel
            opts.indel_min_var_ratio = std::atof(optarg);
            break;
        //case 'l': // all_variant
        //    opts.all_variant = true;
        //    break;
        case 'h': // help
            usage();
            return EXIT_SUCCESS;
        default:
            return EXIT_FAILURE;
        }
    }

    if (bam_fn != nullptr) {
        std::cerr << "Input alignment file: " << bam_fn << std::endl;
    } else {
        std::cerr << "No input alignment file given" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (ref != nullptr) {
        std::cerr << "Reference file: " << ref << std::endl;
    } else {
        std::cerr << "No reference file given" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (sample_name != nullptr) {
        std::cerr << "Sample name: " << sample_name << std::endl;
    } else {
        std::cerr << "No sample name given" << std::endl;
        exit(EXIT_FAILURE);
    }
    if( seq_type !=nullptr ){
        if( strcmp(seq_type, "fwd") ==0 || strcmp(seq_type, "rev")==0 || strcmp(seq_type,"unknown")==0 ){
            std::cerr << "Sequencing type: " << seq_type <<std::endl;
        }else{
            std::cerr<< "Sequencing type: "<< seq_type <<" is not right, must be one of fwd,rev,unknown."<<std::endl;
            exit(EXIT_FAILURE);
        }
    } else {
        seq_type = new char[strlen("rev") + 1];
        strcpy(seq_type, "rev");
        std::cerr<<"The sequencing type is not set, auto-set to rev by default."<<std::endl;
    }
    if( opts.extend < 0 ){
        std::cerr << "The extend size must >= 0 " << std::endl;
        exit(EXIT_FAILURE);
    }
    if (tmp_min_mapq < 0 || tmp_min_mapq > 60) {
        std::cerr << "Minimum mapping quality must be in the range [0, 60]" << std::endl;
        exit(EXIT_FAILURE);
    }
    opts.min_mapq = tmp_min_mapq;
    if (outdir == nullptr) {
        std::cerr << "No output directory given" << std::endl;
        exit(EXIT_FAILURE);
    }else{
        size_t dir_len = strlen(outdir)+12;
        char* command=(char*)calloc(dir_len,sizeof(char));
	   	sprintf(command, "mkdir -p %s/", outdir);
        system(command);
    }
    if( opts.filter_primer ){
        std::string tmp_file(primer_fn);
        if( !std::ifstream(tmp_file)){
            std::cerr<<"Primer file:"<<primer_fn<<" is not exist."<<std::endl;
            exit(EXIT_FAILURE);
        }
    }

    Bam bam(bam_fn, ref);
    size_t fn_len = strlen(outdir) + strlen(sample_name) + 10;
    char* out_fn=(char*)calloc(fn_len, sizeof(char));
    sprintf(out_fn, "%s/%s.vcf.gz", outdir, sample_name);

    vcfFile *vcf_fp;
    const char *write_mode = "wz";
    if ((vcf_fp = hts_open(out_fn, write_mode)) == nullptr) {
        std::cerr << "Unable to open output VCF file" << std::endl;
        exit(EXIT_FAILURE);
    }
    bcf_hdr_t *vcf_hdr;
    if ((vcf_hdr = bcf_hdr_init("w")) == nullptr) {
        std::cerr << "Unable to create variant header" << std::endl;
        exit(EXIT_FAILURE);
    }

    bcf_hdr_set_version(vcf_hdr, "VCFv4.3");
    //regions_list_t regions;
    //regions.reserve(bam._hdr->n_targets);
    std::cerr << "Found " << std::to_string(bam._hdr->n_targets) << " regions in bam" << std::endl;
    int max_len = 16569; // chrM size
    for (int i = 0; i < bam._hdr->n_targets; ++i) {
        std::cerr<<"Chr:"<<std::string(bam._hdr->target_name[i])<<", Length:"<<bam._hdr->target_len[i]<<std::endl;
        //regions.push_back(std::make_pair(std::string(bam._hdr->target_name[i]), bam._hdr->target_len[i]));
        //if (bam._hdr->target_len[i] > max_len ) max_len = bam._hdr->target_len[i];
    }


    ReferenceSequence refseq(ref);
    EventScanner events(opts.read_max_mis_ratio,opts.read_max_gap_ratio,opts.min_base_score,opts.read_len_min,opts.filter_primer);
    CoverageCounter coverages((size_t)max_len, opts.max_cov);
    Primer primer(opts.distance, opts.read_len_max, seq_type);
    VcfWriter writer(coverages, refseq, events, vcf_hdr, vcf_fp, &opts);

    char* prefix=(char*)calloc(fn_len, sizeof(char));
    sprintf(prefix, "%s/%s", outdir, sample_name);
    MapStat mapstat(coverages, &opts, prefix);

    samFile *out;
    char omode[3];
    memset(omode,'\0', 3);
    strncpy(omode, "wb", 2);
    if( out_bam != nullptr ){
        std::string outBam(out_bam);
        size_t fn_len = strlen(out_bam);
        if(fn_len >= 4 && strcmp((out_bam + fn_len - 4),".bam") !=0 ){
            outBam +=".bam";
        }
        out = sam_open(outBam.c_str(), omode);
        if ( !out ) {
            std::cerr<<"failed to open \""<<outBam<<" for output"<<std::endl;
            return 1;
        }
        if (sam_hdr_write(out, bam._hdr) < 0) {
            std::cerr<<"failed to write header"<<std::endl;
            return 1;
        }
    }

    std::string chr_name="chrM";
    refseq.set_region(chr_name.c_str());
    std::cerr<<"Reading Bam for variant calling..."<<std::endl;
    if( opts.filter_primer ){
        primer.read_chr_primer(primer_fn,chr_name);
    }

    const int chr_idx = sam_hdr_name2tid(bam._hdr,chr_name.c_str());
    const int16_t filter_flags = BAM_FUNMAP | BAM_FQCFAIL;
    const int16_t second_flags = BAM_FSECONDARY | BAM_FSUPPLEMENTARY;
    bam1_t *curr_rec = bam_init1();
    while( sam_read1(bam._sf, bam._hdr, curr_rec)>=0 ){
        bool supplementary = false;
        //std::string read_name = bam_get_qname(curr_rec);
        if( curr_rec->core.flag  & second_flags ) supplementary = true;
        if( !supplementary ) ++opts.total_reads;
        if( curr_rec->core.flag & filter_flags ){
            ++opts.fail_unmap_reads;
        }else if( curr_rec->core.flag & BAM_FSECONDARY ){
            ++opts.no_primary;
        }else{
            if( !supplementary ) ++opts.map_reads;
            if( !supplementary && (curr_rec->core.flag  & BAM_FPROPER_PAIR) ) ++opts.paired_reads;
            if( curr_rec->core.qual > opts.min_mapq ){
                if( !supplementary ){
                    ++opts.mapqual_reads;
                    if (bam_is_rev(curr_rec)) ++opts.rev_reads;
                    else ++opts.fwd_reads;
                }
            }
        }
        //if( !supplementary && curr_rec->core.l_qseq < opts.read_len_min ) ++opts.read_short;

        // for chrM SNP and Indel
        if( curr_rec->core.flag & filter_flags ) continue; // unmap
        if( curr_rec->core.flag & BAM_FSECONDARY ) continue; // secondary alignment
        if( curr_rec->core.qual < opts.min_mapq ) continue; // low mapping quality
        if( chr_idx != curr_rec->core.tid ) continue;

        uint8_t *p;
        if ((p = bam_aux_get(curr_rec,"AS")) != NULL) {
	        // char* xa_tag;
	        // xa_tag = bam_aux2Z(p);
	        // std::string xa_tag_str(xa_tag);
            //std::cerr<<"read name:"<<bam_get_qname(curr_rec)<<",Tag:"<<xa_tag<<std::endl;
            if( bam_aux2i(p) < 20 ) continue;
        }
        //std::cerr<<"read name:"<<bam_get_qname(curr_rec)<<",pos:"<<std::to_string(curr_rec->core.pos)<<std::endl;
        events.collect_variant(curr_rec, refseq, coverages, primer);
        int idx = events.get_filter_flag();
        if( !idx && out_bam !=nullptr ){
            int i = sam_write1(out, bam._hdr, curr_rec);
        }else{
            //std::cerr<<"read_name:"<<bam_get_qname(curr_rec)<<",idx:"<<std::to_string(idx)<<std::endl;
        }
    }

    std::cerr<<"Writing VCF...\n";
    //for median depth
    opts.read_short = events.get_short_read_count();
    mapstat.calCoverage(chr_name, max_len);
    mapstat.printDepth(chr_name, max_len, false);
    opts.min_low_depth = 2 * mapstat.get_median_depth() / 10;
    if( opts.min_low_depth <  opts.min_site_dp ) opts.min_low_depth = opts.min_site_dp;

    writer.setup_vcf(sample_name, argc, argv, MTVARIANT_VERSION, ref);
    if (bcf_hdr_write(vcf_fp, vcf_hdr)) {
        std::cerr << "Error: Failed to write VCF header" << std::endl;
        exit(EXIT_FAILURE);
    }

    writer.print_variant_buffer(chr_name, max_len);
    if( out_bam !=nullptr ){
        if (sam_close(out) < 0) {
            fprintf(stderr, "error closing output file\n");
            exit(EXIT_FAILURE);
        }
    }
    if( opts.filter_primer && out_primer_file !=nullptr ) {
        FILE *fout = fopen(out_primer_file,"w");
        if( !fout ) std::cerr<<"Fail to open primer depth file:"<<out_primer_file<<std::endl;
        int size = primer._primer_list.size();
        int total = 0;
        for(int i=0;i<size;++i){
            total += coverages.get_coverage(i, COVSIDX_PRIMER);
        }
        fprintf(fout,"F_Primer_Start\tR_Primer_End\tReadCount\tPercentage(%)\n");
        for(int i=0; i<size; ++i ){
            fprintf(fout,"%d\t%d\t%d\t%.4f\n", primer._primer_list[i].fwd_start,primer._primer_list[i].rev_end,coverages.get_coverage(i, COVSIDX_PRIMER),100.0*coverages.get_coverage(i, COVSIDX_PRIMER)/(1.0*total));
        }
        fclose(fout);
    }
    coverages.reset();
    events.reset();

    mapstat.printStat();

    hts_close(vcf_fp);
    bcf_hdr_destroy(vcf_hdr);

    std::cerr << "Finished" << std::endl;
    return EXIT_SUCCESS;
}

} // namespace mtVariant

int main(int argc, char **argv)
{
    return mtVariant::mtVariant_main(argc, argv);
}
