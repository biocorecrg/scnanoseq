
get_fastqc_counts()
{
    fastqc_file=$1
    counts=$(unzip -p ${fastqc_file} $(basename ${fastqc_file} .zip)/fastqc_data.txt | \
        grep 'Total Sequences' | \
        cut -f2 -d$'\t')
    echo $counts

}

get_nanoplot_counts()
{
    nanoplot_file=$1
    counts=$(grep 'Number of reads' $nanoplot_file | awk '{print $NF}' | cut -f1 -d'.' | sed 's/,//g')
    echo $counts
}

get_flagstat_counts()
{
    flagstat_file=$1
    counts=$(grep 'primary mapped' $flagstat_file | cut -d ' ' -f1)
    echo $counts
}

output=""
input=""

while [[ $# -gt 0 ]]
do
    flag=$1

    case "${flag}" in
        --input) input=$2; shift;;
        --output) output=$2; shift;;
        --type) type=$2; shift;;

        *) echo "Unknown option $1 ${reset}" && exit 1
    esac
    shift
done

header=""
data=""

header="Sample,Input_reads,Tagged_reads,Mapped_reads,Dedup_reads"
echo "$header" > $output

for sample_name in $(for file in $(readlink -f $input)/*.zip; do basename $file; done | cut -f1 -d'.' | sort -u)
do

    ###############
    # INPUT_FILES #
    ###############

    raw_fastqc="${sample_name}.raw_fastqc.zip"
    raw_nanoplot="${sample_name}.raw_NanoStats.txt"

    extract_fastqc="${sample_name}.extracted_fastqc.zip"
    extract_nanoplot="${sample_name}.extracted_NanoStats.txt"

    correct_csv="${sample_name}.corrected_bc_umi.tsv"
    data="$(basename $sample_name)"
    
    mapped_stats="${sample_name}.${type}.tagged.flagstat"
    dedup_stats="${sample_name}.${type}.umitools_dedup.flagstat"
    
    ####################
    # RAW FASTQ COUNTS #
    ####################
    if [[ -s "$raw_fastqc" ]]
    then
        fastqc_counts=$(get_fastqc_counts "$raw_fastqc")
        data="$data,$fastqc_counts"
    elif [[ -s "$raw_nanoplot" ]]
    then
        nanoplot_counts=$(get_nanoplot_counts "$raw_nanoplot")
        data="$data,$nanoplot_counts"
    else
        data="$data,"
    fi

    #################
    # TAGGED COUNTS #
    #################
    if [[ -s $correct_csv ]]
    then
        correct_counts=$(cut -f6 $correct_csv | awk '{if ($0 != "") {print $0}}' | wc -l)
        data="$data,$correct_counts"
    elif [[ -s "$extract_fastqc" ]]
    then
        extract_counts=$(get_fastqc_counts "$extract_fastqc")
        data="$data,$extract_counts"
    elif [[ -s "$extract_nanoplot" ]]
    then
        nanoplot_counts=$(get_nanoplot_counts "$extract_nanoplot")
        data="$data,$nanoplot_counts"
    else
        data="$data,"
    fi

    #################
    # MAPPED COUNTS #
    #################
    if [[ -s $mapped_stats ]]
    then
        mapped_counts=$(get_flagstat_counts "$mapped_stats")
        data="$data,$mapped_counts"
    else
        data="$data,"
    fi

    #################
    # DEDUP COUNTS #
    #################
    if [[ -s $dedup_stats ]]
    then
        dedup_counts=$(get_flagstat_counts "$dedup_stats")
        data="$data,$dedup_counts"
    else
        data="$data,"
    fi
    
    echo "$data" >> $output
done
