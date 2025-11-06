#!/usr/bin/env python

##########################################################################
### Script to extract CB and UMI from read's header for Parse datasets ###
##########################################################################

#Import the required libraries:
import argparse
import pysam

def extract_tags(read_name, delimiter):
    """Extract the cell barcode from read name using the given delimiter."""
    try:
        return read_name.split(delimiter)[-2], read_name.split(delimiter)[-1]

    except Exception:
        return None

def main():
    parser = argparse.ArgumentParser(description="Extract CB tag from read name and add to BAM")
    parser.add_argument("-i", "--input", required=True, help="Input BAM file")
    parser.add_argument("-o", "--output", required=True, help="Output BAM file")
    parser.add_argument("-d", "--delimiter", default="__", help="Delimiter used to split read name [default: '__']")

    args = parser.parse_args()

    in_bam = pysam.AlignmentFile(args.input, "rb")
    out_bam = pysam.AlignmentFile(args.output, "wb", header=in_bam.header)

    for read in in_bam:

        #Extract info from header:
        cb, ub = extract_tags(read.query_name, args.delimiter)

        #Add it as tags:
        if cb:
            read.set_tag("CB", cb, value_type='Z')

        if ub:
            read.set_tag("UB", ub, value_type='Z')

        out_bam.write(read)

    in_bam.close()
    out_bam.close()
    print(f"Done. Written output to {args.output}")

if __name__ == "__main__":
    main()
