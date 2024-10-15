#!/usr/bin/env python

# Import required modules
import argparse
import gzip  # For reading gzipped files


def rrBLUP_hapmap(vcf_file, outfile):
    print("Writing the hapmap file using rrBLUP encoding.")

    # Create filename
    output_filename = outfile + '_hmp.txt'

    # Open input gzipped VCF file and output file
    with gzip.open(vcf_file, 'rt') as vcf, open(output_filename, 'w') as handle:

        # Lists for handling chromosome names
        chrom_l = []

        # Start reading the VCF file
        for line in vcf:
            # Skip meta-information lines
            if line.startswith('##'):
                continue

            # Header line with sample information
            elif line.startswith('#CHROM'):
                # Split the line on tabs
                header = line.strip().split('\t')
                format_field = header.index('FORMAT')
                samples = header[format_field + 1:]

                # Write header for the hapmap output
                handle.write('rs#\tallele\tchrom\tpos\t' + '\t'.join(samples) + '\n')

            # Process genotype data
            else:
                # Create a list to store the data for this line
                toprint = []

                # Split the variant data
                fields = line.strip().split('\t')
                chrom = fields[0]
                position = fields[1]
                ref_allele = fields[3]
                alt_allele = fields[4]
                genotypes = fields[9:]

                # Handle chromosome names
                if chrom not in chrom_l:
                    chrom_l.append(chrom)
                chrom_name = str(chrom_l.index(chrom) + 1)

                # Create a name for the SNP and combine ref/alt alleles
                snp_id = 'S' + chrom_name + '_' + position
                alleles = ref_allele + '/' + alt_allele

                # Add SNP ID, alleles, chromosome, and position to output list
                toprint.extend([snp_id, alleles, chrom_name, position])

                # Process each genotype and encode in rrBLUP format
                for g in genotypes:
                    call = g.split(':')[0]  # Extract the genotype call
                    if call == '0/0':
                        toprint.append('-1')  # Homozygous reference
                    elif call == '0/1' or call == '1/0':
                        toprint.append('0')  # Heterozygous
                    elif call == '1/1':
                        toprint.append('1')  # Homozygous alternate
                    else:
                        toprint.append('NA')  # Missing data or unknown case

                # Write the data line to the output file
                handle.write('\t'.join(toprint) + '\n')

    print(f"File was written as {output_filename}")


# Argument parser setup
parser = argparse.ArgumentParser(description="Convert a gzipped VCF file to a hapmap file in rrBLUP format.")
parser.add_argument('-i', '--vcf_in', help='Input gzipped VCF file', required=True)
parser.add_argument('-o', '--outfile', help='Output file basename (without extension)', required=True)
args = parser.parse_args()

# Execute the rrBLUP conversion
rrBLUP_hapmap(args.vcf_in, args.outfile)

