#!/usr/bin/perl -w

use Data::Dumper;
use File::Path qw(make_path);
use File::Basename;
use Getopt::Long;

my $data_dir;
my $writeMappings = 0;

GetOptions("data_dir=s" => \$data_dir,
		   "writeMappings" => \$writeMappings) or die;

die "Please provide a single parameter, a folder with sequencing data." if not $data_dir;

my @sequencing_files = <$data_dir/sequencing/*>;

make_path("$data_dir/salmon_output");

for (@sequencing_files) {
	my $file_name = basename($_);
	my $salmon_quant = "$data_dir/salmon_output/$file_name";
	my $alignment_file = "$data_dir/alignments/$file_name.sam";
	my $alignment_file_bam = "$data_dir/alignments/$file_name.bam";

	my $salmon_command = "salmon quant -i salmon_indexes/DK_CRISPR_salmon_index -l U -r $_ -o $salmon_quant";
	if ($writeMappings) {
		make_path("$data_dir/alignments");
		$salmon_command = "$salmon_command --writeMappings=$alignment_file";
	}
	system($salmon_command);
	
	if ($writeMappings) {
		my $samtools_command = "samtools view -@ 12 -S -b $alignment_file > $alignment_file_bam";
		system($samtools_command);
		unlink($alignment_file);
	}
}
