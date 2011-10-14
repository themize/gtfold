#!/usr/bin/perl
package SuboptNumStructures;
use strict;
use warnings;
use File::Basename;

sub test()
{
  my(%Config) = %{$_[1]};
  my(%Sequences) = %{$_[2]};
  my(%local_sequences) = %{$_[3]};
  my $logger = $_[4];

  my $gtdir = $Config{"G_GTFOLD_DIR"};
  my $rnadir = $Config{"G_RNAFOLD_DIR"};
  my $workdir = $Config{"G_WORK_DIR"};

  my $key;
  my $value;
  my %new_hash = (%local_sequences);
  # For testing pseudoknot detection logic
  # we may only test the sequences specified in local_sequences
  while (($key, $value) = each(%new_hash)) {

	  my $seqname=$key;
	  my $path;
	  my $suffix;

	  my $seqfile = $value;
    my $dirname = dirname($seqfile);
	  my $gtout  = "$workdir$seqname-gt";
    my $rnaout = "$workdir$seqname-rna";
	  my $gtoutfilename  = $workdir."$seqname-gt.ct";

    my $energy;
    for ($energy = 1; $energy <=10; ++$energy) {

      my $gtfile = $gtout.$energy;
      my $rnafile = $rnaout.$energy."_ss.txt";
  	  my $gtcmd;
      my $rnacmd;
      $gtcmd  = "$gtdir/gtfold --subopt $energy $seqfile -o $gtfile > /dev/null 2>&1";
      $rnacmd  = "$rnadir/RNAsubopt -s DAT -e $energy < $seqfile > $rnafile";

	    system("$gtcmd");
	    system("$rnacmd");

      my $gtsubopt_file = $gtfile."_ss.txt";
      my $rnasubopt_file = $rnafile;

      my $gtsorted = $gtsubopt_file."_sorted";
      my $rnasubopt_sorted = $rnafile."_sorted";
      `cat $gtsubopt_file | sed '/[A-Z]/d' | sort -o $gtsorted`;
      `cat $rnasubopt_file | sed '/[A-Z]/d' | sort -o $rnasubopt_sorted`;
      #`sort $gtsubopt_file -o $gtsorted`;
      #`sort $rnasubopt_file -o $rnasubopt_sorted`;

      my $diff_str = `diff $gtsorted $rnasubopt_sorted`;
      #print $diff_str;

      if ($diff_str eq "") {
        $logger->info("TEST PASSED: $seqname: energy delta = $energy: Suboptimal Structures matched");
      }
      else {
        $logger->error("TEST FAILED: $seqname: energy delta = $energy: Suboptimal Structures not matched\n $diff_str");
      }
    }
  }
}
1;
