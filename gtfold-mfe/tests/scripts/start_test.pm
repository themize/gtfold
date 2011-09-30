#!/usr/bin/perl
use strict;
use warnings;
use Module::Load;

use File::Basename;
use File::Path;
use Log::Log4perl qw(:easy);

require 'test_utils.pl';

my $configdir="../config";
Log::Log4perl::init( "$configdir/root-logger.conf" );
my $logger = Log::Log4perl->get_logger;

$logger->info("Starting Tests...");

# Read Parameter File
my $paramfile = "$configdir/test-params.conf";

# Create Work Directory
my $workdir = "../work";
$logger->info("Deleting directory $workdir");
rmtree($workdir, 0, 1);
$logger->info("Creating work directory $workdir");
mkdir $workdir;

my %Config;

%Config = load_config_file($paramfile);

###### Prepare a Hashmap of Sequences ######

my %Sequences;
my @seqdir_arr = @{$Config{"G_SEQUENCE_DIR"}};

my $seq_include_regex = $Config{"G_INCLUDE_SEQUENCES"};
my $seq_exclude_regex = $Config{"G_EXCLUDE_SEQUENCES"};

foreach (@seqdir_arr) {

  opendir(DIR, $_) || die $!;
  while (my $seqfile = readdir(DIR)) {

	  my $seqname;
	  my $path;
	  my $suffix;
	  ($seqname,$path,$suffix) = fileparse($seqfile);

	  if (-d $seqfile) {
		  next;
	  }

    $seqfile = "$_$seqname";
    my $seq_include = (not defined($seq_include_regex)) || ($seqname =~ /$seq_include_regex/);

    my $seq_exclude = (not defined($seq_exclude_regex)) || ($seqname !~ /$seq_exclude_regex/);
    if ( $seq_include && $seq_exclude ) {
      $Sequences{$seqname} = $seqfile;
      $logger->info("Selected Sequence ... $seqname");
    }
  }
}

my $test_include_regex = $Config{"G_INCLUDE_TESTS"};
my $test_exclude_regex = $Config{"G_EXCLUDE_TESTS"};

my $test_list_file = $Config{"G_TEST_LIST_FILE"};

open(TESTLISTFILE, $test_list_file) || die("Could not open file: $test_list_file");

  while (<TESTLISTFILE>) {

      my $testname = $_;
      chomp($testname);

      $testname =~ s/^\s*//;     # Remove spaces at the start of the line
      $testname =~ s/\s*$//;     # Remove spaces at the end of the line
      if ( ($testname !~ /^#/) && ($testname ne "") ) {    # Ignore lines starting with # and blank lines

      my $test_include = (not defined($test_include_regex)) || ($testname =~ /$seq_include_regex/);

      my $test_exclude = (not defined($test_exclude_regex)) || ($testname !~ /$test_exclude_regex/);

      if ( $test_include && $test_exclude ) {
        my $module = $testname;
        load($module);
        $module->test(\%Config, \%Sequences, $logger);
      }
    
    }

  }

