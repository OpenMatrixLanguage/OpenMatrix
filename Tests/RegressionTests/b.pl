#!/usr/bin/perl 

use strict;
use warnings;


my @files = glob("*.hml");

my $hml_home = $ENV{HML_HOME};
my $os = $^O;
my $hml = $hml_home . '/VisSimMATH/x64/Release';
use File::Compare;

if ($os =~ /darwin/i) {
	$hml = $hml . 'macos' . '/hml2 ';
} elsif ($os =~ /linux/i) {
 	$hml = $hml . $os . '/hml2 ';
} else {
	$hml = $hml . '/VisSimMATH.exe ';
}
my $int = 0;
my @failedTests = qw();
foreach my $hmlfile (@files) {
	print "Running $hmlfile ... :: ";
	my $file = substr($hmlfile, 0, index($hmlfile , '.'));
	my $log = $file . '.log';
	my $cmd =  $hml . $hmlfile;
	
	if ($os =~ /linux/i) {
		print $cmd;
		$cmd = $cmd . ' >& ' . $log ;
	} elsif ($os =~ /darwin/i) {
		$cmd = $cmd . ' >& ' . $log ;
	} else {
		$cmd = $cmd . ' >  ' . $log . ' 2>&1';
	}
	system ($cmd); 
	my $reflog = 'RefLogs/' . $file . '.log'; 
	if (compare($log, $reflog) == 0) {
        	print "$hmlfile PASSED\n";
		unlink $log;
	} else {
		$int++;
		push(@failedTests, "$hmlfile");
		print "$hmlfile *********** FAILED ***********\n";
	}
	
}
print "$int tests failed\n";

if (scalar(@failedTests) == 0) {
	print "All tests passed";
	
}
else{
	print "The following tests failed: \n";
	foreach (@failedTests){
		print $_, "\n";
	}
}

exit 0;
