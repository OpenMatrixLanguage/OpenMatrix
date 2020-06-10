#!/usr/bin/perl 



use strict;

use warnings;





my @files = glob("*.hml");



my $hml_home = $ENV{HML_HOME};

my $os = $^O;

my $hml = $hml_home;

use File::Compare;



if ($os =~ /darwin/i) {

	$hml = $hml . 'macos' . '/hml2 '

} elsif ($os =~ /linux/i) {

 	$hml = $hml . $os . '/hml2 ';

} else {

	$hml = $hml . '/stCompose.exe ';

}



foreach my $hmlfile (@files) {

	print "Running $hml $hmlfile ... :: ";

	my $file = substr($hmlfile, 0, index($hmlfile , '.'));

	my $log = $file . '.log';

	my $cmd =  $hml . $hmlfile;

	if ($os =~ /linux/i) {

		print $cmd;

		$cmd = $cmd . ' >& ' . $log ;

	} elsif ($os =~ /darwin/i) {

		$cmd = $cmd . ' >& ' . $log ;

	} else {

		$cmd = $cmd . ' >>  ' . $log . '2>&1';

	}


	system ($cmd); 

	my $reflog = 'RefLogs/' . $file . '.log'; 

	if (compare($log, $reflog) == 0) {

        	print "$hmlfile PASSED\n";

		unlink $log;

	} else {

		print "$hmlfile FAILED\n";

	}

}



exit 0;

