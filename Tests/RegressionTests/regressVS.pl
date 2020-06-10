#!/usr/bin/perl 

use strict;
use warnings;

# capture the time the script started, used for displaying elapsed time
my $t0;
BEGIN
{
    $t0 = time;
}

use File::Copy;
use File::Path;
use File::Glob qw(:globally :nocase);

my @files = glob("*.oml");
my @DIR = glob ("*/*.oml");
my @logFiles = glob ("*.log");
unlink @logFiles;
my@filtered;
my @array = ('', grep -d, glob '*');

my $argsucess = 0;
my $index = 0;
my $TestSpace = 'TestSpace';
#make the test directory if it doesn't already exist
mkdir $TestSpace unless -d $TestSpace;

foreach my $el (@ARGV){
    
    if (index($ARGV[$index], "-include") != -1){
        
       
        my @incDirs = split(/ /, $ARGV[$index]);
        #reset the files directory since I don't want to use it
        undef(@files);
        undef(@DIR);      
        foreach my $dirToInclude (@incDirs){
           
            push (@DIR, glob($dirToInclude . '/' . '*.oml'));
        }
       
        
    }
    else{
         foreach my $element (@DIR){
            @filtered = grep !m{^$el/}, @DIR;        
            $argsucess++;      
        }
    }
   
    $index++;
}
my $filteredLen = scalar @filtered;
if ($filteredLen == 0)
{
    push(@filtered, @DIR);
}

if ($argsucess <= 0 && $#ARGV + 1 > 0 && index($ARGV[0], "-include") == -1){
    print "one or more arguments was not found check spelling\n";
}


my $os = $^O;
use Cwd;
my $oml = cwd() . '/../stCompose/x64/Debug';
#my $oml = cwd() . '/../stCompose/x64/Release';

use File::Compare;
if ($os =~ /darwin/i) {
    $oml = $oml . 'macos' . '/stCompose ';
} elsif ($os =~ /linux/i) {
    $oml = $oml . $os . '/stCompose ';
} else {
    $oml = $oml . '/stCompose.exe ';
}

my $int = 0;
my $goodtests = 0;
my @failedTests = qw();

#just run it regularly
foreach my $omlfile (@files) {
    print "Running $omlfile ... :: ";
    my $file = substr($omlfile, 0, index($omlfile , '.'));
    my $log = $file . '.log';
    my $cmd =  $oml . $omlfile;
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
            $goodtests++;
            print "$omlfile PASSED\n";
        unlink $log;
    } else {
        $int++;
        push(@failedTests, "$omlfile");
        
        print "$omlfile *********** FAILED ***********\n";
    }
   
}

#if it's in a directory we have to do something special (a couple things actually)
#tried some stuff 
foreach my $omlfile (@filtered) {
    #copy and move the file into the parent directory
    copy ($omlfile, 'TestSpace');
    #now that it's moved we have to strip out the directory stuff before the file name in the string
    $omlfile =~ s:.*/([^/]+):$1:;
    #now we can run it as if it was in the regular directory
    
    print "Running $omlfile ... :: ";
    #$omlfile =~ s:.*/([^/]+):$1:;
    
    my $file = substr($omlfile, 0, index($omlfile , '.'));
    my $log = $file . '.log';
    my $cmd =  $oml . 'TestSpace/' . $omlfile;
           

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
            $goodtests++;
            print "$omlfile PASSED\n";
        unlink $log;
    } else {
        $int++;
        push(@failedTests, "$omlfile");
       
        print "$omlfile *********** FAILED ***********\n";
    }
    #delete file 
    unlink 'TestSpace/' . $omlfile;
}
print "$int tests failed\n";

if (scalar(@failedTests) == 0) {
    print "All tests passed, $goodtests were performed\n";
    
}
else{
    print "Number of tests that passed: $goodtests\n";
    print "The following tests failed: \n";
    foreach (@failedTests){
        print $_, "\n";
    }
}
print "Program: $oml\n";

# capture the time the script completed.  Calculate and format elapsed time.
    my $d = time() - $t0;
    my @int = (
        [ 'second', 1                ],
        [ 'minute', 60               ],
        [ 'hour',   60*60            ],
        [ 'day',    60*60*24         ],
        [ 'week',   60*60*24*7       ],
        [ 'month',  60*60*24*30.5    ],
        [ 'year',   60*60*24*30.5*12 ]
    );
    my $i = $#int;
    my @r;
    while ( ($i>=0) && ($d) )
    {
        if ($d / $int[$i] -> [1] >= 1)
        {
            push @r, sprintf "%d %s%s",
                         $d / $int[$i] -> [1],
                         $int[$i]->[0],
                         ( sprintf "%d", $d / $int[$i] -> [1] ) > 1
                             ? 's'
                             : '';
        }
        $d %= $int[$i] -> [1];
        $i--;
    }

    my $runtime = join ", ", @r if @r;
    # warn sprintf "RUNTIME %s\n", $runtime;
    # print elapsed time
    printf "RUNTIME %s\n", $runtime;
    

exit 0;
