use strict;
use warnings;
use File::Copy;
use File::Path;
use File::Glob qw(:globally :nocase);

my @logFiles = glob ("*.log");
my @refLogFiles = glob ("RefLogs/*.log");

foreach my $log (@ARGV)
{
	system("notepad.exe $log");
		
}
exit 0;