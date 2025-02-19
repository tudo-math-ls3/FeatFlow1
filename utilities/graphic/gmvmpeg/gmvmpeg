#!/usr/bin/env perl
#
# by Sven H.M. Buijssen

#
# gmvmpeg - script to automate generation of mpeg movies with GMV.
#
# History:
# 2005/10/14   v3.0.3  - Allow paths that start with an hash (#)
# 2005/10/11   v3.0.2  - Removed the deprecated MPEG-1 encoder. Replaced it with
#                        settings for the MPEG-2 encoder to create MPEG-1 format.
#                      - Create movie even without write permission in current 
#                        working directory
# 2005/05/30   v3.0.1  - Switched to Getopt::Long.
# 2005/05/23   v3.0    - Introduction of parallel mode
#                      - Introduction of flexible pattern to match files against
#                      - Dynamically determine the lines where the helper 
#                        programs are set for the error messages when one of 
#                        these programs can't be found.
#                      - Added failback mechanism: 
#                        Try to determine the path of some helper programs when
#                        the hard-coded path fails
# 2005/05/17   v2.5.5  - Release version for FeatFlow CD 1.3
# 2005/01/14   v2.5.4  - Check whether horizontal and vertical resolution are even.
#                      - New option '--checkonly'
# 2004/05/28   v2.5.3  - When decompressing on-the-fly, show file name before 
#                        decompressing
#                      - Do not hard-code location of perl in first line, but use 
#                        perl from $PATH
# 2004/05/28   v2.5.2  - Fixed bugs concerning automatic file deletion and MPEG-2 
#                        bit rate setting
# 2004/05/28   v2.5.1  - Added support for MPEG-2 movies with bigger resolutions 
#                        than 720x576
# 2004/01/28   v2.5.0  - Added support for MPEG-2 encoding
# ...
# ...                    (Undocumented changes)
# ...
# 2002/../..   v1.?.?  - Added decompressing on-the-fly
# 2002/../..   v1.?.?  - Changes to script to reflect the changes in GMV
#                        (availability of standalone gmvbatch, no snapshots
#                         with gmv itself possible any more)
# 2001/09/01   v1.0.0  - First release
#

# Allow proper programming only
# (a bit like Fortran90's "implicit none")
use strict;
use warnings;

# being able to use tmpname to create temporary files
use POSIX;

# handling of command line options
use Getopt::Long qw(:config no_ignore_case gnu_getopt no_bundling);

# Portably parse a pathname into directory, basename and suffixes
use File::Basename;

# Supply object methods for I/O handles
use IO::Handle;

# Methods for manipulating file specifications.
use File::Spec;

# Enable Parallel invocations of GMV on machines with multiple cpus
my $enableParallelMode = 0;
if (eval "require Parallel::ForkManager;") {
    $enableParallelMode = 1;
}


#
# Location of programs needed
#
my $gmvGL        ="/usr/local/bin/gmvgl";          # OpenGL version of gmv
my $gmvMesa      ="/usr/local/bin/gmv";            # Mesa version of gmv
my $gmvBatch     ="/usr/local/bin/gmvbatch";       # Mesa Batch version of gmv
my $gmvToUse     =$gmvMesa;

my $sgitopnm     ="/usr/local/bin/sgitopnm";       # Convert gmv rgb screenshots to NetPBM PPM format
my $convert      ="/usr/local/bin/convert";        # Convert gmv rgb screenshots to NetPBM PPM format
my $mpegEncoder  ="/usr/local/bin/mpeg2encode";    # MPEG-1/2 encoder (google for mpeg2vidcodec_v12.tar.gz)

# Some utilities
my $uncompress   ="/usr/bin/uncompress";
my $gunzip       ="/usr/bin/gunzip";
my $bunzip2      ="/usr/bin/bunzip2";


#
# Some default values
#
(my $progname=$0) =~ s/^.*\/(.*)/$1/;
my $version       = "3.0.3";
my $author        = "sven.buijssen\@uni-dortmund.de";
my $wildcard      = "%";
my $pathdelimiter = "/";


# Catch ^C
$SIG{INT} = \&cleanExit;
STDOUT->autoflush(1);


# Infrastructure for parallel invocations of GMV
my $pm;
if ($enableParallelMode) {
    $pm = new Parallel::ForkManager(1);
    # Setup a callback for when a child finishes up so we can
    # get it's exit code
    $pm->run_on_finish(
        sub {
           my ($pid) = @_;
           print "*** Process $pid has finished.\n";
        }
    );

    $pm->run_on_start(
        sub {
	    my ($pid, $ident) = @_;
	    print "*** Waiting for forked process $pid to finish.\n";
        }
    );

    $pm->run_on_wait(
        sub {
	    print ".";
        }, 0.5
    );
}



my %options;
&parseargv2();

# Debug
# foreach my $entry (sort keys %options) {
#     print $entry . " " x (30 - length($entry)) . $options{$entry} . "\n";
# }

# Help requested?
if ($options{"help"}) {
    &usage();
    exit 0;
}

# Version requested?
die "Version: $progname $version\nCopyright: Copyright (C) 2001-2005, $author\n\n" if ($options{"version"});

# Greeting
print "$progname version $version\n";

# Check if all programs needed are available
&check_file_existence();
exit if ($options{"checkonly"});


# Set number of parallel processes
if ($enableParallelMode) {
    if ($options{"jobs"} > 1) {
	$pm->set_max_procs($options{"jobs"});
    } else {
	$enableParallelMode = 0;
    }
} else {
    if ($options{"jobs"} > 1) {
        warn "\n$progname: WARNING: perl module Parallel::ForkManager not found.\n" .
	     "Get it from www.cpan.org.\n" .
	     "Parallel mode disabled!\n\n";
    }
}



# Determine all gmv files to include in movie
my $numOfpadWithZeros = 0;
if ($options{"inputfilePattern"} =~ m/($wildcard+)/) {
    $numOfpadWithZeros = length($1);
} else {
    $options{"inputfilePattern"} =~ s/$options{"workingDirectory"}$pathdelimiter//;
    die "\n$progname: Invalid pattern <" . $options{"inputfilePattern"} . ">\n" .
	"specified. No wildcard character $wildcard within found.\n";
}

my @fileList = ();
print "*** Creating file list of gmv files to visualise... ";
print "\n" if ($enableParallelMode || $options{"verbose"} ne "");
for (my $i = $options{"first"};  $i <= $options{"last"};  $i += $options{"stride"}) {
    my $filename = $options{"inputfilePattern"};
    $filename =~ s/($wildcard+)/"0" x ($numOfpadWithZeros - length($i)) . $i/e;

    # Test existence of files
    if ($options{"verbose"} eq "") {
	print "Testing existence of $filename";
	print "["    if ($options{"useCompress"} || $options{"useGzip"} || $options{"useBzip2"});
	print ".Z"   if ($options{"useCompress"});
	print ","    if ($options{"useCompress"} && ($options{"useGzip"} || $options{"useBzip2"}));
	print ".gz"  if ($options{"useGzip"});
	print ","    if ($options{"useGzip"} && $options{"useBzip2"});
	print ".bz2" if ($options{"useBzip2"});
	print "]"    if ($options{"useCompress"} || $options{"useGzip"} || $options{"useBzip2"});
	print "\n";
    }

    if (-r $filename) {
        push @fileList, $filename;

    } elsif ($options{"useCompress"}  &&  -r $filename . ".Z") {
        push @fileList, $filename . ".Z";

    } elsif ($options{"useGzip"}  &&  -r $filename . ".gz") {
        push @fileList, $filename . ".gz";

    } elsif ($options{"useBzip2"}  &&  -r $filename . ".bz2") {
        push @fileList, $filename . ".bz2";
    }

    # If file does not exist and we don't have an adaptive time stepping
    # print error message and exit.
    elsif (! $options{"adaptive"}) {
        die "$progname: $filename: No such file or directory\n";
    }
}

if ($#fileList == -1) {
    die "$progname: No input files found.\n\n" .
	"Use the command line option '--verbose' to debug the file search.\n";
}


my $j=0;
#(my $outbasename = $options{"outputfile"}) =~ s/^(.*)\.[^\.]+$/$1/;
my $outbasename = $options{"outputfile"};
foreach my $file (@fileList) {
    if (! -e $file) {
	warn "$progname: WARNING:\n<$file> has disappeared. Ignoring.\n";
	next;
    }

    $j++;
    my $pid = $pm->start($file) and next if ($enableParallelMode);
    my $tmpname = "";
    my $command;

    #
    # START
    # Code for child process
    #
    my $isTempfile = 0;
    if ($enableParallelMode || $options{"verbose"} eq "") {
	print "*** Processing $file...\n";
    } else {
	print "*** Processing $file...";
    }

    # Decompression needed?
    #
    # compressed?
    if ($file =~ m/\.Z$/) {
	print "      Uncompressing $file...\n" if ($options{"verbose"} eq "");
        $tmpname = POSIX::tmpnam();
	$command = "$uncompress -c $file > $tmpname";
	$command =~ s/\#/\\\#/g;    # Escape all hashs 

        print "\n$command" if ($options{"verbose"} eq "");
        system("$command");
	if ($?) {
	    die "\n\n$progname: $gmvToUse has terminated unexpectedly.\n" .
		"Decompression on-the-fly using <$uncompress> failed: $!\n";
	}
        $file = $tmpname;
	$isTempfile = 1;
    }

    # gzipped?
    elsif ($file =~ m/\.gz$/) {
	print "      Uncompressing $file...\n" if ($options{"verbose"} eq "");
        $tmpname = POSIX::tmpnam();
	$command = "$gunzip -c $file > $tmpname";
	$command =~ s/\#/\\\#/g;    # Escape all hashs 

        print "\n$command" if ($options{"verbose"} eq "");
        system("$command");
	if ($?) {
	    die "\n\n$progname: $gmvToUse has terminated unexpectedly.\n" .
		"Decompression on-the-fly using <$gunzip> failed: $!\n";
	}
        $file = $tmpname;
	$isTempfile = 1;
    }

    # bzip2-ed?
    elsif ($file =~ m/\.bz2$/) {
	print "      Uncompressing $file...\n" if ($options{"verbose"} eq "");
        $tmpname = POSIX::tmpnam();
	$command = "$bunzip2 -c $file > $tmpname";
	$command =~ s/\#/\\\#/g;    # Escape all hashs 

        print "\n$command" if ($options{"verbose"} eq "");
        system("$command");
	if ($?) {
	    die "\n\n$progname: $gmvToUse has terminated unexpectedly.\n" .
		"Decompression on-the-fly using <$bunzip2> failed: $!\n";
	}
        $file = $tmpname;
	$isTempfile = 1;
    }

    &gmv2rgb($gmvToUse, $options{"attributefile"},
	     $options{"resX"}, $options{"resY"},
	     $file, $outbasename, $j, $options{"verbose"});

    # Delete on-the-fly decompressed files.
    if ($isTempfile) {
	print "Removing $file\n" if ($options{"verbose"} eq "");
	unlink "$file";
    }

    print "\n" if ($options{"verbose"} eq "");

    # pass an exit code to finish
    $pm->finish($?) if ($enableParallelMode);

    #
    # END
    # Code for child process
    #

    print " done.\n" unless ($enableParallelMode);
}

$pm->wait_all_children if ($enableParallelMode);

&generate_movie($options{"resX"}, $options{"resY"},
		$options{"outputfile"}, $outbasename, $j, $options{"verbose"});
exit 0;


sub gmv2rgb {
    my ($gmvToUse, $attrfile, $xres, $yres, $filename, $outfile, $j, $verbose) = (@_);

    my $command = "$gmvToUse -m -a $attrfile -w 0 0 $xres $yres -i $filename ".
                 "-s $outfile$j.rgb";
    # Escape all hashs 
    # (especially FeatFlow users might specify attribute or input files with
    # a relative path starting with # - which would be interpreted by /bin/sh as
    # comment! A command 'gmv -a #some comment' simply crashes, however!)
    $command =~ s/\#/\\\#/g;

    print $command . "\n" if ($verbose eq "");
    system("$command $verbose")
	 && die "\n\n$progname: $gmvToUse has terminated unexpectedly.\n" .
	    "There is no image data that can be converted into a movie.\n\n" .
	    "A few guesses (order according to decreasing probability:)\n" .
	    "* You have pressed Ctrl-C.\n" .
	    "* $gmvToUse crashed. Possibly because of an attribute file\n" .
	    "  that was created with a newer version of gmv than the one\n" .
	    "  used in $progname. Check version numbers!\n";

    # Convert gmv output to mpeg input format
    my $convertCall = "";
    my $convertProg = "";
    if ( -X $sgitopnm ) {
	$convertCall = "$sgitopnm $outfile$j.rgb > $outfile$j.ppm";
	$convertProg = $sgitopnm;
    } elsif ( -X $convert ) {
	$convertCall = "$convert $outfile$j.rgb $outfile$j.ppm";
	$convertProg = $convert;
    } else {
	die "$progname: Neither one of the convert programs could be found.";
    }
    # Escape all hashs 
    $convertCall =~ s/\#/\\\#/g;

    if ($verbose eq "") {
	print  "$convertCall\n";
	system("$convertCall")
           && die "\n\n$progname: $convertProg has reported an error.\n" .
                  "Use option --verbose to scrutinise the problem.\n\n";
	print "File triple <$outfile$j.ppm> has been created.\n";
	print "*** ";
    } else {
	system("$convertCall 2>/dev/null")
           && die "\n\n$progname: $convertProg has reported an error.\n" .
                  "Use option --verbose to scrutinise the problem.\n\n";
    }

    return;
}


# Function which turns a sequence of $num files into a MPEG-1/2 movie
sub generate_movie {
    my ($xres, $yres, $outfile, $outbasename, $num, $verbose) = (@_);

    my $tmpname = POSIX::tmpnam();
    my $movieFilename;
    my $statFilename = $tmpname . "_stat.out";

    # Generate mpeg
    if ($options{"mpeg1"}) {
	# MPEG-1 generation
        print "*** Creating parameter file for MPEG-1 encoder... ";
        &create_mpeg_config_file(1, $tmpname, $statFilename, $outfile . "%d", $num, $xres, $yres, $options{"bitrate"});
        $movieFilename = "$outfile.mpeg";

    } else {
	# MPEG-2 generation
        print "*** Creating parameter file for MPEG-2 encoder... ";
        &create_mpeg_config_file(2, $tmpname, $statFilename, $outfile . "%d", $num, $xres, $yres, $options{"bitrate"});
        $movieFilename = "$outfile.m2v";
    }

    print "*** \n*** Parameter file <$tmpname> created.\n" if ($verbose eq "");
    print "done.\n" unless ($verbose eq "");
    print "*** Generating movie $movieFilename ... ";
    my $command = "$mpegEncoder $tmpname $movieFilename";
    # Escape all hashs 
    $command =~ s/\#/\\\#/g;

    if ($verbose eq "") {
        print "\n$command\n";
        system("$command")
            && die "\n\n$progname: $mpegEncoder has reported an error.\n" .
                   "Use option --verbose to scrutinise the problem.\n\n";
        print "*** ";
    } else {
        system("$command 2>/dev/null")
            && die "\n\n$progname: $mpegEncoder has reported an error.\n" .
                   "Use option --verbose to scrutinise the problem.\n\n";
    }
    print "done.\n";

    # Clean up working directory
    if (! $options{"keepfiles"}) {
        print "*** Removing temporary files... ";
        for (my $jj = 0;  $jj <= $num;  $jj++) {
            unlink "$outbasename$jj.rgb", "$outbasename$jj.ppm";
        }
	unlink "$tmpname", $statFilename;
        print "done.\n";
    }

#    print "\nNote: You might consider creating movies in MPEG-2 format which\n" .
#	  "are usually smaller and of better quality.\n" .
#	  "Just use the command line flag '--mpeg2'.\n" if ($options{"mpeg1"});

    print "\nNote: You might consider using the 'invisible mode' for the rendering\n" .
	  "process. There would be no more GMV windows popping up all the time!\n" .
	  "Just use the command line flag '-I' or '--invisible'.\n" if (! $options{"invisible"});

    return;
}


sub usage {
    print <<EOF
Usage:
$progname [ options ]

    -a <file>, --attributes <file>:
           path for gmv attribute file to use
           (default: default.attr)
    -b <number>, --bitrate <number>:
           Specifies bitrate for MPEG movie.
           (default: 5000000 for MPEG-1, 3500000 for MPEG-2)
    --fls <number>,<number>[,<number>]:
           two or three comma separated digits that specify first and
           last index of gmv input file as well as the stride
           If stride is omitted or set to 0, time step is adaptive (default).
    -h, --help:
           this help screen
    -i <pattern>, --input <pattern>:
           pattern for gmv input files, use $wildcard as wildcard character.
           Several $wildcard mean padding with zeros.
           e.g. u.$wildcard.gmv   for u.2.gmv,   u.3.gmv,   ..., u.100.gmv, ...
                u.$wildcard$wildcard.gmv  for u.02.gmv,  u.03.gmv,  ..., u.100.gmv, ...
                u.$wildcard$wildcard$wildcard.gmv for u.002.gmv, u.003.gmv, ..., u.100.gmv, ...
           (default: u.$wildcard.gmv)
    -j <jobs>, --jobs <jobs>:
           Specifies the number of jobs (commands) to run simultaneously.
           If there is more than one -j option, the last one is effective.
           (default: 1)
    -I, --invisible:
           use gmvBatch for rendering process (invisible mode)
    -k, --keep-files:
           don't delete temporary ppm files
    -m, --max:
           Obsolete. Only available for compatibility reasons
           Has been used to specify maximum size of mpeg movie.
           Use -b/--bitrate instead now.
    --mpeg1:
           Use MPEG-1 encoder (default).
    --mpeg2:
           Use MPEG-2 encoder instead of MPEG-1.
    -o <string>, --output <string>:
           name of output file
           (default: movie.mpeg)
    -V, --verbose:
           verbose mode
    --version:
           print version information
    --checkonly:
           verify that all program paths are correct and exit
    --wd <string>:
           set working directory for attribute, input and output files
           (only applied when not given with absolute path)
    -x <number>:
           horizontal resolution (default: 800 for MPEG-1, 720 for MPEG-2)
    -y <number>:
           vertical resolution (default: 600 for MPEG-1, 576 for MPEG-2)
    -Z:
           decompress input files on-the-fly using compress
    -z, -gzip:
           decompress input files on-the-fly using gzip
    --bzip2:
           decompress input files on-the-fly using bzip2

Example:
   $progname -a gmv_example.attr -o example -i u.$wildcard.gmv --fls 1,99,2

EOF
;
   #'
}


sub parseargv2 {
    GetOptions ("attributes|a=s" => \$options{"attributefile"},
		"bitrate|b=s"   => \$options{"bitrate"},
		"bzip2"         => \$options{"useBzip2"},
		"checkonly"     => \$options{"checkonly"},
		"fls=s"         => \$options{"firstLastStride"},
		"gzip|z"        => \$options{"useGzip"},
		"help|h"        => \$options{"help"},
		"input|i=s"     => \$options{"inputfilePattern"},
		"invisible|I"   => \$options{"invisible"},
		"jobs|j=s"      => \$options{"jobs"},
		"keep-files|k"  => \$options{"keepfiles"},
		"max|m=s"       => \$options{"maxfilesize"},
		"mpeg1"         => \$options{"mpeg1"},
		"mpeg2"         => \$options{"mpeg2"},
		"output|o=s"    => \$options{"outputfile"},
		"verbose|V"     => \$options{"verbose"},
		"version|v"     => \$options{"version"},
		"wd=s"          => \$options{"workingDirectory"},
		"x=s"           => \$options{"resX"},
		"y=s"           => \$options{"resY"},
		"Z"             => \$options{"useCompress"}
    ) ||
	die $progname.": ".
	    "Unknown option found.\n" .
	    (" " x (length($progname)+2)) .
	    "Check your options. Use -h for help.\n";

    # Set default values
    $options{"attributefile"}    ||= "default.attr";
    $options{"checkonly"}        ||= 0;
    $options{"firstLastStride"}  ||= "1,100000,1";
    $options{"help"}             ||= 0;
    $options{"inputfilePattern"} ||= "u.$wildcard.gmv";
    $options{"invisible"}        ||= 0;
    $options{"jobs"}             ||= 1;
    $options{"keepfiles"}        ||= 0;
    $options{"maxfilesize"}      ||= 0;
    $options{"outputfile"}       ||= "movie.mpeg";
    $options{"useBzip2"}         ||= 0;
    $options{"useCompress"}      ||= 0;
    $options{"useGzip"}          ||= 0;
    $options{"version"}          ||= 0;
    chomp(my $pwd=`pwd`);
    $options{"workingDirectory"} ||= $pwd;

    $options{"verbose"} ||= "> /dev/null";
    if ($options{"verbose"} =~ m/^1$/) {
	$options{"verbose"} = "";
    }

    if (! defined($options{"mpeg1"})  &&
	! defined($options{"mpeg2"})) {
	$options{"mpeg1"} = 1;
	$options{"mpeg2"} = 0;
    }

    if ($options{"mpeg1"}) {
	$options{"bitrate"} ||= 5000000;
	$options{"resX"}    ||= 800;
	$options{"resY"}    ||= 600;

	if (defined($options{"mpeg2"})  &&  $options{"mpeg2"}) {
	    die "$progname: ERROR in command line arguments:\n" .
		"You cannot state --mpeg1 and --mpeg2 at the same time.\n\n";
	}
	$options{"mpeg2"} = 0;
    }
    elsif ($options{"mpeg2"}) {
	$options{"bitrate"} ||= 3500000;
	$options{"resX"}    ||= 720;
	$options{"resY"}    ||= 576;
	$options{"mpeg1"}     = 0;
    }

    # GMV file sequence:
    # adaptive stepping? valid numbers given?
    $options{"adaptive"} = 0;
    my @list = split(',', $options{"firstLastStride"});
    if ($#list < 1  ||  $#list > 3) {
	die "$progname: ERROR in command line argument -fls\n" .
	    "Syntax has to be:\n   first,last,stride\n\n";
    }
    ($options{"first"}, $options{"last"}, $options{"stride"}) = @list;
    if (! defined($options{"stride"})  ||  $options{"stride"} == 0 ) {
	$options{"adaptive"} = 1;
	$options{"stride"}   = 1;
    }


    # Check if the x- and y-resolution are conforming to MPEG-1 and -2 standard.
    # Otherwise
    # * the MPEG-1 encoder will produce a 19bit sized movie and
    # * the MPEG-2 encoder will die with an error message
    # So, in any case, all the gmv visualisation has been in vain.
    if ($options{"resX"} % 2) {
	die "$progname: ERROR in command line argument -x\n" .
	    "Horizontal size must be an even (4:2:0 / 4:2:2) for MPEG-1 and -2 movies.\n";
    } elsif ($options{"resY"} % 2) {
	die "$progname: ERROR in command line argument -y\n" .
	    "Vertical size must be an even (4:2:0 / 4:2:2) for MPEG-1 and -2 movies.\n";
    }

    # Number of processes reasonable?
    if ($options{"jobs"} < 1  ||  $options{"jobs"} > 8) {
	die "$progname: ERROR in command line argument -j\n" .
	    "Please specify a reasonable value for the number of gmv processes\n" .
	    "to invoke simultaneously.\n\n";
    }

    # Invisible mode?
    $gmvToUse = $gmvBatch if ($options{"invisible"});

    # Make paths absolute
    foreach my $entry ($options{"attributefile"},
		       $options{"inputfilePattern"},
		       $options{"outputfile"}) {
	# Prepend given working directory (current directory by default)
	# if no directory information is contained in $entry
	if (&dirname($entry) !~ m/^$pathdelimiter/  &&
	             $entry  !~ m/$pathdelimiter/) {
	    $entry = $options{"workingDirectory"} . $pathdelimiter . $entry;
	}
    }

    return;
}


# sub parseargv {
# # I should really convert this routine to use Getopt::Long!!
#     chomp(my $pwd=`pwd`);

#     my ($wd) = ($pwd);
#     my ($attrfile, $inbasename, $outfile) = ("default.attr", "u.$wildcard.gmv", "movie");
#     my ($first, $last, $stride, $adaptive) = (1, 100000, 1, 1);
#     my ($invisibleMode, $use_mpeg2, $options{"useCompress"}, $use_gunzip, $use_bunzip2) = (0, 0, 0, 0, 0);
#     my ($keepfiles, $maxmb, $max_bits, $bitrate) = (0, 0, 0, 0);
#     my ($procs) = (1);
#     my ($xres, $yres) = (0, 0);
#     my $verbose = "> /dev/null";
#     my $checkonly = 0;
#     my @list;

# ARGS: for (my $i=0; $i<=$#ARGV; $i++) {
#         my $arg=$ARGV[$i];

#         if ($arg eq "-a") {
#             $attrfile = $ARGV[$i+1];
#         } elsif ($arg eq "-fls") {
#             @list = split(/,/, $ARGV[$i+1]);
#             if ($#list < 1  ||  $#list > 3) {
#                 die "$progname: Syntax error in argument -fls, specifying\n".
#                     "first and last index as well as stride.\n\n".
#                     "Syntax has to be:\n   first,last,stride\n\n";
#             }

#             $first=$list[0];
#             $last=$list[1];
#             if ($#list==2) {
#                 $stride=$list[2];
#             } else {
#                 $stride=0;
#             }

#             # adaptive time stepping if stride=0
#             if ($stride==0) {
#                 $adaptive=1;
#                 $stride=1;
#             }

# 	    $i++;
#         } elsif ($arg eq "-h" || $arg eq "--help") {
#             usage();
#             exit 0;
#         } elsif ($arg eq "-i") {
#             $inbasename=$ARGV[$i+1];
#         } elsif ($arg eq "-I" || $arg eq "--invisible") {
#             $invisibleMode=1;
#             $gmvToUse=$gmvBatch;
#         } elsif ($arg eq "-k" || $arg eq "--keep-files") {
#             $keepfiles=1;
#         } elsif ($arg eq "-m" || $arg eq "--max") {
# #           $maxmb=$ARGV[$i+1];
# 	    $i++;
#         } elsif ($arg eq "-b" || $arg eq "--bitrate") {
# 	    $bitrate=$ARGV[$i+1];
# 	    $i++;
# 	} elsif ($arg eq "--mpeg1") {
# 	    $use_mpeg2=0;
# 	} elsif ($arg eq "--mpeg2") {
# 	    $use_mpeg2=1;
#         } elsif ($arg eq "-o") {
#             $outfile=$ARGV[$i+1];
# 	    $i++;
#         } elsif ($arg eq "-V" || $arg eq "--verbose") {
#             $verbose="";
#         } elsif ($arg eq "--version") {
#             die "$progname version $version\n" .
#                 "use '$progname -h' for a list of options\n";
#         } elsif ($arg eq "-wd") {
#             $wd=$ARGV[$i+1];
#         } elsif ($arg eq "-j") {
# 	    # parallel visualisation mode is experimental, not working yet
#             $procs=int($ARGV[$i+1]);
#             if ($procs < 1  ||  $procs > 8) {
#                 die "$progname: Error in argument -j\n".
#                     "Please specify a reasonable value for the number of gmv processes\n\n";
#             }
# 	    $i++;
#         } elsif ($arg eq "-x") {
#             $xres=$ARGV[$i+1];
# 	    $i++;
#         } elsif ($arg eq "-y") {
#             $yres=$ARGV[$i+1];
# 	    $i++;
#         } elsif ($arg eq "-Z") {
#             $options{"useCompress"}=1;
#         } elsif ($arg eq "-gz") {
#             $use_gunzip=1;
#         } elsif ($arg eq "-bz2") {
#             $use_bunzip2=1;
#         } elsif ($arg eq "--checkonly") {
#             $checkonly = 1;
#         }
#     }

#     # Check if the x- and y-resolution are conforming to MPEG-1 and -2 standard.
#     # Otherwise
#     # * the MPEG-1 encoder will produce a 19bit sized movie and
#     # * the MPEG-2 encoder will die with an error message
#     # So, in any case, all the gmv visualisation has been in vain.
#     if ($xres % 2) {
# 	die "$progname: Horizontal size must be an even (4:2:0 / 4:2:2) for MPEG-1 and -2 movies.\n";
#     } elsif ($yres % 2) {
# 	die "$progname: Vertical size must be an even (4:2:0 / 4:2:2) for MPEG-1 and -2 movies.\n";
#     }

# # 	if ($xres > 720) {
# # 	    print "$progname: Warning: Horizontal size has been set to 720 to conform to MPEG-2 standard.\n";
# # 	    $xres = 720;
# # 	}
# # 	if ($yres > 576) {
# # 	    print "$progname: Warning: Vertical size has been set to 576 to conform to MPEG-2 standard.\n";
# # 	    $yres = 576;
# # 	}

#     $xres    = 800     if ($xres == 0     &&  ! $use_mpeg2);
#     $yres    = 600     if ($yres == 0     &&  ! $use_mpeg2);
#     $bitrate = 5000000 if ($bitrate == 0  &&  ! $use_mpeg2);
#     $xres    = 720     if ($xres == 0     &&    $use_mpeg2);
#     $yres    = 576     if ($yres == 0     &&    $use_mpeg2);
#     $bitrate = 3500000 if ($bitrate == 0  &&    $use_mpeg2);
#     $attrfile   = "$wd/$attrfile"   unless (dirname($attrfile)   =~ m|^/|);
#     $inbasename = "$wd/$inbasename" unless (dirname($inbasename) =~ m|^/|);
#     $outfile    = "$wd/$outfile"    unless (dirname($outfile)    =~ m|^/|);

#     return ($xres, $yres, $procs, $attrfile, $inbasename, $outfile,
#             $first, $last, $stride, $adaptive,
#             $invisibleMode, $use_mpeg2, $options{"useCompress"}, $use_gunzip, $use_bunzip2,
#             $keepfiles, $maxmb, $max_bits, $bitrate, $verbose, $checkonly);
# }


sub check_file_existence {
    my $line = 0;
    my $failBack = "";

    # Check for gmv program
    if (! -X $gmvToUse) {
	$line = &get_line_no("^my .gmvToUse");
	my $line1 = $line - 4;
	die "$progname: GMV program not found or not executable.\n".
	    "$gmvToUse: No such file or directory/Not executable.\n".
	    "Please check the values in line " . $line1 . "-" . $line . " of $progname\n".
	    "and the permissions of <$gmvToUse>.\n" .
	    "Nothing done.\n";
    }

    # Check for convert programs
    if (! -X $sgitopnm) {
	$line = &get_line_no("^my .sgitopnm");

	# Try to determine program location with built-in which
	$failBack = &which("sgitopnm");
	if (defined($failBack)) {
	    warn "$progname: WARNING:\nOne of the convert programs could not be found.\n".
   	         "$sgitopnm: No such file or directory/Not executable.\n".
		 "Please check the value in line " . $line . " of $progname.\n".
		 "Using <$failBack> instead.\n\n";
	    $sgitopnm = $failBack;
	}
    }
    if (! -X $sgitopnm  &&  ! -X $convert) {
	$line = &get_line_no("^my .convert");
	my $line2 = &get_line_no("^my .sgitopnm");

	# Try to determine program location with built-in which
	$failBack = &which("convert");
	if (defined($failBack)) {
	    warn "$progname: WARNING:\nOne of the convert programs could not be found.\n".
   	         "$convert: No such file or directory/Not executable.\n".
		 "Please check the value in line " . $line . " of $progname.\n".
		 "Using <$failBack> instead.\n\n";
	    $convert = $failBack;
	} else {
	    die "$progname: Neither one of the convert programs could be found.\n".
		"$sgitopnm: No such file or directory/Not executable.\n".
		"$convert: No such file or directory/Not executable.\n".
		"Please check the values in line " . $line2 . " and " . $line . " of $progname.\n".
		"Nothing done.\n";
	}
    }

    # Check for mpeg encoder
    if (! -X $mpegEncoder) {
	$line = &get_line_no("^my .mpegEncoder");

	# Try to determine program location with built-in which
	$failBack = &which("mpeg2encode");
	if (defined($failBack)) {
	    warn "$progname: WARNING:\nOne of the convert programs could not be found.\n".
		"$mpegEncoder: No such file or directory/Not executable\n".
		"Please check the value in line " . $line . " of $progname.\n".
		"Using <$failBack> instead.\n\n";
	    $mpegEncoder = $failBack;
	} else {
	    die "$progname: MPEG encoding program not found or not executable.\n".
		"$mpegEncoder: No such file or directory/Not executable\n".
		"Please check the value in line " . $line . " of $progname.\n".
		"Nothing done.\n";
	}
    }

    # Check whether attrib file exists and is readable
    die "$progname: ERROR in command line arguments:\n" .
	"Attribute file <" . $options{"attributefile"} . "> does not exist\n" .
	"or is not readable. Nothing done.\n" if (! -r $options{"attributefile"});

    # Check for utility programs when we are supposed to use them
    if ($options{"useCompress"}  &&  ! -X $uncompress) {
	$line = &get_line_no("^my .uncompress");

	# Try to determine program location with built-in which
	$failBack = &which("uncompress");
	if (defined($failBack)) {
	    warn "$progname: WARNING:\nOne of the convert programs could not be found.\n".
		"$uncompress: No such file or directory/Not executable\n".
		"Please check the value in line " . $line . " of $progname.\n".
		"Using <$failBack> instead.\n\n";
	    $uncompress = $failBack;
	} else {
	    die "$progname: uncompress program not found or not executable.\n".
		"$uncompress: No such file or directory/Not executable\n".
		"Please check the value in line " . $line . " of $progname.\n".
		"Nothing done.\n";
	}
    }
    if ($options{"useGzip"}  &&  ! -X $gunzip) {
	$line = &get_line_no("^my .gunzip");

	# Try to determine program location with built-in which
	$failBack = &which("gunzip");
	if (defined($failBack)) {
	    warn "$progname: WARNING:\nOne of the convert programs could not be found.\n".
		"$gunzip: No such file or directory/Not executable\n".
		"Please check the value in line " . $line . " of $progname.\n".
		"Using <$failBack> instead.\n\n";
	    $gunzip = $failBack;
	} else {
	    die "$progname: gunzip program not found or not executable.\n".
		"$gunzip: No such file or directory/Not executable\n".
		"Please check the value in line " . $line . " of $progname.\n".
		"Nothing done.\n";
	}
    }
    if ($options{"useBzip2"}  &&  ! -X $bunzip2) {
	$line = &get_line_no("^my .bunzip2");

	# Try to determine program location with built-in which
	$failBack = &which("bunzip2");
	if (defined($failBack)) {
	    warn "$progname: WARNING:\nOne of the convert programs could not be found.\n".
		"$bunzip2: No such file or directory/Not executable\n".
		"Please check the value in line " . $line . " of $progname.\n".
		"Using <$failBack> instead.\n\n";
	    $bunzip2 = $failBack;
	} else {
	    die "$progname: uncompress program not found or not executable.\n".
		"$bunzip2: No such file or directory/Not executable\n".
		"Please check the value in line " . $line . " of $progname.\n".
		"Nothing done.\n";
	}
    }

    return;
}


sub create_mpeg_config_file()
{
    my ($mpegVersion, $filename, $statFilename, $name_source_files, $number_frames, $xres, $yres, $bitrate) = (@_);

    open (PARAMFILE, ">", $filename) or
	die "\n\n$progname: Can not open <$filename>: $!\n" .
	    "This file is needed for invocation of MPEG-2 encoder.\n";

    print PARAMFILE
	qq{parameter file for MPEG-2 movie (PAL) from GMV files - created by gmvmpeg $version\n} .
	qq{$name_source_files /* name of source files */\n} .
	qq{-               /* name of reconstructed images ("-": do not store) */\n} .
	qq{-               /* name of intra quant matrix file     ("-": default matrix) */\n} .
	qq{-               /* name of non intra quant matrix file ("-": default matrix) */\n} .
        qq{$statFilename } . " "x(15 - length($statFilename)) . qq{/* name of statistics file ("-": stdout ) */\n} .
	qq{2               /* input picture file format: 0=*.Y,*.U,*.V, 1=*.yuv, 2=*.ppm */\n} .
        qq{$number_frames } . " "x(15 - length($number_frames)) . qq{/* number of frames */\n} .
	qq{1               /* number of first frame */\n} .
	qq{00:00:00:00     /* timecode of first frame */\n} .
	qq{12              /* N (# of frames in GOP) */\n} .
	qq{3               /* M (I/P frame distance) */\n} .
        ($mpegVersion eq 1 ? 1 : 0) .
        qq{               /* ISO/IEC 11172-2 stream */\n} .
	qq{0               /* 0:frame pictures, 1:field pictures */\n} .
        qq{$xres } . " "x(15 - length($xres)) . qq{/* horizontal_size */\n} .
        qq{$yres } . " "x(15 - length($yres)) . qq{/* vertical_size */\n} .
	qq{1               /* aspect_ratio_information 1=square pel, 2=4:3, 3=16:9, 4=2.11:1 */\n} .
	qq{3               /* frame_rate_code 1=23.976, 2=24, 3=25, 4=29.97, 5=30 frames/sec. */\n} .
        qq{$bitrate } . " "x(15 - length($bitrate)) . qq{/* bit_rate (bits/s) */\n} .
	qq{112             /* vbv_buffer_size (in multiples of 16 kbit) */\n} .
	qq{0               /* low_delay  */\n} .
	qq{0               /* constrained_parameters_flag */\n} .
	qq{4               /* Profile ID: Simple = 5, Main = 4, SNR = 3, Spatial = 2, High = 1 */\n} .
	get_level_id($xres, $yres) .
	qq{               /* Level ID:   Low = 10, Main = 8, High 1440 = 6, High = 4          */\n} .
	qq{0               /* progressive_sequence */\n} .
	qq{1               /* chroma_format: 1=4:2:0, 2=4:2:2, 3=4:4:4 */\n} .
	qq{1               /* video_format: 0=comp., 1=PAL, 2=NTSC, 3=SECAM, 4=MAC, 5=unspec. */\n} .
	qq{5               /* color_primaries */\n} .
	qq{5               /* transfer_characteristics */\n} .
	qq{5               /* matrix_coefficients */\n} .
        qq{$xres } . " "x(15 - length($xres)) . qq{/* display_horizontal_size */\n} .
        qq{$yres } . " "x(15 - length($yres)) . qq{/* display_vertical_size */\n} .
	qq{0               /* intra_dc_precision (0: 8 bit, 1: 9 bit, 2: 10 bit, 3: 11 bit */\n} .
	qq{1               /* top_field_first */\n} .
	qq{0 0 0           /* frame_pred_frame_dct (I P B) */\n} .
	qq{0 0 0           /* concealment_motion_vectors (I P B) */\n} .
	qq{1 1 1           /* q_scale_type  (I P B) */\n} .
	qq{1 0 0           /* intra_vlc_format (I P B)*/\n} .
	qq{0 0 0           /* alternate_scan (I P B) */\n} .
	qq{0               /* repeat_first_field */\n} .
	qq{0               /* progressive_frame */\n} .
	qq{0               /* P distance between complete intra slice refresh */\n} .
	qq{0               /* rate control: r (reaction parameter) */\n} .
	qq{0               /* rate control: avg_act (initial average activity) */\n} .
	qq{0               /* rate control: Xi (initial I frame global complexity measure) */\n} .
	qq{0               /* rate control: Xp (initial P frame global complexity measure) */\n} .
	qq{0               /* rate control: Xb (initial B frame global complexity measure) */\n} .
	qq{0               /* rate control: d0i (initial I frame virtual buffer fullness) */\n} .
	qq{0               /* rate control: d0p (initial P frame virtual buffer fullness) */\n} .
	qq{0               /* rate control: d0b (initial B frame virtual buffer fullness) */\n} .
	qq{2 2 11 11       /* P:  forw_hor_f_code forw_vert_f_code search_width/height */\n} .
	qq{1 1 3  3        /* B1: forw_hor_f_code forw_vert_f_code search_width/height */\n} .
	qq{1 1 7  7        /* B1: back_hor_f_code back_vert_f_code search_width/height */\n} .
	qq{1 1 7  7        /* B2: forw_hor_f_code forw_vert_f_code search_width/height */\n} .
	qq{1 1 3  3        /* B2: back_hor_f_code back_vert_f_code search_width/height */\n};

    close (PARAMFILE);

    return;
}


# This function finds the appropriate "level id" (mpeg2 encoding parameter)
# for a given resolution. In case the recommended resolution for
# the highest possible level (coded as '4') is exceeded, still '4'
# is returned.
#
# The resolution limits are coded in separate arrays and can be
# adjusted easily if necessary.
sub get_level_id {
    my ($X, $Y)=(@_);

    my @Level_ID = (10, 8, 6, 4);
    my @X_Bound  = (352, 720, 1440, 1920);
    my @Y_Bound  = (240, 480, 960, 1080);

    my $i=0;
    $i++ while ($X > $X_Bound[$i]) and $i < $#X_Bound;
    $i++ while ($Y > $Y_Bound[$i]) and $i < $#Y_Bound;
    return $Level_ID[$i];
}


# Determine version of this script
# (from the CVS ID in the header of this script or
#  from the hard-coded VERSION constant)
sub get_line_no {
    my ($string) = $_[0];
    my $line = 0;
    my $foundIn = "";

    # Open this script for reading
    open(FILE, "<", $0);
    if (! eof(FILE)) {
        while (<FILE>) {
	    $line++;
            if (m/$string/) {
                $foundIn = $line;
                last;
            }
        }
    }
    close(FILE);

    # Fall back to hard-coded version number if version number unset
    $foundIn = "<unknown>" if ($foundIn eq "");

    return $foundIn;
}


sub cleanExit
{
  my $message;

  $message = shift || "Cancelled";
  print STDERR "$progname: $message.\n";

  exit 1;
}


# Function to have a Unix-portable which within perl
# (inspired from File::Which 0.05)
sub which {
    my ($exec) = @_;

    return undef unless $exec;

    my $all = wantarray;
    my @results = ();
    my @path_ext = ('');

    my @path = File::Spec->path();
    unshift @path, File::Spec->curdir;

    for my $base (map { File::Spec->catfile($_, $exec) } @path) {
	for my $ext (@path_ext) {
            my $file = $base.$ext;
# print STDERR "$file\n";

            if ((-x $file)      # executable, normal case
		and !-d _)      # and we don't want dirs to pass (as they are -x)
            {
# print STDERR "-x: ", -x $file, " -e: ", -e _, " -d: ", -d _, "\n";
		return $file unless $all;
		push @results, $file;       # Make list to return later
            }
        }
    }

    if ($all) {
        return @results;
    } else {
        return undef;
    }
}
# End of gmvmpeg
