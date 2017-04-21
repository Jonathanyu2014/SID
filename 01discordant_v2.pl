#!/usr/bin/perl
#use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use Statistics::Descriptive;
use IO::File;
use threads;
use threads::shared;
my ($Help,$bamlist,$outDir,$refTEs,$nosort,$temp,$run,$nosplit,$Name);
GetOptions(
    "help!"     => \$Help,      ##Help information
    "bamlist=s" => \$bamlist,   #file include all bamfile names prepared to beed processed
    "outDir=s"  => \$outDir,    ##Appoint the output directory.Default is the path of input file.
    "refTEs=s"  => \$refTEs,    ##TE reference file
    "nosort!"   => \$nosort,    ##Input bam file is nosorted.Default is sorted.
    "temp!"     => \$temp,      ##Do not remove intermediate output files. Default is to cleanup.
#   "run=s"     => \$run,       ##The ways to run blast,qsub or local.Default is local.
    "nosplit!"  => \$nosplit,   ##Split the bam file by RG library.Default is no.
    "name=s"    => \$Name,      ##If you input mulit files, you can appoint the name of result file.Default is unknown
);
my $USAGE = <<USAGE;
Usage: perl $0 -bamlist <string> -refTEs <string> [-outDir <string>] [-nosort] [-temp] [-split] [-name]

       -bamlist  A file include all bamfile names prepared to beed processed
       -refTEs   TE reference file
       [-outDir  Appoint the output directory]
       [-nosort  Input bam file is nosorted.Default is sorted]
       [-temp    Do not remove intermediate output files. Default is to cleanup]
       [-nosplit If the input bamfiles haven't  splited by RG library,please choose it.If you choose -nosplit,please ensure the file in bamlist is single.]
       [-name    If you input mulit files, you can appoint the name of result file.Default is unknown.]

Note: $0 requires samtools (v0.1.18), blast to be in the default path
About software:

Name: Specific insertions detector(Sid)
Version: V2.0
Author: ZengYongli Wangyeming
Email: yuqichao\@genomics.cn zengyongli\@genomics.cn wangyeming\@genomics.cn

USAGE

die $USAGE if($Help);
($bamlist && $refTEs) or die $USAGE;
die "$bamlist is empty!\n" if(-z $bamlist);
open IN,"<$bamlist" or die "Fail to open the $bamlist\n";
my @hash_bamlist;
while(<IN>){
    chomp;
    push@hash_bamlist,$_;
}
close IN;
print STDERR "picking unique and repeat reads\n";
my @suffixlist = qw(.bam .txt);
my ($name,$path,$suffix) = fileparse($bamlist,@suffixlist);
$outDir ||="$path/TE";
$nosort ||=0;
$temp   ||=0;
#$run    ||="local";
$Name   ||="unknown";
$outDir=~s/\/$//;
`mkdir $outDir` unless(-e "$outDir");
my $chrlen="$name.chrlen";
open COUT,">$outDir/$chrlen"  or die "Fail to open $outDir/$chrlen\n";
foreach my $bam (@hash_bamlist){
    open IN,"samtools view -h $bam|" or die "Fail to open $bam\n";
    while(my $str=<IN>){
        chomp($str);
        if($str=~/^\@/){
            print COUT "$str\n";
            next;
        }
        last;
    }
    close IN;
    last;
}
close COUT;

my $outfile_info="$name.info";
if($nosplit){
    die "Please ensure the file in bamlist is single due to -split\n" if(($#hash_bamlist)!=0);
    my %info;
    my %library_name;
    my %unique;
    my %repeat;
    my $bam=$hash_bamlist[0];
    open IN,"samtools view $bam|" or die "Fail to open $bam\n";
    if(!($nosort)){
        my $str;
        my %hash;
        my @temp;
        while ($str = <IN>){
            chomp($str);
            my $id = (split /\s+/, $str)[0];
            my $flag=(split /\s+/, $str)[1];
            my $chr= (split /\s+/, $str)[2];
###insertsize###
            my $insertsize=(split /\s+/, $str)[8];
            my $library=$1 if ($str =~ /RG:Z:(\w+)/);
            push @{$info{$library}},abs($insertsize)  if($#{$info{$library}}<(256*1024));
################
            unless(exists $library_name{$library}){
                $library_name{$library} = $library;
                my $unique_file = IO::File->new("$outDir/$library.unique", 'w') or die "Fail to open unique file!\n";
                $unique{$library}=$unique_file;
                my $repeat_file = IO::File->new("$outDir/$library.repeat",'w') or die "Fail to open repeat file!\n";
                $repeat{$library}=$repeat_file;
            }
            my $uu=$unique{$library};
            my $rr=$repeat{$library};
            if (exists $hash{$id}){
                my @arr1 = split /\s+/, $str;
                my @arr2 = split /\s+/, $hash{$id};
                my $ur1 = $1 if ($str =~ /(XT:A:[A-Z])/);
                my $ur2 = $1 if ($hash{$id} =~ /(XT:A:[A-Z])/);
                if (($ur1 eq "XT:A:U" || $ur1 eq "XT:A:M") && ($ur2 eq "XT:A:R" || $arr2[1] & 4)){
                    print $uu "$str\n";
                    print $rr "$hash{$id}\n";
                }elsif (($ur1 eq "XT:A:R" || $arr1[1] & 4) && ($ur2 eq "XT:A:U" || $ur2 eq "XT:A:M")){
                    print $uu "$hash{$id}\n";
                    print $rr "$str\n";
                }elsif ($ur1 eq "XT:A:U" && $ur2 eq "XT:A:U"){
                    if ($arr1[4] > 20 && $arr2[4] < 10){
                        print $uu "$str\n";
                        print $rr "$hash{$id}\n";
                    }elsif ($arr1[4] < 10 && $arr2[4] > 20){
                        print $uu "$hash{$id}\n";
                        print $rr "$str\n";
                    }else{
                        my $r1 = $1 if ($str =~ /X1:i:(\d+)/);
                        my $r2 = $1 if ($hash{$id} =~ /X1:i:(\d+)/);
                        if ($r1 == 0 && $r2 > 0){
                            print $uu "$str\n";
                            print $rr "$hash{$id}\n";
                        }elsif ($r1 > 0 && $r2 == 0){
                            print $uu "$hash{$id}\n";
                            print $rr "$str\n";
                        }
                    }
                }elsif ($ur1 eq "XT:A:U" && $ur2 eq "XT:A:M"){
                    print $uu "$str\n";
                    print $rr "$hash{$id}\n";
                }
                delete $hash{$id};
            }else{
                my @arr1 = split /\s+/, $str;
                my $len;
                if(!($arr1[1] & 2) && ($arr1[1] & 1)){
                    $hash{$id} = $str;
                }elsif (($arr1[1] & 1) && ($arr1[1] & 2)){
                    my $ur1 = $1 if ($str =~ /(XT:A:[A-Z])/);
                    if ($ur1 eq "XT:A:M"){
                        if ($arr1[3] < $arr1[7]){
                            $len = $1 if ($arr1[5] =~ /^(\d+)S/);
                            $hash{$arr1[0]} = $str if ($len >= 20);
                        }else{
                            $len = $1 if ($arr1[5] =~ /(\d+)S$/);
                            if ($len >=20){
                                foreach my $aa (@temp){
                                    my @arr = split /\s+/, $aa;
                                    if ($arr[0] eq $arr1[0]){
                                        print $uu "$aa\n";
                                        print $rr "$str\n";
                                        last;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            push @temp,$str;
            if($#temp>1000){ shift @temp;}
        }
        foreach my $i(values %unique){
            $i->close;
        }
        foreach my $i(values %repeat){
            $i->close;
        }
    }else{
        my ($str1,$str2);
        while($str1=<IN>){
            chomp($str1);
            my @temp1=split /\s+/, $str1;
            my $id1=$temp1[0];
            my $flag1=$temp1[1];
            my $chr1=$temp1[2];
            my $ur1 = $1 if ($str1 =~ /(XT:A:[A-Z])/);
            my $pos1=$temp1[3];
###insertize###
            my $insertsize=$temp1[8];
            my $library1=$1 if ($str1 =~ /RG:Z:(\w+)/);
            push @{$info{$library1}},abs($insertsize)  if($#{$info{$library1}}<(256*1024));
###############
            chomp($str2=<IN>);
            my @temp2=split /\s+/, $str2;
            my $id2=$temp2[0];
            my $flag2=$temp2[1];
            my $chr2=$temp2[2];
            my $ur2 = $1 if ($str2 =~ /(XT:A:[A-Z])/);
            my $pos2=$temp2[3];
###insertsize###
            my $library2=$1 if ($str2 =~ /RG:Z:(\w+)/);
            $insertsize=$temp2[8];
            push @{$info{$library2}},abs($insertsize)  if($#{$info{$library2}}<(256*1024));
################
            die "Please check your input bam flies whether is sorted or another error!\n" unless(($id1 eq $id2)&&($library1 eq $library2));
            unless(exists $library_name{$library1}){
                my $unique_file = (IO::File->new("$outDir/$library1.unique", 'w') or die "Fail to open unique file!\n");
                $unique{$library1}=$unique_file;
                my $repeat_file = (IO::File->new("$outDir/$library1.repeat",'w') or die "Fail to open repeat file!\n");
                $repeat{$library1}=$repeat_file;
            }
            my $len;
            my $uu=$unique{$library1};
            my $rr=$repeat{$library1};
            if( ($flag1&1) && ($flag1&2) && ($flag2&1) && ($flag2&2) ){
                if($ur1 eq "XT:A:U" && $ur2 eq "XT:A:M"){
                    if($pos1>$pos2){ $len = $1 if ($temp2[5] =~ /^(\d+)S/);}
                    elsif($pos1<$pos2){ $len = $1 if ($temp2[5] =~ /(\d+)S$/);}
                    next unless($len >=20);
                    print $uu "$str1\n";
                    print $rr "$str2\n";
                }elsif($ur1 eq "XT:A:M" && $ur2 eq "XT:A:U"){
                    if($pos1<$pos2){ $len = $1 if ($temp1[5] =~ /^(\d+)S/);}
                    elsif($pos1>$pos2){ $len = $1 if ($temp1[5] =~ /(\d+)S$/);}
                    next unless($len >=20);
                    print $rr "$str1\n";
                    print $uu "$str2\n";
                }
            }elsif( ($flag1&1) && (!($flag1&2)) && ($flag2&1) && (!($flag2&2)) ){
                if (($ur1 eq "XT:A:U" || $ur1 eq "XT:A:M") && ($ur2 eq "XT:A:R" || $flag2 & 4)){
                    print $uu "$str1\n";
                    print $rr "$str2\n";
                }elsif (($ur1 eq "XT:A:R" || $flag1 & 4) && ($ur2 eq "XT:A:U" || $ur2 eq "XT:A:M")){
                    print $uu "$str2\n";
                    print $rr "$str1\n";
                }elsif ($ur1 eq "XT:A:U" && $ur2 eq "XT:A:U"){
                    if ($temp1[4] > 20 && $temp2[4] < 10){
                        print $uu "$str1\n";
                        print $rr "$str2\n";
                    }elsif ($temp1[4] < 10 && $temp2[4] > 20){
                        print $uu "$str2\n";
                        print $rr "$str1\n";
                    }else{
                        my $r1 = $1 if ($str1 =~ /X1:i:(\d+)/);
                        my $r2 = $1 if ($str2 =~ /X1:i:(\d+)/);
                        if ($r1 == 0 && $r2 > 0){
                            print $uu "$str1\n";
                            print $rr "$str2\n";
                        }elsif ($r1 > 0 && $r2 == 0){
                            print $uu "$str2\n";
                            print $rr "$str1\n";
                        }
                    }
                }
            }
        }
        foreach my $i(values %unique){
            $i->close;
        }
        foreach my $i(values %repeat){
            $i->close;
        }
    }
    close IN;
    print STDERR "Pick reads is over!\n";
    print STDERR "Calculating the insertsize\n";
    open OUT,">$outDir/$outfile_info" or die "Fail to open $outfile_info\n";
    my %new_info;
    foreach my $key(sort keys %info){
        my $stat = Statistics::Descriptive::Full->new();
        $stat->add_data(@{$info{$key}});
        my $q1 = $stat->percentile(25);
        my $q2 = $stat->median();
        my $q3 = $stat->percentile(75);
        my $min = $q1-2*($q3-$q1);
        my $max = $q3+2*($q3-$q1);
        my $num = ($#{$info{$key}});
        foreach my $value(@{$info{$key}}){
            push@{$new_info{$key}},$value if($value>$min && $value<$max);
        }
        $stat->clear();
    }
    foreach my $key(sort keys %new_info){
        my $stat = Statistics::Descriptive::Full->new();
        $stat->add_data(@{$new_info{$key}});
        my $mean = $stat->mean();
        my $standard_deviation = $stat->standard_deviation();
        print OUT "$key\t$mean\t$standard_deviation\n";
    }
    close OUT;
    print STDERR "Calculation is over!\n";
######create multi_run######
    sub comparison{
        my $unique=shift;
        my $repeat=shift;
        open IN,"<$repeat" or die "Fail to open $repeat\n";
        my $outfastq="$repeat.fa";
        open OUT,">$outfastq" or die "Fail to open $outfastq\n";
        while (<IN>){
            next if (/^\@\w+\s+/);
            my @arr = split;
            print OUT ">$arr[0]\n$arr[9]\n";
        }
        close IN;
        close OUT;
        open OUT,">$repeat.blast.sh" or die "Fail to open $repeat.blast.sh\n";
        print OUT "blastall -p blastn -F F -d $refTEs -m 8 -e 0.2 -W 8 -i $repeat.fa -o $repeat.fa.m8 \n";
        close OUT;
#        if($run eq "qsub"){
#            `perl qsub-sge.pl --ProjectID XXXX --resource vf=1g --jobprefix blast.sh  "$repeat.blast.sh"`;
#        }else{
            `bash "$repeat.blast.sh"`;
#        }
        my $outfilem8="$repeat.fa.m8";
        open IN_M8,"<$outfilem8" or die "Fail to open $outfilem8\n";
        my %read;
        my %read2;
        while (<IN_M8>){
            my @arr = split;
            next if ($arr[3] < 20);
            $read{$arr[0]} = $arr[1] unless (exists $read{$arr[0]});
            $read2{$arr[0]} = $_ unless (exists $read2{$arr[0]});
        }
        close IN_M8;
        open IN_CHR,"<$outDir/$chrlen" or die "Fail to open $outDir/$chrlen\n";
        my $header = "";
        while ( <IN_CHR>){
            $header.=$_;
        }
        close IN_CHR;
        my $outsam="$unique.sam";
        open OUT,">$outsam" or die "Fail to open $outsam\n";
        print OUT "$header";
        open IN,"<$unique" or die "Fail to open $unique\n";
        my $outsam2 ="$repeat.sam";
        open OUT_2,">$outsam2" or die "Fail to open $outsam2\n";
        while(<IN>){
            my @arr = split;
            if (exists $read{$arr[0]}){
                #$_=~s/^([^\s]+)/$1*$read{$arr[0]}/;
                print OUT "$_";
            }
            if(exists $read2{$arr[0]}){
                print OUT_2 "$read2{$arr[0]}\n";
            }
        }
        close IN;
        close OUT;
        close OUT_2;
        my $outfile_bam="$unique.bam";
        my $outfile_srt="$unique.srt";
        my $outfile_srtbam="$unique.srt.bam";
        my $outfile_srtrmdupbam="$unique.srt.rmdup.bam";
        `samtools view -bhS $outsam -o $outfile_bam`;
        `samtools sort $outfile_bam $outfile_srt`;
        `samtools rmdup $outfile_srtbam $outfile_srtrmdupbam`;
        print STDERR "Blasting is over\n";
        if(!($temp)){
            `rm $outfastq $outfilem8`;
            `rm $outsam`;
            `rm $outfile_bam`;
            `rm $outfile_srtbam`;
            `rm $unique`;
            `rm $repeat`;
        #    `rm $repeat.blast.sh`;
        }
    }
    my @runs;
    foreach my $ln (keys %library_name){
        my $u="$outDir/$ln.unique";
        my $r="$outDir/$ln.repeat";
        my $t= threads->create(\&comparison,$u,$r);
        push@runs,$t;
    }
    print STDERR "Blasting\n";
    foreach my $aa(@runs){
        $aa->join();
    }
    print STDERR "Blast is over\n";
    print STDERR "Merging bam file\n";
    my @split_unique;
    my @split_repeat;
    foreach my $ln (keys %library_name){
        my $u="$outDir/$ln.unique.srt.rmdup.bam";
        push@split_unique,$u;
        my $r="$outDir/$ln.repeat.sam";
        push@split_repeat,$r;
    }
	print "@split_repeat\n"; ## Test
    my $merge_file="$outDir/merge.sh";
    my $merge_file2="$outDir/merge2.sh";
    open OUT,">$merge_file" or die "Fail to open $merge_file\n";
    open OUT_2,">$merge_file2" or die "Fail to open $merge_file2\n";
    print OUT "samtools merge -r $outDir/$name.unique.bam";
    print OUT_2 "cat ";
    foreach my $aa(@split_unique){
        print OUT "\t$aa";
    }
    foreach my $aa(@split_repeat){
        print OUT_2 "\t$aa";
    }
    print OUT "\n";
    print OUT_2 ">$outDir/$name.repeat.m8\n";
    close OUT;
    close OUT_2;
    `echo "nosplit uses first cycly!"`; ## TEST
    `bash $merge_file`;
    `bash $merge_file2`;
    `samtools sort -m 2000000000 $outDir/$name.unique.bam $outDir/$name.unique.srt`;
    `samtools rmdup $outDir/$name.unique.srt.bam $outDir/$name.unique.srt.rmdup.bam`;
    if(!($temp)){
        foreach my $aa(@split_unique){
            `rm $aa`;
        }
        foreach my $aa(@split_repeat){
            `rm $aa`;
        }
        `rm $outDir/$name.unique.bam $outDir/$name.unique.srt.bam `; ## $merge_file $merge_file2
    }

}else{
#######create multi_run#######
    sub pick_bam{
        my %info;
        my $bam=shift;
        my ($name1,$path1,$suffix1) = fileparse($bam,@suffixlist);
        my $unique_outfile="$name1.unique";
        my $repeat_outfile="$name1.repeat";
        open IN,"samtools view $bam|" or die "Fail to open $bam\n";
        open ROUT,">$outDir/$repeat_outfile" or die "Fail to open $outDir/$repeat_outfile\n";
        open UOUT,">$outDir/$unique_outfile" or die "Fail to open $outDir/$unique_outfile\n";
        if(!($nosort)){
            my $str;
            my %hash;
            my @temp;
            while ($str = <IN>){
                chomp($str);
                my $id = (split /\s+/, $str)[0];
                my $flag=(split /\s+/, $str)[1];
                my $chr= (split /\s+/, $str)[2];
###insertsize###
                my $insertsize=(split /\s+/, $str)[8];
                my $library=$1 if ($str =~ /RG:Z:(\w+)/);
                push @{$info{$library}},abs($insertsize)  if($#{$info{$library}}<(256*1024));
################
                if (exists $hash{$id}){
                    my @arr1 = split /\s+/, $str;
                    my @arr2 = split /\s+/, $hash{$id};
                    my $ur1 = $1 if ($str =~ /(XT:A:[A-Z])/);
                    my $ur2 = $1 if ($hash{$id} =~ /(XT:A:[A-Z])/);
                    if (($ur1 eq "XT:A:U" || $ur1 eq "XT:A:M") && ($ur2 eq "XT:A:R" || $arr2[1] & 4)){
                        print UOUT "$str\n";
                        print ROUT "$hash{$id}\n";
                    }elsif (($ur1 eq "XT:A:R" || $arr1[1] & 4) && ($ur2 eq "XT:A:U" || $ur2 eq "XT:A:M")){
                        print UOUT "$hash{$id}\n";
                        print ROUT "$str\n";
                    }elsif ($ur1 eq "XT:A:U" && $ur2 eq "XT:A:U"){
                        if ($arr1[4] > 20 && $arr2[4] < 10){
                            print UOUT "$str\n";
                            print ROUT "$hash{$id}\n";
                        }elsif ($arr1[4] < 10 && $arr2[4] > 20){
                            print UOUT "$hash{$id}\n";
                            print ROUT "$str\n";
                        }else{
                            my $r1 = $1 if ($str =~ /X1:i:(\d+)/);
                            my $r2 = $1 if ($hash{$id} =~ /X1:i:(\d+)/);
                            if ($r1 == 0 && $r2 > 0){
                                print UOUT "$str\n";
                                print ROUT "$hash{$id}\n";
                            }elsif ($r1 > 0 && $r2 == 0){
                                print UOUT "$hash{$id}\n";
                                print ROUT "$str\n";
                            }
                        }
                    }elsif ($ur1 eq "XT:A:U" && $ur2 eq "XT:A:M"){
                        print UOUT "$str\n";
                        print ROUT "$hash{$id}\n";
                    }
                    delete $hash{$id};
                }else{
                    my @arr1 = split /\s+/, $str;
                    my $len;
                    if(!($arr1[1] & 2) && ($arr1[1] & 1)){
                        $hash{$id} = $str;
                    }elsif (($arr1[1] & 1) && ($arr1[1] & 2)){
                        my $ur1 = $1 if ($str =~ /(XT:A:[A-Z])/);
                        if ($ur1 eq "XT:A:M"){
                            if ($arr1[3] < $arr1[7]){
                                $len = $1 if ($arr1[5] =~ /^(\d+)S/);
                                $hash{$arr1[0]} = $str if ($len >= 20);
                            }else{
                                $len = $1 if ($arr1[5] =~ /(\d+)S$/);
                                if ($len >=20){
                                    foreach my $aa (@temp){
                                        my @arr = split /\s+/, $aa;
                                        if ($arr[0] eq $arr1[0]){
                                            print UOUT "$aa\n";
                                            print ROUT "$str\n";
                                            last;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                push @temp,$str;
                if($#temp>1000){ shift @temp;}
            }
        }
        else{
            my ($str1,$str2);
            while($str1=<IN>){
                chomp($str1);
                my @temp1=split /\s+/, $str1;
                my $id1=$temp1[0];
                my $flag1=$temp1[1];
                my $chr1=$temp1[2];
                my $ur1 = $1 if ($str1 =~ /(XT:A:[A-Z])/);
                my $pos1=$temp1[3];
###insertize###
                my $insertsize=$temp1[8];
                my $library1=$1 if ($str1 =~ /RG:Z:(\w+)/);
                push @{$info{$library1}},abs($insertsize)  if($#{$info{$library1}}<(256*1024));
###############
                chomp($str2=<IN>);
                my @temp2=split /\s+/, $str2;
                my $id2=$temp2[0];
                my $flag2=$temp2[1];
                my $chr2=$temp2[2];
                my $ur2 = $1 if ($str2 =~ /(XT:A:[A-Z])/);
                my $pos2=$temp2[3];
###insertsize###
                my $library2=$1 if ($str2 =~ /RG:Z:(\w+)/);
                $insertsize=$temp2[8];
                push @{$info{$library2}},abs($insertsize)  if($#{$info{$library2}}<(256*1024));
################
                die "please check your input bam flies whether is sorted or another error!\n" unless(($id1 eq $id2)&&($library1 eq $library2));
                my $len;
                if( ($flag1&1) && ($flag1&2) && ($flag2&1) && ($flag2&2) ){
                    if($ur1 eq "XT:A:U" && $ur2 eq "XT:A:M"){
                        if($pos1>$pos2){ $len = $1 if ($temp2[5] =~ /^(\d+)S/);}
                        elsif($pos1<$pos2){ $len = $1 if ($temp2[5] =~ /(\d+)S$/);}
                        next unless($len >=20);
                        print UOUT "$str1\n";
                        print ROUT "$str2\n";
                    }elsif($ur1 eq "XT:A:M" && $ur2 eq "XT:A:U"){
                        if($pos1<$pos2){ $len = $1 if ($temp1[5] =~ /^(\d+)S/);}
                        elsif($pos1>$pos2){ $len = $1 if ($temp1[5] =~ /(\d+)S$/);}
                        next unless($len >=20);
                        print ROUT "$str1\n";
                        print UOUT "$str2\n";
                    }
                }elsif( ($flag1&1) && (!($flag1&2)) && ($flag2&1) && (!($flag2&2)) ){
                    if (($ur1 eq "XT:A:U" || $ur1 eq "XT:A:M") && ($ur2 eq "XT:A:R" || $flag2 & 4)){
                        print UOUT "$str1\n";
                        print ROUT "$str2\n";
                    }elsif (($ur1 eq "XT:A:R" || $flag1 & 4) && ($ur2 eq "XT:A:U" || $ur2 eq "XT:A:M")){
                        print UOUT "$str2\n";
                        print ROUT "$str1\n";
                    }elsif ($ur1 eq "XT:A:U" && $ur2 eq "XT:A:U"){
                        if ($temp1[4] > 20 && $temp2[4] < 10){
                            print UOUT "$str1\n";
                            print ROUT "$str2\n";
                        }elsif ($temp1[4] < 10 && $temp2[4] > 20){
                            print UOUT "$str2\n";
                            print ROUT "$str1\n";
                        }else{
                            my $r1 = $1 if ($str1 =~ /X1:i:(\d+)/);
                            my $r2 = $1 if ($str2 =~ /X1:i:(\d+)/);
                            if ($r1 == 0 && $r2 > 0){
                                print UOUT "$str1\n";
                                print ROUT "$str2\n";
                            }elsif ($r1 > 0 && $r2 == 0){
                                print UOUT "$str2\n";
                                print ROUT "$str1\n";
                            }
                        }
                    }
                }
            }
        }
        close IN;
        close UOUT;
        close ROUT;
        my $infofile="$name1.info";
        open OUT,">$outDir/$infofile" or die "Fail to open $infofile\n";
        my %new_info;
        foreach my $key(sort keys %info){
            my $stat = Statistics::Descriptive::Full->new();
            $stat->add_data(@{$info{$key}});
            my $q1 = $stat->percentile(25);
            my $q2 = $stat->median();
            my $q3 = $stat->percentile(75);
            my $min = $q1-2*($q3-$q1);
            my $max = $q3+2*($q3-$q1);
            my $num = ($#{$info{$key}});
            foreach my $value(@{$info{$key}}){
                push@{$new_info{$key}},$value if($value>$min && $value<$max);
            }
            $stat->clear();
        }
        foreach my $key(sort keys %new_info){
            my $stat = Statistics::Descriptive::Full->new();
            $stat->add_data(@{$new_info{$key}});
            my $mean = $stat->mean();
            my $standard_deviation = $stat->standard_deviation();
            print OUT "$key\t$mean\t$standard_deviation\n";
        }
        close OUT;
        open IN,"<$outDir/$repeat_outfile" or die "Fail to open $repeat_outfile\n";
        my $outfastq="$name1.repeat.fa";
        open OUT,">$outDir/$outfastq" or die "Fail to open $outfastq\n";
        while (<IN>){
            next if (/^\@\w+\s+/);
            my @arr = split;
            print OUT ">$arr[0]\n$arr[9]\n";
        }
        close IN;
        close OUT;
        my $outfilem8="$name1.repeat.m8";
        my $outfile_blast="$name1.blast.sh";
        open OUT,">$outDir/$outfile_blast" or die "Fail to open $outfile_blast\n";
        print OUT "blastall -p blastn -F F -d $refTEs -m 8 -e 0.2 -W 8 -i $outDir/$outfastq -o $outDir/$outfilem8 \n";
        close OUT;
#        if($run eq "qsub"){
#            `perl qsub-sge.pl --ProjectID XXXXX --resource vf=1g --jobprefix blast.sh  $outDir/$outfile_blast`;
#        }else{
            `bash $outDir/$outfile_blast`;
#        }
        open IN_M8,"<$outDir/$outfilem8" or die "Fail to open $outDir/$outfilem8\n";
        my %read;
        my %read2;
        while (<IN_M8>){
            my @arr = split;
            next if ($arr[3] < 20);
            $read{$arr[0]} = $arr[1] unless (exists $read{$arr[0]});
            $read2{$arr[0]} = $_ unless (exists $read2{$arr[0]});
        }
        close IN_M8;
        open IN_CHR,"<$outDir/$chrlen" or die "Fail to open $outDir/$chrlen\n";
        my $header = "";
        while ( <IN_CHR>){
            $header.=$_;
        }
        close IN_CHR;
        my $outsam="$name1.unique.sam";
        my $outsam2 ="$name1.repeat.sam";
        open OUT,">$outDir/$outsam" or die "Fail to open $outDir/$outsam\n";
        print OUT "$header";
        open IN,"<$outDir/$unique_outfile" or die "Fail to open $unique_outfile\n";
        open OUT_2,">$outDir/$outsam2" or die "Fail to open $outsam2\n";
        while(<IN>){
            my @arr = split;
            if (exists $read{$arr[0]}){
                #$_=~s/^([^\s]+)/$1*$read{$arr[0]}/;
                print OUT "$_";
            }
            if(exists $read2{$arr[0]}){
                print OUT_2 "$read2{$arr[0]}";
            }
        }
        close IN;
        close OUT;
        close OUT_2;
        my $outfile_bam="$name1.unique.bam";
        my $outfile_srt="$name1.unique.srt";
        my $outfile_srtbam="$name1.unique.srt.bam";
        my $outfile_srtrmdupbam="$name1.unique.srt.rmdup.bam";
        `samtools view -bhS $outDir/$outsam -o $outDir/$outfile_bam`;
        `samtools sort $outDir/$outfile_bam $outDir/$outfile_srt`;
        `samtools rmdup $outDir/$outfile_srtbam $outDir/$outfile_srtrmdupbam`;
        print STDERR "Blasting is over\n";
        if(!($temp)){
            `rm $outDir/$outfastq $outDir/$outfilem8`;
            `rm $outDir/$outsam`;
            #`rm $outDir/$outsam2`;
            `rm $outDir/$outfile_bam`;
            `rm $outDir/$outfile_srtbam`;
            `rm $outDir/$repeat_outfile`;
            `rm $outDir/$unique_outfile`;
            `rm $outDir/$outfile_blast`;
        }
    }
    my @run;
    foreach my $bamfile (@hash_bamlist){
        my $t= threads->create(\&pick_bam,$bamfile);
        push@run,$t;
    }
    foreach my $aa(@run){
        $aa->join();
    }
    print STDERR "Pick reads and blast is over\n";
    print STDERR "Merging bam file\n";
    my @split_unique;
    my @split_repeat;
    my @info_hash;
    my $catfile="cat.sh";
    open INFO_OUT,">$outDir/$catfile" or die "Fail to open $catfile\n";
    print INFO_OUT "cat\t";
    foreach my $bamfile (@hash_bamlist){
        my ($name2,$path2,$suffix2) = fileparse($bamfile,@suffixlist);
        my $u="$outDir/$name2.unique.srt.rmdup.bam";
        push@split_unique,$u;
        my $r="$outDir/$name2.repeat.sam";
        push@split_repeat,$r;
        my $info="$outDir/$name2.info";
        push@info_hash,$info;
        print INFO_OUT "$info\t";
    }
    print INFO_OUT ">$outDir/$outfile_info\n";
    close INFO_OUT;
    `bash $outDir/$catfile`;
    my $merge_file="$outDir/merge.sh";
    my $merge_file2="$outDir/merge2.sh";
    open OUT,">$merge_file" or die "Fail to open $merge_file\n";
    open OUT_2,">$merge_file2" or die "Fail to open $merge_file2\n";
    print OUT "samtools merge -r $outDir/$Name.unique.bam";
    print OUT_2 "cat\t";
    foreach my $aa(@split_unique){
        print OUT "\t$aa";
    }
    foreach my $aa(@split_repeat){
        print OUT_2 "\t$aa";
    }
    print OUT "\n";
    print OUT_2 ">$outDir/$Name.repeat.m8\n";
    close OUT;
    close OUT_2;
    `bash $merge_file`;
    `bash $merge_file2`;

    `samtools sort -m 2000000000 $outDir/$Name.unique.bam $outDir/$Name.unique.srt`;
    `samtools rmdup $outDir/$Name.unique.srt.bam $outDir/$Name.unique.srt.rmdup.bam`;
    if(!($temp)){
        foreach my $aa(@split_unique){
            `rm $aa`;
        }
        foreach my $aa(@split_repeat){
            `rm $aa`;
        }
        foreach my $aa(@info_hash){
            `rm $aa`;
        }
        `rm $outDir/$Name.unique.bam $outDir/$Name.unique.srt.bam $merge_file $merge_file2 $outDir/$catfile`;
    }
    `echo "nosplit uses second cycly!"`;
}
