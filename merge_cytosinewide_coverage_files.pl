#!/usr/bin/perl
use warnings;
use strict;
my @files = <*v2.txt.gz>;
my %groups;

# lane4_Y2_TA_Hi_9_ACTTGATG_L004_18814_4_R2_trimmed_bismark_bt2.deduplicated.cov.gz
foreach my $file (@files) {
    my $group = $file;
    $group =~ s/CpG_GC-CR-(\d+-(\w\d+))_.+CpG_report_v2.txt.gz//; ### This is the regex you will need to modify to suit your needs
    $group = $2;
    warn "$1\t$2\tgroup: $group\tfile: $file\n";;
    push @{$groups{$group}},$file;    
}


foreach my $group (keys %groups) {

    print "Group: $group\n";
    $group =~/^\d+-(\w\d+)$/;  ### This is the regex you will need to modify to suit your needs

    my $sample = $group;
    warn "$sample\n";;
  
    warn "sample $sample\n";
    my $outfile = "${sample}_merged.CpG_report_v2.txt.gz";
    warn "Outfile will ne named:\n$outfile\n\n"; sleep(1);
    
    open (OUT,"| gzip -c > $outfile") or die $!;
    my %cov;

    foreach my $file (@{$groups{$group}}) {
      	print "Now reading file $file\n";
	
		open (IN,"zcat $file |") or die $!;
		while (<IN>){
			chomp;
			my ($chr,$pos,$strand,$m,$u,$cont1,$cont2) = (split "\t")[0,1,2,3,4,5,6];
			$cov{$chr}->{$pos}->{m} += $m;
			$cov{$chr}->{$pos}->{u} += $u;
      $cov{$chr}->{$pos}->{str} = $strand;
      $cov{$chr}->{$pos}->{c1} = $cont1;
      $cov{$chr}->{$pos}->{c2} = $cont2;
			#print  $cov{$chr}->{$pos}->{m},"\n";
			#print  $cov{$chr}->{$pos}->{u},"\n"; sleep(1);
			
			# print join ("\t",$chr,$pos,$m,$u),"\n"; sleep(1); 
		}
		close IN;
	}
    warn "Now printing a new, merged coverage file\n";
    foreach my $chr(sort keys %cov){
		foreach my $pos(sort {$a<=>$b} keys %{$cov{$chr}}){
			my $perc = sprintf ("%.2f", $cov{$chr}->{$pos}->{m} / ($cov{$chr}->{$pos}->{m} + $cov{$chr}->{$pos}->{u}) *100 );
			print OUT "$chr\t$pos\t$cov{$chr}->{$pos}->{str}\t$cov{$chr}->{$pos}->{m}\t$cov{$chr}->{$pos}->{u}\t$cov{$chr}->{$pos}->{c1}\t$cov{$chr}->{$pos}->{c2}\n";
		}
	}
    close OUT or die $!;

    print "\n"; sleep(1);

}
