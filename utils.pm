#!/usr/bin/perl
#	Copyright 2011,2012 Dmitri Pervouchine (dp@crg.eu)
#	This file is a part of the IRBIS package.
#	IRBIS package is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#	
#	IRBIS package is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#	
#	You should have received a copy of the GNU General Public License
#	along with IRBIS package.  If not, see <http://www.gnu.org/licenses/>.


use setup;

##########################################################################################
# This package contains subroutines
# 1) load_global_configuration_file
# 2) read_configuration
# 3) read_annotation
# 4) read_introns
# 5) read_signatures
# 6) read_sequences
# 7) read_expression
#
# 8) make_download
# 9) make_clade
###########################################################################################
    %ABBR = %UCSC = %DATATYPE = %REFDB = %CLADE = %NICKNAME = %FULLNAME = ();
    load_global_configuration_file();

#    $READKEY  = 'yes' if(`perl -e 'use Term::ReadKey;print 1;' 2>>/dev/null`);
#    if($READKEY) {
#        use Term::ReadKey;
#    }
    return(1);
##########################################################################################

sub progressbar {
    my ($current, $last, $message) = @_;
    my $width = 2**int(log($WCHAR)/log(2));
    my $i;
    if(int(($width*($current-1))/$last) < int(($width*$current)/$last)) {
        my $k = int(($width*$current)/$last);
        print STDERR "\r$message\[";
        for($i=0;$i<$k;$i++) {print STDERR "=";}
        print STDERR ">" if($k<$width);
        for($i++;$i<$width;$i++) { print STDERR " ";}
        print STDERR "\] ",($current/$last < 0.1 ? " " : ""), int(100*$current/$last),"%";
    }
    print STDERR "\n" if($current==$last);
}

sub assure {
    (my $a, my $b) = @_;
    die("Conflict: $a vs $b\n") if($a && $a ne $b);
}

##########################################################################################

sub load_global_configuration_file {
    my $cfg_file = "config.dat";
    open FILE, $cfg_file || die("Config file ($cfg_file) does not exist, exiting\n");
    while(my $line=<FILE>) {
	next if($line=~/^#/);
	chomp $line;
	(my $id) = split /\s+/, $line;
	(my $id, $ABBR{$id}, $UCSC{$id}, $DATATYPE{$id}, $REFDB{$id}, $CLADE{$id}, $NICKNAME{$id}, $FULLNAME{$id}) = split /\s+/, $line;
	push @{$CLADELIST{$CLADE{$id}}}, $id;
	$basespecies{$CLADE{$id}} = $id unless($basespecies{$CLADE{$id}});
    }
    close FILE;
}

sub make_download {
    open FILE, ">download";
    print FILE "include makefile\n";
    foreach my $z(keys(%UCSC)) {
	next unless($UCSC{$z});
        my $u = $UCSC{$z};	# ucsc name
	my $w = $UCSC{$z};	# ucsc name first character capital (for vsMm9, for instance)
	substr($w,0,1) =~ tr/[a-z]/[A-Z]/;
	$name    = "\$(SOURCE)$z/md5sum.txt ";
	$command = "\$(SOURCE)$z/md5sum.txt:\n";

	if($DATATYPE{$z} eq "fa") {
	    $command.= "\tmkdir -p \$(SOURCE)$z\n";
	    $command.= "\twget http://hgdownload.cse.ucsc.edu/goldenPath/$u/bigZips/$u.fa.gz -O \$(SOURCE)$z/$u.fa.gz\n";
	    $command.= "\tgunzip -f \$(SOURCE)$z/$u.fa.gz\n";
	    $command.= "\twget http://hgdownload.cse.ucsc.edu/goldenPath/$u/bigZips/$u.fa.masked.gz -O \$(SOURCE)$z/$u.fa.masked.gz\n";
	    $command.= "\tgunzip -f \$(SOURCE)$z/$u.fa.masked.gz\n";
	}
	if($DATATYPE{$z} eq "scaffold") {
	    $command.= "\tmkdir -p \$(SOURCE)$z\n";
            $command.= "\twget http://hgdownload.cse.ucsc.edu/goldenPath/$u/bigZips/scaffoldFa.gz -O \$(SOURCE)$z/scaffold.fa.gz\n";
	    $command.= "\tgunzip -f \$(SOURCE)$z/scaffold.fa.gz\n";
            $command.= "\twget http://hgdownload.cse.ucsc.edu/goldenPath/$u/bigZips/scaffoldFaMasked.gz -O \$(SOURCE)$z/scaffold.fa.masked.gz\n";
            $command.= "\tgunzip -f \$(SOURCE)$z/scaffold.fa.masked.gz\n";
	}
	if($DATATYPE{$z} eq "scaffold-zip") {
	    $command.= "\tmkdir -p \$(SOURCE)$z\n";
            $command.= "\twget http://hgdownload.cse.ucsc.edu/goldenPath/$u/bigZips/scaffoldFa.zip -O \$(SOURCE)$z/scaffold.fa.zip\n";
            $command.= "\tunzip -f \$(SOURCE)$z/scaffold.fa.zip -d\$(SOURCE)$z/\n";
            $command.= "\twget http://hgdownload.cse.ucsc.edu/goldenPath/$u/bigZips/scaffoldFaMasked.zip -O \$(SOURCE)$z/scaffold.fa.masked.zip\n";
            $command.= "\tunzip -f \$(SOURCE)$z/scaffold.fa.masked.zip -d \$(SOURCE)$z/\n";
	}
	if($DATATYPE{$z} eq "chrom-zip") {
	    $command.= "\tmkdir -p \$(SOURCE)$z\n";
	    $command.= "\twget http://hgdownload.cse.ucsc.edu/goldenPath/$u/bigZips/chromFa.zip -O \$(SOURCE)$z/chromFa.zip\n";
	    $command.= "\tunzip \$(SOURCE)$z/chromFa.zip -d \$(SOURCE)$z/\n";
	    $command.= "\twget http://hgdownload.cse.ucsc.edu/goldenPath/$u/bigZips/chromFaMasked.zip -O \$(SOURCE)$z/chromFaMasked.zip\n";
	    $command.= "\tunzip \$(SOURCE)$z/chromFaMasked.zip -d \$(SOURCE)$z/\n";
	}
	if($DATATYPE{$z} eq "chrom" || $DATATYPE{$z} eq "chrom-sub" || $DATATYPE{$z} eq "chrom-soft") {
            $command.= "\tmkdir -p \$(SOURCE)$z\n";
	    $command.= "\twget http://hgdownload.cse.ucsc.edu/goldenPath/$u/bigZips/chromFa.tar.gz -O \$(SOURCE)$z/chromFa.tar.gz\n";
            $command.= "\tgunzip -f \$(SOURCE)$z/chromFa.tar.gz\n";
	    $command.= "\ttar -xf \$(SOURCE)$z/chromFa.tar -C \$(SOURCE)$z/\n";
            $command.= "\twget http://hgdownload.cse.ucsc.edu/goldenPath/$u/bigZips/chromFaMasked.tar.gz -O \$(SOURCE)$z/chromFaMasked.tar.gz\n";
            $command.= "\tgunzip -f \$(SOURCE)$z/chromFaMasked.tar.gz\n";
            $command.= "\ttar -xf \$(SOURCE)$z/chromFaMasked.tar -C \$(SOURCE)$z/\n";
	    if($DATATYPE{$z} eq "chrom-sub") {
            	$command.= "\tmv \$(SOURCE)$z/?/*  \$(SOURCE)$z/\n";
            	$command.= "\tmv \$(SOURCE)$z/??/* \$(SOURCE)$z/\n";
            	$command.= "\trm -f -r \$(SOURCE)$z/? \$(SOURCE)$z/??\n";
	    }
            if($DATATYPE{$z} eq "chrom-soft") {
                $command.= "\tmv \$(SOURCE)$z/softMask/* \$(SOURCE)$z/\n";
                $command.= "\tmv \$(SOURCE)$z/hardMask/* \$(SOURCE)$z/\n";
                $command.= "\trm -f -r \$(SOURCE)$z/softMask/ \$(SOURCE)$z/hardMask/\n";
            }

	    $command.= "\trm -f \$(SOURCE)$z/chromFa.tar \$(SOURCE)$z/chromFaMasked.tar\n";
	}
	$command.= "\twget http://hgdownload.cse.ucsc.edu/goldenPath/$u/bigZips/md5sum.txt -O \$(SOURCE)$z/md5sum.txt\n";

	unless($z eq $basespecies{$CLADE{$z}}) {	
	    $name.=    "\$(CHAIN)$REFDB{$z}.$z.all.chain ";
	    $command.= "\$(CHAIN)$REFDB{$z}.$z.all.chain:\n";
	    $command.= "\tmkdir -p \$(CHAIN)\n";
	    $command.= "\twget http://hgdownload.cse.ucsc.edu/goldenPath/$REFDB{$z}/vs$w/$REFDB{$z}.$u.all.chain.gz -O \$(CHAIN)$REFDB{$z}.$z.all.chain.gz\n";
	    $command.= "\tgunzip -f \$(CHAIN)$REFDB{$z}.$z.all.chain.gz\n";
	}

	$command.= "remove-$z\::\n\trm -f -r \$(SOURCE)$z \$(CHAIN)$REFDB{$z}.$z.all.chain\n";
	print FILE "$command\n\n";

    }
    close FILE;
}

sub make_clade {
    my $clade = @_[0];
    my @genes = @{@_[1]};
    my @biotypes = @{@_[2]};

    my $WINDOW = 150;
    my $DOUBLE = 2*$WINDOW;
    my $MARGIN = 5;

    my $aln_command = $out_command = $bed_command = $sus_command = $suw_command = $aln_name = $out_name = $bed_name = $sus_name = $suw_name = $idx_command = undef;
    my $cps_command = $mus_command = $muw_command = $cps_name = $mus_name = $muw_name = $rem_command = undef;
    foreach my $z(@{$CLADELIST{$clade}}) {
	my $f = ($z eq $basespecies{$clade}) ? "\$(METADATA)$clade/$REFDB{$z}.cps" : "\$(METADATA)$clade/$z.out";
	my $g = ($z eq $basespecies{$clade}) ? "-cis" : undef;
        $idx_command.= "\$(SEQ)$clade/$z.idx \$(SEQ)$clade/$z.dbx:\t\$(SOURCE)$z/md5sum.txt\n\tmkdir -p \$(SEQ)$clade/\n\t\$(TRANSF) -i \$(SOURCE)$z -d \$(SEQ)$clade/$z -r\n";

	unless($z eq $basespecies{$clade}) {
    	    $aln_command.= "\$(METADATA)$clade/$z.aln:\t\$(CHAIN)$REFDB{$z}.$z.all.chain \$(METADATA)$clade/$REFDB{$z}.cps.srt\n\tmkdir -p \$(METADATA)$clade/\n";
	    $aln_command.= "\t\$(MAPPER) -i \$(METADATA)$clade/$REFDB{$z}.cps.srt -d \$(CHAIN)$REFDB{$z}.$z.all.chain -W 0 -o \$(METADATA)$clade/$z.aln\n";
    	    $out_command.= "\$(METADATA)$clade/$z.out:\t\$(METADATA)$clade/$z.aln \$(METADATA)$clade/$REFDB{$z}.cps\n";
	    $out_command.= "\t\$(MATCHER) -i \$(METADATA)$clade/$REFDB{$z}.cps -d \$(METADATA)$clade/$z.aln -o \$(METADATA)$clade/$z.out\n";
	    $aln_name.= "\$(METADATA)$clade/$z.aln ";
	    $out_name.= "\$(METADATA)$clade/$z.out ";
	}

	$bed_command.= "\$(METADATA)$clade/$z.bed:\t$f\n\t\$(BEDDER) $g -i $f -o \$(METADATA)$clade/$z.bed\n";
    	$sus_command.= "\$(METADATA)$clade/$z.sus:\t\$(METADATA)$clade/$z.bed \$(SEQ)$clade/$z.idx \$(SEQ)$clade/$z.dbx \$(SEGMENT)\n";
	$sus_command.= "\t\$(SEGMENT) -i \$(METADATA)$clade/$z.bed -d \$(SEQ)$clade/$z -o \$(METADATA)$clade/$z.sus.uns \$(SEGPAR)\n";
	$sus_command.= "\tsort -n -k 1 \$(METADATA)$clade/$z.sus.uns > \$(METADATA)$clade/$z.sus\n\trm -f \$(METADATA)$clade/$z.sus.uns\n";
	$suw_command.= "\$(METADATA)$clade/$z.sus.ind:\t\$(METADATA)$clade/$z.sus\n\t\$(INDEXING) -dir \$(METADATA)$clade/ -file $z.sus -suffix .ind\n";

	if($z eq $basespecies{$clade}) {
       	    $suw_command.= "\$(METADATA)$clade/$REFDB{$z}.sgn:\t\$(METADATA)$clade/$REFDB{$z}.cps \$(SEQ)$clade/$z.idx \$(SEQ)$clade/$z.dbx\n";
    	    $suw_command.= "\t\$(WINDOW) -i \$(METADATA)$clade/$REFDB{$z}.cps -we 10 -wi 10 -d \$(SEQ)$clade/$z -o \$(METADATA)$clade/$REFDB{$z}.sgn -c -t\n";
	}

    	$bed_name.= "\$(METADATA)$clade/$z.bed ";
    	$sus_name.= "\$(METADATA)$clade/$z.sus ";
	$suw_name.= "\$(METADATA)$clade/$z.sus.ind ";

	$rem_command.= "remove-$z\::\n\trm -f \$(METADATA)$clade/$z.* \$(SEQ)$clade/$z.idx \$(SEQ)$clade/$z.dbx\n";
    }

    my $z = $basespecies{$clade};

    foreach $gene(@genes) {
	$cps_command.= "\$(METADATA)$clade/$gene.cps:\t\$(METADATA)$clade/$REFDB{$z}.cps\n\tawk '\$\$11==\"$gene\"' \$(METADATA)$clade/$REFDB{$z}.cps > \$(METADATA)$clade/$gene.cps\n";
	$mus_command.= "\$(METADATA)$clade/$gene.mus:\t\$(METADATA)$clade/$gene.cps $sus_name\n\t\$(MUFFER) -i $sus_name -o \$(METADATA)$clade/$gene.mus < \$(METADATA)$clade/$gene.cps\n";
	$cps_name.= "\$(METADATA)$clade/$gene.cps ";
	$mus_name.= "\$(METADATA)$clade/$gene.mus ";
    }

    $subsets_command = $subsets_name = undef;
    foreach $biotype("protein_coding", "snoRNA", "snRNA", @biotypes) {
	$subsets_command.= "\$(METADATA)$clade/$biotype.cps: \$(METADATA)$clade/$REFDB{$z}.cps\n";
	$subsets_command.= "\tawk '\$\$10==\"$biotype\"' \$(METADATA)$clade/$REFDB{$z}.cps > \$(METADATA)$clade/$biotype.cps\n";
	$subsets_name.= "\$(METADATA)$clade/$biotype.cps ";
    }

    $subsets_command.= "\$(METADATA)$clade/ncpcg.cps: \$(METADATA)$clade/protein\_coding.cps\n";
    $subsets_command.= "\tawk '\$\$15~/[EAI][5N3]/' \$(METADATA)$clade/protein\_coding.cps > \$(METADATA)$clade/ncpcg.cps\n";
    $subsets_command.= "\$(METADATA)$clade/cspcg.cps: \$(METADATA)$clade/protein\_coding.cps\n";
    $subsets_command.= "\tawk '\$\$15~/[EAI]C/' \$(METADATA)$clade/protein\_coding.cps > \$(METADATA)$clade/cspcg.cps\n";
    $subsets_name.=    "\$(METADATA)$clade/ncpcg.cps \$(METADATA)$clade/cspcg.cps ";

    $all_name.="$clade ";
    $clean_name.="$clade-clean ";
    $remove_name.="$clade-remove ";

    return("$clade\::\t$sus_name $suw_name $mus_name $muw_name $subsets_name \$(METADATA)$clade/$REFDB{$z}.sgn \$(METADATA)$clade/$REFDB{$z}.a2a\n\n".
     "$idx_command\n$aln_command\n$out_command\n$bed_command\n$sus_command\n$suw_command\n$cps_command\n$mus_command\n$muw_command\n$subsets_command\n$rem_command\n".
     "$clade-clean\::\n\trm -f $sus_name $suw_name $mus_name $muw_name\n".
     "$clade-remove\::\t$clade-clean\n\trm -f $aln_name $out_name $bed_name $cps_name $subsets_name \$(METADATA)$clade/$REFDB{$z}.sgn \$(METADATA)$clade/$REFDB{$z}.cps ".
     "\$(METADATA)$clade/$REFDB{$z}.cps.srt \$(METADATA)$clade/$REFDB{$z}.int \$(METADATA)$clade/$REFDB{$z}.a2a\n".
     "\$(METADATA)$clade/$REFDB{$z}.cps.srt:\t\$(METADATA)$clade/$REFDB{$z}.cps\n\tsort -k 1,1 -k 2,2n \$(METADATA)$clade/$REFDB{$z}.cps > \$(METADATA)$clade/$REFDB{$z}.cps.srt\n".
     "\$(METADATA)$clade/$REFDB{$z}.a2a: \$(METADATA)$clade/$REFDB{$z}.cps\n\t./_all2all -i \$(METADATA)$clade/$REFDB{$z}.cps -o \$(METADATA)$clade/$REFDB{$z}.a2a -b\n".
     "index::\n\t./indexing -dir \$(METADATA)$clade/ -extension .sus -suffix .ind\n\n");
}


###########################################################################################################
sub read_configuration {
# Read configuration file specified in $config_file_name
# Populate variables $annotation_file_name and $introns_file_name and $path
    $config_file_name = @_[0] if(@_[0]);
    return unless($config_file_name);
    print STDERR "[Configuration:\t$config_file_name";
    unless($filewasread{$config_file_name}) {
	$path = $annotation_file_name = $introns_file_name = $signatures_file_name = $extention = undef;
	%species = %index = ();
    	$nsp = 0;
    	open FILE, $config_file_name || die(" - can't open '$config_file_name'\n");
    	while($line=<FILE>) {
	    chomp $line;
            ($term, $val, $weight) = split /\s+/, $line;
	    $halfsize		  = $val if($term eq "halfsize");
	    $gapsize		  = $val if($term eq "gap");
	    if($term eq "path") {
		$path = $val;
		$path = "$ENV{HOME}$1" if($path=~/^\~(.*)$/);
	    }
            $annotation_file_name = $path.$val if($term eq "annotation");
            $introns_file_name    = $path.$val if($term eq "introns");
	    $signatures_file_name = $path.$val if($term eq "signatures");
            $extention		  = $val if($term eq "extention");
            $index{$val}	  = $nsp if($term eq "species");
	    $weigh{$val}	  = $weight if($term eq "species");
            $species[$nsp++]	  = $val if($term eq "species");
	}
	close FILE;
	$filewasread{$config_file_name}=1;
    }
    print STDERR "]\n";
}

sub read_annotation {
    return unless($annotation_file_name);
    print STDERR "[Annotation:\t$annotation_file_name";
    unless($filewasread{$annotation_file_name}) {
    	open FILE, $annotation_file_name || die("- can't open '$annotation_file_name'\n");
    	while($line=<FILE>) {
            chomp $line;
	    ($chr, $pos, $str, $gene, $site, $type, $ups, $dws, $ensg, $biotype, $name, $use, $tot, $tup, $tdw, $sup, $sdw) = split /\t/, $line;
	    $DATA{$site}     = $line;
            $GENE{$site}     = $gene;
	    $SEGM{$site}     = $tdw;
	    $CHR{$site}	     = $chr;
	    $POSITION{$site} = $pos;
	    $STRAND{$site}   = ($str eq "+") ? 1 : -1;
	    $TYPE{$site}     = $type;
	    $DESC{$site}     = "gene=$ensg name=$name segment=$sdw type=$tdw";
	    $DESC{$site} =~ s/\_/\\\_/g; 
	    push @{$SITES{$gene}}, $site;
	    $NAME{$gene}     = $name;
        }
        close FILE;
	print STDERR " ",0+keys(%DATA);
	$filewasread{$annotation_file_name} = 1;
    }
    print STDERR "]\n";
}

sub getbox {
    my $site = @_[0];
    my $x = $POSITION{$site} + $STRAND{$site}*@_[1];
    my $y = $POSITION{$site} + $STRAND{$site}*@_[2];
    ($x, $y) = sort {$a<=>$b} ($x, $y);
    return("$CHR{$site}\_$x\_$y\_$STRAND{$site}");
}

sub read_introns {
    return unless($introns_file_name);
    print STDERR "[Introns:\t$introns_file_name";
    unless($filewasread{$introns_file_name}) {
    	open FILE, $introns_file_name || die(" - can't open '$introns_file_name'\n");
    	while($line=<FILE>) {
            chomp $line;
	    ($left, $right) = split /\t/, $line;
	    my $gene = $GENE{$left};
   	    next unless($gene);
	    push @{$INTRONS{$gene}}, [$left, $right];
        }
        close FILE;
	print STDERR " ",0+keys(%INTRONS);
	$filewasread{$introns_file_name}=1;
    }	
    print STDERR "]\n";
}

sub read_signatures {
    return unless($signatures_file_name);
    print STDERR "[Signatures:\t$signatures_file_name";
    unless($filewasread{$signatures_file_name}) {
    	open FILE, $signatures_file_name || die(" can't open '$signatures_file_name'\n");
    	while($line=<FILE>) {
            chomp $line;
	    ($site, $chr, $pos, $len, $str, $full, $seq) = split /\t/, $line;
            $seq =~ tr/[a-z]/[A-Z]/;
            $LSGN{$site} = substr($seq, 0,10);
            $RSGN{$site} = substr($seq,10,10);
        }
        close FILE;
	print STDERR " ",0+keys(%LSGN);
	$filewasread{$signatures_file_name}=1;
    }
    print STDERR "]\n";
}

sub read_sequences {
    my %selected = @_;
    if($maf_file_name) {
  	print STDERR "[MAF input $maf_file_name";
  	open FILE, $maf_file_name || die(" can't open '$maf_file_name'\n");
  	while($line=<FILE>) {
    	    if($line=~/^a id=(\d+)/) {
	    	$id = $1;
	    	next unless($selected{$id});
	    	while($line=<FILE>) {
	    	    ($s, $chr, $pos, $len, $str, $full, $seq) = split /\s+/, $line;
	    	    last unless($s eq "s");
	    	    ($name) = split /\./,$chr;
	    	    $org = $index{$name};
	    	    next if($org eq undef);
            	    $seq =~ tr/[A-Z]/[a-z]/;
	    	    $attr[$org][$id] = "$chr\t$pos\t$len\t$str\t$full";
		    $seq = revcomp($seq) if($rc);
	    	    $data[$org][$id] = $seq;
	    	}
    	    }
  	}
  	close FILE;
	print STDERR "(rc)" if($rc);
	print STDERR "]\n";
    }
    else {
  	for($i=0;$i<$nsp;$i++) {
    	    $filename = "$path$species[$i]$extention";
    	    print STDERR "[Reading $filename";
    	    open FILE, "$filename" || die(" can't open '$filename'\n");
	    if(-e "$filename.ind") {
		%ind = split /[\t\n]/, `cat $filename.ind`;
		print STDERR "(ind)";
		foreach $id(keys(%selected)) {
		    seek FILE, $ind{$id}, 0;
		    $line=<FILE>;
		    chomp($line);
		    ($id, $chr, $pos, $len, $str, $full, $seq) = split /\t/, $line;
		    $seq =~ tr/[A-Z]/[a-z]/;
		    $attr[$i][$id] = "$species[$i].$chr\t$pos\t$len\t$str\t$full";
                    $seq = revcomp($seq) if($rc);
                    $data[$i][$id] = $seq;
		}
	    }
	    else {
    	    	while($line=<FILE>) {
 		    chomp($line);
		    ($id, $chr, $pos, $len, $str, $full, $seq) = split /\t/, $line;
		    if($selected{$id}) {
            	    	$seq =~ tr/[A-Z]/[a-z]/;
	    	    	$attr[$i][$id] = "$species[$i].$chr\t$pos\t$len\t$str\t$full";
		    	$seq = revcomp($seq) if($rc);
	    	    	$data[$i][$id] = $seq;
		    }
    	        }
	    }
	    close FILE;
	    print STDERR "(rc)" if($rc);
    	    print STDERR "]\n";
  	}
    }
}

sub revcomp {
    my $res = join("", reverse split //,@_[0]);
    $res =~ tr/[acgt]/[tgca]/;
    return($res);
}

####################################################################################################################################################

sub read_expression {
    my $expression_file_name = @_[0];
    return unless($expression_file_name);
    open FILE, $expression_file_name || die("Can't open '$expression_file_name'\n");
    while($line=<FILE>) {
    	chomp $line;
	my @arr = split /\s+/, $line;
	my $key = shift(@arr);
	if($key eq "data") {
	    my @files=();
	    foreach my $tag(@arr) {
		$id = $1 if($tag =~ /id=(.+)/);
		push @files, $1 if($tag =~ /file=(.+)/);
	    }
	    next unless($id);
	    foreach $file(@files) {
		read_datafile($id, $file);
	    }
	}
	if($key eq "plot") {
	    push @PLOTS, $line 
	}
    }
    close FILE;
}

sub read_datafile {
    my $id = @_[0];
    my $file = @_[1];
    $file = "$ENV{HOME}$1" if($file=~/^\~(.*)$/);
    print STDERR "[Data: $file -> $id ";
    open INPUT, $file || die("Can't open '$file'\n");
    my @header = ();
    while($line=<INPUT>) {
        chomp $line;
        my @arr = split /\s+/, $line;
        if(@header) {
            my $key = shift(@arr);
            next unless($key);
            for($i=0;$i<@header;$i++) {
            	$VALUE{$id}{$key}{$header[$i]} = $arr[$i];
            }
	} 
        else {
            @header = @arr;
        }
    }
    close INPUT;
    print STDERR 0+keys(%{$VALUE{$id}}),"]\n";
}
