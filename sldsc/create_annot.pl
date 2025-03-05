$CT = $ARGV[0];

print "1/ read ref files\n";
for($CHR = 1; $CHR <= 22; $CHR++){ 
    open(IN,"/project/gazal_569/DATA/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC.$CHR.bim") || die;
    while (<IN>) {
        chomp $_; @line=split; 
        $pip{$line[0]}{$line[3]} = 0;
        $cs{$line[0]}{$line[3]} = 0;
        $pip_cs{$line[0]}{$line[3]} = 0;
    }
    close IN;
}

print "2/ read $CT/$CT.txt\n";
open(IN,"$CT/$CT.txt") || die;
while (<IN>) {
    chomp $_; @line=split;
    if (defined($pip{$line[1]}{$line[2]})){
        $thispip = $line[3];
        $incs = 0; for($i = 4; $i <= $#line; $i++){$incs = $incs + $line[$i]}
        #
        if ($thispip>$pip{$line[1]}{$line[2]}){$pip{$line[1]}{$line[2]} = $thispip} 
        if ($incs > 0){
            $cs{$line[1]}{$line[2]} = 1;
            if ($thispip>$pip_cs{$line[1]}{$line[2]}){$pip_cs{$line[1]}{$line[2]} = $thispip} 
        }
    }
}
close IN;

print "3/ create annotations\n";
for($CHR = 1; $CHR <= 22; $CHR++){ 
    open(IN,"/project/gazal_569/DATA/ldsc/reference_files/1000G_EUR_Phase3/plink_files/1000G.EUR.QC.$CHR.bim") || die;
    open(OUT1,"> $CT/pip.$CHR.annot")    || die; printf OUT1 "pip\n";
    open(OUT2,"> $CT/cs.$CHR.annot")     || die; printf OUT2 "cs\n";
    open(OUT3,"> $CT/pip_cs.$CHR.annot") || die; printf OUT3 "pip_cs\n";
    while (<IN>) {
        chomp $_; @line=split; 
        printf OUT1 "$pip{$line[0]}{$line[3]}\n";
        printf OUT2 "$cs{$line[0]}{$line[3]}\n";
        printf OUT3 "$pip_cs{$line[0]}{$line[3]}\n";
    }
    close IN;
    close OUT1;
    close OUT2;
    close OUT3;
}
system ("gzip $CT/*.annot");
