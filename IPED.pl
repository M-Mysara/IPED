##......................Packages Used..........................
use Getopt::Std;
use Cwd;
use File::Basename;
use strict;
#opening fasta file (trimmed)
open FH, $ARGV[0]; 
my @fasta=<FH>;
close FH;
#opening name file (trimmed)
open FH, $ARGV[1]; 
my @name=<FH>;
close FH;
#opening fasta file after contig formation
open FH, $ARGV[2]; 
my @contig_Fasta=<FH>;
close FH;
#opening qual file after contig formation (modified by M.Mysara)
open FH, $ARGV[3]; 
my @contig_Qual=<FH>;
close FH;
my $oo=$ARGV[4];
mkdir ("temp_split");
unlink("./temp_split/$oo.Test.arff");
open FH, ">>./temp_split/$oo.Test.arff";
print FH '@relation sirnahabal

@attribute FPos numeric
@attribute Fhomo numeric
@attribute FQual numeric
@attribute RQual numeric
@attribute RPos numeric
@attribute RFType {E,M,T,D,I,U}
@attribute Error {T,N}

@data
';
my @Main_counter; 								#In case we need not to treat very frequent reads
my @seq_start;
												#To mark the length of pure sequences "." and "-"
for(my $i=0;$i<scalar(@fasta);$i=$i+2){
	my $freq=scalar(split(",",$name[$i/2]));
	if($freq>1000000000000000000){push(@Main_counter,1);}
	else{
		push(@Main_counter,0);
#Getting information from the inPut fasta file 
		my $id=$fasta[$i]; 											#ID
		my $read=uc($fasta[$i+1]);										#Seq
		my $track=1; 	
		my $num=$i;
		#qFlag after extracting ID from contig fasta and contigqual

#Chomp the lines
		chomp($contig_Fasta[$num+1]);
		chomp($contig_Qual[$num+1]);
		chomp($read);
		my @C_fasta;												#Total Contig Nucleotide Positions
		my @F_fasta;												#Total Forward Nucleotide Positions
		my @R_fasta;												#Total Reverse Nucleotide Positions
		my @F_Qual;													#Total Forward Quality value
		my @R_Qual;													#Total Reverse quality value
		my @Contig_Arr=split(" ",$contig_Qual[$num+1]);				#Splited Contig Nucleotide Quality to be concatinated into @F_fasta @R_fasta @F_Qual @R_Qual
		my @Contig_Fas=split("",$contig_Fasta[$num+1]);				#Splited Contig Nucleotide Positions to be concatinated into @C_fasta
		my @F_homo;													#Total homopolymer status for Forward read
		my @R_homo;													#Total homopolymer status for Reverse read
		my $temp_R_homo="x";										#Buffer for former Homopolymer
		my @HOMO=("0","1","2","3","4","5","6","7","8","9","10","11",
		"11","11","11","11","11","11","11","11","11","11");			#Different Homopolymer status;
		my $HOMO_tracker_r=0;										#To count the accumulating Homopolymers
		my @F_Pos;													#Total Nucleotide position according to Forward read
		my @R_Pos; 													#Total Nucleotide position according to Reverse read
		my @Nuc_type;												#Track the type of relationship between the forward and reverse nucleotide
		my $F_Homo_Buffer="";										#to store last nucleotide for forward
		my $F_Homo_Buffer_ind=0;									#to store last nucleotide index for forward
		my $F_Homo_Buffer_status=0;									#to store last nucleotide index for forward
		my $F_Homo_Buffer_status_P=0;								#to store one before the last nucleotide index for forward
		my $F_count=0;												#To count Forward positions
		my $R_count=0;												#To count Reverse positions
		for(my $j=0;$j<scalar(@Contig_Arr);$j++){
			$Contig_Arr[$j]=~/(\d+)([A,G,T,C,N,\.,-])([A,G,T,C,N,\.,-])(\d+)/;
			push(@F_Qual, $1);
			push(@F_fasta, $2);
			push(@R_fasta, $3);
			push(@R_Qual, $4);
			push(@C_fasta, $Contig_Fas[$j]);
			my $temp_R_fasta="-";									#get the Next non gap nucleotide for Reverse Homopolymer identification
			my $qflag=$j+1;											#Qflag
			while($temp_R_fasta eq "-"){$Contig_Arr[$qflag]=~/(\d+)([A,G,T,C,N,\.,-])([A,G,T,C,N,\.,-])(\d+)/;
				$temp_R_fasta=$3;
				$qflag++;
			}
#Delect homo for Forward Read
			if($F_fasta[$j] eq "-"){								#In case of Gap
				push(@F_homo,"0");push (@F_Pos,$F_count);			#O for Opening "Gap", keep the same count
			}	
			elsif($F_fasta[$j] eq "."){								#In case of "."
				push(@F_homo,"0");push (@F_Pos,0);					#M for missing, position count rested (=0)
			}
			elsif($F_fasta[$j] eq "N"){								#In case of N
				push(@F_homo,"-1");									#U for undetermind
				$F_count++;push (@F_Pos,$F_count);					#Increase the count									
				$F_Homo_Buffer="N";									#Update tHE Bufffer for Nucleotide	
				$F_Homo_Buffer_ind=$j;								#Update tHE Bufffer for Nucleotide Position
				$F_Homo_Buffer_status=0;							#Update tHE Bufffer for Nucleotide Homopolymer
				$F_Homo_Buffer_status_P=$F_Homo_Buffer_status;		#Update tHE Former Bufffer for Nucleotide Homopolymer						
			}
			else{
				$F_count++;push (@F_Pos,$F_count);					#Increase the count							
				if($F_fasta[$j] eq $F_Homo_Buffer){					#same as previous position
					if($HOMO[$F_Homo_Buffer_status] eq "11"){
						$F_homo[$F_Homo_Buffer_ind]=$HOMO[$F_Homo_Buffer_status_P+1];
						$F_Homo_Buffer_status_P=$F_Homo_Buffer_status_P+1;
					}
					else{$F_homo[$F_Homo_Buffer_ind]=$HOMO[1];$F_Homo_Buffer_status_P=1;}
					push(@F_homo,"11");
					$F_Homo_Buffer_status=11;
				}
				else{push(@F_homo,$HOMO[0]);$F_Homo_Buffer_status_P=$F_Homo_Buffer_status;$F_Homo_Buffer_status=0;}
					$F_Homo_Buffer=$F_fasta[$j];
					$F_Homo_Buffer_ind=$j;
			}

#Delect homo for Reverse Read
			if($R_fasta[$j] eq "-"){push(@R_homo,"0");push (@R_Pos,$R_count);$HOMO_tracker_r++;}
			elsif($R_fasta[$j] eq "."){push(@R_homo,"0");push (@R_Pos,0);$temp_R_homo=$R_fasta[$j];}
			elsif($R_fasta[$j] eq "N"){push(@R_homo,"-1");$R_count++;push (@R_Pos,$R_count);$temp_R_homo=$R_fasta[$j];}#-1: Unknown
			else{
				$R_count++;push (@R_Pos,$R_count);					#To calculate the position on the reverse read
				if($R_fasta[$j] eq $temp_R_homo){
					$HOMO_tracker_r++;
					push(@R_homo,$HOMO[1]);
					my $gap_fixer=0;
					for(my $o=0; $o<$HOMO_tracker_r;$o++){
						if($R_homo[$j-$o] eq "0"){$gap_fixer++;}
						else{$R_homo[$j-$o]=$HOMO[$o+1-$gap_fixer];}
					}
				}
				elsif($R_fasta[$j] eq $temp_R_fasta){$HOMO_tracker_r=0;push(@R_homo,"11");}#If similar to next one
				else{$HOMO_tracker_r=0;push(@R_homo,$HOMO[0]);}
				$temp_R_homo=$R_fasta[$j];
			}
#Determining the relationship between both called based in forward and reverese							
			if($F_fasta[$j] eq "." || $R_fasta[$j] eq "." ){push(@Nuc_type,"E");}
			elsif($F_fasta[$j] eq "-"){push(@Nuc_type,"D");}
			elsif($R_fasta[$j] eq "-" ){push(@Nuc_type,"I");}
			elsif($F_fasta[$j] eq "N" || $R_fasta[$j] eq "N" ){push(@Nuc_type,"U");}
			elsif($R_fasta[$j] eq $F_fasta[$j]){push(@Nuc_type,"T");}		
			else{push(@Nuc_type,"M");}
		}
#Trim to the input length
		my $start_position=0;
		my $end_position=0;
		my $Flip=0;
		my $seq_=$read;
		$seq_=~s/\.//g;
		$seq_=~s/\-//g;
		my $starter=length($seq_);
		push(@seq_start,$starter);
		$seq_=~s/n/./gi; #########################################################################Temp!! remove without thinking
		if($contig_Fasta[$num+1]=~/^([\w]*)$seq_/){					#Read matches in Farward mode
			$start_position=length($1);
			$end_position=$start_position+length($seq_);
		}
		else{														#try the reverse one
			my $read_=reverse($seq_);
			$read_=~s/A/H/gi;
			$read_=~s/T/A/gi;
			$read_=~s/H/T/gi;
			$read_=~s/C/H/gi;
			$read_=~s/G/C/gi;
			$read_=~s/H/G/gi;
			if($contig_Fasta[$num+1]=~/^([\w]*)$read_/i){			#Read matches in reverse mode
				$start_position=length($1);
				$end_position=$start_position+length($seq_);
				$Flip=1;
			}
			else{													#cannot match Fasta read with Contig read
					print "Error: Something wrong with ".$id." and ".$contig_Fasta[$num]."\n";	
				}
		}
		

		#Printing out the extracted info.
		if($Flip==0){chomp($id);	$id=~s/>//;
			for(my $j=$start_position;$j<$end_position;$j++){	
				my $temp_R_Pos=0;
				my $pre_FH=$F_homo[$j-1];my $pre_FQ=$F_Qual[$j-1];my $pre_FF=$F_fasta[$j-1];my $pre_RF=$R_fasta[$j-1];my $pre_RQ=$R_Qual[$j-1];my $pre_RH=$R_homo[$j-1];
				my $pos_FH=$F_homo[$j+1];my $pos_FQ=$F_Qual[$j+1];my $pos_FF=$F_fasta[$j+1];my $pos_RF=$R_fasta[$j+1];my $pos_RQ=$R_Qual[$j+1];my $pos_RH=$R_homo[$j+1];
				my $F_H=$F_homo[$j];my $R_H=$R_homo[$j];

				#Fixing the start and end missing values
				if($R_Pos[$j]==0){}else{$temp_R_Pos=$R_count-$R_Pos[$j]+1;}
				if($j==$start_position && $start_position == 0){$pre_FH="0";$pre_FQ=0;$pre_FF=".";$pre_RF=".";$pre_RQ=0;$pre_RH="0";}
				elsif($j==$start_position){if(!$pre_FH){$pre_FH=0;}if(!$pre_FQ){$pre_FQ=0;}if(!$pre_FF){$pre_FF=".";}if(!$pre_RF){$pre_RF=".";}if(!$pre_RQ){$pre_RQ=0;}if(!$pre_RH){$pre_RH=0;}}
				elsif($j==$end_position-1 &&$end_position == length($contig_Fasta[$num+1])){$pos_FH="0";$pos_FQ=0;$pos_FF=".";$pos_RF=".";$pos_RQ=0;$pos_RH="0";}
				elsif($j==$end_position-1){if(!$pos_FH){$pos_FH=0;}if(!$pos_FQ){$pos_FQ=0;}if(!$pos_FF){$pos_FF=".";}if(!$pos_RF){$pos_RF=".";}if(!$pos_RQ){$pos_RQ=0;}if(!$pos_RH){$pos_RH=0;}}
				else{}
				print FH $F_Pos[$j].",";
				if(length($F_H)>0){print FH $F_H.",";}else{print FH $id.",";}
				print FH $F_Qual[$j].",";
				print FH $R_Qual[$j].",";
				print FH $temp_R_Pos.",".$Nuc_type[$j].",N\n";
			}
		}
		else{
			chomp($id);$id=~s/>//;
			for(my $j=$end_position-1;$j>=$start_position;$j--){
				my $temp_R_Pos=0;

				#Fixing the start and end missing values
				my $pre_FH=$F_homo[$j-1];my $pre_FQ=$F_Qual[$j-1];my $pre_FF=$F_fasta[$j-1];my $pre_RF=$R_fasta[$j-1];my $pre_RQ=$R_Qual[$j-1];my $pre_RH=$R_homo[$j-1];
				my $pos_FH=$F_homo[$j+1];my $pos_FQ=$F_Qual[$j+1];my $pos_FF=$F_fasta[$j+1];my $pos_RF=$R_fasta[$j+1];my $pos_RQ=$R_Qual[$j+1];my $pos_RH=$R_homo[$j+1];
				if($R_Pos[$j]==0){}else{$temp_R_Pos=$R_count-$R_Pos[$j]+1;}		
				if($j==$start_position && $start_position == 0){$pre_FH="0";$pre_FQ=0;$pre_FF=".";$pre_RF=".";$pre_RQ=0;$pre_RH="0";}
				elsif($j==$start_position){if(!$pre_FH){$pre_FH=0;}if(!$pre_FQ){$pre_FQ=0;}if(!$pre_FF){$pre_FF=".";}if(!$pre_RF){$pre_RF=".";}if(!$pre_RQ){$pre_RQ=0;}if(!$pre_RH){$pre_RH=0;}}
				elsif($j==$end_position-1 &&$end_position == length($contig_Fasta[$num+1])){$pos_FH="0";$pos_FQ=0;$pos_FF=".";$pos_RF=".";$pos_RQ=0;$pos_RH="0";}
				elsif($j==$end_position-1){if(!$pos_FH){$pos_FH=0;}if(!$pos_FQ){$pos_FQ=0;}if(!$pos_FF){$pos_FF=".";}if(!$pos_RF){$pos_RF=".";}if(!$pos_RQ){$pos_RQ=0;}if(!$pos_RH){$pos_RH=0;}}
				else{}	
				print FH $F_Pos[$j].",";
				print FH $F_homo[$j].",";
				print FH $F_Qual[$j].",";
				print FH $R_Qual[$j].",";
				print FH $temp_R_Pos.",".$Nuc_type[$j].",N\n";
			}
		}
	}
}
close FH;


my $model='IPED.model';my $weka='weka.classifiers.meta.Vote';


system"java -Xmx3000M -classpath ./weka.jar $weka -l ./$model -T ./temp_split/$oo.Test.arff -p 0 > ./temp_split/$oo.Test.Final";


system "cut -d \":\" -f3 ./temp_split/$oo.Test.Final | cut -d \" \" -f1 > ./temp_split/$oo.Test.Cut";
open WEKA,"./temp_split/$oo.Test.Cut";
my @res_arr =<WEKA>;
close WEKA;

splice(@res_arr,0,5);

for(my $i=0;$i<scalar(@fasta);$i=$i+2){
	print $fasta[$i];
	if($Main_counter[$i/2]==1){print $fasta[$i+1];}
	else{
		my $seq=$fasta[$i+1];
		chomp($seq);
		my @arr_seq=split("",$seq);
		my @temp_seq=splice(@res_arr,0,$seq_start[$i/2]);
		my $h=0;
		for(my $l=0;$l<scalar(@arr_seq);$l++){
			chomp($temp_seq[$h]);
			if($arr_seq[$l]eq'.'){print '.';}
			elsif($arr_seq[$l]eq'-'){print '-';}
			elsif($temp_seq[$h]eq'T'){$h++;print$arr_seq[$l];}
			elsif($temp_seq[$h]eq'N'){$h++;print"N";}
		}
		print "\n";
	}
}
