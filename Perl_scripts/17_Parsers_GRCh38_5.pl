# Manuel Tardaguila, PhD. Sanger Institute
  
# BEWARE!. This is a script to turn an fmt3 blastn nt alignment into a vcf like file
# it needs further development
# 1. Strings of the alignment shouldn't be analyzed in series of 60 nts, but as a continuous sequence. THIS is serious, indels distributed acroos two lines or more will appear as different indels
# 2. All the deletions part needs testing

#!/usr/bin/perl


use strict;
use warnings;
use Time::localtime;

my $input1=$ARGV[0];
my $output1=$ARGV[1];


print "$input1\n";
print "$output1\n";
my %hash1=();
my %hash1_rev=();
my %hash2=();
my %hash3=();





my %hash_alignment_parsing=();
my $query="NaN";
my $length_query="NaN";
my $subject="NaN";
my $length_subject="NaN";
my $FLAG_line=0;
my $Flag_hit=1;
my @seq_query=();
my @seq_subject=();
 
my $uncover_N_TER=0;
my $uncover_C_TER=0;
my $coverage=0;
my $mismatch=0;
my $identity=0;
my $gaps_query=0;
my $gaps_subject=0;
 
my $counter=5;
 
if (open(INPUT1, $input1))
{
                     
#  Database: /lustre/scratch115/teams/soranzo/projects/MPRA/reference_files/ensem
#bl/GRCh38_r90.all.fa
#    Posted date:  Feb 2, 2018  3:06 PM
#  Number of letters in database: 3,088,286,401
 # Number of sequences in database:  25



#Matrix: blastn matrix 1 -2
#Gap Penalties: Existence: 0, Extension: 2.5
#BLASTN 2.7.1+


#Reference: Zheng Zhang, Scott Schwartz, Lukas Wagner, and Webb
#Miller (2000), "A greedy algorithm for aligning DNA sequences", J
#Comput Biol 2000; 7(1-2):203-14.



#Database: /lustre/scratch115/teams/soranzo/projects/MPRA/reference_files/ensem
#bl/GRCh38_r90.all.fa
#25 sequences; 3,088,286,401 total letters



#Query= rs9374081_6_109310329_2

#Length=501
#                                                                     Score     E
#Sequences producing significant alignments:                          (Bits)  Value

#6  dna:chromosome chromosome:GRCh38:6:1:170805979:1 REF               920     0.0



#Query_1  1          CCAATCTTGCTACTGAAACCTCTCTCATCAAGAGATGGTCTTTGTCACTAAATCCTCTGG  60
#6        109310329  ............................................................  109310388

#Query_1  61         GTACTTTCAGTCACCTTTCTACCTAACCCCTTGGCAGCATTTGACACAGAGTCCTTCTCA  120
#6        109310389  ............................................................  109310448

#Query_1  121        AAGCCTCTCTGCCCTTGGCTTCTACAGCATGAGTCTCTCCTAGCTTTCCTCCAACATTGC  180
#6        109310449  ............................................................  109310508

#Query_1  181        TGAGCCCCTCTTCTCAGACCCTTCTGTAGGGTCCTCTACCTCTGCCAACCCTAAAAtttt  240
#6        109310509  ............................................................  109310568

#Query_1  241        tttATAATCTGTCCCACCACTGCTACCAGGGCCACATTTTGACTTCCGCAGTCCCTAACA  300
#6        109310569  ..........A.................................................  109310628

#Query_1  301        CTTTTACCTTTGTGGGCTCTTCCCTCCTCCATGAAAAACCATTAAAAATTATATTTTAAC  360
#6        109310629  ............................................................  109310688


#Lambda      K        H
 #  1.33    0.621     1.12

#Gapped
#Lambda      K        H
  # 1.28    0.460    0.850

#Effective search space used: 1457670839072
                     
                    while(my $line=<INPUT1>)
		    {
			chomp($line);
			$counter++;
                                unless($line !~ /\w/)
                                {
                                 #  print "LINE:$line:LINE\n";
                                    # Last line condition
                                     
                                    if($line=~/^Effective search space used/)
                                    {
                                        if($Flag_hit == 1)
                                        {
					    # Define the POS arrays

					    my @positions_query=();
					    my @positions_subject=();

					    my @Line_tmp = sort{$a<=>$b} keys %{$hash_alignment_parsing{'QUERY'}};
					    #print "The_array_of_line_indexes_is\t".join('~~',@Line_tmp)."\n";
					    
					    # Start moving through the array of line indexes

					    for(my $i=0;$i<scalar(@Line_tmp);$i++)
					    {
						# My query variables

						my $query="NaN";
						my $query_seq="NaN";
						my $query_begin="NaN";
						my $query_end="NaN";
						
						my @query_seq_tmp_faux=sort keys %{$hash_alignment_parsing{'QUERY'}{$Line_tmp[$i]}{'seq'}};
						
						$query_seq=join('',@query_seq_tmp_faux);
						
						foreach my $query_tok(sort keys %{$hash_alignment_parsing{'QUERY'}{$Line_tmp[$i]}{'coords'}})
						{
						    $query=$query_tok;
						    
						    foreach my $query_begin_tok(sort{$a<=>$b} keys %{$hash_alignment_parsing{'QUERY'}{$Line_tmp[$i]}{'coords'}{$query_tok}})
						    {
							$query_begin=$query_begin_tok;
							foreach my $query_end_tok(sort{$a<=>$b} keys %{$hash_alignment_parsing{'QUERY'}{$Line_tmp[$i]}{'coords'}{$query_tok}{$query_begin_tok}})
							{
							    $query_end=$query_end_tok;
							}
						    }

						}
						
#						print "----------------------------->$query\t$query_begin\t$query_end\n";
#						print "----------------------------->$query_seq\n";

						# My subject variables

						my $subject="NaN";
						my $subject_seq="NaN";
						my $subject_begin="NaN";
						my $subject_end="NaN";
						
						my @subject_seq_tmp_faux=sort keys %{$hash_alignment_parsing{'SUBJECT'}{$Line_tmp[$i]}{'seq'}};
						
						$subject_seq=join('',@subject_seq_tmp_faux);
						
						foreach my $subject_tok(sort keys %{$hash_alignment_parsing{'SUBJECT'}{$Line_tmp[$i]}{'coords'}})
						{
						    $subject=$subject_tok;
						    
						    foreach my $subject_begin_tok(sort{$a<=>$b} keys %{$hash_alignment_parsing{'SUBJECT'}{$Line_tmp[$i]}{'coords'}{$subject_tok}})
						    {
							$subject_begin=$subject_begin_tok;
							foreach my $subject_end_tok(sort{$a<=>$b} keys %{$hash_alignment_parsing{'SUBJECT'}{$Line_tmp[$i]}{'coords'}{$subject_tok}{$subject_begin_tok}})
							{
							    $subject_end=$subject_end_tok;
							}
						    }

						}
						
#						print "$subject\t$subject_begin\t$subject_end\n";
#						print "$subject_seq\n";

						######################################## Here we start going thorigh the sequences to see where the differences are
						my @nucleotide_query_tmp=split("",$query_seq);
						my @nucleotide_query_tmp_shifted=@nucleotide_query_tmp;

						my @nucleotide_subject_tmp=split("",$subject_seq);
						my @nucleotide_subject_tmp_shifted=@nucleotide_subject_tmp;

						my $length_query_alignment=scalar(@nucleotide_query_tmp);
						my $length_subject_alignment=scalar(@nucleotide_subject_tmp);

						my @passed_query_array=();
						my @passed_subject_array=();

						# Subject POS
						my $sbj_intraline_ABS_POS=$subject_begin;
#						print "INITIAL\t$sbj_intraline_ABS_POS\n";
						
						 while(scalar(@nucleotide_subject_tmp_shifted) !=0)
						 {
						     my $query_contender=shift(@nucleotide_query_tmp_shifted);
						     my $subject_contender=shift(@nucleotide_subject_tmp_shifted);
						     $subject_contender=uc($subject_contender);
						     $query_contender=uc($query_contender);

						     push(@passed_query_array,$query_contender);
#						     push(@passed_subject_array,$subject_contender);

						     # We extract a compare nucleotides from each hash

						     if($subject_contender eq '.')
						     {
							 push(@passed_subject_array,$query_contender);
							 # LAST THING TO DO add one for the next line

							 $sbj_intraline_ABS_POS++;
						     }
						     elsif($subject_contender eq '-' && $query_contender eq '-')
						     {
							 # Do nothing, these postions are the result of Multiple Alignment and do not exist in query and subject
						     }
						     elsif($subject_contender ne '-' && $query_contender eq '-')
						     {
							 # GAP in my ELEMENT, record position
							 
							 push(@passed_subject_array,$subject_contender);

							 my $adjusted_POS=$sbj_intraline_ABS_POS-1;
                                                         my $subject_first_NT_GAP=$passed_query_array[scalar(@passed_query_array)-2];
							 push(@{$hash1{$subject}{'del'}{$adjusted_POS}{$query}{'REF'}},$subject_first_NT_GAP);
                                                         push(@{$hash1{$subject}{'del'}{$adjusted_POS}{$query}{'ALT'}},$query_contender);

#							 print "hash1.1{$subject\t'del'\t$adjusted_POS\t$query\t'REF'\t".','."$subject_first_NT_GAP\n";
#							 print "hash1.1{$subject\t'del'\t$adjusted_POS\t$query\t'ALT'\t".','."$query_contender\n";
							 
							 # LAST THING TO DO add one for the next line
							 $sbj_intraline_ABS_POS++;

						     }
						     elsif($subject_contender eq '-' && $query_contender ne '-')
						     {
							 # GAP in GRCh38, record position
							 push(@passed_subject_array,$subject_contender);
                                         
							 my $adjusted_POS=$sbj_intraline_ABS_POS-1;
							 my $subject_first_NT_GAP=$passed_query_array[scalar(@passed_query_array)-2];
							 
							 # Array no hash, take into account possible repetitions

							 push(@{$hash1{$subject}{'in'}{$adjusted_POS}{$query}{'REF'}},$subject_first_NT_GAP);
							 push(@{$hash1{$subject}{'in'}{$adjusted_POS}{$query}{'ALT'}},$query_contender);
#							 $hash1{$subject}{'in'}{$adjusted_POS}{'Element'}{$query}=1;

#							 print "hash1.2\t$subject\t'in'\t$adjusted_POS\t'REF'\t$subject_first_NT_GAP\t$query\n";
#							 print "hash1.2{$subject\t'in'\t$adjusted_POS\t'ALT'\t$query_contender\t$query";

						     }
						     elsif($subject_contender ne '.' && $subject_contender ne '-' && $query_contender ne '-')
						     {
							 push(@passed_subject_array,$subject_contender);
							 # 'Missense' between ENSEMBL and UNIPROT order; subject->query

							 push(@{$hash1{$subject}{'mm'}{$sbj_intraline_ABS_POS}{$query}{'REF'}},$subject_contender);
                                                         push(@{$hash1{$subject}{'mm'}{$sbj_intraline_ABS_POS}{$query}{'ALT'}},$query_contender);

#							 print "hash1.3\t$subject\t'mm'\t$sbj_intraline_ABS_POS\t$query\t'ALT'\t".','."$subject_contender\n";
#							 print "hash1.3\t$subject\t'mm'\t$sbj_intraline_ABS_POS\t$query\t'ALT'\t".','."$query_contender\n";
							 
                                                         # LAST THING TO DO add one for the next line
							 $sbj_intraline_ABS_POS++;
#							 exit;

						     }


						 }# while2

						

					    }# array of line indexes

#					exit;##################################################################################################### EXIT HERE!!

                                        }# If hit ==1
                                        elsif($Flag_hit == 0)
                                       {
                                            print "WARNING_1Flag_hit == 0!!!!$query\t$subject\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\n";
                                            #~ print "$query\t$subject\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\n";
#                                            exit;
                                       }
                                        # Reinitialize variables
                                         
                                        $query="NaN";
                                        $length_query="NaN";
                                        $subject="NaN";
                                        $length_subject="NaN";
                                        $FLAG_line=0;
                                        $Flag_hit=1;
                                        %hash_alignment_parsing=();
                                        @seq_query=();
                                        @seq_subject=();
                                         
                                        $uncover_N_TER=0;
                                        $uncover_C_TER=0;
                                        $coverage=0;
                                        $mismatch=0;
                                        $identity=0;
                                        $gaps_query=0;
                                        $gaps_subject=0;
                                    }   
                                    elsif($line =~ /\*\*\*\*\* No hits found \*\*\*\*\*/)
                                    {
                                        $Flag_hit=0;
                                    }
                                    elsif($line=~/^Query= (.+)/)
                                    {
                                        $query=$1;
                                    }
                                    elsif($line=~/^Length=(.+)/)
                                    {
                                        $length_query=$1;
                                    }
                                    if($Flag_hit != 0)
                                    {
                                        if($line=~/^Query_[^\s]+\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)/)
                                        {
					    $counter=1;
                                            $FLAG_line++;
                                            my $begin=$1;
                                            my $seq=$2;
                                            my $end=$3;

					    $hash_alignment_parsing{'QUERY'}{$FLAG_line}{'seq'}{$seq}=1;
                                            $hash_alignment_parsing{'QUERY'}{$FLAG_line}{'coords'}{$query}{$begin}{$end}=1;

 #                                           print "hash_alignment_parsing{'QUERY'\t$FLAG_line\t'seq'\t$seq\n";
  #                                          print "hash_alignment_parsing{'QUERY'\t$FLAG_line\t'coords'\t$query\t$begin\t$end\n";



                                        }
                                        elsif($counter == 2 && $line=~/^([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)/)
                                        {
                                            my $subject_contender=$1;
                                            my $begin=$2;
                                            my $seq=$3;
                                            my $end=$4;
                                            #~ print "$begin\t$seq\t$end\n";
                                            $hash_alignment_parsing{'SUBJECT'}{$FLAG_line}{'seq'}{$seq}=1;
					    $hash_alignment_parsing{'SUBJECT'}{$FLAG_line}{'coords'}{$subject_contender}{$begin}{$end}=1;
					    
#					    print "hash_alignment_parsing{'SUBJECT'\t$FLAG_line\t'seq'\t$seq\n";
#					    print "hash_alignment_parsing{'SUBJECT'\t$FLAG_line\t'coords'\t$subject_contender\t$begin\t$end\n";

                                        }
                                    }       
                                }
    }# while INPUT1
}else {print "impossible to open INPUT1\n";die;}
 

if(open(OUTPUT1,'>'.$output1))
{
   print OUTPUT1 "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    foreach my $CHROM_tok(sort keys %hash1)
    {
	foreach my $Key_tok(sort keys %{$hash1{$CHROM_tok}})
	{
	    
	    # INSERTIONS

	    if($Key_tok eq 'in')
	    {
		foreach my $POS_tok(sort{$a<=>$b} keys %{$hash1{$CHROM_tok}{$Key_tok}})
		{
		foreach my $Element_tok(sort keys %{$hash1{$CHROM_tok}{$Key_tok}{$POS_tok}})
		{
		 
#		    print "-->$CHROM_tok\t$POS_tok\n";
		    # Get the ALT alleles and the elements that are characterized by this insertion
		    
		    my @REF_tmp=@{$hash1{$CHROM_tok}{$Key_tok}{$POS_tok}{$Element_tok}{'REF'}};
		    my $REF=$REF_tmp[0];
		    my @ALT_tmp=@{$hash1{$CHROM_tok}{$Key_tok}{$POS_tok}{$Element_tok}{'ALT'}};
		    unshift(@ALT_tmp,$REF);
		    my $ALT=join('',@ALT_tmp);
		    
		    $hash2{$CHROM_tok}{$POS_tok}{$REF}{$ALT}{$Element_tok}=1;
		    print OUTPUT1 "$CHROM_tok\t$POS_tok\t\.\t$REF\t$ALT\t100\tPASS\t$Element_tok\n";
		    

		}# Element 
		}# POSITION PHENOMENON

	    }# INSERTIONS

    # MISMATCHS

	    if($Key_tok eq 'mm')
	    {
		foreach my $POS_tok(sort{$a<=>$b} keys %{$hash1{$CHROM_tok}{$Key_tok}})
		{
		foreach my $Element_tok(sort keys %{$hash1{$CHROM_tok}{$Key_tok}{$POS_tok}})
		{
		 
#		    print "-->$CHROM_tok\t$POS_tok\n";
		    # Get the ALT alleles and the elements that are characterized by this insertion
		    
		    my @REF_tmp=@{$hash1{$CHROM_tok}{$Key_tok}{$POS_tok}{$Element_tok}{'REF'}};
		    my $REF=join('',@REF_tmp);
		    my @ALT_tmp=@{$hash1{$CHROM_tok}{$Key_tok}{$POS_tok}{$Element_tok}{'ALT'}};
		    my $ALT=join('',@ALT_tmp);
		    
		    $hash2{$CHROM_tok}{$POS_tok}{$REF}{$ALT}{$Element_tok}=1;
		    print OUTPUT1 "$CHROM_tok\t$POS_tok\t\.\t$REF\t$ALT\t100\tPASS\t$Element_tok\n";
#		    print "$CHROM_tok\t$POS_tok\t$REF\t$ALT\t$Element_tok\n";
		    

		}# Element 
		}# POSITION PHENOMENON

	    }# MISMATCHS

	    if($Key_tok eq 'del')
	    {
		foreach my $POS_tok(sort{$a<=>$b} keys %{$hash1{$CHROM_tok}{$Key_tok}})
		{
		foreach my $Element_tok(sort keys %{$hash1{$CHROM_tok}{$Key_tok}{$POS_tok}})
		{
		 
#		    print "-->$CHROM_tok\t$POS_tok\n";
		    # Get the ALT alleles and the elements that are characterized by this insertion
		    
		    my @REF_tmp=@{$hash1{$CHROM_tok}{$Key_tok}{$POS_tok}{$Element_tok}{'REF'}};
		    
		    my @ALT_tmp=@{$hash1{$CHROM_tok}{$Key_tok}{$POS_tok}{$Element_tok}{'ALT'}};
		    my $ALT=join('',@ALT_tmp);		    
		    
		    push(@REF_tmp,$ALT);
		    my $REF=join('',@REF_tmp);

		    
		    $hash2{$CHROM_tok}{$POS_tok}{$REF}{$ALT}{$Element_tok}=1;
		    print OUTPUT1 "$CHROM_tok\t$POS_tok\t\.\t$REF\t$ALT\t100\tPASS\t$Element_tok\n";
		    

		}# Element 
		}# POSITION PHENOMENON

	    }# DELETIONS

	}

    }

}
                     
sub timestamp {
    my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}
