# Manuel Tardaguila, PhD. University of Florida.

#!/usr/bin/perl


use strict;
use warnings;
use Time::localtime;

my $input1=$ARGV[0];
my $reference=$ARGV[1];
my $output=$ARGV[2];


my %hash0=();
my %hash1=();
my %hash1_rev=();
my %hash2=();
my %hash3=();

 

#~ exit;

my $time='['. timestamp(). ']'."\n";
print "Queries\n";

my $Gene_ID2="NaN";



if (open(INPUT1, $input1))
{

# /nfs/users/nfs_m/mt19/RareVar_2019/MIPRA_design

#>Element_7;TWO_THIRD;REF
#ACTTTGAAGGAAACACACAGTAGCTTAATCATAGCACAGAGATGATAACCTCATGGTTTTCTCTTCTTCCTGGGTCTCAGGAAGCGACGGCTACTCCACCCTTCCCAGTGGTTGCTGGCAGAAAACATTGATGTCGATTAGTATATAGGG
#>Element_7;TWO_THIRD;ALT
#ACTTTGAAGGAAACACACAGTAGCTTAATCATAGCACAGAGATGATAACCTCATGGTTTTCTCTTCTTCCTGGGTCTCAGGAAGCGACGGCTACTCCACCGTTCCCAGTGGTTGCTGGCAGAAAACATTGATGTCGATTAGTATATAGGG
#>Element_7;HALF;REF
#AATCATAGCACAGAGATGATAACCTCATGGTTTTCTCTTCTTCCTGGGTCTCAGGAAGCGACGGCTACTCCACCCTTCCCAGTGGTTGCTGGCAGAAAACATTGATGTCGATTAGTATATAGGGTATATACAGAGAGTGATGGAAAAGAA

	while(my $line=<INPUT1>)
	{
	    chomp($line);
	    $line =~ s/\r|\n//g;
#	    print "LINE:$line:LINE\n";
	    
	    if($line =~/^>(.+)/)
	    {
		$Gene_ID2 = $1;
	    }
	    else
	    {
		my $seq=$line;
		$hash1{$Gene_ID2}{$seq}=1;
		print "hash0\t$Gene_ID2\t$seq\n";
#		exit;
	    }
	}
}else{print "Unable to open files\n";die;}

#exit;

# Once created I am not going to do it more


my @Gene_tmp=reverse sort keys%hash1;

my $global_file=$output; #join($experiment,"global_result.txt","_");
unlink($global_file);
$global_file= $output;



for(my $i =0 ;$i < scalar(@Gene_tmp); $i++)
{

 #   print "Hello_world_2:\t$Gene_tmp[$i]\n";
#exit;

	my $query_file= "query$i.tmp";
	my $result_file= "result$i.tmp";
	
	if(open(QUERY,'>'.$query_file) && open(RESULT,'>'.$result_file))
	{
		foreach my $query_seq_tok(sort keys %{$hash1{$Gene_tmp[$i]}})
		{
			print QUERY ">$Gene_tmp[$i]\n";
			print QUERY "$query_seq_tok\n";

		#	print "Hello_world_3>$Gene_tmp[$i]\n";
		#
               #	print "Hello_world_3\t$query_seq_tok\n";
			#exit;
		}
		#exit;
		
		####
		print "Blasting\t:$Gene_tmp[$i]\t$i\t$time\n";
		
		system("blastn -db $reference -query $query_file -num_alignments 1 -evalue 5e-90  -outfmt 3  > $result_file");		
	}
	 
		
	system("cat $result_file >> $global_file");
	unlink($query_file,$result_file);
}



sub timestamp {
  my $t = localtime;
  return sprintf( "%04d-%02d-%02d_%02d-%02d-%02d",
                  $t->year + 1900, $t->mon + 1, $t->mday,
                  $t->hour, $t->min, $t->sec );
}

