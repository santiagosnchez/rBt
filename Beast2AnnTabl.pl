#!/usr/bin/perl
# Beast2AnnTbl.pl
# ©Santiago Sanchez-Ramirez, University of Toronto

my $usage = "
perl Beast2AnnTbl.pl <your_mcc.tre>

example:
perl Beast2AnnTbl.pl My_Analyses/dataset1.mcc.tre | tee my_annotations.txt

The script will generate a tab-delimited table.
Tips are numbered from 1 to N terminals and nodes 
are numbered from N+1 to N+(N-1). They are ordered
cladewise and should match the node order in APE (R).

email santiago.snchez@gmail.com for details

";

if (!$ARGV[0]){
	die $usage;
}
elsif ($ARGV[0] =~ m/^-*he{0,1}l{0,1}p{0,1}$/){
	die $usage;
}

my @nexus=();
while(<>){
	chomp($_);
	push @nexus, $_;
}

my ($nexTree) = grep { /^\t*tree/ } @nexus;
$nexTree =~ s/tree TREE1 \= \[&R\] //;

my @edges = split(/:/, $nexTree);

my @annot=();
my @phylo=();
my @nodesNtips=();

foreach (@edges){
	$_ =~ m/\[.+?\]/smp;
	my $ann = ${^MATCH};
	my $pre = ${^PREMATCH};
	$pre =~ m/[\(,](\d+)$/;
	my $nt = $1;
	push @annot, $ann;
	push @phylo, $pre;
	if ($nt){
		push @nodesNtips, $nt;
	} else {
		push @nodesNtips, "*";
	}
}

my @backbone1 = grep { /[\(\),]/ } split(//, join(':', @phylo));
my @backbone2 = grep { /[\(\)]/ } split(//, join(':', @phylo));
my $Nnodes = scalar(grep { /\)/ } @backbone1);
my $Ntips = scalar(grep { /,/ } @backbone1) + 1;

my @Nodeind = grep { $nodesNtips[$_] eq '*' } 0 .. $#nodesNtips;

my $ind = 0;
my $currnode = $Ntips;
my @cladeNodes=();
my @open=();
for my $i (0 .. $#backbone2){
	BLOCK: {
		if ($backbone2[$i] eq '('){
			++$currnode;
			if (grep { /^$currnode$/ } (@open, @cladeNodes)){
				#++$currnode;
				redo BLOCK;
			}
			push @open, $currnode;
		}
		if ($backbone2[$i] eq ')'){
			if (grep { /^$currnode$/ } @cladeNodes){
				--$currnode;
				redo BLOCK;
			}
			push @cladeNodes, $currnode;
		}
	}
}

@nodesNtips[@Nodeind] = @cladeNodes;
my %dataord=();
map { $dataord{$nodesNtips[$_] } = $annot[$_] } 0 .. $#annot;

my %annotData=();
for $x (@nodesNtips){
	%annotData = getannot(\$dataord{$x}, \$x, \%annotData);
}

my %getNames=();
map { map { $getNames{$_} = "" } keys %{$annotData{$_}} } keys %annotData;
my @annotNames = sort {$a cmp $b} keys %getNames;

print join("\t", ("node",@annotNames)) . "\n";
for $node (sort {$a <=> $b} keys %annotData){
	print $node . "\t";
	for $ann (@annotNames){
		if ($ann ne $annotNames[$#annotNames]){
			if ($annotData{$node}{$ann}){
				print $annotData{$node}{$ann} . "\t";
			} else {
				print "NA" . "\t";
			}
		} else {
			if ($annotData{$node}{$ann}){
				print $annotData{$node}{$ann} . "\n";
			} else {
				print "NA" . "\n";
			}
		}
	}
}

sub getannot {
	my ($ann,$node,$hash) = @_;
	my @annotTab = split /=/, $$ann;
	my @curr=();
	my @prev=();
	my @next=();
	my %hash=();
	for my $i (0 .. $#annotTab){
		if ($i == 0){
			$annotTab[$i+1] =~ s/[\{\}\]\[]//g;
			@next = split /,/, $annotTab[$i+1];
			$annotTab[$i] =~ s/\[\&//;
			if (scalar(@next) > 2){
				pop @next;
				$$hash{$$node}{$annotTab[$i]} = join(',', @next);
			} else {
				$$hash{$$node}{$annotTab[$i]} = $next[0];
			}
		} elsif ($i != 0 and $i != $#annotTab){
			$annotTab[$i] =~ s/[\{\}\]\[]//g;
			$annotTab[$i+1] =~ s/[\{\}\]\[]//g;
			$annotTab[$i-1] =~ s/[\{\}\]\[]//g;
			@curr = split /,/, $annotTab[$i];
			@next = split /,/, $annotTab[$i+1];
			@prev = split /,/, $annotTab[$i-1];
			if (scalar(@next) > 2){
				pop @next;
				$$hash{$$node}{$curr[$#curr]} = join(',', @next);
			} else {
				$$hash{$$node}{$curr[$#curr]} = $next[0];
			}
		}
	}
	return(%$hash);
}



