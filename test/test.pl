my @queries = [1234, 65, 456];
my %honk = map {; queries => $_ } @queries;
print %honk;
