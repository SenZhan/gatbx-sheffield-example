% This is one simple piece of  code for solving
% De Jong's first test function.
% It aims at practicing gatbx for genetic algorithm with matlab.
% Source: GAToolbox_Documentation
% De Jong's first test function:
% f(x) = Sum(x(i)^2), -512 le x(i) le 512,
% with i goes from 1 all the way to n.
% Here we may as well set n = 20.
% This problem is a little bit trivial, apprantly set all x(i) = 0
% and then f(x) = 0. But it is good to start with.
tic
NIND = 40; %# of individuals
MAXGEN = 5000; % max # of generations
NVAR = 20; % # of variables
PRECI = 20; % precision
GGAP = 0.9; % generation gap

% Build field descriptor
FieldD = [ rep( [ PRECI ], [ 1, NVAR ] ); rep( [ -512; 512 ], [ 1, NVAR ] );...
    rep( [ 1; 0; 1; 1 ], [ 1,NVAR ] ) ];
% Initialize population
Chrom = crtbp( NIND, NVAR * PRECI );
gen = 0; % count generations
ObjV = objfun1 ( bs2rv ( Chrom, FieldD ) ); %objfun1 is offered
% by the gatbx toolbox which is specifically oriented to 
% this question.

% Generation loop
while gen < MAXGEN
    FitnV = ranking ( ObjV ); % rank ObjV for fitness degree(0-2)
    Selch = select ( 'sus', Chrom, FitnV, GGAP ); % select not strongest to mate,
    % rest directly go to next generation
    Selch = recombin( 'xovsp', Selch, 0.7); % recombined
    Selch = mut( Selch ); %mutated
    ObjVSel = objfun1( bs2rv( Selch, FieldD ) );
    [ Chrom, ObjV ] = reins ( Chrom, Selch, 1, 1, ObjV, ObjVSel ); % This is the reforming of Chrom  
    % related ObjV
    gen = gen + 1;
end
Phen = bs2rv (Chrom, FieldD);
toc    
    

