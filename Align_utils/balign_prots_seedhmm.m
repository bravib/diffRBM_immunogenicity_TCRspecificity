function alignc_new = balign_prots_seedhmm(seqs_seed, lea, lepmin, lepmax, seqs_new, weights)

disp(lea)
disp(lepmin)
disp(lepmax)
%xfold=['/home/barbara/Barbara_Bravi/diffrbm_selection/Selection_model/emerson_training_test/HMM_model_' sprintf('%d', lea)];
xfold=['/home/barbara/Barbara_Bravi/diffrbm_selection/Selection_model/Pgen/HMM_model_' sprintf('%d', lea)];

%%%% Build HMM profile from seqs_seed %%%%
rootf = '/home/barbara/Barbara_Bravi';

yes_weight = 0;
if nargin > 5 & weights
 weights = importdata('weights.txt');
 Weights = weights';
 yes_weight = 1;
end

yes_red = 0;
sca=log(2);
%sca=1;
SM = 'BLOSUM62'; % default BLOSUM50
GO = 8; % default 8 for gap penalty

Xdseq = seqs_seed;
M=length(seqs_seed);
if size(Xdseq) == [1 M],
 Xdseq = reshape(seqs_seed, [M 1]);
end

l=zeros(M,1);
lold=0;
   for i=1:M, 
     l(i)=size(Xdseq{i,1}, 2); 
 end
lmax=max(l); 
lmin=min(l);
disp(lmax)
disp(lmin)


af = cell2mat(Xdseq);
N = length(af(1,:));
alignc = zeros(M,N);
for i=1:M
    for j=1:N
       alignc(i,j) = letter2number(af(i,j));
    end
end

alignc_seed = alignc;
af_seed = af;
a_new = af;
a_new2= af;
a_new_repr = a_new2;

A = importdata([xfold '/MatchEmission.txt']);
B = importdata([xfold '/InsertEmission.txt']);
C = importdata([xfold '/NullEmission.txt']);
D = importdata([xfold '/BeginX.txt']);
E = importdata([xfold '/MatchX.txt']);
F = importdata([xfold '/InsertX.txt']);
G = importdata([xfold '/DeleteX.txt']);
H = importdata([xfold '/FlankingInsertX.txt']);
I = importdata([xfold '/LoopX.txt']);
L = importdata([xfold '/NullX.txt']);
hmmodel = hmmprofstruct(lea,'MatchEmission',A,'InsertEmission',B,'NullEmission',C,'BeginX',D,'MatchX',E,'InsertX',F,'DeleteX',G,'FlankingInsertX',H,'LoopX',I,'NullX',L);

Xdseq = seqs_new;

M=length(seqs_new);
disp(M)
if size(Xdseq) == [1 M],
 Xdseq = reshape(seqs_new, [M 1]);
end
l=zeros(M,1);
lold=0;
   for i=1:M, 
     l(i)=size(Xdseq{i,1}, 2);   % Length distribution of new sequences
 end
lmax=max(l); 
lmin=min(l);


indexl=cell(lmax,1);
dist=zeros(lmax,1);
for ll=1:lmax
 indexl{ll}=find(l==ll); % indices of seqs. in each length group
 dist(ll)=length(indexl{ll}); % population of each length group
end

al=cell(lmax,1);
p=cell(lmax,1);
for ll=1:lmax 
al{ll} = cell2mat(Xdseq(indexl{ll},1)); % converts a cell array into an ordinary array
p{ll}  = seqprofile(Xdseq(indexl{ll},1),'gaps','all','counts',true); % seq profiles (i.e. aa count in each position) for each length
end

% align new sequences to the HMM model
LA = size(a_new2,1);
a_new3 = a_new2;
%a_new_repr = a_new2; % probably the most consistent choice is to align to the full final sample?

%scores_new0 = [];
%scores_new0_cons = [];
if lmin <= lmax % 'align', using the HMM probabilistic model, only if the sequences differ from the average length of the profile

for ll = lmin:lmax
newseq = al{ll,1}; 
 if ll ~= lea
	for ns=1:dist(ll)
		[scores,aligned_seqs]=hmmprofalign(hmmodel,newseq(ns,:));
		[zz,indal]=hmmprofmerge(aligned_seqs);
		ali_seq=aligned_seqs(indal);
                %disp('i am here')

		a_new3=[a_new3;ali_seq]; % new alignment
		%scores_new0 = [scores_new0, scores];
		%scores_new0_cons = [scores_new0_cons, nwalign(seqconsensus(a_new_repr), newseq(ns,:), 'Scale', log(2))*length(newseq(ns,:))];
		%scores_new0_cons = [scores_new0_cons, nwalign(seqconsensus(a_new_repr), newseq(ns,:), 'Scale', sca, 'ScoringMatrix', SM, 'GapOpen', GO)];
%{		
		rf=[];
		for u = 1:size(a_new_repr,1),
			rf = [rf, nwalign(a_new_repr(u,:), newseq(ns,:), 'Scale', sca, 'ScoringMatrix', SM, 'GapOpen', GO)];
		end
		scores_new0_cons = [scores_new0_cons, max(rf)];
%}
	end
else

for ns=1:dist(ll)
	   ali_seq=newseq(ns,:);
  		a_new3=[a_new3;ali_seq]; % new alignment
		%[scores,~]=hmmprofalign(hmmodel,ali_seq);
		%scores_new0 = [scores_new0, scores];
		%scores_new0_cons = [scores_new0_cons, nwalign(seqconsensus(a_new_repr), ali_seq, 'Scale', sca, 'ScoringMatrix', SM, 'GapOpen', GO)];

		end
 end
end

else
ll = lmin;

if ll ~= lea
	for ns=1:dist(ll)
		[scores,aligned_seqs]=hmmprofalign(hmmodel,newseq(ns,:));
		[zz,indal]=hmmprofmerge(aligned_seqs);
		ali_seq=aligned_seqs(indal);
		a_new3=[a_new3;ali_seq]; % new alignment
		%scores_new0 = [scores_new0, scores];
		%disp(seqconsensus(a_new2))
		%scores_new0_cons = [scores_new0_cons, nwalign(seqconsensus(a_new_repr), newseq(ns,:),'Scale', sca, 'ScoringMatrix', SM, 'GapOpen', GO)];
	end
else
for ns=1:dist(ll)
	   ali_seq=newseq(ns,:);
		
	   a_new3=[a_new3;ali_seq]; % new alignment
	   %[scores,~]=hmmprofalign(hmmodel,ali_seq);
	     %scores_new0 = [scores_new0, scores];
	     %scores_new0_cons = [scores_new0_cons, nwalign(seqconsensus(a_new_repr), ali_seq, 'Scale', sca, 'ScoringMatrix', SM, 'GapOpen', GO)];

	end
end
end
a_new4 = a_new3((LA + 1) : (LA + length(seqs_new)),:);
cumdist = cumsum(dist);
af = repmat('', lea, size(a_new4,1));
%scores_new = zeros(size(a_new4,1),1);
%scores_new_cons = zeros(size(a_new4,1),1);
af(indexl{lmin,1},:)=a_new4(1:dist(lmin),:);
%scores_new(indexl{lmin,1}) = scores_new0(1:dist(lmin));
%scores_new_cons(indexl{lmin,1}) = scores_new0_cons(1:dist(lmin));
for ll=lmin+1:lmax
   af(indexl{ll,1},:) = a_new4(cumdist(ll-1)+1:cumdist(ll-1)+dist(ll),:); % final profile re-ordered
   %scores_new(indexl{ll,1}) = scores_new0(cumdist(ll-1)+1:cumdist(ll-1) + dist(ll)); % final HMM optimal scores re-ordered
   %scores_new_cons(indexl{ll,1}) = scores_new0_cons(cumdist(ll-1)+1:cumdist(ll-1) + dist(ll)); % final alignment scores
end

N = lea;
alignc = zeros(M-2,N);
for i=1:M
    for j=1:N
       alignc(i,j) = letter2number(af(i,j)); % convert into numbers to count them and obtain frequencies
    end
end
alignc_new = alignc;
end
