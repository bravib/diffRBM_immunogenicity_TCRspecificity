function pwmat = call_seqprofile_weighted(y)

We = load('weights_cef.txt');
q = y
Xdseq = importdata('cefs.txt');
af = cell2mat(Xdseq);

[pwmat] = seqprofile_rew(af,q,transpose(We))

namef=['/home/barbara/Barbara_Bravi/rbm/data/CEF/pwmR.txt'];
dlmwrite(namef, pwmat, 'delimiter' , ' ');
