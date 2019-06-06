function [grade]=tester(data_filename)

%TESTER Tests student projects and returns score
%   TESTER(DATA_FILENAME) tests the receiver in
%   the current directory.  A file named Rx.m must exist, with
%   all the specifications given in the project description.
%
%   DATA_FILENAME is a string like 'easy' which contains
%   the received signal data to be decoded (i.e. it points
%   to the file 'easy.mat').
%
%   For example, tester('easy') tests the receiver in
%   the current directory on the data file 'easy.mat'.

% v1.0 AGK 11/05/2015 initial version

% number of symbols per user per frame
userDataLength=100;

% if data filename was not passed to function, ask for it
if nargin==0
    data_filename=input('Enter data file: ','s');
end

% load data file
load(data_filename)

% call receiver
[decoded_msg y]=Rx(r, rolloff, desired_user);

% plot constellation of soft decisions
figure
plot(y,'.')

% if student's decoded_msg is a column vector, make it a row vector
[a b]=size(decoded_msg);
if b==userDataLength
    decoded_msg=decoded_msg';
end
decoded_msg=decoded_msg(:)';

% sync up decoded_msg with m
% (but we can toss the first 5 frames since they're not graded on those)
cor_vec=conv(double(fliplr(m(desired_user,5*userDataLength+1:end))),double(decoded_msg));
cor_vec=cor_vec(size(m,2)-5*userDataLength:end);
[junk data_start]=max(cor_vec);

% convert from letters to bits
m2=text2bin(m(desired_user,5*userDataLength+1:end));
decoded=text2bin(decoded_msg(data_start:end));

% append Inf's to the end of student vector if it's too short
N=length(m2);
decoded=[decoded Inf*ones(1,N-length(decoded))];

% convert PAM symbols to bits (assumes Grey encoding)
orig_bits=m2;
student_bits=decoded(1:N);

% calculate number of bit errors
e=sum(abs(orig_bits-student_bits)>0);

% calculate bit error rate
BER=e/N;

% calculate grade
if BER<1e-3
    grade=100;
elseif BER>0.5
    grade=0;
elseif BER>0.25
    grade=round(150-300*BER);
else
    grade=round(69-10*log10(BER));
end

% output results
disp(['BER: ' num2str(BER)])
disp(['Grade: ' num2str(grade)]);
