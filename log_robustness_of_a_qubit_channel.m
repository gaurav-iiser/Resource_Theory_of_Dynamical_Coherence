addpath(genpath('QETLAB-0.9'))

%{
code to calculate the log-robustness of a given channel
define a channel N or its Choi matrix J_N as given in the paper 'Dynamical Resource Theory of Quantum Cohernce' by Gaurav Saxena, Eric Chitambar, and Gilad Gour.
This code is written by Gaurav Saxena.
%}
%{
This is an SDP problem whose primal problem is the log of:
1) Find: min (1/|A0|) Tr[w_A]
2) subject to:
    w_A >= J_N
    diag(diag(w_A) = w_A
    w_{A0} = Tr[w_A]u_{A_0}
    w_A >= 0

and dual problem is the log of:
1) Find: max Tr[n J_N]
2) subject to:
    diag(diag(n)) = diag(diag(n_{A0})) otimes u_{A1}
    diag(diag(n_{A1})) = I_{A1}
    n_A >= 0
    
%}
%{
Note: This code finds the log-robustness of qubit channels only.


Requirements:
1) Matlab
2) Qetlab
3) CVX
%}


clear all
fileID = fopen('newdata.txt','w')

clear B

%%{
for count = 1:1
    
    
    %%%%%%% defining J_N %%%%%%%%%%%%%%%
    %J_N = RandomSuperoperator(2,1,1, 1);
    J_N = 0.5*[1 1 1 -1; 1 1 1 -1; 1 1 1 -1; -1 -1 -1 1]; %Choi of Hadamard
    %trace(J_N)
    %{
       if PartialTrace(J_N,1) ~= eye(2)
        count
        eig(J_N)
        eig(PartialTrace(J_N,1))
    end
    %}
    %J_N = Tensor(eye(2),RandomDensityMatrix(2,1));                         %%%%% Choi matrix of the replacement channel with a random density matrix as output %%%%%%%%%%%
    %{
                                                                            %%%%%          Choi matrix of a replacement channel that outputs a plus state       %%%%%%%%%%%
    plus_state = (1/sqrt(2))*[1; 1];
    plus_state_density_mat = plus_state * plus_state';
    J_N = Tensor(eye(2), plus_state_density_mat ); %Replacement Channel (outputs a plus state)
    %}
    %zero_zero = [1 0; 0 0];
    
    %J_N =2 * MaxEntangled(2) * MaxEntangled(2)';                           %%%%%                      Choi matrix of the identity channel                       %%%%%%%%%%%
    %J_N = RandomSuperoperator(2, 1, 1, 1,1);                             %%%%%%                     Choi matrix of a unitary channel                          %%%%%%%%%%%
    %trace(J_N)
    %{
                                                                            %%%%%%                     Choi matrix of dephasing channel                          %%%%%%%%%%%
    J_N = zeros(4);
    J_N(1,1) = 1;
    J_N(4,4) = 1;
    %}
    %{
                                                                            %%%%%%                    Choi matrix of classical channel                           %%%%%%%%%%
    J_N = zeros(4);
        random_num1 = rand;
        random_num2 = rand;
        J_N(1,1) = random_num1;
        J_N(2,2) = 1 - random_num1;
        J_N(3,3) = random_num2;
        J_N(4,4) = 1 - random_num2;
    %}
    %{
                                                                            %%%%%%                    Choi matrix of depolarizing channel                        %%%%%%%%%%
    p = rand;
    J_N = p * eye(4)/2 + (1-p)*[1 0 0 1; 0 0 0 0 ; 0 0 0 0 ; 1 0 0 1];
    %}
    
    %J_N = RandomSuperoperator([2,2],1,1,0,2);                               %%%%%%                     Choi matrix of a random unital channel with 2 Kraus operations  %%%%%%%%%%%%%%
    %J_N = RandomSuperoperator(2,1,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u_A1 = eye(dim) / trace( eye(dim) );
    u_A0 = eye(dim) / trace( eye(dim) );
    
    
    cvx_begin sdp quiet
    dimension = dim;
    variable w(dimension^2,dimension^2) semidefinite
    minimize (trace(w)/trace( eye(dim)) )
    subject to
    w >= J_N;
    diag(diag(w)) == w;
    PartialTrace(w) == trace(w)*u_A0;
    cvx_end
    
    
    LR_N = log2(cvx_optval);
    if LR_N > 1.8
        LR_N
        J_N
    end
    B(count) = LR_N;
    
end
%}

%%%%%%         dual %%%%%%%%%%%%%
%{
for count = 1:20

    %%%%%%% defining J_N %%%%%%%%%%%%%%%
    %J_N = 0.5*[1 1 1 -1; 1 1 1 -1; 1 1 1 -1; -1 -1 -1 1]; %Choi of Hadamard
    %J_N = RandomSuperoperator(2, 1, 0,  1);
%{
        if PartialTrace(J_N) ~= eye(2)
        count
        eig(J_N), eig(PartialTrace(J_N))
    end
%}
    %J_N = Tensor(eye(2),RandomDensityMatrix(2,1));                         %%%%% Choi matrix of the replacement channel with a random density matrix as output %%%%%%%%%%%
%{
                                                                            %%%%%          Choi matrix of a replacement channel that outputs a plus state       %%%%%%%%%%%
    plus_state = (1/sqrt(2))*[1; 1];
    plus_state_density_mat = plus_state * plus_state';
    J_N = Tensor(eye(2), plus_state_density_mat ); %Replacement Channel (outputs a plus state)
%}
    %zero_zero = [1 0; 0 0];
    
    %J_N =2 * MaxEntangled(2) * MaxEntangled(2)';                           %%%%%                      Choi matrix of the identity channel                       %%%%%%%%%%%
    %J_N = RandomSuperoperator(2, 1, 1, 1,1)/2;                             %%%%%%                     Choi matrix of a unitary channel                          %%%%%%%%%%%
%{
                                                                            %%%%%%                     Choi matrix of dephasing channel                          %%%%%%%%%%%
    J_N = zeros(4);
    J_N(1,1) = 1;
    J_N(4,4) = 1;
%}
%{
                                                                            %%%%%%                    Choi matrix of classical channel                           %%%%%%%%%%
    J_N = zeros(4);
        random_num1 = rand;
        random_num2 = rand;
        J_N(1,1) = random_num1;
        J_N(2,2) = 1 - random_num1;
        J_N(3,3) = random_num2;
        J_N(4,4) = 1 - random_num2;
%}
%{
                                                                            %%%%%%                    Choi matrix of depolarizing channel                        %%%%%%%%%%
    p = rand;
    J_N = p * eye(4)/2 + (1-p)*[1 0 0 1; 0 0 0 0 ; 0 0 0 0 ; 1 0 0 1];
%}
    
    %J_N = RandomSuperoperator([2,2],1,1,0,2);                               %%%%%%                     Choi matrix of a random unital channel with 2 Kraus operations  %%%%%%%%%%%%%%
    %J_N = RandomSuperoperator(2,1,1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u_A1 = eye(dim) / trace( eye(dim) );
    u_A0 = eye(dim) / trace( eye(dim) );
    %count = 0;

    cvx_begin sdp quiet
        dimension = dim;
        variable n(dimension^2,dimension^2) semidefinite
        maximize (trace(n * J_N) )
        subject to
            diag(diag(n)) == Tensor( diag( diag(PartialTrace(n) ) ),u_A1 );
            diag(diag(PartialTrace(n,1))) == eye(dim);
    cvx_end
    
    LR_N = log2(cvx_optval);
    
    B(count) = LR_N
    
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = sort(B,'descend');
for i = 1:count
    fprintf(fileID,'%d \n',  B(i));
end
fclose(fileID);


%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Function Definitions  %
%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = dim
d = 2;
end


