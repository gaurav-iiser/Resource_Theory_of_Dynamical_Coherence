addpath(genpath('QETLAB-0.9'))
%{ 
Code to find the interconversion distance between two channels such that N
-> M as given in the paper 'Dynamical Resource Theory of Quantum Cohernce' by Gaurav Saxena, Eric Chitambar, and Gilad Gour.
->This code is written by Gaurav Saxena.
->This code is designed to handle any kind of channels. The only thing to
note is that for any channel, the number of inputs should be equal to the
number of outputs.
->The dimension of any system is the dimension of the density matrix of that
system.
%}

%{
> Multiple Choi matrices are given below. 
> To use a particular Choi matrix, uncomment that line or section(for
instance, in the case of replacement channel outputting a plus state) and
comment other lines.
> Feel free to add any other Choi matrix for any channel but be cautious of
the dimensions of the Choi matrix.


Requirements:
1) Matlab
2) Qetlab
3) CVX
%}

clear all

%J_N = 0.5*[1 1 1 -1; 1 1 1 -1; 1 1 1 -1; -1 -1 -1 1];        %%%%%% Choi of Hadamard <- Maximally Coherent Channel
%{                                                                         %%%%%          Choi matrix of a replacement channel that outputs a plus state       %%%%%%%%%%%
    plus_state = (1/sqrt(2))*[1; 1];
    plus_state_density_mat = plus_state * plus_state';
    J_N = Tensor(eye(2), plus_state_density_mat ); %Replacement Channel (outputs a plus state)
    J_N = Tensor(J_N,J_N); %Replacement channel that outputs two plus states
 %}
for count = 1:1
%%%%%%% defining J_M %%%%%%%%%%%%%%%
    %J_M = RandomSuperoperator(2,1,1, 1);
    J_M = 0.5*[1 1 1 -1; 1 1 1 -1; 1 1 1 -1; -1 -1 -1 1]; %Choi of Hadamard
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
u_A1 = eye(dim_in) / trace( eye(dim_in) );

cvx_begin sdp quiet
d_in = dim_in;
d_out = dim_out;
variable w(d_out^2,d_out^2) semidefinite
variable alph(d_in^2 * d_out^2, d_in^2 * d_out^2) semidefinite
variable lambda nonnegative
minimize lambda 
subject to
        lambda * eye(dim_out) >= PartialTrace(w);
        w >= PartialTrace( alph * Tensor(transpose(J_N), eye(d_out^2)) ,[1,2],[dim_in, dim_in, dim_out, dim_out]) - J_M;
        PartialTrace(alph, [4],[dim_in, dim_in, dim_out, dim_out]) == Tensor( PartialTrace(alph, [2,4],[dim_in, dim_in, dim_out, dim_out]) , u_A1);
        PartialTrace(alph, [1,4],[dim_in, dim_in, dim_out, dim_out]) == eye(dim_in*dim_out);
        for k = 1:dim_in^dim_in
            f = ones(1, dim_in);
            num = k;
            for l = dim_in:1
                [num,rem] = quorem(sym(num), sym(dim_in));
                f(1,l) = f(1,l) + rem;
            end
            
            for i = 1:dim_out               %%for B_0
                for j = 1:dim_out           %%for B_0
                   for u = 1:dim_out        %%for B_1
                        for v = 1:dim_out   %%for B_1
                            if (i ~= j) | (u ~= v)
                                trace(alph * Tensor(A(f), basis(i,j,'out'), basis(u,v,'out') ) ) == 0;
                            end
                        end
                   end
                end
            end   
        end
cvx_end
cvx_optval
end



%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  Function Definitions  %
%%%%%%%%%%%%%%%%%%%%%%%%%% 

function pauli_basis = Pauli_basis
    I = eye(2);
    sigma_1 = [0, 1; 1  0];
    sigma_2 = [0 -1i; 1i  0];
    sigma_3 = [1  0; 0 -1];
    pauli_basis = I;
    pauli_basis(:,:,2) = sigma_1;
    pauli_basis(:,:,3) = sigma_2;
    pauli_basis(:,:,4) = sigma_3;
end

function a = A(f)
   d = dim_in;
   a = zeros(d^2, d^2);
   for k = 1: d
      a = a + Tensor(basis(k,k,'in'),basis(f(k),f(k),'in')); %function definition needs to go in
   end
end

function basis_elem = basis(i,j,chr)
    if strcmp(chr, 'in')
        basis_elem = Tensor(x(i,'in'), x(j,'in')');
    else
        basis_elem = Tensor(x(i,'out'), x(j,'out')');
    end
end



function basis_vector = x(a, chr)
    if strcmp(chr, 'in')
        basis_vector = zeros(dim_in,1);
        basis_vector(a) = 1;
    else
        basis_vector = zeros(dim_out,1);
        basis_vector(a) = 1;
    end
end

function d = dim_in
    d = 4;
end

function d = dim_out
    d = 2;
end

