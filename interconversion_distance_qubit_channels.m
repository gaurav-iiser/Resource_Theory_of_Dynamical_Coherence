addpath(genpath('QETLAB-0.9'))
%{ 
Code to find the interconversion distance between two channels such that N
-> M as given in the paper 'Dynamical Resource Theory of Quantum Cohernce' by Gaurav Saxena, Eric Chitambar, and Gilad Gour.
This code is written by Gaurav Saxena.
This code is designed specifically for qubit channels. 
%}

%{
> Multiple Choi matrices are given below. 
> To use a particular Choi matrix, uncomment that line or section(for
instance, in the case of replacement channel outputting a plus state) and
comment other lines.
> Feel free to add any other Choi matrix for any channel.


Requirements:
1) Matlab
2) Qetlab
3) CVX
%}


clear all

J_N = 0.5*[1 1 1 -1; 1 1 1 -1; 1 1 1 -1; -1 -1 -1 1];        %%%%%% Choi of Hadamard <- Maximally Coherent Channel
for count = 1:100
%%%%%%% defining J_M %%%%%%%%%%%%%%%
    J_M = RandomSuperoperator(2,1,1, 1);
    %J_M = 0.5*[1 1 1 -1; 1 1 1 -1; 1 1 1 -1; -1 -1 -1 1]; %Choi of Hadamard
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

cvx_begin sdp quiet
dimension = dim;

variable w(dimension^2,dimension^2) semidefinite
variable alph(dimension^4,dimension^4) semidefinite
variable lambda nonnegative
minimize lambda 
subject to
        lambda * eye(dim) >= PartialTrace(w);
        w >= PartialTrace( alph * Tensor(transpose(J_N), eye(dim^2)) ,1) - J_M;
        PartialTrace(alph, [4],[dimension, dimension, dimension, dimension]) == Tensor( PartialTrace(alph, [2,4],[dimension, dimension, dimension, dimension]) , u_A1);
        PartialTrace(alph, [1,4],[dimension, dimension, dimension, dimension]) == eye(4);
        for f_1 = 1:dim
            for f_2 = 1:dim
                f = [ f_1  f_2 ];
                for i = 1:dim               %%for B_0
                    for j = 1:dim           %%for B_0
                       for u = 1:dim        %%for B_1
                            for v = 1:dim   %%for B_1
                                if (i ~= j) | (u ~= v)
                                    trace(alph * Tensor(A(f), basis(i,j), basis(u,v) ) ) == 0;
                                end
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
   d = dim;
   a = zeros(d^2, d^2);
   for k = 1: dim
      a = a + Tensor(basis(k,k),basis(f(k),f(k))); %function definition needs to go in
   end
end

function basis_elem = basis(i,j)
        basis_elem = Tensor(x(i), x(j)');
end

function basis_vector = x(a)
    basis_vector = zeros(dim,1);
    basis_vector(a) = 1;
end

function d = dim
    d = 2;
end

function d = dim_in
    d = 2;
end

function d = dim_out
    d = 2;
end

