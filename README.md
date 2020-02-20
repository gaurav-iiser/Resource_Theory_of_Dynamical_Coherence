# Resource_Theory_of_Dynamical_Coherence
Contains all codes useful for reading/understanding, or producing examples from the paper 'Dynamical Resource Theory of Quantum Coherence' (https://arxiv.org/pdf/1910.00708). This link would be updated as soon as the paper is published in a journal.
Currently there are 4 codes available:
1.) miso_monotones.m - In our work, we provide a complete family of monotones for the case of Maximally Incoherent Superchannels (MISC). This code enables one to find the value of the monotone for any qubit channel.
2.) log_robustness_of_a_qubit_channel.m - Using this code, one can find the log-robustness of coherence of a qubit channel which is an operational monotone for dynamical coherence.
3.) interconversion_distance_with_varying_dim.m - This code calculates the interconversion distance from channel N to channel M. That is, if the result of this code is 0 (i.e., has a very small value of the order of ~10^(-9)) then, we can convert a channel N to M.
4.)interconversion_distance_qubit_channel.m - This code is a specific case for the above general code where the channels are qubit channels.
