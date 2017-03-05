./main.o fg.trr input.inp cgtypes.inp cg.dat

fg.trr: input gromacs fine-grained trajectory

input.inp:  
#1    total_number_of_cg_sites
#2..  ;atom_id weight_in_cg cg_site_id
(both atom_id cg_site_id  start from one)

cgtypes.inp
#1    number_of_cg_types
#2..  ;cg_site_id  cg_type
(cg_site_id starts from 1)
