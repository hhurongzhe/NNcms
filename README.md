# NNcms

This is a C++ 17 code, which calculates the chiral NN interactions under c.m. partial-wave basis.

I provide the functions and parameters for N2LO-EMN500 as an example. You can change for your own use of course!

## Structure

- data-cms/: dir for storing the interaction matrix, partial-waves table and momentum mesh.
- src/basic_math.hpp: basic math functions.
- src/configs.hpp: configuration structure for all parameters.
- src/constants.hpp: constants.
- src/gauss_legendre: computes gauss-legendre mesh points and weights. You can replace it.
- src/infile.hpp: used for read .ini file.
- src/interaction_aPWD.hpp: do partial-wave decomposition(PWD).
- src/interaction_part_contact.hpp: contact terms.
- src/interaction_part_pion_exchange.hpp: pion exchange terms.
- src/interaction_all.hpp: adding contact terms and pion exchange terms.
- src/lib_define.hpp: necessory libs.
- src/main.cpp: main function, calculating and writing to files.
- infile.ini: all parameters.
- Makefile: template makefile.

## Quick use

For a normal run, you should follow:

1. enter the NNcms/ dir
2. compile the codes using Makefile
3. make sure data-cms/ dir exists
4. edit parameters in infile.ini
5. edit table_momentum_mesh.txt
6. edit table_partial_waves.txt for your target partial-waves
7. edit table_tlabs.txt for your target Tlabs
8. run the NNcms.x
9. in the end you can see the result files in data-cms/

## Others

It is open for anyone to use.

If you have any further needs or questions, please contact me: rongzhe_hu@pku.edu.cn
