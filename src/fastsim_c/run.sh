clear


echo "removing old fastsim"
rm ./compiled_fastsim
echo "REMOVED!\n"

echo "compiling fastsim"
gcc fastsim.c -o compiled_fastsim -lm
echo "COMPILED!\n"

echo "running fastsim\n"
./compiled_fastsim
echo "\n\n--\nDONE!\n"