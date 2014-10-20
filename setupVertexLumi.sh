echo "Setting up ATLAS"
setupATLAS
echo "Setting up athena 19.0.3"
asetup 19.0.3
echo "Setting up RootCore and compiling"
cd RootCore
./configure
source scripts/setup.sh
cd ..
rc compile