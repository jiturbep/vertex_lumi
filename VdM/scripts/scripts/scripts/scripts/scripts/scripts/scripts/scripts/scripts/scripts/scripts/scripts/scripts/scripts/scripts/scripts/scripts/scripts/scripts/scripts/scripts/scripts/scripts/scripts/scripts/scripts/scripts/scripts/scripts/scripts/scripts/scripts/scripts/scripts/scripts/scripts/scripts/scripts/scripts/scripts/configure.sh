cd RootCore
./configure
cd ../
source RootCore/scripts/setup.sh

${ROOTCOREBIN}/scripts/find_packages.sh
${ROOTCOREBIN}/scripts/compile.sh

echo "All Done."
