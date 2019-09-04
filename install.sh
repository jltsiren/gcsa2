#!/bin/bash
# ./install.sh [directory [jobs]]
# Install GCSA2 headers and library to the home directory or to the specified
# directory. The number of parallel make jobs is equal to the number of
# logical CPU cores unless specified otherwise.

GCSA_DIR=$(pwd)
INCLUDE_DIR=include
INCLUDE=gcsa
LIBRARY_DIR=lib
LIBRARY=libgcsa2.a

# Determine the directory where we will install GBWT.
if [ $# -ge 1 ]; then
  PREFIX=${1}
else
  PREFIX=${HOME}
fi
mkdir -p "${PREFIX}" 2> /dev/null
cd "${PREFIX}" > /dev/null 2>&1
if [ $? != 0 ]; then
	echo "Error: Directory '${PREFIX}' does not exist and cannot be created."
	exit 1
else
	PREFIX=$(pwd -P)
fi
echo "Installing GCSA2 to '${PREFIX}'."

# Determine the number of jobs.
if [ $# -ge 2 ]; then
  JOBS=${2}
else
  JOBS=$(getconf _NPROCESSORS_ONLN)
fi
echo "Running ${JOBS} parallel make jobs."

# Create the include/library directories.
mkdir -p "${INCLUDE_DIR}" 2> /dev/null
if [ $? != 0 ]; then
	echo "Error: Directory '${PREFIX}/${INCLUDE_DIR}' does not exist and cannot be created."
	exit 1
fi
mkdir -p "${LIBRARY_DIR}" 2> /dev/null
if [ $? != 0 ]; then
	echo "Error: Directory '${PREFIX}/${LIBRARY_DIR}' does not exist and cannot be created."
	exit 1
fi

# Build GCSA2.
cd "${GCSA_DIR}"
make clean > /dev/null
if [ $? != 0 ]; then
	echo "Error: Cleanup failed."
	exit 1
fi
make -j "${JOBS}"
if [ $? != 0 ]; then
	echo "Error: Could not compile GCSA2."
	exit 1
fi

# Install GCSA2.
cp -R "${INCLUDE_DIR}/${INCLUDE}" "${PREFIX}/${INCLUDE_DIR}" 2> /dev/null
if [ $? != 0 ]; then
	echo "Error: Could not install the headers."
	exit 1
fi
cp "${LIBRARY}" "${PREFIX}/${LIBRARY_DIR}" 2> /dev/null
if [ $? != 0 ]; then
	echo "Error: Could not install the headers."
	exit 1
fi

echo "GCSA2 installed to '${PREFIX}'."
