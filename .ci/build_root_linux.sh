pushd $DEPS_DIR 

wget -nv http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $DEPS_DIR/miniconda
export PATH="$DEPS_DIR/miniconda/bin:$PATH"
hash -r
conda config --add channels conda-forge
conda config --set channel_priority strict
# conda install --quiet --yes -c conda-forge/label/gcc8 root zstd 
conda create -n my_root_env root -c conda-forge/label/gcc8

# source "$DEPS_DIR/miniconda/bin/thisroot.sh"
# export CXX="$DEPS_DIR/miniconda/bin/g++"
popd
