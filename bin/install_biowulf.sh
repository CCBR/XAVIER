#!/usr/bin/env bash
# Usage:
#  bin/install_biowulf.sh
set -euxo pipefail

repo_path=/data/CCBR/projects/techDev/XAVIER
version=`cat ${repo_path}/VERSION`
install_path=/data/$USER/XAVIER/.${version}
bin_path=${install_path}/bin/

. "/data/CCBR_Pipeliner/db/PipeDB/Conda/etc/profile.d/conda.sh"
conda activate py311

# remove artifacts from prior builds
pushd ${repo_path}
rm -rf build/ *.egg-info
popd

echo "Installing XAVIER to ${install_path}"
pip install ${repo_path} --target ${install_path} --upgrade
chmod a+rx ${install_path}/xavier/bin/*.*
chmod -R a+r ${install_path}

if [[ ":$PATH:" != *":${bin_path}:"* ]];then
    export PATH="${PATH}:${bin_path}"
fi

if [[ ":$PYTHONPATH:" != *":${install_path}:"* ]];then
    export PYTHONPATH="${PYTHONPATH}:${install_path}"
fi