#!/usr/bin/env bash
#
# Pattern Recognition in Neuroimaging
# Practical Lab
#
# Ahmed Abdulkadir and Jonas Richiardi
# for the Medical Imaging Processing Lab at EPFL
# Prof. Dimitri van de Ville
#
# Downloads data sets and software from the institutional web page.
#
# Following archives will be downloaded:
# - pr_data1.zip
# - pr_data2.zip
# - pr_spm8.zip
# - pr_pronto.zip
#
# All the archives will be extracted into /tmp in order to assure
# compatibility and avoid permission problems.

if [ -d /tmp/AUDIO ]
then
    echo "Folder AUDIO already exists but will be removed!"
    echo ""
    echo "removing /tmp/AUDIO"
    rm -rf /tmp/AUDIO
fi

if [ -d /tmp/VISIO ]
then
    echo "Folder AUDIO already exists but will be removed!"
    echo ""
    echo "removing /tmp/AUDIO"
    rm -rf /tmp/VISIO
fi


echo "Linking MATLAB root into /tmp"
NAS1=/net/icitnass1
NAS2=/net/icitnass2
MATLABDIR=/export/software/MATLAB

if [ -d ${NAS1}${MATLABDIR} ]; then
    rm -rf /tmp/MATLAB
    ln -s ${NAS1}${MATLABDIR}/R2010b /tmp/MATLAB
elif [ -d ${NAS2}${MATLABDIR} ]; then
    rm -rf /tmp/MATLAB
    ln -s ${NAS2}${MATLABDIR}/R2010b /tmp/MATLAB
else
    echo
    echo "Error: Could not find MATLAB."
    echo
fi

echo "Changing to /tmp folder"
cd /tmp
echo "Getting SPM8"
wget http://miplab.epfl.ch/teaching/micro-513/download/glm_spm8.zip
echo "Getting PRoNTo"
wget http://miplab.epfl.ch/teaching/micro-513/download/glmprlab_pronto.zip
wget http://miplab.epfl.ch/teaching/micro-513/download/glmprlab_pronto476.zip
echo "Getting single-subject data set (AUDIO)"
wget http://miplab.epfl.ch/teaching/micro-513/download/glmprlab_data1.zip
echo "Getting multi-subjects data set (VISIO)"
wget http://miplab.epfl.ch/teaching/micro-513/download/glmprlab_data2.zip
echo "Extracting SPM8"
unzip glm_spm8.zip > /dev/null && rm -f glm_spm8.zip
echo "Extracting PRoNTo"
unzip glmprlab_pronto476.zip > /dev/null && rm -f glmprlab_pronto476.zip 
unzip glmprlab_pronto.zip > /dev/null && rm -f  glmprlab_pronto.zip
echo "Extracting AUDIO data set"
unzip glmprlab_data1.zip > /dev/null && rm -f   glmprlab_data1.zip 
echo "Extracting VISIO data set"
unzip glmprlab_data2.zip > /dev/null && rm -f glmprlab_data2.zip

echo ""
echo "Download and Extraction Done!"