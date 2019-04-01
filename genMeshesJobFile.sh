NMESHES=2501
NELEMENTS=256
NEX1=7.8    #number of exclusion parameters
NEX2=0.2
RPARAMSLO=-5.53
RPARAMSHI=0.2
MARGIN_LO=-1
MARGIN_R=0.003
MARGIN_U=-1
MARGIN_LE=0.003
SIGMAGPR=0.4
LENGTHSCALE=0.08
LENGTHSCALER=0.05
SIGMOID=1.2
COORDINATEDIST='tiles'
RADIIDIST='logn'

CORES=1

DATESTR=`date +%m-%d-%H-%M-%N`	#datestring for jobfolder name
#Set up file paths
PROJECTDIR="/home/constantin/python/projects/stokesEquation"
JOBNAME="genMesh_coordinatedist=${COORDINATEDIST}_nEl=${NELEMENTS}N1=${NEX1}N2=${NEX2}margins=${MARGIN_LO}_${MARGIN_R}_${MARGIN_U}_${MARGIN_LE}_radiidist=${RADIIDIST}_r=${RPARAMSLO}_${RPARAMSHI}"
JOBDIR="/home/constantin/python/jobs/${JOBNAME}_${DATESTR}"
JOBSCRIPT="${JOBDIR}/genMesh_cluster.py"

#Create job directory and copy source code
rm -rf $JOBDIR
mkdir $JOBDIR
MESHGENSCRIPT="${PROJECTDIR}/genMesh_cluster.py"
cp  $MESHGENSCRIPT $JOBSCRIPT
#Change directory to job directory; completely independent from project directory
cd $JOBDIR

#construct job file
echo "#!/bin/bash" >> ./job_file.sh
echo "#SBATCH --job-name=${JOBNAME}" >> ./job_file.sh
echo "#SBATCH --partition batch_SNB,batch_SKL" >> ./job_file.sh
echo "#SBATCH --output=/home_eth/constantin/OEfiles/${JOBNAME}.%j.out" >> ./job_file.sh
echo "#SBATCH --error=/home_eth/constantin/OEfiles/${JOBNAME}.%j.err" >> ./job_file.sh
echo "#SBATCH --mincpus=${CORES}" >> ./job_file.sh
echo "#SBATCH --mail-type=ALL" >> ./job_file.sh
echo "#SBATCH --mail-user=mailscluster@gmail.com " >> ./job_file.sh
echo "#SBATCH --time=1000:00:00" >> ./job_file.sh

echo "sed -i \"14s/.*/nMeshes = $NMESHES/\" ./genMesh_cluster.py" >> ./job_file.sh
echo "sed -i \"15s/.*/nElements = $NELEMENTS  # PDE discretization/\" ./genMesh_cluster.py" >> ./job_file.sh
echo "sed -i \"21s/.*/nExclusionParams = [$NEX1, $NEX2]/\" ./genMesh_cluster.py" >> ./job_file.sh
echo "sed -i \"22s/.*/coordinateDist = '${COORDINATEDIST}'/\" ./genMesh_cluster.py" >> ./job_file.sh
echo "sed -i \"25s/.*/margins = ($MARGIN_LO, $MARGIN_R, $MARGIN_U, $MARGIN_LE)/\" ./genMesh_cluster.py" >> ./job_file.sh
echo "sed -i \"28s/.*/radiiDist = '${RADIIDIST}'/\" ./genMesh_cluster.py" >> ./job_file.sh
echo "sed -i \"29s/.*/r_params = [$RPARAMSLO, $RPARAMSHI]/\" ./genMesh_cluster.py" >> ./job_file.sh
echo "sed -i \"35s/.*/cov_l = ${LENGTHSCALE}/\" ./genMesh_cluster.py" >> ./job_file.sh
echo "sed -i \"36s/.*/sig_scale = ${SIGMOID}/\" ./genMesh_cluster.py" >> ./job_file.sh
echo "sed -i \"37s/.*/sigmaGP_r = ${SIGMAGPR}/\" ./genMesh_cluster.py" >> ./job_file.sh
echo "sed -i \"38s/.*/lengthScale_r = ${LENGTHSCALER}/\" ./genMesh_cluster.py" >> ./job_file.sh


#Activate fenics environment and run python
echo "source ~/.bashrc" >> ./job_file.sh
echo "conda activate mshr_new" >> ./job_file.sh
echo "ulimit -s 32000" >> ./job_file.sh
#echo "while true; do" >> ./job_file.sh
echo "python ./genMesh_cluster.py" >> ./job_file.sh
#echo "done" >> ./job_file.sh
echo "" >> ./job_file.sh


chmod +x job_file.sh
#directly submit job file
sbatch job_file.sh
#execute job_file.sh in shell
#./job_file.sh



