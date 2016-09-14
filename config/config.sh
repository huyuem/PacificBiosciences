# directories, user specific
declare -x SMRT_HOME='/pbi/analysis/smrtportal/beta/install/smrtanalysis_2.3.0.140936'
# the following directories need write permission of the first R library and genome user.
declare -x RLIBDIR='/home/bhan/Rlib'
declare -x ANNOTATION_DIR='/home/bhan/annotation/'
# declare -x ANNOTATION_DIR='/home/bhan/jobs/tmp'

# SGE related
declare -x SUBMIT_CMD='qsub -V -cwd -pe smp 4 -l h_rt=240:0:0,mem_free=8G'
# common SGE qsub options:
# -V: sync enviroment
# -cwd: use current working directory
# -pe: choose parallel enviroment, use qconf -spl to list all pe
# -q: choose a queue, use qconf -sql to list all queue
# -l: set various limits
