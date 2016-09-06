declare -x SMRT_HOME='/pbi/analysis/smrtportal/beta/install/smrtanalysis_2.3.0.140936'
declare -x SUBMIT_CMD='qsub -V -cwd -pe smp 8 -l h_rt=240:0:0,mem_free=8G'

# common SGE qsub options:
# -V: sync enviroment
# -cwd: use current working directory
# -pe: choose parallel enviroment, use qconf -spl to list all pe
# -q: choose a queue, use qconf -sql to list all queue
# -l: set various limits