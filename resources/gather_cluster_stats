#!/bin/bash
## AUTHOR : Vishal N. Koparde, Ph.D., CCBR, NCI
## DATE : Mar 2021
## This scripts gathers cluster related statistics for jobs run on Biowulf using Snakemake by:
##   > extracting "external" jobids from snakemake.log files
##   > gather cluster stats for each job using sacct command
##   > sorts output by job submission time
##   > output TSV file
##


if [ "$#" != "1" ];then
	echo " bash $0 <snakemake log file> "
	exit 1
fi


function get_sacct_info {
	jobid=$1
	attribute=$2
	x=$(sacct -j $jobid --noheader --format="${attribute}%500"|head -n1|awk '{print $1}')
	echo $x
}

function displaytime {
  local T=$1
  local D=$((T/60/60/24))
  local H=$((T/60/60%24))
  local M=$((T/60%60))
  local S=$((T%60))
  (( $D > 0 )) && printf '%d-' $D
  printf '%02d:%02d:%02d' $H $M $S
}

function get_batchline_variable {
	variable=$1
	grep -m1 "^$1" /dev/shm/${jobid}.sacct.batchline|awk '{print $NF}'
}

function get_secondline_variable {
	variable=$1
	grep -m1 "^$1" /dev/shm/${jobid}.sacct.secondline|awk '{print $NF}'
}

function print_jobid_stats {
jobid=$1
declare -A jobdataarray
export SLURM_TIME_FORMAT="%s" # epoch time
sacct -P -j $jobid --units=G --format="JobID,AveRSS,MaxRSS,ReqMem,NCPUS,AllocCPUS,ReqCPUS,Start,End,Submit,Elapsed,Timelimit,Account,State,Group,GID,JobName,NodeList,Partition,QOS,ReqTRES,User,UID,WorkDir" > /dev/shm/${jobid}.sacct
nlines=$(cat /dev/shm/${jobid}.sacct|wc -l)
jobdataarray["jobid"]="$jobid"
if [ "$nlines" == "2" ];then
	jobdataarray["submit_time"]="JOBNOTACCOUNTABLE"
else
	head -n2 /dev/shm/${jobid}.sacct | \
	awk -F "|" '
{
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' > /dev/shm/${jobid}.sacct.secondline
	grep -m2 "JobID\|batch" /dev/shm/${jobid}.sacct | \
	awk -F "|" '
{
    for (i=1; i<=NF; i++)  {
        a[NR,i] = $i
    }
}
NF>p { p = NF }
END {
    for(j=1; j<=p; j++) {
        str=a[1,j]
        for(i=2; i<=NR; i++){
            str=str" "a[i,j];
        }
        print str
    }
}' > /dev/shm/${jobid}.sacct.batchline
#batch line variables

	jobdataarray["elapsed"]=$(get_batchline_variable "Elapsed")
	jobdataarray["reqcpus"]=$(get_batchline_variable "ReqCPUS")
	jobdataarray["ncpus"]=$(get_batchline_variable "NCPUS")
	jobdataarray["alloccpus"]=$(get_batchline_variable "AllocCPUS")
	jobdataarray["avg_mem"]=$(get_batchline_variable "AveRSS")
	jobdataarray["max_mem"]=$(get_batchline_variable "MaxRSS")
	jobdataarray["req_mem"]=$(get_batchline_variable "ReqMem")

#second line variables

	jobdataarray["submit_time"]=$(get_secondline_variable "Submit")
	jobdataarray["start"]=$(get_secondline_variable "Start")
	jobdataarray["end"]=$(get_secondline_variable "End")
	st=${jobdataarray["submit_time"]}
	jobdataarray["human_submit_time"]=$(date -d @$st|sed "s/ /_/g")
	jobdataarray["state"]=$(get_secondline_variable "State")
	qt=$(echo ${jobdataarray["start"]} ${jobdataarray["submit_time"]}|awk '{print $1-$2}')
	rt=$(echo ${jobdataarray["end"]} ${jobdataarray["submit_time"]}|awk '{print $1-$2}')
	# echo "TESTING"
	# echo jobid:${jobdataarray["jobid"]}
	# echo submit:${jobdataarray["submit_time"]}
	# echo start:${jobdataarray["start"]}
	# echo end:${jobdataarray["end"]}
	# echo elapsed:${jobdataarray["elapsed"]}
	# echo jobid:${jobdataarray["jobid"]}
	# echo "TESTING"
	# exit
	jobdataarray["queued"]=$(displaytime $qt)
	jobdataarray["runtime"]=$(displaytime $rt)
	jobdataarray["job_name"]=$(get_secondline_variable "JobName")
	jobdataarray["time_limit"]=$(get_secondline_variable "Timelimit")
	jobdataarray["node_list"]=$(get_secondline_variable "NodeList")
	jobdataarray["run_node_partition"]=$(get_secondline_variable "Partition")
	jobdataarray["qos"]=$(get_secondline_variable "QOS")
	jobdataarray["username"]=$(get_secondline_variable "User")
	jobdataarray["uid"]=$(get_secondline_variable "UID")
	jobdataarray["groupname"]=$(get_secondline_variable "Group")
	jobdataarray["gid"]=$(get_secondline_variable "GID")
	jobdataarray["account"]=$(get_secondline_variable "Account")
	jobdataarray["workdir"]=$(get_secondline_variable "WorkDir")
	jobdataarray["ReqTRES"]=$(get_secondline_variable "ReqTRES")

fi
echo -ne "${jobdataarray["submit_time"]}\t"
echo -ne "${jobdataarray["human_submit_time"]}\t"
echo -ne "${jobdataarray["jobid"]};${jobdataarray["state"]};${jobdataarray["job_name"]}\t"
echo -ne "${jobdataarray["node_list"]};${jobdataarray["run_node_partition"]};${jobdataarray["qos"]}\t"
echo -ne "${jobdataarray["queued"]};${jobdataarray["elapsed"]};${jobdataarray["runtime"]};${jobdataarray["time_limit"]}\t"
echo -ne "${jobdataarray["reqcpus"]};${jobdataarray["alloccpus"]}\t"
echo -ne "${jobdataarray["req_mem"]};${jobdataarray["max_mem"]}\t"
echo -ne "${jobdataarray["username"]};${jobdataarray["uid"]};${jobdataarray["groupname"]};${jobdataarray["gid"]};${jobdataarray["account"]}\t"
echo -ne "${jobdataarray["ReqTRES"]}\t"
echo -ne "${jobdataarray["workdir"]}\n"
}


snakemakelogfile=$1
externalidslst="/dev/shm/$(basename ${snakemakelogfile}).jobids.lst"
grep "with external jobid" $snakemakelogfile | awk '{print $NF}' | sed "s/['.]//g" | sort | uniq > $externalidslst
echo -ne "##SubmitTime\tHumanSubmitTime\tJobID:JobState:JobName\tNode;Partition:QOS\tQueueTime;RunTime;Q+R;TimeLimit\tReqCPU;AllocCPU\tReqMEM:MaxMem\tUsername:UID:Group:GID:Account\tReqTRES\tWorkdir\n"
while read jid;do
	print_jobid_stats $jid
done < $externalidslst |sort -k1,1n
rm -f $externalidslst /dev/shm/${jobid}.sacct*
