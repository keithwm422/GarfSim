#!/bin/bash

debug_mode=0  #0=normal mode, 1=debug mode
nlimit=20

root="root -q -l -b"

function execute()
{
  if [ $debug_mode -eq 0 ]; then
    echo -e $C_BLU$*$C_RET
    eval $*
  else
    echo -e $C_GRN$*$C_RET
  fi  

  return $?
}

function error()
{
  echo -e $C_RED"Error: $*"$C_RET
}

function makedir()
{ 
  # If target directory exist already.
  # mode:0 skip
  # mode:1 remove
  # mode:2 move
  
  dir=$1
  mode=$2
  
  if [ ! -d "$dir" ]; then
    execute "mkdir -p $dir"
  elif [ $mode -eq 1 ]; then
    execute "rm -rf $dir"
    execute "mkdir -p $dir"
  elif [ $mode -eq 2 ]; then
    execute "mv $dir $dir.org"
    execute "mkdir -p $dir"
  fi  
}

function check()
{
  ##----- Check No. of Process -----
  while true
  do  
    nproc0=`ps -ef | grep converter.exe | grep -v grep | wc -l`
    nproc1=`ps -ef | grep root.exe | grep proc | grep -v grep | wc -l`
    nproc=$(( $nproc0+$nproc1 ))
    if [ $nproc -le $nlimit ]; then
      echo "Process [$nproc jobs]   `date`"
      echo $stg1f
      execute "sleep 2"
      break
    else
      echo "Process [$nproc jobs]   `date`" 
      execute "sleep 2"
    fi  
  done
}

##### main process

C_RED="\033[1;31m"
C_GRN="\033[1;32m"
C_BLU="\033[1;34m"
C_RET="\033[1;39m"

#target_dir="data/stg0"

check
execute "./isochronGenerator_p_fullEfield 14.616 293 &"
wait

