#!/bin/bash
# Usage ./check-sl-status.sh logs/FILE
display_usage() {
  echo -e "\nUsage: $0 CLUSTERLOG \n"
}

if [ $# -eq 0 ]
then
  display_usage
  exit 1
fi

if [[ ( $1 == "--help" ) || $1 == "-h" ]]
then
  display_usage
  exit 0
fi

grep -P "Estimating \d+..." $1 | tail
