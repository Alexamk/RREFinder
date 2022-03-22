#!/bin/sh

unset HHLIB

# Remove all references to the script folder from the path
NEW_PATH=''
while IFS=':' read -ra ADDR; do
  for i in "${ADDR[@]}"; do
    if ! [ "$i" =  "$RRE_SCRIPTS_FOLDER" ]
    then
        NEW_PATH="$NEW_PATH:$i"
    fi
  done
done <<< "$PATH"
# Remove the first colon
export PATH=${NEW_PATH:1}

unset NEW_PATH
unset RRE_SCRIPTS_FOLDER
