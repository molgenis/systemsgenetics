#!/bin/sh

MYGROUPS=$(eval 'echo ${GROUPS_'`whoami`'-undefined}')
NOREV=0000000000000000000000000000000000000000

while read old new ref ; do

	# To push to master...
	if test $ref = refs/heads/master ; then
		#.... you need to be the hudson user
		if test `whoami` != hudson ; then
			echo "`whoami` is not authorized to push to master"
			exit 1
		fi

		# Exit now to avoid error message about naming conventions
		exit 0
	fi

	# if this is a branch with a prefixed name...
	if echo $ref | grep -q "^refs/heads/.*-" ; then
		# do not allow the user to push if they do not have a group set
                 MYGROUPS=`whoami`

		# do not allow the user to push if the branch is not in their group set
		if echo $ref | grep -q "^refs/heads/${MYGROUPS}-" ; then
			:
		elif test `whoami` != hudson ; then
			# give hudson the ultimate permission 
			echo "`whoami` is not authorized to push to $ref"
			echo "Or $ref does not exist! "
                        echo "Please contact your adminstrator to fix this."
			exit 1
		fi
	else # branch does not have a prefix on the form 'prefix-*'
	  	echo "$ref is not a valid branch name. Please consult the naming conventions."
	  	exit 1;
	fi

done

exit 0