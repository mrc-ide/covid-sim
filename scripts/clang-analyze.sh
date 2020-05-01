#!/bin/sh -e

for x in "$@";
do
  if [ "${x}" == "-c" ];
  then
    temp_out="`mktemp /tmp/clang-analyze.out.XXXXX`"
    temp_bin="`mktemp /tmp/clang-analyze.bin.XXXXX`"

    all_options="$@"
    if [ "${all_options}" == "-v" ];
    then
      all_options=""
    fi

    all_options=`echo "${all_options}" | sed 's,^-v ,,g'`
    all_options=`echo "${all_options}" | sed 's, -v$,,g'`
    all_options=`echo "${all_options}" | sed 's, -v , ,g'`

    # analyze
    # -w : ignore regular compiler warnings so 'temp_out' only contains
    # messages from analyzer. '-w' should not be a part of toolchain flags
    # since Hunter toolchain-id calculated using '#pragma message' output which
    # is implemented as a warning message (hence will be suppressed by '-w')
    clang --analyze -w ${all_options} -o "${temp_bin}" 2> "${temp_out}"

    RESULT=0
    [ "$?" == 0 ] || RESULT=1
    [ -s "${temp_out}" ] && RESULT=1

    cat "${temp_out}";
    rm -f "${temp_out}"
    rm -f "${temp_bin}"

    if [ "${RESULT}" == "1" ];
    then
      exit 1;
    fi
  fi
done

# compile real code
clang "$@"
