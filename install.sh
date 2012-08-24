#!/bin/bash


if ls $R_LIBS/00LOCK-chifit >| /dev/null 2>&1
then
  rm -r $R_LIBS/00LOCK-chifit
fi

cd ${HOME}/andreas/dev


R --vanilla <<EOF
install.packages("chifit",repos=NULL)
EOF
