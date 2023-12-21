#!/bin/bash
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/software/sse/easybuild/prefix/software/Anaconda3/2020.07-extras-nsc1/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/software/sse/easybuild/prefix/software/Anaconda3/2020.07-extras-nsc1/etc/profile.d/conda.sh" ]; then
        . "/software/sse/easybuild/prefix/software/Anaconda3/2020.07-extras-nsc1/etc/profile.d/conda.sh"
    else
        export PATH="/software/sse/easybuild/prefix/software/Anaconda3/2020.07-extras-nsc1/bin:$PATH"
    fi
fi
unset __conda_setup
 #<<< conda initialize <<<
