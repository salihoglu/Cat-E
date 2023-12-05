#!/bin/bash

if [ ! -f /tmp/.X1-lock ]; then
    Xvfb :1 -screen 0 1920x1080x24+32 -fbdir /var/tmp > /tmp/Xvfb.log 2>&1 &
fi

export DISPLAY=:1
echo $DISPLAY


R -e "shiny::runApp('/app', host='0.0.0.0', port=3838)"

