#!/bin/sh
env -i 'DYLD_LIBRARY_PATH'='' \
'PATH'='/usr/bin' \
'XAUTHORITY'='/var/run/lightdm/mmarcano/xauthority' \
'PYTHONPATH'='/usr/share/esoreflex-2.9.1/esoreflex/python' \
'ESOREX_CONFIG'='/etc/esoreflex-esorex.rc' \
'ESOREX_RECIPE_CONFIG'='/etc/esoreflex_default_recipe_config.rc' \
'USER'='mmarcano' \
'LANG'='en_US.utf8' \
'DISPLAY'=':0' \
'HOSTNAME'='trantor.ttu.edu' \
'LD_LIBRARY_PATH'='' \
'LOGNAME'='mmarcano' \
'PWD'='/home/mmarcano/Documents/VIMOS/NGC6652/data/NGC6652B1/rawesodata' \
'HOME'='/home/mmarcano' \
'SHLVL'='1' \
'_'='/usr/bin/java' \
esorex --suppress-prefix=TRUE  --output-dir=/home/mmarcano/Documents/VIMOS/NGC6652/data/NGC6652B1/rawesodata/reflex_tmp_products/vimos-ifu/vmifustandard_1/2018-08-31T13:24:56.233 --log-dir=/home/mmarcano/Documents/VIMOS/NGC6652/data/NGC6652B1/rawesodata/reflex_logs/vimos-ifu/vmifustandard_1/2018-08-31T13:24:56.233 --recipe-config=/home/mmarcano/Documents/VIMOS/NGC6652/data/NGC6652B1/rawesodata/reflex_book_keeping/vimos-ifu/vmifustandard_1/2018-08-31T13:24:56.233/vmifustandard.rc --products-sof=products_sof.json vmifustandard /home/mmarcano/Documents/VIMOS/NGC6652/data/NGC6652B1/rawesodata/reflex_book_keeping/vimos-ifu/vmifustandard_1/2018-08-31T13:24:56.233/data.sof

