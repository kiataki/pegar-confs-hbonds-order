#!/bin/csh

set exec = /home1/kiataki/meusprogramasfortran/pegar-conf-hbonds-novo_v2teste/pegarconf-hb-novo_v2teste.x

$exec << fim > confs_hb.com
XYZFILE = out_1m4n_all.xyz 
GAUFILE = cabecalho 
NSOLUTO = 14
NWATERS = 3 
NCONNFS = 250 
NPCHARG = 1500 
fim

