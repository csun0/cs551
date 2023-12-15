N=6
(
for f in *.pdb
do
	((i=i%N)); ((i++==0)) && wait
	/mnt/f/linux/rosetta/main/source/bin/rosetta_scripts.static.linuxgccrelease \
		-parser:protocol relax_struct.xml \
		-in:file:s $f \
		-nstruct 3 \
		-relax:jump_move true \
		-overwrite &
done
)
