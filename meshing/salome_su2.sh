#!/bin/bash

#----------Takes a dumped python script from salome and spits out an su2 mesh by using the cdfmsh.py library-------

if [ $# -lt 1 ]; then
    echo "Provide a salome python file as input"
    exit 1
fi

salome_py_file=$1
output_name=${salome_py_file%".py"}
tmp_file="generate_su2.py"

echo "#!/bin/python3" > $tmp_file

tail -n +2 $salome_py_file | cat >> $tmp_file

echo "from cfdmsh import *" >> $tmp_file
echo "ExportSU2File(mesh=Mesh_1,file=\"$output_name\")" >> $tmp_file

chmod u+x $tmp_file
./$tmp_file #running the script

rm $tmp_file