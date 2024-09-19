**module load ..**
shovill/1.1.0--0
seroba/1.0.0=py36_1-c1
themisto/2.1.0
mgems/1.3.1--h468198e_0
msweep/2.0.0
alignment-writer/0.4.0
poppunk/2.6.0

**Step 1**: Create mash sketch (i.e. reference) for each GPSC isolate assembly (uses all lanes)
./run_mash_sketch.sh \
    ./lanes.txt \
    ./assemblyfind_lanes.txt \
    ./mash_sketches_ref

**Step 2**:  Run mash paste (i.e. combine mash references)
./run_mash_paste.sh ./mash_sketches_ref ./combined_ref.msh

**Step 3**:  Run mash screen for study plate sweep samples (i.e. align deep sequences from study plate sweep samples against mash reference)
./run_mash_screen.sh ./lanes_study.txt ./combined_ref.msh ./mash_output

**Step 4**: Combine isolate assemblies to make one reference FASTA (only lanes that have passed QC*)
./run_combine_fastas.sh ./lanes.txt ./assemblyfind_lanes.txt ./combined_GPS.fna

**Step 5**: Build Themisto index on this reference
./run_themisto_build.sh \
    ./combined_GPS.fna \
    ./themisto_index

**Step 5.1**: Running PopPUNK2
Copy the reference database and cluster file to your directory
./PopPUNK2/GPS_v7 and ./GPS_v7_external_clusters.csv

Copy seroba database to your directory
./seroba/database 

Prepare a 2-column tab-delimited file: 1st sample name 2nd path to the assembly file
query.txt
pf assembly -i lanes.txt -t file > assembly.txt
Create query.txt in excel using the assembly.txt

poppunk_assign --db GPS_v6 --distances GPS_v6/GPS_v6.dists --query query.txt --output GPSC_assignment --external-clustering GPS_v6_external_clusters.csv

**Step 6**: Run themisto align and mSWEEP for:
./run_msweep_pipeline.sh ./lanex_study.txt ./lanes.txt ./combined_GPS.fna ./themisto_index ./GPSC_assignment/GPSC_assignment_external_clusters_new.csv ./seroba/database ./msweep_output/
