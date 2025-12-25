# Graph Perturbation Analysis for Subgraph Counting

## Requirement
Cmake >= 2.8 GCC = 8.4.0

## Datasets
2 small datasets (DBLP and Human) are provided in this repo. Full datasets can be downloaded [here.](https://hkustconnect-my.sharepoint.com/:u:/g/personal/ssunah_connect_ust_hk/EWwS7ixh4NBHriiPHNpUMAkBu8vbH1f37Ug8CPWQdUXj4w?e=0GXEMg)

## Compile
Under the root directory of the project, execute the following commands to compile the source code.

```zsh
mkdir build
cd build
cmake ..
make
```
## Execute
After compiling the source code, you can find the binary file 'GPESCO.out' under the 'build/matching' directory.

-q, The query graph folder

-d, The data graph file path

-crq, Algorithm for CRQ (4 choices: JointEnum, LookAheadEnum, RWC, RWCH)

-drq, Algorithm for DRQ (4 choices: ProbeEnum, IEEnum, RWD, RWDW)

-ql, The query list file path (including all filenames of queries to execute)

-time_limit, Time limit for execution

-k, The k value for TOPKSUB.


Execute the binary with the following command (example):
```zsh
cd build/matching
./GPESCO.out -d ../../dataset/dblp/data_graph/dblp.graph  -q ../../dataset/dblp/query_graph -crq JointEnum -drq IEEnum -ql ../../dataset/qlist.txt -time_limit 600 -k 5
```
Dataset: DBLP, CRQ: JointEnum, DRQ: IEEnum, Time Limit: 600 seconds

