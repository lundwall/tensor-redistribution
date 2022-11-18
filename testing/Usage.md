### Modules
```
env2lmod
module load gcc/8.2.0 openmpi/4.0.2 cmake/3.20.3 python/3.8.5
```

### Single block tests
First use
```
./euler_run_all_dim.sh
```
which will generate a new testing result directory. Then
```
python handle_single_block_results.py [generated directory]
```