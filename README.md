# PyRAMA_multi
Python3 implementation of the Ramachandran plot. 

This is a new version of old PyRAMA which could process complexs with multiple chains, for my own ICA work. You could check each chain, with outliear report now.

Example: one of the chains of a complex:

![image](https://github.com/user-attachments/assets/f8e54cde-a9ca-4bfb-a2d6-17306d022f0c)



Current fork is still under development, so it can not be downloaded from Pypi. It may have some bug, so please contact me if you meet some. If you want to use, please direct clone the repository and use the script `core.py`:

```python
python pyrama/core.py target.pdb
```

Please change the DEFAULT_OUTPUT_DIR in pyrama/config.py to your own work dictory.

