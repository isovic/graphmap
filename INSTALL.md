## 1. Installation

You will need a recent GCC/G++ version (>=4.7) to compile the source.  

To override the default compiler choice you can set GCC (or GCC_MAC on Mac), e.g.:  

```
GCC=/usr/local/bin/g++ make  
```

### 1.1 Initialize submodules
This will automatically initialize/pull the latest version of submodules.  
```
make modules  
```

Submodules are used as source files, so there is no need to pre-compile them in any way.  


### 1.2 Linux
For a Linux release version type:
```
make  
```

To clean, type:
```
make clean  
```

One can also rebuild, which will cause clean and make to be ran sequentially:
```
make rebuild  
```

### 1.3 Mac
```
make mac  

make cleanmac  
make rebuildmac  
```

### 1.4. Compiling the debug version
```
make debug  

make cleandebug  
make rebuilddebug  
```
