Installing PyCUTEst
===================

**Disclaimer: I wrote the instructions according to my own experience installing this thrice, I can't guarantee it will work out of the box for you**

Setup
-----

Choose a new directory to store all the files using:
```bash
export OPTDIR=$HOME/Documents/Optimization
mkdir $OPTDIR
```

Download CUTEst
---------------

Go into the new created directory:
```bash
cd $OPTDIR
```

Download ARCHDefs with the command

```bash
svn checkout --username anonymous \
  http://ccpforge.cse.rl.ac.uk/svn/cutest/archdefs/trunk ./archdefs
```

Download SIFDecode with the command
```
svn checkout --username anonymous \
  http://ccpforge.cse.rl.ac.uk/svn/cutest/sifdecode/trunk ./sifdecode
```

Download CUTEst, with the command
```bash
svn checkout --username anonymous \
  http://ccpforge.cse.rl.ac.uk/svn/cutest/cutest/trunk ./cutest
``` 
all in the same directory.

This will produce three sub-directories:

    ./archdefs
    ./sifdecode
    ./cutest
    
Add Files
---------
Start by adding the following enviromental variables:
```bash
export ARCHDEFS="$OPTDIR/archdefs"
export SIFDECODE="$OPTDIR/sifdecode"
export CUTEST="$OPTDIR/cutest"
```

You shoud posses two different files: ``pycutestitf.c`` ([here](https://gist.github.com/antonior92/df4ce3edcadfedc196fea04c488881ca)) and ``pycutestmgr.py`` ([here](https://gist.github.com/antonior92/b630fee6d98fdb54ebe6bde19847e1b5)). And you should move them to the 
following appropriate locations:
```bash
mv pycutestitf.c $CUTEST/src/tools
mkdir $CUTEST/src/python
mv pycutestmgr.py $CUTEST/src/python
```
    
Install CUTEst
---------------

Install cutest with the command:
```bash
cd $CUTEST
$ARCHDEFS/install_optsuite
```
and follow the instructions

Download Problems Set
---------------------

Download problem set from: http://www.cuter.rl.ac.uk/Problems/mastsif.shtml and save them on ``$OPTDIR/sif``.

Create Cache
------------

Create a cache file

```bash
export PYCUTEST_CACHE="$OPTDIR/cache"
mkdir $PYCUTEST_CACHE
```
    
Modify ``.bashrc``
------------------

Add the following lines to the file "~/.bashrc":
```bash
# *** OPTIMIZATION FILES *** #
export OPTDIR=$HOME/Documents/Optimization  # Change this line accordingly to the path you used

export ARCHDEFS="$OPTDIR/archdefs"
export SIFDECODE="$OPTDIR/sifdecode"
export CUTEST="$OPTDIR/cutest"
export PATH="${SIFDECODE}/bin:${PATH}"
export PATH="${CUTEST}/bin:${PATH}"
export MANPATH="${SIFDECODE}/man:${MANPATH}"
export MANPATH="${CUTEST}/man:${MANPATH}"
export MYARCH="pc.lnx.gfo"  # modify accordingly to your architecture
export MASTSIF="$OPTDIR/sif"

export PYCUTEST="$CUTEST/src/tools"
export PYCUTEST_CACHE="$OPTDIR/cache"
export PYTHONPATH="$PYCUTEST_CACHE:$CUTEST/src/python/"
```

Reference
---------
CUTEst instalation was done acordingly to the guidelines on https://github.com/optimizers/cutest-mirror
   

