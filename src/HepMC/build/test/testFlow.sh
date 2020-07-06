#! /bin/bash
# Generated by cmake

rm -f testFlow.out 
./testFlow

if [ "$?" -ne 0 ]
then
  exit 1;
fi

if [ ! -e testFlow.out ]
then
   echo "ERROR: testFlow.out does not exist"
   exit 1;
fi

OS=`uname`
case "$OS" in
CYGWIN*)
  cmd1=`sed 's/e+0/e+/g'  testFlow.out  | \
     sed 's/e-0/e-/g' | \
	diff -q -b - /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/test/testFlow.output`
  cmd2=`diff -q -b testFlow.out4 testFlow.out5`
  ;;
*)
  cmd1=`diff -q -b testFlow.out /scratch/blixen/plumberg/Pythia_and_HBT/src/HepMC/HepMC-2.06.09/test/testFlow.output`
  cmd2=`diff -q -b testFlow.out4 testFlow.out5`
esac

if [ -n "$cmd1" ]
then
  echo $cmd1
  exit 1;
fi

if [ -n "$cmd2" ]
then
  echo $cmd2
  exit 1;
fi

exit 0;
