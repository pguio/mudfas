#/bin/bash

#
# $Id: bench-solver.sh,v 1.2 2011/03/26 12:56:40 patrick Exp $
#
# Copyright (c) 2000-2011 Patrick Guio <patrick.guio@gmail.com>
# All Rights Reserved.
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


if [ $# -lt 1 ] ; then
	echo $"Usage: $0 {2|3}d-{gridsize|bdyc|kcycle}"
	exit 0
fi

case "$1" in

	2d-gridsize)

	xgr=(103 101 105 113 97)
	ygr=(103 101 105 113 97)

	echo "******************************"
	echo "Testing 2d grid size algorithm"
	echo "******************************"
	for prog in PoissonBoltzmann2d Poisson2d ; do
		if [ -x $prog ] ; then
	  	for i in 0 1 2 3 4 ; do
		  	for j in 0 1 2 3 4 ; do
			  	cmd="./$prog gridsize=${xgr[$i]},${ygr[$j]}"
					echo $cmd 
					$cmd >> $1.log
				done
			done
		fi
	done
	;;

	3d-gridsize)

	xgr=(39 37 25 33)
	ygr=(39 37 25 33)
	zgr=(39 37 25 33)

	echo "******************************"
	echo "Testing 3d grid size algorithm"
	echo "******************************"
	for prog in PoissonBoltzmann3d Poisson3d ; do
		if [ -x $prog ] ; then
	  	for i in 0 1 2 3 ; do
		  	for j in 0 1 2 3 ; do
			  	cmd="./$prog gridsize=${xgr[$i]},${ygr[$j]},${zgr[$j]}"
					echo $cmd 
					$cmd >> $1.log
				done
			done
		fi
	done
	;;

	2d-bdyc)

	bxmn=(0 1 1 2 2)
	bxmx=(0 1 2 1 2)
	bymn=(0 1 1 2 2)
	bymx=(0 1 2 1 2)

	echo "******************************"
	echo "Testing 2d boundary conditions"
	echo "******************************"
	for prog in PoissonBoltzmann2d Poisson2d ; do
		if [ -x $prog ] ; then
	  	for i in 0 1 2 3 4 ; do
		  	for j in 0 1 2 3 4 ; do
			  	cmd="./$prog bcmin=${bxmn[$i]},${bymn[$j]} \
				  bcmax=${bxmx[$i]},${bymx[$j]}"
					echo $cmd 
					$cmd >> $1.log
				done
			done
		fi
	done
	;;

	3d-bdyc)

	bxmn=(0 1 1 2 2)
	bxmx=(0 1 2 1 2)
	bymn=(0 1 1 2 2)
	bymx=(0 1 2 1 2)
	bzmn=(0 1 1 2 2)
	bzmx=(0 1 2 1 2)

	echo "******************************"
	echo "Testing 3d boundary conditions"
	echo "******************************"
	for prog in PoissonBoltzmann3d Poisson3d ; do
		if [ -x $prog ] ; then
			for i in 0 1 2 3 4 ; do
				for j in 0 1 2 3 4 ; do
					for k in 0 1 2 3 4 ; do
						cmd="./$prog bcmin=${bxmn[$i]},${bymn[$j]},${bzmn[$k]} \
						 bcmax=${bxmx[$i]},${bymx[$j]},${bzmx[$k]}"
						echo $cmd
						$cmd >> $1.log
			  	done
				done
			done
		fi
	done
	;;

	kcycle)

	echo "******************"
	echo "Testing cycle type"
	echo "******************"

	for prog in PoissonBoltzmann2d Poisson2d PoissonBoltzmann3d Poisson3d ; do
		if [ -x $prog ] ; then
	  	for fun in 1 2 3 4 ; do
		  	for kcycle in 1 2 ; do
					for intpol in 1 3 ; do
	  				cmd="./$prog fun $fun kcycle $kcycle intpol $intpol"
						echo $cmd
						$cmd >> $1.log
					done
				done
			done
		fi
	done
	;;

esac


