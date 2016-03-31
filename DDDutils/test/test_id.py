#!/usr/bin/env python

import DDDutils.readin.utils as read


p = '/raid/markus/simulations/multimech/frsources/AR3_4/100/config1'

print p

id_init, id_size, id_dir = read.get_ids_from_simulations(p)


print 'id_init, id_size, id_dir = ', id_init, id_size, id_dir
