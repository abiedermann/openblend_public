#!/usr/bin/env python

import sys
sys.path.append("/data/user_storage")
import opentrons.execute
from opentrons import protocol_api
from plate_builder import plate_builder

def run(protocol: protocol_api.ProtocolContext()):
    foo = plate_builder(protocol,sys.argv[1])
    foo.make_plates()

protocol = opentrons.execute.get_protocol_api('2.0')
run(protocol)

