#!/usr/bin/env python3.7
import Adafruit_GPIO
import Adafruit_GPIO.I2C as I2C
import time
import sys
import argparse


argparser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="Select I2C channel multiplexed by TCA9548A")
argparser.add_argument('ch', help="channel", type=int, default=0)
args = argparser.parse_args()

TCA9548A = I2C.get_i2c_device(0x70)

TCA9548A.write8(0, 1<<args.ch)
