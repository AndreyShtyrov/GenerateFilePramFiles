import numpy as np
import math


def RX(phi: float) -> np.ndarray:
    result = np.zeros((3, 3))
    result[0, 0] = 1
    result[1, 1] = math.cos(phi)
    result[2, 2] = math.cos(phi)
    result[1, 2] = math.sin(phi)
    result[2, 1] = -math.sin(phi)
    return result


def RY(phi:  float) -> np.ndarray:
    result = np.zeros((3, 3))
    result[1, 1] = 1
    result[0, 0] = math.cos(phi)
    result[2, 2] = math.cos(phi)
    result[0, 2] = -math.sin(phi)
    result[2, 0] = math.sin(phi)
    return result


def RZ(phi: float) -> np.ndarray:
    result = np.zeros((3, 3))
    result[2, 2] = 1
    result[0, 0] = math.cos(phi)
    result[1, 1] = math.cos(phi)
    result[0, 1] = math.sin(phi)
    result[1, 0] = -math.sin(phi)
    return result