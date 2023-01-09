import numpy as np

class Corpus():

    def __init__(self, init_pos: np.array, init_vel: np.array, mass: float):
        self.pos = init_pos.copy()
        self.velocity = init_vel.copy()
        self.abs_vel = np.linalg.norm(self.velocity)
        self.mass = mass
        self.cntp_vel = 0
        self.dist = 0
        self.udist = 0
        self.vvel_cntp = 0
