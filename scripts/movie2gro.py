# Released under MIT License.
# Copyright (c) 2023 Ladislav Bartos

"""
Converts Multivalent movie file into gro file readable by VMD.
"""

import sys

particles = []
first_step = True

for line in open(sys.argv[1]):
    if line.strip() == "" or line[0] == "#":
        continue

    if line[0] == "@":
        if first_step:
            for i in range(int(line.split()[8])):
                particles.append([])
            first_step = False
        
        continue

    split = line.split()
    particle_n = int(split[0])
    particles[particle_n].append((float(split[1]), float(split[2])))


with open("movie.gro", "w") as output:
    for step in range(len(particles[0])):
        output.write(f"Multivalent simulation frame t= {step:10.5f} step= {step}\n")
        output.write(f"{len(particles)}\n")
        for part in range(len(particles)):
            output.write(f"{1:5d}{'MULT':<5s}{str(part):>5s}{part:5d}{particles[part][step][0]:8.3f}{particles[part][step][1]:8.3f}{0:8.3f}\n")
        output.write("10000 10000 10000\n")