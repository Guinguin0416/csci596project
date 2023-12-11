def normalize(v, max_velocity):
    return [abs(vel / max_velocity) for vel in v]

velocities = []
with open('lmd.xyz', 'r') as file:
    for line in file:
        if line.startswith('Ar'):
            _, _, _, _, vx, vy, vz = line.split()
            velocities.append([float(vx), float(vy), float(vz)])

max_velocity = max(max(abs(v) for v in velocity) for velocity in velocities)
colors = [normalize(v, max_velocity) for v in velocities]

with open('colors.dat', 'w') as file:
    for color in colors:
        file.write(f"{color[0]} {color[1]} {color[2]}\n")
