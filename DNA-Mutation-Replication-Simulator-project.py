from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
import math
import random
import time
import sys
import os
WIN_W = 1280
WIN_H = 760
NUM_PAIRS = 22
PAIR_SPACING = 0.55
HELIX_RADIUS = 1.6
TWIST_PER_PAIR = 32.0
BASE_SPHERE_RADIUS = 0.22
UNZIP_MAX = 2.8          
UNZIP_BREAK = 0.15       
BASES = ['A', 'T', 'G', 'C']
COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
RNA_COMPLEMENT = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}
BASE_COLOR = { 'A': (0.95, 0.30, 0.35), 'T': (0.25, 0.75, 0.95), 'G': (0.95, 0.85, 0.25), 'C': (0.40, 0.95, 0.45), 'U': (0.85, 0.45, 0.95)}
strand_left = []
strand_right = []
mutation_counts = []
total_mutations = 0
cam_yaw = 35.0
cam_height = 4.0
cam_distance = 16.0
cam_target_y = (NUM_PAIRS * PAIR_SPACING) * 0.5
rotation_angle = 0.0
speed = 1.0
paused = False
selected_pair = NUM_PAIRS // 2
selected_pulse = 0.0
unzip_target = 0.0
unzip_value = 0.0
replicating = False
replication_progress = 0.0
replication_done = False
heatmap_on = False
heatmap_blend = 0.0
rna_active = False
rna_progress = 0.0
rna_strand = []
rna_done = False
tour_active = False
tour_time = 0.0
TOUR_DURATION = 10.0
tour_start_yaw = 0.0
tour_start_height = 4.0
particles = []
last_time = time.time()
quadric = None
focus_mode = False
simulation_time = 0.0
events = []
mutation_impacts = []
mutation_history = [(0.0, 0)]






def init_dna():
    global strand_left, strand_right, mutation_counts, total_mutations, mutation_impacts, mutation_history, simulation_time
    strand_left = [random.choice(BASES) for _ in range(NUM_PAIRS)]
    strand_right = [COMPLEMENT[b] for b in strand_left]
    mutation_counts = [0] * NUM_PAIRS
    total_mutations = 0
    mutation_impacts = ["NONE"] * NUM_PAIRS
    mutation_history = [(0.0, 0)]
    simulation_time = 0.0





def base_position(i, side):
    base_angle = rotation_angle + i * TWIST_PER_PAIR
    if side == 'L':
        a = math.radians(base_angle)
        ux = -unzip_value
    else:
        a = math.radians(base_angle + 180.0)
        ux = unzip_value
    x = HELIX_RADIUS * math.cos(a) + ux
    y = i * PAIR_SPACING
    z = HELIX_RADIUS * math.sin(a)
    return (x, y, z)





def new_partner(i, side):
    base_angle = rotation_angle + i * TWIST_PER_PAIR
    if side == 'L':
        a = math.radians(base_angle + 180.0)
        cx = -unzip_value
    else:
        a = math.radians(base_angle)
        cx = unzip_value
    x = HELIX_RADIUS * math.cos(a) + cx
    y = i * PAIR_SPACING
    z = HELIX_RADIUS * math.sin(a)
    return (x, y, z)




def heatmap_color(count):
    if count <= 0:
        return (0.20, 0.95, 0.40)
    else:
        return (0.95, 0.20, 0.20)




def draw_line_3d(p1, p2, thickness=0.04):
    vx = p2[0] - p1[0]
    vy = p2[1] - p1[1]
    vz = p2[2] - p1[2]
    length = math.sqrt(vx*vx + vy*vy + vz*vz)
    if length < 0.0001:
        return
    ax = -vy
    ay = vx
    az = 0.0
    dot_val = vz / length
    dot_val = max(-1.0, min(1.0, dot_val))
    angle_rad = math.acos(dot_val)
    angle_deg = math.degrees(angle_rad)
    glPushMatrix()
    glTranslatef(p1[0], p1[1], p1[2])
    if abs(ax) > 0.0001 or abs(ay) > 0.0001:
        glRotatef(angle_deg, ax, ay, az)
    elif vz < 0:
        glRotatef(180, 1, 0, 0)
    gluCylinder(quadric, thickness, thickness, length, 14, 1)
    glPopMatrix()


def step(t):
    if t <= 0.0:
        return 0.0
    if t >= 1.0:
        return 1.0
    return 0.5 - 0.5 * math.cos(math.pi * t)


def draw_line_2d(x1, y1, x2, y2, thickness=1.5):
    dx = x2 - x1
    dy = y2 - y1
    length = math.sqrt(dx * dx + dy * dy)
    if length < 0.0001:
        return
    nx = -dy / length * thickness * 0.5
    ny = dx / length * thickness * 0.5
    glBegin(GL_QUADS)
    glVertex3f(x1 + nx, y1 + ny, 0.0)
    glVertex3f(x2 + nx, y2 + ny, 0.0)
    glVertex3f(x2 - nx, y2 - ny, 0.0)
    glVertex3f(x1 - nx, y1 - ny, 0.0)
    glEnd()





def draw_grid():
    glColor3f(0.14, 0.16, 0.26)
    s = 12
    y = -1.0
    for i in range(-s, s + 1):
        draw_line_3d((float(i), y, float(-s)), (float(i), y, float(s)), 0.02)
        draw_line_3d((float(-s), y, float(i)), (float(s), y, float(i)), 0.02)




def draw_strand(side):
    if side == 'L':
        glColor3f(0.75, 0.80, 0.95)
    else:
        glColor3f(0.95, 0.85, 0.70)
    for i in range(NUM_PAIRS - 1):
        p1 = base_position(i, side)
        p2 = base_position(i + 1, side)
        draw_line_3d(p1, p2, 0.06)






def pair_rungs():
    if unzip_value > UNZIP_BREAK:
        return
    fade = 1.0 - (unzip_value / UNZIP_BREAK)
    for i in range(NUM_PAIRS):
        pl = base_position(i, 'L')
        pr = base_position(i, 'R')
        if heatmap_on:
            c = heatmap_color(mutation_counts[i])
        else:
            c = (0.55, 0.55, 0.70)
        if i == selected_pair:
            pulse = 0.5 + 0.5 * math.sin(selected_pulse)
            c = (1.0, 0.4 + 0.6 * pulse, 0.3)
        glColor3f(c[0] * fade, c[1] * fade, c[2] * fade)
        draw_line_3d(pl, pr, 0.04)




def all_bases():
    for i in range(NUM_PAIRS):
        for side, letter in [('L', strand_left[i]), ('R', strand_right[i])]:
            pos = base_position(i, side)

            base_col = BASE_COLOR.get(letter, (1.0, 1.0, 1.0))
            if heatmap_blend > 0.001:
                h_col = heatmap_color(mutation_counts[i])
                col = [base_col[0] * (1.0 - heatmap_blend) + h_col[0] * heatmap_blend,base_col[1] * (1.0 - heatmap_blend) + h_col[1] * heatmap_blend,base_col[2] * (1.0 - heatmap_blend) + h_col[2] * heatmap_blend ]
            else:
                col = [base_col[0], base_col[1], base_col[2]]

            scale_mul = 1.0
            if i == selected_pair:
                pulse = 0.5 + 0.5 * math.sin(selected_pulse)
                scale_mul = 1.0 + (0.60 if focus_mode else 0.30) * pulse
                col = [min(1.0, col[0] + 0.35 * pulse),min(1.0, col[1] + 0.35 * pulse), min(1.0, col[2] + 0.35 * pulse) ]
            elif focus_mode:
                col = [col[0] * 0.3, col[1] * 0.3, col[2] * 0.3]

            glColor3f(col[0], col[1], col[2])
            glPushMatrix()
            glTranslatef(pos[0], pos[1], pos[2])
            rad = BASE_SPHERE_RADIUS * scale_mul
            glScalef(rad, rad, rad)
            gluSphere(quadric, 1.0, 24, 20)
            glPopMatrix()

            if paused and i == selected_pair:
                glPushMatrix()
                glTranslatef(pos[0], pos[1], pos[2])
                glColor3f(0.2, 1.0, 1.0)
                ring_rad = 0.45
                segments = 16
                for s in range(segments):
                    a1 = (s / segments) * math.pi * 2.0
                    a2 = ((s + 1.0) / segments) * math.pi * 2.0
                    p1 = (math.cos(a1)*ring_rad, math.sin(a1)*ring_rad, 0.0)
                    p2 = (math.cos(a2)*ring_rad, math.sin(a2)*ring_rad, 0.0)
                    draw_line_3d(p1, p2, 0.02)
                glPopMatrix()





def draw_replication():
    if not (replicating or replication_done):
        return
    full = int(replication_progress)
    frac = replication_progress - full
    for i in range(NUM_PAIRS):
        if i < full:
            growth = 1.0
        elif i == full and frac > 0.0:
            growth = step(frac)
        else:
            continue
        np_l = new_partner(i, 'L')
        new_r_letter = COMPLEMENT[strand_left[i]]
        c = BASE_COLOR[new_r_letter]
        glColor3f(c[0] * 0.95, c[1] * 0.95, c[2] * 0.95)
        glPushMatrix()
        glTranslatef(np_l[0], np_l[1], np_l[2])
        s = BASE_SPHERE_RADIUS * growth
        glScalef(s, s, s)
        gluSphere(quadric, 1.0, 20, 16)
        glPopMatrix()
        lp = base_position(i, 'L')
        glColor3f(0.45, 0.85, 0.55)
        draw_line_3d(lp, np_l, 0.04 * growth)
        np_r = new_partner(i, 'R')
        new_l_letter = COMPLEMENT[strand_right[i]]
        c = BASE_COLOR[new_l_letter]
        glColor3f(c[0] * 0.95, c[1] * 0.95, c[2] * 0.95)
        glPushMatrix()
        glTranslatef(np_r[0], np_r[1], np_r[2])
        glScalef(s, s, s)
        gluSphere(quadric, 1.0, 20, 16)
        glPopMatrix()
        rp = base_position(i, 'R')
        glColor3f(0.45, 0.85, 0.55)
        draw_line_3d(rp, np_r, 0.04 * growth)
    if full >= 2 or (full >= 1 and frac > 0.0):
        glColor3f(0.70, 0.95, 0.75)
        last = min(full, NUM_PAIRS) - 1
        for i in range(last):
            p1 = new_partner(i, 'L')
            p2 = new_partner(i + 1, 'L')
            draw_line_3d(p1, p2, 0.06)

            p1 = new_partner(i, 'R')
            p2 = new_partner(i + 1, 'R')
            draw_line_3d(p1, p2, 0.06)







RNA_X_OFFSET = 5.0

def draw_rna():
    if not rna_active and not rna_strand:
        return
    glColor3f(0.80, 0.45, 0.95)
    for i in range(len(rna_strand) - 1):
        y1 = i * PAIR_SPACING
        y2 = (i + 1) * PAIR_SPACING
        draw_line_3d((RNA_X_OFFSET, y1, 0.0), (RNA_X_OFFSET, y2, 0.0), 0.06)
    if rna_active and len(rna_strand) > 0:
        idx = len(rna_strand)
        if idx < NUM_PAIRS:
            y1 = (idx - 1) * PAIR_SPACING
            frac_ext = step(rna_progress - int(rna_progress))
            y2 = y1 + PAIR_SPACING * frac_ext
            draw_line_3d((RNA_X_OFFSET, y1, 0.0), (RNA_X_OFFSET, y2, 0.0), 0.06)

    for i, letter in enumerate(rna_strand):
        c = BASE_COLOR.get(letter, (1.0, 1.0, 1.0))
        glColor3f(c[0], c[1], c[2])
        glPushMatrix()
        glTranslatef(RNA_X_OFFSET, i * PAIR_SPACING, 0.0)
        glScalef(0.20, 0.20, 0.20)
        gluSphere(quadric, 1.0, 18, 14)
        glPopMatrix()
    if rna_active:
        idx = min(int(rna_progress), NUM_PAIRS - 1)
        lp = base_position(idx, 'L')
        frac = rna_progress - int(rna_progress)
        y_pos = idx * PAIR_SPACING
        
        pulse = 0.5 + 0.5 * math.sin(time.time() * 8.0)
        glColor3f(1.0, 1.0, 0.35 + 0.65 * pulse)
        draw_line_3d(lp, (RNA_X_OFFSET, y_pos, 0.0), 0.05)
        
        glPushMatrix()
        glTranslatef(RNA_X_OFFSET, y_pos, 0.0)
        s = 0.25 + 0.1 * pulse
        glScalef(s, s, s)
        gluSphere(quadric, 1.0, 12, 12)
        glPopMatrix()






def spawn_particles(pos, color, count=40):
    for _ in range(count):
        vx = random.uniform(-1.0, 1.0)
        vy = random.uniform(-1.0, 1.0)
        vz = random.uniform(-1.0, 1.0)
        m = math.sqrt(vx * vx + vy * vy + vz * vz) + 1e-5
        sp = random.uniform(1.2, 2.6)
        vx, vy, vz = vx / m * sp, vy / m * sp, vz / m * sp
        cr = min(1.0, max(0.0, color[0] + random.uniform(-0.2, 0.2)))
        cg = min(1.0, max(0.0, color[1] + random.uniform(-0.2, 0.2)))
        cb = min(1.0, max(0.0, color[2] + random.uniform(-0.2, 0.2)))
        particles.append({'pos': [pos[0], pos[1], pos[2]], 'vel': [vx, vy, vz], 'life': 1.0, 'color': (cr, cg, cb), 'size': random.uniform(0.04, 0.1) })






def update_particles(dt):
    for p in particles[:]:
        p['pos'][0] += p['vel'][0] * dt
        p['pos'][1] += p['vel'][1] * dt
        p['pos'][2] += p['vel'][2] * dt
        p['vel'][1] -= 3.5 * dt
        p['life'] -= 0.9 * dt
        p['size'] = max(0.0, p['size'] - 0.05 * dt)
        if p['life'] <= 0.0 or p['size'] <= 0.0:
            particles.remove(p)




def d_particles():
    for p in particles:
        a = max(0.0, p['life'])
        r, g, b = p['color']
        glColor3f(r * a, g * a, b * a)
        glPushMatrix()
        glTranslatef(p['pos'][0], p['pos'][1], p['pos'][2])
        sc = p['size'] * (0.5 + 0.5 * a)
        glScalef(sc, sc, sc)
        gluSphere(quadric, 1.0, 8, 8)
        glPopMatrix()




def draw_text(x, y, text, font=None):
    if font is None:
        font = GLUT_BITMAP_HELVETICA_18
    glRasterPos2f(x, y)
    for ch in text:
        glutBitmapCharacter(font, ord(ch))




def hud():
    glMatrixMode(GL_PROJECTION)
    glPushMatrix()
    glLoadIdentity()
    gluOrtho2D(0.0, float(WIN_W), 0.0, float(WIN_H))
    glMatrixMode(GL_MODELVIEW)
    glPushMatrix()
    glLoadIdentity()
    glColor3f(0.45, 0.92, 1.0)
    draw_text(20, WIN_H - 32, "DNA  MUTATION  &  REPLICATION  SIMULATOR",GLUT_BITMAP_HELVETICA_18)
    sb = strand_left[selected_pair]
    sbc = strand_right[selected_pair]
    glColor3f(1.0, 1.0, 1.0)
    draw_text(20, WIN_H - 60, "MODE : {}".format(label()))
    draw_text(20, WIN_H - 82, "Selected Pair  : {}   ({} - {})   Mutations: {}   Impact: {}".format(
                  selected_pair + 1, sb, sbc, mutation_counts[selected_pair], mutation_impacts[selected_pair]))
    draw_text(20, WIN_H - 102, "Total Mutations: {}".format(total_mutations))
    draw_text(20, WIN_H - 122, "Speed          : {:.2f} x".format(speed))

    glColor3f(0.6 if not heatmap_on else 1.0, 1.0 if heatmap_on else 0.6, 0.5)
    draw_text(20, WIN_H - 142, "Heatmap        : {}".format("ON" if heatmap_on else "OFF"))
    glColor3f(1.0, 0.6, 0.4) if paused else glColor3f(0.6, 1.0, 0.6)
    draw_text(20, WIN_H - 162, "State          : {}".format("PAUSED" if paused else "RUNNING"))

    glColor3f(1.0, 1.0, 1.0)
    if rna_active:
        rna_status = "Transcribing...  {}/{}".format(min(int(rna_progress), NUM_PAIRS), NUM_PAIRS)
    elif rna_done:
        rna_status = "RNA Complete  ({} bases)".format(len(rna_strand))
    else:
        rna_status = "Idle"
    draw_text(20, WIN_H - 182, "RNA            : " + rna_status)

    if replicating:
        rep_status = "Replicating  {}/{}".format(int(replication_progress), NUM_PAIRS)
    elif replication_done:
        rep_status = "Replicated"
    else:
        rep_status = "Idle"
    draw_text(20, WIN_H - 202, "Replication    : " + rep_status)

    pct = int(unzip_value / UNZIP_MAX * 100.0)
    draw_text(20, WIN_H - 222, "Unzip          : {}%".format(pct))

    if tour_active:
        tour_status = "ACTIVE  ({:.1f}s left)".format(max(0.0, TOUR_DURATION - tour_time))
        glColor3f(1.0, 0.85, 0.4)
    else:
        tour_status = "OFF"
        glColor3f(1.0, 1.0, 1.0)
    draw_text(20, WIN_H - 242, "Auto Tour      : " + tour_status)

    graph_x = 20
    graph_y = 40
    graph_w = 240
    graph_h = 80

    glColor3f(0.55, 0.85, 0.55)
    draw_text(graph_x, graph_y + graph_h + 14, "MUTATION TIMELINE")

    glColor3f(0.30, 0.30, 0.40)
    draw_line_2d(graph_x, graph_y, graph_x + graph_w, graph_y, 1.6)
    draw_line_2d(graph_x, graph_y, graph_x, graph_y + graph_h, 1.6)
    draw_line_2d(graph_x, graph_y + graph_h, graph_x + graph_w, graph_y + graph_h, 1.6)
    draw_line_2d(graph_x + graph_w, graph_y, graph_x + graph_w, graph_y + graph_h, 1.6)

    glColor3f(0.18, 0.20, 0.28)
    for gi in range(1, 4):
        gy = graph_y + (graph_h * gi / 4.0)
        draw_line_2d(graph_x + 1, gy, graph_x + graph_w - 1, gy, 1.0)

    if len(mutation_history) > 0:
        max_time = max(10.0, simulation_time)
        max_muts = max(5, total_mutations)
        points = []
        for t, m in mutation_history:
            px = graph_x + (t / max_time) * graph_w
            py = graph_y + (m / max_muts) * graph_h
            points.append((px, py))
        points.append((graph_x + (simulation_time / max_time) * graph_w, graph_y + (total_mutations / max_muts) * graph_h ))

        glColor3f(1.0, 0.45, 0.45)
        for i in range(len(points) - 1):
            draw_line_2d(points[i][0], points[i][1],
                         points[i + 1][0], points[i + 1][1], 2.0)

        glColor3f(1.0, 0.85, 0.30)
        for px, py in points:
            glPushMatrix()
            glTranslatef(px, py, 0.0)
            glScalef(2.5, 2.5, 1.0)
            gluSphere(quadric, 1.0, 6, 6)
            glPopMatrix()

    glColor3f(0.65, 0.65, 0.75)
    draw_text(graph_x + graph_w + 8, graph_y + graph_h - 8,"{} muts".format(total_mutations))
    draw_text(graph_x + graph_w + 8, graph_y - 4,"{:.0f}s".format(simulation_time))

    legend_x = WIN_W - 220
    legend_y = WIN_H - 32
    glColor3f(0.7, 0.7, 0.9)
    draw_text(legend_x, legend_y, "Bases")
    glColor3f(*BASE_COLOR['A'])
    draw_text(legend_x, legend_y - 22, "A  Adenine")
    glColor3f(*BASE_COLOR['T'])
    draw_text(legend_x, legend_y - 44, "T  Thymine")
    glColor3f(*BASE_COLOR['G'])
    draw_text(legend_x, legend_y - 66, "G  Guanine")
    glColor3f(*BASE_COLOR['C'])
    draw_text(legend_x, legend_y - 88, "C  Cytosine")
    glColor3f(*BASE_COLOR['U'])
    draw_text(legend_x, legend_y - 110, "U  Uracil (RNA)")

    glPopMatrix()
    glMatrixMode(GL_PROJECTION)
    glPopMatrix()
    glMatrixMode(GL_MODELVIEW)



def label():
    if tour_active: 
        return "AUTO CAMERA TOUR"
    if rna_active: 
        return "RNA TRANSCRIPTION"
    if replicating: 
        return "REPLICATION (forming 2 helices)"
    if replication_done:
        return "TWO DAUGHTER HELICES"
    if unzip_value > UNZIP_BREAK:
        return "UNZIPPED (severed)"
    if heatmap_on: 
        return "HEATMAP VIEW"
    if paused:
        return "ANALYZER (PAUSED)"
    return "ROTATING HELIX"





def display():
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(60.0, float(WIN_W) / float(WIN_H), 0.1, 200.0)

    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()

    cam_x = cam_distance * math.cos(math.radians(cam_yaw))
    cam_z = cam_distance * math.sin(math.radians(cam_yaw))
    gluLookAt(cam_x, cam_height + cam_target_y, cam_z,
              0.0, cam_target_y, 0.0,
              0.0, 1.0, 0.0)

    draw_grid()
    draw_strand('L')
    draw_strand('R')
    pair_rungs()
    all_bases()
    draw_replication()
    draw_rna()
    d_particles()

    hud()

    glutSwapBuffers()




def update():
    global rotation_angle, selected_pulse, unzip_value, replicating, replication_progress, replication_done, rna_active, rna_progress, rna_done, tour_active, tour_time, cam_yaw, cam_height, cam_distance, last_time,heatmap_blend, simulation_time, cam_target_y
    now = time.time()
    dt = now - last_time
    last_time = now
    if dt > 0.1:
        dt = 0.1
        
    simulation_time += dt
        
    for e in events:
        e['time_left'] -= dt
        
    target_blend = 1.0 if heatmap_on else 0.0
    heatmap_blend += (target_blend - heatmap_blend) * min(1.0, 3.5 * dt)

    speed_mult = speed
    cam_lerp = min(1.0, 2.2 * dt)
    if focus_mode:
        speed_mult = speed * 0.2
        target_dist = 6.0
        target_y = selected_pair * PAIR_SPACING
        cam_distance += (target_dist - cam_distance) * cam_lerp
        cam_target_y += (target_y - cam_target_y) * cam_lerp
        cam_height += (3.0 - cam_height) * cam_lerp
    else:
        default_target = (NUM_PAIRS * PAIR_SPACING) * 0.5
        cam_target_y += (default_target - cam_target_y) * min(1.0, 1.5 * dt)

    if not paused:
        rotation_angle += 28.0 * speed_mult * dt
    selected_pulse += 4.5 * dt
    diff = unzip_target - unzip_value
    if abs(diff) > 0.001:
        unzip_value += diff * min(1.0, 1.4 * speed_mult * dt)
    else:
        unzip_value = unzip_target
    if replicating and not paused:
        replication_progress += 1.0 * speed_mult * dt
        if replication_progress >= NUM_PAIRS:
            replication_progress = float(NUM_PAIRS)
            replicating = False
            replication_done = True

    if rna_active and not paused:
        rna_progress += 1.0 * speed_mult * dt
        target = min(int(rna_progress), NUM_PAIRS)
        while len(rna_strand) < target:
            i = len(rna_strand)
            rna_strand.append(RNA_COMPLEMENT[strand_left[i]])
        if rna_progress >= NUM_PAIRS:
            rna_progress = float(NUM_PAIRS)
            rna_active = False
            rna_done = True
            
            while len(rna_strand) < NUM_PAIRS:
                i = len(rna_strand)
                rna_strand.append(RNA_COMPLEMENT[strand_left[i]])
    if tour_active:
        tour_time += dt
        t = min(1.0, tour_time / TOUR_DURATION)
        ease = 0.5 - 0.5 * math.cos(math.pi * t)
        cam_yaw = tour_start_yaw + 360.0 * ease
        cam_height = tour_start_height + math.sin(math.pi * t) * 12.0
        if tour_time >= TOUR_DURATION:
            tour_active = False
            tour_time = 0.0
            cam_yaw = tour_start_yaw + 360.0
            cam_height = tour_start_height
    update_particles(dt)
    glutPostRedisplay()





LOW_TRANSITIONS = {('A', 'G'), ('T', 'C')}
MEDIUM_TRANSITIONS = {('A', 'C'), ('G', 'T')}


def impact(old_base, new_base):
    pair = (old_base, new_base)
    if pair in LOW_TRANSITIONS:
        return "LOW"
    if pair in MEDIUM_TRANSITIONS:
        return "MEDIUM"
    return "HIGH"


def do_mutation():
    global total_mutations
    old = strand_left[selected_pair]
    choices = [b for b in BASES if b != old]
    new_base = random.choice(choices)
    strand_left[selected_pair] = new_base
    strand_right[selected_pair] = COMPLEMENT[new_base]
    mutation_counts[selected_pair] += 1
    total_mutations += 1

    mutation_impacts[selected_pair] = impact(old, new_base)
        
    mutation_history.append((simulation_time, total_mutations))
    
    pl = base_position(selected_pair, 'L')
    pr = base_position(selected_pair, 'R')
    mid = ((pl[0] + pr[0]) * 0.5, (pl[1] + pr[1]) * 0.5, (pl[2] + pr[2]) * 0.5)
    spawn_particles(mid, BASE_COLOR[new_base], 45)






def unzip():
    global unzip_target, replicating, replication_progress, replication_done
    if unzip_target < UNZIP_MAX * 0.5:
        unzip_target = UNZIP_MAX
    else:
        unzip_target = 0.0
        if replicating:
            replicating = False
            replication_progress = 0.0
            replication_done = False





def keyboard(key, x, y):
    global paused, heatmap_on, replicating, replication_progress, replication_done, unzip_target, rna_active, rna_progress, rna_strand, rna_done, tour_active, tour_time, speed, cam_distance, selected_pair, focus_mode, tour_start_yaw, tour_start_height
    if key == b'\x1b': 
        os._exit(0)
    elif key == b'm' or key == b'M': 
        do_mutation()
    elif key == b'u' or key == b'U': 
        unzip()
    elif key == b'r' or key == b'R':
        if replication_done or replicating:
            replicating = False
            replication_done = False
            replication_progress = 0.0
            unzip_target = 0.0
        elif unzip_value > UNZIP_MAX * 0.6:
            replicating = True
            replication_progress = 0.0
            replication_done = False
    elif key == b'h' or key == b'H': 
        heatmap_on = not heatmap_on
    elif key == b'p' or key == b'P':
        paused = not paused
    elif key == b'f' or key == b'F':
        focus_mode = not focus_mode
    elif key == b't' or key == b'T':
        if not rna_active:
            rna_active = True
            rna_done = False
            rna_progress = 0.0
            rna_strand = []
    elif key == b'c' or key == b'C':
        tour_active = not tour_active
        tour_time = 0.0
        if tour_active:
            tour_start_yaw = cam_yaw
            tour_start_height = cam_height
    elif key == b'w' or key == b'W': 
        cam_distance = max(5.0, cam_distance - 0.7)
    elif key == b's' or key == b'S':
        cam_distance = min(40.0, cam_distance + 0.7)
    elif key == b'[': 
        selected_pair = (selected_pair - 1) % NUM_PAIRS
    elif key == b']': 
        selected_pair = (selected_pair + 1) % NUM_PAIRS
    elif key == b'+' or key == b'=': 
        speed = min(10.0, speed + 0.2)
    elif key == b'-' or key == b'_': 
        speed = max(0.1, speed - 0.2)

    glutPostRedisplay()





def special(key, x, y):
    global cam_yaw, cam_height
    if key == GLUT_KEY_LEFT:
        cam_yaw += 4.5
    elif key == GLUT_KEY_RIGHT:
        cam_yaw -= 4.5
    elif key == GLUT_KEY_UP:
        cam_height = min(20.0, cam_height + 0.5)
    elif key == GLUT_KEY_DOWN:
        cam_height = max(-10.0, cam_height - 0.5)
    glutPostRedisplay()





def mouse(button, state, x, y):
    if state != GLUT_DOWN:
        return
    if button == GLUT_LEFT_BUTTON:
        do_mutation()
    elif button == GLUT_RIGHT_BUTTON:
        unzip()
    glutPostRedisplay()





def main():
    global quadric, last_time

    init_dna()

    glutInit()
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
    glutInitWindowSize(WIN_W, WIN_H)
    glutInitWindowPosition(60, 40)
    glutCreateWindow(b"DNA Mutation & Replication Simulator")

    quadric = gluNewQuadric()

    glutDisplayFunc(display)
    glutKeyboardFunc(keyboard)
    glutSpecialFunc(special)
    glutMouseFunc(mouse)
    glutIdleFunc(update)

    last_time = time.time()
    glutMainLoop()

if __name__ == '__main__':
    main()
