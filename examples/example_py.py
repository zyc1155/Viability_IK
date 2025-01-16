import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Button
import numpy as np

import pyviability_ik

dt = 0.005
ratio_simulation_control = 10  # simulation_interval / control_interval
simulation_time = 0

# A = np.array([[-1.0, 0.0], [1.0, 0.0], [0.0, -1.0], [0.0, 1.0]])
# t_q_max = np.array([0.0, np.pi / 2, -0.0, np.pi * 2 / 3])
A = np.array(
    [
        [1.38637, 1],
        [-2.89855, -1.21212],
        [-1, -0],
        [0.0, 1.0],
        [0.0, -1.0],
        [2.98413, 1],
    ]
)
t_q_max = np.array([3.73353, -1, 0, np.pi / 2, -0.01, 6.95302])

v_lim = np.array([3.0, 4.0])
a_lim = np.array([10, 15])

p_ref = np.array([-0.5, 1.5])
q_ref = np.zeros(2)
q0 = np.array([0, np.pi / 2])
q1 = np.zeros(2)
dq0 = np.zeros(2)
dq1 = np.zeros(2)
ddq1 = np.zeros(2)
position1, position2 = np.zeros(2)

l1, l2 = 1.0, 1.0  # Lengths of the robot arm segments
base_x = [-1.0, -1.0, 1.5, 1.5]
base_y = [1, 0, 0, 0.5]

fig = plt.figure(figsize=(16, 8))
gs = fig.add_gridspec(2, 2, height_ratios=[2, 1])

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])

ax_button = plt.axes([0.45, 0.01, 0.1, 0.05])  # Adjust the position and size as needed
btn = Button(ax_button, "Pause/Resume")
is_paused = False

line = []
robot_plot = []
line2 = []

output_data = [[], [], [], [], []]

viability_ik = pyviability_ik.VIABILITY_IK(A, t_q_max, v_lim, a_lim, dt)


def obtain_poly(pt):
    point_x = []
    point_y = []
    for p in pt:
        point_x.append(p[0])
        point_y.append(p[1])

    point_x.append(pt[0][0])
    point_y.append(pt[0][1])
    return point_x, point_y


def divide_polygon_edges(pt, n):
    """
    Divide each edge of the polygon into n equal parts.

    Parameters:
    pt (list of list): List containing vertices of the polygon [x, y].
    n (int): Number of equal parts to divide each edge into.

    Returns:
    list of list: List of points dividing each edge into n equal parts.
    """
    all_points = []
    num_vertices = len(pt)

    for i in range(num_vertices):
        # Get the start and end point of the current edge
        x1, y1 = pt[i]
        x2, y2 = pt[(i + 1) % num_vertices]  # Wrap around to the first vertex

        # Generate n+1 points from x1 to x2 and y1 to y2
        x_points = np.linspace(x1, x2, n + 1)
        y_points = np.linspace(y1, y2, n + 1)

        # Combine x and y coordinates into points
        edge_points = np.vstack((x_points, y_points)).T

        # Exclude the last point to avoid duplication with the next edge's first point
        if i < num_vertices - 1:
            edge_points = edge_points[:-1]

        all_points.extend(edge_points)

    # Convert to list of lists and return
    return [list(map(float, point)) for point in all_points]


pt = [
    [0.0, np.pi / 2],
    [0.0, 0.825],
    [0.345, 0.0],
    [2.33, 0.0],
    [2.015, 0.94],
    [1.56, np.pi / 2],
]
point_x, point_y = obtain_poly(pt)
divided_joint_points = divide_polygon_edges(pt, 10)
divided_task_points = divided_joint_points


def jacobian_matrix(q0):
    global l1, l2
    J11 = -l1 * np.sin(q0[0]) - l2 * np.sin(q0[0] + q0[1])
    J12 = -l2 * np.sin(q0[0] + q0[1])
    J21 = l1 * np.cos(q0[0]) + l2 * np.cos(q0[0] + q0[1])
    J22 = l2 * np.cos(q0[0] + q0[1])
    J = np.array([[J11, J12], [J21, J22]])
    return J


def differential_inverse_kinematics(p_ref, q0, dq0, p0):
    dp_d = 5.0 * (p_ref - p0)
    J = jacobian_matrix(q0)
    success = viability_ik.solve(J, dp_d, q0, dq0)
    return success, viability_ik.result()


def forward_kinematics(q1):
    position1 = np.array([l1 * np.cos(q1[0]), l1 * np.sin(q1[0])])
    position2 = np.array(
        [
            position1[0] + l2 * np.cos(q1[0] + q1[1]),
            position1[1] + l2 * np.sin(q1[0] + q1[1]),
        ]
    )
    return position1, position2


def inverse_kinematics(r):
    r_square_sum = np.sum(r**2)
    if r_square_sum < (l1 - l2) ** 2:
        r_square_sum = (l1 - l2) ** 2
        r = r / np.linalg.norm(r) * np.fabs(l1 - l2)
    elif r_square_sum > (l1 + l2) ** 2:
        r_square_sum = (l1 + l2) ** 2
        r = r / np.linalg.norm(r) * np.fabs(l1 + l2)
    alpha = np.sqrt(
        (r_square_sum + l1**2 + l2**2) ** 2 - 2 * ((r_square_sum) ** 2 + l1**4 + l2**4)
    )
    theta1 = np.arctan2(r[1], r[0]) - np.arctan2(alpha, r_square_sum + l1**2 - l2**2)
    theta2 = np.arctan2(alpha, r_square_sum - l1**2 - l2**2)
    return np.array([theta1, theta2])


def draw_robot(position1, position2):
    robot_plot[0].set_data(base_x, base_y)
    robot_plot[1].set_data(
        [0, position1[0], position2[0]], [0, position1[1], position2[1]]
    )
    robot_plot[2].set_offsets([np.zeros(2), position1])


def init_plot():
    global q_ref, divided_task_points, divided_plot_point_x, divided_plot_point_y
    q_ref = inverse_kinematics(p_ref)
    index = 0
    for point in divided_joint_points:
        position1, position2 = forward_kinematics(point)
        divided_task_points[index][0], divided_task_points[index][1] = (
            position2[0],
            position2[1],
        )
        index = index + 1
    divided_plot_point_x, divided_plot_point_y = obtain_poly(divided_task_points)

    ax1.set_xlim(-2, 2)
    ax1.set_ylim(-1.5, 2.5)
    ax2.set_xlim(-0.5, np.pi)
    ax2.set_ylim(-0.5, np.pi)

    ax1.set_aspect("equal")
    ax2.set_aspect("equal")

    line.append(ax1.scatter([], [], color="red"))
    line.append(ax2.scatter([], [], color="red"))
    line.append(ax2.scatter([], [], color="green"))
    line.append(ax2.plot([], [], "black")[0])
    line.append(
        ax1.text(
            1.0,
            1.0,
            "",
            horizontalalignment="right",
            verticalalignment="top",
            transform=ax1.transAxes,
        )
    )
    line.append(ax1.plot([], [], "g--")[0])

    line2.append(ax3.plot([], [], "blue")[0])
    line2.append(ax3.plot([], [], "green")[0])
    line2.append(ax4.plot([], [], "blue")[0])
    line2.append(ax4.plot([], [], "green")[0])

    robot_plot.append(ax1.plot([], [], "black")[0])
    robot_plot.append(ax1.plot([], [], "purple")[0])
    robot_plot.append(
        ax1.scatter([], [], facecolors="none", edgecolors="green")
    )  # facecolors='none', edgecolors='blue' // color="green"
    return (
        line[0],
        line[1],
        line[2],
        line[3],
        line[4],
        line[5],
        line2[0],
        line2[1],
        line2[2],
        line2[3],
        robot_plot[0],
        robot_plot[1],
        robot_plot[2],
    )


def update(frame):
    global q0, dq0, q1, dq1, ddq1, q_ref, p_ref, position1, position2, simulation_time
    output_data[0].append(simulation_time)
    output_data[1].append(dq0[0])
    output_data[2].append(dq0[1])
    output_data[3].append(ddq1[0])
    output_data[4].append(ddq1[1])

    success = True
    for i in range(ratio_simulation_control):
        success, dq1 = differential_inverse_kinematics(p_ref, q0, dq0, position2)
        ddq1 = (dq1 - dq0) / dt
        q1 = q0 + dt * (dq0 + dq1) / 2
        position1, position2 = forward_kinematics(q1)
        q0, dq0 = q1, dq1
        simulation_time += dt
        if not success:
            print("Failed at", simulation_time, "!")

    line[0].set_offsets([p_ref])
    line[1].set_offsets([q_ref])
    line[2].set_offsets([q1])
    line[3].set_data(point_x, point_y)
    line[4].set_text(f"Time: {simulation_time:.1f} s")
    line[5].set_data(divided_plot_point_x, divided_plot_point_y)

    line2[0].set_data(output_data[0], output_data[1])
    line2[1].set_data(output_data[0], output_data[2])
    line2[2].set_data(output_data[0], output_data[3])
    line2[3].set_data(output_data[0], output_data[4])

    ax3.relim()
    ax3.autoscale()
    ax4.relim()
    ax4.autoscale()

    draw_robot(position1, position2)
    return (
        line[0],
        line[1],
        line[2],
        line[3],
        line[4],
        line[5],
        line2[0],
        line2[1],
        line2[2],
        line2[3],
        robot_plot[0],
        robot_plot[1],
        robot_plot[2],
    )


def on_click(event):
    global p_ref, q_ref, is_paused
    if event.inaxes == ax1 and not is_paused:
        p_ref[0], p_ref[1] = event.xdata, event.ydata
        q_ref = inverse_kinematics(p_ref)


def pause(event):
    global is_paused
    is_paused = not is_paused
    if is_paused:
        ani.event_source.stop()
    else:
        ani.event_source.start()


fig.canvas.mpl_connect("button_press_event", on_click)
btn.on_clicked(pause)

ani = FuncAnimation(
    fig,
    update,
    init_func=init_plot,
    interval=dt * ratio_simulation_control * 1000,
    blit=False,
    cache_frame_data=False,
)

plt.show()
