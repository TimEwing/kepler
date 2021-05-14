# standard
import math
# third party
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
# local
from geo import pixel_area
from geo import circular_intersect

def get_args():
    args = {}
    args['ecc'] = 0.00 # eccentricity - unitless
    args['sma'] = 9_000_000 # semi major axis - km
    args['smi'] = (1-args['ecc']) * args['sma'] # semi minor axis - km
    args['foc'] = np.sqrt(args['sma']**2 * (1 - args['ecc']**2)) # focal distance - km
    args['ape'] = 0.0 # argument of perihelion - rad
    args['inc'] = 0.01 # inclination - rad

    args['solar_rad'] = 696_000 # solar radius - km
    args['planet_rad'] = 680_000 # planet radius - km

    args['lightcurve_noise'] = 0.001 # std dev of noise; noise is a normal multiplied by total solar flux
    args['pixels'] = (8,8)

    return args

def x(t, args):
    # a*cos(t) + c
    return args['sma'] * np.cos(t) + args['ecc']*args['sma']
def y(t, args):
    # b*sin(t)
    return args['smi'] * np.sin(t)
# Rotate x and y by argument of perihelion
def x_rot(t, args):
    # x*cos(w) - y*sin(w)
    return x(t, args) * np.cos(args['ape']) + y(t, args) * np.sin(args['ape'])
def y_rot(t, args):
    # x*sin(w) + y*cos(w)
    return -x(t, args) * np.sin(args['ape']) + y(t, args) * np.cos(args['ape'])

# Projection, assuming earth is along the +x direction
def projection_y(t, args):
    # y * sin(inc)
    return x_rot(t, args) * np.sin(args['inc'])
def projection_x(t, args):
    # x
    return y_rot(t, args)

# Numerically solve for transit start and end
def solve(t0, t1, f, target, precision=0.001, max_itr=100):
    if f(t0) > f(t1):
        t0, t1 = t1, t0
    for _ in range(max_itr):
        t = (t0 + t1) / 2
        v0, v1, v = f(t0), f(t1), f(t)
        if v0 <= target <= v:
            t1 = t
        elif v <= target <= v1:
            t0 = t
        else:
            break

        if abs(f(t) - target) <= precision:
            return t
    raise ValueError("No solution for target {}".format(target))
def find_transit(args, precision=0.001):
    def f(t):
        return projection_x(t, args)
    t1 = solve(
        t0=0, 
        t1=np.pi - 0.001, 
        f=f, 
        target=0,
    )
    t2 = solve(
        t0=np.pi,
        t1=2*np.pi - 0.001,
        f=f, 
        target=0,
    )
    if x_rot(t1, args) > 0:
        t = t1
    elif x_rot(t2) > 0:
        t = t2
    else:
        raise ValueError("No transit: {} {}".format(t1, t2))

    window_start = solve(
        t0=t - np.pi/4,
        t1=t + np.pi/4,
        f=f,
        target=-args['solar_rad'] - args['planet_rad'],
    )
    window_end = solve(
        t0=t - np.pi/4,
        t1=t + np.pi/4,
        f=f,
        target=args['solar_rad'] + args['planet_rad'],
    )
    if window_start > window_end:
        window_start, window_end = window_end, window_start
    return window_start, window_end

# Solar radii = sr
def sr_to_km(x):
    return 695_510 * x
def km_to_sr(x):
    return x / 695_510

def show_projection_plot():

    plt.rcParams["axes.edgecolor"] = "0.15"
    plt.rcParams["axes.linewidth"]  = 1.25

    args = get_args()
    window_start, window_end = find_transit(args)
    times = np.arange(window_start, window_end, (window_end - window_start) * 0.005)

    # Setup plots
    fig = plt.figure(figsize=(10, 5))
    ax_data = plt.subplot2grid((2,3), (0,0))
    ax_transit = plt.subplot2grid((2,3), (0,1))
    ax_topdown = plt.subplot2grid((2,3), (0,2))
    ax_lightcurve = plt.subplot2grid((2,3,), (1,0), colspan=3)

    # Pixel data (sensor simulation)
    pix_width = 2 * args['solar_rad'] / args['pixels'][0]
    pix_height = 2 * args['solar_rad'] / args['pixels'][1]
    pix_max = pix_width * pix_height
    data_time_series = []
    for t in times:
        data = np.zeros(args['pixels'])
        for x in range(args['pixels'][0]):
            for y in range(args['pixels'][1]):
                area = pixel_area(
                    pix_x0=x * pix_width - args['solar_rad'],
                    pix_y0=y * pix_width - args['solar_rad'],
                    pix_x1=(x+1) * pix_width - args['solar_rad'],
                    pix_y1=(y+1) * pix_height - args['solar_rad'],
                    sol_x=0,
                    sol_y=0,
                    sol_r=args['solar_rad'],
                    obj_x=projection_x(t, args),
                    obj_y=projection_y(t, args),
                    obj_r=args['planet_rad'],
                )
                data[y,x] = int((area / pix_max) * 255)
        data_time_series.append(data)
    data_image = ax_data.imshow(data_time_series[0], cmap='hot')
    ax_data.set_xticks([])
    ax_data.set_yticks([])

    # Transit
    # Plot styling stuff
    ax_transit.set_facecolor('black')
    ax_transit.set_xlim(-args['solar_rad'] * 1.1, args['solar_rad'] * 1.1)
    ax_transit.set_ylim(-args['solar_rad'] * 1.1, args['solar_rad'] * 1.1)
    ax_transit.set_aspect(1)
    # Planet and stun
    transit_planet = plt.Circle((0,0), args['planet_rad'], color='black', zorder=1)
    transit_sun = plt.Circle((0, 0), args['solar_rad'], color='w', zorder=0)
    ax_transit.add_artist(transit_planet)
    ax_transit.add_artist(transit_sun)
    # Axis tick stuff
    ax_transit.set_xticks([])
    ax_transit.set_yticks([])

    # Topdown plot stuff
    topdown_times = np.arange(0, 2*np.pi, 0.01)
    ax_topdown.plot(
        x_rot(topdown_times, args), 
        y_rot(topdown_times, args),
        linewidth=1,
        color='blue',
    )
    topdown_planet = plt.Circle(
        (0, 0), 
        args['planet_rad'], 
        color='black'
    )
    topdown_sun = plt.Circle(
        (0, 0), 
        args['solar_rad'], 
        color='#ffee44'
    )
    topdown_rect = plt.Rectangle(
        (0, -args['solar_rad']), 
        20*args['sma'], # Some crazy big number
        2*args['solar_rad'],
        color='lightgrey'
    )
    ax_topdown.add_artist(topdown_rect)
    ax_topdown.add_artist(topdown_planet)
    ax_topdown.add_artist(topdown_sun)
    ax_topdown.set_aspect(1)
    ax_topdown.set_xticks([])
    ax_topdown.set_yticks([])

    # lightcurve stuff
    lightcurve_window_start = window_start - (window_end - window_start) * 0.1
    lightcurve_window_end = window_end + (window_end - window_start) * 0.1
    lightcurve_time = np.arange(
        lightcurve_window_start, 
        lightcurve_window_end, 
        (window_end - window_start) * 0.005
    )
    lightcurve_data = []
    for t in lightcurve_time:
        solar_area = pixel_area(
            pix_x0=-args['solar_rad'],
            pix_x1=args['solar_rad'],
            pix_y0=-args['solar_rad'],
            pix_y1=args['solar_rad'],
            sol_x=0,
            sol_y=0,
            sol_r=args['solar_rad'],
            obj_x=projection_x(t, args),
            obj_y=projection_y(t, args),
            obj_r=args['planet_rad'],
        )
        lightcurve_data.append(solar_area)
    lightcurve_data = np.array(lightcurve_data)
    lightcurve_noise = np.random.normal(0, 1, lightcurve_data.shape)
    lightcurve_noise *= np.pi * args['solar_rad']**2
    lightcurve_noise *= args['lightcurve_noise']
    lightcurve_data += lightcurve_noise
    ax_lightcurve.plot(lightcurve_time, lightcurve_data)
    lightcurve_line = ax_lightcurve.axvline(0, color='red')
    lightcurve_line.set_ydata([0, 0.985])
    ax_lightcurve.set_xticks([])
    ax_lightcurve.set_yticks([])
    ax_lightcurve.set_xlim([lightcurve_window_start, lightcurve_window_end])


    def init():
        return transit_planet, transit_sun, topdown_sun, topdown_planet, data_image, lightcurve_line

    def animate(t):
        i,t = t
        transit_planet.center = projection_x(t, args), -projection_y(t, args) # Flip transit about x axis to match data
        topdown_planet.center = x_rot(t, args), y_rot(t, args)
        data_image.set_data(data_time_series[i])
        lightcurve_line.set_xdata(t)
        return transit_planet, transit_sun, topdown_sun, topdown_planet, data_image, lightcurve_line

    ani = animation.FuncAnimation(
        fig, animate, 
        frames=list(enumerate(times)), 
        blit=True, 
        init_func=init, 
        interval=25
    )
    # fig.suptitle("Transit of Jupiter-sized planet in Mercury orbit", fontsize=24)
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=25, codec='mpeg4')
    ani.save('output.mp4', writer=writer)
    plt.show()


if __name__ == '__main__':
    show_projection_plot()
    # show_topdown_plot()