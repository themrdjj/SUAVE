# Propeller_Plots.py
#
# Created: Mar 2021, R. Erhard
# Modified:

from SUAVE.Core import Units
import pylab as plt
import numpy as np
import matplotlib


def plot_propeller_performance(prop,outputs,conditions):
    # Plots local velocities, blade angles, and blade loading
    
    # Setting Latex Font style
    font = {'family' : 'normal',
            'weight' : 'normal',
            'size'   :  22}
    
    matplotlib.rc('font', **font)
    matplotlib.rcParams['mathtext.fontset'] = 'stix'
    matplotlib.rcParams['font.family'] = 'STIXGeneral'
    matplotlib.rc('lines', lw=3)
    
    # plot results corresponding to the ith control point
    i = 0
    
    # extracting outputs for plotting
    psi         = outputs.disc_azimuthal_distribution[i,:,:]
    r           = outputs.disc_radial_distribution[i,:,:]
    Va_V_disc   = outputs.disc_axial_velocity[i]/outputs.velocity[i][0]
    Vt_V_disc   = outputs.disc_tangential_velocity[i]/outputs.velocity[i][0]
    thrust_disc = outputs.disc_thrust_distribution[i]
    torque_disc = outputs.disc_torque_distribution[i]
    
    delta_alpha = (outputs.disc_local_angle_of_attack[i] - 
                   conditions.aerodynamics.angle_of_attack[i,0])/Units.deg  
    
    # completing the full revolution for the disc
    psi          = np.append(psi,np.array([np.ones_like(psi[0])*2*np.pi]),axis=0)
    r            = np.append(r,np.array([r[0]]),axis=0)
    Va_V_disc    = np.append(Va_V_disc,np.array([Va_V_disc[0]]),axis=0)
    Vt_V_disc    = np.append(Vt_V_disc,np.array([Vt_V_disc[0]]),axis=0)
    thrust_disc  = np.append(thrust_disc,np.array([thrust_disc[0]]),axis=0)
    torque_disc  = np.append(torque_disc,np.array([torque_disc[0]]),axis=0)
    delta_alpha  = np.append(delta_alpha,np.array([delta_alpha[0]]),axis=0)
    
    # adjusting so that the hub is included in the plot:
    rh = prop.hub_radius
    
    fig0, axis0 = plt.subplots(subplot_kw=dict(projection='polar'))
    CS_0 = axis0.contourf(psi, r, Va_V_disc,20,cmap=plt.cm.jet)
    cbar0 = plt.colorbar(CS_0, ax=axis0, format=matplotlib.ticker.FormatStrFormatter('%.2f'))
    cbar0.ax.set_ylabel('$\dfrac{V_a}{V_\infty}$',rotation=0,labelpad=25)
    axis0.set_title('Axial Velocity of Propeller',pad=15) 
    thetaticks = np.arange(0,360,45)
    axis0.set_thetagrids(thetaticks)   
    axis0.set_rorigin(-rh)
    axis0.set_yticklabels([])
    
    fig0, axis0 = plt.subplots(subplot_kw=dict(projection='polar'))
    CS_0 = axis0.contourf(psi, r, Vt_V_disc,20,cmap=plt.cm.jet)
    cbar0 = plt.colorbar(CS_0, ax=axis0,format=matplotlib.ticker.FormatStrFormatter('%.2f'))
    cbar0.ax.set_ylabel('$\dfrac{V_t}{V_\infty}$',rotation=0,labelpad=25)
    axis0.set_title('Tangential Velocity of Propeller',pad=15)   
    thetaticks = np.arange(0,360,45)
    axis0.set_thetagrids(thetaticks)
    axis0.set_rorigin(-rh)
    axis0.set_yticklabels([])
    
    fig0, axis0 = plt.subplots(subplot_kw=dict(projection='polar'))
    CS_0 = axis0.contourf(psi, r, thrust_disc,20,cmap=plt.cm.jet)
    cbar0 = plt.colorbar(CS_0, ax=axis0, format=matplotlib.ticker.FormatStrFormatter('%.2f'))
    cbar0.ax.set_ylabel('Thrust (N)',labelpad=25)
    axis0.set_title('Thrust Distribution of Propeller',pad=15)  
    axis0.set_rorigin(-rh)
    axis0.set_yticklabels([])
    
    fig0, axis0 = plt.subplots(subplot_kw=dict(projection='polar'))
    CS_0 = axis0.contourf(psi, r, torque_disc,20,cmap=plt.cm.jet)
    cbar0 = plt.colorbar(CS_0, ax=axis0, format=matplotlib.ticker.FormatStrFormatter('%.2f'))
    cbar0.ax.set_ylabel('Torque (Nm)',labelpad=25)
    axis0.set_title('Torque Distribution of Propeller',pad=15) 
    axis0.set_rorigin(-rh)
    axis0.set_yticklabels([])
    
    fig0, axis0 = plt.subplots(subplot_kw=dict(projection='polar'))
    CS_0 = axis0.contourf(psi, r, delta_alpha,10,cmap=plt.cm.jet)
    cbar0 = plt.colorbar(CS_0, ax=axis0, format=matplotlib.ticker.FormatStrFormatter('%.2f'))
    cbar0.ax.set_ylabel('Angle of Attack (deg)',labelpad=25)
    axis0.set_title('Blade Local Angle of Attack',pad=15) 
    axis0.set_rorigin(-rh)
    axis0.set_yticklabels([])    
    
    return    


def plot_propeller_geometry(prop, save_figure = False, save_filename = "Propeller Geometry, " , file_type = ".png"):
    # Propellers contains the propellers we wish to plot the blade geometry for.
    
    source = prop.tag
    r_R = prop.radius_distribution / prop.tip_radius
    t_b = prop.max_thickness_distribution
    b_D = prop.chord_distribution / (2*prop.tip_radius)
    beta = prop.twist_distribution / Units.deg
    
    fig, ax1 = plt.subplots()
    fig.set_size_inches(10, 6)
    ax2 = ax1.twinx()  # instantiate a second axis that shares the same x-axis

    ax1.plot( r_R, t_b, "--", color = "tab:gray", label=source + ', t_b')
    ax1.plot( r_R, b_D, "--", color = "tab:blue", label=source + ', b_D')
    ax2.plot( r_R, beta, "--", color= "tab:red", label=source + ', $\\beta$')
    
    ax1.set_xlabel("Blade Fractional Radius (r/R)")
    ax1.set_ylabel("Chord and Thickness Distributions ($b/D$, $t/b$)")
    ax1.set_title("Blade Geometry") 
    
    ax2.set_ylabel("Twist Distribution ($\\beta$)")
    
    ax1.grid()
    fig.legend(loc="upper right", bbox_to_anchor=(1, 1), bbox_transform=ax1.transAxes)
    fig.tight_layout()
    #plt.show()
    
    if save_figure:
        save_filename = save_filename + source
        plt.savefig(save_filename + file_type) 
        
    return
