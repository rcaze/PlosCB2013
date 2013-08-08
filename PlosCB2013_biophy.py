# Script mostly written by Romain Caze some part are updated from Andrew Davidson code snippets
# last modification:26/03/13.
from neuron import nrn, h, hclass, run, init
import matplotlib.pyplot as plt
import numpy as np

# Wrappers from Andrew Davidson neuronpy modified and commented by Romain Caze
# They enable to interface NEURON and Python more easily.
class Mechanism(object):
    """
    Create mechanism which will be inserted in the membrane of a cell

    Examples
    --------
    >>> leak = Mechanism('pas', {'e': -65, 'g': 0.0002})
    >>> hh = Mechanism('hh')
    """
    def __init__(self, name, parameters={}):
        """
        Parameters
        ----------
        name: a char
            the label of the mechanism
        parameters: a dictionary
            contains the different parameter of a mechanism
        """
        self.name = name
        self.parameters = parameters

    def insert_into(self, section):
        """
        Method used to insert a mechanism into a section

        Parameters
        ----------
        section: a NEURON section
            the section where the mechanism needs to be inserted

        """
        section.insert(self.name)
        for name, value in self.parameters.items():
            for segment in section:
                mech = getattr(segment, self.name)
                setattr(mech, name, value)

class Section(nrn.Section):
    """
    Create a NEURON section with certain mechanism inserted in it

    Examples
    --------
    >>> soma = Section(L=30, diam=30, mechanisms=[hh, leak])
    >>> apical = Section(L=600, diam=2, nseg=5, mechanisms=[leak],
    ...                  parent=soma, connection_point=0)
    """
    def __init__(self, L, diam, nseg=1, Ra=100, cm=1, mechanisms=[], parent=None, connection_point=1):
        """
        Parameters
        ----------
        L: a float
            length in micrometers
        diam: a float
            diameter in micrometers
        nseg: an int
            number of segements
        Ra: a float
            surface axial resistance Ohm/micrometer square
        cm: a float
            capacitance in F/micrometer square
        mechanisms: a list
            mechanisms to be inserted (need to be created beforehand)
        parent: a NEURON section
            section to which this section is coming from
        connection_point: a float between 0 and 1
            where this section is connected to its parent

        """

        nrn.Section.__init__(self)
        # set geometry
        self.L = L
        self.diam = diam
        self.nseg = nseg
        # set cable properties
        self.Ra = Ra
        self.cm = cm
        # connect to parent section
        if parent:
            self.connect(parent, connection_point, 0)
        # add the mechanisms
        for mechanism in mechanisms:
            mechanism.insert_into(self)

    def record_spikes(self, threshold=-30):
        """
        Record the number of spikes produced in this section,
        which is the number of time a voltage is crossed in
        the middle of a section

        Parameters
        ----------
        threshold: a float
            voltage determining the presence or not of a spike

        Returns
        -------
        nothing, but change the self.spikecount

        """


        self.spiketimes = h.Vector()
        self.spikecount = h.APCount(0.5, sec=self)
        self.spikecount.thresh = threshold
        self.spikecount.record(self.spiketimes)

#My own class inspired by the pydesign tutorial.
#http://www.paedia.info/quickstart/pydesign.html

class BipolarNeuron(object):
    """
    Produce neuron objects with a standard soma
    and two identical dendrites connected on opposite sides of the
    soma. For the dendrites the following parameters can be changed:

    """
    def __init__(self, d_length=50, d_diam=0.4):
        """
        Parameters
        ----------
        d_length: an integer
            length of the dendrites
        d_diameter:
            diameter of each dendrites
        """
        # Creating Mechanisms
        hh = Mechanism('hh')
        pas = Mechanism('pas', {'e':-65,'g':0.0001})

        # Creating the Sections
        self.soma = Section(10, 10, Ra=150, mechanisms=[hh,pas])
        self.Xdend = Section(d_length, d_diam, nseg=10, Ra=150, parent=self.soma, mechanisms=[pas])
        self.Ydend = Section(d_length, d_diam, nseg=10, Ra=150, parent=self.soma, mechanisms=[pas])

    def initialise(self, vrest=-65):
        """
        Initialise the model, to launch before each simulations
        """
        for sec in h.allsec():
            h.finitialize(vrest, sec)
            h.fcurrent(sec)
        h.frecord_init()

    def min_sim(self, TSTOP=100):
        """
        Launch a minimal simulation to test the model and determine its resting potential empirically
        """
        vrec = h.Vector()
        vrec.record(self.soma(0.5)._ref_v)

        for sec in h.allsec():
            h.finitialize(-65, sec)
            h.fcurrent(sec)
        h.frecord_init()

        while h.t < TSTOP: #Launch a simulation
            h.fadvance()

        vrest = np.array(vrec)[-1]

        return vrest

class Simulation(object):
    """
    Create and control a simulation

    Example
    -------
    >>> cell = BipolarNeuron()
    >>> sim = Simulation(cell)
    >>> sim.go()
    >>> sim.show()
    """
    def __init__(self,
                 cell,
                 sim_time=100,
                 dt=0.01):
        """
        Parameters
        ----------
        cell: BipolarNeuron object
            Neuron model to be stimulated
        sim_tim: integer
            Time of the simulation in ms
        protocol_char: list of integers in char form
            Readable description of the stimulation protocol
            'x-y' generate a stimulation episode where site x and y are activated together.
            Site are labelled by positive integers, '0' corresponds to a no-stimulation episode.
        episode_time: an integer
            timelength of an episode in ms
        dt: a float
            the integration timestep (fix)
        """
        self.cell = cell
        self.sim_time = sim_time
        self.dt = dt
        self.syn, self.stim, self.vplay, self.netcon = {}, {}, {}, {}

    def add_ExpSyn(self, section='soma', position=0.5, name='default', tstim=[50], w=.001):
        """
        Create/replace an Expsyn synapse on a given section which is active at the time in tstim

        Comments
        --------
        The sort command is here to make sure that tstim are in the right order. This method
        requires the pre-compiling of vecstim.mod by NEURON.

        """
        self.syn[name] = h.ExpSyn(self.cell.__getattribute__(section)(position))
        self.stim[name] = h.Vector(np.sort(tstim)) # Converting tstim into a NEURON vector (to play in NEURON)
        self.vplay[name] = h.VecStim() # Creating play vectors to interface with NEURON
        self.vplay[name].play(self.stim[name])  # Connecting vector to VecStim object to play them
        self.netcon[name] = h.NetCon(self.vplay[name], self.syn[name]) # Building the netcon object to connect the stims and the synapses
        self.netcon[name].weight[0] = w # Setting the individual weights

    def set_IClamp(self, name='IClamp', delay=1, amp=-1, dur=3):
        """
        Create a current clamp point process on the soma
        """
        stim = h.IClamp(self.cell.soma(0.5))
        stim.delay = delay
        stim.amp = amp
        stim.dur = dur
        self.stim[name] = stim

    def show(self):
        """
        Show the voltage trace after a simulation
        """
        x = np.array(self.rec_t)
        y = np.array(self.rec_v)
        plt.plot(x, y)
        plt.xlabel("Time [ms]")
        plt.ylabel("Voltage [mV]")
        #plt.axis(ymin=-120, ymax=-50)
        plt.show()

    def set_recording(self):
        """
        Set a recording vector to the soma
        """
        # Record Time
        self.rec_t = h.Vector()
        self.rec_t.record(h._ref_t)
        # Record Voltage
        self.rec_v = h.Vector()
        self.rec_v.record(self.cell.soma(0.5)._ref_v)

    def get_recording(self):
        """
        Return the recordings of the somatic voltage and a time axis
        """
        time = np.array(self.rec_t)
        voltage = np.array(self.rec_v)
        return time, voltage

    def go(self, sim_time=None):
        """
        Launch a simulation of a given time

        Parameters
        ----------
        sim_time: an integer
            the time in millisecond of the simulation
            it replaces the self.sim_time if defined

        Comments
        --------
        It seems that when multiple go method are done it does not
        change the output vector.
        """

        h.t = 0
        #Start recording
        self.set_recording()
        h.dt = self.dt
        self.cell.initialise()
        init()
        #while h.t < self.sim_time: #I was using this procedure before do not know which one is better
        #    h.fadvance()
        if sim_time:
            run(sim_time)
        else:
            run(self.sim_time)

    def insert_signal(self,
                       synlocs=[['Xdend',1], ['Xdend',1], ['Ydend',1]],
                       el=0.005,
                       weights=[1,1,1],
                       tstims=[[25,50],[25],[50]]
                      ):
        """
        Add a set of stimulation to a simulation object to reproduce the uncaging protocol we made with Alex.

        Parameters
        ----------
        synlocs: a list of lists made of a char and a number
            sections and positions of the stimulation
        el: a float
            elementary weight value, defining the maximum conductance
        conductance of a synapse in microS
        weights: a list
            multiplicative factors corresponding to the weights
        tstims: a list of lists.
            stimulation times corresponding to each site (one list per site).

        Returns
        -------
        sim: a simulation object
            this simulation is now "decorated" with stimulations
        """
        for i, tstim in enumerate(tstims):
            self.add_ExpSyn(section=synlocs[i][0], position=synlocs[i][1], name='Stream'+str(i), tstim=tstim, w=weights[i]*el)

if __name__ == '__main__':
    mysim = Simulation(BipolarNeuron())
    mysim.insert_signal()
    mysim.go()
    mysim.show()
