import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import tkinter as tk
from tkinter import ttk

# Signal generation functions
def generate_step_signal(length, step_point):
    signal = np.zeros(length)
    signal[step_point:] = 1
    return signal

def generate_triangular_signal(length, peak_position, peak_height):
    signal = np.zeros(length)
    for i in range(peak_position + 1):
        signal[i] = (peak_height / peak_position) * i
    for i in range(peak_position, length):
        signal[i] = peak_height - ((peak_height / (length - peak_position - 1)) * (i - peak_position))
    return signal

def generate_dirac_impulse(length, impulse_position):
    signal = np.zeros(length)
    signal[impulse_position] = 1
    return signal

def generate_sinusoid_signal(length, amplitude, frequency, phase):
    t = np.arange(length)
    signal = amplitude * np.sin(2 * np.pi * frequency * t + phase)
    return signal

def generate_phasor_signal(length, amplitude, frequency, phase):
    t = np.arange(length)
    signal = amplitude * np.exp(1j * (2 * np.pi * frequency * t + phase))
    return signal

def generate_rectangular_signal(length, start, width, height):
    signal = np.zeros(length)
    end = min(start + width, length)
    signal[start:end] = height
    return signal

def generate_signum_signal(length):
    t = np.arange(-length // 2, length // 2)
    signal = np.sign(t)
    return signal

def generate_cosinusoid_signal(length, amplitude, frequency, phase):
    t = np.arange(length)
    signal = amplitude * np.cos(2 * np.pi * frequency * t + phase)
    return signal

def generate_monolateral_exponential_signal(length, rate):
    t = np.arange(length)
    signal = np.exp(-rate * t)
    return signal

# Function to perform convolution and animate
def animate_convolution(signal1, signal2, save_as_gif=False, gif_filename="convolution_animation.gif"):
    conv_result = np.convolve(signal1, signal2, mode='full')
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 10))
    fig.suptitle("Convolution of Two Signals")

    ax1.plot(np.real(signal1), label="Signal 1", color="blue")
    ax1.legend()
    ax1.set_ylabel("Amplitude")
    ax1.grid()

    ax2.plot(np.real(signal2), label="Signal 2", color="orange")
    ax2.legend()
    ax2.set_ylabel("Amplitude")
    ax2.grid()

    ax3.set_xlim(0, len(conv_result))
    ax3.set_ylim(min(conv_result) - 1, max(conv_result) + 1)
    ax3.set_ylabel("Convolution Result")
    ax3.set_xlabel("Time")
    line, = ax3.plot([], [], color="purple")
    annotation = ax3.text(0.8, 0.9, "", transform=ax3.transAxes)

    expression_text = ax3.text(0.05, 0.9, "", transform=ax3.transAxes, fontsize=12)

    def animate(i):
        if i < len(conv_result):
            line.set_data(range(i), conv_result[:i])
            annotation.set_text(f"y[{i}] = {conv_result[i]:.2f}")
            expression_text.set_text(f"y[n] = sum(x[k] * h[n-k])")
        return line, annotation, expression_text

    ani = FuncAnimation(fig, animate, frames=len(conv_result), interval=50, blit=True)

    if save_as_gif:
        ani.save(gif_filename, writer=PillowWriter(fps=20))
        print(f"Animation saved as {gif_filename}")

    plt.tight_layout()
    plt.show()

# GUI Setup
class SignalApp:
    def __init__(self, root):
        self.root = root
        root.title("Signal Convolution GUI")

        self.signal_types = {
            "Step (Gradino)": self.configure_step,
            "Triangular": self.configure_triangular,
            "Dirac Impulse": self.configure_dirac,
            "Sinusoid": self.configure_sinusoid,
            "Phasor": self.configure_phasor,
            "Rectangular": self.configure_rectangular,
            "Signum": self.configure_signum,
            "Cosinusoid": self.configure_cosinusoid,
            "Monolateral Exponential": self.configure_monolateral_exponential
        }

        # Signal selection and parameters
        self.selected_signal1 = tk.StringVar(value="Step (Gradino)")
        self.selected_signal2 = tk.StringVar(value="Step (Gradino)")
        self.params1 = {}
        self.params2 = {}

        # GUI Elements
        tk.Label(root, text="Select Signal 1:").grid(row=0, column=0)
        tk.Label(root, text="Select Signal 2:").grid(row=0, column=1)
        
        self.dropdown1 = ttk.Combobox(root, textvariable=self.selected_signal1, values=list(self.signal_types.keys()))
        self.dropdown2 = ttk.Combobox(root, textvariable=self.selected_signal2, values=list(self.signal_types.keys()))
        
        self.dropdown1.grid(row=1, column=0)
        self.dropdown2.grid(row=1, column=1)
        
        self.param_frame1 = tk.Frame(root)
        self.param_frame1.grid(row=2, column=0, pady=10)
        
        self.param_frame2 = tk.Frame(root)
        self.param_frame2.grid(row=2, column=1, pady=10)

        self.dropdown1.bind("<<ComboboxSelected>>", self.update_params1)
        self.dropdown2.bind("<<ComboboxSelected>>", self.update_params2)

        # Convolution button
        tk.Button(root, text="Run Convolution", command=self.run_convolution).grid(row=3, column=0, columnspan=2, pady=20)

        # Initial update of parameters
        self.update_params1()
        self.update_params2()

    def update_params1(self, *args):
        for widget in self.param_frame1.winfo_children():
            widget.destroy()
        signal_func = self.signal_types[self.selected_signal1.get()]
        self.params1 = signal_func(self.param_frame1)

    def update_params2(self, *args):
        for widget in self.param_frame2.winfo_children():
            widget.destroy()
        signal_func = self.signal_types[self.selected_signal2.get()]
        self.params2 = signal_func(self.param_frame2)

    # Functions to configure each signal's parameters
    def configure_step(self, frame):
        params = {"length": 20, "step_point": 5}
        tk.Label(frame, text="Length").grid(row=0, column=0)
        tk.Entry(frame, textvariable=tk.IntVar(value=params["length"])).grid(row=0, column=1)
        tk.Label(frame, text="Step Point").grid(row=1, column=0)
        tk.Entry(frame, textvariable=tk.IntVar(value=params["step_point"])).grid(row=1, column=1)
        return params

    def configure_triangular(self, frame):
        params = {"length": 20, "peak_position": 10, "peak_height": 1}
        tk.Label(frame, text="Length").grid(row=0, column=0)
        tk.Entry(frame, textvariable=tk.IntVar(value=params["length"])).grid(row=0, column=1)
        tk.Label(frame, text="Peak Position").grid(row=1, column=0)
        tk.Entry(frame, textvariable=tk.IntVar(value=params["peak_position"])).grid(row=1, column=1)
        tk.Label(frame, text="Peak Height").grid(row=2, column=0)
        tk.Entry(frame, textvariable=tk.DoubleVar(value=params["peak_height"])).grid(row=2, column=1)
        return params

    def configure_dirac(self, frame):
        params = {"length": 20, "impulse_position": 8}
        tk.Label(frame, text="Length").grid(row=0, column=0)
        tk.Entry(frame, textvariable=tk.IntVar(value=params["length"])).grid(row=0, column=1)
        tk.Label(frame, text="Impulse Position").grid(row=1, column=0)
        tk.Entry(frame, textvariable=tk.IntVar(value=params["impulse_position"])).grid(row=1, column=1)
        return params

    def configure_sinusoid(self, frame):
        params = {"length": 20, "amplitude": 1, "frequency": 0.1, "phase": 0}
        tk.Label(frame, text="Length").grid(row=0, column=0)
        tk.Entry(frame, textvariable=tk.IntVar(value=params["length"])).grid(row=0, column=1)
        tk.Label(frame, text="Amplitude").grid(row=1, column=0)
        tk.Entry(frame, textvariable=tk.DoubleVar(value=params["amplitude"])).grid(row=1, column=1)
        tk.Label(frame, text="Frequency").grid(row=2, column=0)
        tk.Entry(frame, textvariable=tk.DoubleVar(value=params["frequency"])).grid(row=2, column=1)
        tk.Label(frame, text="Phase").grid(row=3, column=0)
        tk.Entry(frame, textvariable=tk.DoubleVar(value=params["phase"])).grid(row=3, column=1)
        return params

    def configure_phasor(self, frame):
        return self.configure_sinusoid(frame) # Reuse parameters

    def configure_rectangular(self, frame):
        params = {"length": 20, "start": 5, "width": 5, "height": 1}
        tk.Label(frame, text="Length").grid(row=0, column=0)
        tk.Entry(frame, textvariable=tk.IntVar(value=params["length"])).grid(row=0, column=1)
        tk.Label(frame, text="Start").grid(row=1, column=0)
        tk.Entry(frame, textvariable=tk.IntVar(value=params["start"])).grid(row=1, column=1)
        tk.Label(frame, text="Width").grid(row=2, column=0)
        tk.Entry(frame, textvariable=tk.IntVar(value=params["width"])).grid(row=2, column=1)
        tk.Label(frame, text="Height").grid(row=3, column=0)
        tk.Entry(frame, textvariable=tk.DoubleVar(value=params["height"])).grid(row=3, column=1)
        return params

    def configure_signum(self, frame):
        params = {"length": 20}
        tk.Label(frame, text="Length").grid(row=0, column=0)
        tk.Entry(frame, textvariable=tk.IntVar(value=params["length"])).grid(row=0, column=1)
        return params

    def configure_cosinusoid(self, frame):
        return self.configure_sinusoid(frame)

    def configure_monolateral_exponential(self, frame):
        params = {"length": 20, "rate": 0.1}
        tk.Label(frame, text="Length").grid(row=0, column=0)
        tk.Entry(frame, textvariable=tk.IntVar(value=params["length"])).grid(row=0, column=1)
        tk.Label(frame, text="Rate").grid(row=1, column=0)
        tk.Entry(frame, textvariable=tk.DoubleVar(value=params["rate"])).grid(row=1, column=1)
        return params

    def generate_signal(self, signal_type, params):
        if signal_type == "Step (Gradino)":
            return generate_step_signal(params["length"], params["step_point"])
        elif signal_type == "Triangular":
            return generate_triangular_signal(params["length"], params["peak_position"], params["peak_height"])
        elif signal_type == "Dirac Impulse":
            return generate_dirac_impulse(params["length"], params["impulse_position"])
        elif signal_type == "Sinusoid":
            return generate_sinusoid_signal(params["length"], params["amplitude"], params["frequency"], params["phase"])
        elif signal_type == "Phasor":
            return generate_phasor_signal(params["length"], params["amplitude"], params["frequency"], params["phase"])
        elif signal_type == "Rectangular":
            return generate_rectangular_signal(params["length"], params["start"], params["width"], params["height"])
        elif signal_type == "Signum":
            return generate_signum_signal(params["length"])
        elif signal_type == "Cosinusoid":
            return generate_cosinusoid_signal(params["length"], params["amplitude"], params["frequency"], params["phase"])
        elif signal_type == "Monolateral Exponential":
            return generate_monolateral_exponential_signal(params["length"], params["rate"])

    def run_convolution(self):
        signal1 = self.generate_signal(self.selected_signal1.get(), self.params1)
        signal2 = self.generate_signal(self.selected_signal2.get(), self.params2)
        animate_convolution(signal1, signal2)

# Main application loop
root = tk.Tk()
app = SignalApp(root)
root.mainloop()
