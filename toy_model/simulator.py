"""
Toy model for peak evolution as a function of time
"""

import numpy as np
import random


class Parameter:
    def __init__(
        self, value: float, rate: float = 0, onset_time: float = 0, end_time: float = 0
    ):
        self.value = value
        self.initial_value = value
        self.rate = rate
        self.onset_time = onset_time
        self.end_time = end_time

    def update(self, time: float):
        """
        TODO: may want to replace this by an error function
        """
        if time >= self.onset_time and time < self.end_time:
            self.value = self.initial_value + self.rate * (time - self.onset_time)
        elif time >= self.end_time:
            self.value = self.initial_value + self.rate * (
                self.end_time - self.onset_time
            )
        else:
            self.value = self.initial_value
        return max(0., self.value)


class Peak:
    def __init__(self, position: Parameter, amplitude: Parameter, width: Parameter):
        # Peak parameters
        self.amplitude = amplitude
        self.width = width
        self.position = position

    def evaluate(self, x):
        return self.amplitude.value * np.exp(
            -0.5 * ((x - self.position.value) / self.width.value) ** 2
        )

    def update(self, time):
        return (
            self.amplitude.update(time),
            self.width.update(time),
            self.position.update(time),
        )


class Pattern:
    def __init__(self, peaks: list, x: np.ndarray):
        self.peaks = peaks
        self.x = x

    def evaluate(self):
        y = np.zeros_like(self.x, dtype=float)
        for peak in self.peaks:
            y += peak.evaluate(self.x)
        return y

    def update(self, time):
        return [peak.update(time) for peak in self.peaks]


def create_pattern(max_time=100, n_peaks=5, n_points=1000):
    MAX_POSITION = 100

    # Pick a transition time
    transition_time = random.uniform(10, max_time)

    # Create peaks
    peaks = []
    for _ in range(n_peaks):
        # Is this peak involved in the phase transision?
        # If so, pass along the transition time
        _transition_time = transition_time if random.random() < 0.5 else 0
        peak = create_peak(max_time, MAX_POSITION, _transition_time)
        peaks.append(peak)

    # Create x values
    x = np.linspace(0, MAX_POSITION, n_points)
    pattern = Pattern(peaks, x)
    return pattern, transition_time




def create_peak(max_time=100, max_position=100, transition_time=0):
    MAX_AMPLITUDE = 20
    MAX_WIDTH = 2
    MAX_EXPANSION = 2 * MAX_WIDTH

    # Amplitude
    amplitude_rate = 0
    end_time = 0
    if transition_time > 0 and random.random() < 0.5:
        delta_t = (max_time - transition_time) * random.uniform(0.1, 0.9)
        amplitude_rate = MAX_AMPLITUDE * random.uniform(-1, 1) / delta_t + 1
        end_time = transition_time + delta_t

    print("TRANSITION TIME", transition_time, end_time, amplitude_rate)

    amplitude = Parameter(random.uniform(0, MAX_AMPLITUDE), amplitude_rate, transition_time, end_time)
    width = Parameter(random.uniform(1, MAX_WIDTH))

    if random.random() < 0.3:
        max_rate = MAX_EXPANSION / max_time
        position = Parameter(
            random.uniform(0, max_position),
            max_rate * random.uniform(-1, 1),
            0,
            max_time,
        )
    else:
        position = Parameter(random.uniform(0, max_position))

    return Peak(position, amplitude, width)