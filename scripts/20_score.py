"""Original light background score for the explainer video, synthesised from scratch with numpy (no samples).

  pads     slow progression F maj9 - D m9 - Bb maj7 - C sus2 (two cycles), soft detuned voices
  bass     chord roots, sine with a little warmth
  plucks   sparse F-major-pentatonic plucks, a gentle "heartbeat" figure
  swell    a soft rise into the closing card
A simple stereo reverb (feedback combs + allpass) glues it together.

Writes video/score.wav (44.1 kHz stereo). Usage: python scripts/20_score.py [DURATION_S]
"""
import sys
import wave
from pathlib import Path

import numpy as np
from scipy.signal import lfilter

SR = 44100
RNG = np.random.default_rng(50)
A4 = 440.0
F3 = -16  # semitones from A4
# chord voicings, semitones from A4
CHORDS = [
    [-16, -9, -5, -1, 3],    # F maj9 (F C E G ... )
    [-19, -12, -7, -2, 0],   # D m9
    [-14, -7, -3, 1, 5],     # Bb maj7 add9
    [-21, -14, -7, -2, 2],   # C sus2
]
ROOTS = [-28, -31, -26, -33]
PENT = [0, 2, 4, 7, 9, 12, 14]  # F major pentatonic, semitones above F4


def hz(s):
    return A4 * 2 ** (s / 12)


def env(n, a, r):
    e = np.ones(n)
    a, r = min(int(a * SR), n // 2), min(int(r * SR), n // 2)
    e[:a] = np.linspace(0, 1, a) ** 2
    e[n - r:] *= np.linspace(1, 0, r) ** 2
    return e


def onepole(x, fc):
    a = 1 - np.exp(-2 * np.pi * fc / SR)
    return lfilter([a], [1, a - 1], x)


def pad_note(f, n):
    t = np.arange(n) / SR
    out = np.zeros(n)
    for det in (-0.08, 0.0, 0.07):
        ff = f * 2 ** (det / 12)
        ph = RNG.uniform(0, 2 * np.pi)
        # a few harmonics with falling weights: soft, string-pad-like
        for k, w in ((1, 1.0), (2, 0.35), (3, 0.15), (4, 0.06)):
            out += w * np.sin(2 * np.pi * ff * k * t + ph * k)
    vib = 1 + 0.15 * np.sin(2 * np.pi * 0.2 * t + RNG.uniform(0, 6))
    return out * vib / 3


def pluck(f, n):
    t = np.arange(n) / SR
    e = np.exp(-t * 3.2) * (1 - np.exp(-t * 400))
    s = np.sin(2 * np.pi * f * t) + 0.25 * np.sin(2 * np.pi * 2 * f * t) * np.exp(-t * 6)
    return s * e


def reverb(x, mix=0.35):
    def comb(x, d, g):  # y[n] = x[n] + g*y[n-d], computed one delay-length block at a time
        y = x.copy()
        for i in range(d, len(x), d):
            j = min(i + d, len(x))
            y[i:j] += g * y[i - d:j - d]
        return y

    def allpass(x, d, g):  # y[n] = -g*x[n] + x[n-d] + g*y[n-d]
        y = -g * x
        y[d:] += x[:-d]
        for i in range(d, len(x), d):
            j = min(i + d, len(x))
            y[i:j] += g * y[i - d:j - d]
        return y

    out = np.zeros_like(x)
    for d, g in ((1557, 0.84), (1617, 0.83), (1491, 0.85), (1422, 0.84), (1277, 0.82), (1356, 0.83)):
        out += comb(x, d, g)
    out /= 6
    for d, g in ((225, 0.5), (556, 0.5)):
        out = allpass(out, d, g)
    out = onepole(out, 5000)
    return (1 - mix) * x + mix * out


def render(duration):
    n = int(duration * SR)
    L = np.zeros(n); R = np.zeros(n)
    chord_s = duration / 8  # two passes of four chords
    for i in range(8):
        c = i % 4
        s0 = int(i * chord_s * SR)
        ln = int((chord_s + 2.4) * SR)  # overlap for smooth changes
        ln = min(ln, n - s0)
        e = env(ln, 2.0, 2.4)
        for j, semi in enumerate(CHORDS[c]):
            v = onepole(pad_note(hz(semi), ln), 1400) * e * 0.11
            pan = 0.3 + 0.4 * (j / (len(CHORDS[c]) - 1))
            L[s0:s0 + ln] += v * (1 - pan); R[s0:s0 + ln] += v * pan
        t = np.arange(ln) / SR
        b = (np.sin(2 * np.pi * hz(ROOTS[c]) * t) + 0.2 * np.sin(4 * np.pi * hz(ROOTS[c]) * t)) * e * 0.16
        L[s0:s0 + ln] += b; R[s0:s0 + ln] += b
    # plucks: a sparse pattern from 6 s, thinning during the "alternative" section, back for the ending
    beat = 60 / 84
    tt = 6.0
    k = 0
    while tt < duration - 2.5:
        dens = 0.55 if 34 < tt < 46 else 0.8
        if RNG.random() < dens:
            c = min(int(tt / chord_s), 7) % 4
            # prefer chord-friendly pentatonic notes
            semi = -4 + PENT[RNG.integers(0, len(PENT))]
            ln = int(2.0 * SR)
            s0 = int(tt * SR)
            if s0 + ln < n:
                p = pluck(hz(semi), ln) * (0.15 + 0.05 * RNG.random())
                pan = RNG.uniform(0.25, 0.75)
                L[s0:s0 + ln] += p * (1 - pan); R[s0:s0 + ln] += p * pan
        tt += beat * (1 if k % 3 else 2)
        k += 1
    # soft noise swell into the end card
    s0 = int((duration - 11) * SR); ln = int(6 * SR)
    sw = onepole(onepole(RNG.standard_normal(ln), 600), 600) * np.linspace(0, 1, ln) ** 2 * 0.03
    sw *= env(ln, 0.1, 1.5)
    L[s0:s0 + ln] += sw; R[s0:s0 + ln] += sw
    # final resolving chord
    s0 = int((duration - 5.5) * SR); ln = n - s0
    e = env(ln, 0.8, 3.5)
    for semi in (-16, -9, -4, 3, 8):
        v = onepole(pad_note(hz(semi), ln), 1800) * e * 0.09
        L[s0:s0 + ln] += v; R[s0:s0 + ln] += v
    L, R = reverb(L), reverb(R * 0.98)
    st = np.stack([L, R], 1)
    st *= env(n, 1.5, 2.5)[:, None]
    st = np.tanh(st * 1.2) / 1.2
    st *= 10 ** (-5 / 20) / np.max(np.abs(st))  # peak at -5 dBFS (~-19 LUFS): gentle but audible on a phone
    return st


if __name__ == "__main__":
    dur = float(sys.argv[1]) if len(sys.argv) > 1 else 56.0
    out = Path(__file__).resolve().parent.parent / "video" / "score.wav"
    out.parent.mkdir(exist_ok=True)
    pcm = (render(dur) * 32767).astype("<i2")
    with wave.open(str(out), "wb") as w:
        w.setnchannels(2); w.setsampwidth(2); w.setframerate(SR)
        w.writeframes(pcm.tobytes())
    print("wrote", out, f"{dur:.1f}s")
