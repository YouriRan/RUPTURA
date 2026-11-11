"""Render the column viewer's curves as a fixed-axis H.264 animation."""
from __future__ import annotations

import html
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import threading
import uuid

import numpy as np

from .plots import figure, metadata

_RENDER_LOCK = threading.Lock()


def encoder():
    executable = shutil.which("ffmpeg")
    if executable:
        return executable
    try:
        from imageio_ffmpeg import get_ffmpeg_exe
        return get_ffmpeg_exe()
    except (ImportError, RuntimeError) as error:
        raise ValueError("MP4 export needs FFmpeg; install ruptura[web] or add ffmpeg to PATH.") from error


def plain(value):
    return html.unescape(re.sub(r"<[^>]+>", "", value or ""))


def export_mp4(directory, options):
    fps = options.get("fps", 20)
    if isinstance(fps, bool) or not isinstance(fps, (int, float)) or not np.isfinite(fps) or fps != int(fps) or not 1 <= fps <= 60:
        raise ValueError("Frames per second must be a whole number from 1 to 60.")
    axis = options.get("xAxis", "position")
    if axis not in ("position", "time"):
        raise ValueError("Choose time or position for the column x-axis.")
    directory = Path(directory)
    meta = metadata(directory)
    if meta["kind"] == "MixturePrediction":
        raise ValueError("MP4 export is available for column profiles.")
    executable = encoder()
    try:
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from matplotlib.figure import Figure
    except ImportError as error:
        raise ValueError("MP4 export needs Matplotlib; install ruptura[web] in the Lab Python environment.") from error
    count = meta["frames"] if axis == "position" else len(meta["positions"])
    base = {**options, "plot": "profile", "xAxis": axis}

    def chart(frame):
        return figure(directory, {**base, "frame" if axis == "position" else "grid": frame})

    # Build every frame once and keep only what the renderer draws; fix both axes across the video so changing
    # frames cannot rescale the data.
    frames = []
    limits = {"x": [np.inf, -np.inf], "y": [np.inf, -np.inf]}
    for frame in range(count):
        plot = chart(frame)
        traces = [(np.asarray(trace.x, dtype=float), np.asarray(trace.y, dtype=float), plain(trace.name))
                  for trace in plot.data]
        for trace in traces:
            for key, values in zip(("x", "y"), trace):
                values = values[np.isfinite(values)]
                if values.size:
                    limits[key][0] = min(limits[key][0], float(values.min()))
                    limits[key][1] = max(limits[key][1], float(values.max()))
        vertical_lines = [shape.x0 for shape in plot.layout.shapes or () if shape.type == "line" and shape.x0 == shape.x1]
        frames.append((traces, vertical_lines, plain(plot.layout.xaxis.title.text),
                       plain(plot.layout.yaxis.title.text), plain(plot.layout.title.text)))
    if not all(np.isfinite(value) for pair in limits.values() for value in pair):
        raise ValueError("Select at least one component with finite data before exporting an MP4.")
    layout = plot.layout  # theme colours are the same in every frame
    for key, pair in limits.items():
        padding = (pair[1] - pair[0]) * (0.05 if key == "y" else 0)
        if pair[0] == pair[1]:
            padding = max(abs(pair[0]) * 0.05, 0.5)
        pair[0] -= padding
        pair[1] += padding
    target = directory / f"column-profiles-{int(fps)}fps-{uuid.uuid4().hex[:8]}.mp4"
    with tempfile.TemporaryDirectory(prefix=".video-", dir=directory) as temporary, _RENDER_LOCK:
        output = Path(temporary) / "animation.mp4"
        command = [executable, "-hide_banner", "-loglevel", "error", "-y", "-f", "rawvideo", "-pixel_format", "rgba",
                   "-video_size", "1280x720", "-framerate", str(int(fps)), "-i", "pipe:0", "-an", "-c:v", "libx264",
                   "-preset", "fast", "-crf", "20", "-pix_fmt", "yuv420p", "-movflags", "+faststart", str(output)]
        with tempfile.TemporaryFile() as errors:
            process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, stderr=errors)
            try:
                for frame, (traces, vertical_lines, xlabel, ylabel, title) in enumerate(frames):
                    light = options.get("theme") == "light"
                    ink = "#182e3b" if light else "#e6edf3"
                    canvas_figure = Figure(figsize=(12.8, 7.2), dpi=100, facecolor=layout.paper_bgcolor)
                    canvas = FigureCanvasAgg(canvas_figure)
                    axes = canvas_figure.add_subplot(111)
                    axes.set_facecolor(layout.plot_bgcolor)
                    colors = layout.colorway
                    for i, (x, y, name) in enumerate(traces):
                        axes.plot(x, y, label=name, color=colors[i % len(colors)], linewidth=2,
                                  marker="o" if options.get("markers") else None, markersize=3)
                    for x0 in vertical_lines:
                        axes.axvline(x0, color=ink, linestyle=":", alpha=0.5)
                    axes.set(xlim=limits["x"], ylim=limits["y"], xlabel=xlabel, ylabel=ylabel, title=title)
                    axes.tick_params(colors=ink)
                    axes.xaxis.label.set_color(ink)
                    axes.yaxis.label.set_color(ink)
                    axes.title.set_color(ink)
                    for spine in axes.spines.values():
                        spine.set_color(ink)
                    axes.grid(alpha=0.2, color=ink)
                    axes.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=min(4, len(traces)),
                                facecolor=layout.paper_bgcolor, labelcolor=ink)
                    units = options.get("xUnits", "min")
                    factor = 60 if units == "s" else 1 / 60 if units == "hr" else 1
                    label = (f"Time: {meta['times'][frame] * factor:.6g} {units}" if axis == "position"
                             else f"Position: {meta['positions'][frame]:.6g} m")
                    canvas_figure.text(0.96, 0.96, label, ha="right", color=ink)
                    canvas_figure.subplots_adjust(left=0.12, right=0.96, top=0.85, bottom=0.22)
                    canvas.draw()
                    process.stdin.write(canvas.buffer_rgba().tobytes())
                    canvas_figure.clear()
                process.stdin.close()
                if process.wait(timeout=60):
                    errors.seek(0)
                    raise ValueError("MP4 encoding failed: " + errors.read().decode(errors="replace")[-2000:])
            except BaseException:
                if process.poll() is None:
                    process.kill()
                process.wait()
                if not process.stdin.closed:
                    try:
                        process.stdin.close()
                    except OSError:
                        pass
                raise
        output.replace(target)
    return target
