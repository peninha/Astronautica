from graphviz import Digraph
import os

# Define o caminho de saída
output_dir = "../docs/diagrams/"
os.makedirs(output_dir, exist_ok=True)  # Garante que a pasta exista

# Cria o diagrama UML
uml = Digraph(format="png", name="uml_diagram")
uml.attr(rankdir="LR", size="8,5")

# Define classes
uml.node("Body", "Body\n- mass\n- radius\n...")
uml.node("Frame", "Frame\n- origin\n- orientation\n- velocity_offset\n...\n+ to_bodycentric()\n+ from_bodycentric()")
uml.node("Orbit", "Orbit\n- get_points()\n+ state_vectors_at_time()\n...")
uml.node("Trajectory", "Trajectory\n- orbits\n- maneuvers\n- get_points()\n...")
uml.node("Maneuver", "Maneuver\n- delta_v\n- time\n...")
uml.node("Plotter", "Plotter\n- frame\n+ plot_orbit()\n+ plot_trajectory()")

# Define relações
uml.edge("Body", "Orbit", label="defines center for")
uml.edge("Frame", "Orbit", label="transforms frame")
uml.edge("Frame", "Trajectory", label="transforms frame")
uml.edge("Orbit", "Trajectory", label="contains")
uml.edge("Maneuver", "Trajectory", label="contains")
uml.edge("Frame", "Plotter", label="uses for")
uml.edge("Orbit", "Plotter", label="plots")
uml.edge("Trajectory", "Plotter", label="plots")

# Salva o diagrama na pasta definida
uml.render(os.path.join(output_dir, "uml_diagram"), format="png", cleanup=True)
