
def write_vert(file, nVert, dimensions, shift = [0, 0, 0]):
    stepX = (dimensions["x"]/(nVert["x"] - 1.0))
    stepY = (dimensions["y"]/(nVert["y"] - 1.0))
    stepZ = (dimensions["z"]/(nVert["z"] - 1.0))
    for z in range(nVert["z"]):
        for y in range(nVert["y"]):
            for x in range(nVert["x"]):
                file.write(str(x * stepX - shift[0]) + " " + str(y * stepY - shift[1]) + " " + str(z * stepZ - shift[2]) + "\n")


def write_cells(file, nVert):
    for i in range(int(nVert["x"]) - 1):
        for j in range(int(nVert["y"]) - 1):
            for k in range(int(nVert["z"]) - 1):

                file.write(
                    "8 " +
                    str(i + j * nVert["x"] + k * (nVert["x"] * nVert["y"])) + " " +
                    str(i+1 + j * nVert["x"] + k * (nVert["x"] * nVert["y"])) + " " +
                    str(i+1 + (j + 1) * nVert["x"] + k * (nVert["x"] * nVert["y"])) + " " +
                    str(i + (j + 1) * nVert["x"] + k * (nVert["x"] * nVert["y"])) + " " +
                    str(i + j * nVert["x"] + (k+1) * (nVert["x"] * nVert["y"])) + " " +
                    str(i + 1 + j * nVert["x"] + (k+1) * (nVert["x"] * nVert["y"])) + " " +
                    str(i + (j + 1) * nVert["x"] + 1 + (k+1) * (nVert["x"] * nVert["y"])) + " " +
                    str(i + (j + 1) * nVert["x"] + (k+1) * (nVert["x"] * nVert["y"])) + "\n"
                )

def generate_vtk_box(name, dimension, number_of_vertices):
    file = open(name, "w");

    file.write(
        """# vtk DataFile Version 2.0
3D test mesh
ASCII
DATASET UNSTRUCTURED_GRID
POINTS """)

    file.write(str(numberOfVertices["x"] *
                   numberOfVertices["y"] *
                   numberOfVertices["z"]) + " double\n")
    write_vert(file, numberOfVertices, dimensions)

    file.write("\nCELLS ")
    nCells = (numberOfVertices["x"] - 1) * (numberOfVertices["y"] - 1) * (numberOfVertices["z"] - 1)
    file.write(str(nCells) + " " + str(nCells * 9) + "\n")

    write_cells(file, numberOfVertices)

    file.write("\nCELL_TYPES " + str(nCells) + "\n")
    for i in range(nCells):
        file.write("12\n")

if __name__ == "__main__":
    dimensions = {"x":0.5, "y":0.5, "z":1}

    for num_cells in [9]:
        numberOfVertices = {"x":num_cells + 1, "y":num_cells + 1, "z":num_cells+1}
        generate_vtk_box(f"mesh3D.vtk", dimensions, numberOfVertices)
