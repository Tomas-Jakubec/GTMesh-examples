# A simple script generating block or rectangle
def write_vert(file, nVert, dimensions, shift = [0, 0, 0]):
    stepX = (dimensions["x"]/(nVert["x"] - 1.0))
    stepY = (dimensions["y"]/(nVert["y"] - 1.0))
    stepZ = (dimensions["z"]/(nVert["z"] - 1.0))
    for z in range(nVert["z"]):
        for y in range(nVert["y"]):
            for x in range(nVert["x"]):
                file.write(str(x * stepX - shift[0]) + " " + str(y * stepY - shift[1]) + " " + str(z * stepZ - shift[2]) + "\n")

def write_vert_2D(file, nVert, dimensions, shift = [0, 0, 0]):
    stepX = (dimensions["x"]/(nVert["x"] - 1.0))
    stepY = (dimensions["y"]/(nVert["y"] - 1.0))
    for y in range(nVert["y"]):
        for x in range(nVert["x"]):
            file.write(str(x * stepX - shift[0]) + " " + str(y * stepY - shift[1]) + " 0.0 \n")

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

def write_cells_2D(file, nVert):
    for i in range(int(nVert["x"]) - 1):
        for j in range(int(nVert["y"]) - 1):
            file.write(
                "4 " +
                str(i + j * nVert["x"]) + " " +
                str(i+1 + j * nVert["x"]) + " " +
                str(i+1 + (j + 1) * nVert["x"]) + " " +
                str(i + (j + 1) * nVert["x"]) + "\n"
            )

def generate_vtk_block(name, dimensions, number_of_vertices):
    with open(name, "w") as file:

        file.write(
            """# vtk DataFile Version 2.0
3D test mesh
ASCII
DATASET UNSTRUCTURED_GRID
POINTS """)

        file.write(str(number_of_vertices["x"] *
                    number_of_vertices["y"] *
                    number_of_vertices["z"]) + " double\n")
        write_vert(file, number_of_vertices, dimensions)

        file.write("\nCELLS ")
        nCells = (number_of_vertices["x"] - 1) * (number_of_vertices["y"] - 1) * (number_of_vertices["z"] - 1)
        file.write(str(nCells) + " " + str(nCells * 9) + "\n")

        write_cells(file, number_of_vertices)

        file.write("\nCELL_TYPES " + str(nCells) + "\n")
        for _ in range(nCells):
            file.write("12\n")

def generate_vtk_rectangle(name, dimensions, number_of_vertices):
    with open(name, "w") as file:

        file.write(
            """# vtk DataFile Version 2.0
2D test mesh
ASCII
DATASET UNSTRUCTURED_GRID
POINTS """)

        file.write(str(number_of_vertices["x"] *
                    number_of_vertices["y"]) + " double\n")
        write_vert_2D(file, number_of_vertices, dimensions)

        file.write("\nCELLS ")
        nCells = (number_of_vertices["x"] - 1) * (number_of_vertices["y"] - 1)
        file.write(str(nCells) + " " + str(nCells * 5) + "\n")

        write_cells_2D(file, number_of_vertices)

        file.write("\nCELL_TYPES " + str(nCells) + "\n")
        for _ in range(nCells):
            file.write("9\n")

def main():
    """
    Main function of the simple mesh generator. 
    It produces two files with meshes. One 2D and one 3D.
    Modify the parameters of the functions to produce finer meshes.
    """
    dimensions_2D = {"x":0.5, "y":1}
    num_cells_2D = 20
    number_of_vertices_2D = {"x":num_cells_2D + 1, "y":num_cells_2D * 2 + 1}
    generate_vtk_rectangle(f"mesh2D.vtk", dimensions_2D, number_of_vertices_2D)

    dimensions_3D = {"x":0.5, "y":0.5, "z":1}
    num_cells_3D = 9
    number_of_vertices_3D = {"x":num_cells_3D + 1, "y":num_cells_3D + 1, "z":num_cells_3D * 2 + 1}
    generate_vtk_block(f"mesh3D.vtk", dimensions_3D, number_of_vertices_3D)

if __name__ == "__main__":
    main()