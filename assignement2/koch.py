import numpy as np
import matplotlib.pyplot as plt

def koch_generator(start: np.array, end: np.array) -> np.array:
    start = np.array(start)
    end = np.array(end)
    
    dx, dy = (end - start) / 4  # steps

    p1 = start
    p2 = p1 + np.array([dx, dy]) 
    p3 = p2 + np.array([-dy, dx]) # raise
    p4 = p3 + np.array([dx, dy])  
    p5 = p4 + np.array([dy, -dx])   # lower
    p6 = p5 + np.array([dy, -dx]) # useful later when higher levels
    p7 = p6 + np.array([dx, dy])
    p8 = p7 + np.array([-dy, dx]) # raise
    p9 = end

    return np.array([p1, p2, p3, p4, p5, p6, p7, p8, p9])

def generate_koch_curve(start: np.array, end: np.array, level: int) -> np.array:
    if level == 0:
        return np.array([start, end])    
    
    points = koch_generator(start, end)
    result = []

    for i in range(len(points) - 1):
        result.extend(generate_koch_curve(points[i], points[i + 1], level - 1))

    result.append(end) #add end since len(points)-1
    
    return np.array(result)

def visualize_koch(points: np.array):
    fig, ax = plt.subplots()
    plt.plot(points[:, 0], points[:, 1])
    plt.show()

# Example usage
start_point = np.array([0, 0])
end_point = np.array([1, 0])
level = 2  # level of koch curve

koch_points = generate_koch_curve(start_point, end_point, level)
visualize_koch(koch_points)

