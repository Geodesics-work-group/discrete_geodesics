# discrete_geodesics
A python script to calculate gedesics on different surfaces.

## Getting started

The script uses the libraries [NumPy](http://www.numpy.org/) and [Plotly](https://plot.ly/python/) for the calculation and graphic representation of surfaces.

For information about the installation of the libraries visit the following pages:

* [NumPy](https://www.scipy.org/install.html)
* [Plotly](https://plot.ly/python/getting-started/)

## Changing parameters

```python
h = ? # size of the intervals of t

distance = ? # total t

u = [?] # initial x position

v = [?] # initial y position

p = [?]  # initial x velocity

q = [?]  # initial y velocity
```

To change any of the values of the parametes of the Initial Value Problem, as well as the step taken and the total distance travelled, replace the "?" symbols with the desired values.

Note that in the sphere, the initial velocity vector (p,q) is normalized, so as to set with d the length of the arc generated by the algorithm.

