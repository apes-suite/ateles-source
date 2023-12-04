title: Heat Conduction in 3D

The 3D heat conduction is configured by choosing the equation name to be
"heat_3d".
It takes one parameter, the heat conductivity "k":

```lua
equation = {
  name = 'heat_3d',
  k = 1
}
```

The initial condition has to be provided for the temperature.
An example is provided with a [sinus distribution](sinus_temperature).

