title: Heat Conduction in 2D

The 2D heat conduction is configured by choosing the equation name to be
"heat_2d".
It takes one parameter, the heat conductivity "k":

```lua
equation = {
  name = 'heat_2d',
  k = 1
}
```

The initial condition has to be provided for the temperature.
An example is provided with a [sinus distribution](sinus_temperature).
