/**
# Spark Model

This module defines a spark (hot spot) which is useful for the
ingition of fuel-oxidizer mixtures in combustion simulations.
The idea, based on the implementation of [Cuoci et al. 2013](#cuoci2013numerical),
is to increase the temperature in a small spherical region of
domain defined by the user.
$$
T = T_s - (T_s - 300.)e^{-30(t-t_s)}
$$
where $T_s$ is the peak spark temperature, while $t_s$ is the
time of activation of the spark.
*/

#define SPARK

/**
## User data

We define a struct that gathers all the data required by the
spark model. The spark must be applied on a temperature field
which must be provided as *spark.T*.
*/

struct SparkModel {
  coord position;
  double time;
  double duration;
  double temperature;
  double diameter;
  double baseline;
  scalar T;
  bool linear;
  bool constant;
};

struct SparkModel spark = {
  .position = {0.,0.,0.}, // Coordinate of the ignition point
  .time = 0.,             // Time at which the ignition starts
  .duration = 0.,         // Duration of the spark
  .temperature = 0.,      // Maximum temperature of the spark
  .diameter = 0.,         // Diameter of the spark
  .baseline = 300.,       // Baseline/Minimum temperature value
  .linear = false,        // Linear temperature profile
  .constant = false,      // Constant temperature profile
};

/**
## Implementation

Main spark event, where the cells included by the
spark diameter are identified, and the temperature
increase is applied on these cells.
*/

#if dimension == 1
# define sparkd(x, y, R) (sq(R) - sq(x - spark.position.x))
#elif dimension == 2
# define sparkd(x, y, R) (sq(R) - sq(x - spark.position.x) - sq(y - spark.position.y))
#elif dimension == 3
# define sparkd(x, y, z, R) (sq(R) - sq(x - spark.position.x) - sq(y - spark.position.y) - sq(z - spark.position.z))
#endif

event set_spark (i++) {
#if TREE
  if (t <= spark.time) {
    extern int maxlevel;
    refine (sparkd(x, y, 0.5*spark.diameter) > 0. && level < maxlevel);
  }
#endif
  scalar spark_T = spark.T;
  if (t >= spark.time && t <= (spark.time + spark.duration)) {
    foreach() {
      if (sparkd(x,y,0.5*spark.diameter) > 0.) {
        if (spark.linear)
          spark_T[] = spark.baseline + (spark.temperature - spark.baseline)*(t - spark.time)/spark.duration;
        else if (spark.constant)
          spark_T[] = spark.temperature;
        else
          spark_T[] = spark.temperature -
            (spark.temperature - spark.baseline)*exp (-30.*(t - spark.time));
      }
    }
  }
}

/**
## References

~~~bib
@article{cuoci2013numerical,
  title={Numerical modeling of laminar flames with detailed kinetics based on the operator-splitting method},
  author={Cuoci, Alberto and Frassoldati, Alessio and Faravelli, Tiziano and Ranzi, Eliseo},
  journal={Energy \& fuels},
  volume={27},
  number={12},
  pages={7730--7753},
  year={2013},
  publisher={ACS Publications}
}
~~~
*/

