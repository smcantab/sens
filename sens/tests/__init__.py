from _test_sa_sampler import *
from _test_weighted_pick import *
from _test_ns_lj import *
from _test_sens_lj import *
from _test_build_database import *
from _test_minima_searcher import *
from _test_sens_exact_lj import *
from _test_sa_sampler_ns import *
from _test_sens_exact_compute_energy_HSA import *

long_test = False
if long_test:
    from _test_ns_lj_long import *
    from _test_sens_lj_long import *

if __name__ == "__main__":
    unittest.main()  