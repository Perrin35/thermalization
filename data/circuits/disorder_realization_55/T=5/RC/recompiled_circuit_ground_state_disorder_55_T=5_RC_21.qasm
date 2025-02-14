OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.9094698) q[0];
sx q[0];
rz(-2.1535518) q[0];
sx q[0];
rz(2.7890132) q[0];
rz(2.4352788) q[1];
sx q[1];
rz(-2.3568454) q[1];
sx q[1];
rz(3.1142601) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4135983) q[0];
sx q[0];
rz(-1.2354654) q[0];
sx q[0];
rz(-1.3616614) q[0];
rz(-pi) q[1];
rz(0.17345239) q[2];
sx q[2];
rz(-1.4372184) q[2];
sx q[2];
rz(2.0138149) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.24522745) q[1];
sx q[1];
rz(-1.7385529) q[1];
sx q[1];
rz(-1.6025402) q[1];
rz(-pi) q[2];
rz(0.54630791) q[3];
sx q[3];
rz(-0.27537333) q[3];
sx q[3];
rz(2.7978123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0693543) q[2];
sx q[2];
rz(-1.3212997) q[2];
sx q[2];
rz(1.5864774) q[2];
rz(2.41411) q[3];
sx q[3];
rz(-0.95557094) q[3];
sx q[3];
rz(0.47310841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6721866) q[0];
sx q[0];
rz(-1.7962026) q[0];
sx q[0];
rz(2.192002) q[0];
rz(-2.8166215) q[1];
sx q[1];
rz(-1.3409216) q[1];
sx q[1];
rz(-0.51345888) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0605704) q[0];
sx q[0];
rz(-0.96310421) q[0];
sx q[0];
rz(-2.6430703) q[0];
rz(1.2113234) q[2];
sx q[2];
rz(-0.7604593) q[2];
sx q[2];
rz(-0.20010848) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1850441) q[1];
sx q[1];
rz(-2.0868851) q[1];
sx q[1];
rz(-2.1538877) q[1];
x q[2];
rz(-2.2545037) q[3];
sx q[3];
rz(-0.70264951) q[3];
sx q[3];
rz(1.4969373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.010633858) q[2];
sx q[2];
rz(-1.8124688) q[2];
sx q[2];
rz(1.0884253) q[2];
rz(-0.66793495) q[3];
sx q[3];
rz(-0.1440983) q[3];
sx q[3];
rz(-1.7264504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82069355) q[0];
sx q[0];
rz(-0.90105337) q[0];
sx q[0];
rz(-1.8841085) q[0];
rz(-3.056774) q[1];
sx q[1];
rz(-0.091400472) q[1];
sx q[1];
rz(0.65748293) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5705868) q[0];
sx q[0];
rz(-1.1234385) q[0];
sx q[0];
rz(-1.7156832) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59254139) q[2];
sx q[2];
rz(-2.0503902) q[2];
sx q[2];
rz(-1.3035989) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7967826) q[1];
sx q[1];
rz(-1.0940897) q[1];
sx q[1];
rz(1.3784798) q[1];
rz(-pi) q[2];
rz(-0.90963962) q[3];
sx q[3];
rz(-1.8170905) q[3];
sx q[3];
rz(-1.9734188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76454863) q[2];
sx q[2];
rz(-2.720764) q[2];
sx q[2];
rz(-1.8910889) q[2];
rz(1.4416384) q[3];
sx q[3];
rz(-1.3911894) q[3];
sx q[3];
rz(0.82956782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4568951) q[0];
sx q[0];
rz(-2.1853515) q[0];
sx q[0];
rz(-0.60234219) q[0];
rz(-0.96386987) q[1];
sx q[1];
rz(-0.78097051) q[1];
sx q[1];
rz(-2.9197555) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.202163) q[0];
sx q[0];
rz(-1.5630645) q[0];
sx q[0];
rz(-3.1303694) q[0];
x q[1];
rz(1.1603951) q[2];
sx q[2];
rz(-0.56804915) q[2];
sx q[2];
rz(-0.33524738) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1138933) q[1];
sx q[1];
rz(-0.71281071) q[1];
sx q[1];
rz(0.63473397) q[1];
rz(-pi) q[2];
rz(0.45000123) q[3];
sx q[3];
rz(-1.5061989) q[3];
sx q[3];
rz(0.76417506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7092789) q[2];
sx q[2];
rz(-1.1608492) q[2];
sx q[2];
rz(-0.088689001) q[2];
rz(2.6426219) q[3];
sx q[3];
rz(-1.5191398) q[3];
sx q[3];
rz(0.078908198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0093507) q[0];
sx q[0];
rz(-3.0920588) q[0];
sx q[0];
rz(-0.90231878) q[0];
rz(-0.29014507) q[1];
sx q[1];
rz(-0.90466181) q[1];
sx q[1];
rz(-1.1660928) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0058188) q[0];
sx q[0];
rz(-2.5516833) q[0];
sx q[0];
rz(2.8481872) q[0];
x q[1];
rz(0.27767105) q[2];
sx q[2];
rz(-1.3364961) q[2];
sx q[2];
rz(-1.9512343) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0048827) q[1];
sx q[1];
rz(-0.45082475) q[1];
sx q[1];
rz(-2.2656127) q[1];
rz(1.5968634) q[3];
sx q[3];
rz(-1.8223624) q[3];
sx q[3];
rz(-2.9575916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0151998) q[2];
sx q[2];
rz(-0.79605278) q[2];
sx q[2];
rz(1.2882721) q[2];
rz(-1.0334233) q[3];
sx q[3];
rz(-1.1181592) q[3];
sx q[3];
rz(-1.6667295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.425659) q[0];
sx q[0];
rz(-0.10361828) q[0];
sx q[0];
rz(-0.15956751) q[0];
rz(-2.5345934) q[1];
sx q[1];
rz(-1.1530575) q[1];
sx q[1];
rz(-2.0807696) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2552065) q[0];
sx q[0];
rz(-1.7632503) q[0];
sx q[0];
rz(0.9631971) q[0];
x q[1];
rz(0.69691806) q[2];
sx q[2];
rz(-0.3437416) q[2];
sx q[2];
rz(-2.4375107) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.80172864) q[1];
sx q[1];
rz(-0.9802981) q[1];
sx q[1];
rz(-0.16335674) q[1];
rz(-pi) q[2];
rz(-1.4727598) q[3];
sx q[3];
rz(-2.1475002) q[3];
sx q[3];
rz(-0.85395472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1242421) q[2];
sx q[2];
rz(-0.81581798) q[2];
sx q[2];
rz(-2.7537277) q[2];
rz(-1.806949) q[3];
sx q[3];
rz(-2.4692061) q[3];
sx q[3];
rz(2.1147125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89011985) q[0];
sx q[0];
rz(-0.9820745) q[0];
sx q[0];
rz(2.2482596) q[0];
rz(-2.6003301) q[1];
sx q[1];
rz(-0.90507871) q[1];
sx q[1];
rz(2.0792686) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73570255) q[0];
sx q[0];
rz(-2.2828809) q[0];
sx q[0];
rz(1.9183561) q[0];
x q[1];
rz(-0.063060305) q[2];
sx q[2];
rz(-0.82411924) q[2];
sx q[2];
rz(0.39651878) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7819643) q[1];
sx q[1];
rz(-0.55622314) q[1];
sx q[1];
rz(0.044597711) q[1];
x q[2];
rz(-1.3618062) q[3];
sx q[3];
rz(-0.7721484) q[3];
sx q[3];
rz(1.0812372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.58747753) q[2];
sx q[2];
rz(-2.1671961) q[2];
sx q[2];
rz(1.973773) q[2];
rz(0.81985146) q[3];
sx q[3];
rz(-0.76698118) q[3];
sx q[3];
rz(2.4988373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.576936) q[0];
sx q[0];
rz(-2.7695203) q[0];
sx q[0];
rz(-1.8336953) q[0];
rz(-1.4564184) q[1];
sx q[1];
rz(-1.6273727) q[1];
sx q[1];
rz(-3.0110722) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7537788) q[0];
sx q[0];
rz(-2.8243833) q[0];
sx q[0];
rz(-1.3340201) q[0];
rz(1.0840262) q[2];
sx q[2];
rz(-2.3064838) q[2];
sx q[2];
rz(-1.6693008) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3106892) q[1];
sx q[1];
rz(-0.92055087) q[1];
sx q[1];
rz(1.6199153) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4198615) q[3];
sx q[3];
rz(-0.85866195) q[3];
sx q[3];
rz(-0.16068383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.1859583) q[2];
sx q[2];
rz(-2.3027577) q[2];
sx q[2];
rz(-1.4107417) q[2];
rz(0.12061067) q[3];
sx q[3];
rz(-1.4863622) q[3];
sx q[3];
rz(1.6284778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8677583) q[0];
sx q[0];
rz(-2.5262316) q[0];
sx q[0];
rz(-0.14061418) q[0];
rz(-2.4070542) q[1];
sx q[1];
rz(-1.4334375) q[1];
sx q[1];
rz(-0.75070423) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0561309) q[0];
sx q[0];
rz(-1.805725) q[0];
sx q[0];
rz(-0.30130193) q[0];
x q[1];
rz(1.2241988) q[2];
sx q[2];
rz(-1.7945514) q[2];
sx q[2];
rz(-2.4181197) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0486601) q[1];
sx q[1];
rz(-1.6322929) q[1];
sx q[1];
rz(1.7503121) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4846913) q[3];
sx q[3];
rz(-0.70843452) q[3];
sx q[3];
rz(-1.3902067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.004403) q[2];
sx q[2];
rz(-2.4418094) q[2];
sx q[2];
rz(-2.3889551) q[2];
rz(2.802012) q[3];
sx q[3];
rz(-1.7759674) q[3];
sx q[3];
rz(0.021421758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.705377) q[0];
sx q[0];
rz(-0.74420539) q[0];
sx q[0];
rz(0.63802737) q[0];
rz(0.72884196) q[1];
sx q[1];
rz(-1.809027) q[1];
sx q[1];
rz(-2.3400838) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4214695) q[0];
sx q[0];
rz(-2.3178708) q[0];
sx q[0];
rz(-0.27702866) q[0];
x q[1];
rz(-0.24157816) q[2];
sx q[2];
rz(-2.2890586) q[2];
sx q[2];
rz(-3.0294543) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4364797) q[1];
sx q[1];
rz(-1.6098238) q[1];
sx q[1];
rz(-2.3062333) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9800945) q[3];
sx q[3];
rz(-0.93220369) q[3];
sx q[3];
rz(-2.0929008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8375497) q[2];
sx q[2];
rz(-1.7767228) q[2];
sx q[2];
rz(-2.020828) q[2];
rz(-0.70901999) q[3];
sx q[3];
rz(-0.46375436) q[3];
sx q[3];
rz(-1.9161061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40502248) q[0];
sx q[0];
rz(-0.46157349) q[0];
sx q[0];
rz(2.9622958) q[0];
rz(0.90514056) q[1];
sx q[1];
rz(-2.1642579) q[1];
sx q[1];
rz(-0.49107818) q[1];
rz(-2.7837743) q[2];
sx q[2];
rz(-2.2651548) q[2];
sx q[2];
rz(1.6969778) q[2];
rz(-1.0488679) q[3];
sx q[3];
rz(-1.8199788) q[3];
sx q[3];
rz(2.882973) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
