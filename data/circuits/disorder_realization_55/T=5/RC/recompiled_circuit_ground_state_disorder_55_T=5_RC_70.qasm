OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2321229) q[0];
sx q[0];
rz(-0.98804086) q[0];
sx q[0];
rz(0.35257941) q[0];
rz(-0.70631385) q[1];
sx q[1];
rz(-3.92634) q[1];
sx q[1];
rz(9.3974455) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8407878) q[0];
sx q[0];
rz(-2.7484924) q[0];
sx q[0];
rz(-0.53728763) q[0];
rz(-pi) q[1];
rz(2.9681403) q[2];
sx q[2];
rz(-1.4372184) q[2];
sx q[2];
rz(-2.0138149) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3308709) q[1];
sx q[1];
rz(-1.5394982) q[1];
sx q[1];
rz(-2.9737531) q[1];
x q[2];
rz(-2.9046999) q[3];
sx q[3];
rz(-1.4290571) q[3];
sx q[3];
rz(-2.4439892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0693543) q[2];
sx q[2];
rz(-1.3212997) q[2];
sx q[2];
rz(-1.5551152) q[2];
rz(-0.72748264) q[3];
sx q[3];
rz(-0.95557094) q[3];
sx q[3];
rz(0.47310841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6721866) q[0];
sx q[0];
rz(-1.7962026) q[0];
sx q[0];
rz(2.192002) q[0];
rz(-0.32497111) q[1];
sx q[1];
rz(-1.3409216) q[1];
sx q[1];
rz(0.51345888) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8425861) q[0];
sx q[0];
rz(-2.3761099) q[0];
sx q[0];
rz(0.96852901) q[0];
rz(-pi) q[1];
rz(-2.2983488) q[2];
sx q[2];
rz(-1.3258889) q[2];
sx q[2];
rz(1.6365192) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.70043515) q[1];
sx q[1];
rz(-1.0714021) q[1];
sx q[1];
rz(0.59696377) q[1];
rz(-pi) q[2];
x q[2];
rz(0.98987119) q[3];
sx q[3];
rz(-1.1503011) q[3];
sx q[3];
rz(-2.6591501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.1309588) q[2];
sx q[2];
rz(-1.8124688) q[2];
sx q[2];
rz(-2.0531674) q[2];
rz(2.4736577) q[3];
sx q[3];
rz(-0.1440983) q[3];
sx q[3];
rz(-1.7264504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3208991) q[0];
sx q[0];
rz(-0.90105337) q[0];
sx q[0];
rz(-1.8841085) q[0];
rz(0.084818689) q[1];
sx q[1];
rz(-0.091400472) q[1];
sx q[1];
rz(-2.4841097) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8958998) q[0];
sx q[0];
rz(-2.6728711) q[0];
sx q[0];
rz(2.8492941) q[0];
rz(-pi) q[1];
rz(2.3917603) q[2];
sx q[2];
rz(-2.3978491) q[2];
sx q[2];
rz(-2.2735571) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.3448101) q[1];
sx q[1];
rz(-2.0475029) q[1];
sx q[1];
rz(1.3784798) q[1];
rz(-pi) q[2];
rz(0.90963962) q[3];
sx q[3];
rz(-1.8170905) q[3];
sx q[3];
rz(1.9734188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.377044) q[2];
sx q[2];
rz(-2.720764) q[2];
sx q[2];
rz(1.8910889) q[2];
rz(-1.4416384) q[3];
sx q[3];
rz(-1.3911894) q[3];
sx q[3];
rz(-0.82956782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
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
rz(2.1777228) q[1];
sx q[1];
rz(-2.3606221) q[1];
sx q[1];
rz(2.9197555) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6312799) q[0];
sx q[0];
rz(-1.5820192) q[0];
sx q[0];
rz(-1.5785286) q[0];
rz(1.9811976) q[2];
sx q[2];
rz(-2.5735435) q[2];
sx q[2];
rz(2.8063453) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3419822) q[1];
sx q[1];
rz(-2.1253808) q[1];
sx q[1];
rz(-1.0971402) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1476361) q[3];
sx q[3];
rz(-2.6872928) q[3];
sx q[3];
rz(0.67377485) q[3];
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
rz(3.0529037) q[2];
rz(-2.6426219) q[3];
sx q[3];
rz(-1.5191398) q[3];
sx q[3];
rz(3.0626845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0093507) q[0];
sx q[0];
rz(-3.0920588) q[0];
sx q[0];
rz(-2.2392739) q[0];
rz(-0.29014507) q[1];
sx q[1];
rz(-0.90466181) q[1];
sx q[1];
rz(-1.1660928) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48448823) q[0];
sx q[0];
rz(-1.00923) q[0];
sx q[0];
rz(-1.3795555) q[0];
rz(-pi) q[1];
rz(0.27767105) q[2];
sx q[2];
rz(-1.3364961) q[2];
sx q[2];
rz(1.1903583) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.20953791) q[1];
sx q[1];
rz(-1.2880858) q[1];
sx q[1];
rz(-1.2147895) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8899447) q[3];
sx q[3];
rz(-1.5960428) q[3];
sx q[3];
rz(-1.3803052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.12639283) q[2];
sx q[2];
rz(-2.3455399) q[2];
sx q[2];
rz(1.2882721) q[2];
rz(2.1081693) q[3];
sx q[3];
rz(-1.1181592) q[3];
sx q[3];
rz(-1.6667295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.425659) q[0];
sx q[0];
rz(-0.10361828) q[0];
sx q[0];
rz(2.9820251) q[0];
rz(-2.5345934) q[1];
sx q[1];
rz(-1.9885352) q[1];
sx q[1];
rz(2.0807696) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2552065) q[0];
sx q[0];
rz(-1.7632503) q[0];
sx q[0];
rz(2.1783955) q[0];
rz(-pi) q[1];
rz(-2.873704) q[2];
sx q[2];
rz(-1.7888336) q[2];
sx q[2];
rz(-1.6074558) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2810106) q[1];
sx q[1];
rz(-1.4352903) q[1];
sx q[1];
rz(-2.1675088) q[1];
rz(-1.6688329) q[3];
sx q[3];
rz(-2.1475002) q[3];
sx q[3];
rz(0.85395472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1242421) q[2];
sx q[2];
rz(-2.3257747) q[2];
sx q[2];
rz(-0.38786495) q[2];
rz(1.806949) q[3];
sx q[3];
rz(-0.67238656) q[3];
sx q[3];
rz(2.1147125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89011985) q[0];
sx q[0];
rz(-2.1595182) q[0];
sx q[0];
rz(2.2482596) q[0];
rz(-0.54126254) q[1];
sx q[1];
rz(-0.90507871) q[1];
sx q[1];
rz(-2.0792686) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8996595) q[0];
sx q[0];
rz(-2.3627653) q[0];
sx q[0];
rz(-2.7657484) q[0];
rz(-pi) q[1];
rz(-1.5028033) q[2];
sx q[2];
rz(-2.3927702) q[2];
sx q[2];
rz(2.6523726) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4121288) q[1];
sx q[1];
rz(-1.0151912) q[1];
sx q[1];
rz(1.5985065) q[1];
rz(-1.7797864) q[3];
sx q[3];
rz(-0.7721484) q[3];
sx q[3];
rz(2.0603555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.58747753) q[2];
sx q[2];
rz(-2.1671961) q[2];
sx q[2];
rz(-1.973773) q[2];
rz(2.3217412) q[3];
sx q[3];
rz(-2.3746115) q[3];
sx q[3];
rz(2.4988373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5646566) q[0];
sx q[0];
rz(-0.37207237) q[0];
sx q[0];
rz(1.8336953) q[0];
rz(-1.4564184) q[1];
sx q[1];
rz(-1.6273727) q[1];
sx q[1];
rz(-3.0110722) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.505064) q[0];
sx q[0];
rz(-1.8788596) q[0];
sx q[0];
rz(0.076856514) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0840262) q[2];
sx q[2];
rz(-2.3064838) q[2];
sx q[2];
rz(-1.6693008) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3106892) q[1];
sx q[1];
rz(-0.92055087) q[1];
sx q[1];
rz(-1.6199153) q[1];
rz(-0.71575266) q[3];
sx q[3];
rz(-1.0470611) q[3];
sx q[3];
rz(-0.88818404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.1859583) q[2];
sx q[2];
rz(-2.3027577) q[2];
sx q[2];
rz(1.7308509) q[2];
rz(-0.12061067) q[3];
sx q[3];
rz(-1.6552304) q[3];
sx q[3];
rz(-1.5131148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8677583) q[0];
sx q[0];
rz(-0.61536106) q[0];
sx q[0];
rz(-0.14061418) q[0];
rz(0.73453844) q[1];
sx q[1];
rz(-1.4334375) q[1];
sx q[1];
rz(-0.75070423) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.984042) q[0];
sx q[0];
rz(-2.7617402) q[0];
sx q[0];
rz(-0.67870875) q[0];
x q[1];
rz(-2.9042013) q[2];
sx q[2];
rz(-1.9084045) q[2];
sx q[2];
rz(-0.92729502) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.80440758) q[1];
sx q[1];
rz(-2.951943) q[1];
sx q[1];
rz(1.9028765) q[1];
rz(0.86419515) q[3];
sx q[3];
rz(-1.5148124) q[3];
sx q[3];
rz(-2.8955481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1371896) q[2];
sx q[2];
rz(-0.69978324) q[2];
sx q[2];
rz(0.75263754) q[2];
rz(-0.33958069) q[3];
sx q[3];
rz(-1.3656253) q[3];
sx q[3];
rz(3.1201709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.705377) q[0];
sx q[0];
rz(-2.3973873) q[0];
sx q[0];
rz(0.63802737) q[0];
rz(2.4127507) q[1];
sx q[1];
rz(-1.809027) q[1];
sx q[1];
rz(-0.80150882) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0251533) q[0];
sx q[0];
rz(-0.7873767) q[0];
sx q[0];
rz(1.8579432) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24157816) q[2];
sx q[2];
rz(-2.2890586) q[2];
sx q[2];
rz(3.0294543) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2406225) q[1];
sx q[1];
rz(-2.3055446) q[1];
sx q[1];
rz(0.052608629) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3574886) q[3];
sx q[3];
rz(-0.65591988) q[3];
sx q[3];
rz(1.3154958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.30404299) q[2];
sx q[2];
rz(-1.3648698) q[2];
sx q[2];
rz(-2.020828) q[2];
rz(2.4325727) q[3];
sx q[3];
rz(-0.46375436) q[3];
sx q[3];
rz(-1.9161061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7365702) q[0];
sx q[0];
rz(-0.46157349) q[0];
sx q[0];
rz(2.9622958) q[0];
rz(0.90514056) q[1];
sx q[1];
rz(-2.1642579) q[1];
sx q[1];
rz(-0.49107818) q[1];
rz(-1.9689365) q[2];
sx q[2];
rz(-2.3742994) q[2];
sx q[2];
rz(-0.91579043) q[2];
rz(-1.0988476) q[3];
sx q[3];
rz(-2.5682378) q[3];
sx q[3];
rz(0.90709472) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
