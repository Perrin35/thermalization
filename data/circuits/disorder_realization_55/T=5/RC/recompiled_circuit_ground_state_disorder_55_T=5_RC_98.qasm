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
rz(-0.70631385) q[1];
sx q[1];
rz(-0.7847473) q[1];
sx q[1];
rz(0.027332505) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3008049) q[0];
sx q[0];
rz(-2.7484924) q[0];
sx q[0];
rz(-0.53728763) q[0];
rz(-0.66157056) q[2];
sx q[2];
rz(-2.9230766) q[2];
sx q[2];
rz(-2.0486346) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8107218) q[1];
sx q[1];
rz(-1.6020944) q[1];
sx q[1];
rz(-0.16783952) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23689278) q[3];
sx q[3];
rz(-1.4290571) q[3];
sx q[3];
rz(0.6976034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0722384) q[2];
sx q[2];
rz(-1.8202929) q[2];
sx q[2];
rz(-1.5551152) q[2];
rz(2.41411) q[3];
sx q[3];
rz(-2.1860217) q[3];
sx q[3];
rz(-0.47310841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6721866) q[0];
sx q[0];
rz(-1.3453901) q[0];
sx q[0];
rz(2.192002) q[0];
rz(0.32497111) q[1];
sx q[1];
rz(-1.3409216) q[1];
sx q[1];
rz(-0.51345888) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0810222) q[0];
sx q[0];
rz(-0.96310421) q[0];
sx q[0];
rz(-2.6430703) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9302692) q[2];
sx q[2];
rz(-0.7604593) q[2];
sx q[2];
rz(2.9414842) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9565485) q[1];
sx q[1];
rz(-2.0868851) q[1];
sx q[1];
rz(-0.98770492) q[1];
x q[2];
rz(-0.98987119) q[3];
sx q[3];
rz(-1.9912915) q[3];
sx q[3];
rz(-2.6591501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.010633858) q[2];
sx q[2];
rz(-1.3291239) q[2];
sx q[2];
rz(2.0531674) q[2];
rz(2.4736577) q[3];
sx q[3];
rz(-2.9974944) q[3];
sx q[3];
rz(1.7264504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82069355) q[0];
sx q[0];
rz(-0.90105337) q[0];
sx q[0];
rz(-1.2574842) q[0];
rz(0.084818689) q[1];
sx q[1];
rz(-0.091400472) q[1];
sx q[1];
rz(-2.4841097) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5705868) q[0];
sx q[0];
rz(-2.0181542) q[0];
sx q[0];
rz(1.4259095) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.74983238) q[2];
sx q[2];
rz(-2.3978491) q[2];
sx q[2];
rz(0.86803555) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8264933) q[1];
sx q[1];
rz(-1.7414474) q[1];
sx q[1];
rz(-0.48433372) q[1];
rz(2.231953) q[3];
sx q[3];
rz(-1.3245021) q[3];
sx q[3];
rz(1.9734188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.377044) q[2];
sx q[2];
rz(-2.720764) q[2];
sx q[2];
rz(-1.2505038) q[2];
rz(-1.4416384) q[3];
sx q[3];
rz(-1.7504033) q[3];
sx q[3];
rz(-2.3120248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68469754) q[0];
sx q[0];
rz(-2.1853515) q[0];
sx q[0];
rz(-2.5392505) q[0];
rz(-2.1777228) q[1];
sx q[1];
rz(-0.78097051) q[1];
sx q[1];
rz(2.9197555) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5103127) q[0];
sx q[0];
rz(-1.5820192) q[0];
sx q[0];
rz(1.5785286) q[0];
rz(-pi) q[1];
rz(1.9811976) q[2];
sx q[2];
rz(-2.5735435) q[2];
sx q[2];
rz(-0.33524738) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.027699359) q[1];
sx q[1];
rz(-2.4287819) q[1];
sx q[1];
rz(-0.63473397) q[1];
rz(-pi) q[2];
rz(-0.1476361) q[3];
sx q[3];
rz(-2.6872928) q[3];
sx q[3];
rz(2.4678178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.43231371) q[2];
sx q[2];
rz(-1.1608492) q[2];
sx q[2];
rz(3.0529037) q[2];
rz(2.6426219) q[3];
sx q[3];
rz(-1.5191398) q[3];
sx q[3];
rz(0.078908198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0093507) q[0];
sx q[0];
rz(-3.0920588) q[0];
sx q[0];
rz(-0.90231878) q[0];
rz(0.29014507) q[1];
sx q[1];
rz(-0.90466181) q[1];
sx q[1];
rz(-1.9754999) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13577382) q[0];
sx q[0];
rz(-2.5516833) q[0];
sx q[0];
rz(-2.8481872) q[0];
rz(1.3275215) q[2];
sx q[2];
rz(-1.3009047) q[2];
sx q[2];
rz(2.6950767) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.20953791) q[1];
sx q[1];
rz(-1.2880858) q[1];
sx q[1];
rz(-1.9268031) q[1];
rz(-3.0405255) q[3];
sx q[3];
rz(-2.8887081) q[3];
sx q[3];
rz(-0.28836461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0151998) q[2];
sx q[2];
rz(-2.3455399) q[2];
sx q[2];
rz(-1.8533206) q[2];
rz(-2.1081693) q[3];
sx q[3];
rz(-2.0234334) q[3];
sx q[3];
rz(1.4748632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7159336) q[0];
sx q[0];
rz(-0.10361828) q[0];
sx q[0];
rz(-0.15956751) q[0];
rz(0.60699925) q[1];
sx q[1];
rz(-1.9885352) q[1];
sx q[1];
rz(-1.0608231) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55218766) q[0];
sx q[0];
rz(-0.97597835) q[0];
sx q[0];
rz(-2.908559) q[0];
rz(0.69691806) q[2];
sx q[2];
rz(-2.7978511) q[2];
sx q[2];
rz(2.4375107) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.80172864) q[1];
sx q[1];
rz(-0.9802981) q[1];
sx q[1];
rz(-2.9782359) q[1];
rz(-2.9922375) q[3];
sx q[3];
rz(-2.5575482) q[3];
sx q[3];
rz(0.67549878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1242421) q[2];
sx q[2];
rz(-0.81581798) q[2];
sx q[2];
rz(-2.7537277) q[2];
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
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2514728) q[0];
sx q[0];
rz(-2.1595182) q[0];
sx q[0];
rz(0.89333308) q[0];
rz(0.54126254) q[1];
sx q[1];
rz(-0.90507871) q[1];
sx q[1];
rz(-1.062324) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73570255) q[0];
sx q[0];
rz(-2.2828809) q[0];
sx q[0];
rz(1.2232365) q[0];
x q[1];
rz(-1.5028033) q[2];
sx q[2];
rz(-0.74882245) q[2];
sx q[2];
rz(0.48922005) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9683065) q[1];
sx q[1];
rz(-1.5943375) q[1];
sx q[1];
rz(2.5858155) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3618062) q[3];
sx q[3];
rz(-0.7721484) q[3];
sx q[3];
rz(2.0603555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.58747753) q[2];
sx q[2];
rz(-0.97439659) q[2];
sx q[2];
rz(-1.1678196) q[2];
rz(2.3217412) q[3];
sx q[3];
rz(-0.76698118) q[3];
sx q[3];
rz(-2.4988373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5646566) q[0];
sx q[0];
rz(-2.7695203) q[0];
sx q[0];
rz(1.8336953) q[0];
rz(-1.4564184) q[1];
sx q[1];
rz(-1.5142199) q[1];
sx q[1];
rz(3.0110722) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.099206) q[0];
sx q[0];
rz(-1.644028) q[0];
sx q[0];
rz(1.8797148) q[0];
rz(-0.79733915) q[2];
sx q[2];
rz(-1.9249462) q[2];
sx q[2];
rz(0.43978271) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4314506) q[1];
sx q[1];
rz(-1.5317066) q[1];
sx q[1];
rz(2.4907656) q[1];
rz(2.2240488) q[3];
sx q[3];
rz(-2.1751479) q[3];
sx q[3];
rz(2.0487599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.1859583) q[2];
sx q[2];
rz(-2.3027577) q[2];
sx q[2];
rz(1.7308509) q[2];
rz(0.12061067) q[3];
sx q[3];
rz(-1.6552304) q[3];
sx q[3];
rz(1.5131148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8677583) q[0];
sx q[0];
rz(-0.61536106) q[0];
sx q[0];
rz(0.14061418) q[0];
rz(-2.4070542) q[1];
sx q[1];
rz(-1.4334375) q[1];
sx q[1];
rz(-0.75070423) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1575506) q[0];
sx q[0];
rz(-0.37985248) q[0];
sx q[0];
rz(-0.67870875) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2241988) q[2];
sx q[2];
rz(-1.3470413) q[2];
sx q[2];
rz(-0.72347298) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.80440758) q[1];
sx q[1];
rz(-0.18964968) q[1];
sx q[1];
rz(-1.2387161) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4846913) q[3];
sx q[3];
rz(-2.4331581) q[3];
sx q[3];
rz(1.3902067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.004403) q[2];
sx q[2];
rz(-2.4418094) q[2];
sx q[2];
rz(-2.3889551) q[2];
rz(-2.802012) q[3];
sx q[3];
rz(-1.3656253) q[3];
sx q[3];
rz(-3.1201709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.705377) q[0];
sx q[0];
rz(-0.74420539) q[0];
sx q[0];
rz(2.5035653) q[0];
rz(-0.72884196) q[1];
sx q[1];
rz(-1.3325656) q[1];
sx q[1];
rz(0.80150882) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7201231) q[0];
sx q[0];
rz(-2.3178708) q[0];
sx q[0];
rz(-2.864564) q[0];
x q[1];
rz(-2.9000145) q[2];
sx q[2];
rz(-0.8525341) q[2];
sx q[2];
rz(0.11213839) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3190075) q[1];
sx q[1];
rz(-0.73627824) q[1];
sx q[1];
rz(-1.5126615) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9800945) q[3];
sx q[3];
rz(-2.209389) q[3];
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
rz(-2.6778383) q[3];
sx q[3];
rz(-1.2254865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7365702) q[0];
sx q[0];
rz(-2.6800192) q[0];
sx q[0];
rz(-0.17929684) q[0];
rz(0.90514056) q[1];
sx q[1];
rz(-2.1642579) q[1];
sx q[1];
rz(-0.49107818) q[1];
rz(1.1726562) q[2];
sx q[2];
rz(-2.3742994) q[2];
sx q[2];
rz(-0.91579043) q[2];
rz(-0.2855339) q[3];
sx q[3];
rz(-2.0750506) q[3];
sx q[3];
rz(-1.688523) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
