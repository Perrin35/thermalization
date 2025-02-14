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
rz(-2.7890132) q[0];
rz(-0.70631385) q[1];
sx q[1];
rz(-0.7847473) q[1];
sx q[1];
rz(0.027332505) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72799435) q[0];
sx q[0];
rz(-1.9061273) q[0];
sx q[0];
rz(-1.7799313) q[0];
rz(1.4352082) q[2];
sx q[2];
rz(-1.3989046) q[2];
sx q[2];
rz(0.41968771) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.057291659) q[1];
sx q[1];
rz(-2.9708869) q[1];
sx q[1];
rz(-0.18526669) q[1];
rz(-0.54630791) q[3];
sx q[3];
rz(-2.8662193) q[3];
sx q[3];
rz(-0.34378036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0693543) q[2];
sx q[2];
rz(-1.8202929) q[2];
sx q[2];
rz(1.5551152) q[2];
rz(2.41411) q[3];
sx q[3];
rz(-2.1860217) q[3];
sx q[3];
rz(2.6684842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
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
rz(-1.800671) q[1];
sx q[1];
rz(-0.51345888) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18841132) q[0];
sx q[0];
rz(-1.1674178) q[0];
sx q[0];
rz(0.90103407) q[0];
rz(-0.32294257) q[2];
sx q[2];
rz(-2.2720798) q[2];
sx q[2];
rz(-2.8632134) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4411575) q[1];
sx q[1];
rz(-2.0701906) q[1];
sx q[1];
rz(0.59696377) q[1];
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
rz(3.1309588) q[2];
sx q[2];
rz(-1.8124688) q[2];
sx q[2];
rz(2.0531674) q[2];
rz(2.4736577) q[3];
sx q[3];
rz(-0.1440983) q[3];
sx q[3];
rz(1.4151423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3208991) q[0];
sx q[0];
rz(-2.2405393) q[0];
sx q[0];
rz(-1.2574842) q[0];
rz(-0.084818689) q[1];
sx q[1];
rz(-0.091400472) q[1];
sx q[1];
rz(2.4841097) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5705868) q[0];
sx q[0];
rz(-2.0181542) q[0];
sx q[0];
rz(1.7156832) q[0];
x q[1];
rz(-0.59254139) q[2];
sx q[2];
rz(-2.0503902) q[2];
sx q[2];
rz(1.8379938) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8264933) q[1];
sx q[1];
rz(-1.7414474) q[1];
sx q[1];
rz(-2.6572589) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.182193) q[3];
sx q[3];
rz(-0.69903421) q[3];
sx q[3];
rz(-2.4352668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.377044) q[2];
sx q[2];
rz(-2.720764) q[2];
sx q[2];
rz(1.2505038) q[2];
rz(1.6999543) q[3];
sx q[3];
rz(-1.7504033) q[3];
sx q[3];
rz(-2.3120248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68469754) q[0];
sx q[0];
rz(-0.95624113) q[0];
sx q[0];
rz(2.5392505) q[0];
rz(2.1777228) q[1];
sx q[1];
rz(-0.78097051) q[1];
sx q[1];
rz(0.22183713) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6312799) q[0];
sx q[0];
rz(-1.5595734) q[0];
sx q[0];
rz(-1.5785286) q[0];
x q[1];
rz(-1.0413136) q[2];
sx q[2];
rz(-1.7871249) q[2];
sx q[2];
rz(-1.554503) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.79961046) q[1];
sx q[1];
rz(-2.1253808) q[1];
sx q[1];
rz(-1.0971402) q[1];
rz(1.6425124) q[3];
sx q[3];
rz(-1.1218024) q[3];
sx q[3];
rz(0.83779369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7092789) q[2];
sx q[2];
rz(-1.9807434) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0093507) q[0];
sx q[0];
rz(-3.0920588) q[0];
sx q[0];
rz(2.2392739) q[0];
rz(0.29014507) q[1];
sx q[1];
rz(-0.90466181) q[1];
sx q[1];
rz(-1.9754999) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0058188) q[0];
sx q[0];
rz(-2.5516833) q[0];
sx q[0];
rz(-2.8481872) q[0];
rz(0.71641123) q[2];
sx q[2];
rz(-0.36135095) q[2];
sx q[2];
rz(-2.8385069) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2578967) q[1];
sx q[1];
rz(-1.2295112) q[1];
sx q[1];
rz(-2.8410556) q[1];
rz(-pi) q[2];
rz(-1.5447292) q[3];
sx q[3];
rz(-1.3192303) q[3];
sx q[3];
rz(2.9575916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0151998) q[2];
sx q[2];
rz(-0.79605278) q[2];
sx q[2];
rz(-1.2882721) q[2];
rz(1.0334233) q[3];
sx q[3];
rz(-1.1181592) q[3];
sx q[3];
rz(1.6667295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.425659) q[0];
sx q[0];
rz(-3.0379744) q[0];
sx q[0];
rz(2.9820251) q[0];
rz(0.60699925) q[1];
sx q[1];
rz(-1.9885352) q[1];
sx q[1];
rz(2.0807696) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2552065) q[0];
sx q[0];
rz(-1.3783424) q[0];
sx q[0];
rz(0.9631971) q[0];
rz(-pi) q[1];
rz(-0.26788864) q[2];
sx q[2];
rz(-1.3527591) q[2];
sx q[2];
rz(-1.6074558) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.80172864) q[1];
sx q[1];
rz(-0.9802981) q[1];
sx q[1];
rz(-2.9782359) q[1];
rz(-0.57890602) q[3];
sx q[3];
rz(-1.6529377) q[3];
sx q[3];
rz(-2.3711747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0173505) q[2];
sx q[2];
rz(-2.3257747) q[2];
sx q[2];
rz(-2.7537277) q[2];
rz(-1.806949) q[3];
sx q[3];
rz(-2.4692061) q[3];
sx q[3];
rz(-1.0268802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2514728) q[0];
sx q[0];
rz(-0.9820745) q[0];
sx q[0];
rz(-2.2482596) q[0];
rz(0.54126254) q[1];
sx q[1];
rz(-2.2365139) q[1];
sx q[1];
rz(1.062324) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8996595) q[0];
sx q[0];
rz(-0.77882732) q[0];
sx q[0];
rz(0.3758442) q[0];
x q[1];
rz(-2.3184653) q[2];
sx q[2];
rz(-1.6170653) q[2];
sx q[2];
rz(-1.1314162) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7294638) q[1];
sx q[1];
rz(-1.0151912) q[1];
sx q[1];
rz(1.5985065) q[1];
x q[2];
rz(2.9422308) q[3];
sx q[3];
rz(-2.3219675) q[3];
sx q[3];
rz(2.3481751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5541151) q[2];
sx q[2];
rz(-2.1671961) q[2];
sx q[2];
rz(1.1678196) q[2];
rz(-0.81985146) q[3];
sx q[3];
rz(-2.3746115) q[3];
sx q[3];
rz(-0.64275536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.576936) q[0];
sx q[0];
rz(-0.37207237) q[0];
sx q[0];
rz(1.3078974) q[0];
rz(1.4564184) q[1];
sx q[1];
rz(-1.5142199) q[1];
sx q[1];
rz(-3.0110722) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042386656) q[0];
sx q[0];
rz(-1.644028) q[0];
sx q[0];
rz(-1.2618778) q[0];
rz(-pi) q[1];
rz(2.0575664) q[2];
sx q[2];
rz(-2.3064838) q[2];
sx q[2];
rz(1.6693008) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.91192818) q[1];
sx q[1];
rz(-2.4897631) q[1];
sx q[1];
rz(3.0771281) q[1];
rz(-pi) q[2];
rz(2.42584) q[3];
sx q[3];
rz(-1.0470611) q[3];
sx q[3];
rz(2.2534086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9556344) q[2];
sx q[2];
rz(-0.83883494) q[2];
sx q[2];
rz(-1.4107417) q[2];
rz(-0.12061067) q[3];
sx q[3];
rz(-1.6552304) q[3];
sx q[3];
rz(1.6284778) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2738344) q[0];
sx q[0];
rz(-2.5262316) q[0];
sx q[0];
rz(-0.14061418) q[0];
rz(2.4070542) q[1];
sx q[1];
rz(-1.7081552) q[1];
sx q[1];
rz(2.3908884) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44245369) q[0];
sx q[0];
rz(-1.2780196) q[0];
sx q[0];
rz(-1.8163766) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23739135) q[2];
sx q[2];
rz(-1.9084045) q[2];
sx q[2];
rz(0.92729502) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3371851) q[1];
sx q[1];
rz(-2.951943) q[1];
sx q[1];
rz(-1.9028765) q[1];
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
rz(-2.1371896) q[2];
sx q[2];
rz(-2.4418094) q[2];
sx q[2];
rz(-2.3889551) q[2];
rz(-2.802012) q[3];
sx q[3];
rz(-1.7759674) q[3];
sx q[3];
rz(-0.021421758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4362157) q[0];
sx q[0];
rz(-2.3973873) q[0];
sx q[0];
rz(-2.5035653) q[0];
rz(-2.4127507) q[1];
sx q[1];
rz(-1.3325656) q[1];
sx q[1];
rz(-0.80150882) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4214695) q[0];
sx q[0];
rz(-0.82372181) q[0];
sx q[0];
rz(-0.27702866) q[0];
rz(-pi) q[1];
rz(-2.9000145) q[2];
sx q[2];
rz(-0.8525341) q[2];
sx q[2];
rz(0.11213839) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.90097016) q[1];
sx q[1];
rz(-0.83604807) q[1];
sx q[1];
rz(3.088984) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.784104) q[3];
sx q[3];
rz(-2.4856728) q[3];
sx q[3];
rz(-1.8260969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8375497) q[2];
sx q[2];
rz(-1.7767228) q[2];
sx q[2];
rz(-2.020828) q[2];
rz(2.4325727) q[3];
sx q[3];
rz(-2.6778383) q[3];
sx q[3];
rz(-1.2254865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(2.7837743) q[2];
sx q[2];
rz(-0.87643788) q[2];
sx q[2];
rz(-1.4446148) q[2];
rz(-2.8560588) q[3];
sx q[3];
rz(-1.0665421) q[3];
sx q[3];
rz(1.4530696) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
