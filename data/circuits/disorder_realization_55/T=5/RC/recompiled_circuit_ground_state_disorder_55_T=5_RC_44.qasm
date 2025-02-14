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
rz(2.4352788) q[1];
sx q[1];
rz(-2.3568454) q[1];
sx q[1];
rz(3.1142601) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77307149) q[0];
sx q[0];
rz(-1.3734682) q[0];
sx q[0];
rz(-2.7993566) q[0];
rz(-pi) q[1];
rz(-0.66157056) q[2];
sx q[2];
rz(-2.9230766) q[2];
sx q[2];
rz(-2.0486346) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8107218) q[1];
sx q[1];
rz(-1.6020944) q[1];
sx q[1];
rz(0.16783952) q[1];
rz(-pi) q[2];
rz(-0.23689278) q[3];
sx q[3];
rz(-1.7125355) q[3];
sx q[3];
rz(0.6976034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.469406) q[0];
sx q[0];
rz(-1.3453901) q[0];
sx q[0];
rz(2.192002) q[0];
rz(2.8166215) q[1];
sx q[1];
rz(-1.3409216) q[1];
sx q[1];
rz(0.51345888) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0810222) q[0];
sx q[0];
rz(-2.1784884) q[0];
sx q[0];
rz(-2.6430703) q[0];
x q[1];
rz(-1.2113234) q[2];
sx q[2];
rz(-0.7604593) q[2];
sx q[2];
rz(-2.9414842) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1850441) q[1];
sx q[1];
rz(-1.0547075) q[1];
sx q[1];
rz(-2.1538877) q[1];
x q[2];
rz(-0.88708892) q[3];
sx q[3];
rz(-2.4389431) q[3];
sx q[3];
rz(1.4969373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.010633858) q[2];
sx q[2];
rz(-1.3291239) q[2];
sx q[2];
rz(1.0884253) q[2];
rz(-0.66793495) q[3];
sx q[3];
rz(-2.9974944) q[3];
sx q[3];
rz(1.7264504) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3208991) q[0];
sx q[0];
rz(-0.90105337) q[0];
sx q[0];
rz(-1.2574842) q[0];
rz(-3.056774) q[1];
sx q[1];
rz(-0.091400472) q[1];
sx q[1];
rz(0.65748293) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8958998) q[0];
sx q[0];
rz(-0.4687216) q[0];
sx q[0];
rz(-0.29229852) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5490513) q[2];
sx q[2];
rz(-1.0912025) q[2];
sx q[2];
rz(1.3035989) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.056524371) q[1];
sx q[1];
rz(-2.6303362) q[1];
sx q[1];
rz(-2.7871217) q[1];
x q[2];
rz(-2.231953) q[3];
sx q[3];
rz(-1.3245021) q[3];
sx q[3];
rz(-1.9734188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.377044) q[2];
sx q[2];
rz(-2.720764) q[2];
sx q[2];
rz(-1.8910889) q[2];
rz(-1.6999543) q[3];
sx q[3];
rz(-1.7504033) q[3];
sx q[3];
rz(-0.82956782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4568951) q[0];
sx q[0];
rz(-0.95624113) q[0];
sx q[0];
rz(-2.5392505) q[0];
rz(2.1777228) q[1];
sx q[1];
rz(-2.3606221) q[1];
sx q[1];
rz(2.9197555) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9070061) q[0];
sx q[0];
rz(-3.127964) q[0];
sx q[0];
rz(-0.6032633) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8922563) q[2];
sx q[2];
rz(-1.0549003) q[2];
sx q[2];
rz(0.14125401) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.1138933) q[1];
sx q[1];
rz(-2.4287819) q[1];
sx q[1];
rz(0.63473397) q[1];
x q[2];
rz(1.6425124) q[3];
sx q[3];
rz(-1.1218024) q[3];
sx q[3];
rz(-2.303799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.43231371) q[2];
sx q[2];
rz(-1.9807434) q[2];
sx q[2];
rz(-3.0529037) q[2];
rz(0.49897075) q[3];
sx q[3];
rz(-1.5191398) q[3];
sx q[3];
rz(3.0626845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13224193) q[0];
sx q[0];
rz(-0.049533822) q[0];
sx q[0];
rz(2.2392739) q[0];
rz(-0.29014507) q[1];
sx q[1];
rz(-2.2369308) q[1];
sx q[1];
rz(-1.9754999) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48448823) q[0];
sx q[0];
rz(-1.00923) q[0];
sx q[0];
rz(1.7620371) q[0];
rz(2.8639216) q[2];
sx q[2];
rz(-1.8050965) q[2];
sx q[2];
rz(1.1903583) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2578967) q[1];
sx q[1];
rz(-1.2295112) q[1];
sx q[1];
rz(-2.8410556) q[1];
rz(3.0405255) q[3];
sx q[3];
rz(-2.8887081) q[3];
sx q[3];
rz(-2.853228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0151998) q[2];
sx q[2];
rz(-2.3455399) q[2];
sx q[2];
rz(-1.2882721) q[2];
rz(-1.0334233) q[3];
sx q[3];
rz(-1.1181592) q[3];
sx q[3];
rz(-1.6667295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7159336) q[0];
sx q[0];
rz(-3.0379744) q[0];
sx q[0];
rz(-0.15956751) q[0];
rz(2.5345934) q[1];
sx q[1];
rz(-1.9885352) q[1];
sx q[1];
rz(1.0608231) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2552065) q[0];
sx q[0];
rz(-1.3783424) q[0];
sx q[0];
rz(-0.9631971) q[0];
rz(-2.873704) q[2];
sx q[2];
rz(-1.3527591) q[2];
sx q[2];
rz(1.6074558) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80172864) q[1];
sx q[1];
rz(-0.9802981) q[1];
sx q[1];
rz(0.16335674) q[1];
rz(-pi) q[2];
rz(-1.4727598) q[3];
sx q[3];
rz(-0.99409249) q[3];
sx q[3];
rz(0.85395472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0173505) q[2];
sx q[2];
rz(-2.3257747) q[2];
sx q[2];
rz(-2.7537277) q[2];
rz(-1.3346437) q[3];
sx q[3];
rz(-0.67238656) q[3];
sx q[3];
rz(-1.0268802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2514728) q[0];
sx q[0];
rz(-2.1595182) q[0];
sx q[0];
rz(-2.2482596) q[0];
rz(2.6003301) q[1];
sx q[1];
rz(-0.90507871) q[1];
sx q[1];
rz(1.062324) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5389299) q[0];
sx q[0];
rz(-1.8315803) q[0];
sx q[0];
rz(-0.74270026) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.063060305) q[2];
sx q[2];
rz(-0.82411924) q[2];
sx q[2];
rz(-2.7450739) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7294638) q[1];
sx q[1];
rz(-2.1264014) q[1];
sx q[1];
rz(-1.5985065) q[1];
rz(-2.331953) q[3];
sx q[3];
rz(-1.4255377) q[3];
sx q[3];
rz(0.64034772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.58747753) q[2];
sx q[2];
rz(-2.1671961) q[2];
sx q[2];
rz(1.973773) q[2];
rz(-0.81985146) q[3];
sx q[3];
rz(-0.76698118) q[3];
sx q[3];
rz(0.64275536) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5646566) q[0];
sx q[0];
rz(-0.37207237) q[0];
sx q[0];
rz(-1.3078974) q[0];
rz(-1.6851743) q[1];
sx q[1];
rz(-1.5142199) q[1];
sx q[1];
rz(-3.0110722) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6365287) q[0];
sx q[0];
rz(-1.8788596) q[0];
sx q[0];
rz(-3.0647361) q[0];
rz(-pi) q[1];
rz(-1.0840262) q[2];
sx q[2];
rz(-0.83510885) q[2];
sx q[2];
rz(1.4722919) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.91192818) q[1];
sx q[1];
rz(-2.4897631) q[1];
sx q[1];
rz(-3.0771281) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71575266) q[3];
sx q[3];
rz(-2.0945315) q[3];
sx q[3];
rz(0.88818404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.1859583) q[2];
sx q[2];
rz(-0.83883494) q[2];
sx q[2];
rz(1.4107417) q[2];
rz(0.12061067) q[3];
sx q[3];
rz(-1.6552304) q[3];
sx q[3];
rz(-1.6284778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8677583) q[0];
sx q[0];
rz(-2.5262316) q[0];
sx q[0];
rz(-3.0009785) q[0];
rz(0.73453844) q[1];
sx q[1];
rz(-1.4334375) q[1];
sx q[1];
rz(2.3908884) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44245369) q[0];
sx q[0];
rz(-1.8635731) q[0];
sx q[0];
rz(1.325216) q[0];
rz(-pi) q[1];
rz(-1.9173938) q[2];
sx q[2];
rz(-1.7945514) q[2];
sx q[2];
rz(0.72347298) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0929326) q[1];
sx q[1];
rz(-1.6322929) q[1];
sx q[1];
rz(-1.7503121) q[1];
rz(-pi) q[2];
rz(-2.2773975) q[3];
sx q[3];
rz(-1.5148124) q[3];
sx q[3];
rz(0.2460446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.004403) q[2];
sx q[2];
rz(-2.4418094) q[2];
sx q[2];
rz(-0.75263754) q[2];
rz(2.802012) q[3];
sx q[3];
rz(-1.3656253) q[3];
sx q[3];
rz(3.1201709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4362157) q[0];
sx q[0];
rz(-2.3973873) q[0];
sx q[0];
rz(0.63802737) q[0];
rz(2.4127507) q[1];
sx q[1];
rz(-1.3325656) q[1];
sx q[1];
rz(0.80150882) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8014098) q[0];
sx q[0];
rz(-1.7728285) q[0];
sx q[0];
rz(-0.8043181) q[0];
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
rz(pi/2) q[0];
rz(0.82258517) q[1];
sx q[1];
rz(-0.73627824) q[1];
sx q[1];
rz(1.5126615) q[1];
rz(-pi) q[2];
rz(1.3574886) q[3];
sx q[3];
rz(-2.4856728) q[3];
sx q[3];
rz(1.3154958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8375497) q[2];
sx q[2];
rz(-1.3648698) q[2];
sx q[2];
rz(-2.020828) q[2];
rz(2.4325727) q[3];
sx q[3];
rz(-0.46375436) q[3];
sx q[3];
rz(1.2254865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-2.7837743) q[2];
sx q[2];
rz(-2.2651548) q[2];
sx q[2];
rz(1.6969778) q[2];
rz(-2.0927248) q[3];
sx q[3];
rz(-1.3216139) q[3];
sx q[3];
rz(-0.25861964) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
