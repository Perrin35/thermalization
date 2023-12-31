OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5821563) q[0];
sx q[0];
rz(-0.59098935) q[0];
sx q[0];
rz(0.58340573) q[0];
rz(2.9572339) q[1];
sx q[1];
rz(-0.98362041) q[1];
sx q[1];
rz(-0.89259994) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28461449) q[0];
sx q[0];
rz(-1.0191139) q[0];
sx q[0];
rz(-2.7906228) q[0];
rz(-pi) q[1];
rz(-2.5298654) q[2];
sx q[2];
rz(-2.3731542) q[2];
sx q[2];
rz(-2.7035463) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.8436369) q[1];
sx q[1];
rz(-1.8814109) q[1];
sx q[1];
rz(-0.85308869) q[1];
rz(-2.3905972) q[3];
sx q[3];
rz(-1.5392116) q[3];
sx q[3];
rz(-2.3834474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9154174) q[2];
sx q[2];
rz(-1.5390652) q[2];
sx q[2];
rz(1.9809451) q[2];
rz(2.9246269) q[3];
sx q[3];
rz(-0.522885) q[3];
sx q[3];
rz(-1.0552361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0579257) q[0];
sx q[0];
rz(-2.1766429) q[0];
sx q[0];
rz(0.57587409) q[0];
rz(-1.8946164) q[1];
sx q[1];
rz(-1.2966825) q[1];
sx q[1];
rz(-1.974568) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.089433) q[0];
sx q[0];
rz(-0.26425996) q[0];
sx q[0];
rz(1.3454076) q[0];
x q[1];
rz(-1.9677656) q[2];
sx q[2];
rz(-2.5644828) q[2];
sx q[2];
rz(-1.4156262) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.963672) q[1];
sx q[1];
rz(-0.76621395) q[1];
sx q[1];
rz(-1.6597762) q[1];
rz(-pi) q[2];
rz(1.0807651) q[3];
sx q[3];
rz(-0.96066517) q[3];
sx q[3];
rz(2.3572363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1077659) q[2];
sx q[2];
rz(-1.9627389) q[2];
sx q[2];
rz(-2.9272184) q[2];
rz(-0.073444627) q[3];
sx q[3];
rz(-0.44973222) q[3];
sx q[3];
rz(-0.28081056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4784933) q[0];
sx q[0];
rz(-0.61613023) q[0];
sx q[0];
rz(1.9146772) q[0];
rz(0.40027174) q[1];
sx q[1];
rz(-1.8881256) q[1];
sx q[1];
rz(2.1267557) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0600216) q[0];
sx q[0];
rz(-2.2131753) q[0];
sx q[0];
rz(1.5006256) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8995908) q[2];
sx q[2];
rz(-1.9654462) q[2];
sx q[2];
rz(1.637527) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6940569) q[1];
sx q[1];
rz(-0.60255614) q[1];
sx q[1];
rz(-1.7130997) q[1];
rz(-pi) q[2];
rz(-1.8869927) q[3];
sx q[3];
rz(-0.71101515) q[3];
sx q[3];
rz(0.3871813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0456475) q[2];
sx q[2];
rz(-2.084705) q[2];
sx q[2];
rz(-2.2423559) q[2];
rz(-2.4441161) q[3];
sx q[3];
rz(-1.8557502) q[3];
sx q[3];
rz(-0.60788679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65524453) q[0];
sx q[0];
rz(-2.0589893) q[0];
sx q[0];
rz(-2.7096601) q[0];
rz(-2.5090384) q[1];
sx q[1];
rz(-2.7245941) q[1];
sx q[1];
rz(2.5057709) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5084002) q[0];
sx q[0];
rz(-2.0661372) q[0];
sx q[0];
rz(0.20423996) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7556778) q[2];
sx q[2];
rz(-0.73542483) q[2];
sx q[2];
rz(-2.3787969) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7446049) q[1];
sx q[1];
rz(-2.3481391) q[1];
sx q[1];
rz(-0.26054392) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5084247) q[3];
sx q[3];
rz(-2.79106) q[3];
sx q[3];
rz(-1.7115953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2137961) q[2];
sx q[2];
rz(-0.14363229) q[2];
sx q[2];
rz(2.4528743) q[2];
rz(0.33411807) q[3];
sx q[3];
rz(-1.170661) q[3];
sx q[3];
rz(-2.9978602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14389811) q[0];
sx q[0];
rz(-2.4425638) q[0];
sx q[0];
rz(2.0671663) q[0];
rz(0.74514666) q[1];
sx q[1];
rz(-1.6586168) q[1];
sx q[1];
rz(0.27854663) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51019788) q[0];
sx q[0];
rz(-2.3017831) q[0];
sx q[0];
rz(0.98548074) q[0];
rz(-pi) q[1];
rz(-1.6804382) q[2];
sx q[2];
rz(-2.1796558) q[2];
sx q[2];
rz(-1.7313752) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3814195) q[1];
sx q[1];
rz(-2.1790494) q[1];
sx q[1];
rz(-0.45254032) q[1];
rz(-2.1513125) q[3];
sx q[3];
rz(-2.2970082) q[3];
sx q[3];
rz(-0.57675225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4429861) q[2];
sx q[2];
rz(-1.7437982) q[2];
sx q[2];
rz(2.8473575) q[2];
rz(0.081929835) q[3];
sx q[3];
rz(-0.51920813) q[3];
sx q[3];
rz(-3.1053655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0329523) q[0];
sx q[0];
rz(-2.2886798) q[0];
sx q[0];
rz(-0.0090573514) q[0];
rz(-2.5065705) q[1];
sx q[1];
rz(-0.68990866) q[1];
sx q[1];
rz(-3.0335398) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6022588) q[0];
sx q[0];
rz(-2.5784011) q[0];
sx q[0];
rz(-0.62298933) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.33424218) q[2];
sx q[2];
rz(-0.40194449) q[2];
sx q[2];
rz(-2.3964756) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5080155) q[1];
sx q[1];
rz(-1.0883691) q[1];
sx q[1];
rz(-2.0111994) q[1];
rz(-1.990854) q[3];
sx q[3];
rz(-0.58745158) q[3];
sx q[3];
rz(-2.1849039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.99284995) q[2];
sx q[2];
rz(-1.381258) q[2];
sx q[2];
rz(-1.8704869) q[2];
rz(0.078401119) q[3];
sx q[3];
rz(-1.6631118) q[3];
sx q[3];
rz(1.1318077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25061297) q[0];
sx q[0];
rz(-1.5654634) q[0];
sx q[0];
rz(-0.65761956) q[0];
rz(-1.3972067) q[1];
sx q[1];
rz(-1.1746635) q[1];
sx q[1];
rz(-2.2479642) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052505715) q[0];
sx q[0];
rz(-2.1930709) q[0];
sx q[0];
rz(-1.1629348) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1003175) q[2];
sx q[2];
rz(-0.5404226) q[2];
sx q[2];
rz(-1.4735917) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.212008) q[1];
sx q[1];
rz(-1.3699023) q[1];
sx q[1];
rz(0.68540539) q[1];
x q[2];
rz(3.0393533) q[3];
sx q[3];
rz(-0.7623626) q[3];
sx q[3];
rz(1.8765212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.39067337) q[2];
sx q[2];
rz(-0.94736391) q[2];
sx q[2];
rz(-2.612109) q[2];
rz(0.47618619) q[3];
sx q[3];
rz(-1.3323077) q[3];
sx q[3];
rz(0.34255323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3787518) q[0];
sx q[0];
rz(-1.6289926) q[0];
sx q[0];
rz(1.09028) q[0];
rz(0.11225637) q[1];
sx q[1];
rz(-1.1039762) q[1];
sx q[1];
rz(-1.1539248) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.740828) q[0];
sx q[0];
rz(-2.4908998) q[0];
sx q[0];
rz(-3.0853737) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23837337) q[2];
sx q[2];
rz(-2.5775238) q[2];
sx q[2];
rz(0.43018815) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6653319) q[1];
sx q[1];
rz(-1.6325103) q[1];
sx q[1];
rz(0.76121059) q[1];
rz(-pi) q[2];
rz(1.0941986) q[3];
sx q[3];
rz(-2.4123441) q[3];
sx q[3];
rz(0.13343982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5899137) q[2];
sx q[2];
rz(-2.7605197) q[2];
sx q[2];
rz(-0.38044688) q[2];
rz(1.1278661) q[3];
sx q[3];
rz(-1.7815855) q[3];
sx q[3];
rz(-1.5765223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8001051) q[0];
sx q[0];
rz(-1.8999758) q[0];
sx q[0];
rz(-2.5019116) q[0];
rz(1.2387964) q[1];
sx q[1];
rz(-1.7265373) q[1];
sx q[1];
rz(-1.9715086) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026924883) q[0];
sx q[0];
rz(-1.9244734) q[0];
sx q[0];
rz(-0.35004079) q[0];
rz(-pi) q[1];
rz(-0.77634546) q[2];
sx q[2];
rz(-1.7661957) q[2];
sx q[2];
rz(0.59567829) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.942109) q[1];
sx q[1];
rz(-1.0892727) q[1];
sx q[1];
rz(1.7755309) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89948489) q[3];
sx q[3];
rz(-1.6791108) q[3];
sx q[3];
rz(1.2182106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3141979) q[2];
sx q[2];
rz(-0.65076995) q[2];
sx q[2];
rz(-0.38273746) q[2];
rz(-0.9283723) q[3];
sx q[3];
rz(-1.174077) q[3];
sx q[3];
rz(-0.66463566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2255573) q[0];
sx q[0];
rz(-1.501843) q[0];
sx q[0];
rz(2.8826707) q[0];
rz(2.4312773) q[1];
sx q[1];
rz(-1.0537035) q[1];
sx q[1];
rz(-2.6616667) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5537162) q[0];
sx q[0];
rz(-1.237861) q[0];
sx q[0];
rz(-2.2399708) q[0];
rz(-pi) q[1];
rz(2.9051022) q[2];
sx q[2];
rz(-0.91954008) q[2];
sx q[2];
rz(2.1361534) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5727947) q[1];
sx q[1];
rz(-2.5017782) q[1];
sx q[1];
rz(0.6154284) q[1];
x q[2];
rz(0.5273401) q[3];
sx q[3];
rz(-2.1747327) q[3];
sx q[3];
rz(0.37249836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2991128) q[2];
sx q[2];
rz(-2.2128236) q[2];
sx q[2];
rz(-1.3170362) q[2];
rz(-1.8995829) q[3];
sx q[3];
rz(-0.95359355) q[3];
sx q[3];
rz(2.4035113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9823572) q[0];
sx q[0];
rz(-3.1095105) q[0];
sx q[0];
rz(1.6788917) q[0];
rz(-2.1622529) q[1];
sx q[1];
rz(-1.0995438) q[1];
sx q[1];
rz(-0.88811036) q[1];
rz(0.32817763) q[2];
sx q[2];
rz(-0.50445088) q[2];
sx q[2];
rz(2.9403461) q[2];
rz(-0.76673037) q[3];
sx q[3];
rz(-0.37692108) q[3];
sx q[3];
rz(2.2212096) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
