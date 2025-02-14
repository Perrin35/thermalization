OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.390653) q[0];
sx q[0];
rz(-0.57354623) q[0];
sx q[0];
rz(-0.030800495) q[0];
rz(-2.7544694) q[1];
sx q[1];
rz(-0.20063278) q[1];
sx q[1];
rz(0.90955847) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5873588) q[0];
sx q[0];
rz(-0.92992126) q[0];
sx q[0];
rz(-0.72312042) q[0];
rz(2.4417905) q[2];
sx q[2];
rz(-2.0462817) q[2];
sx q[2];
rz(0.83719992) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8768898) q[1];
sx q[1];
rz(-2.8976421) q[1];
sx q[1];
rz(0.96744858) q[1];
rz(-pi) q[2];
rz(-3.0372871) q[3];
sx q[3];
rz(-1.3005195) q[3];
sx q[3];
rz(-0.53846151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51547852) q[2];
sx q[2];
rz(-1.4592183) q[2];
sx q[2];
rz(-2.1558351) q[2];
rz(-1.23729) q[3];
sx q[3];
rz(-2.3111549) q[3];
sx q[3];
rz(1.5498836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3378231) q[0];
sx q[0];
rz(-1.060312) q[0];
sx q[0];
rz(2.7983303) q[0];
rz(-2.901851) q[1];
sx q[1];
rz(-2.7570351) q[1];
sx q[1];
rz(-3.1125617) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.507258) q[0];
sx q[0];
rz(-1.0936948) q[0];
sx q[0];
rz(1.3647337) q[0];
x q[1];
rz(3.0986039) q[2];
sx q[2];
rz(-2.413297) q[2];
sx q[2];
rz(0.16253072) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5933725) q[1];
sx q[1];
rz(-0.97785972) q[1];
sx q[1];
rz(1.5895459) q[1];
rz(-1.9594479) q[3];
sx q[3];
rz(-0.93224469) q[3];
sx q[3];
rz(-0.030872542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3451738) q[2];
sx q[2];
rz(-2.5721305) q[2];
sx q[2];
rz(0.84863895) q[2];
rz(0.38550115) q[3];
sx q[3];
rz(-1.4717088) q[3];
sx q[3];
rz(-2.8482385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8141396) q[0];
sx q[0];
rz(-1.3804945) q[0];
sx q[0];
rz(-3.0019548) q[0];
rz(-2.8756554) q[1];
sx q[1];
rz(-1.9640924) q[1];
sx q[1];
rz(1.0136484) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0767727) q[0];
sx q[0];
rz(-1.0014604) q[0];
sx q[0];
rz(0.4607309) q[0];
rz(-pi) q[1];
rz(1.1590454) q[2];
sx q[2];
rz(-1.9728185) q[2];
sx q[2];
rz(-2.2184203) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2340548) q[1];
sx q[1];
rz(-2.6352596) q[1];
sx q[1];
rz(2.7975026) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6740213) q[3];
sx q[3];
rz(-1.6397392) q[3];
sx q[3];
rz(-0.063223039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.93918148) q[2];
sx q[2];
rz(-0.95042578) q[2];
sx q[2];
rz(0.44516304) q[2];
rz(-2.4116481) q[3];
sx q[3];
rz(-1.4533307) q[3];
sx q[3];
rz(1.2030407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79778033) q[0];
sx q[0];
rz(-0.36425632) q[0];
sx q[0];
rz(-1.9601747) q[0];
rz(-1.6254609) q[1];
sx q[1];
rz(-1.6226035) q[1];
sx q[1];
rz(2.6536062) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51633126) q[0];
sx q[0];
rz(-0.30540403) q[0];
sx q[0];
rz(-2.5918079) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71731848) q[2];
sx q[2];
rz(-1.0861168) q[2];
sx q[2];
rz(-0.48004638) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9727461) q[1];
sx q[1];
rz(-0.98727814) q[1];
sx q[1];
rz(2.8295687) q[1];
rz(0.41792385) q[3];
sx q[3];
rz(-1.2439787) q[3];
sx q[3];
rz(0.55909294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2130012) q[2];
sx q[2];
rz(-1.4139621) q[2];
sx q[2];
rz(-0.50814) q[2];
rz(0.43934923) q[3];
sx q[3];
rz(-2.5033247) q[3];
sx q[3];
rz(0.58258575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64684922) q[0];
sx q[0];
rz(-1.2127533) q[0];
sx q[0];
rz(-0.61432046) q[0];
rz(2.1454504) q[1];
sx q[1];
rz(-2.3042945) q[1];
sx q[1];
rz(0.070929758) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3230266) q[0];
sx q[0];
rz(-2.1686318) q[0];
sx q[0];
rz(-2.7672309) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0148914) q[2];
sx q[2];
rz(-1.8307865) q[2];
sx q[2];
rz(-0.6636493) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.48551895) q[1];
sx q[1];
rz(-0.89886256) q[1];
sx q[1];
rz(0.80826064) q[1];
rz(-2.2493241) q[3];
sx q[3];
rz(-2.1315977) q[3];
sx q[3];
rz(2.4099007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7190242) q[2];
sx q[2];
rz(-1.4113798) q[2];
sx q[2];
rz(-0.076449797) q[2];
rz(-0.092175305) q[3];
sx q[3];
rz(-1.1276827) q[3];
sx q[3];
rz(0.60747147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80716625) q[0];
sx q[0];
rz(-0.18944117) q[0];
sx q[0];
rz(0.46519753) q[0];
rz(-0.094296433) q[1];
sx q[1];
rz(-1.0153208) q[1];
sx q[1];
rz(-1.2634855) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7955728) q[0];
sx q[0];
rz(-1.4515522) q[0];
sx q[0];
rz(-1.9247965) q[0];
rz(-pi) q[1];
rz(-1.5689079) q[2];
sx q[2];
rz(-0.93853653) q[2];
sx q[2];
rz(2.8304932) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0537529) q[1];
sx q[1];
rz(-2.1422763) q[1];
sx q[1];
rz(-2.0896795) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1912548) q[3];
sx q[3];
rz(-1.1357682) q[3];
sx q[3];
rz(-1.3635927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5491526) q[2];
sx q[2];
rz(-2.0024039) q[2];
sx q[2];
rz(0.69063866) q[2];
rz(-1.6820827) q[3];
sx q[3];
rz(-3.1091318) q[3];
sx q[3];
rz(2.8800817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9948147) q[0];
sx q[0];
rz(-1.0228461) q[0];
sx q[0];
rz(-2.6652375) q[0];
rz(-1.9626544) q[1];
sx q[1];
rz(-0.6911239) q[1];
sx q[1];
rz(1.3271416) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19366143) q[0];
sx q[0];
rz(-1.1817939) q[0];
sx q[0];
rz(-0.48402985) q[0];
rz(-0.05169241) q[2];
sx q[2];
rz(-1.0982656) q[2];
sx q[2];
rz(2.56531) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33846291) q[1];
sx q[1];
rz(-1.5005365) q[1];
sx q[1];
rz(0.52609013) q[1];
rz(-pi) q[2];
rz(1.5487212) q[3];
sx q[3];
rz(-0.95609236) q[3];
sx q[3];
rz(-2.2382977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2771161) q[2];
sx q[2];
rz(-1.6552552) q[2];
sx q[2];
rz(-2.3790996) q[2];
rz(-1.3660733) q[3];
sx q[3];
rz(-1.6906747) q[3];
sx q[3];
rz(2.8097927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45258006) q[0];
sx q[0];
rz(-0.068050139) q[0];
sx q[0];
rz(-2.3633603) q[0];
rz(-1.7315841) q[1];
sx q[1];
rz(-0.86763132) q[1];
sx q[1];
rz(-0.57077879) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1335412) q[0];
sx q[0];
rz(-0.92317177) q[0];
sx q[0];
rz(-0.53093853) q[0];
rz(-pi) q[1];
rz(-1.8536751) q[2];
sx q[2];
rz(-1.6215417) q[2];
sx q[2];
rz(0.03374781) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.82680271) q[1];
sx q[1];
rz(-0.52779454) q[1];
sx q[1];
rz(-1.7262682) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.93511411) q[3];
sx q[3];
rz(-1.550727) q[3];
sx q[3];
rz(0.76616132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1125696) q[2];
sx q[2];
rz(-1.3444291) q[2];
sx q[2];
rz(2.6018108) q[2];
rz(-2.6953186) q[3];
sx q[3];
rz(-0.73791426) q[3];
sx q[3];
rz(2.3257183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096263252) q[0];
sx q[0];
rz(-2.6427866) q[0];
sx q[0];
rz(-0.86736429) q[0];
rz(1.7034886) q[1];
sx q[1];
rz(-2.4080364) q[1];
sx q[1];
rz(-0.47384706) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4854914) q[0];
sx q[0];
rz(-2.7655619) q[0];
sx q[0];
rz(-0.72304) q[0];
x q[1];
rz(2.176193) q[2];
sx q[2];
rz(-1.2393481) q[2];
sx q[2];
rz(1.6939577) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6240166) q[1];
sx q[1];
rz(-1.4291457) q[1];
sx q[1];
rz(-2.4017548) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2603277) q[3];
sx q[3];
rz(-1.9816958) q[3];
sx q[3];
rz(-1.9104065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.69502407) q[2];
sx q[2];
rz(-0.5646255) q[2];
sx q[2];
rz(2.0107644) q[2];
rz(-1.1787777) q[3];
sx q[3];
rz(-1.3935573) q[3];
sx q[3];
rz(-2.668837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9216565) q[0];
sx q[0];
rz(-0.49880767) q[0];
sx q[0];
rz(-2.1529799) q[0];
rz(1.0722718) q[1];
sx q[1];
rz(-1.6055454) q[1];
sx q[1];
rz(1.4354717) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0611183) q[0];
sx q[0];
rz(-2.1466564) q[0];
sx q[0];
rz(-0.030035069) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9133116) q[2];
sx q[2];
rz(-1.5633226) q[2];
sx q[2];
rz(2.8808336) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2769673) q[1];
sx q[1];
rz(-1.5709779) q[1];
sx q[1];
rz(-1.6305883) q[1];
rz(-pi) q[2];
rz(-1.2352474) q[3];
sx q[3];
rz(-2.23508) q[3];
sx q[3];
rz(1.2227525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.84206218) q[2];
sx q[2];
rz(-2.5008423) q[2];
sx q[2];
rz(-0.93471849) q[2];
rz(-2.0958021) q[3];
sx q[3];
rz(-2.9691634) q[3];
sx q[3];
rz(0.24833965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0228148) q[0];
sx q[0];
rz(-1.3640484) q[0];
sx q[0];
rz(1.7400297) q[0];
rz(-0.10764311) q[1];
sx q[1];
rz(-1.5168774) q[1];
sx q[1];
rz(0.37227896) q[1];
rz(2.8714928) q[2];
sx q[2];
rz(-2.3348689) q[2];
sx q[2];
rz(0.058090799) q[2];
rz(2.1984062) q[3];
sx q[3];
rz(-1.9780157) q[3];
sx q[3];
rz(-2.1383022) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
