OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.20733362) q[0];
sx q[0];
rz(3.7319558) q[0];
sx q[0];
rz(9.0537602) q[0];
rz(-0.38129216) q[1];
sx q[1];
rz(-0.59950221) q[1];
sx q[1];
rz(1.376027) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57173079) q[0];
sx q[0];
rz(-0.96643448) q[0];
sx q[0];
rz(2.5992924) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4260169) q[2];
sx q[2];
rz(-1.1872429) q[2];
sx q[2];
rz(0.54563145) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0378117) q[1];
sx q[1];
rz(-1.722464) q[1];
sx q[1];
rz(-1.9019466) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0891221) q[3];
sx q[3];
rz(-2.3215508) q[3];
sx q[3];
rz(-0.47580556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0573037) q[2];
sx q[2];
rz(-0.40033445) q[2];
sx q[2];
rz(-2.1526745) q[2];
rz(2.3890498) q[3];
sx q[3];
rz(-1.9957333) q[3];
sx q[3];
rz(-0.83077103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9343524) q[0];
sx q[0];
rz(-3.0293284) q[0];
sx q[0];
rz(-1.9616615) q[0];
rz(2.143899) q[1];
sx q[1];
rz(-1.2832063) q[1];
sx q[1];
rz(-0.72431272) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2253101) q[0];
sx q[0];
rz(-0.75964576) q[0];
sx q[0];
rz(-1.1067953) q[0];
x q[1];
rz(-0.19947796) q[2];
sx q[2];
rz(-1.6530767) q[2];
sx q[2];
rz(0.85180887) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7608632) q[1];
sx q[1];
rz(-0.41178307) q[1];
sx q[1];
rz(2.1058583) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23726666) q[3];
sx q[3];
rz(-1.7627343) q[3];
sx q[3];
rz(-0.81258472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.90536845) q[2];
sx q[2];
rz(-1.3304109) q[2];
sx q[2];
rz(2.779707) q[2];
rz(3.0055255) q[3];
sx q[3];
rz(-0.55570221) q[3];
sx q[3];
rz(-3.0959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9996465) q[0];
sx q[0];
rz(-2.0479585) q[0];
sx q[0];
rz(-1.3954337) q[0];
rz(2.6793001) q[1];
sx q[1];
rz(-0.42458436) q[1];
sx q[1];
rz(-1.9225072) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8273979) q[0];
sx q[0];
rz(-2.6959246) q[0];
sx q[0];
rz(0.92339869) q[0];
x q[1];
rz(0.8692603) q[2];
sx q[2];
rz(-1.3420891) q[2];
sx q[2];
rz(1.9321835) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3238941) q[1];
sx q[1];
rz(-0.59191275) q[1];
sx q[1];
rz(-1.1281668) q[1];
rz(-pi) q[2];
rz(-0.015467042) q[3];
sx q[3];
rz(-1.8387715) q[3];
sx q[3];
rz(-1.0166575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.42157713) q[2];
sx q[2];
rz(-0.74024671) q[2];
sx q[2];
rz(1.6960309) q[2];
rz(-0.56882632) q[3];
sx q[3];
rz(-2.2918662) q[3];
sx q[3];
rz(3.0310757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9179984) q[0];
sx q[0];
rz(-2.693394) q[0];
sx q[0];
rz(-2.6233327) q[0];
rz(-0.7154243) q[1];
sx q[1];
rz(-2.0253069) q[1];
sx q[1];
rz(2.3148361) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2402156) q[0];
sx q[0];
rz(-0.6745406) q[0];
sx q[0];
rz(2.4504689) q[0];
rz(-2.3739359) q[2];
sx q[2];
rz(-1.4738238) q[2];
sx q[2];
rz(1.3354288) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0747448) q[1];
sx q[1];
rz(-1.2859935) q[1];
sx q[1];
rz(2.7436101) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4171962) q[3];
sx q[3];
rz(-1.3339692) q[3];
sx q[3];
rz(-2.608992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.15239079) q[2];
sx q[2];
rz(-0.19897904) q[2];
sx q[2];
rz(-1.7626804) q[2];
rz(0.072323024) q[3];
sx q[3];
rz(-2.3291589) q[3];
sx q[3];
rz(1.6453843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3146661) q[0];
sx q[0];
rz(-2.5214654) q[0];
sx q[0];
rz(-1.1258874) q[0];
rz(2.2391438) q[1];
sx q[1];
rz(-0.97389644) q[1];
sx q[1];
rz(-2.856423) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42442214) q[0];
sx q[0];
rz(-1.5283913) q[0];
sx q[0];
rz(-1.0958584) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9589013) q[2];
sx q[2];
rz(-0.94411196) q[2];
sx q[2];
rz(-2.2270122) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3030745) q[1];
sx q[1];
rz(-0.31305602) q[1];
sx q[1];
rz(-1.5694373) q[1];
x q[2];
rz(2.7399555) q[3];
sx q[3];
rz(-1.4816227) q[3];
sx q[3];
rz(-0.13084403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.29331648) q[2];
sx q[2];
rz(-2.6360376) q[2];
sx q[2];
rz(-0.53945333) q[2];
rz(-0.30682492) q[3];
sx q[3];
rz(-2.2570733) q[3];
sx q[3];
rz(-2.6873798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8834615) q[0];
sx q[0];
rz(-2.3911609) q[0];
sx q[0];
rz(0.15701292) q[0];
rz(-2.4482588) q[1];
sx q[1];
rz(-0.88070977) q[1];
sx q[1];
rz(1.3670115) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0317504) q[0];
sx q[0];
rz(-1.6037918) q[0];
sx q[0];
rz(-1.6389636) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1742646) q[2];
sx q[2];
rz(-0.8562932) q[2];
sx q[2];
rz(-0.29010233) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4821266) q[1];
sx q[1];
rz(-0.30059338) q[1];
sx q[1];
rz(-1.8468922) q[1];
rz(-2.8664687) q[3];
sx q[3];
rz(-0.36558662) q[3];
sx q[3];
rz(2.5610353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.29629016) q[2];
sx q[2];
rz(-2.2915816) q[2];
sx q[2];
rz(2.7381251) q[2];
rz(-0.48163313) q[3];
sx q[3];
rz(-1.0721595) q[3];
sx q[3];
rz(-2.6223555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8994609) q[0];
sx q[0];
rz(-0.88328981) q[0];
sx q[0];
rz(-0.8738628) q[0];
rz(-2.6938687) q[1];
sx q[1];
rz(-2.402585) q[1];
sx q[1];
rz(1.1706932) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3160352) q[0];
sx q[0];
rz(-0.023840126) q[0];
sx q[0];
rz(0.26160474) q[0];
rz(-pi) q[1];
rz(-3.0253719) q[2];
sx q[2];
rz(-1.5201609) q[2];
sx q[2];
rz(-1.2223787) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.891174) q[1];
sx q[1];
rz(-2.182057) q[1];
sx q[1];
rz(2.3120018) q[1];
rz(1.8740158) q[3];
sx q[3];
rz(-1.056864) q[3];
sx q[3];
rz(-2.5336888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0447023) q[2];
sx q[2];
rz(-0.56920749) q[2];
sx q[2];
rz(-2.5308385) q[2];
rz(-0.47510535) q[3];
sx q[3];
rz(-2.0510309) q[3];
sx q[3];
rz(-0.92774123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24191813) q[0];
sx q[0];
rz(-0.11706676) q[0];
sx q[0];
rz(-2.8444667) q[0];
rz(-1.7469453) q[1];
sx q[1];
rz(-1.9906094) q[1];
sx q[1];
rz(2.4954605) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4532115) q[0];
sx q[0];
rz(-1.4697945) q[0];
sx q[0];
rz(-1.7779011) q[0];
rz(-1.7392776) q[2];
sx q[2];
rz(-2.0373166) q[2];
sx q[2];
rz(-0.70665765) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1505193) q[1];
sx q[1];
rz(-1.0725613) q[1];
sx q[1];
rz(-1.9626161) q[1];
x q[2];
rz(2.2046702) q[3];
sx q[3];
rz(-2.2735032) q[3];
sx q[3];
rz(-0.015451775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0769161) q[2];
sx q[2];
rz(-0.95305324) q[2];
sx q[2];
rz(-2.6867552) q[2];
rz(-0.70139766) q[3];
sx q[3];
rz(-2.111179) q[3];
sx q[3];
rz(1.1340244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69944537) q[0];
sx q[0];
rz(-3*pi/16) q[0];
sx q[0];
rz(0.79750693) q[0];
rz(-0.51756716) q[1];
sx q[1];
rz(-2.3289754) q[1];
sx q[1];
rz(3.033175) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71589564) q[0];
sx q[0];
rz(-1.6285537) q[0];
sx q[0];
rz(-2.2349368) q[0];
rz(-1.8234532) q[2];
sx q[2];
rz(-2.9571819) q[2];
sx q[2];
rz(0.84469634) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0260967) q[1];
sx q[1];
rz(-2.8671088) q[1];
sx q[1];
rz(-2.2309169) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9547144) q[3];
sx q[3];
rz(-2.2664321) q[3];
sx q[3];
rz(-2.7587492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1291528) q[2];
sx q[2];
rz(-1.7947349) q[2];
sx q[2];
rz(2.8016395) q[2];
rz(0.41839504) q[3];
sx q[3];
rz(-2.5451626) q[3];
sx q[3];
rz(-2.4160014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5323935) q[0];
sx q[0];
rz(-2.7476855) q[0];
sx q[0];
rz(-0.6788196) q[0];
rz(0.36418307) q[1];
sx q[1];
rz(-1.697425) q[1];
sx q[1];
rz(3.0864339) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96345761) q[0];
sx q[0];
rz(-1.9473416) q[0];
sx q[0];
rz(-1.5012653) q[0];
rz(1.7807547) q[2];
sx q[2];
rz(-1.0210438) q[2];
sx q[2];
rz(-1.9090261) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0810869) q[1];
sx q[1];
rz(-2.2844237) q[1];
sx q[1];
rz(0.75018261) q[1];
x q[2];
rz(1.4961365) q[3];
sx q[3];
rz(-2.6844822) q[3];
sx q[3];
rz(-0.32170579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1577592) q[2];
sx q[2];
rz(-2.1201717) q[2];
sx q[2];
rz(0.49017635) q[2];
rz(3.0040719) q[3];
sx q[3];
rz(-2.0420045) q[3];
sx q[3];
rz(-0.93808758) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4162083) q[0];
sx q[0];
rz(-1.2513456) q[0];
sx q[0];
rz(-0.67847897) q[0];
rz(0.2086808) q[1];
sx q[1];
rz(-2.0188257) q[1];
sx q[1];
rz(1.6001736) q[1];
rz(0.31221496) q[2];
sx q[2];
rz(-1.4785462) q[2];
sx q[2];
rz(-1.2350456) q[2];
rz(-0.18110885) q[3];
sx q[3];
rz(-2.8977179) q[3];
sx q[3];
rz(0.44117622) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];