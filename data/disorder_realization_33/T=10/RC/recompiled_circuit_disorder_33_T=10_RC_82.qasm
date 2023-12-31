OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5053951) q[0];
sx q[0];
rz(-2.8656821) q[0];
sx q[0];
rz(-1.3077868) q[0];
rz(1.1360599) q[1];
sx q[1];
rz(-0.93568957) q[1];
sx q[1];
rz(1.5703262) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64496541) q[0];
sx q[0];
rz(-2.6500406) q[0];
sx q[0];
rz(-2.3636723) q[0];
x q[1];
rz(0.88408006) q[2];
sx q[2];
rz(-2.2097217) q[2];
sx q[2];
rz(-1.0654578) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.03304122) q[1];
sx q[1];
rz(-1.812495) q[1];
sx q[1];
rz(0.34717314) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9340431) q[3];
sx q[3];
rz(-1.3875291) q[3];
sx q[3];
rz(-0.43915425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2661665) q[2];
sx q[2];
rz(-2.8484919) q[2];
sx q[2];
rz(-1.1323294) q[2];
rz(1.6752361) q[3];
sx q[3];
rz(-1.3365859) q[3];
sx q[3];
rz(1.0124538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19673008) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(2.9557513) q[0];
rz(-0.56022412) q[1];
sx q[1];
rz(-1.2954243) q[1];
sx q[1];
rz(-2.9247608) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29138716) q[0];
sx q[0];
rz(-2.4225525) q[0];
sx q[0];
rz(-1.1262116) q[0];
rz(-pi) q[1];
rz(2.1255323) q[2];
sx q[2];
rz(-1.9811355) q[2];
sx q[2];
rz(1.0085269) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.617802) q[1];
sx q[1];
rz(-0.70980598) q[1];
sx q[1];
rz(1.8989423) q[1];
x q[2];
rz(0.84083765) q[3];
sx q[3];
rz(-2.3008122) q[3];
sx q[3];
rz(-0.18883146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8314787) q[2];
sx q[2];
rz(-2.3159413) q[2];
sx q[2];
rz(-1.2878093) q[2];
rz(2.3790322) q[3];
sx q[3];
rz(-1.972714) q[3];
sx q[3];
rz(0.30502239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4644311) q[0];
sx q[0];
rz(-2.7966249) q[0];
sx q[0];
rz(-2.537354) q[0];
rz(1.8151981) q[1];
sx q[1];
rz(-1.7809968) q[1];
sx q[1];
rz(-2.2089829) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7639097) q[0];
sx q[0];
rz(-1.8881067) q[0];
sx q[0];
rz(-3.0850287) q[0];
x q[1];
rz(2.9595397) q[2];
sx q[2];
rz(-1.3472054) q[2];
sx q[2];
rz(-1.455866) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.052913594) q[1];
sx q[1];
rz(-1.0148078) q[1];
sx q[1];
rz(1.4312137) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2266085) q[3];
sx q[3];
rz(-2.4768156) q[3];
sx q[3];
rz(1.8613601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9937667) q[2];
sx q[2];
rz(-2.0596762) q[2];
sx q[2];
rz(-2.0489342) q[2];
rz(2.5993733) q[3];
sx q[3];
rz(-1.0850302) q[3];
sx q[3];
rz(-0.96737635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7820691) q[0];
sx q[0];
rz(-0.096465915) q[0];
sx q[0];
rz(2.6413667) q[0];
rz(2.3362828) q[1];
sx q[1];
rz(-1.9814682) q[1];
sx q[1];
rz(1.6436228) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.362975) q[0];
sx q[0];
rz(-2.551429) q[0];
sx q[0];
rz(1.0233364) q[0];
x q[1];
rz(1.8565606) q[2];
sx q[2];
rz(-0.19837241) q[2];
sx q[2];
rz(0.49516585) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3264309) q[1];
sx q[1];
rz(-1.3109428) q[1];
sx q[1];
rz(1.3304779) q[1];
rz(-2.6763776) q[3];
sx q[3];
rz(-1.9808931) q[3];
sx q[3];
rz(-2.0801534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3952289) q[2];
sx q[2];
rz(-2.5791898) q[2];
sx q[2];
rz(2.4397819) q[2];
rz(0.83135215) q[3];
sx q[3];
rz(-2.1777007) q[3];
sx q[3];
rz(2.519616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2410626) q[0];
sx q[0];
rz(-2.5456972) q[0];
sx q[0];
rz(0.81533122) q[0];
rz(-1.5218081) q[1];
sx q[1];
rz(-2.3074469) q[1];
sx q[1];
rz(-2.0933847) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4002776) q[0];
sx q[0];
rz(-1.321723) q[0];
sx q[0];
rz(-1.8427909) q[0];
rz(-pi) q[1];
rz(1.5585209) q[2];
sx q[2];
rz(-2.2222812) q[2];
sx q[2];
rz(2.2001681) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8480307) q[1];
sx q[1];
rz(-1.3333496) q[1];
sx q[1];
rz(2.8192987) q[1];
rz(-2.1861595) q[3];
sx q[3];
rz(-1.2293929) q[3];
sx q[3];
rz(-2.0379025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.52577019) q[2];
sx q[2];
rz(-2.5746391) q[2];
sx q[2];
rz(2.0416416) q[2];
rz(-0.82529092) q[3];
sx q[3];
rz(-1.0422948) q[3];
sx q[3];
rz(2.2560789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(2.0681756) q[0];
sx q[0];
rz(-2.5475579) q[0];
sx q[0];
rz(2.2391879) q[0];
rz(2.1249318) q[1];
sx q[1];
rz(-2.0817751) q[1];
sx q[1];
rz(3.0117603) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5141402) q[0];
sx q[0];
rz(-2.1462626) q[0];
sx q[0];
rz(1.1260202) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1461357) q[2];
sx q[2];
rz(-1.8838922) q[2];
sx q[2];
rz(0.42195937) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9784769) q[1];
sx q[1];
rz(-0.69677959) q[1];
sx q[1];
rz(-2.5549868) q[1];
rz(-pi) q[2];
rz(2.8270709) q[3];
sx q[3];
rz(-0.57146996) q[3];
sx q[3];
rz(0.86626569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8292024) q[2];
sx q[2];
rz(-0.94909334) q[2];
sx q[2];
rz(0.20425805) q[2];
rz(-1.9355109) q[3];
sx q[3];
rz(-1.6198502) q[3];
sx q[3];
rz(0.23541418) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4181353) q[0];
sx q[0];
rz(-1.8122939) q[0];
sx q[0];
rz(-1.4468505) q[0];
rz(-1.8824668) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(0.68626219) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2080363) q[0];
sx q[0];
rz(-0.66970034) q[0];
sx q[0];
rz(2.3963388) q[0];
rz(-2.3115736) q[2];
sx q[2];
rz(-1.1168715) q[2];
sx q[2];
rz(1.9759076) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8441019) q[1];
sx q[1];
rz(-0.19980783) q[1];
sx q[1];
rz(-1.1872477) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8920184) q[3];
sx q[3];
rz(-1.1676844) q[3];
sx q[3];
rz(1.5954799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.69616047) q[2];
sx q[2];
rz(-1.3779209) q[2];
sx q[2];
rz(-3.1398204) q[2];
rz(-0.56162515) q[3];
sx q[3];
rz(-0.91149819) q[3];
sx q[3];
rz(1.5047489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.6034265) q[0];
sx q[0];
rz(-0.68646938) q[0];
sx q[0];
rz(-1.6954533) q[0];
rz(-0.7810477) q[1];
sx q[1];
rz(-1.8361517) q[1];
sx q[1];
rz(1.5015645) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38599309) q[0];
sx q[0];
rz(-2.0634683) q[0];
sx q[0];
rz(3.1194869) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.115032) q[2];
sx q[2];
rz(-1.6325258) q[2];
sx q[2];
rz(2.314687) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6540263) q[1];
sx q[1];
rz(-1.8974202) q[1];
sx q[1];
rz(1.4766272) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1498508) q[3];
sx q[3];
rz(-2.1175044) q[3];
sx q[3];
rz(-2.9255097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7897196) q[2];
sx q[2];
rz(-1.7799653) q[2];
sx q[2];
rz(-1.8224576) q[2];
rz(-1.2119279) q[3];
sx q[3];
rz(-1.8550248) q[3];
sx q[3];
rz(-2.8222728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8050352) q[0];
sx q[0];
rz(-0.55258495) q[0];
sx q[0];
rz(1.9375027) q[0];
rz(-0.38326344) q[1];
sx q[1];
rz(-2.6158694) q[1];
sx q[1];
rz(-2.7899172) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29031819) q[0];
sx q[0];
rz(-1.9418678) q[0];
sx q[0];
rz(-1.0822269) q[0];
rz(-0.69182379) q[2];
sx q[2];
rz(-1.3335388) q[2];
sx q[2];
rz(-1.1292063) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5896776) q[1];
sx q[1];
rz(-2.3304686) q[1];
sx q[1];
rz(-3.0685436) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6978108) q[3];
sx q[3];
rz(-0.1212596) q[3];
sx q[3];
rz(-1.6463943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7982771) q[2];
sx q[2];
rz(-1.1077935) q[2];
sx q[2];
rz(1.8593672) q[2];
rz(-1.4964237) q[3];
sx q[3];
rz(-1.5346425) q[3];
sx q[3];
rz(1.055868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4984109) q[0];
sx q[0];
rz(-1.2675985) q[0];
sx q[0];
rz(-0.19432755) q[0];
rz(-2.1037897) q[1];
sx q[1];
rz(-0.56832814) q[1];
sx q[1];
rz(2.1077572) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72458306) q[0];
sx q[0];
rz(-1.7010744) q[0];
sx q[0];
rz(-0.91086046) q[0];
rz(-pi) q[1];
rz(-2.2718272) q[2];
sx q[2];
rz(-1.1681721) q[2];
sx q[2];
rz(1.0588156) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.918805) q[1];
sx q[1];
rz(-0.58090392) q[1];
sx q[1];
rz(1.3851628) q[1];
rz(-pi) q[2];
rz(-1.3528353) q[3];
sx q[3];
rz(-2.2721604) q[3];
sx q[3];
rz(-0.71818128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0795435) q[2];
sx q[2];
rz(-2.1958308) q[2];
sx q[2];
rz(-0.6357843) q[2];
rz(2.87129) q[3];
sx q[3];
rz(-0.79939866) q[3];
sx q[3];
rz(1.5283782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6939659) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(1.4355961) q[1];
sx q[1];
rz(-1.5789079) q[1];
sx q[1];
rz(0.78067738) q[1];
rz(1.6171261) q[2];
sx q[2];
rz(-2.5382858) q[2];
sx q[2];
rz(-0.45679191) q[2];
rz(1.900832) q[3];
sx q[3];
rz(-1.5041372) q[3];
sx q[3];
rz(2.1048673) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
