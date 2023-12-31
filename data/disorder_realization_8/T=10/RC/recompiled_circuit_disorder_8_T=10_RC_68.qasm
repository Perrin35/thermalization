OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8060057) q[0];
sx q[0];
rz(-0.94526362) q[0];
sx q[0];
rz(-0.52559108) q[0];
rz(-2.8984012) q[1];
sx q[1];
rz(-1.2326198) q[1];
sx q[1];
rz(-0.90484172) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.605699) q[0];
sx q[0];
rz(-0.68471691) q[0];
sx q[0];
rz(-2.2936355) q[0];
rz(-pi) q[1];
rz(-1.7674255) q[2];
sx q[2];
rz(-1.7927205) q[2];
sx q[2];
rz(-1.0151528) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9368254) q[1];
sx q[1];
rz(-1.966404) q[1];
sx q[1];
rz(-0.25524615) q[1];
rz(1.2600793) q[3];
sx q[3];
rz(-1.3907392) q[3];
sx q[3];
rz(0.77375274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.20377542) q[2];
sx q[2];
rz(-1.4062466) q[2];
sx q[2];
rz(-0.096244372) q[2];
rz(2.105666) q[3];
sx q[3];
rz(-2.7544498) q[3];
sx q[3];
rz(2.9878785) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44089833) q[0];
sx q[0];
rz(-0.39114025) q[0];
sx q[0];
rz(-0.76517117) q[0];
rz(1.2922497) q[1];
sx q[1];
rz(-0.48520979) q[1];
sx q[1];
rz(-2.4786425) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63755858) q[0];
sx q[0];
rz(-0.090888977) q[0];
sx q[0];
rz(0.81764098) q[0];
rz(-0.40132482) q[2];
sx q[2];
rz(-0.95762816) q[2];
sx q[2];
rz(-2.3344628) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6389097) q[1];
sx q[1];
rz(-2.0109634) q[1];
sx q[1];
rz(-1.6619976) q[1];
x q[2];
rz(-2.6547673) q[3];
sx q[3];
rz(-1.7142222) q[3];
sx q[3];
rz(1.5407345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.092676) q[2];
sx q[2];
rz(-1.0485704) q[2];
sx q[2];
rz(-2.8125787) q[2];
rz(0.66550231) q[3];
sx q[3];
rz(-0.21829675) q[3];
sx q[3];
rz(1.3177419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5927785) q[0];
sx q[0];
rz(-0.22286649) q[0];
sx q[0];
rz(2.9192525) q[0];
rz(2.1242583) q[1];
sx q[1];
rz(-0.72128123) q[1];
sx q[1];
rz(0.51868784) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0761622) q[0];
sx q[0];
rz(-0.87834529) q[0];
sx q[0];
rz(2.6718219) q[0];
rz(-pi) q[1];
rz(1.6665886) q[2];
sx q[2];
rz(-1.3057858) q[2];
sx q[2];
rz(1.4694627) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.31836244) q[1];
sx q[1];
rz(-1.743411) q[1];
sx q[1];
rz(0.77962064) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5716001) q[3];
sx q[3];
rz(-1.5161848) q[3];
sx q[3];
rz(-2.2664203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2805933) q[2];
sx q[2];
rz(-2.7797647) q[2];
sx q[2];
rz(0.068543531) q[2];
rz(0.60244256) q[3];
sx q[3];
rz(-2.3790363) q[3];
sx q[3];
rz(-3.0025735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.65072) q[0];
sx q[0];
rz(-2.3644709) q[0];
sx q[0];
rz(-2.9673476) q[0];
rz(-0.53025591) q[1];
sx q[1];
rz(-1.4825876) q[1];
sx q[1];
rz(-2.6285016) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95272321) q[0];
sx q[0];
rz(-0.9523305) q[0];
sx q[0];
rz(1.4737211) q[0];
rz(-pi) q[1];
rz(-0.80715837) q[2];
sx q[2];
rz(-1.1045189) q[2];
sx q[2];
rz(1.6023139) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.063555524) q[1];
sx q[1];
rz(-2.6647898) q[1];
sx q[1];
rz(-0.30492353) q[1];
rz(-1.7170834) q[3];
sx q[3];
rz(-0.20855599) q[3];
sx q[3];
rz(-1.8640765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.47485581) q[2];
sx q[2];
rz(-2.7754144) q[2];
sx q[2];
rz(0.22988698) q[2];
rz(-2.722548) q[3];
sx q[3];
rz(-1.34904) q[3];
sx q[3];
rz(0.45927799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4693562) q[0];
sx q[0];
rz(-0.75119632) q[0];
sx q[0];
rz(-0.67681926) q[0];
rz(0.49304402) q[1];
sx q[1];
rz(-0.9489916) q[1];
sx q[1];
rz(2.5255323) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13000935) q[0];
sx q[0];
rz(-1.6066215) q[0];
sx q[0];
rz(1.623276) q[0];
x q[1];
rz(-0.30454238) q[2];
sx q[2];
rz(-2.7603622) q[2];
sx q[2];
rz(2.9074557) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7580326) q[1];
sx q[1];
rz(-2.5639503) q[1];
sx q[1];
rz(0.010096117) q[1];
x q[2];
rz(-2.2474399) q[3];
sx q[3];
rz(-1.3291306) q[3];
sx q[3];
rz(1.579293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.73413509) q[2];
sx q[2];
rz(-0.55441684) q[2];
sx q[2];
rz(-0.25536728) q[2];
rz(1.5363961) q[3];
sx q[3];
rz(-0.95241773) q[3];
sx q[3];
rz(-2.3674964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7191294) q[0];
sx q[0];
rz(-2.1753949) q[0];
sx q[0];
rz(0.47250026) q[0];
rz(2.6155112) q[1];
sx q[1];
rz(-0.20985797) q[1];
sx q[1];
rz(2.2568259) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.103325) q[0];
sx q[0];
rz(-1.190435) q[0];
sx q[0];
rz(-3.0618219) q[0];
x q[1];
rz(-2.9315345) q[2];
sx q[2];
rz(-1.5846328) q[2];
sx q[2];
rz(-2.7163598) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.51984519) q[1];
sx q[1];
rz(-1.2832844) q[1];
sx q[1];
rz(1.1430172) q[1];
rz(-pi) q[2];
rz(2.8956036) q[3];
sx q[3];
rz(-1.0503328) q[3];
sx q[3];
rz(-2.1465079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.74449599) q[2];
sx q[2];
rz(-0.8049736) q[2];
sx q[2];
rz(-2.8302144) q[2];
rz(1.7729676) q[3];
sx q[3];
rz(-2.6840648) q[3];
sx q[3];
rz(-2.6122724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7235274) q[0];
sx q[0];
rz(-2.8630246) q[0];
sx q[0];
rz(-0.061070651) q[0];
rz(-3.1014077) q[1];
sx q[1];
rz(-1.1611074) q[1];
sx q[1];
rz(-0.73289245) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31156763) q[0];
sx q[0];
rz(-0.68651474) q[0];
sx q[0];
rz(0.65450432) q[0];
x q[1];
rz(2.272846) q[2];
sx q[2];
rz(-1.1940496) q[2];
sx q[2];
rz(-1.4993315) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.43436189) q[1];
sx q[1];
rz(-1.5773298) q[1];
sx q[1];
rz(0.9450426) q[1];
rz(-1.2039127) q[3];
sx q[3];
rz(-0.11493472) q[3];
sx q[3];
rz(-0.68655187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0682893) q[2];
sx q[2];
rz(-1.9446334) q[2];
sx q[2];
rz(-0.51458365) q[2];
rz(-1.2375281) q[3];
sx q[3];
rz(-2.5585744) q[3];
sx q[3];
rz(-2.5966743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.085389) q[0];
sx q[0];
rz(-1.6563002) q[0];
sx q[0];
rz(-2.432166) q[0];
rz(1.6363232) q[1];
sx q[1];
rz(-2.067833) q[1];
sx q[1];
rz(2.8628796) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7473135) q[0];
sx q[0];
rz(-2.8907667) q[0];
sx q[0];
rz(-2.5670811) q[0];
rz(-0.3785554) q[2];
sx q[2];
rz(-2.5494529) q[2];
sx q[2];
rz(-0.90422599) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.41234327) q[1];
sx q[1];
rz(-1.7885498) q[1];
sx q[1];
rz(-1.8885683) q[1];
x q[2];
rz(-1.9023499) q[3];
sx q[3];
rz(-0.23570508) q[3];
sx q[3];
rz(1.2329668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8401106) q[2];
sx q[2];
rz(-1.297941) q[2];
sx q[2];
rz(3.0855132) q[2];
rz(-2.2864443) q[3];
sx q[3];
rz(-0.44848281) q[3];
sx q[3];
rz(-0.40518951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1466325) q[0];
sx q[0];
rz(-2.9652847) q[0];
sx q[0];
rz(-0.96889281) q[0];
rz(2.6682207) q[1];
sx q[1];
rz(-2.3469766) q[1];
sx q[1];
rz(1.1425346) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56209598) q[0];
sx q[0];
rz(-1.8651433) q[0];
sx q[0];
rz(-3.0466945) q[0];
x q[1];
rz(1.5732473) q[2];
sx q[2];
rz(-2.2741389) q[2];
sx q[2];
rz(-3.0424812) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7361703) q[1];
sx q[1];
rz(-1.0338963) q[1];
sx q[1];
rz(1.7169397) q[1];
x q[2];
rz(-2.9386018) q[3];
sx q[3];
rz(-0.91068017) q[3];
sx q[3];
rz(0.70860329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.896495) q[2];
sx q[2];
rz(-1.2197887) q[2];
sx q[2];
rz(-1.0207821) q[2];
rz(-2.8178689) q[3];
sx q[3];
rz(-0.75298572) q[3];
sx q[3];
rz(3.1304205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97994119) q[0];
sx q[0];
rz(-0.027898235) q[0];
sx q[0];
rz(-0.7014057) q[0];
rz(-0.91570634) q[1];
sx q[1];
rz(-1.0083102) q[1];
sx q[1];
rz(1.2385626) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9311213) q[0];
sx q[0];
rz(-2.5744372) q[0];
sx q[0];
rz(1.9803067) q[0];
rz(-pi) q[1];
rz(1.5194703) q[2];
sx q[2];
rz(-2.3231299) q[2];
sx q[2];
rz(0.78879702) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8348332) q[1];
sx q[1];
rz(-2.4750437) q[1];
sx q[1];
rz(-1.3798316) q[1];
x q[2];
rz(-2.1632807) q[3];
sx q[3];
rz(-1.6654135) q[3];
sx q[3];
rz(2.832151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.68676585) q[2];
sx q[2];
rz(-2.0451615) q[2];
sx q[2];
rz(0.043838538) q[2];
rz(1.2010126) q[3];
sx q[3];
rz(-0.73533708) q[3];
sx q[3];
rz(2.1380077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6702406) q[0];
sx q[0];
rz(-2.4155407) q[0];
sx q[0];
rz(1.7759905) q[0];
rz(-0.025370601) q[1];
sx q[1];
rz(-1.3062968) q[1];
sx q[1];
rz(1.2702373) q[1];
rz(2.2123443) q[2];
sx q[2];
rz(-0.48454787) q[2];
sx q[2];
rz(1.5018644) q[2];
rz(-0.17037114) q[3];
sx q[3];
rz(-0.80633612) q[3];
sx q[3];
rz(2.7663305) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
