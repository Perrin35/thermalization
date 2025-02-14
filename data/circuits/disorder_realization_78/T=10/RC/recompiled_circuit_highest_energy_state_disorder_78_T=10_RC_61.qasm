OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3820833) q[0];
sx q[0];
rz(-0.74892646) q[0];
sx q[0];
rz(1.7480667) q[0];
rz(-2.3387609) q[1];
sx q[1];
rz(3.9376942) q[1];
sx q[1];
rz(12.035523) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33541629) q[0];
sx q[0];
rz(-0.62549504) q[0];
sx q[0];
rz(2.2124887) q[0];
rz(-0.10239281) q[2];
sx q[2];
rz(-0.16764417) q[2];
sx q[2];
rz(1.8554115) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.894132) q[1];
sx q[1];
rz(-1.0399721) q[1];
sx q[1];
rz(0.08695072) q[1];
rz(-0.10194709) q[3];
sx q[3];
rz(-1.0539712) q[3];
sx q[3];
rz(2.068813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.73534766) q[2];
sx q[2];
rz(-0.9333868) q[2];
sx q[2];
rz(2.1437342) q[2];
rz(-0.88456279) q[3];
sx q[3];
rz(-1.9622784) q[3];
sx q[3];
rz(-2.4295889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74580055) q[0];
sx q[0];
rz(-2.6048248) q[0];
sx q[0];
rz(-1.0774379) q[0];
rz(-2.5224345) q[1];
sx q[1];
rz(-1.2666603) q[1];
sx q[1];
rz(1.223863) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3676476) q[0];
sx q[0];
rz(-1.2804693) q[0];
sx q[0];
rz(-2.9591857) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6009839) q[2];
sx q[2];
rz(-1.5726302) q[2];
sx q[2];
rz(0.48170127) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2457285) q[1];
sx q[1];
rz(-0.99555496) q[1];
sx q[1];
rz(-2.0644581) q[1];
rz(-pi) q[2];
rz(0.84544648) q[3];
sx q[3];
rz(-2.2152112) q[3];
sx q[3];
rz(3.0385142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.18664843) q[2];
sx q[2];
rz(-0.84543219) q[2];
sx q[2];
rz(-1.0364944) q[2];
rz(1.9519818) q[3];
sx q[3];
rz(-1.9018973) q[3];
sx q[3];
rz(3.0724683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4619197) q[0];
sx q[0];
rz(-2.9000977) q[0];
sx q[0];
rz(-1.4659721) q[0];
rz(-1.8849323) q[1];
sx q[1];
rz(-2.4938221) q[1];
sx q[1];
rz(0.57674903) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3925331) q[0];
sx q[0];
rz(-2.1196094) q[0];
sx q[0];
rz(-2.4486923) q[0];
rz(2.6827963) q[2];
sx q[2];
rz(-1.5374628) q[2];
sx q[2];
rz(-1.3622487) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8654792) q[1];
sx q[1];
rz(-0.41655585) q[1];
sx q[1];
rz(1.4426484) q[1];
x q[2];
rz(-1.7347832) q[3];
sx q[3];
rz(-1.6583558) q[3];
sx q[3];
rz(-1.6370186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3999148) q[2];
sx q[2];
rz(-0.33881131) q[2];
sx q[2];
rz(-2.3728288) q[2];
rz(-0.1712884) q[3];
sx q[3];
rz(-1.4408709) q[3];
sx q[3];
rz(-2.7870074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4549558) q[0];
sx q[0];
rz(-2.2577715) q[0];
sx q[0];
rz(2.6052642) q[0];
rz(-2.3388011) q[1];
sx q[1];
rz(-1.2840459) q[1];
sx q[1];
rz(3.0435496) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6199099) q[0];
sx q[0];
rz(-0.3371211) q[0];
sx q[0];
rz(1.328892) q[0];
rz(-pi) q[1];
rz(-0.78196193) q[2];
sx q[2];
rz(-2.0148009) q[2];
sx q[2];
rz(2.9733016) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.72934276) q[1];
sx q[1];
rz(-0.55604425) q[1];
sx q[1];
rz(-1.1792609) q[1];
rz(-pi) q[2];
rz(-2.3290702) q[3];
sx q[3];
rz(-0.56216633) q[3];
sx q[3];
rz(0.49725975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98196491) q[2];
sx q[2];
rz(-2.1106796) q[2];
sx q[2];
rz(-0.31943303) q[2];
rz(-0.56644136) q[3];
sx q[3];
rz(-0.74685493) q[3];
sx q[3];
rz(1.1613065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7403858) q[0];
sx q[0];
rz(-1.8219319) q[0];
sx q[0];
rz(2.5886986) q[0];
rz(-0.7792019) q[1];
sx q[1];
rz(-0.82256493) q[1];
sx q[1];
rz(-0.78831569) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4944276) q[0];
sx q[0];
rz(-0.26965037) q[0];
sx q[0];
rz(-0.2616051) q[0];
rz(-pi) q[1];
rz(-3.0296828) q[2];
sx q[2];
rz(-2.9315278) q[2];
sx q[2];
rz(-0.0029759759) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.32994575) q[1];
sx q[1];
rz(-1.1378985) q[1];
sx q[1];
rz(-1.2654773) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8184999) q[3];
sx q[3];
rz(-1.8616397) q[3];
sx q[3];
rz(1.3210306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.39764443) q[2];
sx q[2];
rz(-1.1279736) q[2];
sx q[2];
rz(2.9359342) q[2];
rz(1.7913943) q[3];
sx q[3];
rz(-0.74068991) q[3];
sx q[3];
rz(3.1309483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3288997) q[0];
sx q[0];
rz(-2.6943272) q[0];
sx q[0];
rz(1.5146259) q[0];
rz(2.336592) q[1];
sx q[1];
rz(-1.6945508) q[1];
sx q[1];
rz(2.8256493) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7132063) q[0];
sx q[0];
rz(-1.9613736) q[0];
sx q[0];
rz(-1.9414563) q[0];
rz(1.1996026) q[2];
sx q[2];
rz(-2.1495869) q[2];
sx q[2];
rz(2.3166192) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.22623569) q[1];
sx q[1];
rz(-1.8457883) q[1];
sx q[1];
rz(2.8746241) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31031761) q[3];
sx q[3];
rz(-1.7162706) q[3];
sx q[3];
rz(3.1273354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.67779764) q[2];
sx q[2];
rz(-0.69910502) q[2];
sx q[2];
rz(-0.15764906) q[2];
rz(2.6080103) q[3];
sx q[3];
rz(-1.6505417) q[3];
sx q[3];
rz(1.6238448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5428298) q[0];
sx q[0];
rz(-2.6451126) q[0];
sx q[0];
rz(0.98094034) q[0];
rz(-2.6988103) q[1];
sx q[1];
rz(-1.9552224) q[1];
sx q[1];
rz(0.47169366) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2344246) q[0];
sx q[0];
rz(-1.2786073) q[0];
sx q[0];
rz(2.0227595) q[0];
rz(0.27916698) q[2];
sx q[2];
rz(-1.4067603) q[2];
sx q[2];
rz(-3.0085473) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.16817465) q[1];
sx q[1];
rz(-1.4354295) q[1];
sx q[1];
rz(-1.584704) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6568233) q[3];
sx q[3];
rz(-0.38082424) q[3];
sx q[3];
rz(-2.4172779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9337351) q[2];
sx q[2];
rz(-2.2722878) q[2];
sx q[2];
rz(-0.11824879) q[2];
rz(0.56784981) q[3];
sx q[3];
rz(-1.7852781) q[3];
sx q[3];
rz(0.80287272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0589013) q[0];
sx q[0];
rz(-1.7401798) q[0];
sx q[0];
rz(2.4543104) q[0];
rz(1.2673238) q[1];
sx q[1];
rz(-1.2143538) q[1];
sx q[1];
rz(-0.6257239) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95803496) q[0];
sx q[0];
rz(-1.5661217) q[0];
sx q[0];
rz(-1.5918533) q[0];
x q[1];
rz(2.2082616) q[2];
sx q[2];
rz(-1.0984761) q[2];
sx q[2];
rz(-2.460091) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-pi/4) q[1];
sx q[1];
rz(-2.56678) q[1];
sx q[1];
rz(1.5479799) q[1];
rz(2.1504907) q[3];
sx q[3];
rz(-1.3091581) q[3];
sx q[3];
rz(-1.8707448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.22831336) q[2];
sx q[2];
rz(-1.8931171) q[2];
sx q[2];
rz(1.0949562) q[2];
rz(2.6878808) q[3];
sx q[3];
rz(-2.3504421) q[3];
sx q[3];
rz(2.9554844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5943282) q[0];
sx q[0];
rz(-0.493395) q[0];
sx q[0];
rz(-0.067597978) q[0];
rz(-1.0293845) q[1];
sx q[1];
rz(-1.8600347) q[1];
sx q[1];
rz(-3.0060815) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1757723) q[0];
sx q[0];
rz(-1.0403353) q[0];
sx q[0];
rz(-0.6483174) q[0];
x q[1];
rz(-2.9642815) q[2];
sx q[2];
rz(-2.3925892) q[2];
sx q[2];
rz(2.2875691) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1780336) q[1];
sx q[1];
rz(-1.8521223) q[1];
sx q[1];
rz(0.60124361) q[1];
x q[2];
rz(2.4289651) q[3];
sx q[3];
rz(-0.36667675) q[3];
sx q[3];
rz(-2.3617876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6873261) q[2];
sx q[2];
rz(-2.8886075) q[2];
sx q[2];
rz(-2.6542286) q[2];
rz(0.38960114) q[3];
sx q[3];
rz(-0.95249683) q[3];
sx q[3];
rz(-3.0756557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.49366632) q[0];
sx q[0];
rz(-2.0559897) q[0];
sx q[0];
rz(-1.3847463) q[0];
rz(-2.0866277) q[1];
sx q[1];
rz(-0.87433785) q[1];
sx q[1];
rz(-0.25996444) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4892762) q[0];
sx q[0];
rz(-0.92149177) q[0];
sx q[0];
rz(0.32353521) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0412156) q[2];
sx q[2];
rz(-1.4401199) q[2];
sx q[2];
rz(1.5575156) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.00947) q[1];
sx q[1];
rz(-0.94317052) q[1];
sx q[1];
rz(-0.043379003) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22180827) q[3];
sx q[3];
rz(-1.3768428) q[3];
sx q[3];
rz(-1.8499595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7197623) q[2];
sx q[2];
rz(-2.7037342) q[2];
sx q[2];
rz(1.7913294) q[2];
rz(-2.0566025) q[3];
sx q[3];
rz(-1.2144054) q[3];
sx q[3];
rz(1.1291645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0706901) q[0];
sx q[0];
rz(-2.4546843) q[0];
sx q[0];
rz(-0.51921459) q[0];
rz(-1.6593973) q[1];
sx q[1];
rz(-1.2980325) q[1];
sx q[1];
rz(1.4549805) q[1];
rz(2.0878592) q[2];
sx q[2];
rz(-1.0521493) q[2];
sx q[2];
rz(-0.041368749) q[2];
rz(2.8419513) q[3];
sx q[3];
rz(-1.089923) q[3];
sx q[3];
rz(-3.0724473) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
