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
rz(-1.4200014) q[0];
sx q[0];
rz(5.3237259) q[0];
sx q[0];
rz(11.525679) q[0];
rz(2.8031082) q[1];
sx q[1];
rz(-1.6522633) q[1];
sx q[1];
rz(-1.7117865) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.734935) q[0];
sx q[0];
rz(-1.8419319) q[0];
sx q[0];
rz(-1.447729) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0645039) q[2];
sx q[2];
rz(-0.66254751) q[2];
sx q[2];
rz(0.27499378) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0104645) q[1];
sx q[1];
rz(-2.2704723) q[1];
sx q[1];
rz(1.9424428) q[1];
x q[2];
rz(2.7316507) q[3];
sx q[3];
rz(-0.60527181) q[3];
sx q[3];
rz(-2.2469478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1374984) q[2];
sx q[2];
rz(-1.4278922) q[2];
sx q[2];
rz(-0.049169866) q[2];
rz(-0.54685012) q[3];
sx q[3];
rz(-2.2873736) q[3];
sx q[3];
rz(-1.2949519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8787254) q[0];
sx q[0];
rz(-0.85064369) q[0];
sx q[0];
rz(-0.039462939) q[0];
rz(0.70611686) q[1];
sx q[1];
rz(-1.9673653) q[1];
sx q[1];
rz(1.2605234) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17968125) q[0];
sx q[0];
rz(-2.7196809) q[0];
sx q[0];
rz(1.4279813) q[0];
rz(-2.4751365) q[2];
sx q[2];
rz(-2.0100231) q[2];
sx q[2];
rz(-1.1117488) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9866922) q[1];
sx q[1];
rz(-3.0873484) q[1];
sx q[1];
rz(-2.2510347) q[1];
rz(-pi) q[2];
rz(-0.19393215) q[3];
sx q[3];
rz(-1.5054879) q[3];
sx q[3];
rz(2.0146305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4827106) q[2];
sx q[2];
rz(-1.005268) q[2];
sx q[2];
rz(-0.91747326) q[2];
rz(2.7190123) q[3];
sx q[3];
rz(-1.249142) q[3];
sx q[3];
rz(-1.0379855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4611886) q[0];
sx q[0];
rz(-2.5867511) q[0];
sx q[0];
rz(-2.19221) q[0];
rz(1.8934911) q[1];
sx q[1];
rz(-0.58369842) q[1];
sx q[1];
rz(-1.6277574) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7010403) q[0];
sx q[0];
rz(-1.3258385) q[0];
sx q[0];
rz(2.7666758) q[0];
rz(-pi) q[1];
rz(0.21612303) q[2];
sx q[2];
rz(-2.0449567) q[2];
sx q[2];
rz(0.40008998) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.8366295) q[1];
sx q[1];
rz(-2.8439972) q[1];
sx q[1];
rz(2.8862428) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.85905) q[3];
sx q[3];
rz(-1.4465528) q[3];
sx q[3];
rz(-1.8831203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.44024399) q[2];
sx q[2];
rz(-1.1260208) q[2];
sx q[2];
rz(1.5234285) q[2];
rz(2.8273888) q[3];
sx q[3];
rz(-2.115963) q[3];
sx q[3];
rz(-1.6395125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0358129) q[0];
sx q[0];
rz(-2.950225) q[0];
sx q[0];
rz(-2.9845003) q[0];
rz(3.0063903) q[1];
sx q[1];
rz(-2.356485) q[1];
sx q[1];
rz(2.8083727) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5557809) q[0];
sx q[0];
rz(-0.69792047) q[0];
sx q[0];
rz(-2.3960953) q[0];
x q[1];
rz(2.2991247) q[2];
sx q[2];
rz(-1.4204475) q[2];
sx q[2];
rz(1.6347234) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.0097447952) q[1];
sx q[1];
rz(-1.3074083) q[1];
sx q[1];
rz(-1.2089648) q[1];
rz(-pi) q[2];
rz(-0.31308324) q[3];
sx q[3];
rz(-0.88290324) q[3];
sx q[3];
rz(0.41447251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.80804431) q[2];
sx q[2];
rz(-2.4925888) q[2];
sx q[2];
rz(-1.6928847) q[2];
rz(0.73271218) q[3];
sx q[3];
rz(-1.2196187) q[3];
sx q[3];
rz(-1.609751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44959679) q[0];
sx q[0];
rz(-1.6549598) q[0];
sx q[0];
rz(-2.9100371) q[0];
rz(2.4202276) q[1];
sx q[1];
rz(-2.2504579) q[1];
sx q[1];
rz(2.2931633) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4130294) q[0];
sx q[0];
rz(-0.60270488) q[0];
sx q[0];
rz(-2.187285) q[0];
rz(-pi) q[1];
rz(0.5428585) q[2];
sx q[2];
rz(-1.2194389) q[2];
sx q[2];
rz(1.5962102) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1047804) q[1];
sx q[1];
rz(-1.1835349) q[1];
sx q[1];
rz(-1.4951863) q[1];
rz(0.23691688) q[3];
sx q[3];
rz(-1.6730783) q[3];
sx q[3];
rz(-1.8558242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1344177) q[2];
sx q[2];
rz(-2.538372) q[2];
sx q[2];
rz(1.940654) q[2];
rz(2.4522771) q[3];
sx q[3];
rz(-1.9715693) q[3];
sx q[3];
rz(0.56572604) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1597964) q[0];
sx q[0];
rz(-1.5317651) q[0];
sx q[0];
rz(1.7933886) q[0];
rz(-0.7181522) q[1];
sx q[1];
rz(-2.5703057) q[1];
sx q[1];
rz(1.2917554) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51707089) q[0];
sx q[0];
rz(-0.41247955) q[0];
sx q[0];
rz(-0.58858354) q[0];
rz(2.767105) q[2];
sx q[2];
rz(-1.0051703) q[2];
sx q[2];
rz(-2.781812) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3499432) q[1];
sx q[1];
rz(-1.4022489) q[1];
sx q[1];
rz(1.6577634) q[1];
rz(0.060543493) q[3];
sx q[3];
rz(-2.745564) q[3];
sx q[3];
rz(-2.372449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9888088) q[2];
sx q[2];
rz(-0.96502105) q[2];
sx q[2];
rz(-0.40210813) q[2];
rz(2.1357338) q[3];
sx q[3];
rz(-2.918225) q[3];
sx q[3];
rz(-1.4373826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-2.724029) q[0];
sx q[0];
rz(-0.14823866) q[0];
sx q[0];
rz(1.8120026) q[0];
rz(-1.3279462) q[1];
sx q[1];
rz(-1.7247533) q[1];
sx q[1];
rz(-0.45164576) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62147776) q[0];
sx q[0];
rz(-1.661944) q[0];
sx q[0];
rz(0.4439179) q[0];
rz(-pi) q[1];
x q[1];
rz(0.30173413) q[2];
sx q[2];
rz(-2.6987751) q[2];
sx q[2];
rz(-2.2254965) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0942119) q[1];
sx q[1];
rz(-1.3806731) q[1];
sx q[1];
rz(-3.090108) q[1];
rz(-pi) q[2];
rz(-2.3450646) q[3];
sx q[3];
rz(-2.3083271) q[3];
sx q[3];
rz(0.11175534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.20246501) q[2];
sx q[2];
rz(-1.7697325) q[2];
sx q[2];
rz(-0.28787127) q[2];
rz(-1.4816083) q[3];
sx q[3];
rz(-0.85166728) q[3];
sx q[3];
rz(1.456858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(2.562029) q[0];
sx q[0];
rz(-0.96051878) q[0];
sx q[0];
rz(-2.5886743) q[0];
rz(1.2506073) q[1];
sx q[1];
rz(-2.7313373) q[1];
sx q[1];
rz(-1.7611354) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.834873) q[0];
sx q[0];
rz(-1.0579029) q[0];
sx q[0];
rz(-1.9901278) q[0];
rz(-3.0199261) q[2];
sx q[2];
rz(-0.98699283) q[2];
sx q[2];
rz(-0.19106272) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.41696527) q[1];
sx q[1];
rz(-1.529602) q[1];
sx q[1];
rz(-0.23708701) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16726027) q[3];
sx q[3];
rz(-2.3338807) q[3];
sx q[3];
rz(2.4204262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5914501) q[2];
sx q[2];
rz(-1.7937135) q[2];
sx q[2];
rz(0.71510092) q[2];
rz(3.0969369) q[3];
sx q[3];
rz(-0.5831334) q[3];
sx q[3];
rz(-2.513212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1735246) q[0];
sx q[0];
rz(-2.3723497) q[0];
sx q[0];
rz(2.5007057) q[0];
rz(1.6454654) q[1];
sx q[1];
rz(-2.75664) q[1];
sx q[1];
rz(2.4900751) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9891883) q[0];
sx q[0];
rz(-1.8274283) q[0];
sx q[0];
rz(2.2874831) q[0];
x q[1];
rz(-1.8046494) q[2];
sx q[2];
rz(-1.0923947) q[2];
sx q[2];
rz(-3.0870147) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99371893) q[1];
sx q[1];
rz(-0.9269971) q[1];
sx q[1];
rz(0.46624513) q[1];
rz(-pi) q[2];
rz(1.9189758) q[3];
sx q[3];
rz(-1.3908252) q[3];
sx q[3];
rz(2.6258385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7351825) q[2];
sx q[2];
rz(-1.0446905) q[2];
sx q[2];
rz(0.45808211) q[2];
rz(3.1297019) q[3];
sx q[3];
rz(-1.0606822) q[3];
sx q[3];
rz(0.021765821) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7387725) q[0];
sx q[0];
rz(-2.9467376) q[0];
sx q[0];
rz(3.1126157) q[0];
rz(0.92652357) q[1];
sx q[1];
rz(-1.4484826) q[1];
sx q[1];
rz(2.7353333) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0800015) q[0];
sx q[0];
rz(-1.5621981) q[0];
sx q[0];
rz(0.28010578) q[0];
rz(-0.92490893) q[2];
sx q[2];
rz(-1.8573833) q[2];
sx q[2];
rz(-2.1429581) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.528842) q[1];
sx q[1];
rz(-1.6019197) q[1];
sx q[1];
rz(1.2791388) q[1];
rz(-pi) q[2];
rz(-0.80244949) q[3];
sx q[3];
rz(-0.93367773) q[3];
sx q[3];
rz(-3.1010026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77068344) q[2];
sx q[2];
rz(-1.4981046) q[2];
sx q[2];
rz(-0.59477273) q[2];
rz(2.4316783) q[3];
sx q[3];
rz(-2.1845332) q[3];
sx q[3];
rz(2.1962568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
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
rz(-0.57472537) q[0];
sx q[0];
rz(-2.18676) q[0];
sx q[0];
rz(-2.2796897) q[0];
rz(0.30543874) q[1];
sx q[1];
rz(-1.8157235) q[1];
sx q[1];
rz(-2.9095412) q[1];
rz(-1.1535063) q[2];
sx q[2];
rz(-1.8883033) q[2];
sx q[2];
rz(0.89486833) q[2];
rz(-1.7954682) q[3];
sx q[3];
rz(-1.9763038) q[3];
sx q[3];
rz(1.4122813) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
