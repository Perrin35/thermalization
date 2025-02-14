OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.1342993) q[0];
sx q[0];
rz(8.2259895) q[0];
sx q[0];
rz(10.209957) q[0];
rz(2.5685318) q[1];
sx q[1];
rz(2.1905724) q[1];
sx q[1];
rz(11.934927) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.310869) q[0];
sx q[0];
rz(-1.1397624) q[0];
sx q[0];
rz(2.2742989) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78978059) q[2];
sx q[2];
rz(-2.1319816) q[2];
sx q[2];
rz(-1.9927858) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8353442) q[1];
sx q[1];
rz(-0.95148309) q[1];
sx q[1];
rz(2.1080416) q[1];
rz(-pi) q[2];
rz(2.9358125) q[3];
sx q[3];
rz(-2.6884086) q[3];
sx q[3];
rz(3.1169716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.50420061) q[2];
sx q[2];
rz(-1.7039508) q[2];
sx q[2];
rz(-0.27153095) q[2];
rz(-3.0951989) q[3];
sx q[3];
rz(-1.7312739) q[3];
sx q[3];
rz(-0.36762777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72535998) q[0];
sx q[0];
rz(-0.027149057) q[0];
sx q[0];
rz(2.695988) q[0];
rz(2.7815869) q[1];
sx q[1];
rz(-1.5841192) q[1];
sx q[1];
rz(0.97703591) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.925219) q[0];
sx q[0];
rz(-2.3686103) q[0];
sx q[0];
rz(-2.3712914) q[0];
rz(-0.71828281) q[2];
sx q[2];
rz(-2.278285) q[2];
sx q[2];
rz(1.2808509) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.313873) q[1];
sx q[1];
rz(-1.4196175) q[1];
sx q[1];
rz(0.97444356) q[1];
rz(-pi) q[2];
rz(-2.4384624) q[3];
sx q[3];
rz(-1.0323914) q[3];
sx q[3];
rz(-2.3634274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.52686024) q[2];
sx q[2];
rz(-1.3440259) q[2];
sx q[2];
rz(-0.942222) q[2];
rz(-0.49643907) q[3];
sx q[3];
rz(-1.6783291) q[3];
sx q[3];
rz(-2.0739323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0577069) q[0];
sx q[0];
rz(-1.3297465) q[0];
sx q[0];
rz(0.65621334) q[0];
rz(-2.6612142) q[1];
sx q[1];
rz(-2.1800133) q[1];
sx q[1];
rz(-2.5286455) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71343173) q[0];
sx q[0];
rz(-1.1640004) q[0];
sx q[0];
rz(-3.0359546) q[0];
rz(-pi) q[1];
rz(-3.1334932) q[2];
sx q[2];
rz(-0.73609422) q[2];
sx q[2];
rz(1.8437633) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7483116) q[1];
sx q[1];
rz(-1.5996115) q[1];
sx q[1];
rz(3.1080212) q[1];
rz(-pi) q[2];
rz(-2.6426389) q[3];
sx q[3];
rz(-1.6315579) q[3];
sx q[3];
rz(2.2419597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2059242) q[2];
sx q[2];
rz(-1.0574477) q[2];
sx q[2];
rz(-2.617188) q[2];
rz(2.2762401) q[3];
sx q[3];
rz(-2.083358) q[3];
sx q[3];
rz(0.66044468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.48786369) q[0];
sx q[0];
rz(-1.8346584) q[0];
sx q[0];
rz(0.096268795) q[0];
rz(-3.1277711) q[1];
sx q[1];
rz(-0.50058573) q[1];
sx q[1];
rz(1.0523419) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6276564) q[0];
sx q[0];
rz(-3.0110783) q[0];
sx q[0];
rz(-1.9399727) q[0];
rz(-pi) q[1];
x q[1];
rz(1.123548) q[2];
sx q[2];
rz(-2.5442225) q[2];
sx q[2];
rz(-1.1505466) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2675002) q[1];
sx q[1];
rz(-0.37768957) q[1];
sx q[1];
rz(-0.51376359) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0963444) q[3];
sx q[3];
rz(-1.6017672) q[3];
sx q[3];
rz(2.3207959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4200165) q[2];
sx q[2];
rz(-0.53109157) q[2];
sx q[2];
rz(2.6225923) q[2];
rz(-2.0429933) q[3];
sx q[3];
rz(-1.4562621) q[3];
sx q[3];
rz(2.4115653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5191583) q[0];
sx q[0];
rz(-2.9485478) q[0];
sx q[0];
rz(0.73424196) q[0];
rz(1.2892067) q[1];
sx q[1];
rz(-1.0451008) q[1];
sx q[1];
rz(0.9085013) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72469456) q[0];
sx q[0];
rz(-2.1229345) q[0];
sx q[0];
rz(0.79663251) q[0];
rz(-0.20573767) q[2];
sx q[2];
rz(-2.5605953) q[2];
sx q[2];
rz(0.28827039) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.1724194) q[1];
sx q[1];
rz(-0.73639077) q[1];
sx q[1];
rz(-1.4672155) q[1];
rz(-pi) q[2];
rz(1.8374422) q[3];
sx q[3];
rz(-1.5943267) q[3];
sx q[3];
rz(1.2465828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.044540731) q[2];
sx q[2];
rz(-0.78623199) q[2];
sx q[2];
rz(-2.457974) q[2];
rz(-2.0475856) q[3];
sx q[3];
rz(-1.8180327) q[3];
sx q[3];
rz(-1.5692284) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.053726824) q[0];
sx q[0];
rz(-1.6931067) q[0];
sx q[0];
rz(-2.7159765) q[0];
rz(-2.5348237) q[1];
sx q[1];
rz(-2.4538222) q[1];
sx q[1];
rz(-1.2026131) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15213293) q[0];
sx q[0];
rz(-2.0589863) q[0];
sx q[0];
rz(-0.91849416) q[0];
x q[1];
rz(1.5853556) q[2];
sx q[2];
rz(-1.636616) q[2];
sx q[2];
rz(-1.751339) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.21638097) q[1];
sx q[1];
rz(-0.93384508) q[1];
sx q[1];
rz(0.94016192) q[1];
rz(-pi) q[2];
rz(-1.3404863) q[3];
sx q[3];
rz(-1.7969683) q[3];
sx q[3];
rz(0.61128547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.46835029) q[2];
sx q[2];
rz(-2.612256) q[2];
sx q[2];
rz(0.62874111) q[2];
rz(-2.8041503) q[3];
sx q[3];
rz(-1.3597666) q[3];
sx q[3];
rz(-0.76702816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5479946) q[0];
sx q[0];
rz(-0.070040528) q[0];
sx q[0];
rz(1.8016169) q[0];
rz(-0.10546671) q[1];
sx q[1];
rz(-2.1601845) q[1];
sx q[1];
rz(0.14428219) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54050416) q[0];
sx q[0];
rz(-1.2799037) q[0];
sx q[0];
rz(0.40368794) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7970798) q[2];
sx q[2];
rz(-2.4206875) q[2];
sx q[2];
rz(-0.048335942) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6444417) q[1];
sx q[1];
rz(-1.2220009) q[1];
sx q[1];
rz(-1.3494874) q[1];
rz(-pi) q[2];
rz(0.65001092) q[3];
sx q[3];
rz(-0.78449215) q[3];
sx q[3];
rz(-0.9264526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0353126) q[2];
sx q[2];
rz(-1.148618) q[2];
sx q[2];
rz(2.3036352) q[2];
rz(-0.49977866) q[3];
sx q[3];
rz(-0.79334799) q[3];
sx q[3];
rz(2.4404073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4226294) q[0];
sx q[0];
rz(-2.0327649) q[0];
sx q[0];
rz(0.75123373) q[0];
rz(1.2535837) q[1];
sx q[1];
rz(-1.3026078) q[1];
sx q[1];
rz(-1.7734897) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4104162) q[0];
sx q[0];
rz(-0.77922693) q[0];
sx q[0];
rz(-0.91996361) q[0];
rz(-pi) q[1];
rz(1.8644731) q[2];
sx q[2];
rz(-1.2674567) q[2];
sx q[2];
rz(-1.9544698) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82097365) q[1];
sx q[1];
rz(-1.6553967) q[1];
sx q[1];
rz(-2.3449347) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7799928) q[3];
sx q[3];
rz(-1.0269093) q[3];
sx q[3];
rz(-1.7392478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.83296835) q[2];
sx q[2];
rz(-0.37165752) q[2];
sx q[2];
rz(0.0011681636) q[2];
rz(-0.33231169) q[3];
sx q[3];
rz(-1.6806335) q[3];
sx q[3];
rz(-0.47962475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7697656) q[0];
sx q[0];
rz(-1.9244939) q[0];
sx q[0];
rz(2.0455072) q[0];
rz(-1.7425884) q[1];
sx q[1];
rz(-1.636248) q[1];
sx q[1];
rz(-2.2217264) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70912659) q[0];
sx q[0];
rz(-1.6133012) q[0];
sx q[0];
rz(2.0086238) q[0];
x q[1];
rz(-0.28025744) q[2];
sx q[2];
rz(-1.1667023) q[2];
sx q[2];
rz(-1.0755444) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.87756094) q[1];
sx q[1];
rz(-1.4204867) q[1];
sx q[1];
rz(0.64964215) q[1];
rz(2.1036622) q[3];
sx q[3];
rz(-1.6954331) q[3];
sx q[3];
rz(0.79844284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2190735) q[2];
sx q[2];
rz(-0.24524958) q[2];
sx q[2];
rz(-1.5101439) q[2];
rz(-2.505693) q[3];
sx q[3];
rz(-0.70521516) q[3];
sx q[3];
rz(-2.6270134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.1849599) q[0];
sx q[0];
rz(-2.2745467) q[0];
sx q[0];
rz(1.6459203) q[0];
rz(-0.11790672) q[1];
sx q[1];
rz(-1.7472183) q[1];
sx q[1];
rz(-2.8848021) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8816329) q[0];
sx q[0];
rz(-1.6266394) q[0];
sx q[0];
rz(1.1373489) q[0];
rz(1.8689972) q[2];
sx q[2];
rz(-2.1352945) q[2];
sx q[2];
rz(0.80882458) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.28370441) q[1];
sx q[1];
rz(-1.8562177) q[1];
sx q[1];
rz(-1.2976451) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2076693) q[3];
sx q[3];
rz(-1.052891) q[3];
sx q[3];
rz(0.3929485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.39766496) q[2];
sx q[2];
rz(-2.3981018) q[2];
sx q[2];
rz(2.1916154) q[2];
rz(2.5395565) q[3];
sx q[3];
rz(-1.9046015) q[3];
sx q[3];
rz(0.33815798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6947185) q[0];
sx q[0];
rz(-0.61114408) q[0];
sx q[0];
rz(-1.1993988) q[0];
rz(-1.0097722) q[1];
sx q[1];
rz(-2.2191681) q[1];
sx q[1];
rz(2.9422599) q[1];
rz(1.8795602) q[2];
sx q[2];
rz(-1.8539351) q[2];
sx q[2];
rz(-1.2507982) q[2];
rz(0.77824577) q[3];
sx q[3];
rz(-1.2846071) q[3];
sx q[3];
rz(-1.0461109) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
