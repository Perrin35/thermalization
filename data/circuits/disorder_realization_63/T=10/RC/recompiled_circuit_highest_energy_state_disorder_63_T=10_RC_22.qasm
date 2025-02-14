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
rz(-1.7173362) q[0];
sx q[0];
rz(-1.9219226) q[0];
sx q[0];
rz(-2.5665459) q[0];
rz(0.14248928) q[1];
sx q[1];
rz(-1.9228851) q[1];
sx q[1];
rz(-0.86427468) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48102114) q[0];
sx q[0];
rz(-2.8459274) q[0];
sx q[0];
rz(0.85938248) q[0];
rz(-2.2021212) q[2];
sx q[2];
rz(-1.4465904) q[2];
sx q[2];
rz(0.4865464) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6576523) q[1];
sx q[1];
rz(-0.53002702) q[1];
sx q[1];
rz(0.0040667314) q[1];
rz(-pi) q[2];
rz(2.5745886) q[3];
sx q[3];
rz(-2.0006913) q[3];
sx q[3];
rz(-0.56741949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.064726) q[2];
sx q[2];
rz(-1.6419502) q[2];
sx q[2];
rz(-2.1323252) q[2];
rz(2.3276954) q[3];
sx q[3];
rz(-2.8129357) q[3];
sx q[3];
rz(0.3391268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086804248) q[0];
sx q[0];
rz(-1.4084933) q[0];
sx q[0];
rz(-2.8993697) q[0];
rz(-1.9107266) q[1];
sx q[1];
rz(-2.5098398) q[1];
sx q[1];
rz(-0.79808527) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0882814) q[0];
sx q[0];
rz(-1.6836731) q[0];
sx q[0];
rz(0.61932024) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2276137) q[2];
sx q[2];
rz(-1.7084885) q[2];
sx q[2];
rz(-2.340449) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4520769) q[1];
sx q[1];
rz(-2.544446) q[1];
sx q[1];
rz(2.4680016) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7834218) q[3];
sx q[3];
rz(-1.4205975) q[3];
sx q[3];
rz(1.9615441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.65608394) q[2];
sx q[2];
rz(-2.0723497) q[2];
sx q[2];
rz(-0.83267027) q[2];
rz(2.5101856) q[3];
sx q[3];
rz(-0.74149817) q[3];
sx q[3];
rz(3.1413063) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0176508) q[0];
sx q[0];
rz(-0.90921062) q[0];
sx q[0];
rz(2.9538474) q[0];
rz(-0.84367696) q[1];
sx q[1];
rz(-1.8673106) q[1];
sx q[1];
rz(-2.5086596) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.303514) q[0];
sx q[0];
rz(-1.9331341) q[0];
sx q[0];
rz(3.0658609) q[0];
rz(-0.68667163) q[2];
sx q[2];
rz(-1.8455659) q[2];
sx q[2];
rz(0.49897721) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5714076) q[1];
sx q[1];
rz(-1.8708036) q[1];
sx q[1];
rz(0.5350929) q[1];
x q[2];
rz(-1.2253148) q[3];
sx q[3];
rz(-0.36697432) q[3];
sx q[3];
rz(-2.2030694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8741499) q[2];
sx q[2];
rz(-0.94620693) q[2];
sx q[2];
rz(1.2713185) q[2];
rz(1.4960131) q[3];
sx q[3];
rz(-1.6451903) q[3];
sx q[3];
rz(-0.97150826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2520168) q[0];
sx q[0];
rz(-2.3887964) q[0];
sx q[0];
rz(2.4821607) q[0];
rz(-1.3176428) q[1];
sx q[1];
rz(-1.3392071) q[1];
sx q[1];
rz(2.3514294) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9371344) q[0];
sx q[0];
rz(-0.35872981) q[0];
sx q[0];
rz(2.782194) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1151562) q[2];
sx q[2];
rz(-0.038264878) q[2];
sx q[2];
rz(0.9963893) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.15406628) q[1];
sx q[1];
rz(-0.48168698) q[1];
sx q[1];
rz(0.14194004) q[1];
rz(-pi) q[2];
rz(-1.2983367) q[3];
sx q[3];
rz(-1.0333697) q[3];
sx q[3];
rz(-2.7481273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.059375199) q[2];
sx q[2];
rz(-1.6148115) q[2];
sx q[2];
rz(0.29602948) q[2];
rz(-2.5791903) q[3];
sx q[3];
rz(-0.86807576) q[3];
sx q[3];
rz(3.0276022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.247308) q[0];
sx q[0];
rz(-1.1881275) q[0];
sx q[0];
rz(-2.9344015) q[0];
rz(-2.1098792) q[1];
sx q[1];
rz(-1.1528015) q[1];
sx q[1];
rz(-1.4854887) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0061917) q[0];
sx q[0];
rz(-1.1814607) q[0];
sx q[0];
rz(0.72442309) q[0];
rz(-1.1274952) q[2];
sx q[2];
rz(-2.5336899) q[2];
sx q[2];
rz(-3.0435002) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.19073701) q[1];
sx q[1];
rz(-1.9190291) q[1];
sx q[1];
rz(1.7067043) q[1];
rz(-2.3014004) q[3];
sx q[3];
rz(-1.6512135) q[3];
sx q[3];
rz(1.059746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.038387211) q[2];
sx q[2];
rz(-0.81192553) q[2];
sx q[2];
rz(-2.0337598) q[2];
rz(0.92249089) q[3];
sx q[3];
rz(-1.731571) q[3];
sx q[3];
rz(-0.24423519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6237727) q[0];
sx q[0];
rz(-0.51384846) q[0];
sx q[0];
rz(-0.22860953) q[0];
rz(-0.27433431) q[1];
sx q[1];
rz(-2.3286596) q[1];
sx q[1];
rz(2.6913604) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73295702) q[0];
sx q[0];
rz(-1.8224323) q[0];
sx q[0];
rz(-2.3711414) q[0];
x q[1];
rz(-3.1338056) q[2];
sx q[2];
rz(-1.3415563) q[2];
sx q[2];
rz(-2.0321329) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6630312) q[1];
sx q[1];
rz(-1.8218464) q[1];
sx q[1];
rz(0.51086225) q[1];
x q[2];
rz(1.1563017) q[3];
sx q[3];
rz(-1.8179781) q[3];
sx q[3];
rz(0.37581623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.36394128) q[2];
sx q[2];
rz(-2.5886017) q[2];
sx q[2];
rz(-3.1138163) q[2];
rz(-0.76644301) q[3];
sx q[3];
rz(-1.981363) q[3];
sx q[3];
rz(-0.69507039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22208333) q[0];
sx q[0];
rz(-1.5064025) q[0];
sx q[0];
rz(-2.5191504) q[0];
rz(2.4355603) q[1];
sx q[1];
rz(-2.249554) q[1];
sx q[1];
rz(0.58696729) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3764899) q[0];
sx q[0];
rz(-0.73552948) q[0];
sx q[0];
rz(-2.2549948) q[0];
x q[1];
rz(-0.47759159) q[2];
sx q[2];
rz(-0.53672635) q[2];
sx q[2];
rz(-2.0411185) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.40499726) q[1];
sx q[1];
rz(-0.51144236) q[1];
sx q[1];
rz(-2.9963125) q[1];
rz(-3.1312084) q[3];
sx q[3];
rz(-2.2828013) q[3];
sx q[3];
rz(0.66085789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.075835) q[2];
sx q[2];
rz(-1.4328052) q[2];
sx q[2];
rz(-2.5471121) q[2];
rz(0.25932702) q[3];
sx q[3];
rz(-1.8766873) q[3];
sx q[3];
rz(-1.2172786) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038430564) q[0];
sx q[0];
rz(-0.90534627) q[0];
sx q[0];
rz(0.578798) q[0];
rz(-1.3102866) q[1];
sx q[1];
rz(-1.2625887) q[1];
sx q[1];
rz(-0.8126747) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7395514) q[0];
sx q[0];
rz(-1.427141) q[0];
sx q[0];
rz(2.6584635) q[0];
x q[1];
rz(1.1720288) q[2];
sx q[2];
rz(-2.5426455) q[2];
sx q[2];
rz(2.2676165) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1393237) q[1];
sx q[1];
rz(-1.9319169) q[1];
sx q[1];
rz(-2.4206672) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2872386) q[3];
sx q[3];
rz(-2.4378447) q[3];
sx q[3];
rz(0.077251369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7909214) q[2];
sx q[2];
rz(-1.2988043) q[2];
sx q[2];
rz(-2.217963) q[2];
rz(2.1212497) q[3];
sx q[3];
rz(-2.2444921) q[3];
sx q[3];
rz(-3.0237107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0608805) q[0];
sx q[0];
rz(-1.7055644) q[0];
sx q[0];
rz(-1.6545779) q[0];
rz(-0.05038536) q[1];
sx q[1];
rz(-2.3245508) q[1];
sx q[1];
rz(-2.494716) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9084307) q[0];
sx q[0];
rz(-1.1457503) q[0];
sx q[0];
rz(0.90954291) q[0];
x q[1];
rz(2.0503309) q[2];
sx q[2];
rz(-2.1796436) q[2];
sx q[2];
rz(-3.0225168) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9641831) q[1];
sx q[1];
rz(-0.46176592) q[1];
sx q[1];
rz(0.48514556) q[1];
rz(-0.96109747) q[3];
sx q[3];
rz(-0.25843474) q[3];
sx q[3];
rz(-2.1375382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6369624) q[2];
sx q[2];
rz(-1.5865822) q[2];
sx q[2];
rz(2.385425) q[2];
rz(1.470083) q[3];
sx q[3];
rz(-2.3556637) q[3];
sx q[3];
rz(-0.80529958) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.037755448) q[0];
sx q[0];
rz(-1.2498195) q[0];
sx q[0];
rz(-2.0334429) q[0];
rz(-2.8773819) q[1];
sx q[1];
rz(-1.4690396) q[1];
sx q[1];
rz(0.60233751) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2767576) q[0];
sx q[0];
rz(-1.2981155) q[0];
sx q[0];
rz(0.14436594) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9649023) q[2];
sx q[2];
rz(-1.4193221) q[2];
sx q[2];
rz(-3.1164411) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8089167) q[1];
sx q[1];
rz(-2.4530468) q[1];
sx q[1];
rz(1.6469514) q[1];
rz(-2.4906795) q[3];
sx q[3];
rz(-1.7554566) q[3];
sx q[3];
rz(1.0024651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5084874) q[2];
sx q[2];
rz(-0.54163951) q[2];
sx q[2];
rz(-0.87149054) q[2];
rz(-1.9320711) q[3];
sx q[3];
rz(-0.28035527) q[3];
sx q[3];
rz(0.59396321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7021983) q[0];
sx q[0];
rz(-1.7560503) q[0];
sx q[0];
rz(-0.51881292) q[0];
rz(2.5595472) q[1];
sx q[1];
rz(-1.1036292) q[1];
sx q[1];
rz(-1.6694952) q[1];
rz(2.4051278) q[2];
sx q[2];
rz(-1.8546551) q[2];
sx q[2];
rz(-1.3649884) q[2];
rz(1.4006079) q[3];
sx q[3];
rz(-0.88930973) q[3];
sx q[3];
rz(1.1607778) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
