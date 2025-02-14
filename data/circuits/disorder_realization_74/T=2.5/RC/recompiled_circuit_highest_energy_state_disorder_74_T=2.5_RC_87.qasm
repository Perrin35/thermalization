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
rz(-2.3251301) q[0];
sx q[0];
rz(-0.10179585) q[0];
sx q[0];
rz(-0.53959227) q[0];
rz(0.49630961) q[1];
sx q[1];
rz(-0.30975431) q[1];
sx q[1];
rz(-2.6024979) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2236299) q[0];
sx q[0];
rz(-1.5453891) q[0];
sx q[0];
rz(2.7001803) q[0];
rz(-pi) q[1];
rz(-1.4194054) q[2];
sx q[2];
rz(-1.2149723) q[2];
sx q[2];
rz(-0.21499888) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4161168) q[1];
sx q[1];
rz(-1.892964) q[1];
sx q[1];
rz(2.9941818) q[1];
rz(-pi) q[2];
rz(1.6786537) q[3];
sx q[3];
rz(-0.65179208) q[3];
sx q[3];
rz(0.45676132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.76217905) q[2];
sx q[2];
rz(-0.87612408) q[2];
sx q[2];
rz(1.1890821) q[2];
rz(1.9573697) q[3];
sx q[3];
rz(-2.2370179) q[3];
sx q[3];
rz(-1.5466461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14520833) q[0];
sx q[0];
rz(-1.5897911) q[0];
sx q[0];
rz(-1.8183964) q[0];
rz(2.6595751) q[1];
sx q[1];
rz(-0.91164416) q[1];
sx q[1];
rz(0.97420305) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27069399) q[0];
sx q[0];
rz(-2.2091731) q[0];
sx q[0];
rz(-2.9186072) q[0];
x q[1];
rz(-1.4575794) q[2];
sx q[2];
rz(-1.4721451) q[2];
sx q[2];
rz(1.9960008) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5017101) q[1];
sx q[1];
rz(-0.5941662) q[1];
sx q[1];
rz(-2.7133248) q[1];
rz(2.5752546) q[3];
sx q[3];
rz(-1.7741331) q[3];
sx q[3];
rz(-1.4161863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1366068) q[2];
sx q[2];
rz(-2.5544781) q[2];
sx q[2];
rz(0.74964398) q[2];
rz(-0.43198112) q[3];
sx q[3];
rz(-2.0260729) q[3];
sx q[3];
rz(1.8384793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5169446) q[0];
sx q[0];
rz(-2.8693146) q[0];
sx q[0];
rz(-2.5352449) q[0];
rz(-1.3336522) q[1];
sx q[1];
rz(-1.9275815) q[1];
sx q[1];
rz(-0.25951728) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57588803) q[0];
sx q[0];
rz(-1.7262926) q[0];
sx q[0];
rz(-1.8207654) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3679996) q[2];
sx q[2];
rz(-1.742387) q[2];
sx q[2];
rz(-0.20538782) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.77713457) q[1];
sx q[1];
rz(-1.196047) q[1];
sx q[1];
rz(-1.115429) q[1];
rz(1.0084148) q[3];
sx q[3];
rz(-1.9451688) q[3];
sx q[3];
rz(-1.238747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2567265) q[2];
sx q[2];
rz(-1.7776411) q[2];
sx q[2];
rz(1.5474896) q[2];
rz(0.78440845) q[3];
sx q[3];
rz(-1.6268077) q[3];
sx q[3];
rz(-1.5756395) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6467658) q[0];
sx q[0];
rz(-2.0443199) q[0];
sx q[0];
rz(-2.8821017) q[0];
rz(-1.9208113) q[1];
sx q[1];
rz(-2.6916598) q[1];
sx q[1];
rz(-1.6993914) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7176162) q[0];
sx q[0];
rz(-2.0956796) q[0];
sx q[0];
rz(-0.54846008) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0109826) q[2];
sx q[2];
rz(-0.97960661) q[2];
sx q[2];
rz(-1.755205) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4220548) q[1];
sx q[1];
rz(-1.9214848) q[1];
sx q[1];
rz(1.8338051) q[1];
rz(-pi) q[2];
rz(-3.0715748) q[3];
sx q[3];
rz(-0.88980955) q[3];
sx q[3];
rz(0.79047608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0356902) q[2];
sx q[2];
rz(-1.8041958) q[2];
sx q[2];
rz(2.7316459) q[2];
rz(-1.9299054) q[3];
sx q[3];
rz(-2.4544921) q[3];
sx q[3];
rz(-1.3338026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7164417) q[0];
sx q[0];
rz(-2.4254159) q[0];
sx q[0];
rz(2.8675365) q[0];
rz(-0.43824276) q[1];
sx q[1];
rz(-1.848105) q[1];
sx q[1];
rz(2.4058707) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5925795) q[0];
sx q[0];
rz(-1.6578478) q[0];
sx q[0];
rz(-0.75496952) q[0];
rz(-pi) q[1];
rz(-0.14938905) q[2];
sx q[2];
rz(-1.427703) q[2];
sx q[2];
rz(-2.9516475) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.055701697) q[1];
sx q[1];
rz(-0.38905479) q[1];
sx q[1];
rz(0.86430734) q[1];
rz(0.93345668) q[3];
sx q[3];
rz(-0.75691716) q[3];
sx q[3];
rz(-2.8035823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.938544) q[2];
sx q[2];
rz(-1.4578578) q[2];
sx q[2];
rz(-2.5277444) q[2];
rz(-2.7326873) q[3];
sx q[3];
rz(-0.96858612) q[3];
sx q[3];
rz(-2.3992505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3530389) q[0];
sx q[0];
rz(-0.63816324) q[0];
sx q[0];
rz(-1.2370538) q[0];
rz(1.4452112) q[1];
sx q[1];
rz(-2.684869) q[1];
sx q[1];
rz(2.9663185) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0171368) q[0];
sx q[0];
rz(-1.5457898) q[0];
sx q[0];
rz(1.7437922) q[0];
rz(-pi) q[1];
rz(-0.66368033) q[2];
sx q[2];
rz(-2.366407) q[2];
sx q[2];
rz(2.7221687) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8126909) q[1];
sx q[1];
rz(-2.0022503) q[1];
sx q[1];
rz(-2.575909) q[1];
rz(-0.82476576) q[3];
sx q[3];
rz(-2.3769393) q[3];
sx q[3];
rz(-2.1869756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.96482977) q[2];
sx q[2];
rz(-1.8076597) q[2];
sx q[2];
rz(-2.1691587) q[2];
rz(1.6644299) q[3];
sx q[3];
rz(-1.1494613) q[3];
sx q[3];
rz(0.76178637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2870188) q[0];
sx q[0];
rz(-3.119097) q[0];
sx q[0];
rz(-0.35312411) q[0];
rz(-1.0278206) q[1];
sx q[1];
rz(-1.9408344) q[1];
sx q[1];
rz(-2.1305398) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38336223) q[0];
sx q[0];
rz(-0.35759059) q[0];
sx q[0];
rz(0.041186853) q[0];
rz(-0.89247668) q[2];
sx q[2];
rz(-0.8304285) q[2];
sx q[2];
rz(2.3811946) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8450635) q[1];
sx q[1];
rz(-1.4140714) q[1];
sx q[1];
rz(1.2178376) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8854463) q[3];
sx q[3];
rz(-1.8487747) q[3];
sx q[3];
rz(1.4643471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8280243) q[2];
sx q[2];
rz(-2.821533) q[2];
sx q[2];
rz(-1.774452) q[2];
rz(1.8732871) q[3];
sx q[3];
rz(-1.6627848) q[3];
sx q[3];
rz(0.84793276) q[3];
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
rz(1.1139514) q[0];
sx q[0];
rz(-2.5340762) q[0];
sx q[0];
rz(-1.129958) q[0];
rz(-0.21895151) q[1];
sx q[1];
rz(-1.7349225) q[1];
sx q[1];
rz(2.2311282) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8559056) q[0];
sx q[0];
rz(-0.98074161) q[0];
sx q[0];
rz(-1.1678215) q[0];
x q[1];
rz(-0.16436188) q[2];
sx q[2];
rz(-0.48381643) q[2];
sx q[2];
rz(-1.9828895) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1840802) q[1];
sx q[1];
rz(-1.0614479) q[1];
sx q[1];
rz(2.8345076) q[1];
rz(-0.72515709) q[3];
sx q[3];
rz(-1.2598383) q[3];
sx q[3];
rz(-2.6127315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.8255446) q[2];
sx q[2];
rz(-0.80222183) q[2];
sx q[2];
rz(0.82553378) q[2];
rz(0.46142203) q[3];
sx q[3];
rz(-1.2666707) q[3];
sx q[3];
rz(-2.6958444) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6398741) q[0];
sx q[0];
rz(-1.8032782) q[0];
sx q[0];
rz(-2.3097532) q[0];
rz(0.72744751) q[1];
sx q[1];
rz(-1.7214382) q[1];
sx q[1];
rz(-0.096253455) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1883586) q[0];
sx q[0];
rz(-1.112845) q[0];
sx q[0];
rz(-2.1403378) q[0];
rz(2.8873994) q[2];
sx q[2];
rz(-1.5113748) q[2];
sx q[2];
rz(0.78320247) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2029496) q[1];
sx q[1];
rz(-1.661282) q[1];
sx q[1];
rz(1.3590779) q[1];
x q[2];
rz(-2.9828137) q[3];
sx q[3];
rz(-1.3769866) q[3];
sx q[3];
rz(-0.69413444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.40621296) q[2];
sx q[2];
rz(-1.7071743) q[2];
sx q[2];
rz(-1.5489138) q[2];
rz(-2.5088572) q[3];
sx q[3];
rz(-1.9097208) q[3];
sx q[3];
rz(0.24937853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(11/(9*pi)) q[0];
sx q[0];
rz(-2.1010375) q[0];
sx q[0];
rz(2.7885875) q[0];
rz(-0.32265916) q[1];
sx q[1];
rz(-2.8162075) q[1];
sx q[1];
rz(-0.37193146) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8648099) q[0];
sx q[0];
rz(-1.7686592) q[0];
sx q[0];
rz(-1.2643306) q[0];
rz(2.811889) q[2];
sx q[2];
rz(-1.4011382) q[2];
sx q[2];
rz(-0.28126954) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.10807762) q[1];
sx q[1];
rz(-1.4107499) q[1];
sx q[1];
rz(-2.8041502) q[1];
rz(-pi) q[2];
rz(3.102469) q[3];
sx q[3];
rz(-0.84661603) q[3];
sx q[3];
rz(0.20871221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.17452621) q[2];
sx q[2];
rz(-1.1897949) q[2];
sx q[2];
rz(-1.2971499) q[2];
rz(0.8477115) q[3];
sx q[3];
rz(-1.984237) q[3];
sx q[3];
rz(-2.1367836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2321155) q[0];
sx q[0];
rz(-1.7625325) q[0];
sx q[0];
rz(0.43176227) q[0];
rz(-1.3879981) q[1];
sx q[1];
rz(-1.2139865) q[1];
sx q[1];
rz(-0.075275631) q[1];
rz(0.83609381) q[2];
sx q[2];
rz(-2.2420364) q[2];
sx q[2];
rz(-1.4902761) q[2];
rz(0.79741464) q[3];
sx q[3];
rz(-1.8731464) q[3];
sx q[3];
rz(-2.1220589) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
