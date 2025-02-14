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
rz(0.49864545) q[0];
sx q[0];
rz(-2.5957624) q[0];
sx q[0];
rz(-0.11830615) q[0];
rz(0.91953295) q[1];
sx q[1];
rz(-1.8383205) q[1];
sx q[1];
rz(-1.2143171) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3312553) q[0];
sx q[0];
rz(-2.6542568) q[0];
sx q[0];
rz(-1.1054071) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5164175) q[2];
sx q[2];
rz(-2.0707651) q[2];
sx q[2];
rz(-1.2631877) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8439629) q[1];
sx q[1];
rz(-1.4342035) q[1];
sx q[1];
rz(2.1928487) q[1];
rz(-3.1200527) q[3];
sx q[3];
rz(-1.2976754) q[3];
sx q[3];
rz(-0.5601976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.97293568) q[2];
sx q[2];
rz(-1.1309705) q[2];
sx q[2];
rz(1.0966148) q[2];
rz(0.24122572) q[3];
sx q[3];
rz(-1.3696407) q[3];
sx q[3];
rz(2.296804) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.089652561) q[0];
sx q[0];
rz(-2.0515428) q[0];
sx q[0];
rz(-2.7570778) q[0];
rz(0.3081201) q[1];
sx q[1];
rz(-2.6607951) q[1];
sx q[1];
rz(0.62517977) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85672545) q[0];
sx q[0];
rz(-1.856904) q[0];
sx q[0];
rz(1.5377783) q[0];
rz(1.4433631) q[2];
sx q[2];
rz(-2.2837167) q[2];
sx q[2];
rz(0.29297239) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5601666) q[1];
sx q[1];
rz(-2.4871965) q[1];
sx q[1];
rz(3.0017716) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6513807) q[3];
sx q[3];
rz(-1.0102838) q[3];
sx q[3];
rz(-1.1309159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.008931) q[2];
sx q[2];
rz(-0.93175685) q[2];
sx q[2];
rz(0.30161944) q[2];
rz(-0.53145069) q[3];
sx q[3];
rz(-2.2558432) q[3];
sx q[3];
rz(-2.0187995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6279491) q[0];
sx q[0];
rz(-1.0114089) q[0];
sx q[0];
rz(1.1847786) q[0];
rz(2.9966677) q[1];
sx q[1];
rz(-1.2624319) q[1];
sx q[1];
rz(-1.5473993) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74565174) q[0];
sx q[0];
rz(-2.3303623) q[0];
sx q[0];
rz(0.70180362) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9504531) q[2];
sx q[2];
rz(-1.7413503) q[2];
sx q[2];
rz(-1.4100687) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.48563924) q[1];
sx q[1];
rz(-1.1879293) q[1];
sx q[1];
rz(-0.98321557) q[1];
x q[2];
rz(-0.66158847) q[3];
sx q[3];
rz(-2.6433655) q[3];
sx q[3];
rz(-3.1116886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3323815) q[2];
sx q[2];
rz(-0.96497649) q[2];
sx q[2];
rz(-1.8154209) q[2];
rz(-2.2281846) q[3];
sx q[3];
rz(-1.2181506) q[3];
sx q[3];
rz(-2.6083045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42295414) q[0];
sx q[0];
rz(-2.664743) q[0];
sx q[0];
rz(1.3612716) q[0];
rz(-0.71296972) q[1];
sx q[1];
rz(-1.1402592) q[1];
sx q[1];
rz(2.4580809) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8806802) q[0];
sx q[0];
rz(-2.144038) q[0];
sx q[0];
rz(3.1295958) q[0];
rz(-1.3472802) q[2];
sx q[2];
rz(-1.3663251) q[2];
sx q[2];
rz(2.7140537) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.57505137) q[1];
sx q[1];
rz(-1.05541) q[1];
sx q[1];
rz(0.75830663) q[1];
x q[2];
rz(-0.24432791) q[3];
sx q[3];
rz(-1.4736946) q[3];
sx q[3];
rz(-1.8116784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5061364) q[2];
sx q[2];
rz(-2.7244302) q[2];
sx q[2];
rz(-0.48640856) q[2];
rz(2.6247315) q[3];
sx q[3];
rz(-1.9603399) q[3];
sx q[3];
rz(2.1080871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.085623398) q[0];
sx q[0];
rz(-1.6317246) q[0];
sx q[0];
rz(-1.7686718) q[0];
rz(-1.4068475) q[1];
sx q[1];
rz(-2.5739248) q[1];
sx q[1];
rz(-0.47703201) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1005362) q[0];
sx q[0];
rz(-1.1541799) q[0];
sx q[0];
rz(-0.21505298) q[0];
x q[1];
rz(2.1978756) q[2];
sx q[2];
rz(-2.5260128) q[2];
sx q[2];
rz(-2.1114388) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7030695) q[1];
sx q[1];
rz(-1.0756936) q[1];
sx q[1];
rz(0.34214051) q[1];
rz(-pi) q[2];
rz(2.2955016) q[3];
sx q[3];
rz(-2.2109004) q[3];
sx q[3];
rz(0.37615955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.10739747) q[2];
sx q[2];
rz(-0.98661462) q[2];
sx q[2];
rz(2.6743215) q[2];
rz(2.9723736) q[3];
sx q[3];
rz(-1.7305814) q[3];
sx q[3];
rz(-0.30538487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6071103) q[0];
sx q[0];
rz(-0.87962532) q[0];
sx q[0];
rz(0.18950263) q[0];
rz(2.8684008) q[1];
sx q[1];
rz(-2.0676282) q[1];
sx q[1];
rz(2.3409519) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9439745) q[0];
sx q[0];
rz(-1.6139784) q[0];
sx q[0];
rz(-1.3373242) q[0];
x q[1];
rz(0.34363644) q[2];
sx q[2];
rz(-0.26310194) q[2];
sx q[2];
rz(-0.084089078) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.26594083) q[1];
sx q[1];
rz(-1.7663301) q[1];
sx q[1];
rz(3.0102848) q[1];
x q[2];
rz(-2.4005049) q[3];
sx q[3];
rz(-2.4588278) q[3];
sx q[3];
rz(-2.6094451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2022986) q[2];
sx q[2];
rz(-0.17387986) q[2];
sx q[2];
rz(2.6046216) q[2];
rz(-1.8592853) q[3];
sx q[3];
rz(-1.3064281) q[3];
sx q[3];
rz(-1.6622701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.031161664) q[0];
sx q[0];
rz(-2.3261676) q[0];
sx q[0];
rz(-1.8102616) q[0];
rz(1.7000465) q[1];
sx q[1];
rz(-1.6621637) q[1];
sx q[1];
rz(-1.9190681) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5957907) q[0];
sx q[0];
rz(-1.8078601) q[0];
sx q[0];
rz(1.14181) q[0];
x q[1];
rz(0.55180727) q[2];
sx q[2];
rz(-1.5071609) q[2];
sx q[2];
rz(0.76630521) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.26204208) q[1];
sx q[1];
rz(-2.4051504) q[1];
sx q[1];
rz(1.2465191) q[1];
rz(-pi) q[2];
rz(0.11925536) q[3];
sx q[3];
rz(-1.1290765) q[3];
sx q[3];
rz(0.91631232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.20087251) q[2];
sx q[2];
rz(-2.6437289) q[2];
sx q[2];
rz(0.3328003) q[2];
rz(-2.8424272) q[3];
sx q[3];
rz(-1.9009512) q[3];
sx q[3];
rz(-0.021473264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12896319) q[0];
sx q[0];
rz(-0.7681995) q[0];
sx q[0];
rz(-0.27498883) q[0];
rz(-3.0601652) q[1];
sx q[1];
rz(-2.3402201) q[1];
sx q[1];
rz(-1.3224695) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0819495) q[0];
sx q[0];
rz(-2.9026051) q[0];
sx q[0];
rz(-1.3168524) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7733388) q[2];
sx q[2];
rz(-2.2241631) q[2];
sx q[2];
rz(2.7453932) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7863955) q[1];
sx q[1];
rz(-1.2010368) q[1];
sx q[1];
rz(-0.74798297) q[1];
x q[2];
rz(2.4149618) q[3];
sx q[3];
rz(-1.0612349) q[3];
sx q[3];
rz(-1.0973615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4203809) q[2];
sx q[2];
rz(-0.12004852) q[2];
sx q[2];
rz(1.4018641) q[2];
rz(1.2841355) q[3];
sx q[3];
rz(-1.1893585) q[3];
sx q[3];
rz(3.106626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5597124) q[0];
sx q[0];
rz(-1.5927097) q[0];
sx q[0];
rz(2.6522719) q[0];
rz(-0.020546546) q[1];
sx q[1];
rz(-1.2362044) q[1];
sx q[1];
rz(0.70125854) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1398292) q[0];
sx q[0];
rz(-1.0047303) q[0];
sx q[0];
rz(-2.8553814) q[0];
x q[1];
rz(1.6270774) q[2];
sx q[2];
rz(-2.4372134) q[2];
sx q[2];
rz(3.0075108) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.78811443) q[1];
sx q[1];
rz(-0.81190049) q[1];
sx q[1];
rz(-0.32737704) q[1];
rz(-pi) q[2];
rz(0.55072983) q[3];
sx q[3];
rz(-2.3756873) q[3];
sx q[3];
rz(-2.285241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7265085) q[2];
sx q[2];
rz(-2.5231611) q[2];
sx q[2];
rz(2.6611967) q[2];
rz(-1.4091617) q[3];
sx q[3];
rz(-1.3536072) q[3];
sx q[3];
rz(0.33717808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3163863) q[0];
sx q[0];
rz(-1.9012863) q[0];
sx q[0];
rz(-0.41148841) q[0];
rz(1.8972137) q[1];
sx q[1];
rz(-2.0020516) q[1];
sx q[1];
rz(-0.32434514) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.024687) q[0];
sx q[0];
rz(-1.9431264) q[0];
sx q[0];
rz(-3.0850379) q[0];
x q[1];
rz(2.3408808) q[2];
sx q[2];
rz(-2.6691648) q[2];
sx q[2];
rz(-2.1683482) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.56411769) q[1];
sx q[1];
rz(-2.1337957) q[1];
sx q[1];
rz(3.0433118) q[1];
x q[2];
rz(-1.5839229) q[3];
sx q[3];
rz(-1.6919976) q[3];
sx q[3];
rz(-1.9399724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9559429) q[2];
sx q[2];
rz(-2.0273384) q[2];
sx q[2];
rz(-3.0246217) q[2];
rz(-2.1864435) q[3];
sx q[3];
rz(-0.17894608) q[3];
sx q[3];
rz(-0.54168934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2627926) q[0];
sx q[0];
rz(-2.265082) q[0];
sx q[0];
rz(-2.2646917) q[0];
rz(2.0211438) q[1];
sx q[1];
rz(-2.3255377) q[1];
sx q[1];
rz(-1.9768523) q[1];
rz(-1.9399613) q[2];
sx q[2];
rz(-0.537048) q[2];
sx q[2];
rz(-0.67759003) q[2];
rz(0.025259947) q[3];
sx q[3];
rz(-1.8413481) q[3];
sx q[3];
rz(-3.1006364) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
