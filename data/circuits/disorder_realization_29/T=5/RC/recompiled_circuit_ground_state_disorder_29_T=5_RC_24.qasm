OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.090341181) q[0];
sx q[0];
rz(2.0988965) q[0];
sx q[0];
rz(12.287696) q[0];
rz(-1.9947808) q[1];
sx q[1];
rz(-2.6169701) q[1];
sx q[1];
rz(1.2710748) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81084736) q[0];
sx q[0];
rz(-2.8804104) q[0];
sx q[0];
rz(0.07075858) q[0];
x q[1];
rz(2.4248883) q[2];
sx q[2];
rz(-0.19858805) q[2];
sx q[2];
rz(-2.8633397) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6587352) q[1];
sx q[1];
rz(-0.71719155) q[1];
sx q[1];
rz(-1.6916005) q[1];
rz(-pi) q[2];
rz(0.75436169) q[3];
sx q[3];
rz(-2.3791557) q[3];
sx q[3];
rz(2.0319394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.110454) q[2];
sx q[2];
rz(-1.9240856) q[2];
sx q[2];
rz(-0.26220599) q[2];
rz(-1.0555222) q[3];
sx q[3];
rz(-1.2473829) q[3];
sx q[3];
rz(-3.054255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6156886) q[0];
sx q[0];
rz(-0.50351244) q[0];
sx q[0];
rz(2.9600034) q[0];
rz(-1.7407821) q[1];
sx q[1];
rz(-1.9488275) q[1];
sx q[1];
rz(0.84091944) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070905237) q[0];
sx q[0];
rz(-1.5531666) q[0];
sx q[0];
rz(-2.8215462) q[0];
rz(-0.96581755) q[2];
sx q[2];
rz(-0.73558319) q[2];
sx q[2];
rz(-3.0515665) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7279198) q[1];
sx q[1];
rz(-2.9693954) q[1];
sx q[1];
rz(1.3581469) q[1];
rz(-pi) q[2];
rz(2.0729468) q[3];
sx q[3];
rz(-2.4663894) q[3];
sx q[3];
rz(-1.5218376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38222868) q[2];
sx q[2];
rz(-2.9266734) q[2];
sx q[2];
rz(2.0767029) q[2];
rz(-3.0801638) q[3];
sx q[3];
rz(-1.3353142) q[3];
sx q[3];
rz(1.6399062) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29838022) q[0];
sx q[0];
rz(-1.1792264) q[0];
sx q[0];
rz(0.19905736) q[0];
rz(0.82636034) q[1];
sx q[1];
rz(-0.91824707) q[1];
sx q[1];
rz(-2.3645649) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13764405) q[0];
sx q[0];
rz(-0.17034082) q[0];
sx q[0];
rz(1.9714703) q[0];
x q[1];
rz(2.1146745) q[2];
sx q[2];
rz(-1.9444398) q[2];
sx q[2];
rz(2.0622562) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9553634) q[1];
sx q[1];
rz(-1.0252153) q[1];
sx q[1];
rz(-2.8377691) q[1];
rz(-pi) q[2];
rz(-0.2774419) q[3];
sx q[3];
rz(-0.99511569) q[3];
sx q[3];
rz(2.4572008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0394773) q[2];
sx q[2];
rz(-1.9934883) q[2];
sx q[2];
rz(-2.6285505) q[2];
rz(-2.3593864) q[3];
sx q[3];
rz(-1.9358044) q[3];
sx q[3];
rz(1.3768844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9272598) q[0];
sx q[0];
rz(-1.5434649) q[0];
sx q[0];
rz(-1.4904892) q[0];
rz(2.6093042) q[1];
sx q[1];
rz(-0.63096255) q[1];
sx q[1];
rz(-1.5375563) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2435115) q[0];
sx q[0];
rz(-1.6771206) q[0];
sx q[0];
rz(2.9778019) q[0];
rz(0.093393383) q[2];
sx q[2];
rz(-1.0265145) q[2];
sx q[2];
rz(-1.613008) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.58504471) q[1];
sx q[1];
rz(-0.93657199) q[1];
sx q[1];
rz(-0.08794959) q[1];
rz(0.52738366) q[3];
sx q[3];
rz(-1.77447) q[3];
sx q[3];
rz(-2.7319997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6060467) q[2];
sx q[2];
rz(-0.54680768) q[2];
sx q[2];
rz(3.1287076) q[2];
rz(-1.2031201) q[3];
sx q[3];
rz(-1.4667526) q[3];
sx q[3];
rz(1.0217246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0032229) q[0];
sx q[0];
rz(-0.62507451) q[0];
sx q[0];
rz(-1.7417997) q[0];
rz(2.7895582) q[1];
sx q[1];
rz(-1.0540009) q[1];
sx q[1];
rz(-2.3037516) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.261917) q[0];
sx q[0];
rz(-0.81955541) q[0];
sx q[0];
rz(-2.9330334) q[0];
x q[1];
rz(2.3532392) q[2];
sx q[2];
rz(-0.78561312) q[2];
sx q[2];
rz(-2.2503302) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6988236) q[1];
sx q[1];
rz(-1.6976446) q[1];
sx q[1];
rz(2.5731099) q[1];
x q[2];
rz(-2.3901108) q[3];
sx q[3];
rz(-2.9183636) q[3];
sx q[3];
rz(0.022836784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1167404) q[2];
sx q[2];
rz(-0.96472538) q[2];
sx q[2];
rz(-1.9579197) q[2];
rz(0.27070326) q[3];
sx q[3];
rz(-2.201122) q[3];
sx q[3];
rz(-1.6089599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.026641332) q[0];
sx q[0];
rz(-2.5178435) q[0];
sx q[0];
rz(-2.5624516) q[0];
rz(1.9193513) q[1];
sx q[1];
rz(-1.8341583) q[1];
sx q[1];
rz(-1.1096257) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1117119) q[0];
sx q[0];
rz(-1.7222026) q[0];
sx q[0];
rz(0.16686186) q[0];
rz(-0.17758421) q[2];
sx q[2];
rz(-2.8320667) q[2];
sx q[2];
rz(0.9847275) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.077495726) q[1];
sx q[1];
rz(-1.9019097) q[1];
sx q[1];
rz(1.8331794) q[1];
rz(-2.6251276) q[3];
sx q[3];
rz(-1.0871917) q[3];
sx q[3];
rz(0.91396871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9455202) q[2];
sx q[2];
rz(-1.2976982) q[2];
sx q[2];
rz(-2.784139) q[2];
rz(1.9516021) q[3];
sx q[3];
rz(-1.9414732) q[3];
sx q[3];
rz(-1.3397217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5424407) q[0];
sx q[0];
rz(-1.0170794) q[0];
sx q[0];
rz(-0.79208148) q[0];
rz(0.7041086) q[1];
sx q[1];
rz(-0.45583615) q[1];
sx q[1];
rz(1.5584996) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7501935) q[0];
sx q[0];
rz(-2.0349841) q[0];
sx q[0];
rz(1.5740921) q[0];
rz(-pi) q[1];
rz(-2.547491) q[2];
sx q[2];
rz(-2.1081446) q[2];
sx q[2];
rz(-2.0955476) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5643107) q[1];
sx q[1];
rz(-2.7140712) q[1];
sx q[1];
rz(3.0583818) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3901863) q[3];
sx q[3];
rz(-1.7370686) q[3];
sx q[3];
rz(0.64054447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.11129561) q[2];
sx q[2];
rz(-2.3534677) q[2];
sx q[2];
rz(-0.0041858717) q[2];
rz(-0.26675102) q[3];
sx q[3];
rz(-1.6175852) q[3];
sx q[3];
rz(2.4193616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5608212) q[0];
sx q[0];
rz(-0.72951356) q[0];
sx q[0];
rz(0.15604493) q[0];
rz(1.5566298) q[1];
sx q[1];
rz(-0.86235756) q[1];
sx q[1];
rz(-2.5468266) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1536008) q[0];
sx q[0];
rz(-2.5770441) q[0];
sx q[0];
rz(2.6507244) q[0];
x q[1];
rz(-1.2807219) q[2];
sx q[2];
rz(-2.2422487) q[2];
sx q[2];
rz(0.69556505) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8957828) q[1];
sx q[1];
rz(-1.0059662) q[1];
sx q[1];
rz(-0.71706949) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.90601633) q[3];
sx q[3];
rz(-0.95905868) q[3];
sx q[3];
rz(-2.8869528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.56738371) q[2];
sx q[2];
rz(-1.5790066) q[2];
sx q[2];
rz(1.0279921) q[2];
rz(1.7993641) q[3];
sx q[3];
rz(-2.4864311) q[3];
sx q[3];
rz(1.8067693) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8780355) q[0];
sx q[0];
rz(-1.6850543) q[0];
sx q[0];
rz(0.79826075) q[0];
rz(2.3410666) q[1];
sx q[1];
rz(-0.48897484) q[1];
sx q[1];
rz(0.54623234) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9437708) q[0];
sx q[0];
rz(-1.6795625) q[0];
sx q[0];
rz(-2.1740211) q[0];
x q[1];
rz(1.9990218) q[2];
sx q[2];
rz(-1.0248803) q[2];
sx q[2];
rz(2.2466618) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.45193397) q[1];
sx q[1];
rz(-1.0946737) q[1];
sx q[1];
rz(2.1228288) q[1];
rz(-2.8343938) q[3];
sx q[3];
rz(-2.9044378) q[3];
sx q[3];
rz(-3.0320771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.74786782) q[2];
sx q[2];
rz(-2.0515029) q[2];
sx q[2];
rz(0.063035034) q[2];
rz(-0.6978327) q[3];
sx q[3];
rz(-0.2581667) q[3];
sx q[3];
rz(0.91309083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6460655) q[0];
sx q[0];
rz(-2.263948) q[0];
sx q[0];
rz(-2.7255507) q[0];
rz(-1.7795732) q[1];
sx q[1];
rz(-2.0174618) q[1];
sx q[1];
rz(1.8488041) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36609367) q[0];
sx q[0];
rz(-1.3028569) q[0];
sx q[0];
rz(-1.8084722) q[0];
rz(-pi) q[1];
rz(-1.850239) q[2];
sx q[2];
rz(-2.0115174) q[2];
sx q[2];
rz(-0.27383495) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.67329183) q[1];
sx q[1];
rz(-2.3929962) q[1];
sx q[1];
rz(-2.2842201) q[1];
rz(-pi) q[2];
rz(-0.26117321) q[3];
sx q[3];
rz(-1.9800485) q[3];
sx q[3];
rz(0.60200426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.68299874) q[2];
sx q[2];
rz(-1.3364044) q[2];
sx q[2];
rz(2.2340753) q[2];
rz(1.8806184) q[3];
sx q[3];
rz(-0.57456273) q[3];
sx q[3];
rz(-1.6921836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4954729) q[0];
sx q[0];
rz(-1.6874122) q[0];
sx q[0];
rz(0.73943403) q[0];
rz(2.6932035) q[1];
sx q[1];
rz(-1.3403475) q[1];
sx q[1];
rz(1.1833804) q[1];
rz(2.7896055) q[2];
sx q[2];
rz(-1.864973) q[2];
sx q[2];
rz(-0.61054338) q[2];
rz(3.0082416) q[3];
sx q[3];
rz(-0.84682087) q[3];
sx q[3];
rz(-3.1414733) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
