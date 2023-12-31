OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.9632602) q[0];
sx q[0];
rz(-1.652521) q[0];
sx q[0];
rz(0.89515495) q[0];
rz(2.826638) q[1];
sx q[1];
rz(2.0575674) q[1];
sx q[1];
rz(10.881012) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2468949) q[0];
sx q[0];
rz(-3.0092735) q[0];
sx q[0];
rz(-1.4207065) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1158754) q[2];
sx q[2];
rz(-1.8946049) q[2];
sx q[2];
rz(0.88534249) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1617042) q[1];
sx q[1];
rz(-2.4474505) q[1];
sx q[1];
rz(-2.90467) q[1];
rz(-pi) q[2];
rz(-1.5159357) q[3];
sx q[3];
rz(-1.8901955) q[3];
sx q[3];
rz(1.1064135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.75498092) q[2];
sx q[2];
rz(-1.3957916) q[2];
sx q[2];
rz(-2.6888729) q[2];
rz(-2.9833941) q[3];
sx q[3];
rz(-2.4499564) q[3];
sx q[3];
rz(0.89481568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2125856) q[0];
sx q[0];
rz(-2.043262) q[0];
sx q[0];
rz(-1.989495) q[0];
rz(-1.903803) q[1];
sx q[1];
rz(-1.5367616) q[1];
sx q[1];
rz(-0.47098413) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1176227) q[0];
sx q[0];
rz(-1.515944) q[0];
sx q[0];
rz(2.0532002) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9085625) q[2];
sx q[2];
rz(-0.62440364) q[2];
sx q[2];
rz(1.0242467) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4091332) q[1];
sx q[1];
rz(-2.0583378) q[1];
sx q[1];
rz(-3.0797466) q[1];
rz(-pi) q[2];
rz(0.095951565) q[3];
sx q[3];
rz(-0.6978242) q[3];
sx q[3];
rz(0.56688353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.058078893) q[2];
sx q[2];
rz(-2.5364272) q[2];
sx q[2];
rz(-2.8857968) q[2];
rz(1.4852218) q[3];
sx q[3];
rz(-1.1896313) q[3];
sx q[3];
rz(2.9799057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26329041) q[0];
sx q[0];
rz(-2.0354164) q[0];
sx q[0];
rz(2.8702452) q[0];
rz(-2.4052606) q[1];
sx q[1];
rz(-1.5356179) q[1];
sx q[1];
rz(-2.7022865) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8205748) q[0];
sx q[0];
rz(-1.4903755) q[0];
sx q[0];
rz(-3.1119425) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13621026) q[2];
sx q[2];
rz(-2.7779707) q[2];
sx q[2];
rz(1.9759535) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15553741) q[1];
sx q[1];
rz(-1.7441161) q[1];
sx q[1];
rz(-0.9539414) q[1];
x q[2];
rz(-1.0147694) q[3];
sx q[3];
rz(-1.1960256) q[3];
sx q[3];
rz(2.8127363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1674041) q[2];
sx q[2];
rz(-1.5524652) q[2];
sx q[2];
rz(2.9411194) q[2];
rz(0.75508562) q[3];
sx q[3];
rz(-1.0198159) q[3];
sx q[3];
rz(0.42373207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.077483594) q[0];
sx q[0];
rz(-1.5492726) q[0];
sx q[0];
rz(1.8970998) q[0];
rz(2.3311133) q[1];
sx q[1];
rz(-1.3296209) q[1];
sx q[1];
rz(-2.2669852) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8133102) q[0];
sx q[0];
rz(-2.5173442) q[0];
sx q[0];
rz(1.3024131) q[0];
rz(1.4250408) q[2];
sx q[2];
rz(-1.1894023) q[2];
sx q[2];
rz(-2.8085453) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3687467) q[1];
sx q[1];
rz(-0.44279848) q[1];
sx q[1];
rz(2.5889791) q[1];
rz(-1.6455669) q[3];
sx q[3];
rz(-1.6755591) q[3];
sx q[3];
rz(1.2391702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0041634) q[2];
sx q[2];
rz(-0.88976294) q[2];
sx q[2];
rz(-1.4271663) q[2];
rz(0.066120474) q[3];
sx q[3];
rz(-2.7757006) q[3];
sx q[3];
rz(1.6920413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78113294) q[0];
sx q[0];
rz(-2.9678678) q[0];
sx q[0];
rz(2.5710035) q[0];
rz(-0.55496201) q[1];
sx q[1];
rz(-0.73736063) q[1];
sx q[1];
rz(-0.74329174) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0750908) q[0];
sx q[0];
rz(-2.4473303) q[0];
sx q[0];
rz(-1.2230722) q[0];
x q[1];
rz(1.4666124) q[2];
sx q[2];
rz(-2.0698692) q[2];
sx q[2];
rz(-0.68945976) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9651282) q[1];
sx q[1];
rz(-1.3511409) q[1];
sx q[1];
rz(-1.7996644) q[1];
rz(-pi) q[2];
rz(-0.80446135) q[3];
sx q[3];
rz(-1.3692229) q[3];
sx q[3];
rz(-1.7427012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9479998) q[2];
sx q[2];
rz(-1.7717382) q[2];
sx q[2];
rz(-0.26838475) q[2];
rz(1.0466446) q[3];
sx q[3];
rz(-0.36473754) q[3];
sx q[3];
rz(-0.31782761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5140117) q[0];
sx q[0];
rz(-1.528897) q[0];
sx q[0];
rz(0.50672379) q[0];
rz(-0.2535893) q[1];
sx q[1];
rz(-1.8702303) q[1];
sx q[1];
rz(-1.0553029) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2282288) q[0];
sx q[0];
rz(-2.1852487) q[0];
sx q[0];
rz(-2.9173304) q[0];
x q[1];
rz(-0.92743404) q[2];
sx q[2];
rz(-2.3775527) q[2];
sx q[2];
rz(1.6709136) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1736974) q[1];
sx q[1];
rz(-2.3579862) q[1];
sx q[1];
rz(2.150279) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.7430195) q[3];
sx q[3];
rz(-1.395874) q[3];
sx q[3];
rz(-2.7878441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.48173299) q[2];
sx q[2];
rz(-1.0446171) q[2];
sx q[2];
rz(0.8824904) q[2];
rz(-2.4957538) q[3];
sx q[3];
rz(-1.1487938) q[3];
sx q[3];
rz(-1.8576436) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6362474) q[0];
sx q[0];
rz(-1.212965) q[0];
sx q[0];
rz(0.77254599) q[0];
rz(1.4121217) q[1];
sx q[1];
rz(-1.2071143) q[1];
sx q[1];
rz(-2.5678182) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44874292) q[0];
sx q[0];
rz(-2.9208555) q[0];
sx q[0];
rz(1.0462532) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3705809) q[2];
sx q[2];
rz(-1.5924615) q[2];
sx q[2];
rz(1.3862762) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7635203) q[1];
sx q[1];
rz(-2.1280648) q[1];
sx q[1];
rz(2.4850363) q[1];
x q[2];
rz(-1.0049099) q[3];
sx q[3];
rz(-2.8661869) q[3];
sx q[3];
rz(-2.3435081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.039375719) q[2];
sx q[2];
rz(-2.6874459) q[2];
sx q[2];
rz(0.77073628) q[2];
rz(2.7052774) q[3];
sx q[3];
rz(-1.8728914) q[3];
sx q[3];
rz(-1.8541981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8687826) q[0];
sx q[0];
rz(-1.9406809) q[0];
sx q[0];
rz(1.1707206) q[0];
rz(-2.6314578) q[1];
sx q[1];
rz(-1.7850103) q[1];
sx q[1];
rz(-1.2957113) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0062795) q[0];
sx q[0];
rz(-1.8989925) q[0];
sx q[0];
rz(2.6562064) q[0];
rz(-2.2175118) q[2];
sx q[2];
rz(-2.7087822) q[2];
sx q[2];
rz(2.2518287) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.012055339) q[1];
sx q[1];
rz(-2.4730198) q[1];
sx q[1];
rz(1.936391) q[1];
rz(3.1136884) q[3];
sx q[3];
rz(-2.0114261) q[3];
sx q[3];
rz(0.84907109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.43508139) q[2];
sx q[2];
rz(-1.4979829) q[2];
sx q[2];
rz(-0.56813017) q[2];
rz(-0.98012296) q[3];
sx q[3];
rz(-2.1025434) q[3];
sx q[3];
rz(0.19395104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66529626) q[0];
sx q[0];
rz(-0.81403533) q[0];
sx q[0];
rz(-2.4639159) q[0];
rz(0.19605818) q[1];
sx q[1];
rz(-2.129107) q[1];
sx q[1];
rz(-0.83818865) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3121376) q[0];
sx q[0];
rz(-1.2078309) q[0];
sx q[0];
rz(-1.3296933) q[0];
rz(-1.499275) q[2];
sx q[2];
rz(-2.3841249) q[2];
sx q[2];
rz(2.1069991) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5327685) q[1];
sx q[1];
rz(-0.84034398) q[1];
sx q[1];
rz(-0.26809147) q[1];
x q[2];
rz(-2.8832199) q[3];
sx q[3];
rz(-1.4900472) q[3];
sx q[3];
rz(-2.0852058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3395386) q[2];
sx q[2];
rz(-2.4513117) q[2];
sx q[2];
rz(-1.8748803) q[2];
rz(0.96261111) q[3];
sx q[3];
rz(-1.5714785) q[3];
sx q[3];
rz(-0.094749711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4203913) q[0];
sx q[0];
rz(-1.9045916) q[0];
sx q[0];
rz(0.23751968) q[0];
rz(2.1233842) q[1];
sx q[1];
rz(-0.84914452) q[1];
sx q[1];
rz(-2.9097897) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7983539) q[0];
sx q[0];
rz(-2.624247) q[0];
sx q[0];
rz(-1.377064) q[0];
x q[1];
rz(-0.78328697) q[2];
sx q[2];
rz(-2.8056393) q[2];
sx q[2];
rz(-0.17620262) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17813645) q[1];
sx q[1];
rz(-1.3061211) q[1];
sx q[1];
rz(-1.4241649) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3123355) q[3];
sx q[3];
rz(-2.4371394) q[3];
sx q[3];
rz(-0.36650141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1101749) q[2];
sx q[2];
rz(-1.2538223) q[2];
sx q[2];
rz(-1.0882264) q[2];
rz(2.7534289) q[3];
sx q[3];
rz(-2.4813014) q[3];
sx q[3];
rz(-2.3378519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5768455) q[0];
sx q[0];
rz(-1.7871465) q[0];
sx q[0];
rz(-0.47252895) q[0];
rz(-2.172773) q[1];
sx q[1];
rz(-0.70822721) q[1];
sx q[1];
rz(0.72475564) q[1];
rz(0.62933915) q[2];
sx q[2];
rz(-0.39471252) q[2];
sx q[2];
rz(-0.82804745) q[2];
rz(2.9536392) q[3];
sx q[3];
rz(-1.5516075) q[3];
sx q[3];
rz(2.8967378) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
