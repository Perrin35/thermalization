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
rz(-1.0840253) q[1];
sx q[1];
rz(1.6853583) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0460912) q[0];
sx q[0];
rz(-1.4399733) q[0];
sx q[0];
rz(-0.019898947) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37402447) q[2];
sx q[2];
rz(-2.0846539) q[2];
sx q[2];
rz(-0.49486578) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.97988843) q[1];
sx q[1];
rz(-0.6941422) q[1];
sx q[1];
rz(0.23692268) q[1];
x q[2];
rz(-0.16430328) q[3];
sx q[3];
rz(-0.32391732) q[3];
sx q[3];
rz(1.8620373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3866117) q[2];
sx q[2];
rz(-1.3957916) q[2];
sx q[2];
rz(0.45271978) q[2];
rz(-0.1581986) q[3];
sx q[3];
rz(-0.69163624) q[3];
sx q[3];
rz(0.89481568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.92900705) q[0];
sx q[0];
rz(-1.0983306) q[0];
sx q[0];
rz(-1.1520977) q[0];
rz(-1.903803) q[1];
sx q[1];
rz(-1.5367616) q[1];
sx q[1];
rz(-0.47098413) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.02397) q[0];
sx q[0];
rz(-1.6256486) q[0];
sx q[0];
rz(-1.0883925) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1678796) q[2];
sx q[2];
rz(-1.7657585) q[2];
sx q[2];
rz(-0.82414579) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.86733782) q[1];
sx q[1];
rz(-1.6254289) q[1];
sx q[1];
rz(1.0824624) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0456411) q[3];
sx q[3];
rz(-0.6978242) q[3];
sx q[3];
rz(-0.56688353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.058078893) q[2];
sx q[2];
rz(-2.5364272) q[2];
sx q[2];
rz(0.2557959) q[2];
rz(-1.4852218) q[3];
sx q[3];
rz(-1.1896313) q[3];
sx q[3];
rz(-2.9799057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.8783022) q[0];
sx q[0];
rz(-2.0354164) q[0];
sx q[0];
rz(-2.8702452) q[0];
rz(2.4052606) q[1];
sx q[1];
rz(-1.5356179) q[1];
sx q[1];
rz(-0.43930611) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1742451) q[0];
sx q[0];
rz(-0.085701533) q[0];
sx q[0];
rz(-1.9232737) q[0];
x q[1];
rz(0.36053948) q[2];
sx q[2];
rz(-1.6191102) q[2];
sx q[2];
rz(-2.609032) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9860552) q[1];
sx q[1];
rz(-1.7441161) q[1];
sx q[1];
rz(0.9539414) q[1];
x q[2];
rz(-0.4337173) q[3];
sx q[3];
rz(-1.0573514) q[3];
sx q[3];
rz(1.0182667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1674041) q[2];
sx q[2];
rz(-1.5524652) q[2];
sx q[2];
rz(-0.20047323) q[2];
rz(0.75508562) q[3];
sx q[3];
rz(-1.0198159) q[3];
sx q[3];
rz(-2.7178606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0641091) q[0];
sx q[0];
rz(-1.5492726) q[0];
sx q[0];
rz(1.2444929) q[0];
rz(0.81047932) q[1];
sx q[1];
rz(-1.8119718) q[1];
sx q[1];
rz(0.8746075) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8133102) q[0];
sx q[0];
rz(-0.62424849) q[0];
sx q[0];
rz(1.3024131) q[0];
rz(-pi) q[1];
rz(1.4250408) q[2];
sx q[2];
rz(-1.1894023) q[2];
sx q[2];
rz(-2.8085453) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9676535) q[1];
sx q[1];
rz(-1.1974918) q[1];
sx q[1];
rz(1.8147545) q[1];
rz(0.10505418) q[3];
sx q[3];
rz(-1.6451562) q[3];
sx q[3];
rz(-0.33945938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0041634) q[2];
sx q[2];
rz(-2.2518297) q[2];
sx q[2];
rz(1.4271663) q[2];
rz(0.066120474) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(-1.6920413) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3604597) q[0];
sx q[0];
rz(-2.9678678) q[0];
sx q[0];
rz(2.5710035) q[0];
rz(0.55496201) q[1];
sx q[1];
rz(-2.404232) q[1];
sx q[1];
rz(2.3983009) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3742204) q[0];
sx q[0];
rz(-1.351007) q[0];
sx q[0];
rz(0.90669294) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18851738) q[2];
sx q[2];
rz(-2.6326615) q[2];
sx q[2];
rz(0.9045507) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1764644) q[1];
sx q[1];
rz(-1.7904518) q[1];
sx q[1];
rz(1.3419282) q[1];
rz(-pi) q[2];
rz(2.3371313) q[3];
sx q[3];
rz(-1.7723697) q[3];
sx q[3];
rz(1.7427012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1935929) q[2];
sx q[2];
rz(-1.7717382) q[2];
sx q[2];
rz(-0.26838475) q[2];
rz(2.0949481) q[3];
sx q[3];
rz(-2.7768551) q[3];
sx q[3];
rz(2.823765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.627581) q[0];
sx q[0];
rz(-1.6126957) q[0];
sx q[0];
rz(-0.50672379) q[0];
rz(-0.2535893) q[1];
sx q[1];
rz(-1.8702303) q[1];
sx q[1];
rz(-1.0553029) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5366093) q[0];
sx q[0];
rz(-0.6491001) q[0];
sx q[0];
rz(1.2654632) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2141586) q[2];
sx q[2];
rz(-0.76403996) q[2];
sx q[2];
rz(1.6709136) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9678952) q[1];
sx q[1];
rz(-2.3579862) q[1];
sx q[1];
rz(-0.99131363) q[1];
rz(-pi) q[2];
rz(1.33527) q[3];
sx q[3];
rz(-0.84170656) q[3];
sx q[3];
rz(-2.083076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.48173299) q[2];
sx q[2];
rz(-1.0446171) q[2];
sx q[2];
rz(0.8824904) q[2];
rz(-0.64583889) q[3];
sx q[3];
rz(-1.1487938) q[3];
sx q[3];
rz(-1.283949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5053453) q[0];
sx q[0];
rz(-1.9286276) q[0];
sx q[0];
rz(-0.77254599) q[0];
rz(1.4121217) q[1];
sx q[1];
rz(-1.9344784) q[1];
sx q[1];
rz(-0.57377446) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.505578) q[0];
sx q[0];
rz(-1.4609219) q[0];
sx q[0];
rz(-1.3789603) q[0];
rz(1.7710118) q[2];
sx q[2];
rz(-1.5924615) q[2];
sx q[2];
rz(-1.3862762) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7942012) q[1];
sx q[1];
rz(-2.3080491) q[1];
sx q[1];
rz(0.79574037) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1366828) q[3];
sx q[3];
rz(-2.8661869) q[3];
sx q[3];
rz(0.79808455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1022169) q[2];
sx q[2];
rz(-2.6874459) q[2];
sx q[2];
rz(2.3708564) q[2];
rz(-0.43631521) q[3];
sx q[3];
rz(-1.8728914) q[3];
sx q[3];
rz(-1.8541981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.27281) q[0];
sx q[0];
rz(-1.9406809) q[0];
sx q[0];
rz(1.1707206) q[0];
rz(-0.51013485) q[1];
sx q[1];
rz(-1.3565823) q[1];
sx q[1];
rz(1.8458813) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1575748) q[0];
sx q[0];
rz(-2.5630953) q[0];
sx q[0];
rz(2.5111141) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8700656) q[2];
sx q[2];
rz(-1.2294793) q[2];
sx q[2];
rz(1.5580387) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8746652) q[1];
sx q[1];
rz(-1.7942567) q[1];
sx q[1];
rz(-2.2063971) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5116974) q[3];
sx q[3];
rz(-2.700138) q[3];
sx q[3];
rz(0.7837226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.43508139) q[2];
sx q[2];
rz(-1.6436098) q[2];
sx q[2];
rz(2.5734625) q[2];
rz(-0.98012296) q[3];
sx q[3];
rz(-1.0390493) q[3];
sx q[3];
rz(2.9476416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66529626) q[0];
sx q[0];
rz(-2.3275573) q[0];
sx q[0];
rz(0.67767674) q[0];
rz(0.19605818) q[1];
sx q[1];
rz(-1.0124857) q[1];
sx q[1];
rz(-2.303404) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4351589) q[0];
sx q[0];
rz(-2.7088232) q[0];
sx q[0];
rz(2.58034) q[0];
rz(3.0741192) q[2];
sx q[2];
rz(-2.3258492) q[2];
sx q[2];
rz(2.2052854) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.14324489) q[1];
sx q[1];
rz(-1.76941) q[1];
sx q[1];
rz(2.3193588) q[1];
rz(0.25837274) q[3];
sx q[3];
rz(-1.4900472) q[3];
sx q[3];
rz(-2.0852058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3395386) q[2];
sx q[2];
rz(-0.69028091) q[2];
sx q[2];
rz(-1.8748803) q[2];
rz(-2.1789815) q[3];
sx q[3];
rz(-1.5701141) q[3];
sx q[3];
rz(0.094749711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.4203913) q[0];
sx q[0];
rz(-1.9045916) q[0];
sx q[0];
rz(-0.23751968) q[0];
rz(-2.1233842) q[1];
sx q[1];
rz(-2.2924481) q[1];
sx q[1];
rz(0.231803) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34323877) q[0];
sx q[0];
rz(-0.51734561) q[0];
sx q[0];
rz(-1.377064) q[0];
rz(-pi) q[1];
rz(-2.3583057) q[2];
sx q[2];
rz(-2.8056393) q[2];
sx q[2];
rz(-2.96539) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4494891) q[1];
sx q[1];
rz(-0.30174258) q[1];
sx q[1];
rz(-0.49441378) q[1];
rz(-pi) q[2];
x q[2];
rz(0.2139123) q[3];
sx q[3];
rz(-2.247346) q[3];
sx q[3];
rz(0.032534508) q[3];
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
rz(2.0533662) q[2];
rz(-2.7534289) q[3];
sx q[3];
rz(-0.66029125) q[3];
sx q[3];
rz(0.80374074) q[3];
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
rz(-1.5768455) q[0];
sx q[0];
rz(-1.3544461) q[0];
sx q[0];
rz(2.6690637) q[0];
rz(-0.96881962) q[1];
sx q[1];
rz(-2.4333654) q[1];
sx q[1];
rz(-2.416837) q[1];
rz(0.3248365) q[2];
sx q[2];
rz(-1.3424716) q[2];
sx q[2];
rz(-2.9906103) q[2];
rz(-0.10235056) q[3];
sx q[3];
rz(-0.1889189) q[3];
sx q[3];
rz(1.2253996) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];