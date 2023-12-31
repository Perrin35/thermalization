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
rz(-1.4562343) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0955015) q[0];
sx q[0];
rz(-1.7016194) q[0];
sx q[0];
rz(-3.1216937) q[0];
x q[1];
rz(0.99631359) q[2];
sx q[2];
rz(-0.62553863) q[2];
sx q[2];
rz(1.9728945) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.97988843) q[1];
sx q[1];
rz(-0.6941422) q[1];
sx q[1];
rz(0.23692268) q[1];
rz(1.625657) q[3];
sx q[3];
rz(-1.2513972) q[3];
sx q[3];
rz(2.0351792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.75498092) q[2];
sx q[2];
rz(-1.745801) q[2];
sx q[2];
rz(0.45271978) q[2];
rz(2.9833941) q[3];
sx q[3];
rz(-2.4499564) q[3];
sx q[3];
rz(-0.89481568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.92900705) q[0];
sx q[0];
rz(-1.0983306) q[0];
sx q[0];
rz(-1.1520977) q[0];
rz(1.903803) q[1];
sx q[1];
rz(-1.5367616) q[1];
sx q[1];
rz(-2.6706085) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4244714) q[0];
sx q[0];
rz(-1.0891799) q[0];
sx q[0];
rz(-0.061901285) q[0];
rz(-pi) q[1];
rz(-2.1678796) q[2];
sx q[2];
rz(-1.3758341) q[2];
sx q[2];
rz(-2.3174469) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4091332) q[1];
sx q[1];
rz(-2.0583378) q[1];
sx q[1];
rz(3.0797466) q[1];
rz(-pi) q[2];
rz(1.4906293) q[3];
sx q[3];
rz(-2.2647694) q[3];
sx q[3];
rz(2.699664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0835138) q[2];
sx q[2];
rz(-2.5364272) q[2];
sx q[2];
rz(0.2557959) q[2];
rz(-1.4852218) q[3];
sx q[3];
rz(-1.1896313) q[3];
sx q[3];
rz(0.16168693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.8783022) q[0];
sx q[0];
rz(-1.1061763) q[0];
sx q[0];
rz(0.27134744) q[0];
rz(-2.4052606) q[1];
sx q[1];
rz(-1.6059748) q[1];
sx q[1];
rz(2.7022865) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8205748) q[0];
sx q[0];
rz(-1.6512172) q[0];
sx q[0];
rz(-0.029650173) q[0];
rz(-pi) q[1];
rz(1.6224242) q[2];
sx q[2];
rz(-1.9308959) q[2];
sx q[2];
rz(1.0200295) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.5369536) q[1];
sx q[1];
rz(-2.1770658) q[1];
sx q[1];
rz(-2.9301675) q[1];
x q[2];
rz(-2.2112591) q[3];
sx q[3];
rz(-0.6593245) q[3];
sx q[3];
rz(-1.7742771) q[3];
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
rz(2.386507) q[3];
sx q[3];
rz(-1.0198159) q[3];
sx q[3];
rz(-0.42373207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-3.0641091) q[0];
sx q[0];
rz(-1.5923201) q[0];
sx q[0];
rz(1.8970998) q[0];
rz(-0.81047932) q[1];
sx q[1];
rz(-1.3296209) q[1];
sx q[1];
rz(0.8746075) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0229605) q[0];
sx q[0];
rz(-1.4151787) q[0];
sx q[0];
rz(-0.96373425) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7941197) q[2];
sx q[2];
rz(-0.407019) q[2];
sx q[2];
rz(-0.70870542) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.306327) q[1];
sx q[1];
rz(-1.3439461) q[1];
sx q[1];
rz(-0.38362417) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6455669) q[3];
sx q[3];
rz(-1.4660335) q[3];
sx q[3];
rz(-1.9024224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.13742927) q[2];
sx q[2];
rz(-0.88976294) q[2];
sx q[2];
rz(-1.4271663) q[2];
rz(-3.0754722) q[3];
sx q[3];
rz(-2.7757006) q[3];
sx q[3];
rz(1.6920413) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3604597) q[0];
sx q[0];
rz(-0.17372486) q[0];
sx q[0];
rz(-0.57058913) q[0];
rz(0.55496201) q[1];
sx q[1];
rz(-2.404232) q[1];
sx q[1];
rz(2.3983009) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0665019) q[0];
sx q[0];
rz(-0.69426232) q[0];
sx q[0];
rz(-1.9185205) q[0];
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
sx q[0];
rz(-pi/2) q[0];
rz(2.78391) q[1];
sx q[1];
rz(-2.8257113) q[1];
sx q[1];
rz(2.3481579) q[1];
x q[2];
rz(0.27637847) q[3];
sx q[3];
rz(-0.82377269) q[3];
sx q[3];
rz(-2.7793022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1935929) q[2];
sx q[2];
rz(-1.3698545) q[2];
sx q[2];
rz(-0.26838475) q[2];
rz(2.0949481) q[3];
sx q[3];
rz(-0.36473754) q[3];
sx q[3];
rz(0.31782761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5140117) q[0];
sx q[0];
rz(-1.6126957) q[0];
sx q[0];
rz(0.50672379) q[0];
rz(0.2535893) q[1];
sx q[1];
rz(-1.8702303) q[1];
sx q[1];
rz(1.0553029) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6049833) q[0];
sx q[0];
rz(-0.6491001) q[0];
sx q[0];
rz(1.2654632) q[0];
rz(-pi) q[1];
rz(-0.91674532) q[2];
sx q[2];
rz(-1.9987717) q[2];
sx q[2];
rz(-2.5452754) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1046022) q[1];
sx q[1];
rz(-1.1739507) q[1];
sx q[1];
rz(-2.2657822) q[1];
rz(2.3985732) q[3];
sx q[3];
rz(-1.395874) q[3];
sx q[3];
rz(0.35374853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6598597) q[2];
sx q[2];
rz(-1.0446171) q[2];
sx q[2];
rz(2.2591023) q[2];
rz(0.64583889) q[3];
sx q[3];
rz(-1.1487938) q[3];
sx q[3];
rz(-1.8576436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5053453) q[0];
sx q[0];
rz(-1.9286276) q[0];
sx q[0];
rz(2.3690467) q[0];
rz(-1.729471) q[1];
sx q[1];
rz(-1.2071143) q[1];
sx q[1];
rz(0.57377446) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44874292) q[0];
sx q[0];
rz(-2.9208555) q[0];
sx q[0];
rz(-1.0462532) q[0];
rz(-pi) q[1];
rz(1.6793208) q[2];
sx q[2];
rz(-2.9402241) q[2];
sx q[2];
rz(-3.0634207) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7942012) q[1];
sx q[1];
rz(-0.83354356) q[1];
sx q[1];
rz(-2.3458523) q[1];
rz(-pi) q[2];
rz(1.8049559) q[3];
sx q[3];
rz(-1.7171211) q[3];
sx q[3];
rz(-0.22406604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.039375719) q[2];
sx q[2];
rz(-2.6874459) q[2];
sx q[2];
rz(-2.3708564) q[2];
rz(2.7052774) q[3];
sx q[3];
rz(-1.8728914) q[3];
sx q[3];
rz(-1.8541981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8687826) q[0];
sx q[0];
rz(-1.2009118) q[0];
sx q[0];
rz(1.9708721) q[0];
rz(0.51013485) q[1];
sx q[1];
rz(-1.3565823) q[1];
sx q[1];
rz(1.2957113) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1575748) q[0];
sx q[0];
rz(-0.57849738) q[0];
sx q[0];
rz(0.63047854) q[0];
rz(-pi) q[1];
rz(-1.2175351) q[2];
sx q[2];
rz(-1.3152939) q[2];
sx q[2];
rz(0.080163408) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.012055339) q[1];
sx q[1];
rz(-0.66857282) q[1];
sx q[1];
rz(1.936391) q[1];
rz(-0.027904228) q[3];
sx q[3];
rz(-2.0114261) q[3];
sx q[3];
rz(0.84907109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7065113) q[2];
sx q[2];
rz(-1.4979829) q[2];
sx q[2];
rz(2.5734625) q[2];
rz(2.1614697) q[3];
sx q[3];
rz(-2.1025434) q[3];
sx q[3];
rz(-2.9476416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4762964) q[0];
sx q[0];
rz(-2.3275573) q[0];
sx q[0];
rz(-0.67767674) q[0];
rz(-0.19605818) q[1];
sx q[1];
rz(-2.129107) q[1];
sx q[1];
rz(-2.303404) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8294551) q[0];
sx q[0];
rz(-1.2078309) q[0];
sx q[0];
rz(-1.3296933) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3269862) q[2];
sx q[2];
rz(-1.6199154) q[2];
sx q[2];
rz(0.58821046) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9233401) q[1];
sx q[1];
rz(-2.372101) q[1];
sx q[1];
rz(-1.8583276) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8348654) q[3];
sx q[3];
rz(-0.2704276) q[3];
sx q[3];
rz(-2.9234147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3395386) q[2];
sx q[2];
rz(-0.69028091) q[2];
sx q[2];
rz(1.2667123) q[2];
rz(2.1789815) q[3];
sx q[3];
rz(-1.5714785) q[3];
sx q[3];
rz(-3.0468429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72120136) q[0];
sx q[0];
rz(-1.9045916) q[0];
sx q[0];
rz(2.904073) q[0];
rz(-1.0182084) q[1];
sx q[1];
rz(-2.2924481) q[1];
sx q[1];
rz(-0.231803) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3964513) q[0];
sx q[0];
rz(-1.475435) q[0];
sx q[0];
rz(1.061529) q[0];
rz(-1.3292153) q[2];
sx q[2];
rz(-1.3350147) q[2];
sx q[2];
rz(0.98824046) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.17813645) q[1];
sx q[1];
rz(-1.3061211) q[1];
sx q[1];
rz(1.4241649) q[1];
rz(-0.2139123) q[3];
sx q[3];
rz(-0.89424664) q[3];
sx q[3];
rz(-3.1090581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0314177) q[2];
sx q[2];
rz(-1.8877703) q[2];
sx q[2];
rz(2.0533662) q[2];
rz(0.38816372) q[3];
sx q[3];
rz(-0.66029125) q[3];
sx q[3];
rz(0.80374074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5768455) q[0];
sx q[0];
rz(-1.7871465) q[0];
sx q[0];
rz(-0.47252895) q[0];
rz(0.96881962) q[1];
sx q[1];
rz(-0.70822721) q[1];
sx q[1];
rz(0.72475564) q[1];
rz(-0.3248365) q[2];
sx q[2];
rz(-1.799121) q[2];
sx q[2];
rz(0.15098235) q[2];
rz(-0.18795342) q[3];
sx q[3];
rz(-1.5516075) q[3];
sx q[3];
rz(2.8967378) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
