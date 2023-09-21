OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.17833248) q[0];
sx q[0];
rz(-1.4890716) q[0];
sx q[0];
rz(-0.89515495) q[0];
rz(2.826638) q[1];
sx q[1];
rz(2.0575674) q[1];
sx q[1];
rz(10.881012) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0460912) q[0];
sx q[0];
rz(-1.4399733) q[0];
sx q[0];
rz(-0.019898947) q[0];
rz(0.37402447) q[2];
sx q[2];
rz(-2.0846539) q[2];
sx q[2];
rz(-2.6467269) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8573157) q[1];
sx q[1];
rz(-2.2419062) q[1];
sx q[1];
rz(1.3778694) q[1];
rz(-0.16430328) q[3];
sx q[3];
rz(-0.32391732) q[3];
sx q[3];
rz(-1.2795554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.75498092) q[2];
sx q[2];
rz(-1.745801) q[2];
sx q[2];
rz(2.6888729) q[2];
rz(-0.1581986) q[3];
sx q[3];
rz(-2.4499564) q[3];
sx q[3];
rz(2.246777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
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
rz(1.989495) q[0];
rz(1.903803) q[1];
sx q[1];
rz(-1.6048311) q[1];
sx q[1];
rz(-0.47098413) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7171213) q[0];
sx q[0];
rz(-2.0524128) q[0];
sx q[0];
rz(-0.061901285) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1678796) q[2];
sx q[2];
rz(-1.7657585) q[2];
sx q[2];
rz(2.3174469) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5405611) q[1];
sx q[1];
rz(-2.6504576) q[1];
sx q[1];
rz(1.4547552) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.095951565) q[3];
sx q[3];
rz(-0.6978242) q[3];
sx q[3];
rz(-0.56688353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0835138) q[2];
sx q[2];
rz(-2.5364272) q[2];
sx q[2];
rz(-2.8857968) q[2];
rz(1.6563709) q[3];
sx q[3];
rz(-1.1896313) q[3];
sx q[3];
rz(-2.9799057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8783022) q[0];
sx q[0];
rz(-2.0354164) q[0];
sx q[0];
rz(-0.27134744) q[0];
rz(-0.73633206) q[1];
sx q[1];
rz(-1.5356179) q[1];
sx q[1];
rz(-0.43930611) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8894316) q[0];
sx q[0];
rz(-1.541242) q[0];
sx q[0];
rz(-1.4903402) q[0];
rz(-pi) q[1];
rz(1.5191684) q[2];
sx q[2];
rz(-1.2106967) q[2];
sx q[2];
rz(-2.1215631) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1767133) q[1];
sx q[1];
rz(-0.63767725) q[1];
sx q[1];
rz(-1.2769075) q[1];
rz(2.1268232) q[3];
sx q[3];
rz(-1.1960256) q[3];
sx q[3];
rz(2.8127363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.97418857) q[2];
sx q[2];
rz(-1.5891275) q[2];
sx q[2];
rz(-2.9411194) q[2];
rz(-2.386507) q[3];
sx q[3];
rz(-2.1217767) q[3];
sx q[3];
rz(2.7178606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.077483594) q[0];
sx q[0];
rz(-1.5923201) q[0];
sx q[0];
rz(1.2444929) q[0];
rz(-2.3311133) q[1];
sx q[1];
rz(-1.3296209) q[1];
sx q[1];
rz(-0.8746075) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8133102) q[0];
sx q[0];
rz(-0.62424849) q[0];
sx q[0];
rz(-1.3024131) q[0];
rz(-pi) q[1];
rz(-1.7165519) q[2];
sx q[2];
rz(-1.1894023) q[2];
sx q[2];
rz(0.33304735) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9676535) q[1];
sx q[1];
rz(-1.1974918) q[1];
sx q[1];
rz(-1.8147545) q[1];
rz(-1.6455669) q[3];
sx q[3];
rz(-1.4660335) q[3];
sx q[3];
rz(-1.2391702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0041634) q[2];
sx q[2];
rz(-0.88976294) q[2];
sx q[2];
rz(-1.7144263) q[2];
rz(3.0754722) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(-1.4495513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3604597) q[0];
sx q[0];
rz(-2.9678678) q[0];
sx q[0];
rz(-0.57058913) q[0];
rz(-2.5866306) q[1];
sx q[1];
rz(-2.404232) q[1];
sx q[1];
rz(-0.74329174) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63440454) q[0];
sx q[0];
rz(-2.2162063) q[0];
sx q[0];
rz(2.8651644) q[0];
x q[1];
rz(0.18851738) q[2];
sx q[2];
rz(-0.50893116) q[2];
sx q[2];
rz(0.9045507) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.78391) q[1];
sx q[1];
rz(-0.3158814) q[1];
sx q[1];
rz(-0.79343474) q[1];
rz(-pi) q[2];
x q[2];
rz(1.284243) q[3];
sx q[3];
rz(-0.78714579) q[3];
sx q[3];
rz(3.1084276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9479998) q[2];
sx q[2];
rz(-1.7717382) q[2];
sx q[2];
rz(-2.8732079) q[2];
rz(-1.0466446) q[3];
sx q[3];
rz(-0.36473754) q[3];
sx q[3];
rz(-2.823765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5140117) q[0];
sx q[0];
rz(-1.6126957) q[0];
sx q[0];
rz(0.50672379) q[0];
rz(-0.2535893) q[1];
sx q[1];
rz(-1.8702303) q[1];
sx q[1];
rz(2.0862897) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21181606) q[0];
sx q[0];
rz(-1.7535216) q[0];
sx q[0];
rz(0.94434785) q[0];
rz(2.2141586) q[2];
sx q[2];
rz(-0.76403996) q[2];
sx q[2];
rz(-1.6709136) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.222059) q[1];
sx q[1];
rz(-2.2026081) q[1];
sx q[1];
rz(-0.49948378) q[1];
rz(-pi) q[2];
rz(-0.7430195) q[3];
sx q[3];
rz(-1.395874) q[3];
sx q[3];
rz(-2.7878441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6598597) q[2];
sx q[2];
rz(-2.0969756) q[2];
sx q[2];
rz(-0.8824904) q[2];
rz(0.64583889) q[3];
sx q[3];
rz(-1.9927988) q[3];
sx q[3];
rz(1.8576436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6362474) q[0];
sx q[0];
rz(-1.9286276) q[0];
sx q[0];
rz(-2.3690467) q[0];
rz(1.729471) q[1];
sx q[1];
rz(-1.9344784) q[1];
sx q[1];
rz(0.57377446) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44874292) q[0];
sx q[0];
rz(-2.9208555) q[0];
sx q[0];
rz(-2.0953395) q[0];
x q[1];
rz(1.7710118) q[2];
sx q[2];
rz(-1.5924615) q[2];
sx q[2];
rz(-1.3862762) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.37807235) q[1];
sx q[1];
rz(-1.0135279) q[1];
sx q[1];
rz(2.4850363) q[1];
rz(-2.1366828) q[3];
sx q[3];
rz(-2.8661869) q[3];
sx q[3];
rz(2.3435081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.039375719) q[2];
sx q[2];
rz(-0.45414671) q[2];
sx q[2];
rz(0.77073628) q[2];
rz(-2.7052774) q[3];
sx q[3];
rz(-1.8728914) q[3];
sx q[3];
rz(1.8541981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.27281) q[0];
sx q[0];
rz(-1.2009118) q[0];
sx q[0];
rz(-1.9708721) q[0];
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
rz(1.8745236) q[0];
sx q[0];
rz(-2.0282312) q[0];
sx q[0];
rz(1.2033071) q[0];
rz(-1.2175351) q[2];
sx q[2];
rz(-1.8262987) q[2];
sx q[2];
rz(-0.080163408) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.012055339) q[1];
sx q[1];
rz(-2.4730198) q[1];
sx q[1];
rz(1.936391) q[1];
x q[2];
rz(-1.6298953) q[3];
sx q[3];
rz(-2.700138) q[3];
sx q[3];
rz(2.3578701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7065113) q[2];
sx q[2];
rz(-1.6436098) q[2];
sx q[2];
rz(-0.56813017) q[2];
rz(-2.1614697) q[3];
sx q[3];
rz(-2.1025434) q[3];
sx q[3];
rz(2.9476416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66529626) q[0];
sx q[0];
rz(-2.3275573) q[0];
sx q[0];
rz(0.67767674) q[0];
rz(2.9455345) q[1];
sx q[1];
rz(-1.0124857) q[1];
sx q[1];
rz(2.303404) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3121376) q[0];
sx q[0];
rz(-1.2078309) q[0];
sx q[0];
rz(1.8118993) q[0];
x q[1];
rz(-1.499275) q[2];
sx q[2];
rz(-0.75746775) q[2];
sx q[2];
rz(-2.1069991) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2182525) q[1];
sx q[1];
rz(-0.76949161) q[1];
sx q[1];
rz(-1.283265) q[1];
rz(-pi) q[2];
rz(0.30672726) q[3];
sx q[3];
rz(-2.8711651) q[3];
sx q[3];
rz(-2.9234147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.80205408) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.4203913) q[0];
sx q[0];
rz(-1.9045916) q[0];
sx q[0];
rz(2.904073) q[0];
rz(-2.1233842) q[1];
sx q[1];
rz(-2.2924481) q[1];
sx q[1];
rz(-2.9097897) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7451413) q[0];
sx q[0];
rz(-1.475435) q[0];
sx q[0];
rz(2.0800637) q[0];
x q[1];
rz(-2.899029) q[2];
sx q[2];
rz(-1.8055658) q[2];
sx q[2];
rz(2.501542) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4494891) q[1];
sx q[1];
rz(-0.30174258) q[1];
sx q[1];
rz(-2.6471789) q[1];
x q[2];
rz(-0.88296367) q[3];
sx q[3];
rz(-1.737088) q[3];
sx q[3];
rz(-1.7385141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0314177) q[2];
sx q[2];
rz(-1.8877703) q[2];
sx q[2];
rz(-1.0882264) q[2];
rz(-0.38816372) q[3];
sx q[3];
rz(-0.66029125) q[3];
sx q[3];
rz(2.3378519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-0.62933915) q[2];
sx q[2];
rz(-2.7468801) q[2];
sx q[2];
rz(2.3135452) q[2];
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