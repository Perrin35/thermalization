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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0955015) q[0];
sx q[0];
rz(-1.7016194) q[0];
sx q[0];
rz(-3.1216937) q[0];
x q[1];
rz(1.0257172) q[2];
sx q[2];
rz(-1.8946049) q[2];
sx q[2];
rz(2.2562502) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.284277) q[1];
sx q[1];
rz(-2.2419062) q[1];
sx q[1];
rz(1.7637232) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5159357) q[3];
sx q[3];
rz(-1.2513972) q[3];
sx q[3];
rz(-2.0351792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.75498092) q[2];
sx q[2];
rz(-1.745801) q[2];
sx q[2];
rz(2.6888729) q[2];
rz(0.1581986) q[3];
sx q[3];
rz(-2.4499564) q[3];
sx q[3];
rz(0.89481568) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2125856) q[0];
sx q[0];
rz(-2.043262) q[0];
sx q[0];
rz(1.1520977) q[0];
rz(1.903803) q[1];
sx q[1];
rz(-1.5367616) q[1];
sx q[1];
rz(0.47098413) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4244714) q[0];
sx q[0];
rz(-2.0524128) q[0];
sx q[0];
rz(-0.061901285) q[0];
x q[1];
rz(-2.9071964) q[2];
sx q[2];
rz(-0.98653754) q[2];
sx q[2];
rz(2.5258979) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4091332) q[1];
sx q[1];
rz(-1.0832548) q[1];
sx q[1];
rz(3.0797466) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4460375) q[3];
sx q[3];
rz(-1.5091981) q[3];
sx q[3];
rz(1.0775281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.058078893) q[2];
sx q[2];
rz(-0.60516548) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.6059748) q[1];
sx q[1];
rz(-2.7022865) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1742451) q[0];
sx q[0];
rz(-0.085701533) q[0];
sx q[0];
rz(-1.9232737) q[0];
x q[1];
rz(-0.36053948) q[2];
sx q[2];
rz(-1.5224824) q[2];
sx q[2];
rz(-2.609032) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5369536) q[1];
sx q[1];
rz(-0.96452689) q[1];
sx q[1];
rz(0.21142516) q[1];
rz(-2.1268232) q[3];
sx q[3];
rz(-1.945567) q[3];
sx q[3];
rz(2.8127363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.97418857) q[2];
sx q[2];
rz(-1.5524652) q[2];
sx q[2];
rz(-2.9411194) q[2];
rz(-2.386507) q[3];
sx q[3];
rz(-1.0198159) q[3];
sx q[3];
rz(0.42373207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-0.077483594) q[0];
sx q[0];
rz(-1.5492726) q[0];
sx q[0];
rz(-1.2444929) q[0];
rz(-0.81047932) q[1];
sx q[1];
rz(-1.3296209) q[1];
sx q[1];
rz(-2.2669852) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8133102) q[0];
sx q[0];
rz(-0.62424849) q[0];
sx q[0];
rz(-1.3024131) q[0];
x q[1];
rz(0.34747296) q[2];
sx q[2];
rz(-0.407019) q[2];
sx q[2];
rz(-2.4328872) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3687467) q[1];
sx q[1];
rz(-0.44279848) q[1];
sx q[1];
rz(0.55261353) q[1];
rz(-pi) q[2];
rz(1.4960257) q[3];
sx q[3];
rz(-1.6755591) q[3];
sx q[3];
rz(-1.9024224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13742927) q[2];
sx q[2];
rz(-0.88976294) q[2];
sx q[2];
rz(-1.4271663) q[2];
rz(-3.0754722) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(1.4495513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3604597) q[0];
sx q[0];
rz(-0.17372486) q[0];
sx q[0];
rz(2.5710035) q[0];
rz(-0.55496201) q[1];
sx q[1];
rz(-2.404232) q[1];
sx q[1];
rz(0.74329174) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0665019) q[0];
sx q[0];
rz(-0.69426232) q[0];
sx q[0];
rz(-1.2230722) q[0];
rz(1.4666124) q[2];
sx q[2];
rz(-1.0717234) q[2];
sx q[2];
rz(-2.4521329) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9651282) q[1];
sx q[1];
rz(-1.3511409) q[1];
sx q[1];
rz(-1.3419282) q[1];
x q[2];
rz(2.8652142) q[3];
sx q[3];
rz(-2.31782) q[3];
sx q[3];
rz(0.36229047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1935929) q[2];
sx q[2];
rz(-1.3698545) q[2];
sx q[2];
rz(2.8732079) q[2];
rz(2.0949481) q[3];
sx q[3];
rz(-0.36473754) q[3];
sx q[3];
rz(0.31782761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.627581) q[0];
sx q[0];
rz(-1.6126957) q[0];
sx q[0];
rz(-2.6348689) q[0];
rz(-0.2535893) q[1];
sx q[1];
rz(-1.2713623) q[1];
sx q[1];
rz(-2.0862897) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2282288) q[0];
sx q[0];
rz(-0.95634395) q[0];
sx q[0];
rz(0.22426228) q[0];
rz(-2.6199117) q[2];
sx q[2];
rz(-2.1573967) q[2];
sx q[2];
rz(-2.4751543) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.036990449) q[1];
sx q[1];
rz(-1.1739507) q[1];
sx q[1];
rz(-0.87581046) q[1];
rz(-pi) q[2];
rz(0.25552337) q[3];
sx q[3];
rz(-2.3821085) q[3];
sx q[3];
rz(-1.4042735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6598597) q[2];
sx q[2];
rz(-1.0446171) q[2];
sx q[2];
rz(0.8824904) q[2];
rz(-2.4957538) q[3];
sx q[3];
rz(-1.9927988) q[3];
sx q[3];
rz(-1.283949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6362474) q[0];
sx q[0];
rz(-1.212965) q[0];
sx q[0];
rz(2.3690467) q[0];
rz(-1.729471) q[1];
sx q[1];
rz(-1.9344784) q[1];
sx q[1];
rz(-0.57377446) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44874292) q[0];
sx q[0];
rz(-2.9208555) q[0];
sx q[0];
rz(-2.0953395) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7710118) q[2];
sx q[2];
rz(-1.5491312) q[2];
sx q[2];
rz(1.7553165) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3473914) q[1];
sx q[1];
rz(-0.83354356) q[1];
sx q[1];
rz(-2.3458523) q[1];
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
rz(-3.1022169) q[2];
sx q[2];
rz(-2.6874459) q[2];
sx q[2];
rz(2.3708564) q[2];
rz(-2.7052774) q[3];
sx q[3];
rz(-1.8728914) q[3];
sx q[3];
rz(1.8541981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.27281) q[0];
sx q[0];
rz(-1.9406809) q[0];
sx q[0];
rz(-1.1707206) q[0];
rz(0.51013485) q[1];
sx q[1];
rz(-1.7850103) q[1];
sx q[1];
rz(1.8458813) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8745236) q[0];
sx q[0];
rz(-1.1133615) q[0];
sx q[0];
rz(1.9382856) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2715271) q[2];
sx q[2];
rz(-1.9121133) q[2];
sx q[2];
rz(-1.5835539) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8746652) q[1];
sx q[1];
rz(-1.7942567) q[1];
sx q[1];
rz(0.93519559) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1136884) q[3];
sx q[3];
rz(-2.0114261) q[3];
sx q[3];
rz(-2.2925216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7065113) q[2];
sx q[2];
rz(-1.4979829) q[2];
sx q[2];
rz(-2.5734625) q[2];
rz(2.1614697) q[3];
sx q[3];
rz(-2.1025434) q[3];
sx q[3];
rz(-2.9476416) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66529626) q[0];
sx q[0];
rz(-2.3275573) q[0];
sx q[0];
rz(-0.67767674) q[0];
rz(-2.9455345) q[1];
sx q[1];
rz(-2.129107) q[1];
sx q[1];
rz(2.303404) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4351589) q[0];
sx q[0];
rz(-0.43276946) q[0];
sx q[0];
rz(0.5612527) q[0];
rz(-0.81460641) q[2];
sx q[2];
rz(-1.5216773) q[2];
sx q[2];
rz(-0.58821046) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2182525) q[1];
sx q[1];
rz(-2.372101) q[1];
sx q[1];
rz(1.8583276) q[1];
rz(-2.8348654) q[3];
sx q[3];
rz(-0.2704276) q[3];
sx q[3];
rz(-0.21817792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.80205408) q[2];
sx q[2];
rz(-2.4513117) q[2];
sx q[2];
rz(-1.8748803) q[2];
rz(0.96261111) q[3];
sx q[3];
rz(-1.5701141) q[3];
sx q[3];
rz(0.094749711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4203913) q[0];
sx q[0];
rz(-1.2370011) q[0];
sx q[0];
rz(2.904073) q[0];
rz(-1.0182084) q[1];
sx q[1];
rz(-0.84914452) q[1];
sx q[1];
rz(-2.9097897) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12122614) q[0];
sx q[0];
rz(-2.0775284) q[0];
sx q[0];
rz(0.10911848) q[0];
rz(-pi) q[1];
rz(2.3583057) q[2];
sx q[2];
rz(-2.8056393) q[2];
sx q[2];
rz(2.96539) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4494891) q[1];
sx q[1];
rz(-2.8398501) q[1];
sx q[1];
rz(-0.49441378) q[1];
rz(-2.9276804) q[3];
sx q[3];
rz(-0.89424664) q[3];
sx q[3];
rz(3.1090581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1101749) q[2];
sx q[2];
rz(-1.8877703) q[2];
sx q[2];
rz(2.0533662) q[2];
rz(2.7534289) q[3];
sx q[3];
rz(-0.66029125) q[3];
sx q[3];
rz(-0.80374074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5768455) q[0];
sx q[0];
rz(-1.7871465) q[0];
sx q[0];
rz(-0.47252895) q[0];
rz(2.172773) q[1];
sx q[1];
rz(-2.4333654) q[1];
sx q[1];
rz(-2.416837) q[1];
rz(-1.3303403) q[2];
sx q[2];
rz(-1.8869055) q[2];
sx q[2];
rz(-1.4958924) q[2];
rz(1.5512636) q[3];
sx q[3];
rz(-1.7587147) q[3];
sx q[3];
rz(1.3295909) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
