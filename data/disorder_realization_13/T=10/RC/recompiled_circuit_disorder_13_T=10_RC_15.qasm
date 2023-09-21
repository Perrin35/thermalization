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
rz(2.2464377) q[0];
rz(-0.31495467) q[1];
sx q[1];
rz(-2.0575674) q[1];
sx q[1];
rz(-1.6853583) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0955015) q[0];
sx q[0];
rz(-1.4399733) q[0];
sx q[0];
rz(-0.019898947) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1452791) q[2];
sx q[2];
rz(-2.516054) q[2];
sx q[2];
rz(1.1686981) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8573157) q[1];
sx q[1];
rz(-2.2419062) q[1];
sx q[1];
rz(-1.3778694) q[1];
rz(-pi) q[2];
rz(0.16430328) q[3];
sx q[3];
rz(-2.8176753) q[3];
sx q[3];
rz(1.8620373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3866117) q[2];
sx q[2];
rz(-1.3957916) q[2];
sx q[2];
rz(0.45271978) q[2];
rz(-0.1581986) q[3];
sx q[3];
rz(-2.4499564) q[3];
sx q[3];
rz(2.246777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2125856) q[0];
sx q[0];
rz(-2.043262) q[0];
sx q[0];
rz(1.989495) q[0];
rz(-1.903803) q[1];
sx q[1];
rz(-1.6048311) q[1];
sx q[1];
rz(-2.6706085) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1176227) q[0];
sx q[0];
rz(-1.6256486) q[0];
sx q[0];
rz(1.0883925) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9085625) q[2];
sx q[2];
rz(-2.517189) q[2];
sx q[2];
rz(1.0242467) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.7324595) q[1];
sx q[1];
rz(-2.0583378) q[1];
sx q[1];
rz(-0.061846102) q[1];
rz(-3.0456411) q[3];
sx q[3];
rz(-2.4437685) q[3];
sx q[3];
rz(-0.56688353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0835138) q[2];
sx q[2];
rz(-2.5364272) q[2];
sx q[2];
rz(-0.2557959) q[2];
rz(1.6563709) q[3];
sx q[3];
rz(-1.1896313) q[3];
sx q[3];
rz(0.16168693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26329041) q[0];
sx q[0];
rz(-2.0354164) q[0];
sx q[0];
rz(-0.27134744) q[0];
rz(-2.4052606) q[1];
sx q[1];
rz(-1.5356179) q[1];
sx q[1];
rz(0.43930611) q[1];
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
rz(-1.6224242) q[2];
sx q[2];
rz(-1.2106967) q[2];
sx q[2];
rz(-2.1215631) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.15553741) q[1];
sx q[1];
rz(-1.3974766) q[1];
sx q[1];
rz(-2.1876513) q[1];
rz(2.1268232) q[3];
sx q[3];
rz(-1.945567) q[3];
sx q[3];
rz(-2.8127363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1674041) q[2];
sx q[2];
rz(-1.5891275) q[2];
sx q[2];
rz(2.9411194) q[2];
rz(2.386507) q[3];
sx q[3];
rz(-1.0198159) q[3];
sx q[3];
rz(2.7178606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0641091) q[0];
sx q[0];
rz(-1.5492726) q[0];
sx q[0];
rz(-1.2444929) q[0];
rz(2.3311133) q[1];
sx q[1];
rz(-1.3296209) q[1];
sx q[1];
rz(0.8746075) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32828242) q[0];
sx q[0];
rz(-2.5173442) q[0];
sx q[0];
rz(-1.8391795) q[0];
x q[1];
rz(0.38509102) q[2];
sx q[2];
rz(-1.706012) q[2];
sx q[2];
rz(-1.9584292) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.7728459) q[1];
sx q[1];
rz(-2.6987942) q[1];
sx q[1];
rz(-2.5889791) q[1];
rz(-pi) q[2];
rz(-0.61769684) q[3];
sx q[3];
rz(-3.0129637) q[3];
sx q[3];
rz(2.5240412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.13742927) q[2];
sx q[2];
rz(-2.2518297) q[2];
sx q[2];
rz(-1.7144263) q[2];
rz(3.0754722) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(1.6920413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3604597) q[0];
sx q[0];
rz(-0.17372486) q[0];
sx q[0];
rz(2.5710035) q[0];
rz(-2.5866306) q[1];
sx q[1];
rz(-2.404232) q[1];
sx q[1];
rz(2.3983009) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0750908) q[0];
sx q[0];
rz(-0.69426232) q[0];
sx q[0];
rz(-1.9185205) q[0];
rz(-2.6402316) q[2];
sx q[2];
rz(-1.4793581) q[2];
sx q[2];
rz(-0.83133343) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.78391) q[1];
sx q[1];
rz(-0.3158814) q[1];
sx q[1];
rz(0.79343474) q[1];
rz(-pi) q[2];
rz(-2.3371313) q[3];
sx q[3];
rz(-1.3692229) q[3];
sx q[3];
rz(-1.3988914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1935929) q[2];
sx q[2];
rz(-1.3698545) q[2];
sx q[2];
rz(0.26838475) q[2];
rz(2.0949481) q[3];
sx q[3];
rz(-0.36473754) q[3];
sx q[3];
rz(0.31782761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5140117) q[0];
sx q[0];
rz(-1.6126957) q[0];
sx q[0];
rz(-2.6348689) q[0];
rz(2.8880033) q[1];
sx q[1];
rz(-1.8702303) q[1];
sx q[1];
rz(2.0862897) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5366093) q[0];
sx q[0];
rz(-2.4924926) q[0];
sx q[0];
rz(-1.8761294) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2141586) q[2];
sx q[2];
rz(-0.76403996) q[2];
sx q[2];
rz(1.6709136) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1046022) q[1];
sx q[1];
rz(-1.9676419) q[1];
sx q[1];
rz(2.2657822) q[1];
x q[2];
rz(0.7430195) q[3];
sx q[3];
rz(-1.395874) q[3];
sx q[3];
rz(2.7878441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.48173299) q[2];
sx q[2];
rz(-2.0969756) q[2];
sx q[2];
rz(-2.2591023) q[2];
rz(-0.64583889) q[3];
sx q[3];
rz(-1.9927988) q[3];
sx q[3];
rz(1.283949) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5053453) q[0];
sx q[0];
rz(-1.9286276) q[0];
sx q[0];
rz(-0.77254599) q[0];
rz(-1.729471) q[1];
sx q[1];
rz(-1.9344784) q[1];
sx q[1];
rz(2.5678182) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6928497) q[0];
sx q[0];
rz(-0.22073711) q[0];
sx q[0];
rz(-1.0462532) q[0];
rz(1.6793208) q[2];
sx q[2];
rz(-0.2013686) q[2];
sx q[2];
rz(3.0634207) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7942012) q[1];
sx q[1];
rz(-2.3080491) q[1];
sx q[1];
rz(-0.79574037) q[1];
rz(-2.9912234) q[3];
sx q[3];
rz(-1.339185) q[3];
sx q[3];
rz(-1.7600972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.039375719) q[2];
sx q[2];
rz(-0.45414671) q[2];
sx q[2];
rz(-2.3708564) q[2];
rz(2.7052774) q[3];
sx q[3];
rz(-1.8728914) q[3];
sx q[3];
rz(1.2873945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.27281) q[0];
sx q[0];
rz(-1.9406809) q[0];
sx q[0];
rz(1.1707206) q[0];
rz(2.6314578) q[1];
sx q[1];
rz(-1.3565823) q[1];
sx q[1];
rz(1.8458813) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2670691) q[0];
sx q[0];
rz(-2.0282312) q[0];
sx q[0];
rz(-1.2033071) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92408085) q[2];
sx q[2];
rz(-0.43281049) q[2];
sx q[2];
rz(0.88976394) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.46591972) q[1];
sx q[1];
rz(-2.188176) q[1];
sx q[1];
rz(-0.27523756) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0115764) q[3];
sx q[3];
rz(-1.545558) q[3];
sx q[3];
rz(2.4079635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7065113) q[2];
sx q[2];
rz(-1.4979829) q[2];
sx q[2];
rz(-0.56813017) q[2];
rz(2.1614697) q[3];
sx q[3];
rz(-2.1025434) q[3];
sx q[3];
rz(-2.9476416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66529626) q[0];
sx q[0];
rz(-2.3275573) q[0];
sx q[0];
rz(-0.67767674) q[0];
rz(0.19605818) q[1];
sx q[1];
rz(-2.129107) q[1];
sx q[1];
rz(-0.83818865) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3121376) q[0];
sx q[0];
rz(-1.2078309) q[0];
sx q[0];
rz(-1.8118993) q[0];
x q[1];
rz(-0.81460641) q[2];
sx q[2];
rz(-1.6199154) q[2];
sx q[2];
rz(-2.5533822) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9983478) q[1];
sx q[1];
rz(-1.3721826) q[1];
sx q[1];
rz(0.82223383) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8832199) q[3];
sx q[3];
rz(-1.4900472) q[3];
sx q[3];
rz(2.0852058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3395386) q[2];
sx q[2];
rz(-0.69028091) q[2];
sx q[2];
rz(-1.8748803) q[2];
rz(2.1789815) q[3];
sx q[3];
rz(-1.5701141) q[3];
sx q[3];
rz(3.0468429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.4203913) q[0];
sx q[0];
rz(-1.2370011) q[0];
sx q[0];
rz(-0.23751968) q[0];
rz(-1.0182084) q[1];
sx q[1];
rz(-0.84914452) q[1];
sx q[1];
rz(0.231803) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12122614) q[0];
sx q[0];
rz(-1.0640642) q[0];
sx q[0];
rz(-0.10911848) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24256369) q[2];
sx q[2];
rz(-1.3360268) q[2];
sx q[2];
rz(0.64005062) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3540436) q[1];
sx q[1];
rz(-1.7122867) q[1];
sx q[1];
rz(-2.8741807) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9276804) q[3];
sx q[3];
rz(-0.89424664) q[3];
sx q[3];
rz(0.032534508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1101749) q[2];
sx q[2];
rz(-1.8877703) q[2];
sx q[2];
rz(2.0533662) q[2];
rz(-2.7534289) q[3];
sx q[3];
rz(-2.4813014) q[3];
sx q[3];
rz(-0.80374074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5768455) q[0];
sx q[0];
rz(-1.3544461) q[0];
sx q[0];
rz(2.6690637) q[0];
rz(-2.172773) q[1];
sx q[1];
rz(-0.70822721) q[1];
sx q[1];
rz(0.72475564) q[1];
rz(-1.8112524) q[2];
sx q[2];
rz(-1.2546872) q[2];
sx q[2];
rz(1.6457002) q[2];
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
