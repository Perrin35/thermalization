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
rz(2.826638) q[1];
sx q[1];
rz(-1.0840253) q[1];
sx q[1];
rz(-1.4562343) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47269868) q[0];
sx q[0];
rz(-1.5510674) q[0];
sx q[0];
rz(-1.701645) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99631359) q[2];
sx q[2];
rz(-2.516054) q[2];
sx q[2];
rz(-1.9728945) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.734182) q[1];
sx q[1];
rz(-1.721518) q[1];
sx q[1];
rz(-2.461344) q[1];
rz(0.31984826) q[3];
sx q[3];
rz(-1.6228798) q[3];
sx q[3];
rz(-2.6944514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3866117) q[2];
sx q[2];
rz(-1.745801) q[2];
sx q[2];
rz(-0.45271978) q[2];
rz(2.9833941) q[3];
sx q[3];
rz(-2.4499564) q[3];
sx q[3];
rz(-0.89481568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2125856) q[0];
sx q[0];
rz(-1.0983306) q[0];
sx q[0];
rz(1.989495) q[0];
rz(-1.903803) q[1];
sx q[1];
rz(-1.6048311) q[1];
sx q[1];
rz(0.47098413) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7171213) q[0];
sx q[0];
rz(-1.0891799) q[0];
sx q[0];
rz(-0.061901285) q[0];
x q[1];
rz(-2.9071964) q[2];
sx q[2];
rz(-0.98653754) q[2];
sx q[2];
rz(2.5258979) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5405611) q[1];
sx q[1];
rz(-2.6504576) q[1];
sx q[1];
rz(1.4547552) q[1];
rz(-pi) q[2];
rz(-2.4460375) q[3];
sx q[3];
rz(-1.5091981) q[3];
sx q[3];
rz(2.0640645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.058078893) q[2];
sx q[2];
rz(-0.60516548) q[2];
sx q[2];
rz(-2.8857968) q[2];
rz(1.4852218) q[3];
sx q[3];
rz(-1.9519613) q[3];
sx q[3];
rz(0.16168693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26329041) q[0];
sx q[0];
rz(-2.0354164) q[0];
sx q[0];
rz(-2.8702452) q[0];
rz(0.73633206) q[1];
sx q[1];
rz(-1.5356179) q[1];
sx q[1];
rz(0.43930611) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8894316) q[0];
sx q[0];
rz(-1.541242) q[0];
sx q[0];
rz(-1.4903402) q[0];
rz(-pi) q[1];
rz(3.0053824) q[2];
sx q[2];
rz(-2.7779707) q[2];
sx q[2];
rz(-1.9759535) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5369536) q[1];
sx q[1];
rz(-2.1770658) q[1];
sx q[1];
rz(2.9301675) q[1];
x q[2];
rz(2.2112591) q[3];
sx q[3];
rz(-0.6593245) q[3];
sx q[3];
rz(1.7742771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.97418857) q[2];
sx q[2];
rz(-1.5891275) q[2];
sx q[2];
rz(2.9411194) q[2];
rz(-2.386507) q[3];
sx q[3];
rz(-2.1217767) q[3];
sx q[3];
rz(2.7178606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(3.0641091) q[0];
sx q[0];
rz(-1.5923201) q[0];
sx q[0];
rz(-1.8970998) q[0];
rz(2.3311133) q[1];
sx q[1];
rz(-1.3296209) q[1];
sx q[1];
rz(0.8746075) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8133102) q[0];
sx q[0];
rz(-0.62424849) q[0];
sx q[0];
rz(1.3024131) q[0];
rz(-pi) q[1];
rz(2.7941197) q[2];
sx q[2];
rz(-0.407019) q[2];
sx q[2];
rz(-0.70870542) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.306327) q[1];
sx q[1];
rz(-1.3439461) q[1];
sx q[1];
rz(0.38362417) q[1];
x q[2];
rz(-0.61769684) q[3];
sx q[3];
rz(-3.0129637) q[3];
sx q[3];
rz(2.5240412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.13742927) q[2];
sx q[2];
rz(-2.2518297) q[2];
sx q[2];
rz(1.7144263) q[2];
rz(3.0754722) q[3];
sx q[3];
rz(-2.7757006) q[3];
sx q[3];
rz(1.4495513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3604597) q[0];
sx q[0];
rz(-0.17372486) q[0];
sx q[0];
rz(0.57058913) q[0];
rz(2.5866306) q[1];
sx q[1];
rz(-2.404232) q[1];
sx q[1];
rz(0.74329174) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0750908) q[0];
sx q[0];
rz(-0.69426232) q[0];
sx q[0];
rz(1.2230722) q[0];
rz(-pi) q[1];
rz(2.6402316) q[2];
sx q[2];
rz(-1.4793581) q[2];
sx q[2];
rz(0.83133343) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1764644) q[1];
sx q[1];
rz(-1.7904518) q[1];
sx q[1];
rz(1.7996644) q[1];
x q[2];
rz(0.80446135) q[3];
sx q[3];
rz(-1.7723697) q[3];
sx q[3];
rz(1.3988914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1935929) q[2];
sx q[2];
rz(-1.3698545) q[2];
sx q[2];
rz(-0.26838475) q[2];
rz(-2.0949481) q[3];
sx q[3];
rz(-0.36473754) q[3];
sx q[3];
rz(2.823765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.627581) q[0];
sx q[0];
rz(-1.528897) q[0];
sx q[0];
rz(-0.50672379) q[0];
rz(-2.8880033) q[1];
sx q[1];
rz(-1.2713623) q[1];
sx q[1];
rz(2.0862897) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5366093) q[0];
sx q[0];
rz(-2.4924926) q[0];
sx q[0];
rz(1.8761294) q[0];
rz(0.91674532) q[2];
sx q[2];
rz(-1.9987717) q[2];
sx q[2];
rz(-0.59631729) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1736974) q[1];
sx q[1];
rz(-0.78360644) q[1];
sx q[1];
rz(0.99131363) q[1];
rz(0.25552337) q[3];
sx q[3];
rz(-0.75948411) q[3];
sx q[3];
rz(-1.7373191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6598597) q[2];
sx q[2];
rz(-1.0446171) q[2];
sx q[2];
rz(-0.8824904) q[2];
rz(-0.64583889) q[3];
sx q[3];
rz(-1.9927988) q[3];
sx q[3];
rz(-1.8576436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5053453) q[0];
sx q[0];
rz(-1.9286276) q[0];
sx q[0];
rz(0.77254599) q[0];
rz(1.4121217) q[1];
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
rz(-3.0550802) q[0];
sx q[0];
rz(-1.7614613) q[0];
sx q[0];
rz(3.029682) q[0];
rz(-pi) q[1];
rz(-1.3705809) q[2];
sx q[2];
rz(-1.5491312) q[2];
sx q[2];
rz(1.3862762) q[2];
rz(-pi) q[3];
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
rz(0.79574037) q[1];
x q[2];
rz(-1.8049559) q[3];
sx q[3];
rz(-1.4244716) q[3];
sx q[3];
rz(2.9175266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.039375719) q[2];
sx q[2];
rz(-0.45414671) q[2];
sx q[2];
rz(-0.77073628) q[2];
rz(2.7052774) q[3];
sx q[3];
rz(-1.2687012) q[3];
sx q[3];
rz(-1.2873945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8687826) q[0];
sx q[0];
rz(-1.9406809) q[0];
sx q[0];
rz(-1.1707206) q[0];
rz(2.6314578) q[1];
sx q[1];
rz(-1.7850103) q[1];
sx q[1];
rz(1.2957113) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13531317) q[0];
sx q[0];
rz(-1.8989925) q[0];
sx q[0];
rz(-0.4853863) q[0];
rz(-pi) q[1];
rz(-2.2175118) q[2];
sx q[2];
rz(-0.43281049) q[2];
sx q[2];
rz(-2.2518287) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2669275) q[1];
sx q[1];
rz(-1.7942567) q[1];
sx q[1];
rz(-0.93519559) q[1];
rz(1.6298953) q[3];
sx q[3];
rz(-0.44145465) q[3];
sx q[3];
rz(2.3578701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.43508139) q[2];
sx q[2];
rz(-1.6436098) q[2];
sx q[2];
rz(-2.5734625) q[2];
rz(0.98012296) q[3];
sx q[3];
rz(-2.1025434) q[3];
sx q[3];
rz(-0.19395104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.4762964) q[0];
sx q[0];
rz(-0.81403533) q[0];
sx q[0];
rz(-2.4639159) q[0];
rz(-2.9455345) q[1];
sx q[1];
rz(-1.0124857) q[1];
sx q[1];
rz(-2.303404) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4351589) q[0];
sx q[0];
rz(-2.7088232) q[0];
sx q[0];
rz(2.58034) q[0];
rz(0.067473472) q[2];
sx q[2];
rz(-0.81574342) q[2];
sx q[2];
rz(-0.93630723) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5327685) q[1];
sx q[1];
rz(-0.84034398) q[1];
sx q[1];
rz(0.26809147) q[1];
rz(-2.8348654) q[3];
sx q[3];
rz(-2.8711651) q[3];
sx q[3];
rz(0.21817792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.80205408) q[2];
sx q[2];
rz(-0.69028091) q[2];
sx q[2];
rz(-1.2667123) q[2];
rz(-0.96261111) q[3];
sx q[3];
rz(-1.5701141) q[3];
sx q[3];
rz(-0.094749711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
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
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0203665) q[0];
sx q[0];
rz(-1.0640642) q[0];
sx q[0];
rz(0.10911848) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8123774) q[2];
sx q[2];
rz(-1.806578) q[2];
sx q[2];
rz(-2.1533522) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.69210359) q[1];
sx q[1];
rz(-2.8398501) q[1];
sx q[1];
rz(0.49441378) q[1];
x q[2];
rz(0.88296367) q[3];
sx q[3];
rz(-1.4045047) q[3];
sx q[3];
rz(-1.7385141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1101749) q[2];
sx q[2];
rz(-1.8877703) q[2];
sx q[2];
rz(-1.0882264) q[2];
rz(2.7534289) q[3];
sx q[3];
rz(-2.4813014) q[3];
sx q[3];
rz(0.80374074) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
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
rz(1.3303403) q[2];
sx q[2];
rz(-1.2546872) q[2];
sx q[2];
rz(1.6457002) q[2];
rz(-1.5512636) q[3];
sx q[3];
rz(-1.3828779) q[3];
sx q[3];
rz(-1.8120017) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
