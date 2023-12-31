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
rz(4.6306643) q[0];
sx q[0];
rz(10.319933) q[0];
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
rz(-1.0955015) q[0];
sx q[0];
rz(-1.7016194) q[0];
sx q[0];
rz(-0.019898947) q[0];
rz(0.37402447) q[2];
sx q[2];
rz(-2.0846539) q[2];
sx q[2];
rz(0.49486578) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.284277) q[1];
sx q[1];
rz(-2.2419062) q[1];
sx q[1];
rz(1.3778694) q[1];
rz(-pi) q[2];
rz(-0.16430328) q[3];
sx q[3];
rz(-0.32391732) q[3];
sx q[3];
rz(-1.2795554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92900705) q[0];
sx q[0];
rz(-1.0983306) q[0];
sx q[0];
rz(-1.1520977) q[0];
rz(1.2377897) q[1];
sx q[1];
rz(-1.5367616) q[1];
sx q[1];
rz(2.6706085) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4244714) q[0];
sx q[0];
rz(-2.0524128) q[0];
sx q[0];
rz(-0.061901285) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9071964) q[2];
sx q[2];
rz(-2.1550551) q[2];
sx q[2];
rz(-2.5258979) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.86733782) q[1];
sx q[1];
rz(-1.5161637) q[1];
sx q[1];
rz(-1.0824624) q[1];
x q[2];
rz(-1.6509634) q[3];
sx q[3];
rz(-0.87682322) q[3];
sx q[3];
rz(0.44192867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0835138) q[2];
sx q[2];
rz(-2.5364272) q[2];
sx q[2];
rz(2.8857968) q[2];
rz(-1.6563709) q[3];
sx q[3];
rz(-1.1896313) q[3];
sx q[3];
rz(2.9799057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(2.7022865) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96734756) q[0];
sx q[0];
rz(-3.0558911) q[0];
sx q[0];
rz(1.9232737) q[0];
rz(-pi) q[1];
rz(1.6224242) q[2];
sx q[2];
rz(-1.9308959) q[2];
sx q[2];
rz(-2.1215631) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9648793) q[1];
sx q[1];
rz(-0.63767725) q[1];
sx q[1];
rz(-1.2769075) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.4337173) q[3];
sx q[3];
rz(-2.0842413) q[3];
sx q[3];
rz(-1.0182667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.97418857) q[2];
sx q[2];
rz(-1.5524652) q[2];
sx q[2];
rz(-0.20047323) q[2];
rz(2.386507) q[3];
sx q[3];
rz(-2.1217767) q[3];
sx q[3];
rz(0.42373207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.077483594) q[0];
sx q[0];
rz(-1.5923201) q[0];
sx q[0];
rz(1.2444929) q[0];
rz(-0.81047932) q[1];
sx q[1];
rz(-1.3296209) q[1];
sx q[1];
rz(-2.2669852) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1186322) q[0];
sx q[0];
rz(-1.4151787) q[0];
sx q[0];
rz(-2.1778584) q[0];
x q[1];
rz(1.4250408) q[2];
sx q[2];
rz(-1.9521904) q[2];
sx q[2];
rz(2.8085453) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.306327) q[1];
sx q[1];
rz(-1.3439461) q[1];
sx q[1];
rz(-0.38362417) q[1];
x q[2];
rz(-3.0365385) q[3];
sx q[3];
rz(-1.6451562) q[3];
sx q[3];
rz(2.8021333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0041634) q[2];
sx q[2];
rz(-2.2518297) q[2];
sx q[2];
rz(1.4271663) q[2];
rz(0.066120474) q[3];
sx q[3];
rz(-2.7757006) q[3];
sx q[3];
rz(-1.4495513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3604597) q[0];
sx q[0];
rz(-0.17372486) q[0];
sx q[0];
rz(-2.5710035) q[0];
rz(0.55496201) q[1];
sx q[1];
rz(-2.404232) q[1];
sx q[1];
rz(-0.74329174) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3742204) q[0];
sx q[0];
rz(-1.7905856) q[0];
sx q[0];
rz(-2.2348997) q[0];
x q[1];
rz(-1.4666124) q[2];
sx q[2];
rz(-2.0698692) q[2];
sx q[2];
rz(-2.4521329) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.35768269) q[1];
sx q[1];
rz(-0.3158814) q[1];
sx q[1];
rz(-0.79343474) q[1];
x q[2];
rz(-2.3371313) q[3];
sx q[3];
rz(-1.7723697) q[3];
sx q[3];
rz(-1.7427012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9479998) q[2];
sx q[2];
rz(-1.7717382) q[2];
sx q[2];
rz(2.8732079) q[2];
rz(-2.0949481) q[3];
sx q[3];
rz(-0.36473754) q[3];
sx q[3];
rz(2.823765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5140117) q[0];
sx q[0];
rz(-1.6126957) q[0];
sx q[0];
rz(2.6348689) q[0];
rz(0.2535893) q[1];
sx q[1];
rz(-1.2713623) q[1];
sx q[1];
rz(-1.0553029) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2282288) q[0];
sx q[0];
rz(-2.1852487) q[0];
sx q[0];
rz(2.9173304) q[0];
x q[1];
rz(-2.6199117) q[2];
sx q[2];
rz(-2.1573967) q[2];
sx q[2];
rz(0.66643836) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.222059) q[1];
sx q[1];
rz(-0.93898458) q[1];
sx q[1];
rz(-0.49948378) q[1];
rz(-pi) q[2];
rz(1.8063227) q[3];
sx q[3];
rz(-0.84170656) q[3];
sx q[3];
rz(2.083076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6598597) q[2];
sx q[2];
rz(-1.0446171) q[2];
sx q[2];
rz(-2.2591023) q[2];
rz(0.64583889) q[3];
sx q[3];
rz(-1.9927988) q[3];
sx q[3];
rz(-1.283949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6362474) q[0];
sx q[0];
rz(-1.9286276) q[0];
sx q[0];
rz(2.3690467) q[0];
rz(1.729471) q[1];
sx q[1];
rz(-1.9344784) q[1];
sx q[1];
rz(-2.5678182) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.505578) q[0];
sx q[0];
rz(-1.4609219) q[0];
sx q[0];
rz(1.3789603) q[0];
rz(-pi) q[1];
rz(-1.7710118) q[2];
sx q[2];
rz(-1.5491312) q[2];
sx q[2];
rz(1.7553165) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.80572762) q[1];
sx q[1];
rz(-2.1153567) q[1];
sx q[1];
rz(-0.90421275) q[1];
rz(-pi) q[2];
rz(-1.3366367) q[3];
sx q[3];
rz(-1.7171211) q[3];
sx q[3];
rz(-0.22406604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.039375719) q[2];
sx q[2];
rz(-2.6874459) q[2];
sx q[2];
rz(-2.3708564) q[2];
rz(-0.43631521) q[3];
sx q[3];
rz(-1.2687012) q[3];
sx q[3];
rz(-1.2873945) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8687826) q[0];
sx q[0];
rz(-1.2009118) q[0];
sx q[0];
rz(-1.9708721) q[0];
rz(-0.51013485) q[1];
sx q[1];
rz(-1.7850103) q[1];
sx q[1];
rz(-1.8458813) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2670691) q[0];
sx q[0];
rz(-1.1133615) q[0];
sx q[0];
rz(-1.2033071) q[0];
x q[1];
rz(0.92408085) q[2];
sx q[2];
rz(-2.7087822) q[2];
sx q[2];
rz(2.2518287) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2669275) q[1];
sx q[1];
rz(-1.3473359) q[1];
sx q[1];
rz(0.93519559) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5116974) q[3];
sx q[3];
rz(-2.700138) q[3];
sx q[3];
rz(-2.3578701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7065113) q[2];
sx q[2];
rz(-1.6436098) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66529626) q[0];
sx q[0];
rz(-0.81403533) q[0];
sx q[0];
rz(2.4639159) q[0];
rz(2.9455345) q[1];
sx q[1];
rz(-2.129107) q[1];
sx q[1];
rz(0.83818865) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7064338) q[0];
sx q[0];
rz(-0.43276946) q[0];
sx q[0];
rz(-0.5612527) q[0];
rz(0.81460641) q[2];
sx q[2];
rz(-1.5216773) q[2];
sx q[2];
rz(0.58821046) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2182525) q[1];
sx q[1];
rz(-2.372101) q[1];
sx q[1];
rz(-1.283265) q[1];
rz(2.8832199) q[3];
sx q[3];
rz(-1.4900472) q[3];
sx q[3];
rz(-1.0563869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.80205408) q[2];
sx q[2];
rz(-0.69028091) q[2];
sx q[2];
rz(-1.2667123) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72120136) q[0];
sx q[0];
rz(-1.2370011) q[0];
sx q[0];
rz(-2.904073) q[0];
rz(2.1233842) q[1];
sx q[1];
rz(-0.84914452) q[1];
sx q[1];
rz(-2.9097897) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12122614) q[0];
sx q[0];
rz(-2.0775284) q[0];
sx q[0];
rz(3.0324742) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.899029) q[2];
sx q[2];
rz(-1.8055658) q[2];
sx q[2];
rz(-0.64005062) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9634562) q[1];
sx q[1];
rz(-1.8354715) q[1];
sx q[1];
rz(1.7174277) q[1];
rz(2.9276804) q[3];
sx q[3];
rz(-0.89424664) q[3];
sx q[3];
rz(-3.1090581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0314177) q[2];
sx q[2];
rz(-1.2538223) q[2];
sx q[2];
rz(-2.0533662) q[2];
rz(0.38816372) q[3];
sx q[3];
rz(-2.4813014) q[3];
sx q[3];
rz(-0.80374074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5647472) q[0];
sx q[0];
rz(-1.7871465) q[0];
sx q[0];
rz(-0.47252895) q[0];
rz(-0.96881962) q[1];
sx q[1];
rz(-2.4333654) q[1];
sx q[1];
rz(-2.416837) q[1];
rz(2.8167562) q[2];
sx q[2];
rz(-1.799121) q[2];
sx q[2];
rz(0.15098235) q[2];
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
