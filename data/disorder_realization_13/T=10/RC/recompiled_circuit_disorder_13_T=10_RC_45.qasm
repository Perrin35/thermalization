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
rz(-0.31495467) q[1];
sx q[1];
rz(-2.0575674) q[1];
sx q[1];
rz(-1.6853583) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2468949) q[0];
sx q[0];
rz(-0.1323192) q[0];
sx q[0];
rz(-1.7208862) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37402447) q[2];
sx q[2];
rz(-1.0569388) q[2];
sx q[2];
rz(-2.6467269) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.97988843) q[1];
sx q[1];
rz(-2.4474505) q[1];
sx q[1];
rz(-0.23692268) q[1];
rz(2.9772894) q[3];
sx q[3];
rz(-2.8176753) q[3];
sx q[3];
rz(-1.8620373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3866117) q[2];
sx q[2];
rz(-1.745801) q[2];
sx q[2];
rz(0.45271978) q[2];
rz(-2.9833941) q[3];
sx q[3];
rz(-0.69163624) q[3];
sx q[3];
rz(2.246777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2125856) q[0];
sx q[0];
rz(-1.0983306) q[0];
sx q[0];
rz(1.1520977) q[0];
rz(1.2377897) q[1];
sx q[1];
rz(-1.5367616) q[1];
sx q[1];
rz(-0.47098413) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5841056) q[0];
sx q[0];
rz(-2.6563245) q[0];
sx q[0];
rz(1.4529865) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9071964) q[2];
sx q[2];
rz(-0.98653754) q[2];
sx q[2];
rz(0.61569475) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4091332) q[1];
sx q[1];
rz(-2.0583378) q[1];
sx q[1];
rz(0.061846102) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0456411) q[3];
sx q[3];
rz(-0.6978242) q[3];
sx q[3];
rz(-2.5747091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0835138) q[2];
sx q[2];
rz(-0.60516548) q[2];
sx q[2];
rz(0.2557959) q[2];
rz(-1.6563709) q[3];
sx q[3];
rz(-1.9519613) q[3];
sx q[3];
rz(-2.9799057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8783022) q[0];
sx q[0];
rz(-1.1061763) q[0];
sx q[0];
rz(-0.27134744) q[0];
rz(-2.4052606) q[1];
sx q[1];
rz(-1.5356179) q[1];
sx q[1];
rz(0.43930611) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1742451) q[0];
sx q[0];
rz(-3.0558911) q[0];
sx q[0];
rz(1.9232737) q[0];
rz(-pi) q[1];
rz(-2.7810532) q[2];
sx q[2];
rz(-1.5224824) q[2];
sx q[2];
rz(-0.53256065) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9648793) q[1];
sx q[1];
rz(-0.63767725) q[1];
sx q[1];
rz(-1.8646851) q[1];
x q[2];
rz(-0.4337173) q[3];
sx q[3];
rz(-1.0573514) q[3];
sx q[3];
rz(1.0182667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.97418857) q[2];
sx q[2];
rz(-1.5524652) q[2];
sx q[2];
rz(-2.9411194) q[2];
rz(0.75508562) q[3];
sx q[3];
rz(-1.0198159) q[3];
sx q[3];
rz(-2.7178606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
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
rz(-1.8970998) q[0];
rz(2.3311133) q[1];
sx q[1];
rz(-1.3296209) q[1];
sx q[1];
rz(0.8746075) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8133102) q[0];
sx q[0];
rz(-0.62424849) q[0];
sx q[0];
rz(1.3024131) q[0];
x q[1];
rz(1.4250408) q[2];
sx q[2];
rz(-1.9521904) q[2];
sx q[2];
rz(-0.33304735) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.7728459) q[1];
sx q[1];
rz(-0.44279848) q[1];
sx q[1];
rz(-0.55261353) q[1];
rz(0.61769684) q[3];
sx q[3];
rz(-3.0129637) q[3];
sx q[3];
rz(-2.5240412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.13742927) q[2];
sx q[2];
rz(-0.88976294) q[2];
sx q[2];
rz(1.4271663) q[2];
rz(-3.0754722) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(1.4495513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5071881) q[0];
sx q[0];
rz(-2.2162063) q[0];
sx q[0];
rz(-2.8651644) q[0];
x q[1];
rz(-1.6749803) q[2];
sx q[2];
rz(-1.0717234) q[2];
sx q[2];
rz(0.68945976) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.44504657) q[1];
sx q[1];
rz(-1.3475218) q[1];
sx q[1];
rz(-2.9162507) q[1];
rz(-0.27637847) q[3];
sx q[3];
rz(-2.31782) q[3];
sx q[3];
rz(-2.7793022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9479998) q[2];
sx q[2];
rz(-1.7717382) q[2];
sx q[2];
rz(2.8732079) q[2];
rz(-1.0466446) q[3];
sx q[3];
rz(-2.7768551) q[3];
sx q[3];
rz(-0.31782761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.627581) q[0];
sx q[0];
rz(-1.528897) q[0];
sx q[0];
rz(-2.6348689) q[0];
rz(2.8880033) q[1];
sx q[1];
rz(-1.8702303) q[1];
sx q[1];
rz(-1.0553029) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9297766) q[0];
sx q[0];
rz(-1.388071) q[0];
sx q[0];
rz(-0.94434785) q[0];
rz(-pi) q[1];
rz(0.91674532) q[2];
sx q[2];
rz(-1.9987717) q[2];
sx q[2];
rz(2.5452754) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9195337) q[1];
sx q[1];
rz(-0.93898458) q[1];
sx q[1];
rz(-2.6421089) q[1];
rz(-2.3985732) q[3];
sx q[3];
rz(-1.395874) q[3];
sx q[3];
rz(2.7878441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6598597) q[2];
sx q[2];
rz(-1.0446171) q[2];
sx q[2];
rz(-0.8824904) q[2];
rz(2.4957538) q[3];
sx q[3];
rz(-1.1487938) q[3];
sx q[3];
rz(1.8576436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
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
rz(-2.3690467) q[0];
rz(-1.729471) q[1];
sx q[1];
rz(-1.9344784) q[1];
sx q[1];
rz(-0.57377446) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6360146) q[0];
sx q[0];
rz(-1.6806707) q[0];
sx q[0];
rz(-1.3789603) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7710118) q[2];
sx q[2];
rz(-1.5491312) q[2];
sx q[2];
rz(1.3862762) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7635203) q[1];
sx q[1];
rz(-1.0135279) q[1];
sx q[1];
rz(2.4850363) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8049559) q[3];
sx q[3];
rz(-1.7171211) q[3];
sx q[3];
rz(2.9175266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.039375719) q[2];
sx q[2];
rz(-2.6874459) q[2];
sx q[2];
rz(-2.3708564) q[2];
rz(-0.43631521) q[3];
sx q[3];
rz(-1.8728914) q[3];
sx q[3];
rz(-1.8541981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.27281) q[0];
sx q[0];
rz(-1.2009118) q[0];
sx q[0];
rz(-1.9708721) q[0];
rz(2.6314578) q[1];
sx q[1];
rz(-1.7850103) q[1];
sx q[1];
rz(-1.8458813) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13531317) q[0];
sx q[0];
rz(-1.8989925) q[0];
sx q[0];
rz(0.4853863) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92408085) q[2];
sx q[2];
rz(-0.43281049) q[2];
sx q[2];
rz(-2.2518287) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2669275) q[1];
sx q[1];
rz(-1.3473359) q[1];
sx q[1];
rz(2.2063971) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.027904228) q[3];
sx q[3];
rz(-1.1301665) q[3];
sx q[3];
rz(-0.84907109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7065113) q[2];
sx q[2];
rz(-1.4979829) q[2];
sx q[2];
rz(-2.5734625) q[2];
rz(0.98012296) q[3];
sx q[3];
rz(-1.0390493) q[3];
sx q[3];
rz(0.19395104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4762964) q[0];
sx q[0];
rz(-2.3275573) q[0];
sx q[0];
rz(2.4639159) q[0];
rz(2.9455345) q[1];
sx q[1];
rz(-2.129107) q[1];
sx q[1];
rz(0.83818865) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3121376) q[0];
sx q[0];
rz(-1.2078309) q[0];
sx q[0];
rz(1.8118993) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6423177) q[2];
sx q[2];
rz(-2.3841249) q[2];
sx q[2];
rz(1.0345936) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.14324489) q[1];
sx q[1];
rz(-1.3721826) q[1];
sx q[1];
rz(0.82223383) q[1];
rz(-pi) q[2];
rz(0.30672726) q[3];
sx q[3];
rz(-2.8711651) q[3];
sx q[3];
rz(0.21817792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3395386) q[2];
sx q[2];
rz(-2.4513117) q[2];
sx q[2];
rz(-1.2667123) q[2];
rz(-0.96261111) q[3];
sx q[3];
rz(-1.5714785) q[3];
sx q[3];
rz(-3.0468429) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4203913) q[0];
sx q[0];
rz(-1.2370011) q[0];
sx q[0];
rz(0.23751968) q[0];
rz(1.0182084) q[1];
sx q[1];
rz(-0.84914452) q[1];
sx q[1];
rz(-0.231803) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0203665) q[0];
sx q[0];
rz(-1.0640642) q[0];
sx q[0];
rz(3.0324742) q[0];
rz(-pi) q[1];
x q[1];
rz(2.899029) q[2];
sx q[2];
rz(-1.8055658) q[2];
sx q[2];
rz(0.64005062) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4494891) q[1];
sx q[1];
rz(-2.8398501) q[1];
sx q[1];
rz(2.6471789) q[1];
x q[2];
rz(0.88296367) q[3];
sx q[3];
rz(-1.4045047) q[3];
sx q[3];
rz(-1.7385141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0314177) q[2];
sx q[2];
rz(-1.8877703) q[2];
sx q[2];
rz(1.0882264) q[2];
rz(-0.38816372) q[3];
sx q[3];
rz(-0.66029125) q[3];
sx q[3];
rz(2.3378519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5647472) q[0];
sx q[0];
rz(-1.3544461) q[0];
sx q[0];
rz(2.6690637) q[0];
rz(-2.172773) q[1];
sx q[1];
rz(-0.70822721) q[1];
sx q[1];
rz(0.72475564) q[1];
rz(1.3303403) q[2];
sx q[2];
rz(-1.2546872) q[2];
sx q[2];
rz(1.6457002) q[2];
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