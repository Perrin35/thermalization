OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8213788) q[0];
sx q[0];
rz(7.0452375) q[0];
sx q[0];
rz(10.308916) q[0];
rz(2.60131) q[1];
sx q[1];
rz(1.5958495) q[1];
sx q[1];
rz(11.933239) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54612371) q[0];
sx q[0];
rz(-1.4065388) q[0];
sx q[0];
rz(1.4558653) q[0];
rz(-2.9146637) q[2];
sx q[2];
rz(-2.6440776) q[2];
sx q[2];
rz(1.1253446) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.84304436) q[1];
sx q[1];
rz(-1.6982276) q[1];
sx q[1];
rz(-2.861388) q[1];
rz(-1.9820561) q[3];
sx q[3];
rz(-1.0602128) q[3];
sx q[3];
rz(3.0298281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4151609) q[2];
sx q[2];
rz(-1.8981322) q[2];
sx q[2];
rz(-2.4746056) q[2];
rz(-0.64217448) q[3];
sx q[3];
rz(-0.29688501) q[3];
sx q[3];
rz(1.4750397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17382962) q[0];
sx q[0];
rz(-0.9742631) q[0];
sx q[0];
rz(0.94671384) q[0];
rz(-3.0831378) q[1];
sx q[1];
rz(-2.7775601) q[1];
sx q[1];
rz(-2.2296925) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5240604) q[0];
sx q[0];
rz(-2.723411) q[0];
sx q[0];
rz(-1.6642182) q[0];
rz(-pi) q[1];
rz(1.6591658) q[2];
sx q[2];
rz(-0.67567935) q[2];
sx q[2];
rz(0.32127221) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2649041) q[1];
sx q[1];
rz(-1.1238681) q[1];
sx q[1];
rz(1.8403711) q[1];
rz(-1.5188317) q[3];
sx q[3];
rz(-1.1295737) q[3];
sx q[3];
rz(-0.085746229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.45137206) q[2];
sx q[2];
rz(-0.25190941) q[2];
sx q[2];
rz(2.8042262) q[2];
rz(0.99006027) q[3];
sx q[3];
rz(-1.8127352) q[3];
sx q[3];
rz(-1.2216074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5628691) q[0];
sx q[0];
rz(-0.27414411) q[0];
sx q[0];
rz(-1.8291871) q[0];
rz(-0.21526543) q[1];
sx q[1];
rz(-1.5039517) q[1];
sx q[1];
rz(-0.85830918) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6570411) q[0];
sx q[0];
rz(-0.27052292) q[0];
sx q[0];
rz(1.7596308) q[0];
rz(-pi) q[1];
rz(-2.1444001) q[2];
sx q[2];
rz(-0.39948398) q[2];
sx q[2];
rz(-1.6459207) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.55187622) q[1];
sx q[1];
rz(-1.7360506) q[1];
sx q[1];
rz(0.44323968) q[1];
x q[2];
rz(-0.25562197) q[3];
sx q[3];
rz(-1.2838253) q[3];
sx q[3];
rz(0.14959344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.12576088) q[2];
sx q[2];
rz(-0.55700004) q[2];
sx q[2];
rz(2.1324943) q[2];
rz(-0.5674181) q[3];
sx q[3];
rz(-1.5083454) q[3];
sx q[3];
rz(1.5198038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1384945) q[0];
sx q[0];
rz(-0.69366413) q[0];
sx q[0];
rz(-0.19288119) q[0];
rz(1.0640954) q[1];
sx q[1];
rz(-2.8654983) q[1];
sx q[1];
rz(-0.9562794) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0378936) q[0];
sx q[0];
rz(-2.4913945) q[0];
sx q[0];
rz(-0.29348394) q[0];
rz(-pi) q[1];
rz(0.85727875) q[2];
sx q[2];
rz(-2.3403717) q[2];
sx q[2];
rz(-0.73357936) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6156441) q[1];
sx q[1];
rz(-0.61405776) q[1];
sx q[1];
rz(-0.55974069) q[1];
rz(-2.2364081) q[3];
sx q[3];
rz(-1.3263316) q[3];
sx q[3];
rz(0.28206952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3095779) q[2];
sx q[2];
rz(-0.39065209) q[2];
sx q[2];
rz(-2.2008956) q[2];
rz(-0.37788033) q[3];
sx q[3];
rz(-0.89972073) q[3];
sx q[3];
rz(0.64062029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2864895) q[0];
sx q[0];
rz(-1.6096977) q[0];
sx q[0];
rz(1.8163649) q[0];
rz(-2.8638966) q[1];
sx q[1];
rz(-2.1546021) q[1];
sx q[1];
rz(-1.206254) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2640058) q[0];
sx q[0];
rz(-2.9336399) q[0];
sx q[0];
rz(1.963574) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5304753) q[2];
sx q[2];
rz(-1.9131994) q[2];
sx q[2];
rz(-0.26027203) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5110276) q[1];
sx q[1];
rz(-1.2113809) q[1];
sx q[1];
rz(-3.1267704) q[1];
x q[2];
rz(-1.4021192) q[3];
sx q[3];
rz(-1.3323829) q[3];
sx q[3];
rz(-1.3937529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9557934) q[2];
sx q[2];
rz(-0.15115559) q[2];
sx q[2];
rz(1.7929057) q[2];
rz(-1.0079314) q[3];
sx q[3];
rz(-1.4069087) q[3];
sx q[3];
rz(0.036788363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0808829) q[0];
sx q[0];
rz(-0.01132975) q[0];
sx q[0];
rz(2.9445373) q[0];
rz(-1.9673678) q[1];
sx q[1];
rz(-2.2815727) q[1];
sx q[1];
rz(0.42541447) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2046609) q[0];
sx q[0];
rz(-2.6789894) q[0];
sx q[0];
rz(-3.0499056) q[0];
x q[1];
rz(-1.8749008) q[2];
sx q[2];
rz(-1.5291404) q[2];
sx q[2];
rz(2.0071941) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1716291) q[1];
sx q[1];
rz(-0.49996129) q[1];
sx q[1];
rz(0.5699711) q[1];
rz(-pi) q[2];
rz(0.69611551) q[3];
sx q[3];
rz(-2.1543739) q[3];
sx q[3];
rz(0.35120526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.5328488) q[2];
sx q[2];
rz(-0.11561919) q[2];
sx q[2];
rz(-2.4382584) q[2];
rz(-0.15130791) q[3];
sx q[3];
rz(-1.9688828) q[3];
sx q[3];
rz(0.54733384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5587826) q[0];
sx q[0];
rz(-2.3334239) q[0];
sx q[0];
rz(2.1852344) q[0];
rz(0.76902485) q[1];
sx q[1];
rz(-1.8509357) q[1];
sx q[1];
rz(2.1020611) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7881986) q[0];
sx q[0];
rz(-1.0958584) q[0];
sx q[0];
rz(0.026508412) q[0];
rz(0.86918594) q[2];
sx q[2];
rz(-1.0381881) q[2];
sx q[2];
rz(0.55583304) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7733147) q[1];
sx q[1];
rz(-1.8551143) q[1];
sx q[1];
rz(-1.2523238) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.124362) q[3];
sx q[3];
rz(-0.38700208) q[3];
sx q[3];
rz(-1.6978881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9397554) q[2];
sx q[2];
rz(-2.4079683) q[2];
sx q[2];
rz(-1.5038097) q[2];
rz(-0.87320295) q[3];
sx q[3];
rz(-2.2538897) q[3];
sx q[3];
rz(2.8100815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37209508) q[0];
sx q[0];
rz(-2.7517419) q[0];
sx q[0];
rz(1.1497585) q[0];
rz(0.86504966) q[1];
sx q[1];
rz(-1.9416092) q[1];
sx q[1];
rz(1.4335776) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2835049) q[0];
sx q[0];
rz(-0.69779396) q[0];
sx q[0];
rz(-2.1608864) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9351129) q[2];
sx q[2];
rz(-2.0215567) q[2];
sx q[2];
rz(-0.70012602) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.149274) q[1];
sx q[1];
rz(-2.092835) q[1];
sx q[1];
rz(-3.0401413) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9870342) q[3];
sx q[3];
rz(-1.7572973) q[3];
sx q[3];
rz(2.4881252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4992497) q[2];
sx q[2];
rz(-1.6982102) q[2];
sx q[2];
rz(-0.7439822) q[2];
rz(-2.3104525) q[3];
sx q[3];
rz(-1.4011413) q[3];
sx q[3];
rz(0.35541117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1019854) q[0];
sx q[0];
rz(-2.2672741) q[0];
sx q[0];
rz(-2.4264477) q[0];
rz(-1.8483509) q[1];
sx q[1];
rz(-1.5903541) q[1];
sx q[1];
rz(-2.901758) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6699246) q[0];
sx q[0];
rz(-1.0201766) q[0];
sx q[0];
rz(-2.1651833) q[0];
rz(-pi) q[1];
rz(-1.6314028) q[2];
sx q[2];
rz(-2.1327634) q[2];
sx q[2];
rz(0.89982239) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.72425408) q[1];
sx q[1];
rz(-2.2431233) q[1];
sx q[1];
rz(0.21604854) q[1];
rz(-0.49107869) q[3];
sx q[3];
rz(-1.5415191) q[3];
sx q[3];
rz(1.0539583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0750601) q[2];
sx q[2];
rz(-1.5873453) q[2];
sx q[2];
rz(-1.6481605) q[2];
rz(-3.1028124) q[3];
sx q[3];
rz(-1.6438234) q[3];
sx q[3];
rz(-1.7356167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1455602) q[0];
sx q[0];
rz(-1.6479489) q[0];
sx q[0];
rz(-0.62553585) q[0];
rz(-1.372555) q[1];
sx q[1];
rz(-0.99634975) q[1];
sx q[1];
rz(0.58403429) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9720219) q[0];
sx q[0];
rz(-1.2002647) q[0];
sx q[0];
rz(-1.2967923) q[0];
rz(-pi) q[1];
rz(2.4768684) q[2];
sx q[2];
rz(-1.692284) q[2];
sx q[2];
rz(-2.6912376) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6713114) q[1];
sx q[1];
rz(-1.352942) q[1];
sx q[1];
rz(2.4163428) q[1];
x q[2];
rz(0.32665334) q[3];
sx q[3];
rz(-2.4295837) q[3];
sx q[3];
rz(-0.53381843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.11223101) q[2];
sx q[2];
rz(-1.6718622) q[2];
sx q[2];
rz(-1.7870129) q[2];
rz(0.77238798) q[3];
sx q[3];
rz(-1.9410917) q[3];
sx q[3];
rz(-1.6105917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1188485) q[0];
sx q[0];
rz(-2.4130555) q[0];
sx q[0];
rz(1.5424315) q[0];
rz(0.45730418) q[1];
sx q[1];
rz(-0.87304579) q[1];
sx q[1];
rz(2.9541448) q[1];
rz(3.0627261) q[2];
sx q[2];
rz(-1.7577684) q[2];
sx q[2];
rz(-2.656542) q[2];
rz(2.0304092) q[3];
sx q[3];
rz(-1.8596953) q[3];
sx q[3];
rz(0.019712899) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
