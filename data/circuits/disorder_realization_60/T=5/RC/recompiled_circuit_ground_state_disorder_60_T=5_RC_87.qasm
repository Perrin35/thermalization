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
rz(-0.54028264) q[1];
sx q[1];
rz(-1.5958495) q[1];
sx q[1];
rz(-0.63313142) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0057982) q[0];
sx q[0];
rz(-1.6841736) q[0];
sx q[0];
rz(2.9762639) q[0];
x q[1];
rz(-0.22692899) q[2];
sx q[2];
rz(-0.49751505) q[2];
sx q[2];
rz(1.1253446) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76430815) q[1];
sx q[1];
rz(-1.8486684) q[1];
sx q[1];
rz(1.7033401) q[1];
rz(1.1595366) q[3];
sx q[3];
rz(-1.0602128) q[3];
sx q[3];
rz(3.0298281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7264318) q[2];
sx q[2];
rz(-1.2434604) q[2];
sx q[2];
rz(2.4746056) q[2];
rz(2.4994182) q[3];
sx q[3];
rz(-2.8447076) q[3];
sx q[3];
rz(1.6665529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17382962) q[0];
sx q[0];
rz(-0.9742631) q[0];
sx q[0];
rz(2.1948788) q[0];
rz(3.0831378) q[1];
sx q[1];
rz(-2.7775601) q[1];
sx q[1];
rz(2.2296925) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0094441) q[0];
sx q[0];
rz(-1.5329038) q[0];
sx q[0];
rz(-1.1542341) q[0];
rz(-pi) q[1];
rz(-0.070621892) q[2];
sx q[2];
rz(-0.89824072) q[2];
sx q[2];
rz(0.20820752) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3286314) q[1];
sx q[1];
rz(-1.3282623) q[1];
sx q[1];
rz(0.46142908) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44174444) q[3];
sx q[3];
rz(-1.6177804) q[3];
sx q[3];
rz(1.4628425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6902206) q[2];
sx q[2];
rz(-2.8896832) q[2];
sx q[2];
rz(2.8042262) q[2];
rz(0.99006027) q[3];
sx q[3];
rz(-1.3288574) q[3];
sx q[3];
rz(1.2216074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5628691) q[0];
sx q[0];
rz(-0.27414411) q[0];
sx q[0];
rz(1.8291871) q[0];
rz(2.9263272) q[1];
sx q[1];
rz(-1.5039517) q[1];
sx q[1];
rz(2.2832835) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2374683) q[0];
sx q[0];
rz(-1.6209813) q[0];
sx q[0];
rz(1.8367358) q[0];
rz(-0.22521727) q[2];
sx q[2];
rz(-1.2379939) q[2];
sx q[2];
rz(1.0343347) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5897164) q[1];
sx q[1];
rz(-1.405542) q[1];
sx q[1];
rz(-2.698353) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8859707) q[3];
sx q[3];
rz(-1.2838253) q[3];
sx q[3];
rz(2.9919992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.12576088) q[2];
sx q[2];
rz(-0.55700004) q[2];
sx q[2];
rz(-1.0090984) q[2];
rz(-2.5741746) q[3];
sx q[3];
rz(-1.5083454) q[3];
sx q[3];
rz(-1.5198038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1384945) q[0];
sx q[0];
rz(-0.69366413) q[0];
sx q[0];
rz(-2.9487115) q[0];
rz(-2.0774972) q[1];
sx q[1];
rz(-2.8654983) q[1];
sx q[1];
rz(-0.9562794) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0378936) q[0];
sx q[0];
rz(-0.65019817) q[0];
sx q[0];
rz(2.8481087) q[0];
rz(-pi) q[1];
x q[1];
rz(0.85727875) q[2];
sx q[2];
rz(-0.801221) q[2];
sx q[2];
rz(-2.4080133) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6156441) q[1];
sx q[1];
rz(-0.61405776) q[1];
sx q[1];
rz(2.581852) q[1];
rz(-pi) q[2];
rz(2.2364081) q[3];
sx q[3];
rz(-1.3263316) q[3];
sx q[3];
rz(2.8595231) q[3];
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
rz(0.37788033) q[3];
sx q[3];
rz(-2.2418719) q[3];
sx q[3];
rz(-2.5009724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2864895) q[0];
sx q[0];
rz(-1.531895) q[0];
sx q[0];
rz(-1.3252277) q[0];
rz(0.27769604) q[1];
sx q[1];
rz(-2.1546021) q[1];
sx q[1];
rz(1.9353386) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4770484) q[0];
sx q[0];
rz(-1.7627075) q[0];
sx q[0];
rz(-0.080587793) q[0];
x q[1];
rz(-0.61111738) q[2];
sx q[2];
rz(-1.9131994) q[2];
sx q[2];
rz(2.8813206) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5531471) q[1];
sx q[1];
rz(-2.781885) q[1];
sx q[1];
rz(-1.5313696) q[1];
x q[2];
rz(0.24171452) q[3];
sx q[3];
rz(-1.4069342) q[3];
sx q[3];
rz(3.0047447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9557934) q[2];
sx q[2];
rz(-0.15115559) q[2];
sx q[2];
rz(1.348687) q[2];
rz(1.0079314) q[3];
sx q[3];
rz(-1.734684) q[3];
sx q[3];
rz(-3.1048043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0808829) q[0];
sx q[0];
rz(-0.01132975) q[0];
sx q[0];
rz(-0.19705535) q[0];
rz(1.9673678) q[1];
sx q[1];
rz(-0.86001992) q[1];
sx q[1];
rz(-2.7161782) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4482314) q[0];
sx q[0];
rz(-1.5299242) q[0];
sx q[0];
rz(-0.46092439) q[0];
rz(-pi) q[1];
rz(1.4324911) q[2];
sx q[2];
rz(-2.8347361) q[2];
sx q[2];
rz(-0.30447659) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1716291) q[1];
sx q[1];
rz(-0.49996129) q[1];
sx q[1];
rz(-2.5716216) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4454771) q[3];
sx q[3];
rz(-0.98721877) q[3];
sx q[3];
rz(-0.35120526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.60874385) q[2];
sx q[2];
rz(-0.11561919) q[2];
sx q[2];
rz(2.4382584) q[2];
rz(-2.9902847) q[3];
sx q[3];
rz(-1.9688828) q[3];
sx q[3];
rz(-0.54733384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5828101) q[0];
sx q[0];
rz(-2.3334239) q[0];
sx q[0];
rz(-0.95635828) q[0];
rz(-2.3725678) q[1];
sx q[1];
rz(-1.8509357) q[1];
sx q[1];
rz(-1.0395315) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2295264) q[0];
sx q[0];
rz(-1.5472224) q[0];
sx q[0];
rz(1.0957155) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4843487) q[2];
sx q[2];
rz(-0.98117706) q[2];
sx q[2];
rz(-0.60962617) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.49756296) q[1];
sx q[1];
rz(-2.7179238) q[1];
sx q[1];
rz(-2.3217141) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9305228) q[3];
sx q[3];
rz(-1.8976334) q[3];
sx q[3];
rz(2.0321996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2018373) q[2];
sx q[2];
rz(-0.73362437) q[2];
sx q[2];
rz(1.5038097) q[2];
rz(2.2683897) q[3];
sx q[3];
rz(-2.2538897) q[3];
sx q[3];
rz(-0.33151117) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37209508) q[0];
sx q[0];
rz(-0.3898507) q[0];
sx q[0];
rz(1.1497585) q[0];
rz(0.86504966) q[1];
sx q[1];
rz(-1.1999835) q[1];
sx q[1];
rz(1.7080151) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5763315) q[0];
sx q[0];
rz(-1.0076242) q[0];
sx q[0];
rz(2.7050326) q[0];
rz(-pi) q[1];
rz(0.63460372) q[2];
sx q[2];
rz(-0.57159492) q[2];
sx q[2];
rz(-1.7224479) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.62919697) q[1];
sx q[1];
rz(-1.6586972) q[1];
sx q[1];
rz(2.0950661) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0074437) q[3];
sx q[3];
rz(-2.6877207) q[3];
sx q[3];
rz(-1.8271104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4992497) q[2];
sx q[2];
rz(-1.4433824) q[2];
sx q[2];
rz(2.3976105) q[2];
rz(-2.3104525) q[3];
sx q[3];
rz(-1.4011413) q[3];
sx q[3];
rz(-2.7861815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0396073) q[0];
sx q[0];
rz(-2.2672741) q[0];
sx q[0];
rz(-2.4264477) q[0];
rz(1.2932418) q[1];
sx q[1];
rz(-1.5903541) q[1];
sx q[1];
rz(-2.901758) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47166804) q[0];
sx q[0];
rz(-2.1214161) q[0];
sx q[0];
rz(2.1651833) q[0];
rz(-pi) q[1];
rz(-1.6314028) q[2];
sx q[2];
rz(-1.0088293) q[2];
sx q[2];
rz(2.2417703) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4173386) q[1];
sx q[1];
rz(-2.2431233) q[1];
sx q[1];
rz(-0.21604854) q[1];
rz(-2.650514) q[3];
sx q[3];
rz(-1.6000736) q[3];
sx q[3];
rz(-2.0876344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0665325) q[2];
sx q[2];
rz(-1.5542474) q[2];
sx q[2];
rz(-1.6481605) q[2];
rz(-3.1028124) q[3];
sx q[3];
rz(-1.4977692) q[3];
sx q[3];
rz(-1.405976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99603242) q[0];
sx q[0];
rz(-1.6479489) q[0];
sx q[0];
rz(-2.5160568) q[0];
rz(-1.7690376) q[1];
sx q[1];
rz(-0.99634975) q[1];
sx q[1];
rz(-0.58403429) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3119572) q[0];
sx q[0];
rz(-2.6845508) q[0];
sx q[0];
rz(-2.5331927) q[0];
rz(1.7246849) q[2];
sx q[2];
rz(-2.2297572) q[2];
sx q[2];
rz(1.9264592) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6713114) q[1];
sx q[1];
rz(-1.7886506) q[1];
sx q[1];
rz(-0.72524986) q[1];
rz(-0.68525641) q[3];
sx q[3];
rz(-1.7820089) q[3];
sx q[3];
rz(1.8535456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.11223101) q[2];
sx q[2];
rz(-1.6718622) q[2];
sx q[2];
rz(1.3545797) q[2];
rz(2.3692047) q[3];
sx q[3];
rz(-1.9410917) q[3];
sx q[3];
rz(-1.5310009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0227441) q[0];
sx q[0];
rz(-0.72853715) q[0];
sx q[0];
rz(-1.5991612) q[0];
rz(0.45730418) q[1];
sx q[1];
rz(-0.87304579) q[1];
sx q[1];
rz(2.9541448) q[1];
rz(1.3832548) q[2];
sx q[2];
rz(-1.493307) q[2];
sx q[2];
rz(-1.1004352) q[2];
rz(-0.32021602) q[3];
sx q[3];
rz(-2.0099985) q[3];
sx q[3];
rz(1.7306001) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
