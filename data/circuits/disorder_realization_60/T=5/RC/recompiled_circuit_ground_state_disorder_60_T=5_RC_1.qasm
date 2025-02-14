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
rz(-2.3795405) q[0];
sx q[0];
rz(2.2574545) q[0];
rz(-0.54028264) q[1];
sx q[1];
rz(-1.5958495) q[1];
sx q[1];
rz(-0.63313142) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0057982) q[0];
sx q[0];
rz(-1.6841736) q[0];
sx q[0];
rz(-2.9762639) q[0];
rz(-pi) q[1];
rz(0.48670003) q[2];
sx q[2];
rz(-1.6783769) q[2];
sx q[2];
rz(-0.24525596) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.76430815) q[1];
sx q[1];
rz(-1.2929243) q[1];
sx q[1];
rz(-1.4382526) q[1];
rz(-pi) q[2];
rz(-0.61986519) q[3];
sx q[3];
rz(-2.4976118) q[3];
sx q[3];
rz(-2.5247273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4151609) q[2];
sx q[2];
rz(-1.2434604) q[2];
sx q[2];
rz(2.4746056) q[2];
rz(-0.64217448) q[3];
sx q[3];
rz(-2.8447076) q[3];
sx q[3];
rz(-1.4750397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.967763) q[0];
sx q[0];
rz(-0.9742631) q[0];
sx q[0];
rz(-2.1948788) q[0];
rz(0.058454839) q[1];
sx q[1];
rz(-0.36403257) q[1];
sx q[1];
rz(2.2296925) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13214853) q[0];
sx q[0];
rz(-1.6086888) q[0];
sx q[0];
rz(1.9873585) q[0];
rz(-2.2445685) q[2];
sx q[2];
rz(-1.5155715) q[2];
sx q[2];
rz(1.3185475) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.30787779) q[1];
sx q[1];
rz(-0.51719292) q[1];
sx q[1];
rz(2.6343995) q[1];
rz(2.6998482) q[3];
sx q[3];
rz(-1.5238122) q[3];
sx q[3];
rz(-1.6787501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.45137206) q[2];
sx q[2];
rz(-0.25190941) q[2];
sx q[2];
rz(-2.8042262) q[2];
rz(-2.1515324) q[3];
sx q[3];
rz(-1.8127352) q[3];
sx q[3];
rz(1.9199853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5628691) q[0];
sx q[0];
rz(-0.27414411) q[0];
sx q[0];
rz(1.3124056) q[0];
rz(2.9263272) q[1];
sx q[1];
rz(-1.637641) q[1];
sx q[1];
rz(-2.2832835) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68033531) q[0];
sx q[0];
rz(-1.3051998) q[0];
sx q[0];
rz(-3.0895825) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1444001) q[2];
sx q[2];
rz(-2.7421087) q[2];
sx q[2];
rz(-1.6459207) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.55187622) q[1];
sx q[1];
rz(-1.7360506) q[1];
sx q[1];
rz(-2.698353) q[1];
rz(-pi) q[2];
rz(-1.8668601) q[3];
sx q[3];
rz(-1.8157457) q[3];
sx q[3];
rz(-1.3473657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0158318) q[2];
sx q[2];
rz(-0.55700004) q[2];
sx q[2];
rz(1.0090984) q[2];
rz(-0.5674181) q[3];
sx q[3];
rz(-1.5083454) q[3];
sx q[3];
rz(-1.6217888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0030982) q[0];
sx q[0];
rz(-2.4479285) q[0];
sx q[0];
rz(-0.19288119) q[0];
rz(-1.0640954) q[1];
sx q[1];
rz(-0.27609438) q[1];
sx q[1];
rz(-0.9562794) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.372626) q[0];
sx q[0];
rz(-1.7468233) q[0];
sx q[0];
rz(2.5123216) q[0];
x q[1];
rz(2.2334571) q[2];
sx q[2];
rz(-1.081433) q[2];
sx q[2];
rz(0.2950677) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2697105) q[1];
sx q[1];
rz(-2.080889) q[1];
sx q[1];
rz(-1.2126232) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2364081) q[3];
sx q[3];
rz(-1.3263316) q[3];
sx q[3];
rz(2.8595231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8320147) q[2];
sx q[2];
rz(-2.7509406) q[2];
sx q[2];
rz(2.2008956) q[2];
rz(-2.7637123) q[3];
sx q[3];
rz(-0.89972073) q[3];
sx q[3];
rz(-0.64062029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
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
rz(-0.9869906) q[1];
sx q[1];
rz(-1.9353386) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6645443) q[0];
sx q[0];
rz(-1.7627075) q[0];
sx q[0];
rz(0.080587793) q[0];
rz(-pi) q[1];
rz(1.1603068) q[2];
sx q[2];
rz(-2.1417981) q[2];
sx q[2];
rz(2.0620907) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5110276) q[1];
sx q[1];
rz(-1.2113809) q[1];
sx q[1];
rz(3.1267704) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24171452) q[3];
sx q[3];
rz(-1.7346584) q[3];
sx q[3];
rz(0.13684798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9557934) q[2];
sx q[2];
rz(-0.15115559) q[2];
sx q[2];
rz(1.348687) q[2];
rz(2.1336613) q[3];
sx q[3];
rz(-1.4069087) q[3];
sx q[3];
rz(0.036788363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0808829) q[0];
sx q[0];
rz(-0.01132975) q[0];
sx q[0];
rz(-2.9445373) q[0];
rz(-1.1742249) q[1];
sx q[1];
rz(-2.2815727) q[1];
sx q[1];
rz(-0.42541447) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0393164) q[0];
sx q[0];
rz(-1.1102866) q[0];
sx q[0];
rz(1.6164245) q[0];
rz(-pi) q[1];
rz(-1.8749008) q[2];
sx q[2];
rz(-1.6124523) q[2];
sx q[2];
rz(-2.0071941) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1716291) q[1];
sx q[1];
rz(-0.49996129) q[1];
sx q[1];
rz(2.5716216) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3415575) q[3];
sx q[3];
rz(-0.87585319) q[3];
sx q[3];
rz(-1.3390697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5328488) q[2];
sx q[2];
rz(-0.11561919) q[2];
sx q[2];
rz(2.4382584) q[2];
rz(-2.9902847) q[3];
sx q[3];
rz(-1.9688828) q[3];
sx q[3];
rz(2.5942588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5828101) q[0];
sx q[0];
rz(-0.8081688) q[0];
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
rz(-1.5943702) q[0];
sx q[0];
rz(2.0458771) q[0];
x q[1];
rz(-0.83073406) q[2];
sx q[2];
rz(-0.85262876) q[2];
sx q[2];
rz(1.5560395) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6440297) q[1];
sx q[1];
rz(-2.7179238) q[1];
sx q[1];
rz(-0.81987857) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9045181) q[3];
sx q[3];
rz(-1.7705373) q[3];
sx q[3];
rz(-0.39272308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2018373) q[2];
sx q[2];
rz(-0.73362437) q[2];
sx q[2];
rz(1.5038097) q[2];
rz(-2.2683897) q[3];
sx q[3];
rz(-2.2538897) q[3];
sx q[3];
rz(-2.8100815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37209508) q[0];
sx q[0];
rz(-0.3898507) q[0];
sx q[0];
rz(-1.1497585) q[0];
rz(2.276543) q[1];
sx q[1];
rz(-1.9416092) q[1];
sx q[1];
rz(1.7080151) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85808774) q[0];
sx q[0];
rz(-2.4437987) q[0];
sx q[0];
rz(2.1608864) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2064798) q[2];
sx q[2];
rz(-2.0215567) q[2];
sx q[2];
rz(-2.4414666) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5123957) q[1];
sx q[1];
rz(-1.6586972) q[1];
sx q[1];
rz(-2.0950661) q[1];
rz(-1.1545584) q[3];
sx q[3];
rz(-1.7572973) q[3];
sx q[3];
rz(0.65346741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4992497) q[2];
sx q[2];
rz(-1.6982102) q[2];
sx q[2];
rz(-2.3976105) q[2];
rz(-0.83114019) q[3];
sx q[3];
rz(-1.7404514) q[3];
sx q[3];
rz(0.35541117) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0396073) q[0];
sx q[0];
rz(-2.2672741) q[0];
sx q[0];
rz(0.71514493) q[0];
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
rz(-2.6699246) q[0];
sx q[0];
rz(-2.1214161) q[0];
sx q[0];
rz(2.1651833) q[0];
x q[1];
rz(-1.6314028) q[2];
sx q[2];
rz(-1.0088293) q[2];
sx q[2];
rz(-0.89982239) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.72425408) q[1];
sx q[1];
rz(-2.2431233) q[1];
sx q[1];
rz(-2.9255441) q[1];
rz(2.650514) q[3];
sx q[3];
rz(-1.5415191) q[3];
sx q[3];
rz(1.0539583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0750601) q[2];
sx q[2];
rz(-1.5542474) q[2];
sx q[2];
rz(-1.6481605) q[2];
rz(-3.1028124) q[3];
sx q[3];
rz(-1.6438234) q[3];
sx q[3];
rz(1.405976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99603242) q[0];
sx q[0];
rz(-1.4936438) q[0];
sx q[0];
rz(0.62553585) q[0];
rz(1.7690376) q[1];
sx q[1];
rz(-0.99634975) q[1];
sx q[1];
rz(0.58403429) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16957078) q[0];
sx q[0];
rz(-1.2002647) q[0];
sx q[0];
rz(1.2967923) q[0];
x q[1];
rz(-1.7246849) q[2];
sx q[2];
rz(-0.91183543) q[2];
sx q[2];
rz(1.9264592) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.088773397) q[1];
sx q[1];
rz(-0.86630834) q[1];
sx q[1];
rz(1.2831844) q[1];
rz(-pi) q[2];
rz(1.8409506) q[3];
sx q[3];
rz(-0.90357257) q[3];
sx q[3];
rz(3.0285578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0293616) q[2];
sx q[2];
rz(-1.4697305) q[2];
sx q[2];
rz(1.3545797) q[2];
rz(-2.3692047) q[3];
sx q[3];
rz(-1.200501) q[3];
sx q[3];
rz(1.6105917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0227441) q[0];
sx q[0];
rz(-2.4130555) q[0];
sx q[0];
rz(1.5424315) q[0];
rz(2.6842885) q[1];
sx q[1];
rz(-2.2685469) q[1];
sx q[1];
rz(-0.1874478) q[1];
rz(-1.3832548) q[2];
sx q[2];
rz(-1.6482856) q[2];
sx q[2];
rz(2.0411574) q[2];
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
