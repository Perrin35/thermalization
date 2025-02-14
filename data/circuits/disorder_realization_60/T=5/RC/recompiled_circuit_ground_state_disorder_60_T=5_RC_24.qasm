OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3202138) q[0];
sx q[0];
rz(-0.76205215) q[0];
sx q[0];
rz(-2.2574545) q[0];
rz(2.60131) q[1];
sx q[1];
rz(-1.5457431) q[1];
sx q[1];
rz(0.63313142) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0057982) q[0];
sx q[0];
rz(-1.6841736) q[0];
sx q[0];
rz(0.16532872) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9146637) q[2];
sx q[2];
rz(-0.49751505) q[2];
sx q[2];
rz(2.016248) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2985483) q[1];
sx q[1];
rz(-1.4433651) q[1];
sx q[1];
rz(-2.861388) q[1];
rz(-pi) q[2];
rz(-1.1595366) q[3];
sx q[3];
rz(-1.0602128) q[3];
sx q[3];
rz(0.11176457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7264318) q[2];
sx q[2];
rz(-1.8981322) q[2];
sx q[2];
rz(-2.4746056) q[2];
rz(-2.4994182) q[3];
sx q[3];
rz(-2.8447076) q[3];
sx q[3];
rz(1.4750397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.967763) q[0];
sx q[0];
rz(-0.9742631) q[0];
sx q[0];
rz(-0.94671384) q[0];
rz(0.058454839) q[1];
sx q[1];
rz(-2.7775601) q[1];
sx q[1];
rz(0.91190016) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4218876) q[0];
sx q[0];
rz(-1.1545517) q[0];
sx q[0];
rz(-3.1001607) q[0];
rz(-pi) q[1];
rz(2.2445685) q[2];
sx q[2];
rz(-1.5155715) q[2];
sx q[2];
rz(1.8230452) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.30787779) q[1];
sx q[1];
rz(-0.51719292) q[1];
sx q[1];
rz(-0.50719313) q[1];
x q[2];
rz(-0.1095406) q[3];
sx q[3];
rz(-2.6975204) q[3];
sx q[3];
rz(0.20694297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.45137206) q[2];
sx q[2];
rz(-2.8896832) q[2];
sx q[2];
rz(0.33736649) q[2];
rz(-0.99006027) q[3];
sx q[3];
rz(-1.3288574) q[3];
sx q[3];
rz(1.9199853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57872355) q[0];
sx q[0];
rz(-2.8674485) q[0];
sx q[0];
rz(1.8291871) q[0];
rz(-0.21526543) q[1];
sx q[1];
rz(-1.637641) q[1];
sx q[1];
rz(0.85830918) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4612573) q[0];
sx q[0];
rz(-1.8363928) q[0];
sx q[0];
rz(-3.0895825) q[0];
rz(-pi) q[1];
rz(-0.99719259) q[2];
sx q[2];
rz(-2.7421087) q[2];
sx q[2];
rz(-1.6459207) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.94098136) q[1];
sx q[1];
rz(-2.0075782) q[1];
sx q[1];
rz(1.753356) q[1];
rz(-pi) q[2];
rz(1.8668601) q[3];
sx q[3];
rz(-1.325847) q[3];
sx q[3];
rz(-1.3473657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.12576088) q[2];
sx q[2];
rz(-2.5845926) q[2];
sx q[2];
rz(-1.0090984) q[2];
rz(-0.5674181) q[3];
sx q[3];
rz(-1.6332473) q[3];
sx q[3];
rz(-1.5198038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1384945) q[0];
sx q[0];
rz(-0.69366413) q[0];
sx q[0];
rz(2.9487115) q[0];
rz(-1.0640954) q[1];
sx q[1];
rz(-2.8654983) q[1];
sx q[1];
rz(0.9562794) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76896669) q[0];
sx q[0];
rz(-1.7468233) q[0];
sx q[0];
rz(-2.5123216) q[0];
rz(2.2334571) q[2];
sx q[2];
rz(-2.0601597) q[2];
sx q[2];
rz(2.846525) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.51814521) q[1];
sx q[1];
rz(-1.259874) q[1];
sx q[1];
rz(2.6030933) q[1];
x q[2];
rz(-1.1868915) q[3];
sx q[3];
rz(-2.4389746) q[3];
sx q[3];
rz(-1.5877569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3095779) q[2];
sx q[2];
rz(-0.39065209) q[2];
sx q[2];
rz(0.94069707) q[2];
rz(-0.37788033) q[3];
sx q[3];
rz(-2.2418719) q[3];
sx q[3];
rz(2.5009724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85510319) q[0];
sx q[0];
rz(-1.6096977) q[0];
sx q[0];
rz(1.3252277) q[0];
rz(-0.27769604) q[1];
sx q[1];
rz(-2.1546021) q[1];
sx q[1];
rz(1.206254) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2640058) q[0];
sx q[0];
rz(-2.9336399) q[0];
sx q[0];
rz(1.963574) q[0];
rz(-pi) q[1];
rz(-2.5857193) q[2];
sx q[2];
rz(-2.4519359) q[2];
sx q[2];
rz(1.3841615) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6305651) q[1];
sx q[1];
rz(-1.9302118) q[1];
sx q[1];
rz(0.014822249) q[1];
x q[2];
rz(-2.8998781) q[3];
sx q[3];
rz(-1.4069342) q[3];
sx q[3];
rz(3.0047447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.18579927) q[2];
sx q[2];
rz(-0.15115559) q[2];
sx q[2];
rz(1.7929057) q[2];
rz(1.0079314) q[3];
sx q[3];
rz(-1.734684) q[3];
sx q[3];
rz(-3.1048043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0808829) q[0];
sx q[0];
rz(-3.1302629) q[0];
sx q[0];
rz(-2.9445373) q[0];
rz(1.1742249) q[1];
sx q[1];
rz(-0.86001992) q[1];
sx q[1];
rz(-0.42541447) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1022762) q[0];
sx q[0];
rz(-1.1102866) q[0];
sx q[0];
rz(-1.5251681) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0979359) q[2];
sx q[2];
rz(-1.8746286) q[2];
sx q[2];
rz(-0.44946656) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1716291) q[1];
sx q[1];
rz(-2.6416314) q[1];
sx q[1];
rz(-0.5699711) q[1];
rz(-pi) q[2];
rz(0.80003512) q[3];
sx q[3];
rz(-2.2657395) q[3];
sx q[3];
rz(-1.8025229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5328488) q[2];
sx q[2];
rz(-3.0259735) q[2];
sx q[2];
rz(-2.4382584) q[2];
rz(2.9902847) q[3];
sx q[3];
rz(-1.1727099) q[3];
sx q[3];
rz(2.5942588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5587826) q[0];
sx q[0];
rz(-0.8081688) q[0];
sx q[0];
rz(2.1852344) q[0];
rz(-2.3725678) q[1];
sx q[1];
rz(-1.8509357) q[1];
sx q[1];
rz(-1.0395315) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9120662) q[0];
sx q[0];
rz(-1.5943702) q[0];
sx q[0];
rz(2.0458771) q[0];
x q[1];
rz(2.4843487) q[2];
sx q[2];
rz(-2.1604156) q[2];
sx q[2];
rz(-0.60962617) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8468561) q[1];
sx q[1];
rz(-1.2655317) q[1];
sx q[1];
rz(-0.29851361) q[1];
rz(2.9305228) q[3];
sx q[3];
rz(-1.8976334) q[3];
sx q[3];
rz(-1.1093931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2018373) q[2];
sx q[2];
rz(-0.73362437) q[2];
sx q[2];
rz(-1.5038097) q[2];
rz(2.2683897) q[3];
sx q[3];
rz(-2.2538897) q[3];
sx q[3];
rz(2.8100815) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7694976) q[0];
sx q[0];
rz(-0.3898507) q[0];
sx q[0];
rz(-1.1497585) q[0];
rz(-2.276543) q[1];
sx q[1];
rz(-1.9416092) q[1];
sx q[1];
rz(-1.7080151) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85808774) q[0];
sx q[0];
rz(-0.69779396) q[0];
sx q[0];
rz(2.1608864) q[0];
rz(-pi) q[1];
rz(1.9351129) q[2];
sx q[2];
rz(-2.0215567) q[2];
sx q[2];
rz(-2.4414666) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5123957) q[1];
sx q[1];
rz(-1.4828955) q[1];
sx q[1];
rz(2.0950661) q[1];
x q[2];
rz(-1.9870342) q[3];
sx q[3];
rz(-1.7572973) q[3];
sx q[3];
rz(-0.65346741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.64234298) q[2];
sx q[2];
rz(-1.6982102) q[2];
sx q[2];
rz(0.7439822) q[2];
rz(0.83114019) q[3];
sx q[3];
rz(-1.4011413) q[3];
sx q[3];
rz(-2.7861815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0396073) q[0];
sx q[0];
rz(-2.2672741) q[0];
sx q[0];
rz(-0.71514493) q[0];
rz(1.8483509) q[1];
sx q[1];
rz(-1.5903541) q[1];
sx q[1];
rz(-0.23983461) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6699246) q[0];
sx q[0];
rz(-1.0201766) q[0];
sx q[0];
rz(-2.1651833) q[0];
x q[1];
rz(0.56279601) q[2];
sx q[2];
rz(-1.6220731) q[2];
sx q[2];
rz(2.4382961) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0630602) q[1];
sx q[1];
rz(-2.4405747) q[1];
sx q[1];
rz(-1.8338507) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6039944) q[3];
sx q[3];
rz(-2.0616459) q[3];
sx q[3];
rz(-2.6404078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0665325) q[2];
sx q[2];
rz(-1.5873453) q[2];
sx q[2];
rz(1.4934322) q[2];
rz(-0.0387803) q[3];
sx q[3];
rz(-1.4977692) q[3];
sx q[3];
rz(1.405976) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1455602) q[0];
sx q[0];
rz(-1.4936438) q[0];
sx q[0];
rz(-2.5160568) q[0];
rz(-1.7690376) q[1];
sx q[1];
rz(-0.99634975) q[1];
sx q[1];
rz(-0.58403429) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9720219) q[0];
sx q[0];
rz(-1.941328) q[0];
sx q[0];
rz(1.2967923) q[0];
x q[1];
rz(2.9461924) q[2];
sx q[2];
rz(-0.6740734) q[2];
sx q[2];
rz(0.96701996) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6713114) q[1];
sx q[1];
rz(-1.352942) q[1];
sx q[1];
rz(-2.4163428) q[1];
rz(-pi) q[2];
rz(-2.4563362) q[3];
sx q[3];
rz(-1.3595837) q[3];
sx q[3];
rz(1.8535456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0293616) q[2];
sx q[2];
rz(-1.6718622) q[2];
sx q[2];
rz(-1.3545797) q[2];
rz(-0.77238798) q[3];
sx q[3];
rz(-1.9410917) q[3];
sx q[3];
rz(1.6105917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1188485) q[0];
sx q[0];
rz(-0.72853715) q[0];
sx q[0];
rz(-1.5991612) q[0];
rz(-0.45730418) q[1];
sx q[1];
rz(-2.2685469) q[1];
sx q[1];
rz(-0.1874478) q[1];
rz(1.3832548) q[2];
sx q[2];
rz(-1.493307) q[2];
sx q[2];
rz(-1.1004352) q[2];
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
