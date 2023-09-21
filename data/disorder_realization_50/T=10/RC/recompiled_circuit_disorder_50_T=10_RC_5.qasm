OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2621736) q[0];
sx q[0];
rz(-1.7466495) q[0];
sx q[0];
rz(3.1403132) q[0];
rz(1.4446422) q[1];
sx q[1];
rz(-1.1029707) q[1];
sx q[1];
rz(-0.77494088) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9026731) q[0];
sx q[0];
rz(-0.15107778) q[0];
sx q[0];
rz(2.8047049) q[0];
x q[1];
rz(0.84471976) q[2];
sx q[2];
rz(-1.8978999) q[2];
sx q[2];
rz(-2.6127882) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.70667627) q[1];
sx q[1];
rz(-1.4286563) q[1];
sx q[1];
rz(-3.1013156) q[1];
rz(-0.28489057) q[3];
sx q[3];
rz(-1.4244411) q[3];
sx q[3];
rz(0.26998664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0960192) q[2];
sx q[2];
rz(-1.1117659) q[2];
sx q[2];
rz(1.9457031) q[2];
rz(-1.1536417) q[3];
sx q[3];
rz(-2.3524645) q[3];
sx q[3];
rz(-1.7529863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4202704) q[0];
sx q[0];
rz(-0.35636154) q[0];
sx q[0];
rz(1.8700245) q[0];
rz(-2.0416073) q[1];
sx q[1];
rz(-1.1106691) q[1];
sx q[1];
rz(-1.3756479) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.401706) q[0];
sx q[0];
rz(-1.4835843) q[0];
sx q[0];
rz(2.7313822) q[0];
rz(-pi) q[1];
rz(-3.0963418) q[2];
sx q[2];
rz(-1.3327193) q[2];
sx q[2];
rz(-2.0766052) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.980629) q[1];
sx q[1];
rz(-1.3914319) q[1];
sx q[1];
rz(-0.37100002) q[1];
x q[2];
rz(-0.35956412) q[3];
sx q[3];
rz(-1.0457977) q[3];
sx q[3];
rz(-1.1121225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9282844) q[2];
sx q[2];
rz(-1.3457315) q[2];
sx q[2];
rz(-0.48970547) q[2];
rz(2.0080163) q[3];
sx q[3];
rz(-2.9462892) q[3];
sx q[3];
rz(-1.0666696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5415444) q[0];
sx q[0];
rz(-1.2788037) q[0];
sx q[0];
rz(-2.9009853) q[0];
rz(2.799017) q[1];
sx q[1];
rz(-0.97476417) q[1];
sx q[1];
rz(-1.906357) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9440207) q[0];
sx q[0];
rz(-0.70957843) q[0];
sx q[0];
rz(0.59528415) q[0];
rz(-2.7690976) q[2];
sx q[2];
rz(-1.8066346) q[2];
sx q[2];
rz(1.7861988) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8894549) q[1];
sx q[1];
rz(-2.7733907) q[1];
sx q[1];
rz(-2.3705269) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0577621) q[3];
sx q[3];
rz(-2.1420797) q[3];
sx q[3];
rz(0.42698241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.76413313) q[2];
sx q[2];
rz(-1.6230134) q[2];
sx q[2];
rz(1.4366478) q[2];
rz(1.7403587) q[3];
sx q[3];
rz(-1.2652206) q[3];
sx q[3];
rz(0.60825545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075994611) q[0];
sx q[0];
rz(-1.6125212) q[0];
sx q[0];
rz(-2.3377989) q[0];
rz(0.94961387) q[1];
sx q[1];
rz(-1.6751553) q[1];
sx q[1];
rz(0.11985699) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0940995) q[0];
sx q[0];
rz(-2.1764628) q[0];
sx q[0];
rz(-3.0981307) q[0];
x q[1];
rz(-0.36977936) q[2];
sx q[2];
rz(-2.7571207) q[2];
sx q[2];
rz(2.0535) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.67112982) q[1];
sx q[1];
rz(-1.5884591) q[1];
sx q[1];
rz(-1.8030241) q[1];
x q[2];
rz(1.0355789) q[3];
sx q[3];
rz(-1.03973) q[3];
sx q[3];
rz(-0.5862743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.63311657) q[2];
sx q[2];
rz(-1.8957596) q[2];
sx q[2];
rz(-3.0692549) q[2];
rz(-0.37483254) q[3];
sx q[3];
rz(-2.4852677) q[3];
sx q[3];
rz(1.9434631) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6937834) q[0];
sx q[0];
rz(-0.57149514) q[0];
sx q[0];
rz(-0.48686349) q[0];
rz(0.72987366) q[1];
sx q[1];
rz(-2.2327773) q[1];
sx q[1];
rz(-1.9015076) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86868984) q[0];
sx q[0];
rz(-1.5575952) q[0];
sx q[0];
rz(0.23916434) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8334332) q[2];
sx q[2];
rz(-1.8034168) q[2];
sx q[2];
rz(2.5831985) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.005849) q[1];
sx q[1];
rz(-1.504335) q[1];
sx q[1];
rz(0.60484109) q[1];
rz(-pi) q[2];
rz(-0.92091839) q[3];
sx q[3];
rz(-0.71483597) q[3];
sx q[3];
rz(2.9361847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1029677) q[2];
sx q[2];
rz(-0.90927783) q[2];
sx q[2];
rz(-2.4385578) q[2];
rz(-1.7317584) q[3];
sx q[3];
rz(-1.3491646) q[3];
sx q[3];
rz(2.9838802) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1773961) q[0];
sx q[0];
rz(-1.8494158) q[0];
sx q[0];
rz(2.1512206) q[0];
rz(-0.05274996) q[1];
sx q[1];
rz(-0.92664781) q[1];
sx q[1];
rz(-1.4809158) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4303362) q[0];
sx q[0];
rz(-2.7529786) q[0];
sx q[0];
rz(1.3769763) q[0];
rz(0.00077498282) q[2];
sx q[2];
rz(-1.1256071) q[2];
sx q[2];
rz(-2.0489401) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.78696886) q[1];
sx q[1];
rz(-1.0458535) q[1];
sx q[1];
rz(-2.900219) q[1];
rz(-3.092993) q[3];
sx q[3];
rz(-0.89266333) q[3];
sx q[3];
rz(2.5130659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1331553) q[2];
sx q[2];
rz(-1.6494273) q[2];
sx q[2];
rz(-2.5793502) q[2];
rz(2.0810614) q[3];
sx q[3];
rz(-0.74354592) q[3];
sx q[3];
rz(2.8806768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9437207) q[0];
sx q[0];
rz(-1.3523538) q[0];
sx q[0];
rz(-1.6725756) q[0];
rz(1.0143657) q[1];
sx q[1];
rz(-2.1140153) q[1];
sx q[1];
rz(-1.8168824) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46885195) q[0];
sx q[0];
rz(-1.6981089) q[0];
sx q[0];
rz(-1.5936113) q[0];
rz(-pi) q[1];
rz(0.84682805) q[2];
sx q[2];
rz(-2.4412051) q[2];
sx q[2];
rz(-2.2221153) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3130256) q[1];
sx q[1];
rz(-1.1638068) q[1];
sx q[1];
rz(1.160497) q[1];
x q[2];
rz(-0.91388254) q[3];
sx q[3];
rz(-1.8808639) q[3];
sx q[3];
rz(-0.1012181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5614732) q[2];
sx q[2];
rz(-1.3886398) q[2];
sx q[2];
rz(2.6980147) q[2];
rz(-2.1929072) q[3];
sx q[3];
rz(-1.1921927) q[3];
sx q[3];
rz(-0.77478066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.401944) q[0];
sx q[0];
rz(-2.1304603) q[0];
sx q[0];
rz(0.46052128) q[0];
rz(0.1000239) q[1];
sx q[1];
rz(-2.1477551) q[1];
sx q[1];
rz(1.8519648) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.078482) q[0];
sx q[0];
rz(-1.0149628) q[0];
sx q[0];
rz(-1.0312992) q[0];
rz(-pi) q[1];
rz(2.3233534) q[2];
sx q[2];
rz(-2.1365676) q[2];
sx q[2];
rz(2.408839) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.21266567) q[1];
sx q[1];
rz(-1.7576828) q[1];
sx q[1];
rz(3.0287663) q[1];
rz(-pi) q[2];
rz(-1.4054221) q[3];
sx q[3];
rz(-2.227042) q[3];
sx q[3];
rz(1.0429494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.9926247) q[2];
sx q[2];
rz(-2.830539) q[2];
sx q[2];
rz(-1.0827433) q[2];
rz(0.096171245) q[3];
sx q[3];
rz(-1.7008737) q[3];
sx q[3];
rz(1.910803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7643395) q[0];
sx q[0];
rz(-0.23582533) q[0];
sx q[0];
rz(-0.33426958) q[0];
rz(-1.224068) q[1];
sx q[1];
rz(-1.5806438) q[1];
sx q[1];
rz(-0.28265488) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3960421) q[0];
sx q[0];
rz(-0.96691416) q[0];
sx q[0];
rz(-2.1349483) q[0];
rz(-pi) q[1];
rz(-2.5731509) q[2];
sx q[2];
rz(-0.38707765) q[2];
sx q[2];
rz(2.8708411) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5987451) q[1];
sx q[1];
rz(-1.7464906) q[1];
sx q[1];
rz(0.41136841) q[1];
rz(-1.4275527) q[3];
sx q[3];
rz(-2.4510265) q[3];
sx q[3];
rz(-2.7018202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.21215542) q[2];
sx q[2];
rz(-1.9781457) q[2];
sx q[2];
rz(-0.7412509) q[2];
rz(-0.50179982) q[3];
sx q[3];
rz(-1.3391677) q[3];
sx q[3];
rz(0.58399502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.042645) q[0];
sx q[0];
rz(-0.89530033) q[0];
sx q[0];
rz(-2.6328971) q[0];
rz(0.11518654) q[1];
sx q[1];
rz(-2.4337264) q[1];
sx q[1];
rz(-0.68181109) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17446974) q[0];
sx q[0];
rz(-0.57708626) q[0];
sx q[0];
rz(-1.2601) q[0];
x q[1];
rz(-0.85068662) q[2];
sx q[2];
rz(-2.2518573) q[2];
sx q[2];
rz(1.016664) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.0971113) q[1];
sx q[1];
rz(-2.0437355) q[1];
sx q[1];
rz(0.60826917) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.82716771) q[3];
sx q[3];
rz(-0.70561545) q[3];
sx q[3];
rz(-3.0527715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.29601413) q[2];
sx q[2];
rz(-1.6342376) q[2];
sx q[2];
rz(-2.8005023) q[2];
rz(-2.0579445) q[3];
sx q[3];
rz(-2.3458979) q[3];
sx q[3];
rz(0.17102374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9733799) q[0];
sx q[0];
rz(-1.4773049) q[0];
sx q[0];
rz(2.1422577) q[0];
rz(2.534261) q[1];
sx q[1];
rz(-2.2317531) q[1];
sx q[1];
rz(-2.9324525) q[1];
rz(0.66424673) q[2];
sx q[2];
rz(-1.1462117) q[2];
sx q[2];
rz(0.29427634) q[2];
rz(-2.7568983) q[3];
sx q[3];
rz(-1.3833429) q[3];
sx q[3];
rz(-0.27237567) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
