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
rz(-1.5457431) q[1];
sx q[1];
rz(0.63313142) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0057982) q[0];
sx q[0];
rz(-1.457419) q[0];
sx q[0];
rz(2.9762639) q[0];
x q[1];
rz(-1.6923793) q[2];
sx q[2];
rz(-2.0544397) q[2];
sx q[2];
rz(1.3822966) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.76430815) q[1];
sx q[1];
rz(-1.2929243) q[1];
sx q[1];
rz(-1.4382526) q[1];
rz(-pi) q[2];
rz(-2.5217275) q[3];
sx q[3];
rz(-0.64398089) q[3];
sx q[3];
rz(0.61686531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7264318) q[2];
sx q[2];
rz(-1.8981322) q[2];
sx q[2];
rz(0.66698709) q[2];
rz(-0.64217448) q[3];
sx q[3];
rz(-2.8447076) q[3];
sx q[3];
rz(-1.4750397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.967763) q[0];
sx q[0];
rz(-2.1673295) q[0];
sx q[0];
rz(0.94671384) q[0];
rz(0.058454839) q[1];
sx q[1];
rz(-2.7775601) q[1];
sx q[1];
rz(0.91190016) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4218876) q[0];
sx q[0];
rz(-1.9870409) q[0];
sx q[0];
rz(-3.1001607) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0709708) q[2];
sx q[2];
rz(-0.89824072) q[2];
sx q[2];
rz(0.20820752) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.87668858) q[1];
sx q[1];
rz(-1.1238681) q[1];
sx q[1];
rz(1.3012215) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44174444) q[3];
sx q[3];
rz(-1.6177804) q[3];
sx q[3];
rz(1.6787501) q[3];
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
rz(2.1515324) q[3];
sx q[3];
rz(-1.8127352) q[3];
sx q[3];
rz(1.2216074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-2.2832835) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90412437) q[0];
sx q[0];
rz(-1.5206114) q[0];
sx q[0];
rz(1.8367358) q[0];
rz(-2.1444001) q[2];
sx q[2];
rz(-2.7421087) q[2];
sx q[2];
rz(1.6459207) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2006113) q[1];
sx q[1];
rz(-2.0075782) q[1];
sx q[1];
rz(-1.753356) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86238548) q[3];
sx q[3];
rz(-0.38194712) q[3];
sx q[3];
rz(0.89513429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0158318) q[2];
sx q[2];
rz(-2.5845926) q[2];
sx q[2];
rz(2.1324943) q[2];
rz(2.5741746) q[3];
sx q[3];
rz(-1.6332473) q[3];
sx q[3];
rz(-1.5198038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(2.1853133) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67502695) q[0];
sx q[0];
rz(-0.95272946) q[0];
sx q[0];
rz(-1.7873554) q[0];
rz(-pi) q[1];
rz(-2.2334571) q[2];
sx q[2];
rz(-2.0601597) q[2];
sx q[2];
rz(-2.846525) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6156441) q[1];
sx q[1];
rz(-2.5275349) q[1];
sx q[1];
rz(-2.581852) q[1];
rz(-pi) q[2];
rz(2.8344735) q[3];
sx q[3];
rz(-2.2132717) q[3];
sx q[3];
rz(-1.100934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85510319) q[0];
sx q[0];
rz(-1.531895) q[0];
sx q[0];
rz(-1.3252277) q[0];
rz(2.8638966) q[1];
sx q[1];
rz(-0.9869906) q[1];
sx q[1];
rz(1.9353386) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6645443) q[0];
sx q[0];
rz(-1.3788851) q[0];
sx q[0];
rz(-3.0610049) q[0];
x q[1];
rz(1.9812859) q[2];
sx q[2];
rz(-2.1417981) q[2];
sx q[2];
rz(1.0795019) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5110276) q[1];
sx q[1];
rz(-1.9302118) q[1];
sx q[1];
rz(3.1267704) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7394735) q[3];
sx q[3];
rz(-1.8092098) q[3];
sx q[3];
rz(-1.7478398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.18579927) q[2];
sx q[2];
rz(-2.9904371) q[2];
sx q[2];
rz(1.7929057) q[2];
rz(1.0079314) q[3];
sx q[3];
rz(-1.734684) q[3];
sx q[3];
rz(0.036788363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0607097) q[0];
sx q[0];
rz(-0.01132975) q[0];
sx q[0];
rz(-2.9445373) q[0];
rz(-1.9673678) q[1];
sx q[1];
rz(-2.2815727) q[1];
sx q[1];
rz(0.42541447) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6933612) q[0];
sx q[0];
rz(-1.5299242) q[0];
sx q[0];
rz(-2.6806683) q[0];
x q[1];
rz(1.4324911) q[2];
sx q[2];
rz(-2.8347361) q[2];
sx q[2];
rz(-0.30447659) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.08845932) q[1];
sx q[1];
rz(-1.3091374) q[1];
sx q[1];
rz(-0.43105521) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80003512) q[3];
sx q[3];
rz(-2.2657395) q[3];
sx q[3];
rz(-1.8025229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5328488) q[2];
sx q[2];
rz(-3.0259735) q[2];
sx q[2];
rz(-0.70333424) q[2];
rz(-0.15130791) q[3];
sx q[3];
rz(-1.1727099) q[3];
sx q[3];
rz(-0.54733384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5587826) q[0];
sx q[0];
rz(-0.8081688) q[0];
sx q[0];
rz(-2.1852344) q[0];
rz(-2.3725678) q[1];
sx q[1];
rz(-1.8509357) q[1];
sx q[1];
rz(-1.0395315) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7881986) q[0];
sx q[0];
rz(-2.0457342) q[0];
sx q[0];
rz(0.026508412) q[0];
rz(-0.86918594) q[2];
sx q[2];
rz(-1.0381881) q[2];
sx q[2];
rz(-0.55583304) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.49756296) q[1];
sx q[1];
rz(-2.7179238) q[1];
sx q[1];
rz(-2.3217141) q[1];
x q[2];
rz(1.0172306) q[3];
sx q[3];
rz(-0.38700208) q[3];
sx q[3];
rz(1.4437046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2018373) q[2];
sx q[2];
rz(-0.73362437) q[2];
sx q[2];
rz(1.6377829) q[2];
rz(2.2683897) q[3];
sx q[3];
rz(-0.88770294) q[3];
sx q[3];
rz(0.33151117) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7694976) q[0];
sx q[0];
rz(-2.7517419) q[0];
sx q[0];
rz(1.9918342) q[0];
rz(-2.276543) q[1];
sx q[1];
rz(-1.1999835) q[1];
sx q[1];
rz(1.7080151) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5763315) q[0];
sx q[0];
rz(-1.0076242) q[0];
sx q[0];
rz(0.43656008) q[0];
rz(-pi) q[1];
rz(-2.5069889) q[2];
sx q[2];
rz(-0.57159492) q[2];
sx q[2];
rz(-1.7224479) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5123957) q[1];
sx q[1];
rz(-1.6586972) q[1];
sx q[1];
rz(2.0950661) q[1];
x q[2];
rz(1.9870342) q[3];
sx q[3];
rz(-1.3842954) q[3];
sx q[3];
rz(-0.65346741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.64234298) q[2];
sx q[2];
rz(-1.6982102) q[2];
sx q[2];
rz(0.7439822) q[2];
rz(-2.3104525) q[3];
sx q[3];
rz(-1.4011413) q[3];
sx q[3];
rz(-2.7861815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0396073) q[0];
sx q[0];
rz(-0.87431854) q[0];
sx q[0];
rz(-0.71514493) q[0];
rz(-1.2932418) q[1];
sx q[1];
rz(-1.5512385) q[1];
sx q[1];
rz(-2.901758) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3837357) q[0];
sx q[0];
rz(-0.78690398) q[0];
sx q[0];
rz(-0.73946865) q[0];
x q[1];
rz(-2.5787966) q[2];
sx q[2];
rz(-1.6220731) q[2];
sx q[2];
rz(2.4382961) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.71069392) q[1];
sx q[1];
rz(-1.4022809) q[1];
sx q[1];
rz(-2.2546143) q[1];
rz(-pi) q[2];
rz(3.0795709) q[3];
sx q[3];
rz(-2.6497132) q[3];
sx q[3];
rz(-2.5700701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0665325) q[2];
sx q[2];
rz(-1.5542474) q[2];
sx q[2];
rz(1.4934322) q[2];
rz(3.1028124) q[3];
sx q[3];
rz(-1.4977692) q[3];
sx q[3];
rz(1.405976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1455602) q[0];
sx q[0];
rz(-1.6479489) q[0];
sx q[0];
rz(2.5160568) q[0];
rz(-1.372555) q[1];
sx q[1];
rz(-0.99634975) q[1];
sx q[1];
rz(-2.5575584) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8417977) q[0];
sx q[0];
rz(-1.3158176) q[0];
sx q[0];
rz(-0.38354446) q[0];
x q[1];
rz(-2.9461924) q[2];
sx q[2];
rz(-2.4675193) q[2];
sx q[2];
rz(0.96701996) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6713114) q[1];
sx q[1];
rz(-1.7886506) q[1];
sx q[1];
rz(0.72524986) q[1];
x q[2];
rz(1.8409506) q[3];
sx q[3];
rz(-2.2380201) q[3];
sx q[3];
rz(-3.0285578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11223101) q[2];
sx q[2];
rz(-1.6718622) q[2];
sx q[2];
rz(-1.7870129) q[2];
rz(2.3692047) q[3];
sx q[3];
rz(-1.200501) q[3];
sx q[3];
rz(-1.6105917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1188485) q[0];
sx q[0];
rz(-0.72853715) q[0];
sx q[0];
rz(-1.5991612) q[0];
rz(2.6842885) q[1];
sx q[1];
rz(-2.2685469) q[1];
sx q[1];
rz(-0.1874478) q[1];
rz(-1.7583379) q[2];
sx q[2];
rz(-1.493307) q[2];
sx q[2];
rz(-1.1004352) q[2];
rz(-1.1111835) q[3];
sx q[3];
rz(-1.8596953) q[3];
sx q[3];
rz(0.019712899) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
