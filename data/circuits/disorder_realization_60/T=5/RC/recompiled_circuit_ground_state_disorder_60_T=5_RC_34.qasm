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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5954689) q[0];
sx q[0];
rz(-1.4065388) q[0];
sx q[0];
rz(1.6857273) q[0];
x q[1];
rz(-1.6923793) q[2];
sx q[2];
rz(-1.087153) q[2];
sx q[2];
rz(-1.3822966) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2985483) q[1];
sx q[1];
rz(-1.4433651) q[1];
sx q[1];
rz(0.28020466) q[1];
rz(0.61986519) q[3];
sx q[3];
rz(-2.4976118) q[3];
sx q[3];
rz(2.5247273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7264318) q[2];
sx q[2];
rz(-1.2434604) q[2];
sx q[2];
rz(2.4746056) q[2];
rz(-2.4994182) q[3];
sx q[3];
rz(-2.8447076) q[3];
sx q[3];
rz(1.4750397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.058454839) q[1];
sx q[1];
rz(-2.7775601) q[1];
sx q[1];
rz(2.2296925) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6175323) q[0];
sx q[0];
rz(-0.41818163) q[0];
sx q[0];
rz(1.4773744) q[0];
rz(2.2445685) q[2];
sx q[2];
rz(-1.6260212) q[2];
sx q[2];
rz(1.3185475) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81296125) q[1];
sx q[1];
rz(-1.3282623) q[1];
sx q[1];
rz(-0.46142908) q[1];
rz(-0.1095406) q[3];
sx q[3];
rz(-0.44407223) q[3];
sx q[3];
rz(-0.20694297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6902206) q[2];
sx q[2];
rz(-2.8896832) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5628691) q[0];
sx q[0];
rz(-0.27414411) q[0];
sx q[0];
rz(-1.3124056) q[0];
rz(-2.9263272) q[1];
sx q[1];
rz(-1.637641) q[1];
sx q[1];
rz(2.2832835) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90412437) q[0];
sx q[0];
rz(-1.6209813) q[0];
sx q[0];
rz(1.8367358) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9163754) q[2];
sx q[2];
rz(-1.2379939) q[2];
sx q[2];
rz(2.107258) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5897164) q[1];
sx q[1];
rz(-1.7360506) q[1];
sx q[1];
rz(2.698353) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8668601) q[3];
sx q[3];
rz(-1.8157457) q[3];
sx q[3];
rz(1.7942269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.12576088) q[2];
sx q[2];
rz(-2.5845926) q[2];
sx q[2];
rz(-2.1324943) q[2];
rz(2.5741746) q[3];
sx q[3];
rz(-1.6332473) q[3];
sx q[3];
rz(1.6217888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.0030982) q[0];
sx q[0];
rz(-0.69366413) q[0];
sx q[0];
rz(2.9487115) q[0];
rz(-2.0774972) q[1];
sx q[1];
rz(-0.27609438) q[1];
sx q[1];
rz(0.9562794) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.103699) q[0];
sx q[0];
rz(-0.65019817) q[0];
sx q[0];
rz(2.8481087) q[0];
rz(2.2334571) q[2];
sx q[2];
rz(-1.081433) q[2];
sx q[2];
rz(0.2950677) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.87188211) q[1];
sx q[1];
rz(-2.080889) q[1];
sx q[1];
rz(-1.2126232) q[1];
rz(-pi) q[2];
rz(-0.30711912) q[3];
sx q[3];
rz(-2.2132717) q[3];
sx q[3];
rz(2.0406587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8320147) q[2];
sx q[2];
rz(-0.39065209) q[2];
sx q[2];
rz(0.94069707) q[2];
rz(2.7637123) q[3];
sx q[3];
rz(-2.2418719) q[3];
sx q[3];
rz(-0.64062029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2864895) q[0];
sx q[0];
rz(-1.531895) q[0];
sx q[0];
rz(1.3252277) q[0];
rz(2.8638966) q[1];
sx q[1];
rz(-2.1546021) q[1];
sx q[1];
rz(-1.9353386) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4770484) q[0];
sx q[0];
rz(-1.3788851) q[0];
sx q[0];
rz(0.080587793) q[0];
rz(-pi) q[1];
rz(2.5304753) q[2];
sx q[2];
rz(-1.9131994) q[2];
sx q[2];
rz(2.8813206) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6305651) q[1];
sx q[1];
rz(-1.2113809) q[1];
sx q[1];
rz(0.014822249) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5370977) q[3];
sx q[3];
rz(-0.2911199) q[3];
sx q[3];
rz(-2.0184984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.18579927) q[2];
sx q[2];
rz(-0.15115559) q[2];
sx q[2];
rz(1.7929057) q[2];
rz(-1.0079314) q[3];
sx q[3];
rz(-1.4069087) q[3];
sx q[3];
rz(-3.1048043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0607097) q[0];
sx q[0];
rz(-3.1302629) q[0];
sx q[0];
rz(-0.19705535) q[0];
rz(-1.9673678) q[1];
sx q[1];
rz(-0.86001992) q[1];
sx q[1];
rz(2.7161782) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4482314) q[0];
sx q[0];
rz(-1.5299242) q[0];
sx q[0];
rz(0.46092439) q[0];
x q[1];
rz(1.4324911) q[2];
sx q[2];
rz(-2.8347361) q[2];
sx q[2];
rz(2.8371161) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0531333) q[1];
sx q[1];
rz(-1.3091374) q[1];
sx q[1];
rz(-2.7105374) q[1];
x q[2];
rz(2.3415575) q[3];
sx q[3];
rz(-2.2657395) q[3];
sx q[3];
rz(1.8025229) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.60874385) q[2];
sx q[2];
rz(-3.0259735) q[2];
sx q[2];
rz(-2.4382584) q[2];
rz(-2.9902847) q[3];
sx q[3];
rz(-1.9688828) q[3];
sx q[3];
rz(2.5942588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5587826) q[0];
sx q[0];
rz(-2.3334239) q[0];
sx q[0];
rz(-0.95635828) q[0];
rz(-0.76902485) q[1];
sx q[1];
rz(-1.290657) q[1];
sx q[1];
rz(-1.0395315) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7881986) q[0];
sx q[0];
rz(-2.0457342) q[0];
sx q[0];
rz(-0.026508412) q[0];
rz(-pi) q[1];
rz(2.2724067) q[2];
sx q[2];
rz(-2.1034046) q[2];
sx q[2];
rz(-2.5857596) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.49756296) q[1];
sx q[1];
rz(-2.7179238) q[1];
sx q[1];
rz(-0.81987857) q[1];
x q[2];
rz(-0.2110699) q[3];
sx q[3];
rz(-1.8976334) q[3];
sx q[3];
rz(-1.1093931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9397554) q[2];
sx q[2];
rz(-0.73362437) q[2];
sx q[2];
rz(1.5038097) q[2];
rz(-0.87320295) q[3];
sx q[3];
rz(-0.88770294) q[3];
sx q[3];
rz(0.33151117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37209508) q[0];
sx q[0];
rz(-0.3898507) q[0];
sx q[0];
rz(1.9918342) q[0];
rz(0.86504966) q[1];
sx q[1];
rz(-1.9416092) q[1];
sx q[1];
rz(-1.7080151) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5652611) q[0];
sx q[0];
rz(-2.1339685) q[0];
sx q[0];
rz(-2.7050326) q[0];
rz(0.47793598) q[2];
sx q[2];
rz(-1.897287) q[2];
sx q[2];
rz(0.70604347) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5123957) q[1];
sx q[1];
rz(-1.6586972) q[1];
sx q[1];
rz(2.0950661) q[1];
x q[2];
rz(0.20345466) q[3];
sx q[3];
rz(-1.1622114) q[3];
sx q[3];
rz(-0.83554283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4992497) q[2];
sx q[2];
rz(-1.6982102) q[2];
sx q[2];
rz(0.7439822) q[2];
rz(2.3104525) q[3];
sx q[3];
rz(-1.4011413) q[3];
sx q[3];
rz(2.7861815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
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
rz(-1.5903541) q[1];
sx q[1];
rz(-0.23983461) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3837357) q[0];
sx q[0];
rz(-2.3546887) q[0];
sx q[0];
rz(2.402124) q[0];
rz(1.5101899) q[2];
sx q[2];
rz(-2.1327634) q[2];
sx q[2];
rz(0.89982239) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0785324) q[1];
sx q[1];
rz(-0.70101798) q[1];
sx q[1];
rz(1.307742) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0795709) q[3];
sx q[3];
rz(-2.6497132) q[3];
sx q[3];
rz(0.57152257) q[3];
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
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1455602) q[0];
sx q[0];
rz(-1.6479489) q[0];
sx q[0];
rz(-2.5160568) q[0];
rz(1.7690376) q[1];
sx q[1];
rz(-2.1452429) q[1];
sx q[1];
rz(-0.58403429) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16957078) q[0];
sx q[0];
rz(-1.941328) q[0];
sx q[0];
rz(-1.8448004) q[0];
rz(-pi) q[1];
rz(-2.9461924) q[2];
sx q[2];
rz(-0.6740734) q[2];
sx q[2];
rz(-0.96701996) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8019097) q[1];
sx q[1];
rz(-2.3900635) q[1];
sx q[1];
rz(-0.32210334) q[1];
rz(-1.8409506) q[3];
sx q[3];
rz(-0.90357257) q[3];
sx q[3];
rz(-3.0285578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.11223101) q[2];
sx q[2];
rz(-1.4697305) q[2];
sx q[2];
rz(-1.3545797) q[2];
rz(2.3692047) q[3];
sx q[3];
rz(-1.200501) q[3];
sx q[3];
rz(-1.6105917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
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
rz(1.7583379) q[2];
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
