OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8024017) q[0];
sx q[0];
rz(-1.985476) q[0];
sx q[0];
rz(-0.51751408) q[0];
rz(1.327688) q[1];
sx q[1];
rz(-2.2496536) q[1];
sx q[1];
rz(0.54203027) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1928312) q[0];
sx q[0];
rz(-1.3997215) q[0];
sx q[0];
rz(-0.77617742) q[0];
rz(-pi) q[1];
rz(1.8250685) q[2];
sx q[2];
rz(-2.3208127) q[2];
sx q[2];
rz(0.17707846) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4229976) q[1];
sx q[1];
rz(-1.7684775) q[1];
sx q[1];
rz(2.343178) q[1];
rz(2.5045218) q[3];
sx q[3];
rz(-0.28709295) q[3];
sx q[3];
rz(1.4305103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1338256) q[2];
sx q[2];
rz(-2.0962891) q[2];
sx q[2];
rz(2.2110151) q[2];
rz(-0.13115701) q[3];
sx q[3];
rz(-2.7635062) q[3];
sx q[3];
rz(1.650943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5432878) q[0];
sx q[0];
rz(-2.2853993) q[0];
sx q[0];
rz(-1.244586) q[0];
rz(-0.92567956) q[1];
sx q[1];
rz(-1.2558179) q[1];
sx q[1];
rz(-1.6474887) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7667395) q[0];
sx q[0];
rz(-1.8962066) q[0];
sx q[0];
rz(-2.297657) q[0];
rz(-pi) q[1];
x q[1];
rz(1.235512) q[2];
sx q[2];
rz(-2.2163894) q[2];
sx q[2];
rz(-0.52926999) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0665613) q[1];
sx q[1];
rz(-1.9157456) q[1];
sx q[1];
rz(-0.48081545) q[1];
rz(-pi) q[2];
rz(-0.30416137) q[3];
sx q[3];
rz(-1.5790325) q[3];
sx q[3];
rz(1.3502163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99574789) q[2];
sx q[2];
rz(-2.6971942) q[2];
sx q[2];
rz(-0.90163976) q[2];
rz(1.3580648) q[3];
sx q[3];
rz(-1.3172904) q[3];
sx q[3];
rz(-2.5942514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77483344) q[0];
sx q[0];
rz(-1.1792553) q[0];
sx q[0];
rz(-1.1460079) q[0];
rz(-0.92578069) q[1];
sx q[1];
rz(-1.0113167) q[1];
sx q[1];
rz(-2.944223) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9512874) q[0];
sx q[0];
rz(-1.4337375) q[0];
sx q[0];
rz(-2.1648315) q[0];
rz(-pi) q[1];
rz(-0.80855753) q[2];
sx q[2];
rz(-0.41038358) q[2];
sx q[2];
rz(1.9661599) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.92835411) q[1];
sx q[1];
rz(-0.23121195) q[1];
sx q[1];
rz(1.3998881) q[1];
rz(-2.8160964) q[3];
sx q[3];
rz(-0.56762689) q[3];
sx q[3];
rz(-2.5095255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7388514) q[2];
sx q[2];
rz(-0.98029843) q[2];
sx q[2];
rz(2.951238) q[2];
rz(3.0800152) q[3];
sx q[3];
rz(-2.172251) q[3];
sx q[3];
rz(-2.0891345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9848118) q[0];
sx q[0];
rz(-1.4182014) q[0];
sx q[0];
rz(2.5132827) q[0];
rz(1.5208987) q[1];
sx q[1];
rz(-1.9338806) q[1];
sx q[1];
rz(2.3627088) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91858868) q[0];
sx q[0];
rz(-0.99705925) q[0];
sx q[0];
rz(2.6911435) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2912108) q[2];
sx q[2];
rz(-1.9182854) q[2];
sx q[2];
rz(-2.8742032) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.35793909) q[1];
sx q[1];
rz(-0.049555819) q[1];
sx q[1];
rz(2.8663584) q[1];
rz(0.67653894) q[3];
sx q[3];
rz(-2.153844) q[3];
sx q[3];
rz(-0.18932115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4817619) q[2];
sx q[2];
rz(-2.063282) q[2];
sx q[2];
rz(-4.3241186e-05) q[2];
rz(2.0569885) q[3];
sx q[3];
rz(-1.9856039) q[3];
sx q[3];
rz(-2.675975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4701009) q[0];
sx q[0];
rz(-1.6974314) q[0];
sx q[0];
rz(-1.2985562) q[0];
rz(0.57688722) q[1];
sx q[1];
rz(-1.8678317) q[1];
sx q[1];
rz(-2.6947122) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80432207) q[0];
sx q[0];
rz(-1.7511586) q[0];
sx q[0];
rz(0.65482636) q[0];
x q[1];
rz(-3.066874) q[2];
sx q[2];
rz(-1.8693035) q[2];
sx q[2];
rz(-3.0252688) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3342377) q[1];
sx q[1];
rz(-0.45354834) q[1];
sx q[1];
rz(2.0475044) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6303667) q[3];
sx q[3];
rz(-1.7559663) q[3];
sx q[3];
rz(0.51051729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.4345066) q[2];
sx q[2];
rz(-2.862317) q[2];
sx q[2];
rz(2.9217829) q[2];
rz(1.48742) q[3];
sx q[3];
rz(-1.6916964) q[3];
sx q[3];
rz(-0.69948227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.952482) q[0];
sx q[0];
rz(-2.7087152) q[0];
sx q[0];
rz(-0.89749807) q[0];
rz(-1.7706361) q[1];
sx q[1];
rz(-2.1015344) q[1];
sx q[1];
rz(-2.2460361) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5286551) q[0];
sx q[0];
rz(-0.85451689) q[0];
sx q[0];
rz(-1.3837293) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3562548) q[2];
sx q[2];
rz(-2.3599652) q[2];
sx q[2];
rz(1.2981894) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8297146) q[1];
sx q[1];
rz(-2.6304881) q[1];
sx q[1];
rz(0.73181499) q[1];
rz(-pi) q[2];
x q[2];
rz(0.47116739) q[3];
sx q[3];
rz(-1.1352603) q[3];
sx q[3];
rz(-2.7903231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.36876496) q[2];
sx q[2];
rz(-1.1026829) q[2];
sx q[2];
rz(2.9844798) q[2];
rz(-0.67240063) q[3];
sx q[3];
rz(-1.7417358) q[3];
sx q[3];
rz(3.0290643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2812578) q[0];
sx q[0];
rz(-2.4761138) q[0];
sx q[0];
rz(-1.459664) q[0];
rz(2.7809987) q[1];
sx q[1];
rz(-1.6723885) q[1];
sx q[1];
rz(1.8416539) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0361745) q[0];
sx q[0];
rz(-1.6929394) q[0];
sx q[0];
rz(-0.41897498) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62076135) q[2];
sx q[2];
rz(-1.1859535) q[2];
sx q[2];
rz(-2.8892725) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2092525) q[1];
sx q[1];
rz(-2.8289218) q[1];
sx q[1];
rz(-0.31950848) q[1];
rz(-pi) q[2];
rz(1.9770558) q[3];
sx q[3];
rz(-0.99525577) q[3];
sx q[3];
rz(1.5062603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.931687) q[2];
sx q[2];
rz(-2.1114712) q[2];
sx q[2];
rz(-0.21200655) q[2];
rz(-2.0906406) q[3];
sx q[3];
rz(-1.7651599) q[3];
sx q[3];
rz(3.0534548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5365005) q[0];
sx q[0];
rz(-2.3516646) q[0];
sx q[0];
rz(2.9686046) q[0];
rz(-2.0269003) q[1];
sx q[1];
rz(-2.098691) q[1];
sx q[1];
rz(0.10985049) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86533029) q[0];
sx q[0];
rz(-1.0181576) q[0];
sx q[0];
rz(-1.7981973) q[0];
rz(-pi) q[1];
rz(0.26066653) q[2];
sx q[2];
rz(-2.3015966) q[2];
sx q[2];
rz(2.276401) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.50625769) q[1];
sx q[1];
rz(-1.3999363) q[1];
sx q[1];
rz(1.3116763) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55446378) q[3];
sx q[3];
rz(-1.672847) q[3];
sx q[3];
rz(-2.3481675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7383808) q[2];
sx q[2];
rz(-1.9387559) q[2];
sx q[2];
rz(-2.476725) q[2];
rz(-0.46857771) q[3];
sx q[3];
rz(-1.6178308) q[3];
sx q[3];
rz(2.328228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68798962) q[0];
sx q[0];
rz(-2.5316694) q[0];
sx q[0];
rz(0.74158057) q[0];
rz(1.954151) q[1];
sx q[1];
rz(-1.5267742) q[1];
sx q[1];
rz(-2.5724519) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4837836) q[0];
sx q[0];
rz(-1.2783733) q[0];
sx q[0];
rz(-2.800709) q[0];
x q[1];
rz(2.9958565) q[2];
sx q[2];
rz(-2.6837641) q[2];
sx q[2];
rz(0.65953461) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.38005623) q[1];
sx q[1];
rz(-1.1498701) q[1];
sx q[1];
rz(-0.50889877) q[1];
x q[2];
rz(0.82707246) q[3];
sx q[3];
rz(-2.2571392) q[3];
sx q[3];
rz(-2.4857869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.49447122) q[2];
sx q[2];
rz(-2.1529866) q[2];
sx q[2];
rz(-0.76623255) q[2];
rz(1.9689485) q[3];
sx q[3];
rz(-1.6021043) q[3];
sx q[3];
rz(-0.22181454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93429339) q[0];
sx q[0];
rz(-2.793029) q[0];
sx q[0];
rz(-0.19185129) q[0];
rz(-0.81709298) q[1];
sx q[1];
rz(-2.8849738) q[1];
sx q[1];
rz(-2.4339035) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.548508) q[0];
sx q[0];
rz(-1.3258717) q[0];
sx q[0];
rz(1.3913416) q[0];
x q[1];
rz(0.97918503) q[2];
sx q[2];
rz(-1.8509097) q[2];
sx q[2];
rz(-1.6620987) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.62218828) q[1];
sx q[1];
rz(-0.71336245) q[1];
sx q[1];
rz(0.98825561) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58587425) q[3];
sx q[3];
rz(-1.0535686) q[3];
sx q[3];
rz(-2.6894803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2231458) q[2];
sx q[2];
rz(-2.196849) q[2];
sx q[2];
rz(2.4165912) q[2];
rz(-2.9265192) q[3];
sx q[3];
rz(-2.8408065) q[3];
sx q[3];
rz(-2.422629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.219915) q[0];
sx q[0];
rz(-2.324993) q[0];
sx q[0];
rz(-2.8391229) q[0];
rz(-2.0401781) q[1];
sx q[1];
rz(-1.9093724) q[1];
sx q[1];
rz(-2.2473635) q[1];
rz(1.3742555) q[2];
sx q[2];
rz(-2.5977618) q[2];
sx q[2];
rz(0.14731461) q[2];
rz(-0.054417944) q[3];
sx q[3];
rz(-2.6734753) q[3];
sx q[3];
rz(-2.5772167) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
