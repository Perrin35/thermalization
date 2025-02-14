OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.87822479) q[0];
sx q[0];
rz(3.8642519) q[0];
sx q[0];
rz(10.470358) q[0];
rz(-3.1388406) q[1];
sx q[1];
rz(-1.6834919) q[1];
sx q[1];
rz(-2.4207065) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3430735) q[0];
sx q[0];
rz(-1.6630173) q[0];
sx q[0];
rz(-2.2416122) q[0];
x q[1];
rz(-0.051250967) q[2];
sx q[2];
rz(-1.7619216) q[2];
sx q[2];
rz(2.4911936) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8869141) q[1];
sx q[1];
rz(-2.6638835) q[1];
sx q[1];
rz(0.63666792) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3985596) q[3];
sx q[3];
rz(-1.6662906) q[3];
sx q[3];
rz(-2.9955289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.8568933) q[2];
sx q[2];
rz(-1.2130986) q[2];
sx q[2];
rz(2.6846679) q[2];
rz(2.4871155) q[3];
sx q[3];
rz(-0.5623397) q[3];
sx q[3];
rz(0.42403179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65861312) q[0];
sx q[0];
rz(-2.787866) q[0];
sx q[0];
rz(0.60802996) q[0];
rz(-1.1990625) q[1];
sx q[1];
rz(-0.46747318) q[1];
sx q[1];
rz(2.3466568) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89657822) q[0];
sx q[0];
rz(-2.6908186) q[0];
sx q[0];
rz(-0.71433492) q[0];
rz(-3.0367467) q[2];
sx q[2];
rz(-2.8730132) q[2];
sx q[2];
rz(-1.2152482) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.39917377) q[1];
sx q[1];
rz(-1.8493506) q[1];
sx q[1];
rz(-2.2225077) q[1];
x q[2];
rz(-0.43585082) q[3];
sx q[3];
rz(-0.66879771) q[3];
sx q[3];
rz(-0.28025249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.46637154) q[2];
sx q[2];
rz(-1.0543062) q[2];
sx q[2];
rz(1.6611453) q[2];
rz(-2.772707) q[3];
sx q[3];
rz(-1.2494913) q[3];
sx q[3];
rz(0.7307542) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9166301) q[0];
sx q[0];
rz(-2.1376305) q[0];
sx q[0];
rz(-2.4555901) q[0];
rz(-1.2432159) q[1];
sx q[1];
rz(-0.22626433) q[1];
sx q[1];
rz(-2.539198) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7334977) q[0];
sx q[0];
rz(-1.0801093) q[0];
sx q[0];
rz(-0.48510929) q[0];
rz(-3.1092293) q[2];
sx q[2];
rz(-1.5302858) q[2];
sx q[2];
rz(-1.4537741) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5205875) q[1];
sx q[1];
rz(-2.5650254) q[1];
sx q[1];
rz(-1.7496878) q[1];
rz(-pi) q[2];
rz(-0.32633324) q[3];
sx q[3];
rz(-2.0179393) q[3];
sx q[3];
rz(-2.9016419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7556222) q[2];
sx q[2];
rz(-1.3034416) q[2];
sx q[2];
rz(1.2516359) q[2];
rz(0.55164591) q[3];
sx q[3];
rz(-2.9988204) q[3];
sx q[3];
rz(0.84757203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30152339) q[0];
sx q[0];
rz(-1.9020377) q[0];
sx q[0];
rz(3.0384592) q[0];
rz(2.3221807) q[1];
sx q[1];
rz(-0.90140072) q[1];
sx q[1];
rz(-2.588396) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55692907) q[0];
sx q[0];
rz(-1.2069261) q[0];
sx q[0];
rz(0.30839021) q[0];
rz(-pi) q[1];
rz(2.1125421) q[2];
sx q[2];
rz(-0.61675393) q[2];
sx q[2];
rz(2.665131) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9576615) q[1];
sx q[1];
rz(-1.4344771) q[1];
sx q[1];
rz(-1.0971954) q[1];
x q[2];
rz(-3.0402667) q[3];
sx q[3];
rz(-2.000371) q[3];
sx q[3];
rz(2.5954036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4210522) q[2];
sx q[2];
rz(-0.74137509) q[2];
sx q[2];
rz(-0.53483024) q[2];
rz(-2.363291) q[3];
sx q[3];
rz(-2.4141267) q[3];
sx q[3];
rz(2.9261869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4541723) q[0];
sx q[0];
rz(-1.146831) q[0];
sx q[0];
rz(2.4046894) q[0];
rz(-3.0829644) q[1];
sx q[1];
rz(-2.3643989) q[1];
sx q[1];
rz(2.5545205) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.042087) q[0];
sx q[0];
rz(-1.6027959) q[0];
sx q[0];
rz(-0.93670292) q[0];
rz(-pi) q[1];
rz(-1.8002073) q[2];
sx q[2];
rz(-1.6322114) q[2];
sx q[2];
rz(-3.1385076) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.88564155) q[1];
sx q[1];
rz(-1.5382861) q[1];
sx q[1];
rz(1.2684275) q[1];
rz(-pi) q[2];
rz(-0.30549099) q[3];
sx q[3];
rz(-1.2661714) q[3];
sx q[3];
rz(-2.6412499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.29207841) q[2];
sx q[2];
rz(-0.92486113) q[2];
sx q[2];
rz(2.4936131) q[2];
rz(-2.5255711) q[3];
sx q[3];
rz(-1.6391305) q[3];
sx q[3];
rz(3.093921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-2.7306526) q[0];
sx q[0];
rz(-2.3766282) q[0];
sx q[0];
rz(-2.1867645) q[0];
rz(-3.094589) q[1];
sx q[1];
rz(-2.4689597) q[1];
sx q[1];
rz(2.1254553) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92378855) q[0];
sx q[0];
rz(-1.0402443) q[0];
sx q[0];
rz(1.0180044) q[0];
x q[1];
rz(-0.72332763) q[2];
sx q[2];
rz(-1.1473341) q[2];
sx q[2];
rz(-0.14191574) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2909866) q[1];
sx q[1];
rz(-2.6358446) q[1];
sx q[1];
rz(-3.1116897) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0623779) q[3];
sx q[3];
rz(-0.62566775) q[3];
sx q[3];
rz(0.46712671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59812349) q[2];
sx q[2];
rz(-1.8757966) q[2];
sx q[2];
rz(0.050696105) q[2];
rz(-3.0905881) q[3];
sx q[3];
rz(-1.8404605) q[3];
sx q[3];
rz(0.72371975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26492587) q[0];
sx q[0];
rz(-0.51487881) q[0];
sx q[0];
rz(-1.6616954) q[0];
rz(1.121608) q[1];
sx q[1];
rz(-1.2622967) q[1];
sx q[1];
rz(0.46953896) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74140841) q[0];
sx q[0];
rz(-0.11708507) q[0];
sx q[0];
rz(-1.1598806) q[0];
rz(-pi) q[1];
rz(-1.2079427) q[2];
sx q[2];
rz(-0.41703014) q[2];
sx q[2];
rz(0.94379683) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4457629) q[1];
sx q[1];
rz(-1.1053021) q[1];
sx q[1];
rz(0.19449921) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.77818762) q[3];
sx q[3];
rz(-2.0594849) q[3];
sx q[3];
rz(-0.60569872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7349714) q[2];
sx q[2];
rz(-1.5903641) q[2];
sx q[2];
rz(2.0242019) q[2];
rz(1.0920852) q[3];
sx q[3];
rz(-1.7376244) q[3];
sx q[3];
rz(-0.4357858) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0206873) q[0];
sx q[0];
rz(-1.1130604) q[0];
sx q[0];
rz(-3.0764965) q[0];
rz(-0.33196017) q[1];
sx q[1];
rz(-1.3061701) q[1];
sx q[1];
rz(-1.8665705) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9443566) q[0];
sx q[0];
rz(-1.5612537) q[0];
sx q[0];
rz(1.3413603) q[0];
rz(-1.8529231) q[2];
sx q[2];
rz(-0.45860386) q[2];
sx q[2];
rz(2.0676072) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3447977) q[1];
sx q[1];
rz(-1.5196504) q[1];
sx q[1];
rz(-0.68294345) q[1];
rz(-pi) q[2];
rz(0.17586768) q[3];
sx q[3];
rz(-0.81070903) q[3];
sx q[3];
rz(0.46739235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.28978213) q[2];
sx q[2];
rz(-1.6921356) q[2];
sx q[2];
rz(0.97274485) q[2];
rz(-2.4408477) q[3];
sx q[3];
rz(-2.2795129) q[3];
sx q[3];
rz(2.3166482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4389316) q[0];
sx q[0];
rz(-2.7362566) q[0];
sx q[0];
rz(2.7324556) q[0];
rz(1.9955955) q[1];
sx q[1];
rz(-1.0640249) q[1];
sx q[1];
rz(-1.2629868) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43009297) q[0];
sx q[0];
rz(-1.5709251) q[0];
sx q[0];
rz(1.5789728) q[0];
rz(-pi) q[1];
rz(-1.1767597) q[2];
sx q[2];
rz(-1.3140643) q[2];
sx q[2];
rz(0.22025157) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2330875) q[1];
sx q[1];
rz(-0.89492866) q[1];
sx q[1];
rz(0.29260537) q[1];
rz(0.13554445) q[3];
sx q[3];
rz(-2.5672847) q[3];
sx q[3];
rz(1.3804466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24141773) q[2];
sx q[2];
rz(-1.2215542) q[2];
sx q[2];
rz(-2.9202666) q[2];
rz(-0.48702249) q[3];
sx q[3];
rz(-0.86921391) q[3];
sx q[3];
rz(1.937449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0814334) q[0];
sx q[0];
rz(-0.24516036) q[0];
sx q[0];
rz(1.2231476) q[0];
rz(-3.0606015) q[1];
sx q[1];
rz(-1.6209737) q[1];
sx q[1];
rz(-1.5222668) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8318849) q[0];
sx q[0];
rz(-1.4060347) q[0];
sx q[0];
rz(-0.33640961) q[0];
rz(2.4967983) q[2];
sx q[2];
rz(-0.83907467) q[2];
sx q[2];
rz(3.0562468) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6586858) q[1];
sx q[1];
rz(-1.3750089) q[1];
sx q[1];
rz(-1.2607323) q[1];
rz(-pi) q[2];
rz(-0.49477784) q[3];
sx q[3];
rz(-1.4162546) q[3];
sx q[3];
rz(-3.0734594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0167375) q[2];
sx q[2];
rz(-1.7502461) q[2];
sx q[2];
rz(-0.69046956) q[2];
rz(0.27103439) q[3];
sx q[3];
rz(-0.46509585) q[3];
sx q[3];
rz(-1.5439699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4616213) q[0];
sx q[0];
rz(-0.93628708) q[0];
sx q[0];
rz(-0.17967889) q[0];
rz(2.8996254) q[1];
sx q[1];
rz(-0.89121834) q[1];
sx q[1];
rz(-1.2527087) q[1];
rz(2.7015637) q[2];
sx q[2];
rz(-0.46092214) q[2];
sx q[2];
rz(3.0546247) q[2];
rz(-0.88884066) q[3];
sx q[3];
rz(-1.325229) q[3];
sx q[3];
rz(-1.6390507) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
