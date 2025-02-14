OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.021304) q[0];
sx q[0];
rz(2.6074183) q[0];
sx q[0];
rz(8.3745126) q[0];
rz(1.8771111) q[1];
sx q[1];
rz(-0.85327947) q[1];
sx q[1];
rz(0.54860151) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34615883) q[0];
sx q[0];
rz(-1.2761269) q[0];
sx q[0];
rz(2.8351889) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5496375) q[2];
sx q[2];
rz(-2.1680764) q[2];
sx q[2];
rz(2.6279272) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9160936) q[1];
sx q[1];
rz(-1.0584944) q[1];
sx q[1];
rz(-2.6698547) q[1];
rz(2.9012783) q[3];
sx q[3];
rz(-1.263947) q[3];
sx q[3];
rz(3.0003529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7301664) q[2];
sx q[2];
rz(-2.7226518) q[2];
sx q[2];
rz(-2.1298998) q[2];
rz(-0.88459477) q[3];
sx q[3];
rz(-1.9832289) q[3];
sx q[3];
rz(1.2456892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0678134) q[0];
sx q[0];
rz(-0.99869204) q[0];
sx q[0];
rz(0.78773898) q[0];
rz(2.1630321) q[1];
sx q[1];
rz(-1.1509044) q[1];
sx q[1];
rz(-1.8914793) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.387991) q[0];
sx q[0];
rz(-1.4630254) q[0];
sx q[0];
rz(-2.0600256) q[0];
rz(0.27898326) q[2];
sx q[2];
rz(-1.9595385) q[2];
sx q[2];
rz(1.1065567) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.7693705) q[1];
sx q[1];
rz(-1.5466855) q[1];
sx q[1];
rz(-0.80609821) q[1];
rz(0.21997178) q[3];
sx q[3];
rz(-1.7860054) q[3];
sx q[3];
rz(0.53088354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1193739) q[2];
sx q[2];
rz(-1.4639857) q[2];
sx q[2];
rz(-3.019943) q[2];
rz(-0.71074784) q[3];
sx q[3];
rz(-0.23356479) q[3];
sx q[3];
rz(-2.5392883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0372593) q[0];
sx q[0];
rz(-2.4132044) q[0];
sx q[0];
rz(0.65993586) q[0];
rz(-0.51689369) q[1];
sx q[1];
rz(-2.4362502) q[1];
sx q[1];
rz(1.0924115) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4683974) q[0];
sx q[0];
rz(-1.930113) q[0];
sx q[0];
rz(0.27984377) q[0];
rz(-pi) q[1];
rz(-1.2159816) q[2];
sx q[2];
rz(-0.72941226) q[2];
sx q[2];
rz(-1.5340005) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99169774) q[1];
sx q[1];
rz(-2.3058165) q[1];
sx q[1];
rz(-0.81945547) q[1];
rz(-pi) q[2];
rz(-1.8275683) q[3];
sx q[3];
rz(-0.93807546) q[3];
sx q[3];
rz(-0.69958052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2163781) q[2];
sx q[2];
rz(-2.7984012) q[2];
sx q[2];
rz(-1.421831) q[2];
rz(0.4839932) q[3];
sx q[3];
rz(-1.657594) q[3];
sx q[3];
rz(-1.7234195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5928818) q[0];
sx q[0];
rz(-0.084267862) q[0];
sx q[0];
rz(1.3053869) q[0];
rz(0.0072172324) q[1];
sx q[1];
rz(-0.22166285) q[1];
sx q[1];
rz(-0.73297393) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.510988) q[0];
sx q[0];
rz(-2.9633491) q[0];
sx q[0];
rz(0.83971881) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.29455955) q[2];
sx q[2];
rz(-1.4617051) q[2];
sx q[2];
rz(-2.8677169) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.32506714) q[1];
sx q[1];
rz(-2.1754334) q[1];
sx q[1];
rz(-0.52765347) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6812612) q[3];
sx q[3];
rz(-2.2909431) q[3];
sx q[3];
rz(0.11321774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2923773) q[2];
sx q[2];
rz(-1.4039618) q[2];
sx q[2];
rz(0.89548573) q[2];
rz(0.53564566) q[3];
sx q[3];
rz(-1.5786542) q[3];
sx q[3];
rz(-3.079788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8931005) q[0];
sx q[0];
rz(-2.94815) q[0];
sx q[0];
rz(-2.6336811) q[0];
rz(-1.4503362) q[1];
sx q[1];
rz(-2.129887) q[1];
sx q[1];
rz(-1.3571665) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67726595) q[0];
sx q[0];
rz(-1.6826165) q[0];
sx q[0];
rz(-1.6202116) q[0];
rz(-pi) q[1];
rz(0.72765784) q[2];
sx q[2];
rz(-0.59077677) q[2];
sx q[2];
rz(2.1075005) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.60496441) q[1];
sx q[1];
rz(-1.9603029) q[1];
sx q[1];
rz(3.0252181) q[1];
rz(2.100222) q[3];
sx q[3];
rz(-0.66877194) q[3];
sx q[3];
rz(-3.1094299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.42673972) q[2];
sx q[2];
rz(-1.4740976) q[2];
sx q[2];
rz(-1.1023785) q[2];
rz(1.1357931) q[3];
sx q[3];
rz(-1.2165242) q[3];
sx q[3];
rz(2.3701325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8847467) q[0];
sx q[0];
rz(-0.9698292) q[0];
sx q[0];
rz(0.92887512) q[0];
rz(2.7595787) q[1];
sx q[1];
rz(-0.99645749) q[1];
sx q[1];
rz(1.6928203) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9858157) q[0];
sx q[0];
rz(-2.563463) q[0];
sx q[0];
rz(0.57684071) q[0];
rz(-0.11547757) q[2];
sx q[2];
rz(-2.653947) q[2];
sx q[2];
rz(-2.6754926) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1898489) q[1];
sx q[1];
rz(-0.41297022) q[1];
sx q[1];
rz(-1.8443405) q[1];
rz(-pi) q[2];
rz(1.153451) q[3];
sx q[3];
rz(-0.40606895) q[3];
sx q[3];
rz(-1.690133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.53943071) q[2];
sx q[2];
rz(-0.346589) q[2];
sx q[2];
rz(0.76751417) q[2];
rz(-1.2933939) q[3];
sx q[3];
rz(-0.69197217) q[3];
sx q[3];
rz(0.20492157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9354189) q[0];
sx q[0];
rz(-2.967301) q[0];
sx q[0];
rz(-0.25099227) q[0];
rz(-0.17768606) q[1];
sx q[1];
rz(-1.6975941) q[1];
sx q[1];
rz(-0.65690717) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85992438) q[0];
sx q[0];
rz(-2.826068) q[0];
sx q[0];
rz(0.84976999) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7288293) q[2];
sx q[2];
rz(-2.2185746) q[2];
sx q[2];
rz(-0.72731804) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6823947) q[1];
sx q[1];
rz(-2.6434745) q[1];
sx q[1];
rz(2.47704) q[1];
x q[2];
rz(-2.5890089) q[3];
sx q[3];
rz(-0.72425084) q[3];
sx q[3];
rz(-2.9720705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.80829197) q[2];
sx q[2];
rz(-0.57922816) q[2];
sx q[2];
rz(-2.2283238) q[2];
rz(-2.463786) q[3];
sx q[3];
rz(-1.6401688) q[3];
sx q[3];
rz(1.0031797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.102757) q[0];
sx q[0];
rz(-1.1433733) q[0];
sx q[0];
rz(-2.5586149) q[0];
rz(1.673117) q[1];
sx q[1];
rz(-0.99248326) q[1];
sx q[1];
rz(-2.800422) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4214246) q[0];
sx q[0];
rz(-2.9377794) q[0];
sx q[0];
rz(-1.782062) q[0];
rz(-pi) q[1];
rz(-1.2874574) q[2];
sx q[2];
rz(-1.8036799) q[2];
sx q[2];
rz(2.8849059) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9205089) q[1];
sx q[1];
rz(-2.608641) q[1];
sx q[1];
rz(-2.3531662) q[1];
rz(-0.40557762) q[3];
sx q[3];
rz(-0.93355191) q[3];
sx q[3];
rz(3.0831856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1950281) q[2];
sx q[2];
rz(-2.3174758) q[2];
sx q[2];
rz(-0.80671802) q[2];
rz(-2.6840456) q[3];
sx q[3];
rz(-0.98743192) q[3];
sx q[3];
rz(-0.88361067) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6954527) q[0];
sx q[0];
rz(-2.0841632) q[0];
sx q[0];
rz(2.9651508) q[0];
rz(2.8314619) q[1];
sx q[1];
rz(-1.4896432) q[1];
sx q[1];
rz(2.2055221) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53545685) q[0];
sx q[0];
rz(-1.8671037) q[0];
sx q[0];
rz(-1.7264051) q[0];
rz(2.5199982) q[2];
sx q[2];
rz(-0.69658579) q[2];
sx q[2];
rz(-2.7736349) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.61221695) q[1];
sx q[1];
rz(-1.9381372) q[1];
sx q[1];
rz(0.099766082) q[1];
x q[2];
rz(-3.0210439) q[3];
sx q[3];
rz(-1.0998187) q[3];
sx q[3];
rz(2.4496743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.04756847) q[2];
sx q[2];
rz(-1.3022283) q[2];
sx q[2];
rz(-0.0085208323) q[2];
rz(-2.408037) q[3];
sx q[3];
rz(-0.80500427) q[3];
sx q[3];
rz(0.12695299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9200639) q[0];
sx q[0];
rz(-2.7286752) q[0];
sx q[0];
rz(-0.76989663) q[0];
rz(0.13433111) q[1];
sx q[1];
rz(-0.7901935) q[1];
sx q[1];
rz(1.2199527) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0095739) q[0];
sx q[0];
rz(-1.2694799) q[0];
sx q[0];
rz(-1.3265885) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6069746) q[2];
sx q[2];
rz(-1.8152945) q[2];
sx q[2];
rz(-0.90842694) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7999477) q[1];
sx q[1];
rz(-1.1468256) q[1];
sx q[1];
rz(1.9862224) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0073529) q[3];
sx q[3];
rz(-0.16664342) q[3];
sx q[3];
rz(1.327001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8613345) q[2];
sx q[2];
rz(-2.4836149) q[2];
sx q[2];
rz(-2.3560143) q[2];
rz(-2.806459) q[3];
sx q[3];
rz(-2.1527055) q[3];
sx q[3];
rz(-0.82211632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32421865) q[0];
sx q[0];
rz(-2.2786409) q[0];
sx q[0];
rz(1.8550158) q[0];
rz(0.68589504) q[1];
sx q[1];
rz(-0.58882014) q[1];
sx q[1];
rz(-1.3480766) q[1];
rz(0.98966148) q[2];
sx q[2];
rz(-1.5753395) q[2];
sx q[2];
rz(0.01103845) q[2];
rz(1.4487063) q[3];
sx q[3];
rz(-0.29480903) q[3];
sx q[3];
rz(2.8041425) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
