OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2296978) q[0];
sx q[0];
rz(-1.8954281) q[0];
sx q[0];
rz(1.5204313) q[0];
rz(-5.9235759) q[1];
sx q[1];
rz(3.4103826) q[1];
sx q[1];
rz(8.9094898) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7069076) q[0];
sx q[0];
rz(-1.3680172) q[0];
sx q[0];
rz(-2.7903656) q[0];
rz(-pi) q[1];
rz(2.6679466) q[2];
sx q[2];
rz(-0.59913991) q[2];
sx q[2];
rz(2.0210217) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.34124247) q[1];
sx q[1];
rz(-1.1746695) q[1];
sx q[1];
rz(1.3413097) q[1];
rz(-pi) q[2];
x q[2];
rz(0.73731976) q[3];
sx q[3];
rz(-2.2529369) q[3];
sx q[3];
rz(1.2350262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9109853) q[2];
sx q[2];
rz(-1.5829986) q[2];
sx q[2];
rz(-2.4674463) q[2];
rz(-2.9437183) q[3];
sx q[3];
rz(-1.2400235) q[3];
sx q[3];
rz(-1.6363293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.823371) q[0];
sx q[0];
rz(-2.8405393) q[0];
sx q[0];
rz(-2.9252885) q[0];
rz(3.0220616) q[1];
sx q[1];
rz(-2.0152338) q[1];
sx q[1];
rz(1.4215887) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5622373) q[0];
sx q[0];
rz(-2.6174712) q[0];
sx q[0];
rz(1.7299132) q[0];
rz(0.40070285) q[2];
sx q[2];
rz(-1.3754015) q[2];
sx q[2];
rz(2.4233482) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2188637) q[1];
sx q[1];
rz(-1.1245755) q[1];
sx q[1];
rz(-2.7925744) q[1];
rz(-pi) q[2];
x q[2];
rz(2.414538) q[3];
sx q[3];
rz(-1.1109354) q[3];
sx q[3];
rz(-1.0136295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8372832) q[2];
sx q[2];
rz(-1.195636) q[2];
sx q[2];
rz(2.1514814) q[2];
rz(2.385251) q[3];
sx q[3];
rz(-2.4939311) q[3];
sx q[3];
rz(-2.9653449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-2.9847617) q[0];
sx q[0];
rz(-0.75060833) q[0];
sx q[0];
rz(2.005715) q[0];
rz(1.0095949) q[1];
sx q[1];
rz(-1.986809) q[1];
sx q[1];
rz(-0.44581595) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2899781) q[0];
sx q[0];
rz(-1.5015242) q[0];
sx q[0];
rz(1.4659856) q[0];
x q[1];
rz(3.1059044) q[2];
sx q[2];
rz(-0.70377398) q[2];
sx q[2];
rz(1.698871) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.517934) q[1];
sx q[1];
rz(-1.9821189) q[1];
sx q[1];
rz(2.2208396) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11196158) q[3];
sx q[3];
rz(-0.44977934) q[3];
sx q[3];
rz(1.2686319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77057114) q[2];
sx q[2];
rz(-1.0022481) q[2];
sx q[2];
rz(-2.501343) q[2];
rz(2.443743) q[3];
sx q[3];
rz(-1.4638487) q[3];
sx q[3];
rz(2.8016134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-2.8401538) q[0];
sx q[0];
rz(-2.2541663) q[0];
sx q[0];
rz(1.3288757) q[0];
rz(-1.9245194) q[1];
sx q[1];
rz(-2.4220059) q[1];
sx q[1];
rz(0.77635366) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76809873) q[0];
sx q[0];
rz(-1.2006309) q[0];
sx q[0];
rz(-0.64565701) q[0];
x q[1];
rz(2.0587875) q[2];
sx q[2];
rz(-2.1052202) q[2];
sx q[2];
rz(-0.39235175) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1963819) q[1];
sx q[1];
rz(-0.27881611) q[1];
sx q[1];
rz(1.5295188) q[1];
rz(-2.4569974) q[3];
sx q[3];
rz(-1.5696313) q[3];
sx q[3];
rz(-0.62546414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.85294574) q[2];
sx q[2];
rz(-2.8198346) q[2];
sx q[2];
rz(1.9256437) q[2];
rz(1.1207885) q[3];
sx q[3];
rz(-1.2633163) q[3];
sx q[3];
rz(-2.2585675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-0.23585606) q[0];
sx q[0];
rz(-2.164542) q[0];
sx q[0];
rz(-2.6781154) q[0];
rz(1.0889168) q[1];
sx q[1];
rz(-1.8810279) q[1];
sx q[1];
rz(-1.3190528) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4850906) q[0];
sx q[0];
rz(-0.35995558) q[0];
sx q[0];
rz(-0.75004049) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4013138) q[2];
sx q[2];
rz(-0.78022829) q[2];
sx q[2];
rz(1.7192769) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0214349) q[1];
sx q[1];
rz(-1.6579227) q[1];
sx q[1];
rz(-2.1539262) q[1];
x q[2];
rz(2.960694) q[3];
sx q[3];
rz(-1.1579063) q[3];
sx q[3];
rz(-2.3933786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9702381) q[2];
sx q[2];
rz(-0.78533185) q[2];
sx q[2];
rz(-0.63924092) q[2];
rz(-0.81131896) q[3];
sx q[3];
rz(-2.6526484) q[3];
sx q[3];
rz(-1.607224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2163579) q[0];
sx q[0];
rz(-0.41330591) q[0];
sx q[0];
rz(-0.21892029) q[0];
rz(-1.2159329) q[1];
sx q[1];
rz(-2.6775807) q[1];
sx q[1];
rz(-1.6485515) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0279044) q[0];
sx q[0];
rz(-2.0311715) q[0];
sx q[0];
rz(-1.7229863) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1347009) q[2];
sx q[2];
rz(-2.7240319) q[2];
sx q[2];
rz(-1.895895) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.56144372) q[1];
sx q[1];
rz(-0.83858788) q[1];
sx q[1];
rz(-1.5931908) q[1];
rz(-pi) q[2];
rz(-1.5676751) q[3];
sx q[3];
rz(-1.6699446) q[3];
sx q[3];
rz(1.8120476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1047989) q[2];
sx q[2];
rz(-1.3497738) q[2];
sx q[2];
rz(-1.8297423) q[2];
rz(-1.5051684) q[3];
sx q[3];
rz(-1.8421831) q[3];
sx q[3];
rz(-0.45984355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63603193) q[0];
sx q[0];
rz(-0.4902896) q[0];
sx q[0];
rz(1.7946515) q[0];
rz(1.3268283) q[1];
sx q[1];
rz(-0.83419269) q[1];
sx q[1];
rz(2.7562275) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6198719) q[0];
sx q[0];
rz(-1.5391897) q[0];
sx q[0];
rz(0.35913976) q[0];
rz(-pi) q[1];
rz(2.1009675) q[2];
sx q[2];
rz(-1.1405924) q[2];
sx q[2];
rz(2.8456147) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3431541) q[1];
sx q[1];
rz(-2.0166409) q[1];
sx q[1];
rz(1.6002161) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1256451) q[3];
sx q[3];
rz(-1.1375918) q[3];
sx q[3];
rz(1.0039745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.68606004) q[2];
sx q[2];
rz(-0.81642381) q[2];
sx q[2];
rz(-1.4875937) q[2];
rz(2.9758596) q[3];
sx q[3];
rz(-2.3159852) q[3];
sx q[3];
rz(0.17416557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97891775) q[0];
sx q[0];
rz(-2.5161777) q[0];
sx q[0];
rz(2.9130574) q[0];
rz(0.31271115) q[1];
sx q[1];
rz(-2.2622908) q[1];
sx q[1];
rz(-1.762134) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1995102) q[0];
sx q[0];
rz(-1.6481912) q[0];
sx q[0];
rz(0.40298526) q[0];
x q[1];
rz(-1.7422471) q[2];
sx q[2];
rz(-1.2510692) q[2];
sx q[2];
rz(-0.4188183) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5818565) q[1];
sx q[1];
rz(-0.35071555) q[1];
sx q[1];
rz(1.4697671) q[1];
rz(-pi) q[2];
rz(-0.58016915) q[3];
sx q[3];
rz(-2.3456367) q[3];
sx q[3];
rz(0.43886504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.032893) q[2];
sx q[2];
rz(-2.097082) q[2];
sx q[2];
rz(-1.0408939) q[2];
rz(-2.6563472) q[3];
sx q[3];
rz(-1.8294561) q[3];
sx q[3];
rz(2.7237039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27337209) q[0];
sx q[0];
rz(-0.98216787) q[0];
sx q[0];
rz(-0.81106538) q[0];
rz(1.9048994) q[1];
sx q[1];
rz(-0.86634723) q[1];
sx q[1];
rz(-0.19485697) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2132414) q[0];
sx q[0];
rz(-1.3143063) q[0];
sx q[0];
rz(-0.24730206) q[0];
x q[1];
rz(2.2213908) q[2];
sx q[2];
rz(-2.0363931) q[2];
sx q[2];
rz(-1.0126142) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.3231655) q[1];
sx q[1];
rz(-2.0173397) q[1];
sx q[1];
rz(0.38231229) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55826346) q[3];
sx q[3];
rz(-0.9953863) q[3];
sx q[3];
rz(-0.093274506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.98086944) q[2];
sx q[2];
rz(-1.8941433) q[2];
sx q[2];
rz(-0.54171872) q[2];
rz(1.3228275) q[3];
sx q[3];
rz(-1.3397237) q[3];
sx q[3];
rz(-0.26255783) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5498891) q[0];
sx q[0];
rz(-1.6240969) q[0];
sx q[0];
rz(-2.8926335) q[0];
rz(-1.2173563) q[1];
sx q[1];
rz(-1.267642) q[1];
sx q[1];
rz(2.731146) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7792203) q[0];
sx q[0];
rz(-1.9598613) q[0];
sx q[0];
rz(-0.44045191) q[0];
rz(-1.6617695) q[2];
sx q[2];
rz(-1.5601279) q[2];
sx q[2];
rz(1.4271229) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.16479334) q[1];
sx q[1];
rz(-1.327997) q[1];
sx q[1];
rz(2.040928) q[1];
x q[2];
rz(2.0759567) q[3];
sx q[3];
rz(-1.0215875) q[3];
sx q[3];
rz(2.4545936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2534788) q[2];
sx q[2];
rz(-1.9733182) q[2];
sx q[2];
rz(-2.0237563) q[2];
rz(-3.0228293) q[3];
sx q[3];
rz(-2.4517086) q[3];
sx q[3];
rz(-1.9411055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4111907) q[0];
sx q[0];
rz(-1.7350736) q[0];
sx q[0];
rz(-0.34179678) q[0];
rz(0.047601184) q[1];
sx q[1];
rz(-1.2536512) q[1];
sx q[1];
rz(1.8258078) q[1];
rz(-1.088935) q[2];
sx q[2];
rz(-0.81648044) q[2];
sx q[2];
rz(-0.73406506) q[2];
rz(-0.66524617) q[3];
sx q[3];
rz(-2.2753297) q[3];
sx q[3];
rz(0.5827502) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
