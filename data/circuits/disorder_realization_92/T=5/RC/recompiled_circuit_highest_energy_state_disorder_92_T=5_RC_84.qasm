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
rz(-1.5795213) q[0];
sx q[0];
rz(-3.002394) q[0];
sx q[0];
rz(-0.31779274) q[0];
rz(0.7572445) q[1];
sx q[1];
rz(-1.5905967) q[1];
sx q[1];
rz(1.6928033) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23148558) q[0];
sx q[0];
rz(-0.39662374) q[0];
sx q[0];
rz(2.7825732) q[0];
x q[1];
rz(1.4124539) q[2];
sx q[2];
rz(-1.9082365) q[2];
sx q[2];
rz(2.1580163) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6054666) q[1];
sx q[1];
rz(-1.2647293) q[1];
sx q[1];
rz(-1.1121967) q[1];
rz(-pi) q[2];
rz(-3.0516226) q[3];
sx q[3];
rz(-1.1655088) q[3];
sx q[3];
rz(2.1818415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.455787) q[2];
sx q[2];
rz(-2.0152342) q[2];
sx q[2];
rz(-2.8796999) q[2];
rz(0.68771466) q[3];
sx q[3];
rz(-0.7656289) q[3];
sx q[3];
rz(0.31622893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9460816) q[0];
sx q[0];
rz(-1.8571778) q[0];
sx q[0];
rz(-2.2971357) q[0];
rz(-0.63956815) q[1];
sx q[1];
rz(-2.8896459) q[1];
sx q[1];
rz(-1.8738939) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4019011) q[0];
sx q[0];
rz(-0.88205719) q[0];
sx q[0];
rz(-2.9852377) q[0];
rz(-0.19494178) q[2];
sx q[2];
rz(-0.6856853) q[2];
sx q[2];
rz(0.99239381) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.87593791) q[1];
sx q[1];
rz(-0.75232435) q[1];
sx q[1];
rz(-0.0068412555) q[1];
x q[2];
rz(-0.92314641) q[3];
sx q[3];
rz(-0.83469838) q[3];
sx q[3];
rz(-1.0506964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0212448) q[2];
sx q[2];
rz(-0.52560386) q[2];
sx q[2];
rz(0.28908238) q[2];
rz(0.31288475) q[3];
sx q[3];
rz(-0.94066921) q[3];
sx q[3];
rz(-2.5488034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6784356) q[0];
sx q[0];
rz(-2.1364697) q[0];
sx q[0];
rz(0.67331719) q[0];
rz(0.24766573) q[1];
sx q[1];
rz(-1.0177871) q[1];
sx q[1];
rz(-2.8376875) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16557753) q[0];
sx q[0];
rz(-1.0212693) q[0];
sx q[0];
rz(-0.59161341) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7627569) q[2];
sx q[2];
rz(-2.1386563) q[2];
sx q[2];
rz(-3.0936846) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1294692) q[1];
sx q[1];
rz(-0.15145603) q[1];
sx q[1];
rz(1.1726332) q[1];
rz(-0.21411333) q[3];
sx q[3];
rz(-1.9823834) q[3];
sx q[3];
rz(-2.0579091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5690545) q[2];
sx q[2];
rz(-1.7819541) q[2];
sx q[2];
rz(0.016810091) q[2];
rz(-2.861764) q[3];
sx q[3];
rz(-0.40585104) q[3];
sx q[3];
rz(0.60389891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0957606) q[0];
sx q[0];
rz(-2.026676) q[0];
sx q[0];
rz(-1.8810077) q[0];
rz(-0.49651217) q[1];
sx q[1];
rz(-1.9910201) q[1];
sx q[1];
rz(-1.5042492) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53816089) q[0];
sx q[0];
rz(-1.0350845) q[0];
sx q[0];
rz(-0.49749048) q[0];
rz(-pi) q[1];
x q[1];
rz(1.331318) q[2];
sx q[2];
rz(-2.1000855) q[2];
sx q[2];
rz(3.0893507) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1671686) q[1];
sx q[1];
rz(-1.4054055) q[1];
sx q[1];
rz(1.0252293) q[1];
rz(-pi) q[2];
rz(1.8674866) q[3];
sx q[3];
rz(-2.0540576) q[3];
sx q[3];
rz(-1.7532574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0513976) q[2];
sx q[2];
rz(-0.53956705) q[2];
sx q[2];
rz(2.9810193) q[2];
rz(1.4422902) q[3];
sx q[3];
rz(-1.620404) q[3];
sx q[3];
rz(-3.1062533) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2546688) q[0];
sx q[0];
rz(-1.897432) q[0];
sx q[0];
rz(2.1723893) q[0];
rz(-1.8907549) q[1];
sx q[1];
rz(-2.4617317) q[1];
sx q[1];
rz(2.9108237) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8849205) q[0];
sx q[0];
rz(-1.6065336) q[0];
sx q[0];
rz(-0.0053832609) q[0];
rz(-1.5838301) q[2];
sx q[2];
rz(-0.37838867) q[2];
sx q[2];
rz(2.3005444) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.88420743) q[1];
sx q[1];
rz(-1.1697959) q[1];
sx q[1];
rz(-2.0805418) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1766566) q[3];
sx q[3];
rz(-1.1232383) q[3];
sx q[3];
rz(0.42060619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1988924) q[2];
sx q[2];
rz(-0.57622272) q[2];
sx q[2];
rz(0.25299859) q[2];
rz(-1.7871855) q[3];
sx q[3];
rz(-2.0372882) q[3];
sx q[3];
rz(1.838292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81786466) q[0];
sx q[0];
rz(-2.4619894) q[0];
sx q[0];
rz(0.060081765) q[0];
rz(2.5443063) q[1];
sx q[1];
rz(-0.92514721) q[1];
sx q[1];
rz(3.0937321) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56592076) q[0];
sx q[0];
rz(-1.4788155) q[0];
sx q[0];
rz(0.015181294) q[0];
rz(-pi) q[1];
rz(0.20345236) q[2];
sx q[2];
rz(-1.5566388) q[2];
sx q[2];
rz(-2.3490273) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.3678494) q[1];
sx q[1];
rz(-1.7715104) q[1];
sx q[1];
rz(1.2307005) q[1];
rz(-pi) q[2];
rz(-0.23714785) q[3];
sx q[3];
rz(-2.811069) q[3];
sx q[3];
rz(2.046409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.73703274) q[2];
sx q[2];
rz(-0.76475778) q[2];
sx q[2];
rz(-2.9743189) q[2];
rz(3.0625999) q[3];
sx q[3];
rz(-1.9588214) q[3];
sx q[3];
rz(2.3110068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1792184) q[0];
sx q[0];
rz(-3.0308864) q[0];
sx q[0];
rz(3.0141444) q[0];
rz(0.92370644) q[1];
sx q[1];
rz(-0.5674924) q[1];
sx q[1];
rz(2.4097402) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8475721) q[0];
sx q[0];
rz(-2.5811727) q[0];
sx q[0];
rz(-2.8655478) q[0];
rz(-pi) q[1];
x q[1];
rz(0.42917128) q[2];
sx q[2];
rz(-1.8597684) q[2];
sx q[2];
rz(1.8565551) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4036847) q[1];
sx q[1];
rz(-2.2353454) q[1];
sx q[1];
rz(-0.43364201) q[1];
x q[2];
rz(-0.71890383) q[3];
sx q[3];
rz(-1.2205924) q[3];
sx q[3];
rz(1.7065136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.38714108) q[2];
sx q[2];
rz(-1.0820505) q[2];
sx q[2];
rz(0.88073909) q[2];
rz(0.292101) q[3];
sx q[3];
rz(-0.62551522) q[3];
sx q[3];
rz(2.1706799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77325118) q[0];
sx q[0];
rz(-1.0057978) q[0];
sx q[0];
rz(-0.7811501) q[0];
rz(-1.7225522) q[1];
sx q[1];
rz(-0.96839372) q[1];
sx q[1];
rz(2.9636545) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9289293) q[0];
sx q[0];
rz(-0.50422445) q[0];
sx q[0];
rz(-2.2480856) q[0];
rz(-pi) q[1];
rz(2.0974085) q[2];
sx q[2];
rz(-1.9810314) q[2];
sx q[2];
rz(0.1178478) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9330703) q[1];
sx q[1];
rz(-2.0561777) q[1];
sx q[1];
rz(0.066047864) q[1];
rz(-0.61011647) q[3];
sx q[3];
rz(-2.8233006) q[3];
sx q[3];
rz(1.3194609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.51131821) q[2];
sx q[2];
rz(-1.7061468) q[2];
sx q[2];
rz(-1.1159631) q[2];
rz(2.337194) q[3];
sx q[3];
rz(-0.65785995) q[3];
sx q[3];
rz(0.72949725) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46111742) q[0];
sx q[0];
rz(-2.456433) q[0];
sx q[0];
rz(-0.53801909) q[0];
rz(2.8021725) q[1];
sx q[1];
rz(-1.5433886) q[1];
sx q[1];
rz(-0.0433878) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2983673) q[0];
sx q[0];
rz(-1.4740254) q[0];
sx q[0];
rz(0.065573905) q[0];
rz(-pi) q[1];
rz(-1.814079) q[2];
sx q[2];
rz(-2.5777811) q[2];
sx q[2];
rz(0.76432894) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0108569) q[1];
sx q[1];
rz(-2.3191575) q[1];
sx q[1];
rz(-1.6793516) q[1];
x q[2];
rz(-2.5424388) q[3];
sx q[3];
rz(-1.8945005) q[3];
sx q[3];
rz(-2.168098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1084719) q[2];
sx q[2];
rz(-2.7816009) q[2];
sx q[2];
rz(-1.6610422) q[2];
rz(-2.8451653) q[3];
sx q[3];
rz(-2.0015643) q[3];
sx q[3];
rz(-0.75291434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4794889) q[0];
sx q[0];
rz(-0.39283735) q[0];
sx q[0];
rz(-1.9774849) q[0];
rz(-2.8083943) q[1];
sx q[1];
rz(-1.6644128) q[1];
sx q[1];
rz(0.747917) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5503905) q[0];
sx q[0];
rz(-1.8699614) q[0];
sx q[0];
rz(-0.62320605) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93846069) q[2];
sx q[2];
rz(-1.399164) q[2];
sx q[2];
rz(-0.29612637) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9063776) q[1];
sx q[1];
rz(-3.0864419) q[1];
sx q[1];
rz(-0.13184594) q[1];
x q[2];
rz(-3.0958067) q[3];
sx q[3];
rz(-0.87821315) q[3];
sx q[3];
rz(-1.847895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.36772874) q[2];
sx q[2];
rz(-1.7325956) q[2];
sx q[2];
rz(1.1665974) q[2];
rz(-1.599132) q[3];
sx q[3];
rz(-1.5365994) q[3];
sx q[3];
rz(-1.0421853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3508956) q[0];
sx q[0];
rz(-0.46115524) q[0];
sx q[0];
rz(-0.3140558) q[0];
rz(3.113204) q[1];
sx q[1];
rz(-2.8969565) q[1];
sx q[1];
rz(1.2895186) q[1];
rz(0.97276982) q[2];
sx q[2];
rz(-1.1373873) q[2];
sx q[2];
rz(1.6349229) q[2];
rz(0.52846626) q[3];
sx q[3];
rz(-2.2615643) q[3];
sx q[3];
rz(-0.27558358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
