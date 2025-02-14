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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7069076) q[0];
sx q[0];
rz(-1.3680172) q[0];
sx q[0];
rz(0.35122706) q[0];
x q[1];
rz(1.2688387) q[2];
sx q[2];
rz(-2.0965323) q[2];
sx q[2];
rz(2.5765004) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2559515) q[1];
sx q[1];
rz(-0.4547387) q[1];
sx q[1];
rz(-0.498147) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87936136) q[3];
sx q[3];
rz(-2.1829343) q[3];
sx q[3];
rz(-2.1990051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9109853) q[2];
sx q[2];
rz(-1.5829986) q[2];
sx q[2];
rz(-0.67414635) q[2];
rz(2.9437183) q[3];
sx q[3];
rz(-1.2400235) q[3];
sx q[3];
rz(-1.5052634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.823371) q[0];
sx q[0];
rz(-2.8405393) q[0];
sx q[0];
rz(2.9252885) q[0];
rz(3.0220616) q[1];
sx q[1];
rz(-2.0152338) q[1];
sx q[1];
rz(1.4215887) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0121032) q[0];
sx q[0];
rz(-1.6501745) q[0];
sx q[0];
rz(1.0521655) q[0];
rz(-pi) q[1];
rz(2.6720409) q[2];
sx q[2];
rz(-0.44347635) q[2];
sx q[2];
rz(-1.2823022) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9227289) q[1];
sx q[1];
rz(-1.1245755) q[1];
sx q[1];
rz(2.7925744) q[1];
rz(-0.72705468) q[3];
sx q[3];
rz(-2.0306572) q[3];
sx q[3];
rz(1.0136295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8372832) q[2];
sx q[2];
rz(-1.195636) q[2];
sx q[2];
rz(-2.1514814) q[2];
rz(-2.385251) q[3];
sx q[3];
rz(-2.4939311) q[3];
sx q[3];
rz(2.9653449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9847617) q[0];
sx q[0];
rz(-2.3909843) q[0];
sx q[0];
rz(-1.1358776) q[0];
rz(1.0095949) q[1];
sx q[1];
rz(-1.986809) q[1];
sx q[1];
rz(-0.44581595) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71190055) q[0];
sx q[0];
rz(-1.4662379) q[0];
sx q[0];
rz(0.069653102) q[0];
rz(-0.70345975) q[2];
sx q[2];
rz(-1.5477053) q[2];
sx q[2];
rz(-2.9863043) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.793632) q[1];
sx q[1];
rz(-2.1588481) q[1];
sx q[1];
rz(0.5012725) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11196158) q[3];
sx q[3];
rz(-2.6918133) q[3];
sx q[3];
rz(-1.2686319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.77057114) q[2];
sx q[2];
rz(-2.1393445) q[2];
sx q[2];
rz(-2.501343) q[2];
rz(0.69784969) q[3];
sx q[3];
rz(-1.4638487) q[3];
sx q[3];
rz(0.33997926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
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
rz(0.3014389) q[0];
sx q[0];
rz(-2.2541663) q[0];
sx q[0];
rz(-1.812717) q[0];
rz(-1.2170732) q[1];
sx q[1];
rz(-0.71958676) q[1];
sx q[1];
rz(-2.365239) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76809873) q[0];
sx q[0];
rz(-1.9409618) q[0];
sx q[0];
rz(0.64565701) q[0];
rz(-pi) q[1];
rz(-0.59036915) q[2];
sx q[2];
rz(-1.1554829) q[2];
sx q[2];
rz(-1.4424974) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1963819) q[1];
sx q[1];
rz(-0.27881611) q[1];
sx q[1];
rz(-1.6120738) q[1];
rz(3.1397503) q[3];
sx q[3];
rz(-0.68459612) q[3];
sx q[3];
rz(2.1948333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.2886469) q[2];
sx q[2];
rz(-0.32175803) q[2];
sx q[2];
rz(1.2159489) q[2];
rz(-1.1207885) q[3];
sx q[3];
rz(-1.8782764) q[3];
sx q[3];
rz(-2.2585675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9057366) q[0];
sx q[0];
rz(-0.97705066) q[0];
sx q[0];
rz(2.6781154) q[0];
rz(-1.0889168) q[1];
sx q[1];
rz(-1.2605647) q[1];
sx q[1];
rz(-1.3190528) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3387793) q[0];
sx q[0];
rz(-1.3283214) q[0];
sx q[0];
rz(-2.8728897) q[0];
x q[1];
rz(-1.7402788) q[2];
sx q[2];
rz(-2.3613644) q[2];
sx q[2];
rz(-1.7192769) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5079759) q[1];
sx q[1];
rz(-2.1514261) q[1];
sx q[1];
rz(3.0373322) q[1];
rz(-pi) q[2];
x q[2];
rz(2.960694) q[3];
sx q[3];
rz(-1.9836863) q[3];
sx q[3];
rz(2.3933786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9702381) q[2];
sx q[2];
rz(-2.3562608) q[2];
sx q[2];
rz(-2.5023517) q[2];
rz(-2.3302737) q[3];
sx q[3];
rz(-0.48894426) q[3];
sx q[3];
rz(-1.607224) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2163579) q[0];
sx q[0];
rz(-2.7282867) q[0];
sx q[0];
rz(0.21892029) q[0];
rz(1.9256598) q[1];
sx q[1];
rz(-0.464012) q[1];
sx q[1];
rz(-1.4930412) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7812778) q[0];
sx q[0];
rz(-2.6584315) q[0];
sx q[0];
rz(0.29668087) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0068917787) q[2];
sx q[2];
rz(-0.41756072) q[2];
sx q[2];
rz(1.895895) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.1472125) q[1];
sx q[1];
rz(-1.5874505) q[1];
sx q[1];
rz(2.4092595) q[1];
rz(-0.099148765) q[3];
sx q[3];
rz(-1.5676904) q[3];
sx q[3];
rz(2.9000324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1047989) q[2];
sx q[2];
rz(-1.3497738) q[2];
sx q[2];
rz(1.3118504) q[2];
rz(-1.6364243) q[3];
sx q[3];
rz(-1.8421831) q[3];
sx q[3];
rz(0.45984355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5055607) q[0];
sx q[0];
rz(-0.4902896) q[0];
sx q[0];
rz(1.7946515) q[0];
rz(-1.3268283) q[1];
sx q[1];
rz(-2.3074) q[1];
sx q[1];
rz(-0.38536513) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1043804) q[0];
sx q[0];
rz(-1.9297486) q[0];
sx q[0];
rz(1.6045553) q[0];
rz(2.1009675) q[2];
sx q[2];
rz(-1.1405924) q[2];
sx q[2];
rz(2.8456147) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7984385) q[1];
sx q[1];
rz(-2.0166409) q[1];
sx q[1];
rz(-1.5413766) q[1];
rz(-pi) q[2];
rz(-0.74968289) q[3];
sx q[3];
rz(-2.5307641) q[3];
sx q[3];
rz(-1.2884097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.68606004) q[2];
sx q[2];
rz(-0.81642381) q[2];
sx q[2];
rz(1.4875937) q[2];
rz(-0.16573302) q[3];
sx q[3];
rz(-0.82560742) q[3];
sx q[3];
rz(-0.17416557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97891775) q[0];
sx q[0];
rz(-0.62541494) q[0];
sx q[0];
rz(2.9130574) q[0];
rz(-0.31271115) q[1];
sx q[1];
rz(-0.87930185) q[1];
sx q[1];
rz(-1.762134) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7373567) q[0];
sx q[0];
rz(-1.9725058) q[0];
sx q[0];
rz(-1.6549003) q[0];
x q[1];
rz(-2.6657728) q[2];
sx q[2];
rz(-0.36140051) q[2];
sx q[2];
rz(0.084712383) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.55973616) q[1];
sx q[1];
rz(-2.7908771) q[1];
sx q[1];
rz(-1.4697671) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4346541) q[3];
sx q[3];
rz(-1.1683373) q[3];
sx q[3];
rz(-1.5618531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.032893) q[2];
sx q[2];
rz(-2.097082) q[2];
sx q[2];
rz(1.0408939) q[2];
rz(-2.6563472) q[3];
sx q[3];
rz(-1.8294561) q[3];
sx q[3];
rz(-0.41788873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27337209) q[0];
sx q[0];
rz(-2.1594248) q[0];
sx q[0];
rz(2.3305273) q[0];
rz(1.2366933) q[1];
sx q[1];
rz(-2.2752454) q[1];
sx q[1];
rz(-0.19485697) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14535739) q[0];
sx q[0];
rz(-0.35439098) q[0];
sx q[0];
rz(2.3217391) q[0];
rz(0.92020184) q[2];
sx q[2];
rz(-2.0363931) q[2];
sx q[2];
rz(-2.1289785) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.58043805) q[1];
sx q[1];
rz(-1.2276137) q[1];
sx q[1];
rz(2.047206) q[1];
rz(-pi) q[2];
rz(0.88597383) q[3];
sx q[3];
rz(-0.77903226) q[3];
sx q[3];
rz(-0.94731915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.98086944) q[2];
sx q[2];
rz(-1.2474493) q[2];
sx q[2];
rz(0.54171872) q[2];
rz(1.3228275) q[3];
sx q[3];
rz(-1.3397237) q[3];
sx q[3];
rz(2.8790348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.5498891) q[0];
sx q[0];
rz(-1.6240969) q[0];
sx q[0];
rz(-0.24895915) q[0];
rz(1.9242363) q[1];
sx q[1];
rz(-1.8739506) q[1];
sx q[1];
rz(-2.731146) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7792203) q[0];
sx q[0];
rz(-1.9598613) q[0];
sx q[0];
rz(-0.44045191) q[0];
rz(-pi) q[1];
rz(-1.4798231) q[2];
sx q[2];
rz(-1.5601279) q[2];
sx q[2];
rz(-1.4271229) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9767993) q[1];
sx q[1];
rz(-1.327997) q[1];
sx q[1];
rz(2.040928) q[1];
rz(-pi) q[2];
rz(2.0759567) q[3];
sx q[3];
rz(-1.0215875) q[3];
sx q[3];
rz(2.4545936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2534788) q[2];
sx q[2];
rz(-1.1682744) q[2];
sx q[2];
rz(2.0237563) q[2];
rz(-3.0228293) q[3];
sx q[3];
rz(-0.6898841) q[3];
sx q[3];
rz(1.9411055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73040199) q[0];
sx q[0];
rz(-1.7350736) q[0];
sx q[0];
rz(-0.34179678) q[0];
rz(-3.0939915) q[1];
sx q[1];
rz(-1.2536512) q[1];
sx q[1];
rz(1.8258078) q[1];
rz(0.45817057) q[2];
sx q[2];
rz(-2.2728163) q[2];
sx q[2];
rz(-1.386281) q[2];
rz(-2.1988395) q[3];
sx q[3];
rz(-0.92798622) q[3];
sx q[3];
rz(-0.29792132) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
