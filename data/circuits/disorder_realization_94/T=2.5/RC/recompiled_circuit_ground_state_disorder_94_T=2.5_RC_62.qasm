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
rz(0.35960943) q[1];
sx q[1];
rz(-2.8728027) q[1];
sx q[1];
rz(-0.51528817) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9318207) q[0];
sx q[0];
rz(-1.2270667) q[0];
sx q[0];
rz(-1.7863669) q[0];
x q[1];
rz(0.54606055) q[2];
sx q[2];
rz(-1.8309497) q[2];
sx q[2];
rz(0.85064155) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0019307) q[1];
sx q[1];
rz(-1.7822305) q[1];
sx q[1];
rz(0.40567186) q[1];
rz(-pi) q[2];
rz(2.2622313) q[3];
sx q[3];
rz(-0.9586584) q[3];
sx q[3];
rz(-0.94258756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.23060736) q[2];
sx q[2];
rz(-1.5829986) q[2];
sx q[2];
rz(2.4674463) q[2];
rz(0.19787431) q[3];
sx q[3];
rz(-1.2400235) q[3];
sx q[3];
rz(1.5052634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.823371) q[0];
sx q[0];
rz(-2.8405393) q[0];
sx q[0];
rz(0.2163042) q[0];
rz(3.0220616) q[1];
sx q[1];
rz(-2.0152338) q[1];
sx q[1];
rz(-1.7200039) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57935539) q[0];
sx q[0];
rz(-2.6174712) q[0];
sx q[0];
rz(1.7299132) q[0];
rz(-2.6720409) q[2];
sx q[2];
rz(-0.44347635) q[2];
sx q[2];
rz(-1.8592905) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6232963) q[1];
sx q[1];
rz(-0.5591679) q[1];
sx q[1];
rz(-2.1914047) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.64038527) q[3];
sx q[3];
rz(-0.83723584) q[3];
sx q[3];
rz(-3.0471826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.30430946) q[2];
sx q[2];
rz(-1.9459566) q[2];
sx q[2];
rz(-2.1514814) q[2];
rz(2.385251) q[3];
sx q[3];
rz(-2.4939311) q[3];
sx q[3];
rz(-2.9653449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9847617) q[0];
sx q[0];
rz(-2.3909843) q[0];
sx q[0];
rz(2.005715) q[0];
rz(2.1319977) q[1];
sx q[1];
rz(-1.1547836) q[1];
sx q[1];
rz(-0.44581595) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85161457) q[0];
sx q[0];
rz(-1.5015242) q[0];
sx q[0];
rz(1.4659856) q[0];
rz(-pi) q[1];
x q[1];
rz(0.03568825) q[2];
sx q[2];
rz(-0.70377398) q[2];
sx q[2];
rz(1.4427217) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7103084) q[1];
sx q[1];
rz(-2.388622) q[1];
sx q[1];
rz(-2.1953039) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0296311) q[3];
sx q[3];
rz(-0.44977934) q[3];
sx q[3];
rz(1.2686319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3710215) q[2];
sx q[2];
rz(-2.1393445) q[2];
sx q[2];
rz(0.64024964) q[2];
rz(-0.69784969) q[3];
sx q[3];
rz(-1.4638487) q[3];
sx q[3];
rz(2.8016134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.8401538) q[0];
sx q[0];
rz(-2.2541663) q[0];
sx q[0];
rz(1.3288757) q[0];
rz(-1.2170732) q[1];
sx q[1];
rz(-2.4220059) q[1];
sx q[1];
rz(2.365239) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6049801) q[0];
sx q[0];
rz(-2.1663499) q[0];
sx q[0];
rz(-1.1185297) q[0];
rz(-pi) q[1];
rz(2.4716581) q[2];
sx q[2];
rz(-2.4342854) q[2];
sx q[2];
rz(-2.7279127) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.94521071) q[1];
sx q[1];
rz(-0.27881611) q[1];
sx q[1];
rz(1.5295188) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4569974) q[3];
sx q[3];
rz(-1.5696313) q[3];
sx q[3];
rz(2.5161285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.85294574) q[2];
sx q[2];
rz(-0.32175803) q[2];
sx q[2];
rz(-1.2159489) q[2];
rz(1.1207885) q[3];
sx q[3];
rz(-1.8782764) q[3];
sx q[3];
rz(-0.88302511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23585606) q[0];
sx q[0];
rz(-0.97705066) q[0];
sx q[0];
rz(-2.6781154) q[0];
rz(-2.0526759) q[1];
sx q[1];
rz(-1.2605647) q[1];
sx q[1];
rz(1.3190528) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4396297) q[0];
sx q[0];
rz(-1.3101398) q[0];
sx q[0];
rz(-1.8219276) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7402788) q[2];
sx q[2];
rz(-2.3613644) q[2];
sx q[2];
rz(-1.4223157) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0214349) q[1];
sx q[1];
rz(-1.6579227) q[1];
sx q[1];
rz(-0.98766649) q[1];
x q[2];
rz(-0.18089862) q[3];
sx q[3];
rz(-1.9836863) q[3];
sx q[3];
rz(-0.74821405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9702381) q[2];
sx q[2];
rz(-0.78533185) q[2];
sx q[2];
rz(0.63924092) q[2];
rz(-2.3302737) q[3];
sx q[3];
rz(-2.6526484) q[3];
sx q[3];
rz(-1.5343687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(0.92523471) q[0];
sx q[0];
rz(-2.7282867) q[0];
sx q[0];
rz(-0.21892029) q[0];
rz(1.9256598) q[1];
sx q[1];
rz(-2.6775807) q[1];
sx q[1];
rz(-1.6485515) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0279044) q[0];
sx q[0];
rz(-1.1104212) q[0];
sx q[0];
rz(1.7229863) q[0];
rz(1.5677388) q[2];
sx q[2];
rz(-1.1532461) q[2];
sx q[2];
rz(-1.9034346) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.99438019) q[1];
sx q[1];
rz(-1.5541422) q[1];
sx q[1];
rz(2.4092595) q[1];
x q[2];
rz(0.031367112) q[3];
sx q[3];
rz(-3.0423954) q[3];
sx q[3];
rz(-1.8435696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1047989) q[2];
sx q[2];
rz(-1.7918189) q[2];
sx q[2];
rz(-1.3118504) q[2];
rz(1.5051684) q[3];
sx q[3];
rz(-1.8421831) q[3];
sx q[3];
rz(-2.6817491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63603193) q[0];
sx q[0];
rz(-0.4902896) q[0];
sx q[0];
rz(1.7946515) q[0];
rz(-1.3268283) q[1];
sx q[1];
rz(-2.3074) q[1];
sx q[1];
rz(2.7562275) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.037212278) q[0];
sx q[0];
rz(-1.2118441) q[0];
sx q[0];
rz(1.5370374) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48882882) q[2];
sx q[2];
rz(-2.0483453) q[2];
sx q[2];
rz(1.0350641) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3431541) q[1];
sx q[1];
rz(-2.0166409) q[1];
sx q[1];
rz(1.5413766) q[1];
x q[2];
rz(-0.74968289) q[3];
sx q[3];
rz(-2.5307641) q[3];
sx q[3];
rz(1.8531829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4555326) q[2];
sx q[2];
rz(-0.81642381) q[2];
sx q[2];
rz(1.653999) q[2];
rz(-0.16573302) q[3];
sx q[3];
rz(-2.3159852) q[3];
sx q[3];
rz(-2.9674271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1626749) q[0];
sx q[0];
rz(-0.62541494) q[0];
sx q[0];
rz(0.22853525) q[0];
rz(0.31271115) q[1];
sx q[1];
rz(-2.2622908) q[1];
sx q[1];
rz(1.3794587) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7373567) q[0];
sx q[0];
rz(-1.1690869) q[0];
sx q[0];
rz(1.4866923) q[0];
rz(-1.3993456) q[2];
sx q[2];
rz(-1.8905235) q[2];
sx q[2];
rz(-0.4188183) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.91614281) q[1];
sx q[1];
rz(-1.6054549) q[1];
sx q[1];
rz(-1.9198656) q[1];
x q[2];
rz(0.58016915) q[3];
sx q[3];
rz(-2.3456367) q[3];
sx q[3];
rz(-0.43886504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.10869965) q[2];
sx q[2];
rz(-1.0445107) q[2];
sx q[2];
rz(-2.1006987) q[2];
rz(-0.4852455) q[3];
sx q[3];
rz(-1.3121366) q[3];
sx q[3];
rz(-0.41788873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8682206) q[0];
sx q[0];
rz(-2.1594248) q[0];
sx q[0];
rz(-0.81106538) q[0];
rz(1.9048994) q[1];
sx q[1];
rz(-2.2752454) q[1];
sx q[1];
rz(-2.9467357) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9962353) q[0];
sx q[0];
rz(-2.7872017) q[0];
sx q[0];
rz(2.3217391) q[0];
x q[1];
rz(-2.2213908) q[2];
sx q[2];
rz(-2.0363931) q[2];
sx q[2];
rz(1.0126142) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.81842717) q[1];
sx q[1];
rz(-1.1242529) q[1];
sx q[1];
rz(2.7592804) q[1];
rz(-pi) q[2];
rz(-2.2556188) q[3];
sx q[3];
rz(-2.3625604) q[3];
sx q[3];
rz(-2.1942735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1607232) q[2];
sx q[2];
rz(-1.2474493) q[2];
sx q[2];
rz(-0.54171872) q[2];
rz(1.8187652) q[3];
sx q[3];
rz(-1.8018689) q[3];
sx q[3];
rz(-0.26255783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5498891) q[0];
sx q[0];
rz(-1.5174958) q[0];
sx q[0];
rz(2.8926335) q[0];
rz(1.2173563) q[1];
sx q[1];
rz(-1.267642) q[1];
sx q[1];
rz(0.41044661) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1730984) q[0];
sx q[0];
rz(-1.1652892) q[0];
sx q[0];
rz(-1.9963229) q[0];
rz(-pi) q[1];
rz(-1.4798231) q[2];
sx q[2];
rz(-1.5814648) q[2];
sx q[2];
rz(1.4271229) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1775248) q[1];
sx q[1];
rz(-0.52492889) q[1];
sx q[1];
rz(-2.0711511) q[1];
rz(2.0759567) q[3];
sx q[3];
rz(-2.1200051) q[3];
sx q[3];
rz(0.68699902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.88811389) q[2];
sx q[2];
rz(-1.1682744) q[2];
sx q[2];
rz(-1.1178364) q[2];
rz(-0.11876336) q[3];
sx q[3];
rz(-0.6898841) q[3];
sx q[3];
rz(1.2004872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73040199) q[0];
sx q[0];
rz(-1.406519) q[0];
sx q[0];
rz(2.7997959) q[0];
rz(0.047601184) q[1];
sx q[1];
rz(-1.2536512) q[1];
sx q[1];
rz(1.8258078) q[1];
rz(2.32687) q[2];
sx q[2];
rz(-1.2263032) q[2];
sx q[2];
rz(-2.6487614) q[2];
rz(2.3948492) q[3];
sx q[3];
rz(-2.0604196) q[3];
sx q[3];
rz(-1.4581791) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
