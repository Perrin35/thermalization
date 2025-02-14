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
rz(0.12198099) q[0];
sx q[0];
rz(-0.044959083) q[0];
sx q[0];
rz(-3.0461351) q[0];
rz(1.8010315) q[1];
sx q[1];
rz(3.0756693) q[1];
sx q[1];
rz(8.8311721) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5779764) q[0];
sx q[0];
rz(-1.9270183) q[0];
sx q[0];
rz(-2.6249159) q[0];
rz(-pi) q[1];
rz(2.4727137) q[2];
sx q[2];
rz(-0.19193412) q[2];
sx q[2];
rz(-0.30489433) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2954457) q[1];
sx q[1];
rz(-0.6653924) q[1];
sx q[1];
rz(2.800368) q[1];
x q[2];
rz(1.0450706) q[3];
sx q[3];
rz(-1.8107256) q[3];
sx q[3];
rz(0.6434427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0815214) q[2];
sx q[2];
rz(-1.8895443) q[2];
sx q[2];
rz(0.56212765) q[2];
rz(-0.33956042) q[3];
sx q[3];
rz(-1.5226676) q[3];
sx q[3];
rz(-2.0558527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18315166) q[0];
sx q[0];
rz(-2.9957132) q[0];
sx q[0];
rz(1.1613783) q[0];
rz(2.9005652) q[1];
sx q[1];
rz(-2.3190505) q[1];
sx q[1];
rz(-0.30311146) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68837092) q[0];
sx q[0];
rz(-1.3728598) q[0];
sx q[0];
rz(3.0730547) q[0];
rz(1.4463521) q[2];
sx q[2];
rz(-1.4258988) q[2];
sx q[2];
rz(3.0595487) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9431601) q[1];
sx q[1];
rz(-1.8311856) q[1];
sx q[1];
rz(-2.5501273) q[1];
rz(-pi) q[2];
rz(2.3645401) q[3];
sx q[3];
rz(-0.62824471) q[3];
sx q[3];
rz(-2.6846882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8384398) q[2];
sx q[2];
rz(-1.8586321) q[2];
sx q[2];
rz(0.34070936) q[2];
rz(0.43863145) q[3];
sx q[3];
rz(-2.3352968) q[3];
sx q[3];
rz(0.058569245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4230147) q[0];
sx q[0];
rz(-0.61841643) q[0];
sx q[0];
rz(-0.83786905) q[0];
rz(-0.92309976) q[1];
sx q[1];
rz(-1.2078614) q[1];
sx q[1];
rz(2.7941678) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23419955) q[0];
sx q[0];
rz(-3.0406018) q[0];
sx q[0];
rz(-1.1924465) q[0];
x q[1];
rz(1.3143888) q[2];
sx q[2];
rz(-1.5359582) q[2];
sx q[2];
rz(-0.23171356) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8953462) q[1];
sx q[1];
rz(-1.0965753) q[1];
sx q[1];
rz(-0.62690027) q[1];
rz(-pi) q[2];
x q[2];
rz(0.68183636) q[3];
sx q[3];
rz(-0.41636514) q[3];
sx q[3];
rz(-1.3282466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.122494) q[2];
sx q[2];
rz(-0.118003) q[2];
sx q[2];
rz(0.72354358) q[2];
rz(-1.8540234) q[3];
sx q[3];
rz(-2.5836594) q[3];
sx q[3];
rz(-3.0961228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7751854) q[0];
sx q[0];
rz(-2.7025096) q[0];
sx q[0];
rz(-1.084569) q[0];
rz(1.9547801) q[1];
sx q[1];
rz(-1.2326406) q[1];
sx q[1];
rz(1.4694227) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5195193) q[0];
sx q[0];
rz(-1.5876829) q[0];
sx q[0];
rz(-3.1114717) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7524322) q[2];
sx q[2];
rz(-1.2090328) q[2];
sx q[2];
rz(-3.0856067) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.96257229) q[1];
sx q[1];
rz(-2.031759) q[1];
sx q[1];
rz(-0.42862062) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9061435) q[3];
sx q[3];
rz(-1.9108539) q[3];
sx q[3];
rz(0.53676134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0161418) q[2];
sx q[2];
rz(-2.2405388) q[2];
sx q[2];
rz(-0.93552843) q[2];
rz(-2.9901796) q[3];
sx q[3];
rz(-0.97984034) q[3];
sx q[3];
rz(2.6305731) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22684111) q[0];
sx q[0];
rz(-1.5837357) q[0];
sx q[0];
rz(2.3874808) q[0];
rz(-0.61112815) q[1];
sx q[1];
rz(-1.1437623) q[1];
sx q[1];
rz(-0.67620826) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0026086) q[0];
sx q[0];
rz(-1.7327285) q[0];
sx q[0];
rz(1.9011628) q[0];
rz(-pi) q[1];
rz(2.7774413) q[2];
sx q[2];
rz(-2.0937348) q[2];
sx q[2];
rz(2.5464818) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3776548) q[1];
sx q[1];
rz(-1.7309524) q[1];
sx q[1];
rz(-2.9679144) q[1];
x q[2];
rz(-1.4884454) q[3];
sx q[3];
rz(-0.98656482) q[3];
sx q[3];
rz(1.0109993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.0040697441) q[2];
sx q[2];
rz(-2.4293032) q[2];
sx q[2];
rz(-0.75312692) q[2];
rz(2.5774041) q[3];
sx q[3];
rz(-2.0337532) q[3];
sx q[3];
rz(-0.75642014) q[3];
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
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18993987) q[0];
sx q[0];
rz(-0.48007444) q[0];
sx q[0];
rz(-0.22115627) q[0];
rz(3.0033374) q[1];
sx q[1];
rz(-1.6860551) q[1];
sx q[1];
rz(2.5855605) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5297896) q[0];
sx q[0];
rz(-1.3767813) q[0];
sx q[0];
rz(-0.24873269) q[0];
rz(-pi) q[1];
rz(-1.959461) q[2];
sx q[2];
rz(-1.920653) q[2];
sx q[2];
rz(1.6528724) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9784402) q[1];
sx q[1];
rz(-1.5915055) q[1];
sx q[1];
rz(0.15343517) q[1];
rz(-1.1193854) q[3];
sx q[3];
rz(-1.6777633) q[3];
sx q[3];
rz(-0.2052923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.14744645) q[2];
sx q[2];
rz(-0.7615971) q[2];
sx q[2];
rz(0.93696326) q[2];
rz(-2.7312036) q[3];
sx q[3];
rz(-0.97987163) q[3];
sx q[3];
rz(-2.4858937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
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
rz(-0.74528247) q[0];
sx q[0];
rz(-2.2672125) q[0];
sx q[0];
rz(-2.1790047) q[0];
rz(-2.4707322) q[1];
sx q[1];
rz(-0.52898359) q[1];
sx q[1];
rz(-0.76505351) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6214692) q[0];
sx q[0];
rz(-1.7304106) q[0];
sx q[0];
rz(1.9748715) q[0];
x q[1];
rz(-0.069326055) q[2];
sx q[2];
rz(-1.6717807) q[2];
sx q[2];
rz(2.947383) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9967186) q[1];
sx q[1];
rz(-1.178374) q[1];
sx q[1];
rz(-1.0189572) q[1];
rz(-pi) q[2];
rz(1.2575862) q[3];
sx q[3];
rz(-2.7566285) q[3];
sx q[3];
rz(-2.0244618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.700231) q[2];
sx q[2];
rz(-2.2141778) q[2];
sx q[2];
rz(-0.78511635) q[2];
rz(2.650812) q[3];
sx q[3];
rz(-1.1622585) q[3];
sx q[3];
rz(2.6665915) q[3];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9908776) q[0];
sx q[0];
rz(-3.0325723) q[0];
sx q[0];
rz(-2.9955067) q[0];
rz(-0.27169216) q[1];
sx q[1];
rz(-2.3482595) q[1];
sx q[1];
rz(2.9300516) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78286029) q[0];
sx q[0];
rz(-1.7365337) q[0];
sx q[0];
rz(0.77605794) q[0];
rz(-pi) q[1];
rz(1.1988099) q[2];
sx q[2];
rz(-1.6787207) q[2];
sx q[2];
rz(-1.4733009) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.817343) q[1];
sx q[1];
rz(-1.5012263) q[1];
sx q[1];
rz(-2.9051498) q[1];
rz(-pi) q[2];
x q[2];
rz(0.11084307) q[3];
sx q[3];
rz(-0.97358429) q[3];
sx q[3];
rz(-0.66652121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3286256) q[2];
sx q[2];
rz(-0.93572891) q[2];
sx q[2];
rz(-0.82175559) q[2];
rz(1.8030608) q[3];
sx q[3];
rz(-2.0526363) q[3];
sx q[3];
rz(-0.72430044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(2.5271135) q[0];
sx q[0];
rz(-2.4452657) q[0];
sx q[0];
rz(1.3592199) q[0];
rz(0.95787734) q[1];
sx q[1];
rz(-0.33708894) q[1];
sx q[1];
rz(2.8639796) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2045317) q[0];
sx q[0];
rz(-0.87420428) q[0];
sx q[0];
rz(-2.2995528) q[0];
rz(-2.8238966) q[2];
sx q[2];
rz(-0.80149273) q[2];
sx q[2];
rz(1.7048175) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6093919) q[1];
sx q[1];
rz(-1.4436296) q[1];
sx q[1];
rz(2.723477) q[1];
x q[2];
rz(-1.7429535) q[3];
sx q[3];
rz(-1.1973698) q[3];
sx q[3];
rz(2.7368054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9684888) q[2];
sx q[2];
rz(-0.59440458) q[2];
sx q[2];
rz(-1.8572726) q[2];
rz(2.3641018) q[3];
sx q[3];
rz(-2.5992664) q[3];
sx q[3];
rz(-2.8418181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1182627) q[0];
sx q[0];
rz(-1.4813923) q[0];
sx q[0];
rz(3.135664) q[0];
rz(-1.802035) q[1];
sx q[1];
rz(-2.1052994) q[1];
sx q[1];
rz(0.48066995) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.04881415) q[0];
sx q[0];
rz(-1.5546106) q[0];
sx q[0];
rz(1.5850787) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1176012) q[2];
sx q[2];
rz(-0.725774) q[2];
sx q[2];
rz(0.99601907) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8816068) q[1];
sx q[1];
rz(-2.1602958) q[1];
sx q[1];
rz(-2.0719253) q[1];
rz(-pi) q[2];
rz(-3.0044286) q[3];
sx q[3];
rz(-1.8347632) q[3];
sx q[3];
rz(0.19318737) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9504451) q[2];
sx q[2];
rz(-0.55926776) q[2];
sx q[2];
rz(2.9648103) q[2];
rz(-2.2117129) q[3];
sx q[3];
rz(-1.9527438) q[3];
sx q[3];
rz(-2.1117579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69242351) q[0];
sx q[0];
rz(-1.7208736) q[0];
sx q[0];
rz(2.0509913) q[0];
rz(-1.0724267) q[1];
sx q[1];
rz(-1.6630395) q[1];
sx q[1];
rz(1.6987775) q[1];
rz(-0.49184797) q[2];
sx q[2];
rz(-1.576423) q[2];
sx q[2];
rz(1.5058422) q[2];
rz(1.4789875) q[3];
sx q[3];
rz(-2.1740365) q[3];
sx q[3];
rz(-2.4400644) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
