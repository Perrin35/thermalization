OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.92575443) q[0];
sx q[0];
rz(-0.32045066) q[0];
sx q[0];
rz(-0.12838636) q[0];
rz(-2.1865891) q[1];
sx q[1];
rz(-0.65989143) q[1];
sx q[1];
rz(-0.58712062) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40872657) q[0];
sx q[0];
rz(-0.48983296) q[0];
sx q[0];
rz(-0.31087713) q[0];
rz(-2.0287747) q[2];
sx q[2];
rz(-1.0931778) q[2];
sx q[2];
rz(2.7736563) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94226664) q[1];
sx q[1];
rz(-1.1360511) q[1];
sx q[1];
rz(2.4995657) q[1];
rz(-pi) q[2];
rz(-1.5641603) q[3];
sx q[3];
rz(-0.38449461) q[3];
sx q[3];
rz(-0.16669434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8762274) q[2];
sx q[2];
rz(-0.23398016) q[2];
sx q[2];
rz(2.2205676) q[2];
rz(-1.5246576) q[3];
sx q[3];
rz(-1.2911258) q[3];
sx q[3];
rz(0.18572148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8601473) q[0];
sx q[0];
rz(-2.7466819) q[0];
sx q[0];
rz(1.7922147) q[0];
rz(-2.7941864) q[1];
sx q[1];
rz(-1.9527304) q[1];
sx q[1];
rz(1.4030392) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31579298) q[0];
sx q[0];
rz(-2.1496338) q[0];
sx q[0];
rz(-1.3114503) q[0];
rz(-pi) q[1];
rz(-2.5591056) q[2];
sx q[2];
rz(-2.2299924) q[2];
sx q[2];
rz(0.90345736) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.35866061) q[1];
sx q[1];
rz(-2.4110849) q[1];
sx q[1];
rz(-0.23083861) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8946793) q[3];
sx q[3];
rz(-0.71454701) q[3];
sx q[3];
rz(-1.1199878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8189524) q[2];
sx q[2];
rz(-1.3078657) q[2];
sx q[2];
rz(2.9827706) q[2];
rz(-0.21765503) q[3];
sx q[3];
rz(-1.7526151) q[3];
sx q[3];
rz(1.3505107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7749216) q[0];
sx q[0];
rz(-1.8407624) q[0];
sx q[0];
rz(-1.8721254) q[0];
rz(0.41995755) q[1];
sx q[1];
rz(-2.1407514) q[1];
sx q[1];
rz(-1.5392083) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8458662) q[0];
sx q[0];
rz(-0.55385607) q[0];
sx q[0];
rz(-0.94828301) q[0];
rz(1.3871269) q[2];
sx q[2];
rz(-1.78698) q[2];
sx q[2];
rz(0.64858299) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6091385) q[1];
sx q[1];
rz(-1.6241606) q[1];
sx q[1];
rz(1.1751925) q[1];
rz(-pi) q[2];
rz(1.9004542) q[3];
sx q[3];
rz(-1.0287549) q[3];
sx q[3];
rz(-0.84621669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0594844) q[2];
sx q[2];
rz(-0.84222811) q[2];
sx q[2];
rz(-2.9498937) q[2];
rz(-2.5890403) q[3];
sx q[3];
rz(-1.6165761) q[3];
sx q[3];
rz(1.0244757) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4100274) q[0];
sx q[0];
rz(-2.4931694) q[0];
sx q[0];
rz(-0.60436526) q[0];
rz(-1.8961689) q[1];
sx q[1];
rz(-0.75029293) q[1];
sx q[1];
rz(-0.85974685) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.788279) q[0];
sx q[0];
rz(-1.9082359) q[0];
sx q[0];
rz(0.91633015) q[0];
rz(-pi) q[1];
rz(-3.1273975) q[2];
sx q[2];
rz(-1.4687339) q[2];
sx q[2];
rz(1.0400187) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5936232) q[1];
sx q[1];
rz(-1.3892349) q[1];
sx q[1];
rz(-1.5201925) q[1];
rz(-pi) q[2];
rz(1.3172011) q[3];
sx q[3];
rz(-0.86652866) q[3];
sx q[3];
rz(0.21237838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2802281) q[2];
sx q[2];
rz(-1.0140398) q[2];
sx q[2];
rz(1.4385983) q[2];
rz(1.2922618) q[3];
sx q[3];
rz(-0.82787138) q[3];
sx q[3];
rz(-1.0546225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.0786521) q[0];
sx q[0];
rz(-2.8684454) q[0];
sx q[0];
rz(-2.8237421) q[0];
rz(1.3392797) q[1];
sx q[1];
rz(-0.41792089) q[1];
sx q[1];
rz(0.97871614) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1237549) q[0];
sx q[0];
rz(-2.7394501) q[0];
sx q[0];
rz(-0.077362424) q[0];
x q[1];
rz(0.42338223) q[2];
sx q[2];
rz(-2.9278594) q[2];
sx q[2];
rz(-0.20704421) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5086245) q[1];
sx q[1];
rz(-0.76286722) q[1];
sx q[1];
rz(2.9588225) q[1];
rz(2.3581402) q[3];
sx q[3];
rz(-0.26973596) q[3];
sx q[3];
rz(1.3470284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.19834441) q[2];
sx q[2];
rz(-2.1965616) q[2];
sx q[2];
rz(-1.1172969) q[2];
rz(-2.9605588) q[3];
sx q[3];
rz(-1.2369316) q[3];
sx q[3];
rz(1.9443289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70500526) q[0];
sx q[0];
rz(-2.8745108) q[0];
sx q[0];
rz(-1.8319112) q[0];
rz(1.0218703) q[1];
sx q[1];
rz(-2.3049057) q[1];
sx q[1];
rz(-1.6530564) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2383136) q[0];
sx q[0];
rz(-2.2198703) q[0];
sx q[0];
rz(-2.7191976) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4632881) q[2];
sx q[2];
rz(-1.3875752) q[2];
sx q[2];
rz(1.9187294) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.19410876) q[1];
sx q[1];
rz(-2.6257537) q[1];
sx q[1];
rz(-0.2920002) q[1];
rz(-1.7506787) q[3];
sx q[3];
rz(-2.5733786) q[3];
sx q[3];
rz(3.1364322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66270193) q[2];
sx q[2];
rz(-1.7544489) q[2];
sx q[2];
rz(-0.26408163) q[2];
rz(1.4763907) q[3];
sx q[3];
rz(-2.4614406) q[3];
sx q[3];
rz(-2.7217854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5298117) q[0];
sx q[0];
rz(-1.6407069) q[0];
sx q[0];
rz(1.06426) q[0];
rz(-3.0423959) q[1];
sx q[1];
rz(-1.3490889) q[1];
sx q[1];
rz(1.4605716) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6457466) q[0];
sx q[0];
rz(-0.48862132) q[0];
sx q[0];
rz(-3.0054128) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6978092) q[2];
sx q[2];
rz(-2.840503) q[2];
sx q[2];
rz(-1.3236698) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8122892) q[1];
sx q[1];
rz(-1.9871622) q[1];
sx q[1];
rz(1.8805481) q[1];
rz(2.0123294) q[3];
sx q[3];
rz(-2.0534424) q[3];
sx q[3];
rz(-0.03015524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3508241) q[2];
sx q[2];
rz(-2.0227573) q[2];
sx q[2];
rz(2.5420945) q[2];
rz(-2.2439469) q[3];
sx q[3];
rz(-1.41956) q[3];
sx q[3];
rz(3.1035778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.81383234) q[0];
sx q[0];
rz(-1.0603511) q[0];
sx q[0];
rz(-1.4014442) q[0];
rz(0.071488149) q[1];
sx q[1];
rz(-1.6782327) q[1];
sx q[1];
rz(1.00114) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5073231) q[0];
sx q[0];
rz(-1.8232913) q[0];
sx q[0];
rz(2.6802313) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4779448) q[2];
sx q[2];
rz(-1.8816684) q[2];
sx q[2];
rz(0.03229095) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.82187958) q[1];
sx q[1];
rz(-0.94187573) q[1];
sx q[1];
rz(-1.4075539) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4720932) q[3];
sx q[3];
rz(-2.0560798) q[3];
sx q[3];
rz(2.8593331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1772168) q[2];
sx q[2];
rz(-1.5528677) q[2];
sx q[2];
rz(-0.71844086) q[2];
rz(-0.046772379) q[3];
sx q[3];
rz(-2.434157) q[3];
sx q[3];
rz(2.5717946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8319594) q[0];
sx q[0];
rz(-2.9871812) q[0];
sx q[0];
rz(-0.69044789) q[0];
rz(-0.18028232) q[1];
sx q[1];
rz(-1.0612265) q[1];
sx q[1];
rz(2.8188425) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0092516) q[0];
sx q[0];
rz(-1.4173495) q[0];
sx q[0];
rz(-1.4645534) q[0];
x q[1];
rz(2.1229136) q[2];
sx q[2];
rz(-1.5763057) q[2];
sx q[2];
rz(-1.6679279) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1806644) q[1];
sx q[1];
rz(-0.81392787) q[1];
sx q[1];
rz(0.77427534) q[1];
x q[2];
rz(2.2261103) q[3];
sx q[3];
rz(-0.99952664) q[3];
sx q[3];
rz(2.9397428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47101578) q[2];
sx q[2];
rz(-1.3479503) q[2];
sx q[2];
rz(-1.0178817) q[2];
rz(0.0035303591) q[3];
sx q[3];
rz(-0.96960932) q[3];
sx q[3];
rz(-1.2118916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8135391) q[0];
sx q[0];
rz(-0.73754755) q[0];
sx q[0];
rz(-2.2156583) q[0];
rz(0.13735859) q[1];
sx q[1];
rz(-0.67247144) q[1];
sx q[1];
rz(0.59991178) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1178499) q[0];
sx q[0];
rz(-1.5735007) q[0];
sx q[0];
rz(0.54393025) q[0];
rz(-pi) q[1];
rz(2.4358168) q[2];
sx q[2];
rz(-2.3786491) q[2];
sx q[2];
rz(-2.0313789) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.39420262) q[1];
sx q[1];
rz(-2.5521186) q[1];
sx q[1];
rz(1.516045) q[1];
rz(-pi) q[2];
rz(2.0447015) q[3];
sx q[3];
rz(-1.2792493) q[3];
sx q[3];
rz(2.4959559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4124734) q[2];
sx q[2];
rz(-2.1226661) q[2];
sx q[2];
rz(2.4998383) q[2];
rz(1.9852091) q[3];
sx q[3];
rz(-2.3803847) q[3];
sx q[3];
rz(1.523264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625576) q[0];
sx q[0];
rz(-0.87974822) q[0];
sx q[0];
rz(0.77597822) q[0];
rz(-1.4398126) q[1];
sx q[1];
rz(-1.6701313) q[1];
sx q[1];
rz(-3.0727542) q[1];
rz(-0.65563249) q[2];
sx q[2];
rz(-2.8218695) q[2];
sx q[2];
rz(-0.8994076) q[2];
rz(-0.84550459) q[3];
sx q[3];
rz(-0.68544023) q[3];
sx q[3];
rz(-1.1123085) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
