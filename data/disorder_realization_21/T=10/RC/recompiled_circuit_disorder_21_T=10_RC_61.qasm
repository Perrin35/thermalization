OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.055846) q[0];
sx q[0];
rz(-3.0598109) q[0];
sx q[0];
rz(2.6401289) q[0];
rz(1.4986829) q[1];
sx q[1];
rz(3.5377503) q[1];
sx q[1];
rz(9.1023394) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6591588) q[0];
sx q[0];
rz(-1.6560439) q[0];
sx q[0];
rz(1.2486588) q[0];
rz(-pi) q[1];
rz(-0.85876043) q[2];
sx q[2];
rz(-0.81958629) q[2];
sx q[2];
rz(0.31529266) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0281801) q[1];
sx q[1];
rz(-0.98787381) q[1];
sx q[1];
rz(-0.41863538) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93482165) q[3];
sx q[3];
rz(-1.9536195) q[3];
sx q[3];
rz(2.8179907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.50513187) q[2];
sx q[2];
rz(-2.5487066) q[2];
sx q[2];
rz(-0.55603975) q[2];
rz(-0.83267823) q[3];
sx q[3];
rz(-1.6502389) q[3];
sx q[3];
rz(-0.94579831) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44822025) q[0];
sx q[0];
rz(-1.4602666) q[0];
sx q[0];
rz(2.9843176) q[0];
rz(-0.26113025) q[1];
sx q[1];
rz(-1.7938679) q[1];
sx q[1];
rz(0.10903407) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029023829) q[0];
sx q[0];
rz(-1.7984263) q[0];
sx q[0];
rz(2.9665222) q[0];
rz(-pi) q[1];
rz(2.878506) q[2];
sx q[2];
rz(-1.3795985) q[2];
sx q[2];
rz(-2.4441602) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1034531) q[1];
sx q[1];
rz(-0.98587576) q[1];
sx q[1];
rz(2.1409047) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29626366) q[3];
sx q[3];
rz(-0.54013541) q[3];
sx q[3];
rz(1.6903282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2382425) q[2];
sx q[2];
rz(-1.976333) q[2];
sx q[2];
rz(1.8781352) q[2];
rz(-2.8144828) q[3];
sx q[3];
rz(-1.5771022) q[3];
sx q[3];
rz(-1.9272778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-0.064421244) q[0];
sx q[0];
rz(-3.0922958) q[0];
sx q[0];
rz(1.7984614) q[0];
rz(0.24761565) q[1];
sx q[1];
rz(-0.74664465) q[1];
sx q[1];
rz(-0.48167357) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1789078) q[0];
sx q[0];
rz(-0.6324397) q[0];
sx q[0];
rz(-0.90895598) q[0];
rz(-pi) q[1];
rz(2.0273676) q[2];
sx q[2];
rz(-2.4764428) q[2];
sx q[2];
rz(-2.3014724) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8183221) q[1];
sx q[1];
rz(-2.1774877) q[1];
sx q[1];
rz(-2.2491749) q[1];
x q[2];
rz(0.62044575) q[3];
sx q[3];
rz(-0.9471604) q[3];
sx q[3];
rz(0.40943957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8032916) q[2];
sx q[2];
rz(-0.81739601) q[2];
sx q[2];
rz(2.6417007) q[2];
rz(-0.56097427) q[3];
sx q[3];
rz(-1.8818972) q[3];
sx q[3];
rz(1.4612173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19701476) q[0];
sx q[0];
rz(-0.16600969) q[0];
sx q[0];
rz(2.5894077) q[0];
rz(-1.588297) q[1];
sx q[1];
rz(-0.89893666) q[1];
sx q[1];
rz(1.2447371) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3361928) q[0];
sx q[0];
rz(-1.5691783) q[0];
sx q[0];
rz(2.8021115) q[0];
rz(-1.7719901) q[2];
sx q[2];
rz(-1.6987213) q[2];
sx q[2];
rz(2.163137) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.0075477608) q[1];
sx q[1];
rz(-1.3589077) q[1];
sx q[1];
rz(-2.120554) q[1];
rz(1.3454868) q[3];
sx q[3];
rz(-2.1300737) q[3];
sx q[3];
rz(-2.6252281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.84919471) q[2];
sx q[2];
rz(-1.2595824) q[2];
sx q[2];
rz(-1.1506895) q[2];
rz(1.4771279) q[3];
sx q[3];
rz(-1.632558) q[3];
sx q[3];
rz(2.6586444) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0634336) q[0];
sx q[0];
rz(-0.76197356) q[0];
sx q[0];
rz(3.0601236) q[0];
rz(-3.07913) q[1];
sx q[1];
rz(-1.1413347) q[1];
sx q[1];
rz(1.6385471) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3969288) q[0];
sx q[0];
rz(-1.443112) q[0];
sx q[0];
rz(0.98608195) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4625174) q[2];
sx q[2];
rz(-1.6533972) q[2];
sx q[2];
rz(-2.8464707) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0283708) q[1];
sx q[1];
rz(-1.2591259) q[1];
sx q[1];
rz(3.034364) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2522069) q[3];
sx q[3];
rz(-2.0584403) q[3];
sx q[3];
rz(0.82061758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.5082671) q[2];
sx q[2];
rz(-0.94255629) q[2];
sx q[2];
rz(-1.903669) q[2];
rz(-1.1226908) q[3];
sx q[3];
rz(-2.4653547) q[3];
sx q[3];
rz(-2.6200263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6102819) q[0];
sx q[0];
rz(-1.0045445) q[0];
sx q[0];
rz(-2.8748728) q[0];
rz(2.5807014) q[1];
sx q[1];
rz(-1.8436878) q[1];
sx q[1];
rz(2.3430603) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079433867) q[0];
sx q[0];
rz(-1.2957797) q[0];
sx q[0];
rz(-2.9955203) q[0];
x q[1];
rz(-1.919826) q[2];
sx q[2];
rz(-1.0997699) q[2];
sx q[2];
rz(2.7407321) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.16520271) q[1];
sx q[1];
rz(-0.573728) q[1];
sx q[1];
rz(1.4742673) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.9540326) q[3];
sx q[3];
rz(-0.87152374) q[3];
sx q[3];
rz(1.1773393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.16053998) q[2];
sx q[2];
rz(-1.2489698) q[2];
sx q[2];
rz(-0.36671656) q[2];
rz(1.2612873) q[3];
sx q[3];
rz(-1.4533318) q[3];
sx q[3];
rz(-0.096207531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7325608) q[0];
sx q[0];
rz(-2.2213187) q[0];
sx q[0];
rz(-0.60638705) q[0];
rz(2.9442893) q[1];
sx q[1];
rz(-2.0154672) q[1];
sx q[1];
rz(-0.46404776) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7921917) q[0];
sx q[0];
rz(-2.1878625) q[0];
sx q[0];
rz(-2.6130136) q[0];
x q[1];
rz(-1.3299994) q[2];
sx q[2];
rz(-1.5511998) q[2];
sx q[2];
rz(2.7407676) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.53941599) q[1];
sx q[1];
rz(-2.5578941) q[1];
sx q[1];
rz(-1.0023487) q[1];
rz(-pi) q[2];
rz(-0.23630996) q[3];
sx q[3];
rz(-1.4910306) q[3];
sx q[3];
rz(-0.41111708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.93418926) q[2];
sx q[2];
rz(-2.1384017) q[2];
sx q[2];
rz(2.8835473) q[2];
rz(1.9559654) q[3];
sx q[3];
rz(-1.5153171) q[3];
sx q[3];
rz(-3.1380838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.946452) q[0];
sx q[0];
rz(-1.8608681) q[0];
sx q[0];
rz(-0.38129693) q[0];
rz(-0.095245846) q[1];
sx q[1];
rz(-0.97243273) q[1];
sx q[1];
rz(-1.4415178) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38948108) q[0];
sx q[0];
rz(-1.507326) q[0];
sx q[0];
rz(-1.585929) q[0];
x q[1];
rz(-1.5324138) q[2];
sx q[2];
rz(-0.95853171) q[2];
sx q[2];
rz(-1.6910764) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0563593) q[1];
sx q[1];
rz(-1.4409522) q[1];
sx q[1];
rz(1.4443881) q[1];
rz(-pi) q[2];
rz(1.5278682) q[3];
sx q[3];
rz(-1.9044442) q[3];
sx q[3];
rz(-1.4294525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1429446) q[2];
sx q[2];
rz(-0.412985) q[2];
sx q[2];
rz(2.9150035) q[2];
rz(-0.44858027) q[3];
sx q[3];
rz(-1.5357163) q[3];
sx q[3];
rz(0.8297689) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3806234) q[0];
sx q[0];
rz(-2.3801104) q[0];
sx q[0];
rz(1.3990078) q[0];
rz(-0.31708583) q[1];
sx q[1];
rz(-1.6665019) q[1];
sx q[1];
rz(-0.98659602) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50842677) q[0];
sx q[0];
rz(-2.0619443) q[0];
sx q[0];
rz(-1.3915865) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1836906) q[2];
sx q[2];
rz(-2.3447678) q[2];
sx q[2];
rz(0.56607841) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1074859) q[1];
sx q[1];
rz(-1.7120546) q[1];
sx q[1];
rz(-0.60955255) q[1];
rz(-pi) q[2];
rz(1.1986198) q[3];
sx q[3];
rz(-0.44789207) q[3];
sx q[3];
rz(-0.4263634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2150779) q[2];
sx q[2];
rz(-2.4145917) q[2];
sx q[2];
rz(-2.731936) q[2];
rz(-2.8783197) q[3];
sx q[3];
rz(-1.8299088) q[3];
sx q[3];
rz(0.51945654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39919329) q[0];
sx q[0];
rz(-3.0629459) q[0];
sx q[0];
rz(-1.7364527) q[0];
rz(-0.82110226) q[1];
sx q[1];
rz(-2.2228873) q[1];
sx q[1];
rz(1.7260889) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4840568) q[0];
sx q[0];
rz(-2.7298379) q[0];
sx q[0];
rz(-2.5176237) q[0];
rz(-pi) q[1];
x q[1];
rz(2.892898) q[2];
sx q[2];
rz(-1.9242312) q[2];
sx q[2];
rz(-0.24662185) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6092097) q[1];
sx q[1];
rz(-1.4192974) q[1];
sx q[1];
rz(-1.3551559) q[1];
rz(1.8652925) q[3];
sx q[3];
rz(-2.6101972) q[3];
sx q[3];
rz(0.48081276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4225509) q[2];
sx q[2];
rz(-2.8406403) q[2];
sx q[2];
rz(-0.12410513) q[2];
rz(2.1758046) q[3];
sx q[3];
rz(-1.6671168) q[3];
sx q[3];
rz(-0.66108274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.538095) q[0];
sx q[0];
rz(-0.24833831) q[0];
sx q[0];
rz(-0.86059358) q[0];
rz(-2.8339236) q[1];
sx q[1];
rz(-1.888702) q[1];
sx q[1];
rz(-1.9370334) q[1];
rz(-2.5406191) q[2];
sx q[2];
rz(-0.31243159) q[2];
sx q[2];
rz(2.5675788) q[2];
rz(0.55660558) q[3];
sx q[3];
rz(-1.6986596) q[3];
sx q[3];
rz(-1.2996775) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];