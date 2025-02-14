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
rz(1.2827058) q[0];
sx q[0];
rz(1.016322) q[0];
sx q[0];
rz(10.680351) q[0];
rz(1.7805055) q[1];
sx q[1];
rz(-2.1828916) q[1];
sx q[1];
rz(-0.60671848) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3253759) q[0];
sx q[0];
rz(-2.0097549) q[0];
sx q[0];
rz(-0.75876848) q[0];
rz(-pi) q[1];
rz(-0.032564596) q[2];
sx q[2];
rz(-1.835602) q[2];
sx q[2];
rz(-0.84003583) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.42598924) q[1];
sx q[1];
rz(-1.9349602) q[1];
sx q[1];
rz(-3.0028572) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6912724) q[3];
sx q[3];
rz(-1.4109777) q[3];
sx q[3];
rz(-0.32630703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2029734) q[2];
sx q[2];
rz(-1.9846658) q[2];
sx q[2];
rz(-0.17091664) q[2];
rz(-2.0140698) q[3];
sx q[3];
rz(-2.5850962) q[3];
sx q[3];
rz(0.053430406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.663986) q[0];
sx q[0];
rz(-2.8150616) q[0];
sx q[0];
rz(0.68847454) q[0];
rz(-1.3779878) q[1];
sx q[1];
rz(-2.1207899) q[1];
sx q[1];
rz(-0.14150208) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3162075) q[0];
sx q[0];
rz(-2.0658464) q[0];
sx q[0];
rz(1.113167) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7121041) q[2];
sx q[2];
rz(-1.4802684) q[2];
sx q[2];
rz(-1.8877826) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8162839) q[1];
sx q[1];
rz(-0.46093309) q[1];
sx q[1];
rz(0.67616762) q[1];
rz(-pi) q[2];
rz(3.0612437) q[3];
sx q[3];
rz(-1.1871205) q[3];
sx q[3];
rz(1.9398882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3839533) q[2];
sx q[2];
rz(-0.44718224) q[2];
sx q[2];
rz(3.1128856) q[2];
rz(0.38255295) q[3];
sx q[3];
rz(-1.5111204) q[3];
sx q[3];
rz(0.20146519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89556995) q[0];
sx q[0];
rz(-1.0423132) q[0];
sx q[0];
rz(-0.37164715) q[0];
rz(2.9314575) q[1];
sx q[1];
rz(-0.34966436) q[1];
sx q[1];
rz(-1.9336611) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9555617) q[0];
sx q[0];
rz(-0.81350858) q[0];
sx q[0];
rz(-0.71788089) q[0];
rz(-pi) q[1];
x q[1];
rz(0.082449989) q[2];
sx q[2];
rz(-1.5918613) q[2];
sx q[2];
rz(2.9631071) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4285046) q[1];
sx q[1];
rz(-2.2853973) q[1];
sx q[1];
rz(2.3355161) q[1];
x q[2];
rz(-2.7961416) q[3];
sx q[3];
rz(-1.6092704) q[3];
sx q[3];
rz(-1.0055804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5189884) q[2];
sx q[2];
rz(-3.0486139) q[2];
sx q[2];
rz(0.63817111) q[2];
rz(-2.7291164) q[3];
sx q[3];
rz(-1.6715489) q[3];
sx q[3];
rz(-1.1023869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.07688044) q[0];
sx q[0];
rz(-0.91016722) q[0];
sx q[0];
rz(1.3077211) q[0];
rz(2.4619596) q[1];
sx q[1];
rz(-0.72703397) q[1];
sx q[1];
rz(1.1839428) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7789202) q[0];
sx q[0];
rz(-2.2244172) q[0];
sx q[0];
rz(-1.158276) q[0];
rz(-pi) q[1];
rz(-1.6865191) q[2];
sx q[2];
rz(-1.8050131) q[2];
sx q[2];
rz(2.9701592) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1278504) q[1];
sx q[1];
rz(-0.35687414) q[1];
sx q[1];
rz(-2.9452717) q[1];
rz(-2.3700962) q[3];
sx q[3];
rz(-0.1813387) q[3];
sx q[3];
rz(-0.69252959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.1188941) q[2];
sx q[2];
rz(-2.9939632) q[2];
sx q[2];
rz(-1.0559319) q[2];
rz(2.8583156) q[3];
sx q[3];
rz(-1.8746459) q[3];
sx q[3];
rz(1.383708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088575514) q[0];
sx q[0];
rz(-2.179189) q[0];
sx q[0];
rz(-1.4085294) q[0];
rz(0.22964302) q[1];
sx q[1];
rz(-0.85669986) q[1];
sx q[1];
rz(-1.9502669) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8386779) q[0];
sx q[0];
rz(-1.9698799) q[0];
sx q[0];
rz(-2.9877034) q[0];
rz(-pi) q[1];
rz(1.2651839) q[2];
sx q[2];
rz(-1.4508045) q[2];
sx q[2];
rz(-0.12745276) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.31684418) q[1];
sx q[1];
rz(-1.1543546) q[1];
sx q[1];
rz(-1.29711) q[1];
x q[2];
rz(2.2681885) q[3];
sx q[3];
rz(-2.3583989) q[3];
sx q[3];
rz(2.7627373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4952937) q[2];
sx q[2];
rz(-2.7806492) q[2];
sx q[2];
rz(-1.6634644) q[2];
rz(2.9804001) q[3];
sx q[3];
rz(-1.5321923) q[3];
sx q[3];
rz(-2.351779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6394871) q[0];
sx q[0];
rz(-2.7171071) q[0];
sx q[0];
rz(2.2232527) q[0];
rz(0.25686747) q[1];
sx q[1];
rz(-2.145642) q[1];
sx q[1];
rz(-1.1161944) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6300121) q[0];
sx q[0];
rz(-2.0697547) q[0];
sx q[0];
rz(2.5415238) q[0];
x q[1];
rz(0.29735469) q[2];
sx q[2];
rz(-0.28536826) q[2];
sx q[2];
rz(0.076208027) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5600561) q[1];
sx q[1];
rz(-1.3764842) q[1];
sx q[1];
rz(2.8981853) q[1];
x q[2];
rz(0.61577101) q[3];
sx q[3];
rz(-0.29524657) q[3];
sx q[3];
rz(-0.24850445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.18672289) q[2];
sx q[2];
rz(-0.8064417) q[2];
sx q[2];
rz(0.40697971) q[2];
rz(-2.5278029) q[3];
sx q[3];
rz(-1.611462) q[3];
sx q[3];
rz(-0.49155244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75519049) q[0];
sx q[0];
rz(-0.84119216) q[0];
sx q[0];
rz(2.2156773) q[0];
rz(2.6689802) q[1];
sx q[1];
rz(-2.0874529) q[1];
sx q[1];
rz(2.6714163) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42379728) q[0];
sx q[0];
rz(-1.6542302) q[0];
sx q[0];
rz(0.83441894) q[0];
rz(-pi) q[1];
rz(1.0711925) q[2];
sx q[2];
rz(-2.4322369) q[2];
sx q[2];
rz(1.5200523) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2963841) q[1];
sx q[1];
rz(-1.607158) q[1];
sx q[1];
rz(-0.78486376) q[1];
rz(-1.5298793) q[3];
sx q[3];
rz(-2.4114263) q[3];
sx q[3];
rz(2.1388291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.65972313) q[2];
sx q[2];
rz(-1.0708151) q[2];
sx q[2];
rz(-2.3107963) q[2];
rz(3.0610436) q[3];
sx q[3];
rz(-2.060067) q[3];
sx q[3];
rz(-2.1848333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-1.047121) q[0];
sx q[0];
rz(-1.5249277) q[0];
sx q[0];
rz(2.9834874) q[0];
rz(-0.72599167) q[1];
sx q[1];
rz(-2.7114365) q[1];
sx q[1];
rz(0.37011883) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86240101) q[0];
sx q[0];
rz(-1.4790311) q[0];
sx q[0];
rz(2.807995) q[0];
rz(-pi) q[1];
rz(0.49741715) q[2];
sx q[2];
rz(-2.3193365) q[2];
sx q[2];
rz(-2.9679839) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.98277011) q[1];
sx q[1];
rz(-1.1825605) q[1];
sx q[1];
rz(-2.7169711) q[1];
rz(-0.40970476) q[3];
sx q[3];
rz(-1.0632044) q[3];
sx q[3];
rz(2.0970124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.98081723) q[2];
sx q[2];
rz(-1.2951916) q[2];
sx q[2];
rz(-3.004461) q[2];
rz(1.49617) q[3];
sx q[3];
rz(-0.6936332) q[3];
sx q[3];
rz(-2.6579198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54336035) q[0];
sx q[0];
rz(-1.0022663) q[0];
sx q[0];
rz(0.12635669) q[0];
rz(2.5670746) q[1];
sx q[1];
rz(-0.9205598) q[1];
sx q[1];
rz(-2.394302) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3355646) q[0];
sx q[0];
rz(-1.757041) q[0];
sx q[0];
rz(-1.9646909) q[0];
x q[1];
rz(-2.4392082) q[2];
sx q[2];
rz(-2.0364072) q[2];
sx q[2];
rz(1.7543751) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.14366985) q[1];
sx q[1];
rz(-1.9410681) q[1];
sx q[1];
rz(-2.0475494) q[1];
rz(-pi) q[2];
rz(-2.6797557) q[3];
sx q[3];
rz(-0.29134068) q[3];
sx q[3];
rz(-0.1801462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0390465) q[2];
sx q[2];
rz(-1.1209844) q[2];
sx q[2];
rz(-2.5940564) q[2];
rz(-1.840379) q[3];
sx q[3];
rz(-2.0042714) q[3];
sx q[3];
rz(3.014452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.309677) q[0];
sx q[0];
rz(-1.7597821) q[0];
sx q[0];
rz(-0.58050138) q[0];
rz(-0.25513395) q[1];
sx q[1];
rz(-2.3105123) q[1];
sx q[1];
rz(1.1855804) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0609293) q[0];
sx q[0];
rz(-2.1619045) q[0];
sx q[0];
rz(0.76158686) q[0];
x q[1];
rz(1.2078778) q[2];
sx q[2];
rz(-1.5722154) q[2];
sx q[2];
rz(3.0848173) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.285814) q[1];
sx q[1];
rz(-1.6596982) q[1];
sx q[1];
rz(-2.6707021) q[1];
rz(-0.85082741) q[3];
sx q[3];
rz(-1.730207) q[3];
sx q[3];
rz(-3.0992143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.4127976) q[2];
sx q[2];
rz(-2.2788861) q[2];
sx q[2];
rz(-1.0069138) q[2];
rz(-1.6179196) q[3];
sx q[3];
rz(-2.1606052) q[3];
sx q[3];
rz(-1.9699008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11378743) q[0];
sx q[0];
rz(-1.8814977) q[0];
sx q[0];
rz(0.30497288) q[0];
rz(-1.8995359) q[1];
sx q[1];
rz(-1.2868953) q[1];
sx q[1];
rz(-2.3931265) q[1];
rz(-1.9340877) q[2];
sx q[2];
rz(-0.45798326) q[2];
sx q[2];
rz(1.9275673) q[2];
rz(0.92963582) q[3];
sx q[3];
rz(-2.2916678) q[3];
sx q[3];
rz(-2.5558932) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
