OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.38874415) q[0];
sx q[0];
rz(3.677877) q[0];
sx q[0];
rz(10.372547) q[0];
rz(-1.3287969) q[1];
sx q[1];
rz(4.4089945) q[1];
sx q[1];
rz(10.452527) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9337024) q[0];
sx q[0];
rz(-2.7650802) q[0];
sx q[0];
rz(0.062113751) q[0];
x q[1];
rz(2.0601294) q[2];
sx q[2];
rz(-1.7667734) q[2];
sx q[2];
rz(-2.8010362) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6537522) q[1];
sx q[1];
rz(-2.2334705) q[1];
sx q[1];
rz(-1.0767879) q[1];
rz(-pi) q[2];
x q[2];
rz(0.62698934) q[3];
sx q[3];
rz(-1.9416182) q[3];
sx q[3];
rz(0.8644608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0119005) q[2];
sx q[2];
rz(-1.4346069) q[2];
sx q[2];
rz(2.091308) q[2];
rz(-2.0283279) q[3];
sx q[3];
rz(-0.89171019) q[3];
sx q[3];
rz(-3.0734857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-3.0691836) q[0];
sx q[0];
rz(-1.8962815) q[0];
sx q[0];
rz(-0.29775277) q[0];
rz(2.521926) q[1];
sx q[1];
rz(-2.1344118) q[1];
sx q[1];
rz(1.108095) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71582687) q[0];
sx q[0];
rz(-2.1799488) q[0];
sx q[0];
rz(-2.8858375) q[0];
x q[1];
rz(-0.98845311) q[2];
sx q[2];
rz(-2.4465912) q[2];
sx q[2];
rz(0.37441355) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.121671) q[1];
sx q[1];
rz(-0.53598511) q[1];
sx q[1];
rz(-2.3977445) q[1];
rz(-pi) q[2];
rz(-2.439019) q[3];
sx q[3];
rz(-1.7785903) q[3];
sx q[3];
rz(-2.7489565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6790598) q[2];
sx q[2];
rz(-2.1647537) q[2];
sx q[2];
rz(-2.1662946) q[2];
rz(-2.2235928) q[3];
sx q[3];
rz(-1.8564329) q[3];
sx q[3];
rz(-0.29618922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.179203) q[0];
sx q[0];
rz(-2.2900892) q[0];
sx q[0];
rz(0.54291022) q[0];
rz(2.2593598) q[1];
sx q[1];
rz(-1.135332) q[1];
sx q[1];
rz(-0.96484819) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0438) q[0];
sx q[0];
rz(-1.8505197) q[0];
sx q[0];
rz(-2.8066737) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29693551) q[2];
sx q[2];
rz(-1.6116217) q[2];
sx q[2];
rz(2.9253935) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0286897) q[1];
sx q[1];
rz(-1.5064872) q[1];
sx q[1];
rz(1.6921922) q[1];
x q[2];
rz(-1.987625) q[3];
sx q[3];
rz(-1.3891451) q[3];
sx q[3];
rz(-2.3895398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0470011) q[2];
sx q[2];
rz(-0.61085218) q[2];
sx q[2];
rz(1.1331406) q[2];
rz(-0.23162332) q[3];
sx q[3];
rz(-1.8685721) q[3];
sx q[3];
rz(2.384322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8621181) q[0];
sx q[0];
rz(-0.010443895) q[0];
sx q[0];
rz(1.3765155) q[0];
rz(2.6230985) q[1];
sx q[1];
rz(-1.8771749) q[1];
sx q[1];
rz(-0.24212295) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69232363) q[0];
sx q[0];
rz(-1.4081435) q[0];
sx q[0];
rz(-0.67740324) q[0];
x q[1];
rz(2.7119615) q[2];
sx q[2];
rz(-1.5826844) q[2];
sx q[2];
rz(0.55693835) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.081269216) q[1];
sx q[1];
rz(-2.6594866) q[1];
sx q[1];
rz(0.37270765) q[1];
x q[2];
rz(-0.30947134) q[3];
sx q[3];
rz(-2.7052393) q[3];
sx q[3];
rz(-2.5040124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1233998) q[2];
sx q[2];
rz(-0.96863666) q[2];
sx q[2];
rz(-0.094853178) q[2];
rz(-1.3421966) q[3];
sx q[3];
rz(-1.3972524) q[3];
sx q[3];
rz(-2.7024787) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3577394) q[0];
sx q[0];
rz(-1.8308324) q[0];
sx q[0];
rz(0.36079303) q[0];
rz(-1.3882673) q[1];
sx q[1];
rz(-1.3307945) q[1];
sx q[1];
rz(-2.0070019) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73571262) q[0];
sx q[0];
rz(-2.1942733) q[0];
sx q[0];
rz(-2.1924125) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0005433) q[2];
sx q[2];
rz(-1.2417214) q[2];
sx q[2];
rz(0.82440257) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7736588) q[1];
sx q[1];
rz(-1.415167) q[1];
sx q[1];
rz(2.6917798) q[1];
x q[2];
rz(2.9633425) q[3];
sx q[3];
rz(-2.6986487) q[3];
sx q[3];
rz(0.49803842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0052884) q[2];
sx q[2];
rz(-1.7207928) q[2];
sx q[2];
rz(0.57265442) q[2];
rz(-2.2128361) q[3];
sx q[3];
rz(-2.6199665) q[3];
sx q[3];
rz(-1.1675534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3271493) q[0];
sx q[0];
rz(-1.0961908) q[0];
sx q[0];
rz(-1.2493398) q[0];
rz(-1.918474) q[1];
sx q[1];
rz(-1.616281) q[1];
sx q[1];
rz(1.9893601) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.457068) q[0];
sx q[0];
rz(-3.0355434) q[0];
sx q[0];
rz(1.690879) q[0];
rz(0.67201519) q[2];
sx q[2];
rz(-2.5294371) q[2];
sx q[2];
rz(1.9415346) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2684979) q[1];
sx q[1];
rz(-0.22646204) q[1];
sx q[1];
rz(-0.48392673) q[1];
x q[2];
rz(0.67752083) q[3];
sx q[3];
rz(-1.0761677) q[3];
sx q[3];
rz(-0.90466162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.46889177) q[2];
sx q[2];
rz(-1.9273309) q[2];
sx q[2];
rz(-2.1506298) q[2];
rz(-2.4937566) q[3];
sx q[3];
rz(-2.1760553) q[3];
sx q[3];
rz(2.794054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25794849) q[0];
sx q[0];
rz(-2.914496) q[0];
sx q[0];
rz(-0.062285475) q[0];
rz(0.1858055) q[1];
sx q[1];
rz(-1.4567016) q[1];
sx q[1];
rz(0.3947765) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1226574) q[0];
sx q[0];
rz(-1.1766953) q[0];
sx q[0];
rz(-1.8271441) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4418418) q[2];
sx q[2];
rz(-2.664898) q[2];
sx q[2];
rz(1.9112019) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.77183206) q[1];
sx q[1];
rz(-2.6409915) q[1];
sx q[1];
rz(-2.6066149) q[1];
x q[2];
rz(2.9069101) q[3];
sx q[3];
rz(-0.15771401) q[3];
sx q[3];
rz(-1.7705256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4930967) q[2];
sx q[2];
rz(-0.33005565) q[2];
sx q[2];
rz(-2.8685692) q[2];
rz(-1.8388883) q[3];
sx q[3];
rz(-1.8282993) q[3];
sx q[3];
rz(-0.31204143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.7438695) q[0];
sx q[0];
rz(-1.1345154) q[0];
sx q[0];
rz(1.2851108) q[0];
rz(-1.5015191) q[1];
sx q[1];
rz(-1.7506426) q[1];
sx q[1];
rz(-1.8008908) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1788951) q[0];
sx q[0];
rz(-0.51998752) q[0];
sx q[0];
rz(-2.2682701) q[0];
x q[1];
rz(-3.0807207) q[2];
sx q[2];
rz(-0.69673046) q[2];
sx q[2];
rz(0.75887647) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2015842) q[1];
sx q[1];
rz(-1.9218947) q[1];
sx q[1];
rz(-1.1854978) q[1];
x q[2];
rz(-2.826346) q[3];
sx q[3];
rz(-1.9936221) q[3];
sx q[3];
rz(-1.8002312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1061219) q[2];
sx q[2];
rz(-0.80646986) q[2];
sx q[2];
rz(2.1179874) q[2];
rz(-0.18493955) q[3];
sx q[3];
rz(-2.7513294) q[3];
sx q[3];
rz(-0.24188724) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0446562) q[0];
sx q[0];
rz(-0.99656314) q[0];
sx q[0];
rz(1.5203083) q[0];
rz(-2.8114491) q[1];
sx q[1];
rz(-1.9338927) q[1];
sx q[1];
rz(-2.3044589) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1820113) q[0];
sx q[0];
rz(-2.2416229) q[0];
sx q[0];
rz(1.3120033) q[0];
rz(0.43948549) q[2];
sx q[2];
rz(-0.62289933) q[2];
sx q[2];
rz(-0.57657951) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.21048927) q[1];
sx q[1];
rz(-1.466202) q[1];
sx q[1];
rz(1.4180257) q[1];
rz(-pi) q[2];
rz(-1.0650474) q[3];
sx q[3];
rz(-2.5878083) q[3];
sx q[3];
rz(1.6365285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.60124406) q[2];
sx q[2];
rz(-1.0849755) q[2];
sx q[2];
rz(2.9619651) q[2];
rz(2.1458697) q[3];
sx q[3];
rz(-1.2487753) q[3];
sx q[3];
rz(-1.8306336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3778465) q[0];
sx q[0];
rz(-0.34559956) q[0];
sx q[0];
rz(-2.0843704) q[0];
rz(-3.0341042) q[1];
sx q[1];
rz(-1.2534393) q[1];
sx q[1];
rz(0.9799788) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36469034) q[0];
sx q[0];
rz(-1.5199465) q[0];
sx q[0];
rz(-1.6965673) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0716295) q[2];
sx q[2];
rz(-1.1071148) q[2];
sx q[2];
rz(2.0641363) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.53303888) q[1];
sx q[1];
rz(-1.8095008) q[1];
sx q[1];
rz(-1.1281611) q[1];
rz(-1.8944593) q[3];
sx q[3];
rz(-1.4753046) q[3];
sx q[3];
rz(0.0088866339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.3048627) q[2];
sx q[2];
rz(-1.2048081) q[2];
sx q[2];
rz(-2.1137962) q[2];
rz(1.3868388) q[3];
sx q[3];
rz(-1.8324865) q[3];
sx q[3];
rz(-2.8579779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2733611) q[0];
sx q[0];
rz(-2.1049451) q[0];
sx q[0];
rz(2.0275397) q[0];
rz(-2.3095619) q[1];
sx q[1];
rz(-0.46453005) q[1];
sx q[1];
rz(0.66418905) q[1];
rz(-1.8503415) q[2];
sx q[2];
rz(-2.5646979) q[2];
sx q[2];
rz(-0.23451351) q[2];
rz(1.1888614) q[3];
sx q[3];
rz(-1.100913) q[3];
sx q[3];
rz(2.037896) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
