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
rz(-2.6053083) q[0];
sx q[0];
rz(0.94776881) q[0];
rz(1.8127958) q[1];
sx q[1];
rz(-1.2674018) q[1];
sx q[1];
rz(-1.0277494) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1411113) q[0];
sx q[0];
rz(-1.9465465) q[0];
sx q[0];
rz(1.5953338) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9704291) q[2];
sx q[2];
rz(-0.52414775) q[2];
sx q[2];
rz(-1.5607967) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9103968) q[1];
sx q[1];
rz(-2.3379571) q[1];
sx q[1];
rz(-0.5459783) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0184228) q[3];
sx q[3];
rz(-2.1493704) q[3];
sx q[3];
rz(2.6920126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1296922) q[2];
sx q[2];
rz(-1.7069858) q[2];
sx q[2];
rz(-1.0502846) q[2];
rz(1.1132647) q[3];
sx q[3];
rz(-0.89171019) q[3];
sx q[3];
rz(-3.0734857) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0691836) q[0];
sx q[0];
rz(-1.2453112) q[0];
sx q[0];
rz(0.29775277) q[0];
rz(-0.61966664) q[1];
sx q[1];
rz(-1.0071808) q[1];
sx q[1];
rz(2.0334977) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1381123) q[0];
sx q[0];
rz(-1.3618042) q[0];
sx q[0];
rz(-0.94603993) q[0];
x q[1];
rz(0.98845311) q[2];
sx q[2];
rz(-2.4465912) q[2];
sx q[2];
rz(-0.37441355) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.83924343) q[1];
sx q[1];
rz(-1.9560555) q[1];
sx q[1];
rz(-1.1883931) q[1];
rz(-pi) q[2];
rz(-1.3012582) q[3];
sx q[3];
rz(-2.2552935) q[3];
sx q[3];
rz(-2.136363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.46253282) q[2];
sx q[2];
rz(-0.97683895) q[2];
sx q[2];
rz(-0.97529808) q[2];
rz(2.2235928) q[3];
sx q[3];
rz(-1.2851597) q[3];
sx q[3];
rz(-0.29618922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96238962) q[0];
sx q[0];
rz(-2.2900892) q[0];
sx q[0];
rz(0.54291022) q[0];
rz(-2.2593598) q[1];
sx q[1];
rz(-1.135332) q[1];
sx q[1];
rz(0.96484819) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7643825) q[0];
sx q[0];
rz(-1.8922193) q[0];
sx q[0];
rz(1.8660603) q[0];
rz(-2.8446571) q[2];
sx q[2];
rz(-1.529971) q[2];
sx q[2];
rz(0.21619913) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.112903) q[1];
sx q[1];
rz(-1.5064872) q[1];
sx q[1];
rz(1.4494004) q[1];
rz(-pi) q[2];
x q[2];
rz(1.987625) q[3];
sx q[3];
rz(-1.7524476) q[3];
sx q[3];
rz(0.75205284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.09459153) q[2];
sx q[2];
rz(-0.61085218) q[2];
sx q[2];
rz(-1.1331406) q[2];
rz(-0.23162332) q[3];
sx q[3];
rz(-1.8685721) q[3];
sx q[3];
rz(2.384322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27947458) q[0];
sx q[0];
rz(-3.1311488) q[0];
sx q[0];
rz(1.3765155) q[0];
rz(-0.51849413) q[1];
sx q[1];
rz(-1.2644178) q[1];
sx q[1];
rz(-2.8994697) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3926485) q[0];
sx q[0];
rz(-0.90396515) q[0];
sx q[0];
rz(1.3632266) q[0];
rz(3.1130586) q[2];
sx q[2];
rz(-0.42978537) q[2];
sx q[2];
rz(-2.1536749) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.081269216) q[1];
sx q[1];
rz(-0.4821061) q[1];
sx q[1];
rz(-2.768885) q[1];
x q[2];
rz(-2.7235892) q[3];
sx q[3];
rz(-1.4417218) q[3];
sx q[3];
rz(-1.2152745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.018192856) q[2];
sx q[2];
rz(-0.96863666) q[2];
sx q[2];
rz(-3.0467395) q[2];
rz(1.799396) q[3];
sx q[3];
rz(-1.7443402) q[3];
sx q[3];
rz(-0.43911394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.78385329) q[0];
sx q[0];
rz(-1.8308324) q[0];
sx q[0];
rz(-2.7807996) q[0];
rz(-1.7533253) q[1];
sx q[1];
rz(-1.8107982) q[1];
sx q[1];
rz(1.1345908) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2312233) q[0];
sx q[0];
rz(-1.0783505) q[0];
sx q[0];
rz(-0.72427303) q[0];
rz(-0.88419948) q[2];
sx q[2];
rz(-2.6066385) q[2];
sx q[2];
rz(1.7810437) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7736588) q[1];
sx q[1];
rz(-1.7264257) q[1];
sx q[1];
rz(2.6917798) q[1];
x q[2];
rz(-2.9633425) q[3];
sx q[3];
rz(-2.6986487) q[3];
sx q[3];
rz(2.6435542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0052884) q[2];
sx q[2];
rz(-1.4207999) q[2];
sx q[2];
rz(-2.5689382) q[2];
rz(-0.92875656) q[3];
sx q[3];
rz(-2.6199665) q[3];
sx q[3];
rz(-1.9740392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8144433) q[0];
sx q[0];
rz(-1.0961908) q[0];
sx q[0];
rz(1.8922528) q[0];
rz(-1.2231187) q[1];
sx q[1];
rz(-1.616281) q[1];
sx q[1];
rz(-1.9893601) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5637701) q[0];
sx q[0];
rz(-1.676079) q[0];
sx q[0];
rz(-3.1288414) q[0];
rz(-pi) q[1];
rz(2.639159) q[2];
sx q[2];
rz(-1.9366493) q[2];
sx q[2];
rz(-0.94787129) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7631543) q[1];
sx q[1];
rz(-1.3707146) q[1];
sx q[1];
rz(-1.4640019) q[1];
rz(-0.96529393) q[3];
sx q[3];
rz(-2.1552342) q[3];
sx q[3];
rz(0.30129978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.46889177) q[2];
sx q[2];
rz(-1.2142618) q[2];
sx q[2];
rz(2.1506298) q[2];
rz(0.64783603) q[3];
sx q[3];
rz(-2.1760553) q[3];
sx q[3];
rz(-0.34753862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25794849) q[0];
sx q[0];
rz(-2.914496) q[0];
sx q[0];
rz(-0.062285475) q[0];
rz(2.9557872) q[1];
sx q[1];
rz(-1.6848911) q[1];
sx q[1];
rz(-2.7468162) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6900401) q[0];
sx q[0];
rz(-1.8071022) q[0];
sx q[0];
rz(-0.4060181) q[0];
x q[1];
rz(-2.4418418) q[2];
sx q[2];
rz(-0.47669461) q[2];
sx q[2];
rz(1.9112019) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2784087) q[1];
sx q[1];
rz(-1.3235958) q[1];
sx q[1];
rz(2.7017038) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23468252) q[3];
sx q[3];
rz(-0.15771401) q[3];
sx q[3];
rz(1.371067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4930967) q[2];
sx q[2];
rz(-2.811537) q[2];
sx q[2];
rz(-0.27302343) q[2];
rz(-1.8388883) q[3];
sx q[3];
rz(-1.8282993) q[3];
sx q[3];
rz(-0.31204143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7438695) q[0];
sx q[0];
rz(-1.1345154) q[0];
sx q[0];
rz(-1.8564818) q[0];
rz(1.5015191) q[1];
sx q[1];
rz(-1.7506426) q[1];
sx q[1];
rz(1.8008908) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9467981) q[0];
sx q[0];
rz(-1.9614944) q[0];
sx q[0];
rz(-2.7892053) q[0];
x q[1];
rz(2.4457744) q[2];
sx q[2];
rz(-1.5317481) q[2];
sx q[2];
rz(-2.2829636) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.50780523) q[1];
sx q[1];
rz(-1.9314737) q[1];
sx q[1];
rz(2.7651869) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0128485) q[3];
sx q[3];
rz(-1.2841409) q[3];
sx q[3];
rz(-3.0451881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1061219) q[2];
sx q[2];
rz(-2.3351228) q[2];
sx q[2];
rz(-2.1179874) q[2];
rz(0.18493955) q[3];
sx q[3];
rz(-0.39026323) q[3];
sx q[3];
rz(-0.24188724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0969365) q[0];
sx q[0];
rz(-2.1450295) q[0];
sx q[0];
rz(1.6212844) q[0];
rz(0.3301436) q[1];
sx q[1];
rz(-1.2076999) q[1];
sx q[1];
rz(2.3044589) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9159106) q[0];
sx q[0];
rz(-1.3689694) q[0];
sx q[0];
rz(2.4542144) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2741954) q[2];
sx q[2];
rz(-1.0146078) q[2];
sx q[2];
rz(2.0402758) q[2];
rz(pi/2) q[3];
sx q[3];
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
x q[2];
rz(-2.8505441) q[3];
sx q[3];
rz(-2.0488727) q[3];
sx q[3];
rz(-2.2136798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.60124406) q[2];
sx q[2];
rz(-1.0849755) q[2];
sx q[2];
rz(2.9619651) q[2];
rz(0.99572292) q[3];
sx q[3];
rz(-1.2487753) q[3];
sx q[3];
rz(1.8306336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(0.76374617) q[0];
sx q[0];
rz(-0.34559956) q[0];
sx q[0];
rz(-1.0572222) q[0];
rz(0.10748848) q[1];
sx q[1];
rz(-1.2534393) q[1];
sx q[1];
rz(-2.1616139) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5883334) q[0];
sx q[0];
rz(-3.0059814) q[0];
sx q[0];
rz(1.1853663) q[0];
rz(-pi) q[1];
rz(2.0699632) q[2];
sx q[2];
rz(-1.1071148) q[2];
sx q[2];
rz(-1.0774563) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6411297) q[1];
sx q[1];
rz(-0.49912057) q[1];
sx q[1];
rz(2.0874546) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8633217) q[3];
sx q[3];
rz(-0.33698002) q[3];
sx q[3];
rz(-1.8388336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8367299) q[2];
sx q[2];
rz(-1.9367846) q[2];
sx q[2];
rz(1.0277964) q[2];
rz(-1.7547539) q[3];
sx q[3];
rz(-1.8324865) q[3];
sx q[3];
rz(-2.8579779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86823157) q[0];
sx q[0];
rz(-2.1049451) q[0];
sx q[0];
rz(2.0275397) q[0];
rz(-2.3095619) q[1];
sx q[1];
rz(-0.46453005) q[1];
sx q[1];
rz(0.66418905) q[1];
rz(1.2912512) q[2];
sx q[2];
rz(-2.5646979) q[2];
sx q[2];
rz(-0.23451351) q[2];
rz(-1.9527312) q[3];
sx q[3];
rz(-1.100913) q[3];
sx q[3];
rz(2.037896) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
