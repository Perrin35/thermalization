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
rz(1.8127958) q[1];
sx q[1];
rz(-1.2674018) q[1];
sx q[1];
rz(2.1138432) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2078903) q[0];
sx q[0];
rz(-2.7650802) q[0];
sx q[0];
rz(-0.062113751) q[0];
x q[1];
rz(-0.22123863) q[2];
sx q[2];
rz(-2.049963) q[2];
sx q[2];
rz(2.0146807) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.4878405) q[1];
sx q[1];
rz(-2.2334705) q[1];
sx q[1];
rz(1.0767879) q[1];
x q[2];
rz(0.58524744) q[3];
sx q[3];
rz(-2.4260776) q[3];
sx q[3];
rz(1.9714718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0119005) q[2];
sx q[2];
rz(-1.4346069) q[2];
sx q[2];
rz(-1.0502846) q[2];
rz(-2.0283279) q[3];
sx q[3];
rz(-2.2498825) q[3];
sx q[3];
rz(3.0734857) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0691836) q[0];
sx q[0];
rz(-1.2453112) q[0];
sx q[0];
rz(-2.8438399) q[0];
rz(0.61966664) q[1];
sx q[1];
rz(-1.0071808) q[1];
sx q[1];
rz(1.108095) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0034804) q[0];
sx q[0];
rz(-1.7797884) q[0];
sx q[0];
rz(-0.94603993) q[0];
x q[1];
rz(0.42995288) q[2];
sx q[2];
rz(-1.0064831) q[2];
sx q[2];
rz(-2.058409) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.019921692) q[1];
sx q[1];
rz(-2.6056075) q[1];
sx q[1];
rz(-0.74384816) q[1];
rz(-pi) q[2];
x q[2];
rz(2.439019) q[3];
sx q[3];
rz(-1.7785903) q[3];
sx q[3];
rz(-0.39263615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.46253282) q[2];
sx q[2];
rz(-2.1647537) q[2];
sx q[2];
rz(-2.1662946) q[2];
rz(-0.9179999) q[3];
sx q[3];
rz(-1.2851597) q[3];
sx q[3];
rz(2.8454034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.179203) q[0];
sx q[0];
rz(-0.85150349) q[0];
sx q[0];
rz(2.5986824) q[0];
rz(-0.88223282) q[1];
sx q[1];
rz(-1.135332) q[1];
sx q[1];
rz(-0.96484819) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7643825) q[0];
sx q[0];
rz(-1.2493734) q[0];
sx q[0];
rz(1.2755323) q[0];
x q[1];
rz(2.8446571) q[2];
sx q[2];
rz(-1.6116217) q[2];
sx q[2];
rz(0.21619913) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1988586) q[1];
sx q[1];
rz(-3.004289) q[1];
sx q[1];
rz(2.0595466) q[1];
x q[2];
rz(-1.987625) q[3];
sx q[3];
rz(-1.7524476) q[3];
sx q[3];
rz(2.3895398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.09459153) q[2];
sx q[2];
rz(-2.5307405) q[2];
sx q[2];
rz(2.0084521) q[2];
rz(0.23162332) q[3];
sx q[3];
rz(-1.2730205) q[3];
sx q[3];
rz(2.384322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27947458) q[0];
sx q[0];
rz(-3.1311488) q[0];
sx q[0];
rz(1.7650771) q[0];
rz(-0.51849413) q[1];
sx q[1];
rz(-1.2644178) q[1];
sx q[1];
rz(-2.8994697) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0771368) q[0];
sx q[0];
rz(-2.4479439) q[0];
sx q[0];
rz(-2.8855188) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5838727) q[2];
sx q[2];
rz(-1.1411975) q[2];
sx q[2];
rz(-1.0193046) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.081269216) q[1];
sx q[1];
rz(-2.6594866) q[1];
sx q[1];
rz(2.768885) q[1];
rz(-pi) q[2];
rz(-1.4297156) q[3];
sx q[3];
rz(-1.1564848) q[3];
sx q[3];
rz(-0.29841081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1233998) q[2];
sx q[2];
rz(-2.172956) q[2];
sx q[2];
rz(-3.0467395) q[2];
rz(1.3421966) q[3];
sx q[3];
rz(-1.3972524) q[3];
sx q[3];
rz(2.7024787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78385329) q[0];
sx q[0];
rz(-1.3107603) q[0];
sx q[0];
rz(-0.36079303) q[0];
rz(-1.3882673) q[1];
sx q[1];
rz(-1.3307945) q[1];
sx q[1];
rz(1.1345908) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15123385) q[0];
sx q[0];
rz(-2.2915974) q[0];
sx q[0];
rz(2.4609341) q[0];
rz(-pi) q[1];
rz(-0.35933944) q[2];
sx q[2];
rz(-1.1655072) q[2];
sx q[2];
rz(-2.5422424) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1281631) q[1];
sx q[1];
rz(-2.014782) q[1];
sx q[1];
rz(-1.7432937) q[1];
x q[2];
rz(2.9633425) q[3];
sx q[3];
rz(-2.6986487) q[3];
sx q[3];
rz(0.49803842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1363042) q[2];
sx q[2];
rz(-1.7207928) q[2];
sx q[2];
rz(-0.57265442) q[2];
rz(0.92875656) q[3];
sx q[3];
rz(-0.52162617) q[3];
sx q[3];
rz(1.1675534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3271493) q[0];
sx q[0];
rz(-1.0961908) q[0];
sx q[0];
rz(1.8922528) q[0];
rz(1.918474) q[1];
sx q[1];
rz(-1.5253116) q[1];
sx q[1];
rz(-1.1522326) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6845247) q[0];
sx q[0];
rz(-3.0355434) q[0];
sx q[0];
rz(1.690879) q[0];
rz(-2.639159) q[2];
sx q[2];
rz(-1.9366493) q[2];
sx q[2];
rz(0.94787129) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7631543) q[1];
sx q[1];
rz(-1.3707146) q[1];
sx q[1];
rz(1.6775908) q[1];
rz(-pi) q[2];
rz(2.4640718) q[3];
sx q[3];
rz(-2.065425) q[3];
sx q[3];
rz(2.236931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.46889177) q[2];
sx q[2];
rz(-1.9273309) q[2];
sx q[2];
rz(2.1506298) q[2];
rz(-2.4937566) q[3];
sx q[3];
rz(-0.9655374) q[3];
sx q[3];
rz(-2.794054) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25794849) q[0];
sx q[0];
rz(-2.914496) q[0];
sx q[0];
rz(0.062285475) q[0];
rz(-0.1858055) q[1];
sx q[1];
rz(-1.4567016) q[1];
sx q[1];
rz(2.7468162) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5236854) q[0];
sx q[0];
rz(-0.46645188) q[0];
sx q[0];
rz(-2.594069) q[0];
x q[1];
rz(-0.6997509) q[2];
sx q[2];
rz(-2.664898) q[2];
sx q[2];
rz(1.9112019) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.77183206) q[1];
sx q[1];
rz(-0.50060111) q[1];
sx q[1];
rz(0.53497772) q[1];
rz(-pi) q[2];
rz(1.5338321) q[3];
sx q[3];
rz(-1.4174403) q[3];
sx q[3];
rz(-1.6085898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4930967) q[2];
sx q[2];
rz(-2.811537) q[2];
sx q[2];
rz(-0.27302343) q[2];
rz(-1.3027044) q[3];
sx q[3];
rz(-1.8282993) q[3];
sx q[3];
rz(-2.8295512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7438695) q[0];
sx q[0];
rz(-2.0070772) q[0];
sx q[0];
rz(-1.8564818) q[0];
rz(1.6400736) q[1];
sx q[1];
rz(-1.7506426) q[1];
sx q[1];
rz(-1.8008908) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9467981) q[0];
sx q[0];
rz(-1.9614944) q[0];
sx q[0];
rz(-0.35238738) q[0];
rz(0.060872002) q[2];
sx q[2];
rz(-0.69673046) q[2];
sx q[2];
rz(0.75887647) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.50780523) q[1];
sx q[1];
rz(-1.2101189) q[1];
sx q[1];
rz(-0.37640576) q[1];
x q[2];
rz(-0.31524661) q[3];
sx q[3];
rz(-1.1479706) q[3];
sx q[3];
rz(1.3413615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0354707) q[2];
sx q[2];
rz(-2.3351228) q[2];
sx q[2];
rz(1.0236053) q[2];
rz(-0.18493955) q[3];
sx q[3];
rz(-2.7513294) q[3];
sx q[3];
rz(2.8997054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
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
rz(2.8114491) q[1];
sx q[1];
rz(-1.9338927) q[1];
sx q[1];
rz(2.3044589) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9159106) q[0];
sx q[0];
rz(-1.3689694) q[0];
sx q[0];
rz(2.4542144) q[0];
rz(-1.2741954) q[2];
sx q[2];
rz(-1.0146078) q[2];
sx q[2];
rz(-2.0402758) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7652119) q[1];
sx q[1];
rz(-1.7227255) q[1];
sx q[1];
rz(-0.10581776) q[1];
rz(-pi) q[2];
x q[2];
rz(0.29104851) q[3];
sx q[3];
rz(-2.0488727) q[3];
sx q[3];
rz(0.92791286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5403486) q[2];
sx q[2];
rz(-2.0566172) q[2];
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
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3778465) q[0];
sx q[0];
rz(-0.34559956) q[0];
sx q[0];
rz(-1.0572222) q[0];
rz(3.0341042) q[1];
sx q[1];
rz(-1.2534393) q[1];
sx q[1];
rz(-0.9799788) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5532593) q[0];
sx q[0];
rz(-3.0059814) q[0];
sx q[0];
rz(-1.1853663) q[0];
rz(0.51771848) q[2];
sx q[2];
rz(-1.1283518) q[2];
sx q[2];
rz(2.8874318) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.92614782) q[1];
sx q[1];
rz(-2.0000334) q[1];
sx q[1];
rz(-2.8785359) q[1];
rz(-pi) q[2];
rz(-0.1006871) q[3];
sx q[3];
rz(-1.8929314) q[3];
sx q[3];
rz(-1.5299357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3048627) q[2];
sx q[2];
rz(-1.2048081) q[2];
sx q[2];
rz(2.1137962) q[2];
rz(-1.7547539) q[3];
sx q[3];
rz(-1.8324865) q[3];
sx q[3];
rz(-2.8579779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
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
rz(2.3095619) q[1];
sx q[1];
rz(-2.6770626) q[1];
sx q[1];
rz(-2.4774036) q[1];
rz(-1.8503415) q[2];
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
