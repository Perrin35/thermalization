OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.7528485) q[0];
sx q[0];
rz(-0.53628439) q[0];
sx q[0];
rz(2.1938238) q[0];
rz(1.8127958) q[1];
sx q[1];
rz(-1.2674018) q[1];
sx q[1];
rz(2.1138432) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9337024) q[0];
sx q[0];
rz(-2.7650802) q[0];
sx q[0];
rz(0.062113751) q[0];
rz(-1.9704291) q[2];
sx q[2];
rz(-0.52414775) q[2];
sx q[2];
rz(-1.580796) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6537522) q[1];
sx q[1];
rz(-2.2334705) q[1];
sx q[1];
rz(-1.0767879) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5146033) q[3];
sx q[3];
rz(-1.1999745) q[3];
sx q[3];
rz(-0.8644608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1296922) q[2];
sx q[2];
rz(-1.4346069) q[2];
sx q[2];
rz(1.0502846) q[2];
rz(2.0283279) q[3];
sx q[3];
rz(-2.2498825) q[3];
sx q[3];
rz(-3.0734857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.072409078) q[0];
sx q[0];
rz(-1.2453112) q[0];
sx q[0];
rz(2.8438399) q[0];
rz(0.61966664) q[1];
sx q[1];
rz(-2.1344118) q[1];
sx q[1];
rz(-1.108095) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4257658) q[0];
sx q[0];
rz(-2.1799488) q[0];
sx q[0];
rz(-0.25575511) q[0];
rz(-pi) q[1];
rz(2.1790702) q[2];
sx q[2];
rz(-1.2108742) q[2];
sx q[2];
rz(0.72812176) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3023492) q[1];
sx q[1];
rz(-1.1855372) q[1];
sx q[1];
rz(-1.1883931) q[1];
rz(-pi) q[2];
rz(0.70257367) q[3];
sx q[3];
rz(-1.3630023) q[3];
sx q[3];
rz(2.7489565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46253282) q[2];
sx q[2];
rz(-2.1647537) q[2];
sx q[2];
rz(0.97529808) q[2];
rz(2.2235928) q[3];
sx q[3];
rz(-1.8564329) q[3];
sx q[3];
rz(-2.8454034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96238962) q[0];
sx q[0];
rz(-2.2900892) q[0];
sx q[0];
rz(2.5986824) q[0];
rz(-2.2593598) q[1];
sx q[1];
rz(-2.0062607) q[1];
sx q[1];
rz(2.1767445) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0438) q[0];
sx q[0];
rz(-1.291073) q[0];
sx q[0];
rz(-2.8066737) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.13871128) q[2];
sx q[2];
rz(-2.841946) q[2];
sx q[2];
rz(-1.2219929) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.94273401) q[1];
sx q[1];
rz(-0.13730362) q[1];
sx q[1];
rz(-2.0595466) q[1];
rz(-pi) q[2];
rz(-2.9433555) q[3];
sx q[3];
rz(-1.1612411) q[3];
sx q[3];
rz(-0.73892456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0470011) q[2];
sx q[2];
rz(-0.61085218) q[2];
sx q[2];
rz(-2.0084521) q[2];
rz(-0.23162332) q[3];
sx q[3];
rz(-1.8685721) q[3];
sx q[3];
rz(-0.75727063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8621181) q[0];
sx q[0];
rz(-3.1311488) q[0];
sx q[0];
rz(1.3765155) q[0];
rz(-0.51849413) q[1];
sx q[1];
rz(-1.8771749) q[1];
sx q[1];
rz(2.8994697) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0771368) q[0];
sx q[0];
rz(-2.4479439) q[0];
sx q[0];
rz(0.25607381) q[0];
rz(-1.5838727) q[2];
sx q[2];
rz(-2.0003951) q[2];
sx q[2];
rz(-2.122288) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8072847) q[1];
sx q[1];
rz(-2.0172999) q[1];
sx q[1];
rz(1.3825033) q[1];
rz(2.8321213) q[3];
sx q[3];
rz(-2.7052393) q[3];
sx q[3];
rz(-2.5040124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1233998) q[2];
sx q[2];
rz(-0.96863666) q[2];
sx q[2];
rz(3.0467395) q[2];
rz(-1.3421966) q[3];
sx q[3];
rz(-1.3972524) q[3];
sx q[3];
rz(-2.7024787) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78385329) q[0];
sx q[0];
rz(-1.8308324) q[0];
sx q[0];
rz(-0.36079303) q[0];
rz(1.3882673) q[1];
sx q[1];
rz(-1.8107982) q[1];
sx q[1];
rz(1.1345908) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73571262) q[0];
sx q[0];
rz(-2.1942733) q[0];
sx q[0];
rz(2.1924125) q[0];
rz(-pi) q[1];
rz(2.2573932) q[2];
sx q[2];
rz(-2.6066385) q[2];
sx q[2];
rz(1.7810437) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.36793383) q[1];
sx q[1];
rz(-1.7264257) q[1];
sx q[1];
rz(-0.44981287) q[1];
rz(-pi) q[2];
rz(-1.6547104) q[3];
sx q[3];
rz(-1.1353555) q[3];
sx q[3];
rz(0.30121379) q[3];
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
rz(0.57265442) q[2];
rz(-2.2128361) q[3];
sx q[3];
rz(-2.6199665) q[3];
sx q[3];
rz(-1.1675534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.8144433) q[0];
sx q[0];
rz(-2.0454018) q[0];
sx q[0];
rz(1.2493398) q[0];
rz(-1.918474) q[1];
sx q[1];
rz(-1.5253116) q[1];
sx q[1];
rz(-1.9893601) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5778225) q[0];
sx q[0];
rz(-1.4655136) q[0];
sx q[0];
rz(-3.1288414) q[0];
rz(-pi) q[1];
rz(2.4695775) q[2];
sx q[2];
rz(-2.5294371) q[2];
sx q[2];
rz(-1.9415346) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2684979) q[1];
sx q[1];
rz(-2.9151306) q[1];
sx q[1];
rz(0.48392673) q[1];
rz(0.67752083) q[3];
sx q[3];
rz(-2.065425) q[3];
sx q[3];
rz(-2.236931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.46889177) q[2];
sx q[2];
rz(-1.2142618) q[2];
sx q[2];
rz(0.99096283) q[2];
rz(-0.64783603) q[3];
sx q[3];
rz(-0.9655374) q[3];
sx q[3];
rz(-0.34753862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.25794849) q[0];
sx q[0];
rz(-2.914496) q[0];
sx q[0];
rz(0.062285475) q[0];
rz(2.9557872) q[1];
sx q[1];
rz(-1.4567016) q[1];
sx q[1];
rz(-0.3947765) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61790723) q[0];
sx q[0];
rz(-2.6751408) q[0];
sx q[0];
rz(2.594069) q[0];
x q[1];
rz(2.4418418) q[2];
sx q[2];
rz(-0.47669461) q[2];
sx q[2];
rz(-1.9112019) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.863184) q[1];
sx q[1];
rz(-1.8179968) q[1];
sx q[1];
rz(-2.7017038) q[1];
rz(-pi) q[2];
x q[2];
rz(0.15345927) q[3];
sx q[3];
rz(-1.5342661) q[3];
sx q[3];
rz(3.1094482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64849598) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7438695) q[0];
sx q[0];
rz(-2.0070772) q[0];
sx q[0];
rz(1.8564818) q[0];
rz(1.6400736) q[1];
sx q[1];
rz(-1.39095) q[1];
sx q[1];
rz(1.8008908) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9467981) q[0];
sx q[0];
rz(-1.1800982) q[0];
sx q[0];
rz(0.35238738) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6216535) q[2];
sx q[2];
rz(-0.87561456) q[2];
sx q[2];
rz(2.4620172) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33412877) q[1];
sx q[1];
rz(-0.51528105) q[1];
sx q[1];
rz(-2.3433102) q[1];
rz(-0.31524661) q[3];
sx q[3];
rz(-1.1479706) q[3];
sx q[3];
rz(-1.8002312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1061219) q[2];
sx q[2];
rz(-2.3351228) q[2];
sx q[2];
rz(-2.1179874) q[2];
rz(-0.18493955) q[3];
sx q[3];
rz(-2.7513294) q[3];
sx q[3];
rz(2.8997054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.0969365) q[0];
sx q[0];
rz(-0.99656314) q[0];
sx q[0];
rz(-1.5203083) q[0];
rz(-0.3301436) q[1];
sx q[1];
rz(-1.2076999) q[1];
sx q[1];
rz(0.83713371) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1820113) q[0];
sx q[0];
rz(-2.2416229) q[0];
sx q[0];
rz(1.8295893) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5651921) q[2];
sx q[2];
rz(-1.3199558) q[2];
sx q[2];
rz(-2.5121411) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.76444641) q[1];
sx q[1];
rz(-0.18491491) q[1];
sx q[1];
rz(-0.96692337) q[1];
rz(2.0666396) q[3];
sx q[3];
rz(-1.3131724) q[3];
sx q[3];
rz(2.6356634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.60124406) q[2];
sx q[2];
rz(-1.0849755) q[2];
sx q[2];
rz(-2.9619651) q[2];
rz(0.99572292) q[3];
sx q[3];
rz(-1.8928173) q[3];
sx q[3];
rz(1.3109591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(2.1616139) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1996795) q[0];
sx q[0];
rz(-1.6964039) q[0];
sx q[0];
rz(-0.051254) q[0];
rz(2.3780031) q[2];
sx q[2];
rz(-0.66765235) q[2];
sx q[2];
rz(-1.9612567) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.53303888) q[1];
sx q[1];
rz(-1.3320919) q[1];
sx q[1];
rz(-2.0134316) q[1];
x q[2];
rz(-1.2471334) q[3];
sx q[3];
rz(-1.4753046) q[3];
sx q[3];
rz(-0.0088866339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.3048627) q[2];
sx q[2];
rz(-1.9367846) q[2];
sx q[2];
rz(2.1137962) q[2];
rz(1.7547539) q[3];
sx q[3];
rz(-1.8324865) q[3];
sx q[3];
rz(2.8579779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.1297395) q[2];
sx q[2];
rz(-1.4197299) q[2];
sx q[2];
rz(1.1001669) q[2];
rz(2.5084393) q[3];
sx q[3];
rz(-2.5452151) q[3];
sx q[3];
rz(2.7635318) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];