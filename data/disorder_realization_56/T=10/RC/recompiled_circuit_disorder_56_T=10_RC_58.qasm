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
rz(-0.94776881) q[0];
rz(1.8127958) q[1];
sx q[1];
rz(-1.2674018) q[1];
sx q[1];
rz(-1.0277494) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9337024) q[0];
sx q[0];
rz(-0.37651248) q[0];
sx q[0];
rz(-0.062113751) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22123863) q[2];
sx q[2];
rz(-2.049963) q[2];
sx q[2];
rz(1.1269119) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7386908) q[1];
sx q[1];
rz(-1.1876904) q[1];
sx q[1];
rz(0.72523592) q[1];
rz(-pi) q[2];
rz(2.0184228) q[3];
sx q[3];
rz(-0.99222224) q[3];
sx q[3];
rz(2.6920126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0119005) q[2];
sx q[2];
rz(-1.4346069) q[2];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0691836) q[0];
sx q[0];
rz(-1.8962815) q[0];
sx q[0];
rz(2.8438399) q[0];
rz(0.61966664) q[1];
sx q[1];
rz(-2.1344118) q[1];
sx q[1];
rz(2.0334977) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0034804) q[0];
sx q[0];
rz(-1.7797884) q[0];
sx q[0];
rz(-0.94603993) q[0];
rz(-0.96252243) q[2];
sx q[2];
rz(-1.2108742) q[2];
sx q[2];
rz(0.72812176) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.121671) q[1];
sx q[1];
rz(-0.53598511) q[1];
sx q[1];
rz(-0.74384816) q[1];
x q[2];
rz(0.31538972) q[3];
sx q[3];
rz(-2.4139801) q[3];
sx q[3];
rz(1.7244347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6790598) q[2];
sx q[2];
rz(-0.97683895) q[2];
sx q[2];
rz(2.1662946) q[2];
rz(0.9179999) q[3];
sx q[3];
rz(-1.8564329) q[3];
sx q[3];
rz(2.8454034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.96238962) q[0];
sx q[0];
rz(-0.85150349) q[0];
sx q[0];
rz(-0.54291022) q[0];
rz(-0.88223282) q[1];
sx q[1];
rz(-2.0062607) q[1];
sx q[1];
rz(-2.1767445) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7643825) q[0];
sx q[0];
rz(-1.2493734) q[0];
sx q[0];
rz(-1.2755323) q[0];
rz(-1.5281048) q[2];
sx q[2];
rz(-1.2741158) q[2];
sx q[2];
rz(-1.7745078) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0286897) q[1];
sx q[1];
rz(-1.5064872) q[1];
sx q[1];
rz(1.6921922) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9433555) q[3];
sx q[3];
rz(-1.9803515) q[3];
sx q[3];
rz(2.4026681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27947458) q[0];
sx q[0];
rz(-0.010443895) q[0];
sx q[0];
rz(-1.3765155) q[0];
rz(2.6230985) q[1];
sx q[1];
rz(-1.2644178) q[1];
sx q[1];
rz(-2.8994697) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.449269) q[0];
sx q[0];
rz(-1.4081435) q[0];
sx q[0];
rz(-0.67740324) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7119615) q[2];
sx q[2];
rz(-1.5589082) q[2];
sx q[2];
rz(-0.55693835) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.081269216) q[1];
sx q[1];
rz(-0.4821061) q[1];
sx q[1];
rz(2.768885) q[1];
rz(2.7235892) q[3];
sx q[3];
rz(-1.6998708) q[3];
sx q[3];
rz(-1.2152745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1233998) q[2];
sx q[2];
rz(-2.172956) q[2];
sx q[2];
rz(0.094853178) q[2];
rz(1.799396) q[3];
sx q[3];
rz(-1.7443402) q[3];
sx q[3];
rz(2.7024787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-1.7533253) q[1];
sx q[1];
rz(-1.8107982) q[1];
sx q[1];
rz(-2.0070019) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9103694) q[0];
sx q[0];
rz(-1.0783505) q[0];
sx q[0];
rz(-2.4173196) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7822532) q[2];
sx q[2];
rz(-1.1655072) q[2];
sx q[2];
rz(-0.59935024) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6280917) q[1];
sx q[1];
rz(-0.47422945) q[1];
sx q[1];
rz(2.7952816) q[1];
rz(-1.4868823) q[3];
sx q[3];
rz(-1.1353555) q[3];
sx q[3];
rz(-0.30121379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0052884) q[2];
sx q[2];
rz(-1.4207999) q[2];
sx q[2];
rz(2.5689382) q[2];
rz(2.2128361) q[3];
sx q[3];
rz(-0.52162617) q[3];
sx q[3];
rz(-1.1675534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.8144433) q[0];
sx q[0];
rz(-1.0961908) q[0];
sx q[0];
rz(1.2493398) q[0];
rz(-1.2231187) q[1];
sx q[1];
rz(-1.616281) q[1];
sx q[1];
rz(1.1522326) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5637701) q[0];
sx q[0];
rz(-1.4655136) q[0];
sx q[0];
rz(-3.1288414) q[0];
rz(0.67201519) q[2];
sx q[2];
rz(-2.5294371) q[2];
sx q[2];
rz(1.9415346) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.17105477) q[1];
sx q[1];
rz(-1.6754524) q[1];
sx q[1];
rz(-0.20119757) q[1];
rz(-pi) q[2];
rz(-2.1762987) q[3];
sx q[3];
rz(-2.1552342) q[3];
sx q[3];
rz(2.8402929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6727009) q[2];
sx q[2];
rz(-1.2142618) q[2];
sx q[2];
rz(2.1506298) q[2];
rz(2.4937566) q[3];
sx q[3];
rz(-0.9655374) q[3];
sx q[3];
rz(-0.34753862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25794849) q[0];
sx q[0];
rz(-2.914496) q[0];
sx q[0];
rz(-0.062285475) q[0];
rz(-2.9557872) q[1];
sx q[1];
rz(-1.6848911) q[1];
sx q[1];
rz(-0.3947765) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61790723) q[0];
sx q[0];
rz(-2.6751408) q[0];
sx q[0];
rz(0.54752366) q[0];
x q[1];
rz(2.4418418) q[2];
sx q[2];
rz(-0.47669461) q[2];
sx q[2];
rz(-1.9112019) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3697606) q[1];
sx q[1];
rz(-2.6409915) q[1];
sx q[1];
rz(0.53497772) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9069101) q[3];
sx q[3];
rz(-2.9838786) q[3];
sx q[3];
rz(-1.371067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.64849598) q[2];
sx q[2];
rz(-2.811537) q[2];
sx q[2];
rz(2.8685692) q[2];
rz(-1.3027044) q[3];
sx q[3];
rz(-1.3132934) q[3];
sx q[3];
rz(2.8295512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39772314) q[0];
sx q[0];
rz(-1.1345154) q[0];
sx q[0];
rz(1.2851108) q[0];
rz(1.5015191) q[1];
sx q[1];
rz(-1.39095) q[1];
sx q[1];
rz(1.3407019) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1947945) q[0];
sx q[0];
rz(-1.1800982) q[0];
sx q[0];
rz(0.35238738) q[0];
rz(-pi) q[1];
rz(1.6216535) q[2];
sx q[2];
rz(-2.2659781) q[2];
sx q[2];
rz(0.67957544) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2015842) q[1];
sx q[1];
rz(-1.9218947) q[1];
sx q[1];
rz(-1.9560948) q[1];
rz(2.0128485) q[3];
sx q[3];
rz(-1.2841409) q[3];
sx q[3];
rz(-0.096404508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0354707) q[2];
sx q[2];
rz(-2.3351228) q[2];
sx q[2];
rz(-1.0236053) q[2];
rz(-0.18493955) q[3];
sx q[3];
rz(-2.7513294) q[3];
sx q[3];
rz(-0.24188724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0969365) q[0];
sx q[0];
rz(-0.99656314) q[0];
sx q[0];
rz(1.6212844) q[0];
rz(2.8114491) q[1];
sx q[1];
rz(-1.2076999) q[1];
sx q[1];
rz(-2.3044589) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9595813) q[0];
sx q[0];
rz(-2.2416229) q[0];
sx q[0];
rz(-1.3120033) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2741954) q[2];
sx q[2];
rz(-1.0146078) q[2];
sx q[2];
rz(1.1013168) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3771462) q[1];
sx q[1];
rz(-0.18491491) q[1];
sx q[1];
rz(0.96692337) q[1];
rz(-1.0650474) q[3];
sx q[3];
rz(-0.55378434) q[3];
sx q[3];
rz(-1.6365285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.60124406) q[2];
sx q[2];
rz(-2.0566172) q[2];
sx q[2];
rz(0.17962757) q[2];
rz(0.99572292) q[3];
sx q[3];
rz(-1.8928173) q[3];
sx q[3];
rz(1.3109591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3778465) q[0];
sx q[0];
rz(-2.7959931) q[0];
sx q[0];
rz(-1.0572222) q[0];
rz(-0.10748848) q[1];
sx q[1];
rz(-1.2534393) q[1];
sx q[1];
rz(2.1616139) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1996795) q[0];
sx q[0];
rz(-1.4451888) q[0];
sx q[0];
rz(-3.0903387) q[0];
rz(0.51771848) q[2];
sx q[2];
rz(-1.1283518) q[2];
sx q[2];
rz(-0.25416086) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.500463) q[1];
sx q[1];
rz(-0.49912057) q[1];
sx q[1];
rz(2.0874546) q[1];
rz(-0.1006871) q[3];
sx q[3];
rz(-1.2486613) q[3];
sx q[3];
rz(1.5299357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8367299) q[2];
sx q[2];
rz(-1.2048081) q[2];
sx q[2];
rz(-2.1137962) q[2];
rz(1.7547539) q[3];
sx q[3];
rz(-1.3091062) q[3];
sx q[3];
rz(0.28361472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.1297395) q[2];
sx q[2];
rz(-1.7218628) q[2];
sx q[2];
rz(-2.0414258) q[2];
rz(-2.5084393) q[3];
sx q[3];
rz(-0.59637759) q[3];
sx q[3];
rz(-0.37806088) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
