OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4053722) q[0];
sx q[0];
rz(-0.020981941) q[0];
sx q[0];
rz(1.9606645) q[0];
rz(2.6612072) q[1];
sx q[1];
rz(5.1279629) q[1];
sx q[1];
rz(9.8915107) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8202756) q[0];
sx q[0];
rz(-1.2877687) q[0];
sx q[0];
rz(-2.3664066) q[0];
x q[1];
rz(-2.4619815) q[2];
sx q[2];
rz(-0.70290297) q[2];
sx q[2];
rz(1.694569) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4157652) q[1];
sx q[1];
rz(-1.5752324) q[1];
sx q[1];
rz(-1.8091101) q[1];
x q[2];
rz(-1.0084413) q[3];
sx q[3];
rz(-2.1690024) q[3];
sx q[3];
rz(1.4149208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0878318) q[2];
sx q[2];
rz(-1.4542645) q[2];
sx q[2];
rz(1.9161179) q[2];
rz(2.7671704) q[3];
sx q[3];
rz(-2.008308) q[3];
sx q[3];
rz(2.5882914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9541009) q[0];
sx q[0];
rz(-0.88045374) q[0];
sx q[0];
rz(-2.6075897) q[0];
rz(2.7502637) q[1];
sx q[1];
rz(-2.6983039) q[1];
sx q[1];
rz(-2.2543529) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2132062) q[0];
sx q[0];
rz(-1.9149826) q[0];
sx q[0];
rz(1.8977988) q[0];
rz(-pi) q[1];
rz(0.5440622) q[2];
sx q[2];
rz(-2.3070222) q[2];
sx q[2];
rz(2.6092333) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.45184093) q[1];
sx q[1];
rz(-1.4418169) q[1];
sx q[1];
rz(-2.9366628) q[1];
x q[2];
rz(-0.34167413) q[3];
sx q[3];
rz(-1.0782584) q[3];
sx q[3];
rz(1.1600174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5947764) q[2];
sx q[2];
rz(-1.8211326) q[2];
sx q[2];
rz(-2.7776264) q[2];
rz(-1.4173896) q[3];
sx q[3];
rz(-0.28984362) q[3];
sx q[3];
rz(-2.3072306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.904838) q[0];
sx q[0];
rz(-3.0433488) q[0];
sx q[0];
rz(-0.33016095) q[0];
rz(-2.5579021) q[1];
sx q[1];
rz(-2.3287562) q[1];
sx q[1];
rz(-2.6469753) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30545091) q[0];
sx q[0];
rz(-1.8343121) q[0];
sx q[0];
rz(-1.8515153) q[0];
x q[1];
rz(-1.0801267) q[2];
sx q[2];
rz(-1.3441836) q[2];
sx q[2];
rz(1.8267711) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.45060748) q[1];
sx q[1];
rz(-0.71390426) q[1];
sx q[1];
rz(2.3074478) q[1];
rz(-pi) q[2];
rz(-0.98683896) q[3];
sx q[3];
rz(-1.6641326) q[3];
sx q[3];
rz(-0.10658857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.79564774) q[2];
sx q[2];
rz(-2.4730885) q[2];
sx q[2];
rz(3.0188959) q[2];
rz(0.042938558) q[3];
sx q[3];
rz(-1.0791082) q[3];
sx q[3];
rz(2.9470961) q[3];
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
rz(0.73404679) q[0];
sx q[0];
rz(-2.6341944) q[0];
sx q[0];
rz(-0.88974446) q[0];
rz(2.6906158) q[1];
sx q[1];
rz(-0.86800066) q[1];
sx q[1];
rz(-0.45423347) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.736859) q[0];
sx q[0];
rz(-1.4643719) q[0];
sx q[0];
rz(0.34261301) q[0];
x q[1];
rz(0.42734442) q[2];
sx q[2];
rz(-1.5691535) q[2];
sx q[2];
rz(2.9438007) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.046084) q[1];
sx q[1];
rz(-0.81611505) q[1];
sx q[1];
rz(2.7296394) q[1];
rz(0.20619456) q[3];
sx q[3];
rz(-0.36199328) q[3];
sx q[3];
rz(-2.9235206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.98895994) q[2];
sx q[2];
rz(-2.5071414) q[2];
sx q[2];
rz(-0.75759849) q[2];
rz(0.18975137) q[3];
sx q[3];
rz(-1.8228143) q[3];
sx q[3];
rz(2.5943713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2106638) q[0];
sx q[0];
rz(-2.3543816) q[0];
sx q[0];
rz(-1.8481365) q[0];
rz(-2.5594607) q[1];
sx q[1];
rz(-1.7030092) q[1];
sx q[1];
rz(0.17939803) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1360649) q[0];
sx q[0];
rz(-1.3910157) q[0];
sx q[0];
rz(1.3639569) q[0];
rz(-pi) q[1];
rz(1.8151692) q[2];
sx q[2];
rz(-0.91460157) q[2];
sx q[2];
rz(2.6442621) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0696358) q[1];
sx q[1];
rz(-1.3838125) q[1];
sx q[1];
rz(0.93902875) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6396181) q[3];
sx q[3];
rz(-2.2036512) q[3];
sx q[3];
rz(-2.8090734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3680215) q[2];
sx q[2];
rz(-2.7715235) q[2];
sx q[2];
rz(2.5717112) q[2];
rz(0.44024769) q[3];
sx q[3];
rz(-1.8972242) q[3];
sx q[3];
rz(0.35278916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2920947) q[0];
sx q[0];
rz(-2.8773913) q[0];
sx q[0];
rz(-1.1415035) q[0];
rz(1.3453206) q[1];
sx q[1];
rz(-0.99015403) q[1];
sx q[1];
rz(-0.17123953) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8526588) q[0];
sx q[0];
rz(-1.4698514) q[0];
sx q[0];
rz(0.48152083) q[0];
rz(-pi) q[1];
rz(2.0948579) q[2];
sx q[2];
rz(-1.9643133) q[2];
sx q[2];
rz(1.6167906) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6936612) q[1];
sx q[1];
rz(-2.1235925) q[1];
sx q[1];
rz(-1.3752244) q[1];
x q[2];
rz(-2.1275131) q[3];
sx q[3];
rz(-1.5601592) q[3];
sx q[3];
rz(0.47350804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3693927) q[2];
sx q[2];
rz(-2.2056396) q[2];
sx q[2];
rz(3.0041223) q[2];
rz(1.2482268) q[3];
sx q[3];
rz(-0.56709254) q[3];
sx q[3];
rz(1.221777) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0565979) q[0];
sx q[0];
rz(-2.5039112) q[0];
sx q[0];
rz(1.4884663) q[0];
rz(0.78530637) q[1];
sx q[1];
rz(-1.7410802) q[1];
sx q[1];
rz(1.3135501) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1378325) q[0];
sx q[0];
rz(-1.6562471) q[0];
sx q[0];
rz(-1.6944044) q[0];
rz(-pi) q[1];
rz(2.9503472) q[2];
sx q[2];
rz(-1.4594263) q[2];
sx q[2];
rz(-1.065606) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8911095) q[1];
sx q[1];
rz(-1.4817855) q[1];
sx q[1];
rz(-1.5989892) q[1];
rz(-pi) q[2];
rz(0.87487674) q[3];
sx q[3];
rz(-1.6102487) q[3];
sx q[3];
rz(2.6764398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5588348) q[2];
sx q[2];
rz(-2.1119327) q[2];
sx q[2];
rz(-0.47490698) q[2];
rz(0.20259914) q[3];
sx q[3];
rz(-0.83297268) q[3];
sx q[3];
rz(1.9995662) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1982034) q[0];
sx q[0];
rz(-0.13309637) q[0];
sx q[0];
rz(-2.8357847) q[0];
rz(2.7022779) q[1];
sx q[1];
rz(-0.82249928) q[1];
sx q[1];
rz(1.3154715) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02365774) q[0];
sx q[0];
rz(-2.2542037) q[0];
sx q[0];
rz(1.5562431) q[0];
rz(-pi) q[1];
rz(2.1256746) q[2];
sx q[2];
rz(-1.3279997) q[2];
sx q[2];
rz(1.0958874) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.89375673) q[1];
sx q[1];
rz(-0.79036056) q[1];
sx q[1];
rz(-1.7414581) q[1];
rz(2.1602116) q[3];
sx q[3];
rz(-1.2228612) q[3];
sx q[3];
rz(-0.98435005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4149912) q[2];
sx q[2];
rz(-1.5982268) q[2];
sx q[2];
rz(-0.42465633) q[2];
rz(0.031938227) q[3];
sx q[3];
rz(-1.5141124) q[3];
sx q[3];
rz(2.4922075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.6040819) q[0];
sx q[0];
rz(-2.4896955) q[0];
sx q[0];
rz(0.55484581) q[0];
rz(0.26508731) q[1];
sx q[1];
rz(-1.6731429) q[1];
sx q[1];
rz(1.2010942) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4477312) q[0];
sx q[0];
rz(-0.70708067) q[0];
sx q[0];
rz(-0.62753824) q[0];
rz(-pi) q[1];
rz(-1.7428107) q[2];
sx q[2];
rz(-1.8427249) q[2];
sx q[2];
rz(-3.0957785) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8256067) q[1];
sx q[1];
rz(-0.34994967) q[1];
sx q[1];
rz(3.0341329) q[1];
rz(-pi) q[2];
rz(-2.3917434) q[3];
sx q[3];
rz(-1.5529274) q[3];
sx q[3];
rz(-2.8211879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.889664) q[2];
sx q[2];
rz(-2.1519075) q[2];
sx q[2];
rz(0.84452334) q[2];
rz(-1.1318413) q[3];
sx q[3];
rz(-2.2364538) q[3];
sx q[3];
rz(-1.3593146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1691386) q[0];
sx q[0];
rz(-0.3903946) q[0];
sx q[0];
rz(1.1727232) q[0];
rz(-2.4523465) q[1];
sx q[1];
rz(-1.0968364) q[1];
sx q[1];
rz(-0.26395878) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.496752) q[0];
sx q[0];
rz(-0.94366108) q[0];
sx q[0];
rz(-1.0154614) q[0];
x q[1];
rz(2.2911777) q[2];
sx q[2];
rz(-1.9738832) q[2];
sx q[2];
rz(0.46024761) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.86706272) q[1];
sx q[1];
rz(-1.2559051) q[1];
sx q[1];
rz(0.19537433) q[1];
rz(-pi) q[2];
rz(-3.0194324) q[3];
sx q[3];
rz(-2.2111243) q[3];
sx q[3];
rz(-1.8807008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.76137233) q[2];
sx q[2];
rz(-1.994588) q[2];
sx q[2];
rz(-2.0210361) q[2];
rz(-1.5348966) q[3];
sx q[3];
rz(-0.94529072) q[3];
sx q[3];
rz(1.631261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071038889) q[0];
sx q[0];
rz(-2.4865535) q[0];
sx q[0];
rz(3.0336663) q[0];
rz(0.72036605) q[1];
sx q[1];
rz(-1.9279059) q[1];
sx q[1];
rz(2.7005213) q[1];
rz(-0.82577827) q[2];
sx q[2];
rz(-1.1008395) q[2];
sx q[2];
rz(0.69175082) q[2];
rz(-1.4955487) q[3];
sx q[3];
rz(-2.1863219) q[3];
sx q[3];
rz(-3.0466515) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
