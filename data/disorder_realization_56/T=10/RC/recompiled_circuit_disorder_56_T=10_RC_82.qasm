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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2078903) q[0];
sx q[0];
rz(-2.7650802) q[0];
sx q[0];
rz(-3.0794789) q[0];
rz(-0.22123863) q[2];
sx q[2];
rz(-1.0916296) q[2];
sx q[2];
rz(1.1269119) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6537522) q[1];
sx q[1];
rz(-0.90812212) q[1];
sx q[1];
rz(2.0648048) q[1];
rz(2.5146033) q[3];
sx q[3];
rz(-1.1999745) q[3];
sx q[3];
rz(-2.2771319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1296922) q[2];
sx q[2];
rz(-1.4346069) q[2];
sx q[2];
rz(2.091308) q[2];
rz(1.1132647) q[3];
sx q[3];
rz(-0.89171019) q[3];
sx q[3];
rz(-3.0734857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71582687) q[0];
sx q[0];
rz(-2.1799488) q[0];
sx q[0];
rz(-0.25575511) q[0];
rz(-pi) q[1];
rz(-0.42995288) q[2];
sx q[2];
rz(-1.0064831) q[2];
sx q[2];
rz(-1.0831837) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.88156466) q[1];
sx q[1];
rz(-1.9238872) q[1];
sx q[1];
rz(0.41207037) q[1];
x q[2];
rz(-0.31538972) q[3];
sx q[3];
rz(-2.4139801) q[3];
sx q[3];
rz(1.417158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.46253282) q[2];
sx q[2];
rz(-2.1647537) q[2];
sx q[2];
rz(2.1662946) q[2];
rz(-0.9179999) q[3];
sx q[3];
rz(-1.8564329) q[3];
sx q[3];
rz(0.29618922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96238962) q[0];
sx q[0];
rz(-2.2900892) q[0];
sx q[0];
rz(0.54291022) q[0];
rz(2.2593598) q[1];
sx q[1];
rz(-2.0062607) q[1];
sx q[1];
rz(-2.1767445) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99795139) q[0];
sx q[0];
rz(-0.43295857) q[0];
sx q[0];
rz(0.71822449) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6134878) q[2];
sx q[2];
rz(-1.2741158) q[2];
sx q[2];
rz(-1.7745078) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1988586) q[1];
sx q[1];
rz(-0.13730362) q[1];
sx q[1];
rz(1.0820461) q[1];
rz(-pi) q[2];
rz(0.19823719) q[3];
sx q[3];
rz(-1.9803515) q[3];
sx q[3];
rz(0.73892456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.09459153) q[2];
sx q[2];
rz(-0.61085218) q[2];
sx q[2];
rz(-2.0084521) q[2];
rz(0.23162332) q[3];
sx q[3];
rz(-1.8685721) q[3];
sx q[3];
rz(0.75727063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8621181) q[0];
sx q[0];
rz(-0.010443895) q[0];
sx q[0];
rz(-1.7650771) q[0];
rz(-0.51849413) q[1];
sx q[1];
rz(-1.8771749) q[1];
sx q[1];
rz(2.8994697) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3926485) q[0];
sx q[0];
rz(-2.2376275) q[0];
sx q[0];
rz(-1.7783661) q[0];
rz(-2.7119615) q[2];
sx q[2];
rz(-1.5589082) q[2];
sx q[2];
rz(0.55693835) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.33430797) q[1];
sx q[1];
rz(-2.0172999) q[1];
sx q[1];
rz(1.7590894) q[1];
x q[2];
rz(-0.30947134) q[3];
sx q[3];
rz(-0.43635338) q[3];
sx q[3];
rz(2.5040124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.1233998) q[2];
sx q[2];
rz(-0.96863666) q[2];
sx q[2];
rz(-3.0467395) q[2];
rz(1.3421966) q[3];
sx q[3];
rz(-1.3972524) q[3];
sx q[3];
rz(-0.43911394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3577394) q[0];
sx q[0];
rz(-1.3107603) q[0];
sx q[0];
rz(-0.36079303) q[0];
rz(1.3882673) q[1];
sx q[1];
rz(-1.8107982) q[1];
sx q[1];
rz(1.1345908) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9103694) q[0];
sx q[0];
rz(-2.0632422) q[0];
sx q[0];
rz(2.4173196) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35933944) q[2];
sx q[2];
rz(-1.9760855) q[2];
sx q[2];
rz(0.59935024) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1281631) q[1];
sx q[1];
rz(-2.014782) q[1];
sx q[1];
rz(-1.3982989) q[1];
x q[2];
rz(-1.4868823) q[3];
sx q[3];
rz(-1.1353555) q[3];
sx q[3];
rz(2.8403789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0052884) q[2];
sx q[2];
rz(-1.4207999) q[2];
sx q[2];
rz(2.5689382) q[2];
rz(-0.92875656) q[3];
sx q[3];
rz(-0.52162617) q[3];
sx q[3];
rz(-1.1675534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8144433) q[0];
sx q[0];
rz(-1.0961908) q[0];
sx q[0];
rz(-1.2493398) q[0];
rz(-1.2231187) q[1];
sx q[1];
rz(-1.616281) q[1];
sx q[1];
rz(1.1522326) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0056861176) q[0];
sx q[0];
rz(-1.583477) q[0];
sx q[0];
rz(-1.6760875) q[0];
x q[1];
rz(1.9828898) q[2];
sx q[2];
rz(-1.1043784) q[2];
sx q[2];
rz(-0.42883021) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2684979) q[1];
sx q[1];
rz(-0.22646204) q[1];
sx q[1];
rz(-0.48392673) q[1];
x q[2];
rz(-0.67752083) q[3];
sx q[3];
rz(-1.0761677) q[3];
sx q[3];
rz(-2.236931) q[3];
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
rz(-0.99096283) q[2];
rz(2.4937566) q[3];
sx q[3];
rz(-2.1760553) q[3];
sx q[3];
rz(0.34753862) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8836442) q[0];
sx q[0];
rz(-2.914496) q[0];
sx q[0];
rz(3.0793072) q[0];
rz(2.9557872) q[1];
sx q[1];
rz(-1.6848911) q[1];
sx q[1];
rz(-2.7468162) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5236854) q[0];
sx q[0];
rz(-0.46645188) q[0];
sx q[0];
rz(0.54752366) q[0];
rz(-pi) q[1];
rz(1.2497181) q[2];
sx q[2];
rz(-1.9294538) q[2];
sx q[2];
rz(1.9888339) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9638622) q[1];
sx q[1];
rz(-1.1451671) q[1];
sx q[1];
rz(1.8427986) q[1];
rz(-pi) q[2];
rz(-2.9069101) q[3];
sx q[3];
rz(-0.15771401) q[3];
sx q[3];
rz(1.7705256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.64849598) q[2];
sx q[2];
rz(-2.811537) q[2];
sx q[2];
rz(-2.8685692) q[2];
rz(-1.3027044) q[3];
sx q[3];
rz(-1.8282993) q[3];
sx q[3];
rz(-2.8295512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7438695) q[0];
sx q[0];
rz(-1.1345154) q[0];
sx q[0];
rz(1.8564818) q[0];
rz(-1.6400736) q[1];
sx q[1];
rz(-1.39095) q[1];
sx q[1];
rz(-1.8008908) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2368603) q[0];
sx q[0];
rz(-1.245984) q[0];
sx q[0];
rz(-1.1572641) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4457744) q[2];
sx q[2];
rz(-1.5317481) q[2];
sx q[2];
rz(-2.2829636) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2015842) q[1];
sx q[1];
rz(-1.219698) q[1];
sx q[1];
rz(1.9560948) q[1];
rz(2.826346) q[3];
sx q[3];
rz(-1.1479706) q[3];
sx q[3];
rz(-1.8002312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1061219) q[2];
sx q[2];
rz(-0.80646986) q[2];
sx q[2];
rz(-1.0236053) q[2];
rz(-0.18493955) q[3];
sx q[3];
rz(-0.39026323) q[3];
sx q[3];
rz(-2.8997054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0446562) q[0];
sx q[0];
rz(-2.1450295) q[0];
sx q[0];
rz(1.5203083) q[0];
rz(2.8114491) q[1];
sx q[1];
rz(-1.2076999) q[1];
sx q[1];
rz(0.83713371) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22568208) q[0];
sx q[0];
rz(-1.7726232) q[0];
sx q[0];
rz(0.68737824) q[0];
rz(-pi) q[1];
rz(-1.2741954) q[2];
sx q[2];
rz(-2.1269848) q[2];
sx q[2];
rz(-1.1013168) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9311034) q[1];
sx q[1];
rz(-1.6753907) q[1];
sx q[1];
rz(-1.4180257) q[1];
x q[2];
rz(1.0650474) q[3];
sx q[3];
rz(-0.55378434) q[3];
sx q[3];
rz(-1.5050642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5403486) q[2];
sx q[2];
rz(-2.0566172) q[2];
sx q[2];
rz(-0.17962757) q[2];
rz(-0.99572292) q[3];
sx q[3];
rz(-1.2487753) q[3];
sx q[3];
rz(1.3109591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3778465) q[0];
sx q[0];
rz(-2.7959931) q[0];
sx q[0];
rz(-1.0572222) q[0];
rz(0.10748848) q[1];
sx q[1];
rz(-1.2534393) q[1];
sx q[1];
rz(0.9799788) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5532593) q[0];
sx q[0];
rz(-3.0059814) q[0];
sx q[0];
rz(-1.9562264) q[0];
rz(-pi) q[1];
rz(-2.0699632) q[2];
sx q[2];
rz(-2.0344779) q[2];
sx q[2];
rz(-1.0774563) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6411297) q[1];
sx q[1];
rz(-0.49912057) q[1];
sx q[1];
rz(2.0874546) q[1];
rz(-pi) q[2];
rz(1.8944593) q[3];
sx q[3];
rz(-1.666288) q[3];
sx q[3];
rz(-3.132706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8367299) q[2];
sx q[2];
rz(-1.9367846) q[2];
sx q[2];
rz(-1.0277964) q[2];
rz(1.3868388) q[3];
sx q[3];
rz(-1.8324865) q[3];
sx q[3];
rz(-2.8579779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86823157) q[0];
sx q[0];
rz(-2.1049451) q[0];
sx q[0];
rz(2.0275397) q[0];
rz(-0.83203075) q[1];
sx q[1];
rz(-2.6770626) q[1];
sx q[1];
rz(-2.4774036) q[1];
rz(-2.9639347) q[2];
sx q[2];
rz(-2.1226317) q[2];
sx q[2];
rz(2.5771099) q[2];
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
