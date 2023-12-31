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
rz(-1.3287969) q[1];
sx q[1];
rz(-1.8741908) q[1];
sx q[1];
rz(1.0277494) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0004814) q[0];
sx q[0];
rz(-1.9465465) q[0];
sx q[0];
rz(1.5953338) q[0];
x q[1];
rz(-2.920354) q[2];
sx q[2];
rz(-2.049963) q[2];
sx q[2];
rz(1.1269119) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7386908) q[1];
sx q[1];
rz(-1.9539023) q[1];
sx q[1];
rz(-2.4163567) q[1];
x q[2];
rz(0.58524744) q[3];
sx q[3];
rz(-2.4260776) q[3];
sx q[3];
rz(-1.1701208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0119005) q[2];
sx q[2];
rz(-1.4346069) q[2];
sx q[2];
rz(1.0502846) q[2];
rz(-2.0283279) q[3];
sx q[3];
rz(-2.2498825) q[3];
sx q[3];
rz(-0.068107001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072409078) q[0];
sx q[0];
rz(-1.8962815) q[0];
sx q[0];
rz(2.8438399) q[0];
rz(0.61966664) q[1];
sx q[1];
rz(-1.0071808) q[1];
sx q[1];
rz(1.108095) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71582687) q[0];
sx q[0];
rz(-0.96164385) q[0];
sx q[0];
rz(-0.25575511) q[0];
rz(-pi) q[1];
rz(2.1790702) q[2];
sx q[2];
rz(-1.2108742) q[2];
sx q[2];
rz(0.72812176) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.121671) q[1];
sx q[1];
rz(-2.6056075) q[1];
sx q[1];
rz(2.3977445) q[1];
x q[2];
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
rz(-pi) q[1];
rz(-2.6790598) q[2];
sx q[2];
rz(-2.1647537) q[2];
sx q[2];
rz(2.1662946) q[2];
rz(-0.9179999) q[3];
sx q[3];
rz(-1.2851597) q[3];
sx q[3];
rz(2.8454034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96238962) q[0];
sx q[0];
rz(-0.85150349) q[0];
sx q[0];
rz(-0.54291022) q[0];
rz(0.88223282) q[1];
sx q[1];
rz(-2.0062607) q[1];
sx q[1];
rz(-0.96484819) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99795139) q[0];
sx q[0];
rz(-0.43295857) q[0];
sx q[0];
rz(-2.4233682) q[0];
x q[1];
rz(-1.6134878) q[2];
sx q[2];
rz(-1.2741158) q[2];
sx q[2];
rz(-1.3670849) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.987625) q[3];
sx q[3];
rz(-1.3891451) q[3];
sx q[3];
rz(-2.3895398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0470011) q[2];
sx q[2];
rz(-0.61085218) q[2];
sx q[2];
rz(1.1331406) q[2];
rz(2.9099693) q[3];
sx q[3];
rz(-1.2730205) q[3];
sx q[3];
rz(-2.384322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
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
rz(-2.8621181) q[0];
sx q[0];
rz(-0.010443895) q[0];
sx q[0];
rz(1.3765155) q[0];
rz(-0.51849413) q[1];
sx q[1];
rz(-1.8771749) q[1];
sx q[1];
rz(-0.24212295) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0644558) q[0];
sx q[0];
rz(-2.4479439) q[0];
sx q[0];
rz(0.25607381) q[0];
x q[1];
rz(2.7119615) q[2];
sx q[2];
rz(-1.5826844) q[2];
sx q[2];
rz(-2.5846543) q[2];
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
x q[2];
rz(-1.711877) q[3];
sx q[3];
rz(-1.1564848) q[3];
sx q[3];
rz(0.29841081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1233998) q[2];
sx q[2];
rz(-2.172956) q[2];
sx q[2];
rz(0.094853178) q[2];
rz(-1.3421966) q[3];
sx q[3];
rz(-1.3972524) q[3];
sx q[3];
rz(-2.7024787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3577394) q[0];
sx q[0];
rz(-1.3107603) q[0];
sx q[0];
rz(-0.36079303) q[0];
rz(1.7533253) q[1];
sx q[1];
rz(-1.3307945) q[1];
sx q[1];
rz(1.1345908) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9103694) q[0];
sx q[0];
rz(-2.0632422) q[0];
sx q[0];
rz(-0.72427303) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7822532) q[2];
sx q[2];
rz(-1.1655072) q[2];
sx q[2];
rz(-2.5422424) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1281631) q[1];
sx q[1];
rz(-1.1268106) q[1];
sx q[1];
rz(-1.7432937) q[1];
x q[2];
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
rz(-2.1363042) q[2];
sx q[2];
rz(-1.4207999) q[2];
sx q[2];
rz(2.5689382) q[2];
rz(2.2128361) q[3];
sx q[3];
rz(-2.6199665) q[3];
sx q[3];
rz(1.1675534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8144433) q[0];
sx q[0];
rz(-2.0454018) q[0];
sx q[0];
rz(-1.2493398) q[0];
rz(-1.918474) q[1];
sx q[1];
rz(-1.616281) q[1];
sx q[1];
rz(-1.1522326) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0056861176) q[0];
sx q[0];
rz(-1.583477) q[0];
sx q[0];
rz(-1.4655051) q[0];
x q[1];
rz(-0.50243369) q[2];
sx q[2];
rz(-1.2049434) q[2];
sx q[2];
rz(0.94787129) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.17105477) q[1];
sx q[1];
rz(-1.4661403) q[1];
sx q[1];
rz(-2.9403951) q[1];
rz(-2.4640718) q[3];
sx q[3];
rz(-2.065425) q[3];
sx q[3];
rz(-2.236931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.46889177) q[2];
sx q[2];
rz(-1.9273309) q[2];
sx q[2];
rz(0.99096283) q[2];
rz(2.4937566) q[3];
sx q[3];
rz(-2.1760553) q[3];
sx q[3];
rz(-2.794054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8836442) q[0];
sx q[0];
rz(-0.22709665) q[0];
sx q[0];
rz(-3.0793072) q[0];
rz(-0.1858055) q[1];
sx q[1];
rz(-1.6848911) q[1];
sx q[1];
rz(-2.7468162) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6900401) q[0];
sx q[0];
rz(-1.8071022) q[0];
sx q[0];
rz(-2.7355746) q[0];
rz(-pi) q[1];
x q[1];
rz(2.765352) q[2];
sx q[2];
rz(-1.2708086) q[2];
sx q[2];
rz(0.30181995) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.77183206) q[1];
sx q[1];
rz(-2.6409915) q[1];
sx q[1];
rz(-0.53497772) q[1];
x q[2];
rz(-1.5338321) q[3];
sx q[3];
rz(-1.4174403) q[3];
sx q[3];
rz(-1.5330029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4930967) q[2];
sx q[2];
rz(-0.33005565) q[2];
sx q[2];
rz(-0.27302343) q[2];
rz(1.3027044) q[3];
sx q[3];
rz(-1.3132934) q[3];
sx q[3];
rz(0.31204143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39772314) q[0];
sx q[0];
rz(-2.0070772) q[0];
sx q[0];
rz(-1.8564818) q[0];
rz(-1.6400736) q[1];
sx q[1];
rz(-1.39095) q[1];
sx q[1];
rz(1.3407019) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9047324) q[0];
sx q[0];
rz(-1.245984) q[0];
sx q[0];
rz(-1.9843285) q[0];
x q[1];
rz(3.0807207) q[2];
sx q[2];
rz(-0.69673046) q[2];
sx q[2];
rz(2.3827162) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8074639) q[1];
sx q[1];
rz(-0.51528105) q[1];
sx q[1];
rz(-0.79828243) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96745456) q[3];
sx q[3];
rz(-2.6199322) q[3];
sx q[3];
rz(1.1286917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1061219) q[2];
sx q[2];
rz(-2.3351228) q[2];
sx q[2];
rz(1.0236053) q[2];
rz(-2.9566531) q[3];
sx q[3];
rz(-0.39026323) q[3];
sx q[3];
rz(-0.24188724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0969365) q[0];
sx q[0];
rz(-2.1450295) q[0];
sx q[0];
rz(-1.5203083) q[0];
rz(0.3301436) q[1];
sx q[1];
rz(-1.2076999) q[1];
sx q[1];
rz(-0.83713371) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5846111) q[0];
sx q[0];
rz(-0.71174445) q[0];
sx q[0];
rz(-0.31194375) q[0];
x q[1];
rz(-1.2741954) q[2];
sx q[2];
rz(-2.1269848) q[2];
sx q[2];
rz(-1.1013168) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7652119) q[1];
sx q[1];
rz(-1.4188671) q[1];
sx q[1];
rz(-3.0357749) q[1];
x q[2];
rz(-1.0650474) q[3];
sx q[3];
rz(-2.5878083) q[3];
sx q[3];
rz(-1.5050642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.60124406) q[2];
sx q[2];
rz(-2.0566172) q[2];
sx q[2];
rz(2.9619651) q[2];
rz(-0.99572292) q[3];
sx q[3];
rz(-1.2487753) q[3];
sx q[3];
rz(1.3109591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76374617) q[0];
sx q[0];
rz(-0.34559956) q[0];
sx q[0];
rz(1.0572222) q[0];
rz(-0.10748848) q[1];
sx q[1];
rz(-1.2534393) q[1];
sx q[1];
rz(-0.9799788) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9419132) q[0];
sx q[0];
rz(-1.4451888) q[0];
sx q[0];
rz(3.0903387) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51771848) q[2];
sx q[2];
rz(-1.1283518) q[2];
sx q[2];
rz(-2.8874318) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.6411297) q[1];
sx q[1];
rz(-2.6424721) q[1];
sx q[1];
rz(-2.0874546) q[1];
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
rz(-pi) q[1];
rz(1.8367299) q[2];
sx q[2];
rz(-1.2048081) q[2];
sx q[2];
rz(-2.1137962) q[2];
rz(1.3868388) q[3];
sx q[3];
rz(-1.3091062) q[3];
sx q[3];
rz(2.8579779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(2.1297395) q[2];
sx q[2];
rz(-1.4197299) q[2];
sx q[2];
rz(1.1001669) q[2];
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
