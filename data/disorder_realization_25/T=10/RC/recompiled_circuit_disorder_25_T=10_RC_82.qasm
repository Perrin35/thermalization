OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0796354) q[0];
sx q[0];
rz(-2.893653) q[0];
sx q[0];
rz(0.82013446) q[0];
rz(0.040854383) q[1];
sx q[1];
rz(-2.3526469) q[1];
sx q[1];
rz(0.087892428) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.983842) q[0];
sx q[0];
rz(-2.0294667) q[0];
sx q[0];
rz(1.0755324) q[0];
x q[1];
rz(2.0983661) q[2];
sx q[2];
rz(-1.1967778) q[2];
sx q[2];
rz(2.5803215) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76074782) q[1];
sx q[1];
rz(-1.5686085) q[1];
sx q[1];
rz(2.6373177) q[1];
rz(-pi) q[2];
rz(-1.8294931) q[3];
sx q[3];
rz(-1.5764578) q[3];
sx q[3];
rz(-0.86710801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0007881) q[2];
sx q[2];
rz(-1.6892097) q[2];
sx q[2];
rz(-2.5722356) q[2];
rz(1.6128444) q[3];
sx q[3];
rz(-0.6362392) q[3];
sx q[3];
rz(1.7830085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59812087) q[0];
sx q[0];
rz(-2.9765029) q[0];
sx q[0];
rz(-2.5879481) q[0];
rz(1.9042632) q[1];
sx q[1];
rz(-1.3668704) q[1];
sx q[1];
rz(-1.9083317) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90399088) q[0];
sx q[0];
rz(-1.0407789) q[0];
sx q[0];
rz(0.7941829) q[0];
x q[1];
rz(-1.6987259) q[2];
sx q[2];
rz(-2.2413841) q[2];
sx q[2];
rz(0.48534457) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1439647) q[1];
sx q[1];
rz(-1.3647172) q[1];
sx q[1];
rz(-3.1254083) q[1];
rz(-pi) q[2];
rz(-1.3319098) q[3];
sx q[3];
rz(-0.91730648) q[3];
sx q[3];
rz(1.4373923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1374986) q[2];
sx q[2];
rz(-1.4592905) q[2];
sx q[2];
rz(3.1393576) q[2];
rz(0.8301174) q[3];
sx q[3];
rz(-2.3486962) q[3];
sx q[3];
rz(1.2494276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24901351) q[0];
sx q[0];
rz(-0.61616388) q[0];
sx q[0];
rz(-2.2900443) q[0];
rz(0.7710723) q[1];
sx q[1];
rz(-1.675019) q[1];
sx q[1];
rz(-0.071203701) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.486858) q[0];
sx q[0];
rz(-1.3329957) q[0];
sx q[0];
rz(-0.20240692) q[0];
x q[1];
rz(2.879197) q[2];
sx q[2];
rz(-2.0121644) q[2];
sx q[2];
rz(1.0052094) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.066312) q[1];
sx q[1];
rz(-1.4367141) q[1];
sx q[1];
rz(-2.4877497) q[1];
rz(-pi) q[2];
rz(-0.72910587) q[3];
sx q[3];
rz(-2.2359072) q[3];
sx q[3];
rz(-0.8641181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8292134) q[2];
sx q[2];
rz(-0.92386121) q[2];
sx q[2];
rz(2.6181347) q[2];
rz(0.33189014) q[3];
sx q[3];
rz(-0.060398014) q[3];
sx q[3];
rz(-0.50326842) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4098542) q[0];
sx q[0];
rz(-0.24895746) q[0];
sx q[0];
rz(0.82558924) q[0];
rz(-1.1666974) q[1];
sx q[1];
rz(-2.2749133) q[1];
sx q[1];
rz(1.4473787) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81340504) q[0];
sx q[0];
rz(-1.4077328) q[0];
sx q[0];
rz(-1.8036611) q[0];
rz(1.3924613) q[2];
sx q[2];
rz(-1.5618088) q[2];
sx q[2];
rz(1.0967147) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.96781603) q[1];
sx q[1];
rz(-0.88138352) q[1];
sx q[1];
rz(-2.9348228) q[1];
rz(-pi) q[2];
rz(-0.57595423) q[3];
sx q[3];
rz(-0.42612694) q[3];
sx q[3];
rz(-1.061561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.46431413) q[2];
sx q[2];
rz(-1.364664) q[2];
sx q[2];
rz(0.14492598) q[2];
rz(-1.011301) q[3];
sx q[3];
rz(-2.5144808) q[3];
sx q[3];
rz(-3.1269126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.1458364) q[0];
sx q[0];
rz(-0.053481426) q[0];
sx q[0];
rz(-0.68403912) q[0];
rz(-1.1902635) q[1];
sx q[1];
rz(-0.66527706) q[1];
sx q[1];
rz(-1.2129983) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4764403) q[0];
sx q[0];
rz(-1.9422918) q[0];
sx q[0];
rz(1.8508338) q[0];
rz(-pi) q[1];
rz(2.5118675) q[2];
sx q[2];
rz(-2.1246398) q[2];
sx q[2];
rz(-2.6385006) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.1274384) q[1];
sx q[1];
rz(-1.452983) q[1];
sx q[1];
rz(0.25728667) q[1];
rz(-2.1608859) q[3];
sx q[3];
rz(-0.31168918) q[3];
sx q[3];
rz(2.2463472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6325536) q[2];
sx q[2];
rz(-1.8703439) q[2];
sx q[2];
rz(0.63009134) q[2];
rz(-0.57224327) q[3];
sx q[3];
rz(-0.64544353) q[3];
sx q[3];
rz(0.8852638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0260139) q[0];
sx q[0];
rz(-0.84340874) q[0];
sx q[0];
rz(2.8552326) q[0];
rz(2.5419366) q[1];
sx q[1];
rz(-1.1279794) q[1];
sx q[1];
rz(-2.5568331) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50305784) q[0];
sx q[0];
rz(-0.98031822) q[0];
sx q[0];
rz(-1.9496586) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5969959) q[2];
sx q[2];
rz(-0.95653406) q[2];
sx q[2];
rz(2.7363077) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6905745) q[1];
sx q[1];
rz(-1.283053) q[1];
sx q[1];
rz(2.2698127) q[1];
rz(-pi) q[2];
rz(-1.6213958) q[3];
sx q[3];
rz(-1.7369324) q[3];
sx q[3];
rz(1.43464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.16453234) q[2];
sx q[2];
rz(-0.88447276) q[2];
sx q[2];
rz(1.5131081) q[2];
rz(-2.1785054) q[3];
sx q[3];
rz(-1.3685127) q[3];
sx q[3];
rz(-0.62455463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088783711) q[0];
sx q[0];
rz(-1.7099021) q[0];
sx q[0];
rz(0.54876304) q[0];
rz(-0.051963003) q[1];
sx q[1];
rz(-2.6497662) q[1];
sx q[1];
rz(0.60639492) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43019766) q[0];
sx q[0];
rz(-1.0476307) q[0];
sx q[0];
rz(-2.8682312) q[0];
x q[1];
rz(2.7658471) q[2];
sx q[2];
rz(-2.1990015) q[2];
sx q[2];
rz(1.9919765) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7402621) q[1];
sx q[1];
rz(-1.6832502) q[1];
sx q[1];
rz(-2.3349891) q[1];
rz(-2.2783979) q[3];
sx q[3];
rz(-0.51897012) q[3];
sx q[3];
rz(-0.88245813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2018044) q[2];
sx q[2];
rz(-2.1494892) q[2];
sx q[2];
rz(-0.50841224) q[2];
rz(-1.1172179) q[3];
sx q[3];
rz(-2.9682187) q[3];
sx q[3];
rz(-1.4846876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33928076) q[0];
sx q[0];
rz(-2.7957714) q[0];
sx q[0];
rz(0.75396496) q[0];
rz(2.1915961) q[1];
sx q[1];
rz(-1.6597304) q[1];
sx q[1];
rz(-2.9139013) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3823711) q[0];
sx q[0];
rz(-1.4996487) q[0];
sx q[0];
rz(-1.1573777) q[0];
x q[1];
rz(2.4661602) q[2];
sx q[2];
rz(-2.5917412) q[2];
sx q[2];
rz(0.39381248) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.88340064) q[1];
sx q[1];
rz(-1.0419783) q[1];
sx q[1];
rz(-1.8316815) q[1];
rz(-2.5826107) q[3];
sx q[3];
rz(-1.6056656) q[3];
sx q[3];
rz(-2.1502286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.078538744) q[2];
sx q[2];
rz(-2.9324014) q[2];
sx q[2];
rz(0.83855808) q[2];
rz(-0.36455425) q[3];
sx q[3];
rz(-1.3804599) q[3];
sx q[3];
rz(-0.84028876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39695981) q[0];
sx q[0];
rz(-1.7663706) q[0];
sx q[0];
rz(-2.3378085) q[0];
rz(-1.0546168) q[1];
sx q[1];
rz(-0.63916731) q[1];
sx q[1];
rz(-1.1484336) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.046612) q[0];
sx q[0];
rz(-1.5469095) q[0];
sx q[0];
rz(-1.6427965) q[0];
rz(-pi) q[1];
rz(2.7788413) q[2];
sx q[2];
rz(-2.0667549) q[2];
sx q[2];
rz(-2.9438058) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.941274) q[1];
sx q[1];
rz(-0.66220821) q[1];
sx q[1];
rz(2.3108285) q[1];
rz(-pi) q[2];
rz(3.0461379) q[3];
sx q[3];
rz(-0.3534375) q[3];
sx q[3];
rz(-2.8638338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.049872963) q[2];
sx q[2];
rz(-2.1240978) q[2];
sx q[2];
rz(-2.3633374) q[2];
rz(-0.19566472) q[3];
sx q[3];
rz(-1.0396495) q[3];
sx q[3];
rz(2.1197317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.8675999) q[0];
sx q[0];
rz(-2.1430528) q[0];
sx q[0];
rz(-0.25319779) q[0];
rz(-0.7397488) q[1];
sx q[1];
rz(-2.4052129) q[1];
sx q[1];
rz(0.69828066) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65933278) q[0];
sx q[0];
rz(-1.2685304) q[0];
sx q[0];
rz(1.4955826) q[0];
rz(-3.1168078) q[2];
sx q[2];
rz(-1.6449271) q[2];
sx q[2];
rz(-0.63400808) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0638949) q[1];
sx q[1];
rz(-0.5118891) q[1];
sx q[1];
rz(-2.0236532) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24153696) q[3];
sx q[3];
rz(-2.1250484) q[3];
sx q[3];
rz(0.17061558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.11761052) q[2];
sx q[2];
rz(-2.2587946) q[2];
sx q[2];
rz(2.3013766) q[2];
rz(-0.64030567) q[3];
sx q[3];
rz(-0.87602031) q[3];
sx q[3];
rz(1.1673814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8828076) q[0];
sx q[0];
rz(-2.2866645) q[0];
sx q[0];
rz(-1.4186161) q[0];
rz(-0.029126833) q[1];
sx q[1];
rz(-3.0283785) q[1];
sx q[1];
rz(1.821847) q[1];
rz(2.6559033) q[2];
sx q[2];
rz(-2.4808241) q[2];
sx q[2];
rz(-0.11447699) q[2];
rz(2.9914231) q[3];
sx q[3];
rz(-0.88610813) q[3];
sx q[3];
rz(-0.025972493) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];