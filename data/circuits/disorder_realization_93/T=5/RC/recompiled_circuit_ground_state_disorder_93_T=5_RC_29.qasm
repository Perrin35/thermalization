OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.326958) q[0];
sx q[0];
rz(-3.1231413) q[0];
sx q[0];
rz(2.9676394) q[0];
rz(1.8040909) q[1];
sx q[1];
rz(4.017158) q[1];
sx q[1];
rz(12.234966) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1179292) q[0];
sx q[0];
rz(-2.8011596) q[0];
sx q[0];
rz(-2.6871919) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.69281225) q[2];
sx q[2];
rz(-1.1651784) q[2];
sx q[2];
rz(-1.931153) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1494721) q[1];
sx q[1];
rz(-1.3821342) q[1];
sx q[1];
rz(-0.68575286) q[1];
rz(-pi) q[2];
rz(2.9339497) q[3];
sx q[3];
rz(-0.15842008) q[3];
sx q[3];
rz(1.0760361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6203674) q[2];
sx q[2];
rz(-1.2646893) q[2];
sx q[2];
rz(0.59086665) q[2];
rz(1.2938195) q[3];
sx q[3];
rz(-1.3760682) q[3];
sx q[3];
rz(2.6025492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6227601) q[0];
sx q[0];
rz(-2.9157214) q[0];
sx q[0];
rz(-2.2665005) q[0];
rz(-0.092763364) q[1];
sx q[1];
rz(-1.0567254) q[1];
sx q[1];
rz(-1.0088395) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0419249) q[0];
sx q[0];
rz(-2.0001051) q[0];
sx q[0];
rz(-0.03089182) q[0];
x q[1];
rz(1.5968199) q[2];
sx q[2];
rz(-0.98434005) q[2];
sx q[2];
rz(1.1555188) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.61774594) q[1];
sx q[1];
rz(-2.2154795) q[1];
sx q[1];
rz(-2.9545386) q[1];
x q[2];
rz(1.5126692) q[3];
sx q[3];
rz(-1.584313) q[3];
sx q[3];
rz(2.0292387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2729317) q[2];
sx q[2];
rz(-1.2596143) q[2];
sx q[2];
rz(1.3966365) q[2];
rz(-0.40846387) q[3];
sx q[3];
rz(-2.4558081) q[3];
sx q[3];
rz(0.43533152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8946266) q[0];
sx q[0];
rz(-0.28852391) q[0];
sx q[0];
rz(-2.1911279) q[0];
rz(1.2735927) q[1];
sx q[1];
rz(-1.344695) q[1];
sx q[1];
rz(-1.7705932) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4038611) q[0];
sx q[0];
rz(-1.4329684) q[0];
sx q[0];
rz(0.65562427) q[0];
rz(-pi) q[1];
rz(0.43355219) q[2];
sx q[2];
rz(-1.6759652) q[2];
sx q[2];
rz(2.1702332) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.85670722) q[1];
sx q[1];
rz(-2.6213985) q[1];
sx q[1];
rz(-1.8265282) q[1];
rz(-pi) q[2];
rz(-2.5018243) q[3];
sx q[3];
rz(-1.4832212) q[3];
sx q[3];
rz(-1.4882907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3889968) q[2];
sx q[2];
rz(-1.7277371) q[2];
sx q[2];
rz(-1.5163806) q[2];
rz(-0.73836941) q[3];
sx q[3];
rz(-1.5117896) q[3];
sx q[3];
rz(0.71877688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3589288) q[0];
sx q[0];
rz(-1.6343225) q[0];
sx q[0];
rz(-1.2231109) q[0];
rz(0.75543985) q[1];
sx q[1];
rz(-2.0016045) q[1];
sx q[1];
rz(3.0208407) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31407315) q[0];
sx q[0];
rz(-1.4085007) q[0];
sx q[0];
rz(-1.7364962) q[0];
x q[1];
rz(-2.7677972) q[2];
sx q[2];
rz(-2.5540562) q[2];
sx q[2];
rz(-2.2898506) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.11343458) q[1];
sx q[1];
rz(-2.096281) q[1];
sx q[1];
rz(-2.3281086) q[1];
rz(-pi) q[2];
rz(1.6183245) q[3];
sx q[3];
rz(-1.5521677) q[3];
sx q[3];
rz(-0.19028529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6197023) q[2];
sx q[2];
rz(-1.686124) q[2];
sx q[2];
rz(1.5187368) q[2];
rz(-1.8146727) q[3];
sx q[3];
rz(-0.36553317) q[3];
sx q[3];
rz(-2.0324619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-3.0066852) q[0];
sx q[0];
rz(-1.4445855) q[0];
sx q[0];
rz(-0.69394845) q[0];
rz(-0.42849439) q[1];
sx q[1];
rz(-0.65281147) q[1];
sx q[1];
rz(1.8983715) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4994298) q[0];
sx q[0];
rz(-1.3754396) q[0];
sx q[0];
rz(0.046577297) q[0];
rz(-1.7967729) q[2];
sx q[2];
rz(-1.4046841) q[2];
sx q[2];
rz(1.0122344) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8163838) q[1];
sx q[1];
rz(-1.7933791) q[1];
sx q[1];
rz(1.9680262) q[1];
rz(-2.6596131) q[3];
sx q[3];
rz(-1.478412) q[3];
sx q[3];
rz(0.57541144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1474233) q[2];
sx q[2];
rz(-1.3012412) q[2];
sx q[2];
rz(-2.6541138) q[2];
rz(0.068988919) q[3];
sx q[3];
rz(-2.4853849) q[3];
sx q[3];
rz(0.60855734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2212806) q[0];
sx q[0];
rz(-2.1423036) q[0];
sx q[0];
rz(0.47166237) q[0];
rz(-2.3645875) q[1];
sx q[1];
rz(-1.124137) q[1];
sx q[1];
rz(-1.5583001) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7651833) q[0];
sx q[0];
rz(-1.6003971) q[0];
sx q[0];
rz(-1.5628003) q[0];
rz(-pi) q[1];
rz(-2.7194994) q[2];
sx q[2];
rz(-2.0712086) q[2];
sx q[2];
rz(-1.3304324) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.029249749) q[1];
sx q[1];
rz(-1.7866553) q[1];
sx q[1];
rz(2.6964746) q[1];
rz(-pi) q[2];
rz(2.4151012) q[3];
sx q[3];
rz(-1.3950751) q[3];
sx q[3];
rz(-0.54783152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.16288699) q[2];
sx q[2];
rz(-1.6521613) q[2];
sx q[2];
rz(1.9784652) q[2];
rz(-2.2833521) q[3];
sx q[3];
rz(-1.4505998) q[3];
sx q[3];
rz(0.8391909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6100886) q[0];
sx q[0];
rz(-1.6781582) q[0];
sx q[0];
rz(0.32291821) q[0];
rz(-1.8220176) q[1];
sx q[1];
rz(-2.0009305) q[1];
sx q[1];
rz(-1.742977) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6933029) q[0];
sx q[0];
rz(-2.9041661) q[0];
sx q[0];
rz(-2.3795932) q[0];
rz(-0.9594907) q[2];
sx q[2];
rz(-1.9027019) q[2];
sx q[2];
rz(-1.6417981) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3460878) q[1];
sx q[1];
rz(-0.75199703) q[1];
sx q[1];
rz(3.094207) q[1];
rz(-1.6635464) q[3];
sx q[3];
rz(-1.4418915) q[3];
sx q[3];
rz(3.1069251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4601595) q[2];
sx q[2];
rz(-2.5761649) q[2];
sx q[2];
rz(2.4088755) q[2];
rz(-0.021746246) q[3];
sx q[3];
rz(-1.5214058) q[3];
sx q[3];
rz(-0.17797962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9916423) q[0];
sx q[0];
rz(-3.0198779) q[0];
sx q[0];
rz(0.44418401) q[0];
rz(1.3581879) q[1];
sx q[1];
rz(-1.5779747) q[1];
sx q[1];
rz(-1.1202728) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5356784) q[0];
sx q[0];
rz(-1.4679424) q[0];
sx q[0];
rz(0.82983183) q[0];
rz(-0.6702126) q[2];
sx q[2];
rz(-0.59500098) q[2];
sx q[2];
rz(-0.52192823) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7071735) q[1];
sx q[1];
rz(-1.7203693) q[1];
sx q[1];
rz(-0.24874462) q[1];
rz(1.9838263) q[3];
sx q[3];
rz(-1.1404697) q[3];
sx q[3];
rz(2.4347507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.292395) q[2];
sx q[2];
rz(-2.2870543) q[2];
sx q[2];
rz(1.9067859) q[2];
rz(-2.1644927) q[3];
sx q[3];
rz(-1.2950725) q[3];
sx q[3];
rz(0.60404122) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4567756) q[0];
sx q[0];
rz(-0.026739459) q[0];
sx q[0];
rz(-0.34715125) q[0];
rz(-2.8390362) q[1];
sx q[1];
rz(-2.0358993) q[1];
sx q[1];
rz(0.0060630719) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5232698) q[0];
sx q[0];
rz(-1.5875641) q[0];
sx q[0];
rz(-1.3888593) q[0];
rz(-pi) q[1];
rz(-0.89170189) q[2];
sx q[2];
rz(-1.804934) q[2];
sx q[2];
rz(-0.32227031) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4681983) q[1];
sx q[1];
rz(-1.2691193) q[1];
sx q[1];
rz(-2.2190898) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.047386332) q[3];
sx q[3];
rz(-2.3607254) q[3];
sx q[3];
rz(0.56347195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2431295) q[2];
sx q[2];
rz(-2.1738238) q[2];
sx q[2];
rz(1.751162) q[2];
rz(9/(1*pi)) q[3];
sx q[3];
rz(-1.0172458) q[3];
sx q[3];
rz(1.0677342) q[3];
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
rz(-1.4047563) q[0];
sx q[0];
rz(-2.3911349) q[0];
sx q[0];
rz(-2.2682037) q[0];
rz(-2.5923173) q[1];
sx q[1];
rz(-2.4875689) q[1];
sx q[1];
rz(-2.3347847) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1081558) q[0];
sx q[0];
rz(-1.8315541) q[0];
sx q[0];
rz(-0.64020641) q[0];
rz(0.64555706) q[2];
sx q[2];
rz(-0.19536138) q[2];
sx q[2];
rz(-0.33249172) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5196073) q[1];
sx q[1];
rz(-0.27493536) q[1];
sx q[1];
rz(1.1311929) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96050519) q[3];
sx q[3];
rz(-1.3209855) q[3];
sx q[3];
rz(-0.042674289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.4775382) q[2];
sx q[2];
rz(-1.8287903) q[2];
sx q[2];
rz(0.99267268) q[2];
rz(-0.70332518) q[3];
sx q[3];
rz(-1.1329009) q[3];
sx q[3];
rz(2.4542184) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0032046) q[0];
sx q[0];
rz(-0.78835543) q[0];
sx q[0];
rz(2.0302571) q[0];
rz(-2.5913024) q[1];
sx q[1];
rz(-2.7688409) q[1];
sx q[1];
rz(0.54584835) q[1];
rz(-2.4443378) q[2];
sx q[2];
rz(-1.6959126) q[2];
sx q[2];
rz(-0.87903862) q[2];
rz(-1.1543858) q[3];
sx q[3];
rz(-1.1416804) q[3];
sx q[3];
rz(0.46905372) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
