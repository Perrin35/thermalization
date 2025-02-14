OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5950522) q[0];
sx q[0];
rz(-0.15931436) q[0];
sx q[0];
rz(1.2567318) q[0];
rz(-2.9873084) q[1];
sx q[1];
rz(-1.0839387) q[1];
sx q[1];
rz(-1.3561603) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6124466) q[0];
sx q[0];
rz(-2.6137798) q[0];
sx q[0];
rz(2.0291269) q[0];
x q[1];
rz(-2.9262395) q[2];
sx q[2];
rz(-2.9850508) q[2];
sx q[2];
rz(-2.6555344) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.534733) q[1];
sx q[1];
rz(-2.4280425) q[1];
sx q[1];
rz(-3.0487398) q[1];
rz(1.1514444) q[3];
sx q[3];
rz(-2.0128801) q[3];
sx q[3];
rz(-2.6579554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1821182) q[2];
sx q[2];
rz(-1.139816) q[2];
sx q[2];
rz(-3.0457022) q[2];
rz(0.31283665) q[3];
sx q[3];
rz(-2.0818807) q[3];
sx q[3];
rz(-2.634341) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7392015) q[0];
sx q[0];
rz(-1.6283789) q[0];
sx q[0];
rz(1.2865404) q[0];
rz(-1.3013499) q[1];
sx q[1];
rz(-1.2570612) q[1];
sx q[1];
rz(1.4305065) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1576618) q[0];
sx q[0];
rz(-1.0583911) q[0];
sx q[0];
rz(-0.061467193) q[0];
rz(-pi) q[1];
rz(-0.78708055) q[2];
sx q[2];
rz(-2.5445017) q[2];
sx q[2];
rz(0.40368375) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7315361) q[1];
sx q[1];
rz(-1.6997196) q[1];
sx q[1];
rz(1.6881264) q[1];
rz(-pi) q[2];
rz(-0.28263553) q[3];
sx q[3];
rz(-2.4135597) q[3];
sx q[3];
rz(-1.1580434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8282738) q[2];
sx q[2];
rz(-0.88130772) q[2];
sx q[2];
rz(-2.4271097) q[2];
rz(2.1099527) q[3];
sx q[3];
rz(-2.1561626) q[3];
sx q[3];
rz(-0.45903444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76630509) q[0];
sx q[0];
rz(-2.3630688) q[0];
sx q[0];
rz(-2.9189723) q[0];
rz(-2.8011232) q[1];
sx q[1];
rz(-1.6727996) q[1];
sx q[1];
rz(-0.31390831) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1134046) q[0];
sx q[0];
rz(-2.2566099) q[0];
sx q[0];
rz(-0.96341204) q[0];
rz(2.3000175) q[2];
sx q[2];
rz(-1.1969337) q[2];
sx q[2];
rz(-0.29619994) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7063278) q[1];
sx q[1];
rz(-1.2641331) q[1];
sx q[1];
rz(-2.4198662) q[1];
rz(0.49626428) q[3];
sx q[3];
rz(-1.7682425) q[3];
sx q[3];
rz(2.4112301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8935304) q[2];
sx q[2];
rz(-1.6174822) q[2];
sx q[2];
rz(2.1171872) q[2];
rz(-1.3192568) q[3];
sx q[3];
rz(-0.28194591) q[3];
sx q[3];
rz(-0.2969186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1612741) q[0];
sx q[0];
rz(-1.4264822) q[0];
sx q[0];
rz(1.0262161) q[0];
rz(-0.16464344) q[1];
sx q[1];
rz(-2.7963729) q[1];
sx q[1];
rz(-0.77273291) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29870009) q[0];
sx q[0];
rz(-0.94988454) q[0];
sx q[0];
rz(-2.0218532) q[0];
rz(0.010920694) q[2];
sx q[2];
rz(-0.35775634) q[2];
sx q[2];
rz(-0.27586473) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.31269909) q[1];
sx q[1];
rz(-0.71694198) q[1];
sx q[1];
rz(-1.904105) q[1];
rz(-pi) q[2];
rz(-1.5192612) q[3];
sx q[3];
rz(-1.5750575) q[3];
sx q[3];
rz(-0.86458998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2356977) q[2];
sx q[2];
rz(-0.57196456) q[2];
sx q[2];
rz(2.8254361) q[2];
rz(-2.9302127) q[3];
sx q[3];
rz(-1.4485161) q[3];
sx q[3];
rz(-2.0785418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0036156) q[0];
sx q[0];
rz(-1.5253541) q[0];
sx q[0];
rz(0.77950087) q[0];
rz(-1.4315073) q[1];
sx q[1];
rz(-1.5305488) q[1];
sx q[1];
rz(-0.12758189) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0128827) q[0];
sx q[0];
rz(-1.5361551) q[0];
sx q[0];
rz(-1.4564464) q[0];
x q[1];
rz(-2.50364) q[2];
sx q[2];
rz(-1.049713) q[2];
sx q[2];
rz(-2.6661154) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0605088) q[1];
sx q[1];
rz(-1.9366055) q[1];
sx q[1];
rz(0.27549705) q[1];
rz(-pi) q[2];
rz(0.63739325) q[3];
sx q[3];
rz(-2.0717607) q[3];
sx q[3];
rz(-0.65046179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7369507) q[2];
sx q[2];
rz(-1.5664132) q[2];
sx q[2];
rz(-0.85181063) q[2];
rz(-0.24438721) q[3];
sx q[3];
rz(-0.71072018) q[3];
sx q[3];
rz(-1.3062668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39563018) q[0];
sx q[0];
rz(-1.4816477) q[0];
sx q[0];
rz(3.0356044) q[0];
rz(-1.4251739) q[1];
sx q[1];
rz(-1.1499848) q[1];
sx q[1];
rz(2.0521767) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7558173) q[0];
sx q[0];
rz(-1.5013278) q[0];
sx q[0];
rz(1.8088732) q[0];
rz(-pi) q[1];
rz(-1.127943) q[2];
sx q[2];
rz(-2.2086655) q[2];
sx q[2];
rz(-1.4952212) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2069083) q[1];
sx q[1];
rz(-1.4438153) q[1];
sx q[1];
rz(1.2506574) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7884729) q[3];
sx q[3];
rz(-2.5158415) q[3];
sx q[3];
rz(-1.9754639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0561698) q[2];
sx q[2];
rz(-1.9499754) q[2];
sx q[2];
rz(-0.58758152) q[2];
rz(-1.2299296) q[3];
sx q[3];
rz(-0.39837024) q[3];
sx q[3];
rz(0.57268322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9375482) q[0];
sx q[0];
rz(-1.9119104) q[0];
sx q[0];
rz(-2.6958534) q[0];
rz(-0.84272376) q[1];
sx q[1];
rz(-0.80519599) q[1];
sx q[1];
rz(-2.3420948) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84265316) q[0];
sx q[0];
rz(-1.9439204) q[0];
sx q[0];
rz(-2.0040657) q[0];
x q[1];
rz(-0.92858814) q[2];
sx q[2];
rz(-2.1669496) q[2];
sx q[2];
rz(-0.10732574) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.93175353) q[1];
sx q[1];
rz(-0.60720217) q[1];
sx q[1];
rz(2.3926413) q[1];
rz(-pi) q[2];
rz(-1.6715253) q[3];
sx q[3];
rz(-2.5593649) q[3];
sx q[3];
rz(1.4034206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2004956) q[2];
sx q[2];
rz(-2.7636187) q[2];
sx q[2];
rz(0.7824347) q[2];
rz(-2.2109219) q[3];
sx q[3];
rz(-1.2717671) q[3];
sx q[3];
rz(2.8455287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74872577) q[0];
sx q[0];
rz(-2.1935232) q[0];
sx q[0];
rz(2.0178846) q[0];
rz(0.8404845) q[1];
sx q[1];
rz(-1.9978465) q[1];
sx q[1];
rz(-1.7488272) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9542342) q[0];
sx q[0];
rz(-1.6723044) q[0];
sx q[0];
rz(3.0925173) q[0];
rz(-pi) q[1];
rz(-0.56864777) q[2];
sx q[2];
rz(-2.7629182) q[2];
sx q[2];
rz(2.9616122) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6992174) q[1];
sx q[1];
rz(-1.3997119) q[1];
sx q[1];
rz(-1.2198001) q[1];
rz(-pi) q[2];
rz(-0.61711611) q[3];
sx q[3];
rz(-1.502862) q[3];
sx q[3];
rz(-1.2319433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.30354083) q[2];
sx q[2];
rz(-0.76798648) q[2];
sx q[2];
rz(1.5194019) q[2];
rz(1.0649902) q[3];
sx q[3];
rz(-1.0737373) q[3];
sx q[3];
rz(0.85552335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4143455) q[0];
sx q[0];
rz(-0.94731826) q[0];
sx q[0];
rz(-1.974768) q[0];
rz(-2.836152) q[1];
sx q[1];
rz(-2.0716397) q[1];
sx q[1];
rz(-2.6397612) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4763757) q[0];
sx q[0];
rz(-1.3759651) q[0];
sx q[0];
rz(-1.4754773) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2756349) q[2];
sx q[2];
rz(-0.9532477) q[2];
sx q[2];
rz(-2.7028529) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8055685) q[1];
sx q[1];
rz(-2.7689287) q[1];
sx q[1];
rz(-1.3678846) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5975941) q[3];
sx q[3];
rz(-2.0680313) q[3];
sx q[3];
rz(2.0759283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.66494232) q[2];
sx q[2];
rz(-0.75935894) q[2];
sx q[2];
rz(3.0090028) q[2];
rz(2.8442123) q[3];
sx q[3];
rz(-1.5404276) q[3];
sx q[3];
rz(0.69445777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7731758) q[0];
sx q[0];
rz(-1.4765803) q[0];
sx q[0];
rz(1.6126527) q[0];
rz(1.348314) q[1];
sx q[1];
rz(-1.7650083) q[1];
sx q[1];
rz(-0.37337676) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.901875) q[0];
sx q[0];
rz(-1.7267701) q[0];
sx q[0];
rz(-1.9330935) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48252784) q[2];
sx q[2];
rz(-1.8410701) q[2];
sx q[2];
rz(2.7729098) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8638525) q[1];
sx q[1];
rz(-1.2009632) q[1];
sx q[1];
rz(-0.86994008) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67383234) q[3];
sx q[3];
rz(-0.34075173) q[3];
sx q[3];
rz(1.0125404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0384963) q[2];
sx q[2];
rz(-0.080048397) q[2];
sx q[2];
rz(-2.6888964) q[2];
rz(0.90940851) q[3];
sx q[3];
rz(-1.012864) q[3];
sx q[3];
rz(-0.67210853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4368923) q[0];
sx q[0];
rz(-1.1672651) q[0];
sx q[0];
rz(1.1080909) q[0];
rz(2.3469901) q[1];
sx q[1];
rz(-1.662685) q[1];
sx q[1];
rz(-0.86659238) q[1];
rz(1.5431326) q[2];
sx q[2];
rz(-2.4504708) q[2];
sx q[2];
rz(-2.8219555) q[2];
rz(-0.45171236) q[3];
sx q[3];
rz(-2.3087043) q[3];
sx q[3];
rz(0.63895978) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
