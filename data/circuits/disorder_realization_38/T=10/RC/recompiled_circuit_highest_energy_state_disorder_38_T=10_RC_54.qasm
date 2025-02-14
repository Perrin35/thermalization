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
rz(2.8528557) q[0];
sx q[0];
rz(5.5878162) q[0];
sx q[0];
rz(6.5509808) q[0];
rz(0.42203045) q[1];
sx q[1];
rz(4.0623436) q[1];
sx q[1];
rz(10.698591) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52142329) q[0];
sx q[0];
rz(-1.6156989) q[0];
sx q[0];
rz(-1.3292759) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2620413) q[2];
sx q[2];
rz(-2.3143907) q[2];
sx q[2];
rz(-0.66397053) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5112524) q[1];
sx q[1];
rz(-2.4773438) q[1];
sx q[1];
rz(-3.0306007) q[1];
rz(1.4273604) q[3];
sx q[3];
rz(-1.535245) q[3];
sx q[3];
rz(1.2811023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6196809) q[2];
sx q[2];
rz(-2.5002067) q[2];
sx q[2];
rz(-0.042595159) q[2];
rz(2.8604782) q[3];
sx q[3];
rz(-1.5646076) q[3];
sx q[3];
rz(-0.59578305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6302781) q[0];
sx q[0];
rz(-1.9673286) q[0];
sx q[0];
rz(2.0654772) q[0];
rz(-2.1108744) q[1];
sx q[1];
rz(-1.1117671) q[1];
sx q[1];
rz(-1.3105185) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.977518) q[0];
sx q[0];
rz(-2.4830472) q[0];
sx q[0];
rz(-1.5367299) q[0];
rz(-pi) q[1];
rz(0.4035455) q[2];
sx q[2];
rz(-0.87658823) q[2];
sx q[2];
rz(-0.36118868) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4736997) q[1];
sx q[1];
rz(-1.4034499) q[1];
sx q[1];
rz(-1.6909864) q[1];
rz(1.978068) q[3];
sx q[3];
rz(-2.0680475) q[3];
sx q[3];
rz(2.6866792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.217546) q[2];
sx q[2];
rz(-0.7520389) q[2];
sx q[2];
rz(-0.95620608) q[2];
rz(-1.3077959) q[3];
sx q[3];
rz(-0.81495133) q[3];
sx q[3];
rz(-3.0202878) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2995375) q[0];
sx q[0];
rz(-0.50899035) q[0];
sx q[0];
rz(-0.036238413) q[0];
rz(-2.3402975) q[1];
sx q[1];
rz(-1.5472629) q[1];
sx q[1];
rz(-1.8410199) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56604859) q[0];
sx q[0];
rz(-1.6774208) q[0];
sx q[0];
rz(1.410065) q[0];
rz(0.97359263) q[2];
sx q[2];
rz(-1.6064715) q[2];
sx q[2];
rz(-0.0168456) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6319148) q[1];
sx q[1];
rz(-1.011112) q[1];
sx q[1];
rz(0.51559971) q[1];
rz(-0.43576305) q[3];
sx q[3];
rz(-0.73832694) q[3];
sx q[3];
rz(0.83454013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7654045) q[2];
sx q[2];
rz(-1.3724962) q[2];
sx q[2];
rz(-0.21793951) q[2];
rz(2.229522) q[3];
sx q[3];
rz(-1.4927161) q[3];
sx q[3];
rz(-0.85165858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7243778) q[0];
sx q[0];
rz(-2.0700924) q[0];
sx q[0];
rz(-0.81556129) q[0];
rz(2.0888445) q[1];
sx q[1];
rz(-0.88200724) q[1];
sx q[1];
rz(1.7900593) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6725051) q[0];
sx q[0];
rz(-2.6658305) q[0];
sx q[0];
rz(-0.16938727) q[0];
rz(-pi) q[1];
rz(3.0409524) q[2];
sx q[2];
rz(-0.71906861) q[2];
sx q[2];
rz(-1.4532064) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0114944) q[1];
sx q[1];
rz(-1.5421499) q[1];
sx q[1];
rz(-2.0462034) q[1];
x q[2];
rz(2.5791956) q[3];
sx q[3];
rz(-0.80005433) q[3];
sx q[3];
rz(2.7871015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.012933) q[2];
sx q[2];
rz(-1.1933051) q[2];
sx q[2];
rz(2.7093757) q[2];
rz(-0.84248078) q[3];
sx q[3];
rz(-2.1079) q[3];
sx q[3];
rz(-0.99561083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1603482) q[0];
sx q[0];
rz(-0.46537414) q[0];
sx q[0];
rz(0.55111849) q[0];
rz(-0.61800686) q[1];
sx q[1];
rz(-1.8564686) q[1];
sx q[1];
rz(1.0689703) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2129125) q[0];
sx q[0];
rz(-1.7846264) q[0];
sx q[0];
rz(-2.4545011) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3927116) q[2];
sx q[2];
rz(-0.62070337) q[2];
sx q[2];
rz(-1.701603) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5139937) q[1];
sx q[1];
rz(-0.74666903) q[1];
sx q[1];
rz(-1.7563266) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7154416) q[3];
sx q[3];
rz(-2.0791884) q[3];
sx q[3];
rz(1.4163464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7908287) q[2];
sx q[2];
rz(-2.5574234) q[2];
sx q[2];
rz(-2.186415) q[2];
rz(-0.59349924) q[3];
sx q[3];
rz(-2.1953526) q[3];
sx q[3];
rz(-1.3957297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(1.7668358) q[0];
sx q[0];
rz(-3.0992442) q[0];
sx q[0];
rz(1.7497077) q[0];
rz(-1.1426686) q[1];
sx q[1];
rz(-1.3418158) q[1];
sx q[1];
rz(-1.3124189) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56579933) q[0];
sx q[0];
rz(-0.84747073) q[0];
sx q[0];
rz(-2.2295843) q[0];
rz(-pi) q[1];
rz(0.80760132) q[2];
sx q[2];
rz(-0.86059216) q[2];
sx q[2];
rz(-1.0654895) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3692664) q[1];
sx q[1];
rz(-1.516894) q[1];
sx q[1];
rz(-0.51477706) q[1];
x q[2];
rz(0.057891616) q[3];
sx q[3];
rz(-1.8857919) q[3];
sx q[3];
rz(-2.0712899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7890847) q[2];
sx q[2];
rz(-0.75560537) q[2];
sx q[2];
rz(1.4136723) q[2];
rz(-0.46755725) q[3];
sx q[3];
rz(-2.4831725) q[3];
sx q[3];
rz(2.5889034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5107875) q[0];
sx q[0];
rz(-0.063022114) q[0];
sx q[0];
rz(-2.4600273) q[0];
rz(-1.1850146) q[1];
sx q[1];
rz(-0.76528913) q[1];
sx q[1];
rz(-2.3550745) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1543321) q[0];
sx q[0];
rz(-1.2367147) q[0];
sx q[0];
rz(2.7272237) q[0];
rz(-2.2165197) q[2];
sx q[2];
rz(-2.4074915) q[2];
sx q[2];
rz(2.2112598) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6326846) q[1];
sx q[1];
rz(-0.68636319) q[1];
sx q[1];
rz(-1.7401766) q[1];
rz(1.6166572) q[3];
sx q[3];
rz(-2.4099382) q[3];
sx q[3];
rz(1.5523486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6098392) q[2];
sx q[2];
rz(-2.9015151) q[2];
sx q[2];
rz(1.5698203) q[2];
rz(1.8543367) q[3];
sx q[3];
rz(-1.6183034) q[3];
sx q[3];
rz(-0.94304812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59931961) q[0];
sx q[0];
rz(-2.4241408) q[0];
sx q[0];
rz(1.188311) q[0];
rz(3.1207454) q[1];
sx q[1];
rz(-0.7531082) q[1];
sx q[1];
rz(0.59897024) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1513903) q[0];
sx q[0];
rz(-1.6693475) q[0];
sx q[0];
rz(2.0771785) q[0];
rz(-pi) q[1];
rz(-1.125419) q[2];
sx q[2];
rz(-0.76596224) q[2];
sx q[2];
rz(2.9474023) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7341344) q[1];
sx q[1];
rz(-2.2770513) q[1];
sx q[1];
rz(0.64532537) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7008301) q[3];
sx q[3];
rz(-2.2513933) q[3];
sx q[3];
rz(1.1632533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8591566) q[2];
sx q[2];
rz(-1.2500637) q[2];
sx q[2];
rz(-0.51521987) q[2];
rz(-0.53269261) q[3];
sx q[3];
rz(-1.0849413) q[3];
sx q[3];
rz(-0.93713078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1135947) q[0];
sx q[0];
rz(-2.4585215) q[0];
sx q[0];
rz(-0.37044507) q[0];
rz(-2.0619552) q[1];
sx q[1];
rz(-0.41079435) q[1];
sx q[1];
rz(-3.0830141) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.661299) q[0];
sx q[0];
rz(-0.96503557) q[0];
sx q[0];
rz(1.763271) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7981766) q[2];
sx q[2];
rz(-1.6140964) q[2];
sx q[2];
rz(-0.82466489) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6296719) q[1];
sx q[1];
rz(-1.6663934) q[1];
sx q[1];
rz(-2.3562576) q[1];
rz(-pi) q[2];
rz(-0.58038099) q[3];
sx q[3];
rz(-2.2182052) q[3];
sx q[3];
rz(0.035149487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8687245) q[2];
sx q[2];
rz(-0.803002) q[2];
sx q[2];
rz(-0.33099428) q[2];
rz(-3.1281779) q[3];
sx q[3];
rz(-2.1881073) q[3];
sx q[3];
rz(1.0564055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(0.57824221) q[0];
sx q[0];
rz(-2.712482) q[0];
sx q[0];
rz(-2.6934534) q[0];
rz(-3.0478802) q[1];
sx q[1];
rz(-0.30148503) q[1];
sx q[1];
rz(-1.9261446) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96998065) q[0];
sx q[0];
rz(-1.5413862) q[0];
sx q[0];
rz(2.0038811) q[0];
x q[1];
rz(-2.3289069) q[2];
sx q[2];
rz(-1.2764837) q[2];
sx q[2];
rz(0.049969604) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3642973) q[1];
sx q[1];
rz(-2.5357995) q[1];
sx q[1];
rz(-1.08849) q[1];
x q[2];
rz(-3.0057109) q[3];
sx q[3];
rz(-1.1757506) q[3];
sx q[3];
rz(-2.6108133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.965968) q[2];
sx q[2];
rz(-1.5886687) q[2];
sx q[2];
rz(-0.69941163) q[2];
rz(1.2753963) q[3];
sx q[3];
rz(-1.9857152) q[3];
sx q[3];
rz(-1.8705961) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4324343) q[0];
sx q[0];
rz(-0.86031886) q[0];
sx q[0];
rz(1.1501089) q[0];
rz(-1.3954096) q[1];
sx q[1];
rz(-1.2184873) q[1];
sx q[1];
rz(1.7582735) q[1];
rz(3.0396661) q[2];
sx q[2];
rz(-0.66833767) q[2];
sx q[2];
rz(2.3705033) q[2];
rz(1.6897353) q[3];
sx q[3];
rz(-2.0343067) q[3];
sx q[3];
rz(0.35342356) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
