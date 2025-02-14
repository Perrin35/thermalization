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
rz(-0.91322652) q[0];
sx q[0];
rz(-0.93728596) q[0];
sx q[0];
rz(-0.32567853) q[0];
rz(3.3477793) q[1];
sx q[1];
rz(3.4833796) q[1];
sx q[1];
rz(5.7835328) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2316811) q[0];
sx q[0];
rz(-1.2308111) q[0];
sx q[0];
rz(2.489734) q[0];
x q[1];
rz(1.5743786) q[2];
sx q[2];
rz(-2.56977) q[2];
sx q[2];
rz(-2.6973558) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6898866) q[1];
sx q[1];
rz(-1.8763649) q[1];
sx q[1];
rz(-2.0042104) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1127171) q[3];
sx q[3];
rz(-1.6290974) q[3];
sx q[3];
rz(1.3227303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8832522) q[2];
sx q[2];
rz(-2.3005433) q[2];
sx q[2];
rz(-0.19634761) q[2];
rz(-1.1697065) q[3];
sx q[3];
rz(-1.615808) q[3];
sx q[3];
rz(1.660478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36160463) q[0];
sx q[0];
rz(-1.9444281) q[0];
sx q[0];
rz(-0.25413221) q[0];
rz(-2.3669071) q[1];
sx q[1];
rz(-1.4137555) q[1];
sx q[1];
rz(-0.69001251) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38063506) q[0];
sx q[0];
rz(-1.9467783) q[0];
sx q[0];
rz(1.4248821) q[0];
rz(-2.5326917) q[2];
sx q[2];
rz(-1.5821067) q[2];
sx q[2];
rz(3.0593556) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.17711267) q[1];
sx q[1];
rz(-1.8928253) q[1];
sx q[1];
rz(-2.7142375) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8514107) q[3];
sx q[3];
rz(-2.0737344) q[3];
sx q[3];
rz(0.72301722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.1095573) q[2];
sx q[2];
rz(-3.1041225) q[2];
sx q[2];
rz(-1.0008585) q[2];
rz(1.9199269) q[3];
sx q[3];
rz(-1.6358717) q[3];
sx q[3];
rz(-0.033230335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.074653305) q[0];
sx q[0];
rz(-3.1237488) q[0];
sx q[0];
rz(-2.6113206) q[0];
rz(-0.52337921) q[1];
sx q[1];
rz(-1.4396079) q[1];
sx q[1];
rz(1.9022) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3227993) q[0];
sx q[0];
rz(-1.8086495) q[0];
sx q[0];
rz(0.6799233) q[0];
rz(-1.3472413) q[2];
sx q[2];
rz(-1.3251148) q[2];
sx q[2];
rz(-1.9839191) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0066632) q[1];
sx q[1];
rz(-1.8678209) q[1];
sx q[1];
rz(-3.013738) q[1];
rz(-pi) q[2];
rz(0.40145603) q[3];
sx q[3];
rz(-1.9250416) q[3];
sx q[3];
rz(-0.73520943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8221028) q[2];
sx q[2];
rz(-2.4097061) q[2];
sx q[2];
rz(-0.11719318) q[2];
rz(2.9486837) q[3];
sx q[3];
rz(-1.984805) q[3];
sx q[3];
rz(2.9624654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8837226) q[0];
sx q[0];
rz(-1.9921046) q[0];
sx q[0];
rz(-3.1174739) q[0];
rz(2.3564677) q[1];
sx q[1];
rz(-1.5784011) q[1];
sx q[1];
rz(-1.9804573) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6265921) q[0];
sx q[0];
rz(-0.75567317) q[0];
sx q[0];
rz(2.021778) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9865667) q[2];
sx q[2];
rz(-0.13992684) q[2];
sx q[2];
rz(-2.4292912) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9332446) q[1];
sx q[1];
rz(-2.1365494) q[1];
sx q[1];
rz(1.9557593) q[1];
x q[2];
rz(2.9579453) q[3];
sx q[3];
rz(-2.1602294) q[3];
sx q[3];
rz(-0.69212427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.44196096) q[2];
sx q[2];
rz(-2.4157951) q[2];
sx q[2];
rz(0.3886784) q[2];
rz(-1.9937438) q[3];
sx q[3];
rz(-1.1286705) q[3];
sx q[3];
rz(-1.283099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1175784) q[0];
sx q[0];
rz(-2.3021181) q[0];
sx q[0];
rz(-2.1003387) q[0];
rz(-0.72192347) q[1];
sx q[1];
rz(-1.8903774) q[1];
sx q[1];
rz(1.8183244) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6206317) q[0];
sx q[0];
rz(-2.8537321) q[0];
sx q[0];
rz(2.5858712) q[0];
rz(-pi) q[1];
rz(-1.2385246) q[2];
sx q[2];
rz(-2.4715444) q[2];
sx q[2];
rz(-0.61720309) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1309138) q[1];
sx q[1];
rz(-1.8478573) q[1];
sx q[1];
rz(2.4127736) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3106662) q[3];
sx q[3];
rz(-1.8046265) q[3];
sx q[3];
rz(1.1375858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.84956518) q[2];
sx q[2];
rz(-2.5490856) q[2];
sx q[2];
rz(-2.3805857) q[2];
rz(-1.546953) q[3];
sx q[3];
rz(-1.2706815) q[3];
sx q[3];
rz(2.7775619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4407235) q[0];
sx q[0];
rz(-0.64467621) q[0];
sx q[0];
rz(-1.122129) q[0];
rz(-2.2286277) q[1];
sx q[1];
rz(-2.2367621) q[1];
sx q[1];
rz(-2.1838271) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58272254) q[0];
sx q[0];
rz(-0.98419154) q[0];
sx q[0];
rz(-3.0987958) q[0];
rz(-0.087052778) q[2];
sx q[2];
rz(-1.7886482) q[2];
sx q[2];
rz(-1.717928) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4563728) q[1];
sx q[1];
rz(-1.3351016) q[1];
sx q[1];
rz(-0.46044223) q[1];
x q[2];
rz(0.60868355) q[3];
sx q[3];
rz(-2.9306378) q[3];
sx q[3];
rz(1.8041478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9271348) q[2];
sx q[2];
rz(-2.1467777) q[2];
sx q[2];
rz(-2.4289995) q[2];
rz(-1.7566977) q[3];
sx q[3];
rz(-1.8637928) q[3];
sx q[3];
rz(0.28759292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49823847) q[0];
sx q[0];
rz(-0.21182984) q[0];
sx q[0];
rz(-0.24464749) q[0];
rz(3.0554166) q[1];
sx q[1];
rz(-2.0241604) q[1];
sx q[1];
rz(0.92775956) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.010205) q[0];
sx q[0];
rz(-0.82021111) q[0];
sx q[0];
rz(-1.6774235) q[0];
x q[1];
rz(1.7150925) q[2];
sx q[2];
rz(-1.7789911) q[2];
sx q[2];
rz(-2.352072) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2725506) q[1];
sx q[1];
rz(-2.1107337) q[1];
sx q[1];
rz(-2.115095) q[1];
rz(-pi) q[2];
rz(-0.93337977) q[3];
sx q[3];
rz(-2.1669845) q[3];
sx q[3];
rz(0.0010871452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1722395) q[2];
sx q[2];
rz(-2.4805706) q[2];
sx q[2];
rz(-1.6924525) q[2];
rz(-2.0701854) q[3];
sx q[3];
rz(-1.7457733) q[3];
sx q[3];
rz(1.2034108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0485359) q[0];
sx q[0];
rz(-0.14126784) q[0];
sx q[0];
rz(2.5988044) q[0];
rz(0.12611783) q[1];
sx q[1];
rz(-1.7593971) q[1];
sx q[1];
rz(-2.9468367) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31349438) q[0];
sx q[0];
rz(-0.53088218) q[0];
sx q[0];
rz(-0.33057018) q[0];
rz(-pi) q[1];
rz(-0.82055494) q[2];
sx q[2];
rz(-0.61167292) q[2];
sx q[2];
rz(0.95407971) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4808274) q[1];
sx q[1];
rz(-1.387537) q[1];
sx q[1];
rz(2.6916958) q[1];
x q[2];
rz(0.76169059) q[3];
sx q[3];
rz(-1.7045827) q[3];
sx q[3];
rz(1.9445668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9572218) q[2];
sx q[2];
rz(-1.7879282) q[2];
sx q[2];
rz(1.4380737) q[2];
rz(-0.67972368) q[3];
sx q[3];
rz(-0.34346911) q[3];
sx q[3];
rz(2.9168985) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1723802) q[0];
sx q[0];
rz(-0.42177105) q[0];
sx q[0];
rz(1.0129741) q[0];
rz(0.29533932) q[1];
sx q[1];
rz(-1.4238009) q[1];
sx q[1];
rz(1.0257592) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.140705) q[0];
sx q[0];
rz(-1.6955339) q[0];
sx q[0];
rz(-0.38356218) q[0];
rz(-pi) q[1];
rz(-2.5466333) q[2];
sx q[2];
rz(-2.0987217) q[2];
sx q[2];
rz(2.4558112) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.48183) q[1];
sx q[1];
rz(-2.8952878) q[1];
sx q[1];
rz(-2.984761) q[1];
x q[2];
rz(-1.8904311) q[3];
sx q[3];
rz(-0.87351528) q[3];
sx q[3];
rz(0.13673328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.6369796) q[2];
sx q[2];
rz(-2.8490318) q[2];
sx q[2];
rz(-2.2620849) q[2];
rz(0.28359908) q[3];
sx q[3];
rz(-1.9234383) q[3];
sx q[3];
rz(1.2951736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.857169) q[0];
sx q[0];
rz(-1.5791945) q[0];
sx q[0];
rz(0.6012342) q[0];
rz(-0.55016905) q[1];
sx q[1];
rz(-2.1255122) q[1];
sx q[1];
rz(2.8395431) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6876533) q[0];
sx q[0];
rz(-2.2317356) q[0];
sx q[0];
rz(2.3783408) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43688628) q[2];
sx q[2];
rz(-2.1648975) q[2];
sx q[2];
rz(-0.64615102) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.35827929) q[1];
sx q[1];
rz(-2.9960592) q[1];
sx q[1];
rz(-3.1118433) q[1];
x q[2];
rz(1.0022903) q[3];
sx q[3];
rz(-1.7712815) q[3];
sx q[3];
rz(-2.025163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0740795) q[2];
sx q[2];
rz(-0.16177598) q[2];
sx q[2];
rz(1.8174685) q[2];
rz(0.4829123) q[3];
sx q[3];
rz(-2.3528779) q[3];
sx q[3];
rz(-0.27831349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.0637689) q[0];
sx q[0];
rz(-1.9763197) q[0];
sx q[0];
rz(1.9015953) q[0];
rz(0.77156144) q[1];
sx q[1];
rz(-1.1411219) q[1];
sx q[1];
rz(0.065446767) q[1];
rz(-0.83648079) q[2];
sx q[2];
rz(-1.4118839) q[2];
sx q[2];
rz(-1.6250162) q[2];
rz(2.2167233) q[3];
sx q[3];
rz(-0.51565167) q[3];
sx q[3];
rz(-2.6272163) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
