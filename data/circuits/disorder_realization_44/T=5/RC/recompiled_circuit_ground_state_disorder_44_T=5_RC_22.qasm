OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7527591) q[0];
sx q[0];
rz(1.692481) q[0];
sx q[0];
rz(11.115885) q[0];
rz(0.9736355) q[1];
sx q[1];
rz(-1.7042301) q[1];
sx q[1];
rz(2.2223284) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5333119) q[0];
sx q[0];
rz(-1.601378) q[0];
sx q[0];
rz(1.5536867) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4953142) q[2];
sx q[2];
rz(-2.9913104) q[2];
sx q[2];
rz(0.26963797) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.1303006) q[1];
sx q[1];
rz(-0.71349547) q[1];
sx q[1];
rz(1.43047) q[1];
rz(2.7204442) q[3];
sx q[3];
rz(-2.0055565) q[3];
sx q[3];
rz(-0.3374633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3322525) q[2];
sx q[2];
rz(-1.6925749) q[2];
sx q[2];
rz(-1.4130886) q[2];
rz(0.20279065) q[3];
sx q[3];
rz(-1.3821802) q[3];
sx q[3];
rz(3.0387759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24421144) q[0];
sx q[0];
rz(-1.0845217) q[0];
sx q[0];
rz(0.71075034) q[0];
rz(-2.3731025) q[1];
sx q[1];
rz(-1.0667421) q[1];
sx q[1];
rz(-2.1220727) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33776721) q[0];
sx q[0];
rz(-2.6927184) q[0];
sx q[0];
rz(-0.091218936) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7350082) q[2];
sx q[2];
rz(-0.43638849) q[2];
sx q[2];
rz(-1.9917038) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5078214) q[1];
sx q[1];
rz(-1.941676) q[1];
sx q[1];
rz(-0.033286496) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1250167) q[3];
sx q[3];
rz(-1.6416405) q[3];
sx q[3];
rz(1.4431825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.76356137) q[2];
sx q[2];
rz(-1.2445933) q[2];
sx q[2];
rz(1.2223318) q[2];
rz(1.2126806) q[3];
sx q[3];
rz(-2.1247532) q[3];
sx q[3];
rz(0.71162629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0843622) q[0];
sx q[0];
rz(-2.1570692) q[0];
sx q[0];
rz(1.2179751) q[0];
rz(2.6257264) q[1];
sx q[1];
rz(-2.5936544) q[1];
sx q[1];
rz(2.2191494) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3271752) q[0];
sx q[0];
rz(-1.6285768) q[0];
sx q[0];
rz(-1.6114651) q[0];
rz(-2.416196) q[2];
sx q[2];
rz(-1.9624406) q[2];
sx q[2];
rz(-1.2466516) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6092012) q[1];
sx q[1];
rz(-0.89840404) q[1];
sx q[1];
rz(-0.81873853) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32914583) q[3];
sx q[3];
rz(-1.4107879) q[3];
sx q[3];
rz(-3.126006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.018365232) q[2];
sx q[2];
rz(-2.0570698) q[2];
sx q[2];
rz(1.3398735) q[2];
rz(-0.54723048) q[3];
sx q[3];
rz(-2.7180143) q[3];
sx q[3];
rz(-0.71119285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.055534) q[0];
sx q[0];
rz(-2.9965897) q[0];
sx q[0];
rz(-0.15922971) q[0];
rz(3.1314462) q[1];
sx q[1];
rz(-1.0228913) q[1];
sx q[1];
rz(-2.8841282) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.793593) q[0];
sx q[0];
rz(-1.2158911) q[0];
sx q[0];
rz(-0.76323842) q[0];
rz(-pi) q[1];
rz(-0.79433673) q[2];
sx q[2];
rz(-1.5243425) q[2];
sx q[2];
rz(-0.72512324) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.60002335) q[1];
sx q[1];
rz(-2.5065055) q[1];
sx q[1];
rz(-0.16819185) q[1];
rz(-pi) q[2];
rz(-0.61473989) q[3];
sx q[3];
rz(-2.2343544) q[3];
sx q[3];
rz(3.1082982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3782392) q[2];
sx q[2];
rz(-1.2946318) q[2];
sx q[2];
rz(-0.47284687) q[2];
rz(2.4371448) q[3];
sx q[3];
rz(-1.7770146) q[3];
sx q[3];
rz(-0.64594597) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5064297) q[0];
sx q[0];
rz(-1.1744873) q[0];
sx q[0];
rz(3.0761062) q[0];
rz(-0.40924117) q[1];
sx q[1];
rz(-1.1301273) q[1];
sx q[1];
rz(1.4261036) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62140853) q[0];
sx q[0];
rz(-1.6637633) q[0];
sx q[0];
rz(-0.21958242) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.481856) q[2];
sx q[2];
rz(-1.9581902) q[2];
sx q[2];
rz(-0.26814869) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.1864422) q[1];
sx q[1];
rz(-1.7935866) q[1];
sx q[1];
rz(-2.8563315) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9687443) q[3];
sx q[3];
rz(-1.8573055) q[3];
sx q[3];
rz(2.2423173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.19491974) q[2];
sx q[2];
rz(-2.0544923) q[2];
sx q[2];
rz(1.3067513) q[2];
rz(1.170018) q[3];
sx q[3];
rz(-0.40478671) q[3];
sx q[3];
rz(-3.034333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4829247) q[0];
sx q[0];
rz(-1.985745) q[0];
sx q[0];
rz(2.7401127) q[0];
rz(-2.4688156) q[1];
sx q[1];
rz(-2.3428226) q[1];
sx q[1];
rz(2.2875517) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6333165) q[0];
sx q[0];
rz(-2.9190953) q[0];
sx q[0];
rz(1.7822687) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4638205) q[2];
sx q[2];
rz(-0.23636625) q[2];
sx q[2];
rz(1.4217699) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8484162) q[1];
sx q[1];
rz(-0.59976116) q[1];
sx q[1];
rz(-0.17677115) q[1];
x q[2];
rz(-2.0186508) q[3];
sx q[3];
rz(-1.0077701) q[3];
sx q[3];
rz(-0.26859586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4868769) q[2];
sx q[2];
rz(-2.2701023) q[2];
sx q[2];
rz(1.1775449) q[2];
rz(-2.5901637) q[3];
sx q[3];
rz(-1.5827554) q[3];
sx q[3];
rz(-2.3059755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7198782) q[0];
sx q[0];
rz(-0.28110176) q[0];
sx q[0];
rz(-1.1623435) q[0];
rz(-1.1516736) q[1];
sx q[1];
rz(-1.3538227) q[1];
sx q[1];
rz(2.2231359) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7699444) q[0];
sx q[0];
rz(-1.838883) q[0];
sx q[0];
rz(-1.2275342) q[0];
rz(-pi) q[1];
rz(0.8753885) q[2];
sx q[2];
rz(-1.7729442) q[2];
sx q[2];
rz(-3.1234158) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6035108) q[1];
sx q[1];
rz(-0.0099364837) q[1];
sx q[1];
rz(2.4007829) q[1];
x q[2];
rz(-2.9711668) q[3];
sx q[3];
rz(-1.2316576) q[3];
sx q[3];
rz(0.91811524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1071757) q[2];
sx q[2];
rz(-1.7566661) q[2];
sx q[2];
rz(-2.6386063) q[2];
rz(2.3668187) q[3];
sx q[3];
rz(-2.6796902) q[3];
sx q[3];
rz(-1.7983961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69304943) q[0];
sx q[0];
rz(-0.22668426) q[0];
sx q[0];
rz(1.5013129) q[0];
rz(0.34128183) q[1];
sx q[1];
rz(-1.7106067) q[1];
sx q[1];
rz(2.0223845) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78748465) q[0];
sx q[0];
rz(-1.9409927) q[0];
sx q[0];
rz(-2.5404929) q[0];
rz(2.881024) q[2];
sx q[2];
rz(-2.0379279) q[2];
sx q[2];
rz(1.9595944) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.95204845) q[1];
sx q[1];
rz(-0.91390007) q[1];
sx q[1];
rz(-1.7236962) q[1];
rz(-pi) q[2];
rz(1.4358892) q[3];
sx q[3];
rz(-1.9622318) q[3];
sx q[3];
rz(-1.6329444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.98562733) q[2];
sx q[2];
rz(-2.1817744) q[2];
sx q[2];
rz(2.7723374) q[2];
rz(2.8333832) q[3];
sx q[3];
rz(-0.9674558) q[3];
sx q[3];
rz(1.6109899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.2823328) q[0];
sx q[0];
rz(-1.8349324) q[0];
sx q[0];
rz(0.34307137) q[0];
rz(1.169091) q[1];
sx q[1];
rz(-1.3645376) q[1];
sx q[1];
rz(-1.4287359) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0051826) q[0];
sx q[0];
rz(-0.47629582) q[0];
sx q[0];
rz(2.5695557) q[0];
rz(-2.4506017) q[2];
sx q[2];
rz(-2.1804902) q[2];
sx q[2];
rz(0.26870773) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.27738849) q[1];
sx q[1];
rz(-0.38071796) q[1];
sx q[1];
rz(-1.7362795) q[1];
rz(-2.1200646) q[3];
sx q[3];
rz(-0.77199751) q[3];
sx q[3];
rz(0.87024161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.94065654) q[2];
sx q[2];
rz(-2.9328465) q[2];
sx q[2];
rz(-1.3915871) q[2];
rz(-0.40677795) q[3];
sx q[3];
rz(-1.4429561) q[3];
sx q[3];
rz(-2.2627635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10130356) q[0];
sx q[0];
rz(-2.5387796) q[0];
sx q[0];
rz(1.783675) q[0];
rz(-1.2376002) q[1];
sx q[1];
rz(-2.1289181) q[1];
sx q[1];
rz(2.3416669) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7800723) q[0];
sx q[0];
rz(-1.8006386) q[0];
sx q[0];
rz(-2.337268) q[0];
rz(1.8736035) q[2];
sx q[2];
rz(-1.0410415) q[2];
sx q[2];
rz(-0.15632665) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4373191) q[1];
sx q[1];
rz(-1.4156439) q[1];
sx q[1];
rz(-2.4680733) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.25000817) q[3];
sx q[3];
rz(-1.0984612) q[3];
sx q[3];
rz(1.3046622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.37821975) q[2];
sx q[2];
rz(-1.4241445) q[2];
sx q[2];
rz(0.073089449) q[2];
rz(-0.25589219) q[3];
sx q[3];
rz(-0.81232324) q[3];
sx q[3];
rz(-1.488744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(1.273461) q[0];
sx q[0];
rz(-1.68597) q[0];
sx q[0];
rz(-1.2706533) q[0];
rz(-1.3399667) q[1];
sx q[1];
rz(-1.7908962) q[1];
sx q[1];
rz(-2.4790196) q[1];
rz(0.70941464) q[2];
sx q[2];
rz(-1.5723036) q[2];
sx q[2];
rz(0.2607762) q[2];
rz(1.003391) q[3];
sx q[3];
rz(-2.1324674) q[3];
sx q[3];
rz(2.4218925) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
