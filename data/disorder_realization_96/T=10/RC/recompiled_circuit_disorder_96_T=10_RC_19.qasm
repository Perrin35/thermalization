OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0918026) q[0];
sx q[0];
rz(-3.0135305) q[0];
sx q[0];
rz(-0.81737104) q[0];
rz(-2.1583537) q[1];
sx q[1];
rz(-2.6020738) q[1];
sx q[1];
rz(-1.9411545) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.565633) q[0];
sx q[0];
rz(-0.48086777) q[0];
sx q[0];
rz(-2.5984882) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7898516) q[2];
sx q[2];
rz(-2.5918505) q[2];
sx q[2];
rz(-1.654939) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5719205) q[1];
sx q[1];
rz(-1.306428) q[1];
sx q[1];
rz(2.4155248) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8957542) q[3];
sx q[3];
rz(-0.80245362) q[3];
sx q[3];
rz(1.6368395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1203221) q[2];
sx q[2];
rz(-2.5734084) q[2];
sx q[2];
rz(1.583064) q[2];
rz(0.99672404) q[3];
sx q[3];
rz(-0.45209) q[3];
sx q[3];
rz(0.42580095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5581756) q[0];
sx q[0];
rz(-2.3893864) q[0];
sx q[0];
rz(-0.054071991) q[0];
rz(1.9460829) q[1];
sx q[1];
rz(-1.0369438) q[1];
sx q[1];
rz(2.6057459) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8230096) q[0];
sx q[0];
rz(-2.5479925) q[0];
sx q[0];
rz(-1.7943322) q[0];
x q[1];
rz(-1.8230121) q[2];
sx q[2];
rz(-0.87715845) q[2];
sx q[2];
rz(-1.1064305) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9745001) q[1];
sx q[1];
rz(-1.8247316) q[1];
sx q[1];
rz(-0.42029917) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.348098) q[3];
sx q[3];
rz(-2.9950812) q[3];
sx q[3];
rz(0.033586249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0216996) q[2];
sx q[2];
rz(-2.0844441) q[2];
sx q[2];
rz(0.95834857) q[2];
rz(-0.066453233) q[3];
sx q[3];
rz(-1.5840014) q[3];
sx q[3];
rz(-0.44979969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7217343) q[0];
sx q[0];
rz(-1.9323213) q[0];
sx q[0];
rz(-2.9911175) q[0];
rz(0.45723215) q[1];
sx q[1];
rz(-0.91232863) q[1];
sx q[1];
rz(-3.1157852) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2496101) q[0];
sx q[0];
rz(-1.5907856) q[0];
sx q[0];
rz(-1.32094) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.843156) q[2];
sx q[2];
rz(-2.8189427) q[2];
sx q[2];
rz(1.1622365) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.64297134) q[1];
sx q[1];
rz(-1.1585304) q[1];
sx q[1];
rz(-0.69795124) q[1];
rz(-1.2529536) q[3];
sx q[3];
rz(-2.7405973) q[3];
sx q[3];
rz(2.3272115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0187443) q[2];
sx q[2];
rz(-2.7763425) q[2];
sx q[2];
rz(2.5562111) q[2];
rz(2.9600926) q[3];
sx q[3];
rz(-1.3320965) q[3];
sx q[3];
rz(-1.5766778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.240775) q[0];
sx q[0];
rz(-2.5188991) q[0];
sx q[0];
rz(-0.17661072) q[0];
rz(-0.88090849) q[1];
sx q[1];
rz(-1.0662339) q[1];
sx q[1];
rz(-0.53612971) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49318424) q[0];
sx q[0];
rz(-2.3097561) q[0];
sx q[0];
rz(2.1302845) q[0];
x q[1];
rz(0.23303194) q[2];
sx q[2];
rz(-2.8985902) q[2];
sx q[2];
rz(-2.255893) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0587479) q[1];
sx q[1];
rz(-1.9472329) q[1];
sx q[1];
rz(0.67614268) q[1];
x q[2];
rz(-1.2767775) q[3];
sx q[3];
rz(-2.2259568) q[3];
sx q[3];
rz(-0.40757195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.46999103) q[2];
sx q[2];
rz(-1.7207547) q[2];
sx q[2];
rz(-2.0969351) q[2];
rz(-0.70703834) q[3];
sx q[3];
rz(-0.98393327) q[3];
sx q[3];
rz(-2.7846591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(1.2351284) q[0];
sx q[0];
rz(-1.8179853) q[0];
sx q[0];
rz(-2.0902324) q[0];
rz(1.4936739) q[1];
sx q[1];
rz(-2.5525679) q[1];
sx q[1];
rz(-3.0984745) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9163141) q[0];
sx q[0];
rz(-1.7597223) q[0];
sx q[0];
rz(-2.180045) q[0];
x q[1];
rz(2.0448858) q[2];
sx q[2];
rz(-2.5708963) q[2];
sx q[2];
rz(1.2793465) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9777898) q[1];
sx q[1];
rz(-1.8352574) q[1];
sx q[1];
rz(-2.8603641) q[1];
rz(-2.99302) q[3];
sx q[3];
rz(-0.79509495) q[3];
sx q[3];
rz(-0.038392301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.12895) q[2];
sx q[2];
rz(-0.49762112) q[2];
sx q[2];
rz(-2.7894003) q[2];
rz(0.59018618) q[3];
sx q[3];
rz(-0.47368172) q[3];
sx q[3];
rz(0.56110704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.6234289) q[0];
sx q[0];
rz(-1.9264899) q[0];
sx q[0];
rz(-1.1556926) q[0];
rz(-0.75025264) q[1];
sx q[1];
rz(-2.2022088) q[1];
sx q[1];
rz(-2.0828784) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9583225) q[0];
sx q[0];
rz(-1.1407307) q[0];
sx q[0];
rz(1.3413315) q[0];
x q[1];
rz(-2.5541359) q[2];
sx q[2];
rz(-2.537478) q[2];
sx q[2];
rz(0.47746745) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.42029542) q[1];
sx q[1];
rz(-1.8719721) q[1];
sx q[1];
rz(0.23423127) q[1];
rz(2.433957) q[3];
sx q[3];
rz(-1.71873) q[3];
sx q[3];
rz(-1.1843475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.47026149) q[2];
sx q[2];
rz(-1.7241717) q[2];
sx q[2];
rz(-1.2188101) q[2];
rz(1.1550711) q[3];
sx q[3];
rz(-0.23854908) q[3];
sx q[3];
rz(-1.7003805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(3.0999775) q[0];
sx q[0];
rz(-1.8247373) q[0];
sx q[0];
rz(-1.51145) q[0];
rz(1.7639683) q[1];
sx q[1];
rz(-2.8306077) q[1];
sx q[1];
rz(0.84164936) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3232988) q[0];
sx q[0];
rz(-2.5677486) q[0];
sx q[0];
rz(0.42563514) q[0];
rz(-pi) q[1];
x q[1];
rz(0.26288962) q[2];
sx q[2];
rz(-0.6436231) q[2];
sx q[2];
rz(2.7123244) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.48982606) q[1];
sx q[1];
rz(-0.090432743) q[1];
sx q[1];
rz(-2.138278) q[1];
rz(-pi) q[2];
rz(-0.21861403) q[3];
sx q[3];
rz(-2.300005) q[3];
sx q[3];
rz(-1.5515755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.001361751) q[2];
sx q[2];
rz(-1.3683687) q[2];
sx q[2];
rz(-0.049953071) q[2];
rz(2.4800381) q[3];
sx q[3];
rz(-2.619132) q[3];
sx q[3];
rz(-0.18930999) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2426303) q[0];
sx q[0];
rz(-2.1147418) q[0];
sx q[0];
rz(-1.4021953) q[0];
rz(3.0461123) q[1];
sx q[1];
rz(-1.9752558) q[1];
sx q[1];
rz(2.7239674) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71224371) q[0];
sx q[0];
rz(-2.0040383) q[0];
sx q[0];
rz(-0.54746763) q[0];
rz(-pi) q[1];
rz(-3.0268961) q[2];
sx q[2];
rz(-1.9343997) q[2];
sx q[2];
rz(-2.7483658) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.046125267) q[1];
sx q[1];
rz(-2.1324665) q[1];
sx q[1];
rz(0.559505) q[1];
x q[2];
rz(2.2291064) q[3];
sx q[3];
rz(-2.6819326) q[3];
sx q[3];
rz(2.130079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3461356) q[2];
sx q[2];
rz(-1.4435578) q[2];
sx q[2];
rz(1.0428838) q[2];
rz(-2.4677094) q[3];
sx q[3];
rz(-1.4833114) q[3];
sx q[3];
rz(-0.90014443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3089356) q[0];
sx q[0];
rz(-1.4082264) q[0];
sx q[0];
rz(2.5119264) q[0];
rz(-2.5667403) q[1];
sx q[1];
rz(-1.8300627) q[1];
sx q[1];
rz(-0.94690698) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7787951) q[0];
sx q[0];
rz(-0.4501833) q[0];
sx q[0];
rz(-0.74624004) q[0];
rz(2.7526555) q[2];
sx q[2];
rz(-3.0367594) q[2];
sx q[2];
rz(-2.7742085) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7864101) q[1];
sx q[1];
rz(-1.2372412) q[1];
sx q[1];
rz(1.6598318) q[1];
x q[2];
rz(1.9584719) q[3];
sx q[3];
rz(-2.6689853) q[3];
sx q[3];
rz(-1.2848867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.56069121) q[2];
sx q[2];
rz(-2.0118482) q[2];
sx q[2];
rz(1.2488731) q[2];
rz(-0.71436626) q[3];
sx q[3];
rz(-1.8656105) q[3];
sx q[3];
rz(-2.8760288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(2.0666075) q[0];
sx q[0];
rz(-2.5191436) q[0];
sx q[0];
rz(0.65504909) q[0];
rz(-0.89637268) q[1];
sx q[1];
rz(-0.90677774) q[1];
sx q[1];
rz(0.64430976) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9172168) q[0];
sx q[0];
rz(-1.3972358) q[0];
sx q[0];
rz(1.5149084) q[0];
x q[1];
rz(1.8337433) q[2];
sx q[2];
rz(-0.73336468) q[2];
sx q[2];
rz(-2.3761689) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.555968) q[1];
sx q[1];
rz(-1.854419) q[1];
sx q[1];
rz(-0.3507627) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0993016) q[3];
sx q[3];
rz(-1.1700556) q[3];
sx q[3];
rz(-0.61939643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7908988) q[2];
sx q[2];
rz(-1.6692946) q[2];
sx q[2];
rz(-0.59990668) q[2];
rz(0.89896262) q[3];
sx q[3];
rz(-0.18342429) q[3];
sx q[3];
rz(-1.3658587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8469289) q[0];
sx q[0];
rz(-1.8142721) q[0];
sx q[0];
rz(-0.66692764) q[0];
rz(0.22944336) q[1];
sx q[1];
rz(-0.89090092) q[1];
sx q[1];
rz(0.13577239) q[1];
rz(-3.0145666) q[2];
sx q[2];
rz(-0.95022485) q[2];
sx q[2];
rz(-2.7765204) q[2];
rz(-2.0471845) q[3];
sx q[3];
rz(-0.73275685) q[3];
sx q[3];
rz(-2.3263596) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];