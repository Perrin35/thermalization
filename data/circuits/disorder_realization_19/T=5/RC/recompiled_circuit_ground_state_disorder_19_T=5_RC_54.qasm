OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.26025298) q[0];
sx q[0];
rz(-2.2007359) q[0];
sx q[0];
rz(0.22766222) q[0];
rz(0.28336278) q[1];
sx q[1];
rz(-2.7222187) q[1];
sx q[1];
rz(1.8546606) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57878198) q[0];
sx q[0];
rz(-1.1039886) q[0];
sx q[0];
rz(-0.060725529) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1212767) q[2];
sx q[2];
rz(-1.9971022) q[2];
sx q[2];
rz(-0.21838494) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.097295105) q[1];
sx q[1];
rz(-1.0781647) q[1];
sx q[1];
rz(-0.45483847) q[1];
rz(-pi) q[2];
rz(0.6391949) q[3];
sx q[3];
rz(-1.8248282) q[3];
sx q[3];
rz(2.3237236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.8925573) q[2];
sx q[2];
rz(-2.1302569) q[2];
sx q[2];
rz(2.2300143) q[2];
rz(-2.3859731) q[3];
sx q[3];
rz(-0.32114649) q[3];
sx q[3];
rz(2.6452981) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78254533) q[0];
sx q[0];
rz(-1.8308715) q[0];
sx q[0];
rz(-1.1388592) q[0];
rz(2.5149939) q[1];
sx q[1];
rz(-2.7732924) q[1];
sx q[1];
rz(-2.0638594) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12909938) q[0];
sx q[0];
rz(-1.9097304) q[0];
sx q[0];
rz(-0.61130543) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35524345) q[2];
sx q[2];
rz(-1.1525103) q[2];
sx q[2];
rz(-0.21165161) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0976757) q[1];
sx q[1];
rz(-1.3070135) q[1];
sx q[1];
rz(1.5013807) q[1];
rz(-pi) q[2];
rz(-1.0913565) q[3];
sx q[3];
rz(-0.78244996) q[3];
sx q[3];
rz(-1.2408011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6935912) q[2];
sx q[2];
rz(-0.75746626) q[2];
sx q[2];
rz(2.9288911) q[2];
rz(1.5564144) q[3];
sx q[3];
rz(-2.3978265) q[3];
sx q[3];
rz(2.0785544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0001275) q[0];
sx q[0];
rz(-2.8746222) q[0];
sx q[0];
rz(0.84685999) q[0];
rz(-2.8894539) q[1];
sx q[1];
rz(-2.3138901) q[1];
sx q[1];
rz(-2.491378) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6221049) q[0];
sx q[0];
rz(-1.6148003) q[0];
sx q[0];
rz(-2.8027338) q[0];
rz(-2.479662) q[2];
sx q[2];
rz(-2.227894) q[2];
sx q[2];
rz(-1.3160694) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.099459186) q[1];
sx q[1];
rz(-0.317527) q[1];
sx q[1];
rz(1.2088163) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7103151) q[3];
sx q[3];
rz(-0.85429885) q[3];
sx q[3];
rz(0.34654472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3113159) q[2];
sx q[2];
rz(-2.4429784) q[2];
sx q[2];
rz(2.1777731) q[2];
rz(-1.3612932) q[3];
sx q[3];
rz(-0.45069525) q[3];
sx q[3];
rz(0.042958766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90218246) q[0];
sx q[0];
rz(-0.22265156) q[0];
sx q[0];
rz(-2.1233001) q[0];
rz(2.2564383) q[1];
sx q[1];
rz(-2.8926909) q[1];
sx q[1];
rz(-0.28908602) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8967696) q[0];
sx q[0];
rz(-0.45216783) q[0];
sx q[0];
rz(0.74613476) q[0];
rz(-pi) q[1];
rz(1.9924329) q[2];
sx q[2];
rz(-2.1648063) q[2];
sx q[2];
rz(-2.8624812) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.30384597) q[1];
sx q[1];
rz(-2.0833384) q[1];
sx q[1];
rz(-2.9339177) q[1];
rz(-1.698367) q[3];
sx q[3];
rz(-1.3535548) q[3];
sx q[3];
rz(-2.9590817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3503512) q[2];
sx q[2];
rz(-1.505932) q[2];
sx q[2];
rz(1.852847) q[2];
rz(-0.57981235) q[3];
sx q[3];
rz(-0.47477397) q[3];
sx q[3];
rz(-0.60540664) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8346005) q[0];
sx q[0];
rz(-0.51486105) q[0];
sx q[0];
rz(-2.8173764) q[0];
rz(0.038837198) q[1];
sx q[1];
rz(-2.4034998) q[1];
sx q[1];
rz(2.0447581) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0457458) q[0];
sx q[0];
rz(-0.37558324) q[0];
sx q[0];
rz(0.71806851) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7261502) q[2];
sx q[2];
rz(-2.5553779) q[2];
sx q[2];
rz(1.4582576) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6935743) q[1];
sx q[1];
rz(-1.8724257) q[1];
sx q[1];
rz(-0.043971929) q[1];
rz(-pi) q[2];
rz(1.3665133) q[3];
sx q[3];
rz(-0.42821233) q[3];
sx q[3];
rz(2.7985559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0754806) q[2];
sx q[2];
rz(-2.8754063) q[2];
sx q[2];
rz(0.038662635) q[2];
rz(1.4085116) q[3];
sx q[3];
rz(-1.446529) q[3];
sx q[3];
rz(0.2814289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51327813) q[0];
sx q[0];
rz(-2.3230041) q[0];
sx q[0];
rz(2.1224838) q[0];
rz(2.6844773) q[1];
sx q[1];
rz(-0.26908427) q[1];
sx q[1];
rz(1.7105182) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5810878) q[0];
sx q[0];
rz(-1.771534) q[0];
sx q[0];
rz(-1.597658) q[0];
x q[1];
rz(1.2007295) q[2];
sx q[2];
rz(-0.6218172) q[2];
sx q[2];
rz(-1.3696826) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3574111) q[1];
sx q[1];
rz(-0.37115449) q[1];
sx q[1];
rz(-1.3264912) q[1];
rz(-2.2912628) q[3];
sx q[3];
rz(-1.8763233) q[3];
sx q[3];
rz(-3.1269493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9223601) q[2];
sx q[2];
rz(-1.3693501) q[2];
sx q[2];
rz(-2.5617981) q[2];
rz(-0.29948768) q[3];
sx q[3];
rz(-2.6087285) q[3];
sx q[3];
rz(1.9129491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2806468) q[0];
sx q[0];
rz(-1.4806643) q[0];
sx q[0];
rz(3.0378367) q[0];
rz(-2.5170028) q[1];
sx q[1];
rz(-2.2207405) q[1];
sx q[1];
rz(-1.3821028) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73653664) q[0];
sx q[0];
rz(-2.2782405) q[0];
sx q[0];
rz(0.28178431) q[0];
rz(-1.4598249) q[2];
sx q[2];
rz(-1.8689512) q[2];
sx q[2];
rz(-2.0046749) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.991532) q[1];
sx q[1];
rz(-1.0673017) q[1];
sx q[1];
rz(2.8343384) q[1];
rz(-pi) q[2];
rz(0.66094605) q[3];
sx q[3];
rz(-1.3437643) q[3];
sx q[3];
rz(1.1477136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.78911191) q[2];
sx q[2];
rz(-1.4752957) q[2];
sx q[2];
rz(-2.2769807) q[2];
rz(-0.23884808) q[3];
sx q[3];
rz(-0.75967234) q[3];
sx q[3];
rz(-2.4281832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1777986) q[0];
sx q[0];
rz(-2.5803784) q[0];
sx q[0];
rz(0.23502769) q[0];
rz(0.66559732) q[1];
sx q[1];
rz(-1.1382297) q[1];
sx q[1];
rz(1.6682909) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4427933) q[0];
sx q[0];
rz(-0.59663749) q[0];
sx q[0];
rz(-1.193861) q[0];
x q[1];
rz(0.26745112) q[2];
sx q[2];
rz(-2.6595734) q[2];
sx q[2];
rz(0.4617402) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0966592) q[1];
sx q[1];
rz(-2.9439264) q[1];
sx q[1];
rz(-2.7621619) q[1];
x q[2];
rz(-0.45110945) q[3];
sx q[3];
rz(-1.530297) q[3];
sx q[3];
rz(0.19865741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.16234806) q[2];
sx q[2];
rz(-0.20077106) q[2];
sx q[2];
rz(-2.6806504) q[2];
rz(1.0848684) q[3];
sx q[3];
rz(-0.74551398) q[3];
sx q[3];
rz(-0.78456867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.064086176) q[0];
sx q[0];
rz(-2.6106847) q[0];
sx q[0];
rz(-2.532646) q[0];
rz(2.5755836) q[1];
sx q[1];
rz(-1.9809664) q[1];
sx q[1];
rz(-1.2014679) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.214819) q[0];
sx q[0];
rz(-1.5626217) q[0];
sx q[0];
rz(-1.7063441) q[0];
x q[1];
rz(2.0588834) q[2];
sx q[2];
rz(-1.3859473) q[2];
sx q[2];
rz(2.3019997) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2347327) q[1];
sx q[1];
rz(-0.77488778) q[1];
sx q[1];
rz(3.0084684) q[1];
x q[2];
rz(-1.2668328) q[3];
sx q[3];
rz(-2.0018775) q[3];
sx q[3];
rz(1.4747185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0386049) q[2];
sx q[2];
rz(-1.0706341) q[2];
sx q[2];
rz(1.9496244) q[2];
rz(2.1206756) q[3];
sx q[3];
rz(-0.20191419) q[3];
sx q[3];
rz(1.6010223) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22333071) q[0];
sx q[0];
rz(-2.9409565) q[0];
sx q[0];
rz(-0.9675135) q[0];
rz(1.4406904) q[1];
sx q[1];
rz(-0.43165019) q[1];
sx q[1];
rz(0.38223019) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0456881) q[0];
sx q[0];
rz(-1.75215) q[0];
sx q[0];
rz(-0.99338326) q[0];
rz(-pi) q[1];
x q[1];
rz(1.154243) q[2];
sx q[2];
rz(-1.1692516) q[2];
sx q[2];
rz(-0.43972115) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9542071) q[1];
sx q[1];
rz(-1.6897196) q[1];
sx q[1];
rz(-0.78561659) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.554364) q[3];
sx q[3];
rz(-2.1167618) q[3];
sx q[3];
rz(0.16758238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.091577) q[2];
sx q[2];
rz(-0.86805582) q[2];
sx q[2];
rz(0.78563219) q[2];
rz(3.1146289) q[3];
sx q[3];
rz(-1.5188768) q[3];
sx q[3];
rz(0.54076076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3383269) q[0];
sx q[0];
rz(-1.4563518) q[0];
sx q[0];
rz(-0.8932054) q[0];
rz(-2.4772353) q[1];
sx q[1];
rz(-1.9021481) q[1];
sx q[1];
rz(-0.93217168) q[1];
rz(2.4128466) q[2];
sx q[2];
rz(-1.797429) q[2];
sx q[2];
rz(0.28354473) q[2];
rz(2.5369) q[3];
sx q[3];
rz(-1.5786377) q[3];
sx q[3];
rz(-0.64463401) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
