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
rz(0.34974521) q[0];
sx q[0];
rz(5.3164696) q[0];
sx q[0];
rz(9.7860019) q[0];
rz(2.6063882) q[1];
sx q[1];
rz(-1.9898131) q[1];
sx q[1];
rz(-1.8395543) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0326234) q[0];
sx q[0];
rz(-1.8501256) q[0];
sx q[0];
rz(-0.33786122) q[0];
rz(0.7187219) q[2];
sx q[2];
rz(-1.5437922) q[2];
sx q[2];
rz(-2.570603) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.9141317) q[1];
sx q[1];
rz(-1.1644851) q[1];
sx q[1];
rz(-1.3106457) q[1];
x q[2];
rz(-0.21168904) q[3];
sx q[3];
rz(-1.6394588) q[3];
sx q[3];
rz(0.65955034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.56613049) q[2];
sx q[2];
rz(-0.29511109) q[2];
sx q[2];
rz(-2.5556514) q[2];
rz(-1.640004) q[3];
sx q[3];
rz(-1.8668886) q[3];
sx q[3];
rz(1.2007825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-0.59586278) q[0];
sx q[0];
rz(-2.0450617) q[0];
sx q[0];
rz(1.1281661) q[0];
rz(2.3919171) q[1];
sx q[1];
rz(-1.3355037) q[1];
sx q[1];
rz(-1.4834167) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2735398) q[0];
sx q[0];
rz(-0.85208396) q[0];
sx q[0];
rz(-2.5107288) q[0];
rz(2.5407578) q[2];
sx q[2];
rz(-1.877583) q[2];
sx q[2];
rz(0.50150774) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.0644599) q[1];
sx q[1];
rz(-2.0920447) q[1];
sx q[1];
rz(-1.2724933) q[1];
rz(-pi) q[2];
rz(-0.34134369) q[3];
sx q[3];
rz(-2.0273547) q[3];
sx q[3];
rz(-1.8648636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.38773203) q[2];
sx q[2];
rz(-2.3591177) q[2];
sx q[2];
rz(-1.0046879) q[2];
rz(-3.1333771) q[3];
sx q[3];
rz(-1.3722082) q[3];
sx q[3];
rz(-0.43732244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32366556) q[0];
sx q[0];
rz(-1.8765457) q[0];
sx q[0];
rz(1.0311968) q[0];
rz(2.9464856) q[1];
sx q[1];
rz(-2.52774) q[1];
sx q[1];
rz(0.88417792) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7802641) q[0];
sx q[0];
rz(-1.5705646) q[0];
sx q[0];
rz(-3.1383297) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7102406) q[2];
sx q[2];
rz(-0.49908733) q[2];
sx q[2];
rz(1.2896529) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.40425912) q[1];
sx q[1];
rz(-1.9801723) q[1];
sx q[1];
rz(3.0380556) q[1];
rz(-1.9724794) q[3];
sx q[3];
rz(-2.5163262) q[3];
sx q[3];
rz(0.4180983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.6241793) q[2];
sx q[2];
rz(-1.8884582) q[2];
sx q[2];
rz(0.45428983) q[2];
rz(1.0034466) q[3];
sx q[3];
rz(-2.5277621) q[3];
sx q[3];
rz(-1.412089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91785947) q[0];
sx q[0];
rz(-2.7847325) q[0];
sx q[0];
rz(0.80765635) q[0];
rz(-1.3937021) q[1];
sx q[1];
rz(-1.423577) q[1];
sx q[1];
rz(-1.7234939) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3781698) q[0];
sx q[0];
rz(-2.6900107) q[0];
sx q[0];
rz(-1.2319135) q[0];
rz(-pi) q[1];
rz(2.9190665) q[2];
sx q[2];
rz(-2.4234802) q[2];
sx q[2];
rz(1.1873174) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6114072) q[1];
sx q[1];
rz(-0.3819335) q[1];
sx q[1];
rz(0.032738222) q[1];
rz(-1.0491761) q[3];
sx q[3];
rz(-1.9289013) q[3];
sx q[3];
rz(-1.988609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.35004804) q[2];
sx q[2];
rz(-1.2309265) q[2];
sx q[2];
rz(-0.13354224) q[2];
rz(2.7327025) q[3];
sx q[3];
rz(-0.10433993) q[3];
sx q[3];
rz(1.0094118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-2.908476) q[0];
sx q[0];
rz(-2.3148843) q[0];
sx q[0];
rz(-2.6772461) q[0];
rz(-0.50114477) q[1];
sx q[1];
rz(-2.8906288) q[1];
sx q[1];
rz(0.73890013) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2621612) q[0];
sx q[0];
rz(-1.3152243) q[0];
sx q[0];
rz(-1.4523427) q[0];
rz(-pi) q[1];
x q[1];
rz(0.22245714) q[2];
sx q[2];
rz(-0.55402606) q[2];
sx q[2];
rz(-1.2740327) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6433324) q[1];
sx q[1];
rz(-0.49884379) q[1];
sx q[1];
rz(0.47175924) q[1];
rz(-2.49315) q[3];
sx q[3];
rz(-2.01727) q[3];
sx q[3];
rz(0.97383271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9008987) q[2];
sx q[2];
rz(-1.7186761) q[2];
sx q[2];
rz(0.23441976) q[2];
rz(-0.33469409) q[3];
sx q[3];
rz(-0.60939279) q[3];
sx q[3];
rz(-1.7083683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9944821) q[0];
sx q[0];
rz(-2.2587977) q[0];
sx q[0];
rz(-2.2301646) q[0];
rz(-1.3145187) q[1];
sx q[1];
rz(-2.283137) q[1];
sx q[1];
rz(1.7427157) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11451498) q[0];
sx q[0];
rz(-2.0721779) q[0];
sx q[0];
rz(1.03427) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6009459) q[2];
sx q[2];
rz(-1.7099524) q[2];
sx q[2];
rz(0.88983941) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4928455) q[1];
sx q[1];
rz(-1.4727108) q[1];
sx q[1];
rz(-1.1405888) q[1];
rz(1.9955548) q[3];
sx q[3];
rz(-2.5638608) q[3];
sx q[3];
rz(2.7932624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.92480245) q[2];
sx q[2];
rz(-0.81764597) q[2];
sx q[2];
rz(-1.372288) q[2];
rz(-1.6670082) q[3];
sx q[3];
rz(-2.0788914) q[3];
sx q[3];
rz(-2.8180928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40397662) q[0];
sx q[0];
rz(-1.0219028) q[0];
sx q[0];
rz(-1.2603941) q[0];
rz(-0.34126392) q[1];
sx q[1];
rz(-2.212846) q[1];
sx q[1];
rz(-1.0923247) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72427801) q[0];
sx q[0];
rz(-1.5900185) q[0];
sx q[0];
rz(-2.6249983) q[0];
rz(2.7997362) q[2];
sx q[2];
rz(-1.3876788) q[2];
sx q[2];
rz(1.7464856) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1398692) q[1];
sx q[1];
rz(-0.67808637) q[1];
sx q[1];
rz(1.1820611) q[1];
rz(2.722205) q[3];
sx q[3];
rz(-1.8598781) q[3];
sx q[3];
rz(-1.737864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1616538) q[2];
sx q[2];
rz(-1.3875763) q[2];
sx q[2];
rz(2.2765344) q[2];
rz(3.0480399) q[3];
sx q[3];
rz(-1.716194) q[3];
sx q[3];
rz(0.098492749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.426067) q[0];
sx q[0];
rz(-0.97398296) q[0];
sx q[0];
rz(1.4564212) q[0];
rz(-0.23297019) q[1];
sx q[1];
rz(-1.3825682) q[1];
sx q[1];
rz(1.2196541) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4254974) q[0];
sx q[0];
rz(-2.3566735) q[0];
sx q[0];
rz(-2.7254478) q[0];
rz(-2.4039335) q[2];
sx q[2];
rz(-2.5768777) q[2];
sx q[2];
rz(-1.6961404) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.0387242) q[1];
sx q[1];
rz(-1.7448493) q[1];
sx q[1];
rz(-2.7636488) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4457621) q[3];
sx q[3];
rz(-1.9587818) q[3];
sx q[3];
rz(1.6534905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0435656) q[2];
sx q[2];
rz(-1.7003912) q[2];
sx q[2];
rz(-0.11858693) q[2];
rz(2.3218527) q[3];
sx q[3];
rz(-0.44410646) q[3];
sx q[3];
rz(-0.31258252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1239301) q[0];
sx q[0];
rz(-2.1899962) q[0];
sx q[0];
rz(-0.96624017) q[0];
rz(0.36695925) q[1];
sx q[1];
rz(-2.1789813) q[1];
sx q[1];
rz(0.56698322) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1408932) q[0];
sx q[0];
rz(-0.53017925) q[0];
sx q[0];
rz(2.7383366) q[0];
rz(-pi) q[1];
rz(-2.5091293) q[2];
sx q[2];
rz(-2.2603432) q[2];
sx q[2];
rz(-1.1310215) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2046656) q[1];
sx q[1];
rz(-1.7742397) q[1];
sx q[1];
rz(2.4406747) q[1];
rz(-pi) q[2];
rz(0.35225711) q[3];
sx q[3];
rz(-2.2116419) q[3];
sx q[3];
rz(0.0059222277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.46816599) q[2];
sx q[2];
rz(-1.4665073) q[2];
sx q[2];
rz(2.1231988) q[2];
rz(-0.23032019) q[3];
sx q[3];
rz(-2.6789594) q[3];
sx q[3];
rz(-3.0751626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6224943) q[0];
sx q[0];
rz(-2.1343756) q[0];
sx q[0];
rz(1.8100716) q[0];
rz(1.9271756) q[1];
sx q[1];
rz(-1.6969095) q[1];
sx q[1];
rz(0.4164947) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7365624) q[0];
sx q[0];
rz(-0.89823898) q[0];
sx q[0];
rz(-0.26215464) q[0];
x q[1];
rz(1.3168066) q[2];
sx q[2];
rz(-1.2653627) q[2];
sx q[2];
rz(1.0930201) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0076722) q[1];
sx q[1];
rz(-0.67613542) q[1];
sx q[1];
rz(2.424404) q[1];
rz(-1.1199529) q[3];
sx q[3];
rz(-1.4792253) q[3];
sx q[3];
rz(-1.4705603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.74849621) q[2];
sx q[2];
rz(-1.2747719) q[2];
sx q[2];
rz(-1.3416802) q[2];
rz(1.120535) q[3];
sx q[3];
rz(-1.7399961) q[3];
sx q[3];
rz(0.11693461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31711598) q[0];
sx q[0];
rz(-1.308029) q[0];
sx q[0];
rz(-1.8608004) q[0];
rz(0.9723797) q[1];
sx q[1];
rz(-1.096611) q[1];
sx q[1];
rz(-0.70659804) q[1];
rz(0.62452684) q[2];
sx q[2];
rz(-1.3546158) q[2];
sx q[2];
rz(3.0784567) q[2];
rz(0.032361055) q[3];
sx q[3];
rz(-1.6616884) q[3];
sx q[3];
rz(-1.5189439) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
