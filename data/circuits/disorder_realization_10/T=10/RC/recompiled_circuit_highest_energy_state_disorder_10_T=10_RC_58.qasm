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
rz(-0.53520441) q[1];
sx q[1];
rz(1.9898131) q[1];
sx q[1];
rz(10.726816) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0326234) q[0];
sx q[0];
rz(-1.8501256) q[0];
sx q[0];
rz(-0.33786122) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5349242) q[2];
sx q[2];
rz(-2.2891993) q[2];
sx q[2];
rz(-2.1181698) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3801119) q[1];
sx q[1];
rz(-1.8093361) q[1];
sx q[1];
rz(0.41888361) q[1];
x q[2];
rz(1.5005712) q[3];
sx q[3];
rz(-1.781979) q[3];
sx q[3];
rz(-2.2156032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.56613049) q[2];
sx q[2];
rz(-0.29511109) q[2];
sx q[2];
rz(0.58594123) q[2];
rz(1.640004) q[3];
sx q[3];
rz(-1.274704) q[3];
sx q[3];
rz(-1.9408102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5457299) q[0];
sx q[0];
rz(-1.0965309) q[0];
sx q[0];
rz(2.0134266) q[0];
rz(-0.7496756) q[1];
sx q[1];
rz(-1.3355037) q[1];
sx q[1];
rz(1.658176) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2735398) q[0];
sx q[0];
rz(-0.85208396) q[0];
sx q[0];
rz(-0.63086381) q[0];
x q[1];
rz(-0.60083484) q[2];
sx q[2];
rz(-1.877583) q[2];
sx q[2];
rz(0.50150774) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5239657) q[1];
sx q[1];
rz(-0.59365827) q[1];
sx q[1];
rz(2.6685324) q[1];
rz(2.1690278) q[3];
sx q[3];
rz(-0.56275193) q[3];
sx q[3];
rz(2.5431354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.38773203) q[2];
sx q[2];
rz(-2.3591177) q[2];
sx q[2];
rz(-2.1369047) q[2];
rz(-0.0082155148) q[3];
sx q[3];
rz(-1.3722082) q[3];
sx q[3];
rz(-2.7042702) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32366556) q[0];
sx q[0];
rz(-1.8765457) q[0];
sx q[0];
rz(-1.0311968) q[0];
rz(-2.9464856) q[1];
sx q[1];
rz(-0.61385265) q[1];
sx q[1];
rz(0.88417792) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3510597) q[0];
sx q[0];
rz(-1.5740593) q[0];
sx q[0];
rz(1.5705646) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.45978893) q[2];
sx q[2];
rz(-1.7722691) q[2];
sx q[2];
rz(3.0386864) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7373335) q[1];
sx q[1];
rz(-1.9801723) q[1];
sx q[1];
rz(-3.0380556) q[1];
rz(-pi) q[2];
rz(-2.1572491) q[3];
sx q[3];
rz(-1.3399117) q[3];
sx q[3];
rz(-1.6572052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.6241793) q[2];
sx q[2];
rz(-1.2531345) q[2];
sx q[2];
rz(-0.45428983) q[2];
rz(-1.0034466) q[3];
sx q[3];
rz(-0.61383057) q[3];
sx q[3];
rz(-1.412089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2237332) q[0];
sx q[0];
rz(-2.7847325) q[0];
sx q[0];
rz(-2.3339363) q[0];
rz(-1.3937021) q[1];
sx q[1];
rz(-1.423577) q[1];
sx q[1];
rz(-1.7234939) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0271064) q[0];
sx q[0];
rz(-1.4252121) q[0];
sx q[0];
rz(1.1417861) q[0];
rz(-pi) q[1];
rz(-0.70576422) q[2];
sx q[2];
rz(-1.4250722) q[2];
sx q[2];
rz(0.55225295) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.56546375) q[1];
sx q[1];
rz(-1.189078) q[1];
sx q[1];
rz(1.5839427) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2149736) q[3];
sx q[3];
rz(-0.62320504) q[3];
sx q[3];
rz(0.96545593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7915446) q[2];
sx q[2];
rz(-1.2309265) q[2];
sx q[2];
rz(-3.0080504) q[2];
rz(0.40889016) q[3];
sx q[3];
rz(-3.0372527) q[3];
sx q[3];
rz(1.0094118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23311663) q[0];
sx q[0];
rz(-0.82670832) q[0];
sx q[0];
rz(-0.46434656) q[0];
rz(-0.50114477) q[1];
sx q[1];
rz(-2.8906288) q[1];
sx q[1];
rz(-2.4026925) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.319425) q[0];
sx q[0];
rz(-0.28114265) q[0];
sx q[0];
rz(-2.7168363) q[0];
x q[1];
rz(2.5986669) q[2];
sx q[2];
rz(-1.6871337) q[2];
sx q[2];
rz(0.1067208) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1170144) q[1];
sx q[1];
rz(-1.1305594) q[1];
sx q[1];
rz(1.3280895) q[1];
x q[2];
rz(-0.67023863) q[3];
sx q[3];
rz(-2.3730108) q[3];
sx q[3];
rz(-2.0266899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9008987) q[2];
sx q[2];
rz(-1.4229166) q[2];
sx q[2];
rz(2.9071729) q[2];
rz(-2.8068986) q[3];
sx q[3];
rz(-2.5321999) q[3];
sx q[3];
rz(1.4332244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(2.1471106) q[0];
sx q[0];
rz(-0.88279498) q[0];
sx q[0];
rz(-2.2301646) q[0];
rz(-1.3145187) q[1];
sx q[1];
rz(-2.283137) q[1];
sx q[1];
rz(-1.3988769) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0270777) q[0];
sx q[0];
rz(-2.0721779) q[0];
sx q[0];
rz(-2.1073226) q[0];
rz(0.21199356) q[2];
sx q[2];
rz(-0.14236406) q[2];
sx q[2];
rz(1.1039343) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4928455) q[1];
sx q[1];
rz(-1.4727108) q[1];
sx q[1];
rz(-2.0010038) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1067922) q[3];
sx q[3];
rz(-1.3437931) q[3];
sx q[3];
rz(1.5846399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92480245) q[2];
sx q[2];
rz(-2.3239467) q[2];
sx q[2];
rz(1.7693046) q[2];
rz(-1.6670082) q[3];
sx q[3];
rz(-2.0788914) q[3];
sx q[3];
rz(-2.8180928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.40397662) q[0];
sx q[0];
rz(-1.0219028) q[0];
sx q[0];
rz(-1.8811986) q[0];
rz(0.34126392) q[1];
sx q[1];
rz(-0.9287467) q[1];
sx q[1];
rz(2.0492679) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72427801) q[0];
sx q[0];
rz(-1.5515741) q[0];
sx q[0];
rz(0.5165944) q[0];
x q[1];
rz(-0.34185648) q[2];
sx q[2];
rz(-1.7539138) q[2];
sx q[2];
rz(1.3951071) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.51755136) q[1];
sx q[1];
rz(-0.95145345) q[1];
sx q[1];
rz(0.29631181) q[1];
rz(-pi) q[2];
rz(0.41938765) q[3];
sx q[3];
rz(-1.2817146) q[3];
sx q[3];
rz(-1.737864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9799389) q[2];
sx q[2];
rz(-1.3875763) q[2];
sx q[2];
rz(0.8650583) q[2];
rz(3.0480399) q[3];
sx q[3];
rz(-1.4253987) q[3];
sx q[3];
rz(3.0430999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.426067) q[0];
sx q[0];
rz(-0.97398296) q[0];
sx q[0];
rz(-1.4564212) q[0];
rz(0.23297019) q[1];
sx q[1];
rz(-1.3825682) q[1];
sx q[1];
rz(1.9219386) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15773179) q[0];
sx q[0];
rz(-1.8605352) q[0];
sx q[0];
rz(-2.4012171) q[0];
rz(-2.4039335) q[2];
sx q[2];
rz(-0.56471497) q[2];
sx q[2];
rz(-1.4454522) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.10286843) q[1];
sx q[1];
rz(-1.7448493) q[1];
sx q[1];
rz(0.37794387) q[1];
rz(2.574) q[3];
sx q[3];
rz(-0.78062468) q[3];
sx q[3];
rz(2.7989557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.098027078) q[2];
sx q[2];
rz(-1.7003912) q[2];
sx q[2];
rz(3.0230057) q[2];
rz(-0.81973997) q[3];
sx q[3];
rz(-0.44410646) q[3];
sx q[3];
rz(2.8290101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0176625) q[0];
sx q[0];
rz(-2.1899962) q[0];
sx q[0];
rz(0.96624017) q[0];
rz(-0.36695925) q[1];
sx q[1];
rz(-0.96261135) q[1];
sx q[1];
rz(0.56698322) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1408932) q[0];
sx q[0];
rz(-0.53017925) q[0];
sx q[0];
rz(0.40325608) q[0];
rz(-2.5091293) q[2];
sx q[2];
rz(-2.2603432) q[2];
sx q[2];
rz(2.0105711) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2046656) q[1];
sx q[1];
rz(-1.3673529) q[1];
sx q[1];
rz(0.70091796) q[1];
x q[2];
rz(-0.35225711) q[3];
sx q[3];
rz(-0.92995074) q[3];
sx q[3];
rz(0.0059222277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6734267) q[2];
sx q[2];
rz(-1.6750853) q[2];
sx q[2];
rz(2.1231988) q[2];
rz(-0.23032019) q[3];
sx q[3];
rz(-2.6789594) q[3];
sx q[3];
rz(0.066430062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(2.6224943) q[0];
sx q[0];
rz(-1.007217) q[0];
sx q[0];
rz(1.8100716) q[0];
rz(1.9271756) q[1];
sx q[1];
rz(-1.6969095) q[1];
sx q[1];
rz(0.4164947) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1399222) q[0];
sx q[0];
rz(-0.71437144) q[0];
sx q[0];
rz(-1.2561965) q[0];
rz(-pi) q[1];
rz(-2.8266818) q[2];
sx q[2];
rz(-1.8127894) q[2];
sx q[2];
rz(-2.7417208) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.13392049) q[1];
sx q[1];
rz(-2.4654572) q[1];
sx q[1];
rz(-2.424404) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0216398) q[3];
sx q[3];
rz(-1.6623673) q[3];
sx q[3];
rz(1.4705603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.74849621) q[2];
sx q[2];
rz(-1.2747719) q[2];
sx q[2];
rz(-1.7999125) q[2];
rz(2.0210576) q[3];
sx q[3];
rz(-1.4015965) q[3];
sx q[3];
rz(-3.024658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31711598) q[0];
sx q[0];
rz(-1.8335637) q[0];
sx q[0];
rz(1.2807922) q[0];
rz(2.169213) q[1];
sx q[1];
rz(-2.0449816) q[1];
sx q[1];
rz(2.4349946) q[1];
rz(-0.62452684) q[2];
sx q[2];
rz(-1.7869768) q[2];
sx q[2];
rz(-0.063135978) q[2];
rz(3.1092316) q[3];
sx q[3];
rz(-1.4799043) q[3];
sx q[3];
rz(1.6226488) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
