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
rz(-0.62701464) q[0];
sx q[0];
rz(-2.505317) q[0];
sx q[0];
rz(1.0778435) q[0];
rz(1.0765422) q[1];
sx q[1];
rz(-0.62907469) q[1];
sx q[1];
rz(-1.0318626) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92562308) q[0];
sx q[0];
rz(-1.5711938) q[0];
sx q[0];
rz(0.0041603869) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7912895) q[2];
sx q[2];
rz(-1.7997348) q[2];
sx q[2];
rz(0.13267429) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.87776041) q[1];
sx q[1];
rz(-1.227637) q[1];
sx q[1];
rz(0.43107098) q[1];
rz(-0.21764619) q[3];
sx q[3];
rz(-2.7983694) q[3];
sx q[3];
rz(1.1962593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91118139) q[2];
sx q[2];
rz(-1.1831256) q[2];
sx q[2];
rz(2.6739056) q[2];
rz(-0.45105252) q[3];
sx q[3];
rz(-2.7428198) q[3];
sx q[3];
rz(2.8635645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69538799) q[0];
sx q[0];
rz(-2.9187293) q[0];
sx q[0];
rz(-2.9872802) q[0];
rz(0.34126869) q[1];
sx q[1];
rz(-0.35366615) q[1];
sx q[1];
rz(-0.42010677) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1695009) q[0];
sx q[0];
rz(-1.3779213) q[0];
sx q[0];
rz(2.4986096) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5163563) q[2];
sx q[2];
rz(-0.71085658) q[2];
sx q[2];
rz(2.9756851) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5678011) q[1];
sx q[1];
rz(-2.670799) q[1];
sx q[1];
rz(0.0030229645) q[1];
rz(-2.9734334) q[3];
sx q[3];
rz(-0.75440591) q[3];
sx q[3];
rz(1.8969632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5306065) q[2];
sx q[2];
rz(-0.89908081) q[2];
sx q[2];
rz(2.3685624) q[2];
rz(2.2981339) q[3];
sx q[3];
rz(-1.5777028) q[3];
sx q[3];
rz(0.57190603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17324363) q[0];
sx q[0];
rz(-2.1143715) q[0];
sx q[0];
rz(-2.9395043) q[0];
rz(2.659722) q[1];
sx q[1];
rz(-0.20129573) q[1];
sx q[1];
rz(0.906382) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8608241) q[0];
sx q[0];
rz(-1.1204136) q[0];
sx q[0];
rz(0.76669873) q[0];
x q[1];
rz(0.93744509) q[2];
sx q[2];
rz(-2.5966532) q[2];
sx q[2];
rz(2.3505806) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8123834) q[1];
sx q[1];
rz(-1.0339197) q[1];
sx q[1];
rz(-1.7574134) q[1];
rz(1.2754457) q[3];
sx q[3];
rz(-1.4100299) q[3];
sx q[3];
rz(2.8637342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.071094461) q[2];
sx q[2];
rz(-2.0281894) q[2];
sx q[2];
rz(-0.58007288) q[2];
rz(1.5602559) q[3];
sx q[3];
rz(-1.3908849) q[3];
sx q[3];
rz(-1.4238547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4105014) q[0];
sx q[0];
rz(-0.061256496) q[0];
sx q[0];
rz(3.0938003) q[0];
rz(-1.8675249) q[1];
sx q[1];
rz(-1.2893226) q[1];
sx q[1];
rz(-0.045225708) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83376081) q[0];
sx q[0];
rz(-0.93474301) q[0];
sx q[0];
rz(0.18524203) q[0];
rz(2.8130262) q[2];
sx q[2];
rz(-0.96146482) q[2];
sx q[2];
rz(-0.48664618) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8359591) q[1];
sx q[1];
rz(-0.0063414185) q[1];
sx q[1];
rz(-2.0672227) q[1];
x q[2];
rz(-2.2269656) q[3];
sx q[3];
rz(-2.4883929) q[3];
sx q[3];
rz(-2.7340305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4101326) q[2];
sx q[2];
rz(-1.2612017) q[2];
sx q[2];
rz(0.88978466) q[2];
rz(2.5951923) q[3];
sx q[3];
rz(-2.140464) q[3];
sx q[3];
rz(2.8304097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6379717) q[0];
sx q[0];
rz(-1.4542955) q[0];
sx q[0];
rz(-1.4171492) q[0];
rz(-2.0218938) q[1];
sx q[1];
rz(-1.7599301) q[1];
sx q[1];
rz(-1.959257) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.946163) q[0];
sx q[0];
rz(-0.39349213) q[0];
sx q[0];
rz(2.959743) q[0];
rz(-1.117302) q[2];
sx q[2];
rz(-1.1536479) q[2];
sx q[2];
rz(0.81062775) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.97302283) q[1];
sx q[1];
rz(-0.72271361) q[1];
sx q[1];
rz(-3.0672468) q[1];
rz(2.6695146) q[3];
sx q[3];
rz(-0.80754495) q[3];
sx q[3];
rz(-2.9387752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.917439) q[2];
sx q[2];
rz(-0.76442337) q[2];
sx q[2];
rz(0.42382851) q[2];
rz(2.0297) q[3];
sx q[3];
rz(-0.49911505) q[3];
sx q[3];
rz(2.4901938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79913419) q[0];
sx q[0];
rz(-3.0693711) q[0];
sx q[0];
rz(0.95241958) q[0];
rz(-0.77596387) q[1];
sx q[1];
rz(-2.2064078) q[1];
sx q[1];
rz(1.5663358) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.090388894) q[0];
sx q[0];
rz(-3.0984833) q[0];
sx q[0];
rz(-2.6531522) q[0];
rz(-pi) q[1];
rz(-3.1029019) q[2];
sx q[2];
rz(-2.1771447) q[2];
sx q[2];
rz(-2.4232466) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9989062) q[1];
sx q[1];
rz(-0.64780462) q[1];
sx q[1];
rz(-0.36869375) q[1];
rz(-pi) q[2];
rz(1.2346047) q[3];
sx q[3];
rz(-2.6264274) q[3];
sx q[3];
rz(-3.0807487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7427407) q[2];
sx q[2];
rz(-2.3336918) q[2];
sx q[2];
rz(0.74637949) q[2];
rz(-1.0910723) q[3];
sx q[3];
rz(-1.4336136) q[3];
sx q[3];
rz(2.5892042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.93254507) q[0];
sx q[0];
rz(-3.0952251) q[0];
sx q[0];
rz(2.5982502) q[0];
rz(-2.6047193) q[1];
sx q[1];
rz(-0.32506341) q[1];
sx q[1];
rz(-3.0184025) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42483271) q[0];
sx q[0];
rz(-0.93532978) q[0];
sx q[0];
rz(3.1034971) q[0];
rz(-pi) q[1];
rz(2.6307879) q[2];
sx q[2];
rz(-1.9419365) q[2];
sx q[2];
rz(-2.5696511) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1773273) q[1];
sx q[1];
rz(-0.96879241) q[1];
sx q[1];
rz(-0.89167562) q[1];
rz(-pi) q[2];
rz(0.79027883) q[3];
sx q[3];
rz(-2.1690282) q[3];
sx q[3];
rz(-1.7733396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3049551) q[2];
sx q[2];
rz(-1.7540437) q[2];
sx q[2];
rz(0.48509625) q[2];
rz(-1.1525611) q[3];
sx q[3];
rz(-1.3644812) q[3];
sx q[3];
rz(0.047957234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.7570067) q[0];
sx q[0];
rz(-2.204019) q[0];
sx q[0];
rz(2.6458929) q[0];
rz(-3.0777625) q[1];
sx q[1];
rz(-0.7494691) q[1];
sx q[1];
rz(3.0705423) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9793689) q[0];
sx q[0];
rz(-2.6155229) q[0];
sx q[0];
rz(-2.6499676) q[0];
rz(0.60511968) q[2];
sx q[2];
rz(-1.4265043) q[2];
sx q[2];
rz(0.52478204) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.38842216) q[1];
sx q[1];
rz(-2.8533397) q[1];
sx q[1];
rz(-0.99797319) q[1];
rz(0.17545731) q[3];
sx q[3];
rz(-2.6883295) q[3];
sx q[3];
rz(2.3153265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0157328) q[2];
sx q[2];
rz(-1.4733529) q[2];
sx q[2];
rz(-0.86137548) q[2];
rz(-1.8371948) q[3];
sx q[3];
rz(-0.574489) q[3];
sx q[3];
rz(-1.5931574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9651589) q[0];
sx q[0];
rz(-2.6506944) q[0];
sx q[0];
rz(1.4392256) q[0];
rz(-0.08055117) q[1];
sx q[1];
rz(-1.6702024) q[1];
sx q[1];
rz(-1.7505987) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2918562) q[0];
sx q[0];
rz(-2.448659) q[0];
sx q[0];
rz(-0.87397184) q[0];
rz(-1.5059696) q[2];
sx q[2];
rz(-2.8362084) q[2];
sx q[2];
rz(2.7750654) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.357482) q[1];
sx q[1];
rz(-2.3598089) q[1];
sx q[1];
rz(-2.2428721) q[1];
rz(-2.3426314) q[3];
sx q[3];
rz(-2.6332601) q[3];
sx q[3];
rz(1.5647183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3535658) q[2];
sx q[2];
rz(-1.2949508) q[2];
sx q[2];
rz(-0.12727748) q[2];
rz(2.2150529) q[3];
sx q[3];
rz(-2.8997731) q[3];
sx q[3];
rz(-1.658879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7552898) q[0];
sx q[0];
rz(-2.4918064) q[0];
sx q[0];
rz(1.5007098) q[0];
rz(2.3312148) q[1];
sx q[1];
rz(-2.2893298) q[1];
sx q[1];
rz(2.3694029) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9861575) q[0];
sx q[0];
rz(-0.86830189) q[0];
sx q[0];
rz(-2.2566354) q[0];
rz(-pi) q[1];
rz(-1.7842845) q[2];
sx q[2];
rz(-1.8597892) q[2];
sx q[2];
rz(-1.9910781) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.72550794) q[1];
sx q[1];
rz(-2.5920752) q[1];
sx q[1];
rz(-1.5931507) q[1];
x q[2];
rz(0.89086074) q[3];
sx q[3];
rz(-1.6107585) q[3];
sx q[3];
rz(-3.0639632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3810252) q[2];
sx q[2];
rz(-2.3367391) q[2];
sx q[2];
rz(-2.5123361) q[2];
rz(2.4106846) q[3];
sx q[3];
rz(-2.2980502) q[3];
sx q[3];
rz(-1.6985016) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44169852) q[0];
sx q[0];
rz(-0.48342539) q[0];
sx q[0];
rz(-0.40406686) q[0];
rz(2.7370257) q[1];
sx q[1];
rz(-1.7351983) q[1];
sx q[1];
rz(2.4529967) q[1];
rz(0.2307288) q[2];
sx q[2];
rz(-2.8525272) q[2];
sx q[2];
rz(2.0244103) q[2];
rz(-2.6062905) q[3];
sx q[3];
rz(-1.8284767) q[3];
sx q[3];
rz(-0.60093193) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
