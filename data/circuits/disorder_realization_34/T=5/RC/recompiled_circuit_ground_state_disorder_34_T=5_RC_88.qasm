OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7393957) q[0];
sx q[0];
rz(4.0255044) q[0];
sx q[0];
rz(8.5691353) q[0];
rz(2.9698676) q[1];
sx q[1];
rz(-3.0260234) q[1];
sx q[1];
rz(0.56245437) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5006258) q[0];
sx q[0];
rz(-1.62407) q[0];
sx q[0];
rz(-1.9356273) q[0];
rz(-pi) q[1];
rz(2.9859574) q[2];
sx q[2];
rz(-2.0726207) q[2];
sx q[2];
rz(1.4841472) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.89002883) q[1];
sx q[1];
rz(-1.152413) q[1];
sx q[1];
rz(-2.0189925) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80339949) q[3];
sx q[3];
rz(-2.6028742) q[3];
sx q[3];
rz(-2.5634457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.90613753) q[2];
sx q[2];
rz(-1.9635341) q[2];
sx q[2];
rz(-0.79929024) q[2];
rz(0.47131395) q[3];
sx q[3];
rz(-2.2282232) q[3];
sx q[3];
rz(-0.93588626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1021295) q[0];
sx q[0];
rz(-1.3251745) q[0];
sx q[0];
rz(1.348173) q[0];
rz(-1.7680291) q[1];
sx q[1];
rz(-1.1496239) q[1];
sx q[1];
rz(1.9893533) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3593529) q[0];
sx q[0];
rz(-0.77762654) q[0];
sx q[0];
rz(1.1471143) q[0];
x q[1];
rz(0.94448467) q[2];
sx q[2];
rz(-0.99530333) q[2];
sx q[2];
rz(-1.4552417) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.5785529) q[1];
sx q[1];
rz(-0.7078979) q[1];
sx q[1];
rz(-2.8824174) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3186734) q[3];
sx q[3];
rz(-0.69861087) q[3];
sx q[3];
rz(2.5320092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3422602) q[2];
sx q[2];
rz(-2.7228184) q[2];
sx q[2];
rz(-0.93117923) q[2];
rz(-3.0350507) q[3];
sx q[3];
rz(-1.1139161) q[3];
sx q[3];
rz(0.60025269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0837285) q[0];
sx q[0];
rz(-1.7513542) q[0];
sx q[0];
rz(-0.51783836) q[0];
rz(0.88874108) q[1];
sx q[1];
rz(-2.4338212) q[1];
sx q[1];
rz(-0.4471561) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5260552) q[0];
sx q[0];
rz(-1.3096022) q[0];
sx q[0];
rz(1.5047856) q[0];
rz(-pi) q[1];
rz(-2.1250399) q[2];
sx q[2];
rz(-0.92281658) q[2];
sx q[2];
rz(-2.3584443) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0886317) q[1];
sx q[1];
rz(-1.4354032) q[1];
sx q[1];
rz(0.9779344) q[1];
rz(-1.8497883) q[3];
sx q[3];
rz(-2.0800262) q[3];
sx q[3];
rz(-1.4470456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7670224) q[2];
sx q[2];
rz(-2.3775103) q[2];
sx q[2];
rz(-0.53263295) q[2];
rz(0.24752188) q[3];
sx q[3];
rz(-0.73740021) q[3];
sx q[3];
rz(-2.1029162) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5175051) q[0];
sx q[0];
rz(-0.085260304) q[0];
sx q[0];
rz(3.0349773) q[0];
rz(-0.34635776) q[1];
sx q[1];
rz(-2.296591) q[1];
sx q[1];
rz(1.1603629) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4100944) q[0];
sx q[0];
rz(-2.838352) q[0];
sx q[0];
rz(-0.61654697) q[0];
x q[1];
rz(0.71765064) q[2];
sx q[2];
rz(-2.1776878) q[2];
sx q[2];
rz(-0.90536149) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6667156) q[1];
sx q[1];
rz(-1.2704388) q[1];
sx q[1];
rz(1.9431861) q[1];
x q[2];
rz(2.0959693) q[3];
sx q[3];
rz(-1.9672638) q[3];
sx q[3];
rz(3.0399655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13566636) q[2];
sx q[2];
rz(-0.29834193) q[2];
sx q[2];
rz(0.076210991) q[2];
rz(-0.58586079) q[3];
sx q[3];
rz(-1.1548837) q[3];
sx q[3];
rz(-1.3425672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12014408) q[0];
sx q[0];
rz(-1.3055389) q[0];
sx q[0];
rz(-0.35010499) q[0];
rz(-0.95589751) q[1];
sx q[1];
rz(-1.8616385) q[1];
sx q[1];
rz(1.9116481) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0405925) q[0];
sx q[0];
rz(-0.44389566) q[0];
sx q[0];
rz(-2.1568265) q[0];
rz(-3.1051738) q[2];
sx q[2];
rz(-1.7360188) q[2];
sx q[2];
rz(0.18319229) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1175244) q[1];
sx q[1];
rz(-1.9000016) q[1];
sx q[1];
rz(-0.19964053) q[1];
rz(-1.2754945) q[3];
sx q[3];
rz(-2.5341138) q[3];
sx q[3];
rz(2.5584084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2436287) q[2];
sx q[2];
rz(-2.882759) q[2];
sx q[2];
rz(1.492307) q[2];
rz(2.7367075) q[3];
sx q[3];
rz(-2.3758774) q[3];
sx q[3];
rz(-2.305472) q[3];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1355857) q[0];
sx q[0];
rz(-1.0772935) q[0];
sx q[0];
rz(-0.58240044) q[0];
rz(0.68663418) q[1];
sx q[1];
rz(-1.735894) q[1];
sx q[1];
rz(1.2581717) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95520407) q[0];
sx q[0];
rz(-2.0058647) q[0];
sx q[0];
rz(-3.0032936) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7863726) q[2];
sx q[2];
rz(-1.0632916) q[2];
sx q[2];
rz(1.2025646) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4542089) q[1];
sx q[1];
rz(-1.6404248) q[1];
sx q[1];
rz(1.5147665) q[1];
rz(-pi) q[2];
rz(-1.889713) q[3];
sx q[3];
rz(-0.61781672) q[3];
sx q[3];
rz(-1.8998673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.76576343) q[2];
sx q[2];
rz(-1.9932237) q[2];
sx q[2];
rz(2.1235535) q[2];
rz(-0.24614075) q[3];
sx q[3];
rz(-1.3956416) q[3];
sx q[3];
rz(-1.5181946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70596424) q[0];
sx q[0];
rz(-0.80393296) q[0];
sx q[0];
rz(0.58018082) q[0];
rz(0.14353453) q[1];
sx q[1];
rz(-2.6519471) q[1];
sx q[1];
rz(-0.20763436) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.62791) q[0];
sx q[0];
rz(-1.2259036) q[0];
sx q[0];
rz(-2.3570127) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34822627) q[2];
sx q[2];
rz(-1.019632) q[2];
sx q[2];
rz(1.7334235) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5112111) q[1];
sx q[1];
rz(-2.2464754) q[1];
sx q[1];
rz(0.33501321) q[1];
x q[2];
rz(2.8945699) q[3];
sx q[3];
rz(-2.7678856) q[3];
sx q[3];
rz(1.6173687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84270728) q[2];
sx q[2];
rz(-1.2879813) q[2];
sx q[2];
rz(-2.5353954) q[2];
rz(-0.68743622) q[3];
sx q[3];
rz(-2.0076553) q[3];
sx q[3];
rz(0.70639759) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6523478) q[0];
sx q[0];
rz(-0.027996538) q[0];
sx q[0];
rz(-1.0580753) q[0];
rz(3.0319013) q[1];
sx q[1];
rz(-1.1166162) q[1];
sx q[1];
rz(-1.6995957) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083695166) q[0];
sx q[0];
rz(-2.6287615) q[0];
sx q[0];
rz(-2.4714241) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2953561) q[2];
sx q[2];
rz(-1.0862964) q[2];
sx q[2];
rz(0.4415919) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.40934) q[1];
sx q[1];
rz(-2.2870018) q[1];
sx q[1];
rz(-0.90998896) q[1];
rz(2.9061716) q[3];
sx q[3];
rz(-1.306433) q[3];
sx q[3];
rz(0.89371577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.66120061) q[2];
sx q[2];
rz(-0.19442393) q[2];
sx q[2];
rz(-1.2083758) q[2];
rz(-2.4783573) q[3];
sx q[3];
rz(-1.6845208) q[3];
sx q[3];
rz(-1.9251582) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97063589) q[0];
sx q[0];
rz(-0.96681505) q[0];
sx q[0];
rz(-0.092967689) q[0];
rz(-1.2804821) q[1];
sx q[1];
rz(-2.4130776) q[1];
sx q[1];
rz(-3.0063937) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4358035) q[0];
sx q[0];
rz(-3.0717141) q[0];
sx q[0];
rz(0.70894928) q[0];
rz(-2.8835758) q[2];
sx q[2];
rz(-1.4224367) q[2];
sx q[2];
rz(0.49000636) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1982806) q[1];
sx q[1];
rz(-1.3168646) q[1];
sx q[1];
rz(-1.5487681) q[1];
rz(-pi) q[2];
x q[2];
rz(0.074196176) q[3];
sx q[3];
rz(-1.7230986) q[3];
sx q[3];
rz(1.1049113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4325503) q[2];
sx q[2];
rz(-1.21864) q[2];
sx q[2];
rz(2.3700628) q[2];
rz(-2.6500474) q[3];
sx q[3];
rz(-1.8621657) q[3];
sx q[3];
rz(0.5947203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.58525697) q[0];
sx q[0];
rz(-0.62194967) q[0];
sx q[0];
rz(1.1248032) q[0];
rz(1.8428165) q[1];
sx q[1];
rz(-0.61538428) q[1];
sx q[1];
rz(0.76464701) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1021834) q[0];
sx q[0];
rz(-1.4847857) q[0];
sx q[0];
rz(-0.10552222) q[0];
rz(-2.6420399) q[2];
sx q[2];
rz(-0.95328125) q[2];
sx q[2];
rz(-2.2411186) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5950039) q[1];
sx q[1];
rz(-1.4440795) q[1];
sx q[1];
rz(1.7751883) q[1];
rz(2.5705757) q[3];
sx q[3];
rz(-0.63435508) q[3];
sx q[3];
rz(-2.3650124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3760959) q[2];
sx q[2];
rz(-1.8509879) q[2];
sx q[2];
rz(0.20720227) q[2];
rz(-2.1616705) q[3];
sx q[3];
rz(-1.1332952) q[3];
sx q[3];
rz(-2.9492212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1205263) q[0];
sx q[0];
rz(-1.4561894) q[0];
sx q[0];
rz(-0.86984632) q[0];
rz(-0.5008685) q[1];
sx q[1];
rz(-0.23575467) q[1];
sx q[1];
rz(-2.2101319) q[1];
rz(-1.4516713) q[2];
sx q[2];
rz(-1.7075734) q[2];
sx q[2];
rz(1.3592958) q[2];
rz(2.5475827) q[3];
sx q[3];
rz(-2.8566711) q[3];
sx q[3];
rz(2.2950238) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
