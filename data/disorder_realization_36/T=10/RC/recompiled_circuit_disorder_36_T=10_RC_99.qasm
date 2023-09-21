OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.49178034) q[0];
sx q[0];
rz(3.427504) q[0];
sx q[0];
rz(8.9094845) q[0];
rz(4.4858785) q[1];
sx q[1];
rz(2.9872515) q[1];
sx q[1];
rz(6.8607688) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2429457) q[0];
sx q[0];
rz(-2.4624914) q[0];
sx q[0];
rz(-2.8773017) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.53519997) q[2];
sx q[2];
rz(-2.0173965) q[2];
sx q[2];
rz(1.610178) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9211728) q[1];
sx q[1];
rz(-2.6487659) q[1];
sx q[1];
rz(2.1652031) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3732784) q[3];
sx q[3];
rz(-2.3294805) q[3];
sx q[3];
rz(2.1384359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.704533) q[2];
sx q[2];
rz(-1.5960863) q[2];
sx q[2];
rz(0.68721592) q[2];
rz(1.0152738) q[3];
sx q[3];
rz(-1.3736558) q[3];
sx q[3];
rz(0.12250531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9706443) q[0];
sx q[0];
rz(-2.0630554) q[0];
sx q[0];
rz(1.2600391) q[0];
rz(2.1353703) q[1];
sx q[1];
rz(-0.99199122) q[1];
sx q[1];
rz(2.2959183) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.082367912) q[0];
sx q[0];
rz(-2.0031643) q[0];
sx q[0];
rz(2.5655377) q[0];
rz(-3.0078366) q[2];
sx q[2];
rz(-1.4417366) q[2];
sx q[2];
rz(-1.2226906) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.37296346) q[1];
sx q[1];
rz(-1.0781204) q[1];
sx q[1];
rz(-1.782934) q[1];
rz(-pi) q[2];
rz(-2.1917079) q[3];
sx q[3];
rz(-1.1511027) q[3];
sx q[3];
rz(1.3148395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3530897) q[2];
sx q[2];
rz(-2.916009) q[2];
sx q[2];
rz(0.4804002) q[2];
rz(-1.7885615) q[3];
sx q[3];
rz(-2.0856817) q[3];
sx q[3];
rz(1.9539179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1903494) q[0];
sx q[0];
rz(-0.23878637) q[0];
sx q[0];
rz(2.3685266) q[0];
rz(-3.0103325) q[1];
sx q[1];
rz(-1.2845598) q[1];
sx q[1];
rz(2.0551596) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13585424) q[0];
sx q[0];
rz(-1.5724626) q[0];
sx q[0];
rz(-1.1093596) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92347446) q[2];
sx q[2];
rz(-1.662865) q[2];
sx q[2];
rz(1.1084686) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.31049) q[1];
sx q[1];
rz(-1.5655787) q[1];
sx q[1];
rz(-2.8527841) q[1];
x q[2];
rz(-1.1658737) q[3];
sx q[3];
rz(-1.4110663) q[3];
sx q[3];
rz(0.83042849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1290258) q[2];
sx q[2];
rz(-1.7029224) q[2];
sx q[2];
rz(-1.3712937) q[2];
rz(0.38315547) q[3];
sx q[3];
rz(-1.2569191) q[3];
sx q[3];
rz(0.80254054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4521769) q[0];
sx q[0];
rz(-1.2503662) q[0];
sx q[0];
rz(0.048359811) q[0];
rz(0.16391779) q[1];
sx q[1];
rz(-2.7719031) q[1];
sx q[1];
rz(1.6960467) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4500344) q[0];
sx q[0];
rz(-1.7507179) q[0];
sx q[0];
rz(0.66288373) q[0];
rz(-pi) q[1];
rz(2.7518026) q[2];
sx q[2];
rz(-0.23332694) q[2];
sx q[2];
rz(-1.1167655) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9963035) q[1];
sx q[1];
rz(-2.9128296) q[1];
sx q[1];
rz(2.8537675) q[1];
rz(-pi) q[2];
rz(-1.0501782) q[3];
sx q[3];
rz(-1.8225267) q[3];
sx q[3];
rz(-2.7057735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0754898) q[2];
sx q[2];
rz(-1.3572071) q[2];
sx q[2];
rz(-1.0243105) q[2];
rz(1.6131489) q[3];
sx q[3];
rz(-1.5214835) q[3];
sx q[3];
rz(0.23322341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(1.4399453) q[0];
sx q[0];
rz(-2.3174536) q[0];
sx q[0];
rz(1.2874999) q[0];
rz(0.31907407) q[1];
sx q[1];
rz(-1.5417475) q[1];
sx q[1];
rz(0.85420001) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4695278) q[0];
sx q[0];
rz(-0.85692353) q[0];
sx q[0];
rz(0.010849997) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6536077) q[2];
sx q[2];
rz(-2.2099566) q[2];
sx q[2];
rz(-1.1346863) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2687208) q[1];
sx q[1];
rz(-0.99209058) q[1];
sx q[1];
rz(0.89923664) q[1];
rz(-pi) q[2];
x q[2];
rz(0.79980005) q[3];
sx q[3];
rz(-0.721867) q[3];
sx q[3];
rz(-1.1218027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1896818) q[2];
sx q[2];
rz(-0.59331912) q[2];
sx q[2];
rz(-2.5642776) q[2];
rz(-2.632085) q[3];
sx q[3];
rz(-0.43764344) q[3];
sx q[3];
rz(0.90665162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1381056) q[0];
sx q[0];
rz(-2.0697937) q[0];
sx q[0];
rz(-0.072120897) q[0];
rz(-2.0347319) q[1];
sx q[1];
rz(-2.6289584) q[1];
sx q[1];
rz(-3.0153826) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70996767) q[0];
sx q[0];
rz(-2.9266848) q[0];
sx q[0];
rz(0.14847319) q[0];
x q[1];
rz(2.7563165) q[2];
sx q[2];
rz(-2.3284973) q[2];
sx q[2];
rz(-2.3713881) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5050161) q[1];
sx q[1];
rz(-0.86537433) q[1];
sx q[1];
rz(2.8935863) q[1];
rz(-1.8993127) q[3];
sx q[3];
rz(-1.7528755) q[3];
sx q[3];
rz(1.8329221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2386027) q[2];
sx q[2];
rz(-0.1923407) q[2];
sx q[2];
rz(0.77511707) q[2];
rz(2.3136247) q[3];
sx q[3];
rz(-0.29100806) q[3];
sx q[3];
rz(2.0194139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.5489952) q[0];
sx q[0];
rz(-0.42625517) q[0];
sx q[0];
rz(-0.098408498) q[0];
rz(-1.9495643) q[1];
sx q[1];
rz(-1.8076618) q[1];
sx q[1];
rz(0.55955204) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7647117) q[0];
sx q[0];
rz(-1.6983713) q[0];
sx q[0];
rz(2.5401831) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6227116) q[2];
sx q[2];
rz(-1.8318818) q[2];
sx q[2];
rz(-1.7920997) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0911078) q[1];
sx q[1];
rz(-0.15639601) q[1];
sx q[1];
rz(0.97538235) q[1];
x q[2];
rz(2.5114602) q[3];
sx q[3];
rz(-1.7297941) q[3];
sx q[3];
rz(0.4160479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6825535) q[2];
sx q[2];
rz(-1.8457396) q[2];
sx q[2];
rz(0.34379488) q[2];
rz(0.5665468) q[3];
sx q[3];
rz(-0.44851258) q[3];
sx q[3];
rz(-0.47376537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.375181) q[0];
sx q[0];
rz(-1.3324998) q[0];
sx q[0];
rz(0.73076105) q[0];
rz(0.14239755) q[1];
sx q[1];
rz(-1.2700894) q[1];
sx q[1];
rz(2.2699845) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5430785) q[0];
sx q[0];
rz(-3.0656272) q[0];
sx q[0];
rz(3.024858) q[0];
rz(-pi) q[1];
rz(-1.7224738) q[2];
sx q[2];
rz(-1.1693839) q[2];
sx q[2];
rz(0.36827189) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1069113) q[1];
sx q[1];
rz(-0.92202631) q[1];
sx q[1];
rz(-2.0295063) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1400847) q[3];
sx q[3];
rz(-0.9730556) q[3];
sx q[3];
rz(2.7261581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.72620755) q[2];
sx q[2];
rz(-2.0724847) q[2];
sx q[2];
rz(-2.0020206) q[2];
rz(-1.4987882) q[3];
sx q[3];
rz(-0.39396861) q[3];
sx q[3];
rz(-0.9128226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20294872) q[0];
sx q[0];
rz(-1.4511755) q[0];
sx q[0];
rz(-1.2217481) q[0];
rz(2.9755759) q[1];
sx q[1];
rz(-1.8211726) q[1];
sx q[1];
rz(1.5244012) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0696408) q[0];
sx q[0];
rz(-1.6531634) q[0];
sx q[0];
rz(0.029043341) q[0];
rz(-pi) q[1];
rz(2.7828214) q[2];
sx q[2];
rz(-2.6371187) q[2];
sx q[2];
rz(-2.3141253) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7414654) q[1];
sx q[1];
rz(-1.8997846) q[1];
sx q[1];
rz(-0.16727438) q[1];
rz(2.7306261) q[3];
sx q[3];
rz(-0.68115679) q[3];
sx q[3];
rz(2.4364803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1200072) q[2];
sx q[2];
rz(-1.6759796) q[2];
sx q[2];
rz(-0.35153708) q[2];
rz(-2.0848138) q[3];
sx q[3];
rz(-0.52962279) q[3];
sx q[3];
rz(2.3969011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.4979424) q[0];
sx q[0];
rz(-2.239776) q[0];
sx q[0];
rz(1.3051916) q[0];
rz(2.7611043) q[1];
sx q[1];
rz(-1.0419798) q[1];
sx q[1];
rz(-2.8881853) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8587592) q[0];
sx q[0];
rz(-1.4677605) q[0];
sx q[0];
rz(-1.9673002) q[0];
rz(-pi) q[1];
rz(-2.017574) q[2];
sx q[2];
rz(-0.3728711) q[2];
sx q[2];
rz(0.98137059) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3897755) q[1];
sx q[1];
rz(-2.7174065) q[1];
sx q[1];
rz(0.61702375) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1425584) q[3];
sx q[3];
rz(-1.0034475) q[3];
sx q[3];
rz(2.464307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.59166756) q[2];
sx q[2];
rz(-0.88576907) q[2];
sx q[2];
rz(0.5029451) q[2];
rz(0.89899603) q[3];
sx q[3];
rz(-1.8476202) q[3];
sx q[3];
rz(-1.1635273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9713365) q[0];
sx q[0];
rz(-1.6032871) q[0];
sx q[0];
rz(0.26300318) q[0];
rz(2.4304216) q[1];
sx q[1];
rz(-1.0881337) q[1];
sx q[1];
rz(1.7137391) q[1];
rz(0.74264991) q[2];
sx q[2];
rz(-2.7682318) q[2];
sx q[2];
rz(-2.9329185) q[2];
rz(2.3861804) q[3];
sx q[3];
rz(-2.073954) q[3];
sx q[3];
rz(3.0975773) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];