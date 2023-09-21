OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.5819117) q[0];
sx q[0];
rz(-2.0547325) q[0];
sx q[0];
rz(1.342919) q[0];
rz(1.5496594) q[1];
sx q[1];
rz(-0.078443371) q[1];
sx q[1];
rz(-0.6426386) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0138577) q[0];
sx q[0];
rz(-0.53984261) q[0];
sx q[0];
rz(-0.054245754) q[0];
rz(-pi) q[1];
rz(-0.3354934) q[2];
sx q[2];
rz(-0.53288424) q[2];
sx q[2];
rz(0.83836183) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1145775) q[1];
sx q[1];
rz(-1.6903965) q[1];
sx q[1];
rz(0.12160614) q[1];
x q[2];
rz(-0.23407614) q[3];
sx q[3];
rz(-1.5353068) q[3];
sx q[3];
rz(1.8518098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6671483) q[2];
sx q[2];
rz(-0.64626226) q[2];
sx q[2];
rz(1.2791963) q[2];
rz(0.71875087) q[3];
sx q[3];
rz(-1.5712534) q[3];
sx q[3];
rz(0.32354245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6779125) q[0];
sx q[0];
rz(-1.7872515) q[0];
sx q[0];
rz(-0.98051488) q[0];
rz(-2.9837043) q[1];
sx q[1];
rz(-1.6442559) q[1];
sx q[1];
rz(0.78871361) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65608998) q[0];
sx q[0];
rz(-0.1682818) q[0];
sx q[0];
rz(0.82505723) q[0];
rz(-pi) q[1];
rz(2.8854495) q[2];
sx q[2];
rz(-1.6015341) q[2];
sx q[2];
rz(-0.61402938) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.5309108) q[1];
sx q[1];
rz(-2.1149939) q[1];
sx q[1];
rz(-1.9287841) q[1];
rz(-pi) q[2];
rz(-2.2867145) q[3];
sx q[3];
rz(-2.2023871) q[3];
sx q[3];
rz(-2.7090493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.366189) q[2];
sx q[2];
rz(-2.0663694) q[2];
sx q[2];
rz(0.186084) q[2];
rz(2.4880593) q[3];
sx q[3];
rz(-1.7553522) q[3];
sx q[3];
rz(-0.11793605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2114975) q[0];
sx q[0];
rz(-2.1989172) q[0];
sx q[0];
rz(0.63013664) q[0];
rz(-3.0139626) q[1];
sx q[1];
rz(-0.6487414) q[1];
sx q[1];
rz(-0.72174597) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9741309) q[0];
sx q[0];
rz(-2.2142017) q[0];
sx q[0];
rz(-1.8882303) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0314991) q[2];
sx q[2];
rz(-2.2113872) q[2];
sx q[2];
rz(5.0355807e-05) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0112146) q[1];
sx q[1];
rz(-1.6092669) q[1];
sx q[1];
rz(2.6973666) q[1];
rz(-pi) q[2];
rz(2.8053022) q[3];
sx q[3];
rz(-0.97067562) q[3];
sx q[3];
rz(2.1093413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7599941) q[2];
sx q[2];
rz(-0.95278946) q[2];
sx q[2];
rz(0.90908137) q[2];
rz(-2.6233853) q[3];
sx q[3];
rz(-2.0776694) q[3];
sx q[3];
rz(-0.28809965) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5575314) q[0];
sx q[0];
rz(-2.4052305) q[0];
sx q[0];
rz(-3.0990565) q[0];
rz(2.361239) q[1];
sx q[1];
rz(-2.6413554) q[1];
sx q[1];
rz(-2.3775878) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4064179) q[0];
sx q[0];
rz(-1.5582726) q[0];
sx q[0];
rz(-0.89212117) q[0];
rz(-pi) q[1];
rz(2.9742083) q[2];
sx q[2];
rz(-0.65132729) q[2];
sx q[2];
rz(1.1812047) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2149787) q[1];
sx q[1];
rz(-1.1264631) q[1];
sx q[1];
rz(-0.39919969) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66391151) q[3];
sx q[3];
rz(-1.9466562) q[3];
sx q[3];
rz(-2.1223048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2970695) q[2];
sx q[2];
rz(-1.2244747) q[2];
sx q[2];
rz(-2.6085473) q[2];
rz(2.879203) q[3];
sx q[3];
rz(-0.636594) q[3];
sx q[3];
rz(-2.6791402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2247291) q[0];
sx q[0];
rz(-0.30039766) q[0];
sx q[0];
rz(-1.8977144) q[0];
rz(-0.22661701) q[1];
sx q[1];
rz(-2.367327) q[1];
sx q[1];
rz(-0.40333834) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6809083) q[0];
sx q[0];
rz(-1.4049238) q[0];
sx q[0];
rz(-2.9537863) q[0];
x q[1];
rz(2.7062347) q[2];
sx q[2];
rz(-0.73409664) q[2];
sx q[2];
rz(-1.8045319) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7008608) q[1];
sx q[1];
rz(-1.4167538) q[1];
sx q[1];
rz(-0.53149077) q[1];
rz(0.63286085) q[3];
sx q[3];
rz(-0.78213464) q[3];
sx q[3];
rz(2.6485505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8205745) q[2];
sx q[2];
rz(-2.0901188) q[2];
sx q[2];
rz(-0.78249758) q[2];
rz(2.0292422) q[3];
sx q[3];
rz(-0.4370884) q[3];
sx q[3];
rz(-2.1330244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.430442) q[0];
sx q[0];
rz(-0.40200457) q[0];
sx q[0];
rz(-0.45561403) q[0];
rz(3.0420711) q[1];
sx q[1];
rz(-1.1385304) q[1];
sx q[1];
rz(-0.0064370357) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9329487) q[0];
sx q[0];
rz(-2.3017578) q[0];
sx q[0];
rz(1.2293925) q[0];
rz(-pi) q[1];
rz(-0.10613425) q[2];
sx q[2];
rz(-1.1960293) q[2];
sx q[2];
rz(-0.29667621) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.16033123) q[1];
sx q[1];
rz(-0.48109522) q[1];
sx q[1];
rz(1.0465924) q[1];
rz(-pi) q[2];
rz(0.84900093) q[3];
sx q[3];
rz(-0.95522049) q[3];
sx q[3];
rz(1.4403696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.36879888) q[2];
sx q[2];
rz(-1.6364748) q[2];
sx q[2];
rz(2.2951365) q[2];
rz(2.1438697) q[3];
sx q[3];
rz(-2.3322451) q[3];
sx q[3];
rz(-1.458582) q[3];
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
rz(pi/2) q[3];
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
rz(3.0434175) q[0];
sx q[0];
rz(-2.257453) q[0];
sx q[0];
rz(-3.085882) q[0];
rz(-2.3588691) q[1];
sx q[1];
rz(-2.0109773) q[1];
sx q[1];
rz(-1.1605211) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29138628) q[0];
sx q[0];
rz(-2.0178447) q[0];
sx q[0];
rz(-2.7544751) q[0];
x q[1];
rz(-2.0918526) q[2];
sx q[2];
rz(-1.1546635) q[2];
sx q[2];
rz(0.001948826) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6351663) q[1];
sx q[1];
rz(-2.4392358) q[1];
sx q[1];
rz(0.96794767) q[1];
rz(2.3485687) q[3];
sx q[3];
rz(-2.6245955) q[3];
sx q[3];
rz(1.6263863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0760076) q[2];
sx q[2];
rz(-2.2202754) q[2];
sx q[2];
rz(0.78545061) q[2];
rz(-0.75585946) q[3];
sx q[3];
rz(-1.8463585) q[3];
sx q[3];
rz(2.7261962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-1.8287559) q[0];
sx q[0];
rz(-0.090933196) q[0];
sx q[0];
rz(1.1428517) q[0];
rz(-1.8354592) q[1];
sx q[1];
rz(-1.9530692) q[1];
sx q[1];
rz(-0.41608861) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5160822) q[0];
sx q[0];
rz(-1.7689546) q[0];
sx q[0];
rz(-1.0869736) q[0];
rz(1.847319) q[2];
sx q[2];
rz(-2.6575436) q[2];
sx q[2];
rz(0.36177847) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0033274) q[1];
sx q[1];
rz(-1.7308373) q[1];
sx q[1];
rz(0.40469594) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0307426) q[3];
sx q[3];
rz(-1.1006315) q[3];
sx q[3];
rz(-1.6925616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0351506) q[2];
sx q[2];
rz(-1.687655) q[2];
sx q[2];
rz(2.8184334) q[2];
rz(0.20251003) q[3];
sx q[3];
rz(-1.2812307) q[3];
sx q[3];
rz(-0.57730738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.768196) q[0];
sx q[0];
rz(-0.62921262) q[0];
sx q[0];
rz(1.9966104) q[0];
rz(-1.9454983) q[1];
sx q[1];
rz(-0.15592608) q[1];
sx q[1];
rz(0.51913613) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3648758) q[0];
sx q[0];
rz(-1.4380699) q[0];
sx q[0];
rz(1.9507292) q[0];
rz(1.8066508) q[2];
sx q[2];
rz(-1.7867076) q[2];
sx q[2];
rz(-2.4049135) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0172826) q[1];
sx q[1];
rz(-2.3553684) q[1];
sx q[1];
rz(-3.0148274) q[1];
rz(-1.395123) q[3];
sx q[3];
rz(-1.3364949) q[3];
sx q[3];
rz(-2.7845886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8273932) q[2];
sx q[2];
rz(-1.9284356) q[2];
sx q[2];
rz(0.040977565) q[2];
rz(-0.86769062) q[3];
sx q[3];
rz(-2.4882081) q[3];
sx q[3];
rz(0.51122558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39559078) q[0];
sx q[0];
rz(-2.262291) q[0];
sx q[0];
rz(1.6145153) q[0];
rz(-1.4279667) q[1];
sx q[1];
rz(-0.36179301) q[1];
sx q[1];
rz(1.6428927) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4157383) q[0];
sx q[0];
rz(-3.0634974) q[0];
sx q[0];
rz(-1.8321091) q[0];
rz(1.2730359) q[2];
sx q[2];
rz(-0.16212633) q[2];
sx q[2];
rz(-2.8043384) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5012706) q[1];
sx q[1];
rz(-0.81201279) q[1];
sx q[1];
rz(-2.0444319) q[1];
x q[2];
rz(1.4478217) q[3];
sx q[3];
rz(-1.7719367) q[3];
sx q[3];
rz(-2.5589383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3830118) q[2];
sx q[2];
rz(-2.0643533) q[2];
sx q[2];
rz(-0.14653462) q[2];
rz(0.81418973) q[3];
sx q[3];
rz(-0.79629961) q[3];
sx q[3];
rz(-2.7248785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9183337) q[0];
sx q[0];
rz(-1.2467361) q[0];
sx q[0];
rz(0.99714739) q[0];
rz(-1.2659484) q[1];
sx q[1];
rz(-1.0026149) q[1];
sx q[1];
rz(1.2276585) q[1];
rz(3.0284991) q[2];
sx q[2];
rz(-1.8086776) q[2];
sx q[2];
rz(2.0247731) q[2];
rz(-0.72600611) q[3];
sx q[3];
rz(-0.86961679) q[3];
sx q[3];
rz(1.1179954) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];