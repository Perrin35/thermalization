OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.98439944) q[0];
sx q[0];
rz(2.0318883) q[0];
sx q[0];
rz(11.561031) q[0];
rz(0.73468626) q[1];
sx q[1];
rz(-2.7856196) q[1];
sx q[1];
rz(2.2495143) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5135554) q[0];
sx q[0];
rz(-1.6694067) q[0];
sx q[0];
rz(0.74752083) q[0];
rz(-0.14867013) q[2];
sx q[2];
rz(-0.97236247) q[2];
sx q[2];
rz(2.8658346) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0301117) q[1];
sx q[1];
rz(-1.4731543) q[1];
sx q[1];
rz(0.48508118) q[1];
rz(2.0964811) q[3];
sx q[3];
rz(-1.0307923) q[3];
sx q[3];
rz(1.857855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.72508183) q[2];
sx q[2];
rz(-1.9981013) q[2];
sx q[2];
rz(-1.4408646) q[2];
rz(2.1848988) q[3];
sx q[3];
rz(-2.0243093) q[3];
sx q[3];
rz(2.8982437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44593921) q[0];
sx q[0];
rz(-2.74701) q[0];
sx q[0];
rz(-1.12895) q[0];
rz(-0.24540643) q[1];
sx q[1];
rz(-1.3134198) q[1];
sx q[1];
rz(-2.479877) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45265439) q[0];
sx q[0];
rz(-2.1784349) q[0];
sx q[0];
rz(-1.2217962) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7956338) q[2];
sx q[2];
rz(-1.944146) q[2];
sx q[2];
rz(2.7527347) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7364517) q[1];
sx q[1];
rz(-1.5390522) q[1];
sx q[1];
rz(2.870138) q[1];
x q[2];
rz(1.8850967) q[3];
sx q[3];
rz(-2.8366025) q[3];
sx q[3];
rz(-1.5105607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5976065) q[2];
sx q[2];
rz(-0.49763766) q[2];
sx q[2];
rz(1.0428766) q[2];
rz(1.4526224) q[3];
sx q[3];
rz(-0.85113227) q[3];
sx q[3];
rz(-2.0378621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3747028) q[0];
sx q[0];
rz(-0.68793982) q[0];
sx q[0];
rz(-1.8324628) q[0];
rz(-2.6922928) q[1];
sx q[1];
rz(-0.6895014) q[1];
sx q[1];
rz(-1.4264533) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5912522) q[0];
sx q[0];
rz(-0.36446291) q[0];
sx q[0];
rz(-0.82247199) q[0];
rz(-pi) q[1];
rz(1.6108247) q[2];
sx q[2];
rz(-2.7187901) q[2];
sx q[2];
rz(-1.947248) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.92628463) q[1];
sx q[1];
rz(-2.4325772) q[1];
sx q[1];
rz(0.11911094) q[1];
x q[2];
rz(-3.0578311) q[3];
sx q[3];
rz(-1.3256095) q[3];
sx q[3];
rz(-0.29473588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7766777) q[2];
sx q[2];
rz(-1.1740351) q[2];
sx q[2];
rz(3.0752227) q[2];
rz(-2.5868609) q[3];
sx q[3];
rz(-1.4812255) q[3];
sx q[3];
rz(-1.6248645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1944815) q[0];
sx q[0];
rz(-2.8567061) q[0];
sx q[0];
rz(-0.77199212) q[0];
rz(-2.4783065) q[1];
sx q[1];
rz(-0.74378496) q[1];
sx q[1];
rz(3.0139121) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4467112) q[0];
sx q[0];
rz(-0.90411964) q[0];
sx q[0];
rz(-1.5645909) q[0];
rz(-pi) q[1];
rz(-2.588932) q[2];
sx q[2];
rz(-1.2465806) q[2];
sx q[2];
rz(0.49190258) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0952747) q[1];
sx q[1];
rz(-1.1409586) q[1];
sx q[1];
rz(-0.44989391) q[1];
rz(-pi) q[2];
rz(-1.7288307) q[3];
sx q[3];
rz(-1.1452066) q[3];
sx q[3];
rz(2.0142209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.53106236) q[2];
sx q[2];
rz(-2.1805111) q[2];
sx q[2];
rz(-2.5430211) q[2];
rz(0.5851723) q[3];
sx q[3];
rz(-1.7681311) q[3];
sx q[3];
rz(3.0857871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7912306) q[0];
sx q[0];
rz(-1.5086011) q[0];
sx q[0];
rz(2.1685261) q[0];
rz(0.19634253) q[1];
sx q[1];
rz(-2.0025496) q[1];
sx q[1];
rz(1.6729209) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7356883) q[0];
sx q[0];
rz(-0.87474147) q[0];
sx q[0];
rz(0.52175867) q[0];
x q[1];
rz(-0.25428793) q[2];
sx q[2];
rz(-2.0627705) q[2];
sx q[2];
rz(-0.51153431) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.65485146) q[1];
sx q[1];
rz(-1.2706725) q[1];
sx q[1];
rz(-2.9611118) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1096538) q[3];
sx q[3];
rz(-1.5232289) q[3];
sx q[3];
rz(2.4470234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.68957442) q[2];
sx q[2];
rz(-1.3268027) q[2];
sx q[2];
rz(0.16239521) q[2];
rz(1.1299805) q[3];
sx q[3];
rz(-0.88310784) q[3];
sx q[3];
rz(2.0067154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-0.75055403) q[0];
sx q[0];
rz(-2.6226608) q[0];
sx q[0];
rz(-3.0911875) q[0];
rz(0.16818908) q[1];
sx q[1];
rz(-1.5805809) q[1];
sx q[1];
rz(2.2680297) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45437434) q[0];
sx q[0];
rz(-0.66965398) q[0];
sx q[0];
rz(0.39910103) q[0];
rz(1.8040405) q[2];
sx q[2];
rz(-2.9961259) q[2];
sx q[2];
rz(1.5158397) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.62100121) q[1];
sx q[1];
rz(-1.1258003) q[1];
sx q[1];
rz(0.59153647) q[1];
rz(-pi) q[2];
rz(-1.3662947) q[3];
sx q[3];
rz(-1.2786153) q[3];
sx q[3];
rz(2.146221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0198589) q[2];
sx q[2];
rz(-0.21616082) q[2];
sx q[2];
rz(0.38062322) q[2];
rz(-2.5753283) q[3];
sx q[3];
rz(-1.7411313) q[3];
sx q[3];
rz(2.3316021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-0.64067632) q[0];
sx q[0];
rz(-2.930142) q[0];
sx q[0];
rz(-0.40263116) q[0];
rz(1.7598049) q[1];
sx q[1];
rz(-0.80843061) q[1];
sx q[1];
rz(3.0852539) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0005971) q[0];
sx q[0];
rz(-0.5282481) q[0];
sx q[0];
rz(-2.2597367) q[0];
rz(-2.0414786) q[2];
sx q[2];
rz(-1.6073445) q[2];
sx q[2];
rz(-1.4135897) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.67225918) q[1];
sx q[1];
rz(-1.0067098) q[1];
sx q[1];
rz(-3.1271598) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67752892) q[3];
sx q[3];
rz(-1.9871681) q[3];
sx q[3];
rz(0.29034055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.23400433) q[2];
sx q[2];
rz(-1.3434709) q[2];
sx q[2];
rz(2.6410356) q[2];
rz(-0.060700011) q[3];
sx q[3];
rz(-1.7874291) q[3];
sx q[3];
rz(-2.6589656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57182264) q[0];
sx q[0];
rz(-1.7800542) q[0];
sx q[0];
rz(-1.3609591) q[0];
rz(-0.743615) q[1];
sx q[1];
rz(-1.4185602) q[1];
sx q[1];
rz(-0.13882151) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25538975) q[0];
sx q[0];
rz(-0.74658827) q[0];
sx q[0];
rz(-2.3476178) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2108364) q[2];
sx q[2];
rz(-2.4413707) q[2];
sx q[2];
rz(2.8559358) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.19473361) q[1];
sx q[1];
rz(-1.3689965) q[1];
sx q[1];
rz(-1.9029593) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9337594) q[3];
sx q[3];
rz(-1.1869988) q[3];
sx q[3];
rz(1.919534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4678141) q[2];
sx q[2];
rz(-2.0970586) q[2];
sx q[2];
rz(-0.71411258) q[2];
rz(2.1821187) q[3];
sx q[3];
rz(-1.6474612) q[3];
sx q[3];
rz(2.1625471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4176843) q[0];
sx q[0];
rz(-0.67643106) q[0];
sx q[0];
rz(-0.12829256) q[0];
rz(-0.3872321) q[1];
sx q[1];
rz(-2.5435244) q[1];
sx q[1];
rz(1.3234352) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0664862) q[0];
sx q[0];
rz(-1.7597226) q[0];
sx q[0];
rz(0.099204258) q[0];
rz(-pi) q[1];
rz(0.88302112) q[2];
sx q[2];
rz(-2.3273914) q[2];
sx q[2];
rz(-1.778686) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1897498) q[1];
sx q[1];
rz(-1.7147168) q[1];
sx q[1];
rz(-0.014149498) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3714482) q[3];
sx q[3];
rz(-0.85271612) q[3];
sx q[3];
rz(1.9950641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6748176) q[2];
sx q[2];
rz(-1.520949) q[2];
sx q[2];
rz(-0.20469323) q[2];
rz(1.8681059) q[3];
sx q[3];
rz(-0.73015648) q[3];
sx q[3];
rz(-1.5761121) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2765163) q[0];
sx q[0];
rz(-2.5193494) q[0];
sx q[0];
rz(-0.21328558) q[0];
rz(-0.22008303) q[1];
sx q[1];
rz(-1.2525696) q[1];
sx q[1];
rz(1.3594886) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8187788) q[0];
sx q[0];
rz(-1.5867763) q[0];
sx q[0];
rz(-0.009601618) q[0];
rz(-pi) q[1];
rz(-1.8259117) q[2];
sx q[2];
rz(-1.1099225) q[2];
sx q[2];
rz(-0.67042353) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.90725856) q[1];
sx q[1];
rz(-1.9149019) q[1];
sx q[1];
rz(-1.0306148) q[1];
rz(-pi) q[2];
rz(-2.6892446) q[3];
sx q[3];
rz(-2.1901263) q[3];
sx q[3];
rz(1.0119708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6243817) q[2];
sx q[2];
rz(-1.5456079) q[2];
sx q[2];
rz(-0.33779302) q[2];
rz(-1.4126011) q[3];
sx q[3];
rz(-1.1510886) q[3];
sx q[3];
rz(-0.82144773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.74078858) q[0];
sx q[0];
rz(-2.3176226) q[0];
sx q[0];
rz(1.2299706) q[0];
rz(-2.7578655) q[1];
sx q[1];
rz(-1.4839254) q[1];
sx q[1];
rz(0.76212777) q[1];
rz(2.7376851) q[2];
sx q[2];
rz(-1.8941034) q[2];
sx q[2];
rz(-0.20878172) q[2];
rz(-1.1521793) q[3];
sx q[3];
rz(-0.55057303) q[3];
sx q[3];
rz(0.84921992) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
