OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4797526) q[0];
sx q[0];
rz(-2.2979484) q[0];
sx q[0];
rz(2.9736829) q[0];
rz(1.1711988) q[1];
sx q[1];
rz(3.436915) q[1];
sx q[1];
rz(9.480939) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6457155) q[0];
sx q[0];
rz(-2.0464532) q[0];
sx q[0];
rz(-0.23947421) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.216823) q[2];
sx q[2];
rz(-1.3999108) q[2];
sx q[2];
rz(0.31121635) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15966283) q[1];
sx q[1];
rz(-2.5291981) q[1];
sx q[1];
rz(-2.0484522) q[1];
x q[2];
rz(-2.9763016) q[3];
sx q[3];
rz(-2.878302) q[3];
sx q[3];
rz(0.96112448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.37796676) q[2];
sx q[2];
rz(-0.28181919) q[2];
sx q[2];
rz(0.4326694) q[2];
rz(-1.9487322) q[3];
sx q[3];
rz(-1.2377219) q[3];
sx q[3];
rz(0.38309923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52779657) q[0];
sx q[0];
rz(-0.48848099) q[0];
sx q[0];
rz(-1.8288076) q[0];
rz(2.9361172) q[1];
sx q[1];
rz(-2.165129) q[1];
sx q[1];
rz(-1.1516494) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5713455) q[0];
sx q[0];
rz(-0.76128188) q[0];
sx q[0];
rz(2.0061357) q[0];
rz(2.4685681) q[2];
sx q[2];
rz(-0.7012127) q[2];
sx q[2];
rz(2.6075624) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9559905) q[1];
sx q[1];
rz(-0.62528175) q[1];
sx q[1];
rz(2.37466) q[1];
rz(2.6395256) q[3];
sx q[3];
rz(-1.9126529) q[3];
sx q[3];
rz(1.7934007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.0097222086) q[2];
sx q[2];
rz(-1.6513609) q[2];
sx q[2];
rz(-2.9197664) q[2];
rz(-0.37718537) q[3];
sx q[3];
rz(-2.714034) q[3];
sx q[3];
rz(0.77243531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.31056988) q[0];
sx q[0];
rz(-3.0492058) q[0];
sx q[0];
rz(-3.1047399) q[0];
rz(2.316078) q[1];
sx q[1];
rz(-1.3157536) q[1];
sx q[1];
rz(-3.085014) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7650334) q[0];
sx q[0];
rz(-2.3925836) q[0];
sx q[0];
rz(0.88699938) q[0];
rz(1.4372196) q[2];
sx q[2];
rz(-1.3165511) q[2];
sx q[2];
rz(0.1453407) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.024836) q[1];
sx q[1];
rz(-0.36326888) q[1];
sx q[1];
rz(-0.58961745) q[1];
rz(0.64485456) q[3];
sx q[3];
rz(-1.2554902) q[3];
sx q[3];
rz(-0.23526084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8686707) q[2];
sx q[2];
rz(-1.4977095) q[2];
sx q[2];
rz(2.2154714) q[2];
rz(2.5849294) q[3];
sx q[3];
rz(-0.29354468) q[3];
sx q[3];
rz(1.042897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91519231) q[0];
sx q[0];
rz(-2.4214348) q[0];
sx q[0];
rz(0.74209374) q[0];
rz(2.0023951) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(-2.6779968) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.050042) q[0];
sx q[0];
rz(-0.76489641) q[0];
sx q[0];
rz(-2.0043623) q[0];
rz(-pi) q[1];
rz(-2.8088403) q[2];
sx q[2];
rz(-1.5126265) q[2];
sx q[2];
rz(1.0156877) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0193034) q[1];
sx q[1];
rz(-2.3483854) q[1];
sx q[1];
rz(-1.8122458) q[1];
rz(0.47834088) q[3];
sx q[3];
rz(-2.0739177) q[3];
sx q[3];
rz(2.2174045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3670369) q[2];
sx q[2];
rz(-2.98428) q[2];
sx q[2];
rz(-0.049499361) q[2];
rz(3.0130623) q[3];
sx q[3];
rz(-1.5481719) q[3];
sx q[3];
rz(-3.0226504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0054935) q[0];
sx q[0];
rz(-2.7042784) q[0];
sx q[0];
rz(-2.8438925) q[0];
rz(-0.4822576) q[1];
sx q[1];
rz(-0.75459701) q[1];
sx q[1];
rz(2.1972426) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56086841) q[0];
sx q[0];
rz(-1.7797911) q[0];
sx q[0];
rz(-0.046037721) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3659533) q[2];
sx q[2];
rz(-1.2795942) q[2];
sx q[2];
rz(1.8679384) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9152865) q[1];
sx q[1];
rz(-1.5729135) q[1];
sx q[1];
rz(1.5096942) q[1];
rz(1.827042) q[3];
sx q[3];
rz(-1.1503997) q[3];
sx q[3];
rz(-2.8375569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.258761) q[2];
sx q[2];
rz(-1.0927039) q[2];
sx q[2];
rz(-0.10822254) q[2];
rz(-3.1392858) q[3];
sx q[3];
rz(-1.5294411) q[3];
sx q[3];
rz(2.8172857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7047983) q[0];
sx q[0];
rz(-0.44610766) q[0];
sx q[0];
rz(-0.58445245) q[0];
rz(-2.2553518) q[1];
sx q[1];
rz(-0.61683547) q[1];
sx q[1];
rz(3.086673) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9960105) q[0];
sx q[0];
rz(-1.8390363) q[0];
sx q[0];
rz(0.085573816) q[0];
rz(-pi) q[1];
rz(1.5149649) q[2];
sx q[2];
rz(-2.8179114) q[2];
sx q[2];
rz(-1.8813546) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0262895) q[1];
sx q[1];
rz(-2.4541828) q[1];
sx q[1];
rz(-2.5251212) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6782126) q[3];
sx q[3];
rz(-1.0372835) q[3];
sx q[3];
rz(2.2802071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3871258) q[2];
sx q[2];
rz(-0.12067623) q[2];
sx q[2];
rz(-1.0167271) q[2];
rz(-2.5975442) q[3];
sx q[3];
rz(-0.35990158) q[3];
sx q[3];
rz(1.8090766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6761557) q[0];
sx q[0];
rz(-0.98452079) q[0];
sx q[0];
rz(0.28453919) q[0];
rz(2.1971205) q[1];
sx q[1];
rz(-1.9453134) q[1];
sx q[1];
rz(2.231266) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.401424) q[0];
sx q[0];
rz(-1.5059885) q[0];
sx q[0];
rz(0.054697371) q[0];
x q[1];
rz(0.34142999) q[2];
sx q[2];
rz(-1.176398) q[2];
sx q[2];
rz(2.9801126) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.42156223) q[1];
sx q[1];
rz(-1.7517462) q[1];
sx q[1];
rz(-2.6471495) q[1];
rz(-2.4207553) q[3];
sx q[3];
rz(-2.0257054) q[3];
sx q[3];
rz(2.7618559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3900782) q[2];
sx q[2];
rz(-3.0780767) q[2];
sx q[2];
rz(-0.92203036) q[2];
rz(2.5743124) q[3];
sx q[3];
rz(-1.4138979) q[3];
sx q[3];
rz(-2.1218307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89408016) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(0.12938736) q[0];
rz(0.63240504) q[1];
sx q[1];
rz(-1.0267195) q[1];
sx q[1];
rz(2.8410889) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6209517) q[0];
sx q[0];
rz(-1.2004939) q[0];
sx q[0];
rz(-2.305549) q[0];
rz(0.6408765) q[2];
sx q[2];
rz(-0.87101988) q[2];
sx q[2];
rz(-1.4027632) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.064425163) q[1];
sx q[1];
rz(-2.7556813) q[1];
sx q[1];
rz(-0.47972958) q[1];
rz(-pi) q[2];
rz(1.9037876) q[3];
sx q[3];
rz(-0.78562842) q[3];
sx q[3];
rz(-1.9612519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.58632103) q[2];
sx q[2];
rz(-0.99493146) q[2];
sx q[2];
rz(-0.78197455) q[2];
rz(0.55109763) q[3];
sx q[3];
rz(-1.3827773) q[3];
sx q[3];
rz(0.15792318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5683811) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(0.12776275) q[0];
rz(-2.5993775) q[1];
sx q[1];
rz(-2.1844889) q[1];
sx q[1];
rz(0.75884563) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42620537) q[0];
sx q[0];
rz(-1.9949159) q[0];
sx q[0];
rz(-0.018652648) q[0];
rz(-pi) q[1];
rz(-2.5326469) q[2];
sx q[2];
rz(-1.7296089) q[2];
sx q[2];
rz(-1.2283404) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.54770494) q[1];
sx q[1];
rz(-0.87915671) q[1];
sx q[1];
rz(-1.7734852) q[1];
x q[2];
rz(0.90518732) q[3];
sx q[3];
rz(-1.5593411) q[3];
sx q[3];
rz(-0.36597914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0163429) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(-0.49003595) q[2];
rz(1.4222493) q[3];
sx q[3];
rz(-1.2079206) q[3];
sx q[3];
rz(2.0786044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35995099) q[0];
sx q[0];
rz(-0.61976969) q[0];
sx q[0];
rz(-0.075335659) q[0];
rz(-2.244859) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(2.5316701) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2595554) q[0];
sx q[0];
rz(-1.8543188) q[0];
sx q[0];
rz(-0.77772227) q[0];
rz(-1.2953193) q[2];
sx q[2];
rz(-1.9343978) q[2];
sx q[2];
rz(1.3678577) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6590609) q[1];
sx q[1];
rz(-0.7809124) q[1];
sx q[1];
rz(-0.34300967) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0718594) q[3];
sx q[3];
rz(-2.2483629) q[3];
sx q[3];
rz(0.60615218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.23218368) q[2];
sx q[2];
rz(-0.82513088) q[2];
sx q[2];
rz(-2.4278736) q[2];
rz(0.37832007) q[3];
sx q[3];
rz(-0.49351966) q[3];
sx q[3];
rz(-0.87987125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.338035) q[0];
sx q[0];
rz(-1.150158) q[0];
sx q[0];
rz(-1.5858571) q[0];
rz(0.65080416) q[1];
sx q[1];
rz(-1.6497859) q[1];
sx q[1];
rz(-0.12129687) q[1];
rz(-0.030608721) q[2];
sx q[2];
rz(-1.3749214) q[2];
sx q[2];
rz(2.2236852) q[2];
rz(0.25272947) q[3];
sx q[3];
rz(-1.1128294) q[3];
sx q[3];
rz(2.2313234) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
