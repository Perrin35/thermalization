OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.4322296) q[0];
sx q[0];
rz(-0.95786434) q[0];
sx q[0];
rz(0.14444484) q[0];
rz(0.56675178) q[1];
sx q[1];
rz(2.6161939) q[1];
sx q[1];
rz(8.4470196) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5241961) q[0];
sx q[0];
rz(-0.27271909) q[0];
sx q[0];
rz(-2.8196536) q[0];
rz(0.30774967) q[2];
sx q[2];
rz(-2.0648742) q[2];
sx q[2];
rz(-0.57927629) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6482918) q[1];
sx q[1];
rz(-0.19953218) q[1];
sx q[1];
rz(-0.7944016) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3542777) q[3];
sx q[3];
rz(-1.6978605) q[3];
sx q[3];
rz(-0.78391677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.67291659) q[2];
sx q[2];
rz(-1.1489931) q[2];
sx q[2];
rz(2.2093175) q[2];
rz(2.9428234) q[3];
sx q[3];
rz(-2.0310183) q[3];
sx q[3];
rz(-0.96536243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4345877) q[0];
sx q[0];
rz(-0.90536896) q[0];
sx q[0];
rz(0.36112753) q[0];
rz(-1.7065642) q[1];
sx q[1];
rz(-1.7838493) q[1];
sx q[1];
rz(0.8180058) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5404856) q[0];
sx q[0];
rz(-1.846082) q[0];
sx q[0];
rz(-0.40095532) q[0];
x q[1];
rz(-1.9423219) q[2];
sx q[2];
rz(-1.6852334) q[2];
sx q[2];
rz(2.3105846) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.2172912) q[1];
sx q[1];
rz(-2.0089307) q[1];
sx q[1];
rz(-0.89466722) q[1];
x q[2];
rz(-3.0844968) q[3];
sx q[3];
rz(-1.9340056) q[3];
sx q[3];
rz(2.925194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.25257418) q[2];
sx q[2];
rz(-2.694464) q[2];
sx q[2];
rz(1.7209631) q[2];
rz(-1.8255) q[3];
sx q[3];
rz(-0.75794739) q[3];
sx q[3];
rz(-0.38823286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4017568) q[0];
sx q[0];
rz(-0.49210423) q[0];
sx q[0];
rz(2.2706568) q[0];
rz(2.8254106) q[1];
sx q[1];
rz(-2.8600287) q[1];
sx q[1];
rz(0.2972163) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57097829) q[0];
sx q[0];
rz(-1.9750722) q[0];
sx q[0];
rz(-1.5017609) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77135135) q[2];
sx q[2];
rz(-1.0372247) q[2];
sx q[2];
rz(-0.16650621) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7797459) q[1];
sx q[1];
rz(-1.5231113) q[1];
sx q[1];
rz(-1.7576799) q[1];
rz(-pi) q[2];
x q[2];
rz(0.86621427) q[3];
sx q[3];
rz(-2.2598887) q[3];
sx q[3];
rz(0.44486526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3729942) q[2];
sx q[2];
rz(-2.3186389) q[2];
sx q[2];
rz(-2.7139943) q[2];
rz(-1.9528495) q[3];
sx q[3];
rz(-0.62994981) q[3];
sx q[3];
rz(2.6141613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7052085) q[0];
sx q[0];
rz(-1.5642865) q[0];
sx q[0];
rz(0.77392459) q[0];
rz(-2.4286843) q[1];
sx q[1];
rz(-2.1122825) q[1];
sx q[1];
rz(-0.68177044) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7101689) q[0];
sx q[0];
rz(-1.3623326) q[0];
sx q[0];
rz(-1.3039949) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8800456) q[2];
sx q[2];
rz(-2.5144858) q[2];
sx q[2];
rz(-2.4583465) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.53182787) q[1];
sx q[1];
rz(-1.0351287) q[1];
sx q[1];
rz(-0.51171724) q[1];
x q[2];
rz(-0.34337266) q[3];
sx q[3];
rz(-1.9150754) q[3];
sx q[3];
rz(1.4624649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1327847) q[2];
sx q[2];
rz(-0.71295732) q[2];
sx q[2];
rz(-2.183389) q[2];
rz(2.0751674) q[3];
sx q[3];
rz(-1.2833779) q[3];
sx q[3];
rz(-0.3796033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5565857) q[0];
sx q[0];
rz(-1.7946694) q[0];
sx q[0];
rz(0.56030309) q[0];
rz(2.141748) q[1];
sx q[1];
rz(-2.9381349) q[1];
sx q[1];
rz(-1.5195742) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0409659) q[0];
sx q[0];
rz(-2.0140411) q[0];
sx q[0];
rz(0.53951453) q[0];
rz(-0.6117758) q[2];
sx q[2];
rz(-1.099274) q[2];
sx q[2];
rz(-2.7437291) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6474364) q[1];
sx q[1];
rz(-1.7687706) q[1];
sx q[1];
rz(-2.3784749) q[1];
x q[2];
rz(3.0171379) q[3];
sx q[3];
rz(-0.96154562) q[3];
sx q[3];
rz(2.5973158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4524298) q[2];
sx q[2];
rz(-0.55183691) q[2];
sx q[2];
rz(-2.9186644) q[2];
rz(3.1068504) q[3];
sx q[3];
rz(-1.7545173) q[3];
sx q[3];
rz(0.071578659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.8054304) q[0];
sx q[0];
rz(-2.8283089) q[0];
sx q[0];
rz(1.0700595) q[0];
rz(1.7806212) q[1];
sx q[1];
rz(-2.762251) q[1];
sx q[1];
rz(-1.2840575) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59239913) q[0];
sx q[0];
rz(-2.1201029) q[0];
sx q[0];
rz(2.4248289) q[0];
x q[1];
rz(-0.27988866) q[2];
sx q[2];
rz(-1.1942689) q[2];
sx q[2];
rz(2.1044452) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2388873) q[1];
sx q[1];
rz(-1.1088015) q[1];
sx q[1];
rz(1.2483031) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0809903) q[3];
sx q[3];
rz(-0.59270699) q[3];
sx q[3];
rz(-0.61471516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.47508919) q[2];
sx q[2];
rz(-1.44375) q[2];
sx q[2];
rz(2.5040023) q[2];
rz(0.83667886) q[3];
sx q[3];
rz(-2.0357318) q[3];
sx q[3];
rz(1.7975413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5752983) q[0];
sx q[0];
rz(-2.7116382) q[0];
sx q[0];
rz(-2.5740525) q[0];
rz(-0.42770806) q[1];
sx q[1];
rz(-1.6141012) q[1];
sx q[1];
rz(2.2033851) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.339401) q[0];
sx q[0];
rz(-1.2656478) q[0];
sx q[0];
rz(-1.1423654) q[0];
rz(-pi) q[1];
x q[1];
rz(0.73056716) q[2];
sx q[2];
rz(-1.1267203) q[2];
sx q[2];
rz(0.72545746) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7792556) q[1];
sx q[1];
rz(-1.2149097) q[1];
sx q[1];
rz(-2.2625655) q[1];
rz(-0.70721831) q[3];
sx q[3];
rz(-1.8054609) q[3];
sx q[3];
rz(-2.9747687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.87970916) q[2];
sx q[2];
rz(-1.9677013) q[2];
sx q[2];
rz(-1.3860469) q[2];
rz(-1.322768) q[3];
sx q[3];
rz(-1.1498007) q[3];
sx q[3];
rz(3.1125606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8742074) q[0];
sx q[0];
rz(-0.33339849) q[0];
sx q[0];
rz(-1.7077131) q[0];
rz(1.2738312) q[1];
sx q[1];
rz(-2.0055983) q[1];
sx q[1];
rz(-0.83126718) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87330504) q[0];
sx q[0];
rz(-1.9718861) q[0];
sx q[0];
rz(-1.7307161) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93584658) q[2];
sx q[2];
rz(-1.4668462) q[2];
sx q[2];
rz(-0.39633358) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.89325209) q[1];
sx q[1];
rz(-2.2447531) q[1];
sx q[1];
rz(2.2158951) q[1];
x q[2];
rz(-1.7940815) q[3];
sx q[3];
rz(-1.7956453) q[3];
sx q[3];
rz(-2.4527578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5635809) q[2];
sx q[2];
rz(-1.7687904) q[2];
sx q[2];
rz(-2.1772299) q[2];
rz(-1.1635121) q[3];
sx q[3];
rz(-2.5301299) q[3];
sx q[3];
rz(0.65892974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37373856) q[0];
sx q[0];
rz(-2.4066194) q[0];
sx q[0];
rz(0.97737616) q[0];
rz(-1.3865698) q[1];
sx q[1];
rz(-1.3061378) q[1];
sx q[1];
rz(2.0358553) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39767299) q[0];
sx q[0];
rz(-2.0154675) q[0];
sx q[0];
rz(-1.4414653) q[0];
x q[1];
rz(0.90412037) q[2];
sx q[2];
rz(-1.6534001) q[2];
sx q[2];
rz(0.77043515) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9165009) q[1];
sx q[1];
rz(-1.8019925) q[1];
sx q[1];
rz(-0.1901615) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8589784) q[3];
sx q[3];
rz(-1.9087221) q[3];
sx q[3];
rz(-1.5087138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.608312) q[2];
sx q[2];
rz(-0.38485843) q[2];
sx q[2];
rz(-2.9917955) q[2];
rz(1.7685361) q[3];
sx q[3];
rz(-1.7356197) q[3];
sx q[3];
rz(1.8201374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49366429) q[0];
sx q[0];
rz(-0.51769185) q[0];
sx q[0];
rz(-3.0143484) q[0];
rz(1.6607025) q[1];
sx q[1];
rz(-0.36247411) q[1];
sx q[1];
rz(-0.1677992) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0509335) q[0];
sx q[0];
rz(-1.5899842) q[0];
sx q[0];
rz(-0.030254342) q[0];
rz(-pi) q[1];
rz(2.9211505) q[2];
sx q[2];
rz(-2.5219005) q[2];
sx q[2];
rz(0.1534136) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.64393109) q[1];
sx q[1];
rz(-1.5776458) q[1];
sx q[1];
rz(0.58845206) q[1];
rz(0.12331788) q[3];
sx q[3];
rz(-1.6880369) q[3];
sx q[3];
rz(-0.58837147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.1148791) q[2];
sx q[2];
rz(-0.9388887) q[2];
sx q[2];
rz(-2.3804469) q[2];
rz(-0.090027697) q[3];
sx q[3];
rz(-1.0031676) q[3];
sx q[3];
rz(-2.1910523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5933843) q[0];
sx q[0];
rz(-1.9822639) q[0];
sx q[0];
rz(-0.32250861) q[0];
rz(-2.7535915) q[1];
sx q[1];
rz(-1.7419659) q[1];
sx q[1];
rz(2.3566125) q[1];
rz(-0.8969174) q[2];
sx q[2];
rz(-2.2581836) q[2];
sx q[2];
rz(2.7473292) q[2];
rz(2.9705863) q[3];
sx q[3];
rz(-2.1033559) q[3];
sx q[3];
rz(-0.81567473) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];