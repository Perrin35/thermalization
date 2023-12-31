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
rz(-2.5748409) q[1];
sx q[1];
rz(-2.6161939) q[1];
sx q[1];
rz(2.1638343) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35740556) q[0];
sx q[0];
rz(-1.4854684) q[0];
sx q[0];
rz(-2.8822495) q[0];
rz(-pi) q[1];
rz(1.0563645) q[2];
sx q[2];
rz(-1.8407514) q[2];
sx q[2];
rz(1.1411238) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.830476) q[1];
sx q[1];
rz(-1.7101354) q[1];
sx q[1];
rz(1.4275101) q[1];
rz(-pi) q[2];
rz(-1.3917189) q[3];
sx q[3];
rz(-2.350051) q[3];
sx q[3];
rz(-2.4812428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.67291659) q[2];
sx q[2];
rz(-1.9925995) q[2];
sx q[2];
rz(-0.93227512) q[2];
rz(2.9428234) q[3];
sx q[3];
rz(-2.0310183) q[3];
sx q[3];
rz(-0.96536243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
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
rz(-1.4350285) q[1];
sx q[1];
rz(-1.3577434) q[1];
sx q[1];
rz(-2.3235869) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0571787) q[0];
sx q[0];
rz(-1.955843) q[0];
sx q[0];
rz(1.8684698) q[0];
x q[1];
rz(1.9423219) q[2];
sx q[2];
rz(-1.6852334) q[2];
sx q[2];
rz(0.83100806) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.2172912) q[1];
sx q[1];
rz(-1.1326619) q[1];
sx q[1];
rz(2.2469254) q[1];
rz(-pi) q[2];
rz(0.057095842) q[3];
sx q[3];
rz(-1.207587) q[3];
sx q[3];
rz(0.2163987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8890185) q[2];
sx q[2];
rz(-0.44712862) q[2];
sx q[2];
rz(1.7209631) q[2];
rz(-1.3160926) q[3];
sx q[3];
rz(-0.75794739) q[3];
sx q[3];
rz(0.38823286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7398359) q[0];
sx q[0];
rz(-0.49210423) q[0];
sx q[0];
rz(0.87093583) q[0];
rz(-0.31618205) q[1];
sx q[1];
rz(-2.8600287) q[1];
sx q[1];
rz(0.2972163) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0270099) q[0];
sx q[0];
rz(-1.5073338) q[0];
sx q[0];
rz(-0.40513904) q[0];
rz(-pi) q[1];
rz(-0.77135135) q[2];
sx q[2];
rz(-2.1043679) q[2];
sx q[2];
rz(0.16650621) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.19993648) q[1];
sx q[1];
rz(-1.3841277) q[1];
sx q[1];
rz(0.048528683) q[1];
rz(0.86621427) q[3];
sx q[3];
rz(-0.88170393) q[3];
sx q[3];
rz(-0.44486526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3729942) q[2];
sx q[2];
rz(-2.3186389) q[2];
sx q[2];
rz(2.7139943) q[2];
rz(1.9528495) q[3];
sx q[3];
rz(-0.62994981) q[3];
sx q[3];
rz(-2.6141613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7052085) q[0];
sx q[0];
rz(-1.5642865) q[0];
sx q[0];
rz(0.77392459) q[0];
rz(-0.71290839) q[1];
sx q[1];
rz(-2.1122825) q[1];
sx q[1];
rz(-2.4598222) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9457152) q[0];
sx q[0];
rz(-1.8316852) q[0];
sx q[0];
rz(2.9257141) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21708023) q[2];
sx q[2];
rz(-0.97765572) q[2];
sx q[2];
rz(2.8341688) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.53182787) q[1];
sx q[1];
rz(-1.0351287) q[1];
sx q[1];
rz(2.6298754) q[1];
rz(-pi) q[2];
rz(1.2069615) q[3];
sx q[3];
rz(-1.893265) q[3];
sx q[3];
rz(0.22842562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1327847) q[2];
sx q[2];
rz(-2.4286353) q[2];
sx q[2];
rz(0.95820367) q[2];
rz(-2.0751674) q[3];
sx q[3];
rz(-1.2833779) q[3];
sx q[3];
rz(-2.7619894) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.585007) q[0];
sx q[0];
rz(-1.3469232) q[0];
sx q[0];
rz(2.5812896) q[0];
rz(0.99984461) q[1];
sx q[1];
rz(-2.9381349) q[1];
sx q[1];
rz(-1.6220185) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1006267) q[0];
sx q[0];
rz(-2.0140411) q[0];
sx q[0];
rz(0.53951453) q[0];
x q[1];
rz(2.4155596) q[2];
sx q[2];
rz(-0.75349977) q[2];
sx q[2];
rz(-1.7475278) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.8738185) q[1];
sx q[1];
rz(-0.78332892) q[1];
sx q[1];
rz(2.8591213) q[1];
rz(-pi) q[2];
rz(-2.1836957) q[3];
sx q[3];
rz(-1.6727722) q[3];
sx q[3];
rz(-2.0436055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4524298) q[2];
sx q[2];
rz(-0.55183691) q[2];
sx q[2];
rz(-0.2229283) q[2];
rz(3.1068504) q[3];
sx q[3];
rz(-1.7545173) q[3];
sx q[3];
rz(-3.070014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3361622) q[0];
sx q[0];
rz(-2.8283089) q[0];
sx q[0];
rz(1.0700595) q[0];
rz(1.7806212) q[1];
sx q[1];
rz(-0.37934163) q[1];
sx q[1];
rz(-1.8575352) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5491935) q[0];
sx q[0];
rz(-2.1201029) q[0];
sx q[0];
rz(0.71676371) q[0];
x q[1];
rz(-2.1806296) q[2];
sx q[2];
rz(-0.46513882) q[2];
sx q[2];
rz(-2.7679408) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2594086) q[1];
sx q[1];
rz(-2.5849197) q[1];
sx q[1];
rz(-0.56682079) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0345801) q[3];
sx q[3];
rz(-1.3048733) q[3];
sx q[3];
rz(1.3724316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6665035) q[2];
sx q[2];
rz(-1.6978426) q[2];
sx q[2];
rz(-0.63759032) q[2];
rz(-2.3049138) q[3];
sx q[3];
rz(-2.0357318) q[3];
sx q[3];
rz(-1.3440514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-0.93820757) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81323775) q[0];
sx q[0];
rz(-0.52044808) q[0];
sx q[0];
rz(2.2195199) q[0];
x q[1];
rz(-2.5222048) q[2];
sx q[2];
rz(-2.3084547) q[2];
sx q[2];
rz(-2.7433861) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.36233703) q[1];
sx q[1];
rz(-1.2149097) q[1];
sx q[1];
rz(2.2625655) q[1];
x q[2];
rz(1.8754962) q[3];
sx q[3];
rz(-2.2548171) q[3];
sx q[3];
rz(-1.6001493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.87970916) q[2];
sx q[2];
rz(-1.9677013) q[2];
sx q[2];
rz(1.7555457) q[2];
rz(1.8188247) q[3];
sx q[3];
rz(-1.1498007) q[3];
sx q[3];
rz(-0.02903207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26738527) q[0];
sx q[0];
rz(-0.33339849) q[0];
sx q[0];
rz(-1.4338795) q[0];
rz(1.2738312) q[1];
sx q[1];
rz(-1.1359943) q[1];
sx q[1];
rz(-2.3103255) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48152637) q[0];
sx q[0];
rz(-0.43018451) q[0];
sx q[0];
rz(2.7823886) q[0];
x q[1];
rz(0.12886329) q[2];
sx q[2];
rz(-0.93981758) q[2];
sx q[2];
rz(1.250759) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.23849328) q[1];
sx q[1];
rz(-2.0598754) q[1];
sx q[1];
rz(0.78519435) q[1];
rz(0.76922272) q[3];
sx q[3];
rz(-2.8260494) q[3];
sx q[3];
rz(-0.10570082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5780118) q[2];
sx q[2];
rz(-1.7687904) q[2];
sx q[2];
rz(-0.9643628) q[2];
rz(-1.1635121) q[3];
sx q[3];
rz(-2.5301299) q[3];
sx q[3];
rz(-2.4826629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7678541) q[0];
sx q[0];
rz(-2.4066194) q[0];
sx q[0];
rz(2.1642165) q[0];
rz(1.3865698) q[1];
sx q[1];
rz(-1.3061378) q[1];
sx q[1];
rz(-2.0358553) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2290105) q[0];
sx q[0];
rz(-1.6874896) q[0];
sx q[0];
rz(-0.44793655) q[0];
rz(-pi) q[1];
rz(1.703891) q[2];
sx q[2];
rz(-2.4705952) q[2];
sx q[2];
rz(-2.4457096) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.92377201) q[1];
sx q[1];
rz(-0.2982699) q[1];
sx q[1];
rz(-0.89426269) q[1];
rz(-pi) q[2];
x q[2];
rz(0.35131485) q[3];
sx q[3];
rz(-1.8422541) q[3];
sx q[3];
rz(-0.16004496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.608312) q[2];
sx q[2];
rz(-2.7567342) q[2];
sx q[2];
rz(-0.14979714) q[2];
rz(1.3730565) q[3];
sx q[3];
rz(-1.7356197) q[3];
sx q[3];
rz(-1.8201374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6479284) q[0];
sx q[0];
rz(-2.6239008) q[0];
sx q[0];
rz(-3.0143484) q[0];
rz(1.6607025) q[1];
sx q[1];
rz(-2.7791185) q[1];
sx q[1];
rz(0.1677992) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0566933) q[0];
sx q[0];
rz(-0.035824422) q[0];
sx q[0];
rz(2.5762659) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4160412) q[2];
sx q[2];
rz(-2.1733279) q[2];
sx q[2];
rz(0.11520152) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.93712902) q[1];
sx q[1];
rz(-2.5531054) q[1];
sx q[1];
rz(-0.012339331) q[1];
x q[2];
rz(0.12331788) q[3];
sx q[3];
rz(-1.4535558) q[3];
sx q[3];
rz(0.58837147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.026713513) q[2];
sx q[2];
rz(-2.202704) q[2];
sx q[2];
rz(0.76114571) q[2];
rz(0.090027697) q[3];
sx q[3];
rz(-1.0031676) q[3];
sx q[3];
rz(-0.95054039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5933843) q[0];
sx q[0];
rz(-1.9822639) q[0];
sx q[0];
rz(-0.32250861) q[0];
rz(-0.38800115) q[1];
sx q[1];
rz(-1.3996268) q[1];
sx q[1];
rz(-0.7849801) q[1];
rz(2.491647) q[2];
sx q[2];
rz(-2.2192498) q[2];
sx q[2];
rz(0.50512209) q[2];
rz(-1.0317867) q[3];
sx q[3];
rz(-1.7179334) q[3];
sx q[3];
rz(-2.4739305) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
