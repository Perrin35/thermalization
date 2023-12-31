OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.52656093) q[0];
sx q[0];
rz(-2.5685413) q[0];
sx q[0];
rz(-0.84258643) q[0];
rz(-1.0358345) q[1];
sx q[1];
rz(-2.0422715) q[1];
sx q[1];
rz(1.6834747) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9487171) q[0];
sx q[0];
rz(-2.2964381) q[0];
sx q[0];
rz(-1.2851508) q[0];
rz(-pi) q[1];
rz(0.61383944) q[2];
sx q[2];
rz(-1.5547353) q[2];
sx q[2];
rz(1.8526555) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.38274256) q[1];
sx q[1];
rz(-1.9470125) q[1];
sx q[1];
rz(1.43169) q[1];
rz(-2.9620142) q[3];
sx q[3];
rz(-2.539733) q[3];
sx q[3];
rz(1.4048525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.52790102) q[2];
sx q[2];
rz(-2.1353022) q[2];
sx q[2];
rz(-2.9620985) q[2];
rz(1.9159296) q[3];
sx q[3];
rz(-1.7951199) q[3];
sx q[3];
rz(-0.82204449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3935788) q[0];
sx q[0];
rz(-0.8809692) q[0];
sx q[0];
rz(-2.8161312) q[0];
rz(1.356396) q[1];
sx q[1];
rz(-2.0929095) q[1];
sx q[1];
rz(-1.9869841) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3852859) q[0];
sx q[0];
rz(-3.1161838) q[0];
sx q[0];
rz(-0.73861648) q[0];
rz(1.0500533) q[2];
sx q[2];
rz(-0.69486952) q[2];
sx q[2];
rz(0.75887242) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.11101152) q[1];
sx q[1];
rz(-1.6554553) q[1];
sx q[1];
rz(2.3762523) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96315893) q[3];
sx q[3];
rz(-2.8249486) q[3];
sx q[3];
rz(0.39099993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6894199) q[2];
sx q[2];
rz(-1.2499115) q[2];
sx q[2];
rz(0.88341218) q[2];
rz(2.6702821) q[3];
sx q[3];
rz(-1.4383957) q[3];
sx q[3];
rz(-2.3538891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31323355) q[0];
sx q[0];
rz(-1.6468843) q[0];
sx q[0];
rz(1.5154243) q[0];
rz(0.60107636) q[1];
sx q[1];
rz(-0.54769146) q[1];
sx q[1];
rz(2.0498958) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1356782) q[0];
sx q[0];
rz(-1.0251097) q[0];
sx q[0];
rz(-0.90555993) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0513564) q[2];
sx q[2];
rz(-0.72548496) q[2];
sx q[2];
rz(-1.5067593) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6669238) q[1];
sx q[1];
rz(-0.83819929) q[1];
sx q[1];
rz(1.0679507) q[1];
x q[2];
rz(-3.0136209) q[3];
sx q[3];
rz(-1.5068753) q[3];
sx q[3];
rz(1.909006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8213356) q[2];
sx q[2];
rz(-0.50575033) q[2];
sx q[2];
rz(-0.88095218) q[2];
rz(-1.7679924) q[3];
sx q[3];
rz(-1.526984) q[3];
sx q[3];
rz(-2.1239471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83051935) q[0];
sx q[0];
rz(-1.3922465) q[0];
sx q[0];
rz(-2.7048892) q[0];
rz(0.23315915) q[1];
sx q[1];
rz(-1.2522839) q[1];
sx q[1];
rz(2.8312347) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096075637) q[0];
sx q[0];
rz(-1.1928416) q[0];
sx q[0];
rz(1.4287352) q[0];
x q[1];
rz(-0.68508673) q[2];
sx q[2];
rz(-1.473046) q[2];
sx q[2];
rz(-3.0435261) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0879285) q[1];
sx q[1];
rz(-1.2128608) q[1];
sx q[1];
rz(0.019836516) q[1];
x q[2];
rz(-0.035590812) q[3];
sx q[3];
rz(-1.6944052) q[3];
sx q[3];
rz(1.3328758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0115396) q[2];
sx q[2];
rz(-0.71802846) q[2];
sx q[2];
rz(1.0774353) q[2];
rz(3.0854026) q[3];
sx q[3];
rz(-2.5037933) q[3];
sx q[3];
rz(-1.5475387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9524277) q[0];
sx q[0];
rz(-1.0452894) q[0];
sx q[0];
rz(2.8919343) q[0];
rz(1.5646308) q[1];
sx q[1];
rz(-0.77762929) q[1];
sx q[1];
rz(2.2713984) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10660431) q[0];
sx q[0];
rz(-0.37811324) q[0];
sx q[0];
rz(-0.58017054) q[0];
rz(2.9179847) q[2];
sx q[2];
rz(-2.4325271) q[2];
sx q[2];
rz(0.86415926) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.2493077) q[1];
sx q[1];
rz(-1.2586437) q[1];
sx q[1];
rz(-2.409163) q[1];
x q[2];
rz(1.8364041) q[3];
sx q[3];
rz(-0.57816539) q[3];
sx q[3];
rz(-2.2418914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.27328086) q[2];
sx q[2];
rz(-1.8153278) q[2];
sx q[2];
rz(-0.67374054) q[2];
rz(-0.30361787) q[3];
sx q[3];
rz(-1.9165336) q[3];
sx q[3];
rz(-1.822086) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7917787) q[0];
sx q[0];
rz(-2.204201) q[0];
sx q[0];
rz(0.2579903) q[0];
rz(-2.7164283) q[1];
sx q[1];
rz(-2.185967) q[1];
sx q[1];
rz(-1.4917096) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0817889) q[0];
sx q[0];
rz(-0.44133082) q[0];
sx q[0];
rz(0.23131891) q[0];
rz(0.67955534) q[2];
sx q[2];
rz(-1.2365885) q[2];
sx q[2];
rz(-2.213775) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3867492) q[1];
sx q[1];
rz(-1.0068839) q[1];
sx q[1];
rz(0.90653231) q[1];
rz(-pi) q[2];
rz(-1.8317354) q[3];
sx q[3];
rz(-2.76537) q[3];
sx q[3];
rz(0.46686831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.012718) q[2];
sx q[2];
rz(-2.1779163) q[2];
sx q[2];
rz(-0.091726124) q[2];
rz(0.84364676) q[3];
sx q[3];
rz(-2.1618312) q[3];
sx q[3];
rz(2.2475524) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82350746) q[0];
sx q[0];
rz(-1.894269) q[0];
sx q[0];
rz(-2.7303625) q[0];
rz(2.2757018) q[1];
sx q[1];
rz(-2.829268) q[1];
sx q[1];
rz(-0.033989865) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2798529) q[0];
sx q[0];
rz(-1.4089157) q[0];
sx q[0];
rz(-2.3539761) q[0];
rz(2.2074239) q[2];
sx q[2];
rz(-1.120943) q[2];
sx q[2];
rz(2.5069782) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9719203) q[1];
sx q[1];
rz(-1.6213413) q[1];
sx q[1];
rz(-2.2488942) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61492413) q[3];
sx q[3];
rz(-1.4498386) q[3];
sx q[3];
rz(1.8254335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6035446) q[2];
sx q[2];
rz(-2.5431583) q[2];
sx q[2];
rz(-2.2650488) q[2];
rz(2.792568) q[3];
sx q[3];
rz(-1.2004431) q[3];
sx q[3];
rz(-0.14311895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2440764) q[0];
sx q[0];
rz(-1.709047) q[0];
sx q[0];
rz(-2.7476655) q[0];
rz(0.36755964) q[1];
sx q[1];
rz(-1.3840679) q[1];
sx q[1];
rz(1.4454909) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.53795) q[0];
sx q[0];
rz(-1.1374439) q[0];
sx q[0];
rz(0.005714697) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5567899) q[2];
sx q[2];
rz(-2.1091166) q[2];
sx q[2];
rz(-2.1540097) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6519424) q[1];
sx q[1];
rz(-1.4133269) q[1];
sx q[1];
rz(-0.55649346) q[1];
rz(0.94533841) q[3];
sx q[3];
rz(-1.006554) q[3];
sx q[3];
rz(-2.4953147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2945071) q[2];
sx q[2];
rz(-0.89035788) q[2];
sx q[2];
rz(-0.40714804) q[2];
rz(-1.6242705) q[3];
sx q[3];
rz(-1.9842691) q[3];
sx q[3];
rz(-0.24967641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3354934) q[0];
sx q[0];
rz(-2.6265916) q[0];
sx q[0];
rz(1.2517713) q[0];
rz(0.66954008) q[1];
sx q[1];
rz(-1.1839097) q[1];
sx q[1];
rz(-2.8318185) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71000242) q[0];
sx q[0];
rz(-2.6012523) q[0];
sx q[0];
rz(-0.24775981) q[0];
x q[1];
rz(-1.3938815) q[2];
sx q[2];
rz(-1.5614911) q[2];
sx q[2];
rz(-2.6975346) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.27948353) q[1];
sx q[1];
rz(-2.9276507) q[1];
sx q[1];
rz(1.8474384) q[1];
x q[2];
rz(0.36848948) q[3];
sx q[3];
rz(-2.512305) q[3];
sx q[3];
rz(0.027241782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.70242515) q[2];
sx q[2];
rz(-0.71321407) q[2];
sx q[2];
rz(1.2072198) q[2];
rz(2.1045945) q[3];
sx q[3];
rz(-1.8959277) q[3];
sx q[3];
rz(0.65565482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0913775) q[0];
sx q[0];
rz(-1.8176879) q[0];
sx q[0];
rz(1.2058831) q[0];
rz(-2.5559015) q[1];
sx q[1];
rz(-1.0810477) q[1];
sx q[1];
rz(-1.4996128) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3033894) q[0];
sx q[0];
rz(-1.2801542) q[0];
sx q[0];
rz(-3.035726) q[0];
x q[1];
rz(1.502938) q[2];
sx q[2];
rz(-0.80791622) q[2];
sx q[2];
rz(-2.9615336) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3824532) q[1];
sx q[1];
rz(-2.0331953) q[1];
sx q[1];
rz(0.23666246) q[1];
rz(-pi) q[2];
rz(1.7435944) q[3];
sx q[3];
rz(-1.5732592) q[3];
sx q[3];
rz(0.60009225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2575834) q[2];
sx q[2];
rz(-1.3486226) q[2];
sx q[2];
rz(-2.1949027) q[2];
rz(-2.7729014) q[3];
sx q[3];
rz(-1.576141) q[3];
sx q[3];
rz(2.6856016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7832227) q[0];
sx q[0];
rz(-1.9932278) q[0];
sx q[0];
rz(2.7182462) q[0];
rz(0.070925698) q[1];
sx q[1];
rz(-1.6880886) q[1];
sx q[1];
rz(-0.26500519) q[1];
rz(2.6772066) q[2];
sx q[2];
rz(-2.0225564) q[2];
sx q[2];
rz(0.49973942) q[2];
rz(0.55340135) q[3];
sx q[3];
rz(-2.3407866) q[3];
sx q[3];
rz(1.8368807) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
