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
rz(2.4013588) q[0];
sx q[0];
rz(-1.6594247) q[0];
sx q[0];
rz(-2.8066714) q[0];
rz(-2.6236293) q[1];
sx q[1];
rz(-2.1393175) q[1];
sx q[1];
rz(2.5340773) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018463919) q[0];
sx q[0];
rz(-1.0842609) q[0];
sx q[0];
rz(0.26256613) q[0];
x q[1];
rz(1.363475) q[2];
sx q[2];
rz(-2.6229515) q[2];
sx q[2];
rz(0.27388369) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0820056) q[1];
sx q[1];
rz(-1.4114037) q[1];
sx q[1];
rz(-2.6578147) q[1];
x q[2];
rz(-0.88413864) q[3];
sx q[3];
rz(-0.55301412) q[3];
sx q[3];
rz(2.494604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.5241549) q[2];
sx q[2];
rz(-1.9257156) q[2];
sx q[2];
rz(-2.4784135) q[2];
rz(-0.080862008) q[3];
sx q[3];
rz(-0.20564779) q[3];
sx q[3];
rz(1.1801571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.2422159) q[0];
sx q[0];
rz(-1.3726534) q[0];
sx q[0];
rz(-0.85897613) q[0];
rz(1.8513177) q[1];
sx q[1];
rz(-1.4651508) q[1];
sx q[1];
rz(1.4069517) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039120395) q[0];
sx q[0];
rz(-2.6022096) q[0];
sx q[0];
rz(-2.026019) q[0];
rz(-pi) q[1];
rz(1.2386049) q[2];
sx q[2];
rz(-1.536035) q[2];
sx q[2];
rz(1.3369651) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8415815) q[1];
sx q[1];
rz(-1.4805796) q[1];
sx q[1];
rz(-1.0567699) q[1];
rz(-pi) q[2];
rz(-2.5646016) q[3];
sx q[3];
rz(-1.2170047) q[3];
sx q[3];
rz(-2.4959223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0591639) q[2];
sx q[2];
rz(-2.6291206) q[2];
sx q[2];
rz(1.3272237) q[2];
rz(2.0969157) q[3];
sx q[3];
rz(-2.4278214) q[3];
sx q[3];
rz(-2.7248342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5752207) q[0];
sx q[0];
rz(-0.91082585) q[0];
sx q[0];
rz(0.4253934) q[0];
rz(-1.3771903) q[1];
sx q[1];
rz(-1.4981937) q[1];
sx q[1];
rz(-1.4345217) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4664513) q[0];
sx q[0];
rz(-1.4527307) q[0];
sx q[0];
rz(0.64339126) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6041669) q[2];
sx q[2];
rz(-0.86276189) q[2];
sx q[2];
rz(1.58085) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7075478) q[1];
sx q[1];
rz(-2.2128989) q[1];
sx q[1];
rz(1.0471859) q[1];
rz(-pi) q[2];
rz(0.30607732) q[3];
sx q[3];
rz(-0.97460213) q[3];
sx q[3];
rz(1.5709189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7003358) q[2];
sx q[2];
rz(-1.8852899) q[2];
sx q[2];
rz(-1.8505992) q[2];
rz(1.357632) q[3];
sx q[3];
rz(-1.5057526) q[3];
sx q[3];
rz(1.6443845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.62035471) q[0];
sx q[0];
rz(-2.1339895) q[0];
sx q[0];
rz(-2.3714016) q[0];
rz(2.1417446) q[1];
sx q[1];
rz(-0.60037535) q[1];
sx q[1];
rz(-1.8005449) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4907341) q[0];
sx q[0];
rz(-2.5528209) q[0];
sx q[0];
rz(-0.38932271) q[0];
rz(0.88444986) q[2];
sx q[2];
rz(-1.2053066) q[2];
sx q[2];
rz(-2.3523503) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.51777202) q[1];
sx q[1];
rz(-1.869535) q[1];
sx q[1];
rz(-1.2932528) q[1];
rz(-pi) q[2];
rz(1.8982696) q[3];
sx q[3];
rz(-0.30042111) q[3];
sx q[3];
rz(1.6302366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8009214) q[2];
sx q[2];
rz(-1.069671) q[2];
sx q[2];
rz(0.7684024) q[2];
rz(3.0411804) q[3];
sx q[3];
rz(-1.6008335) q[3];
sx q[3];
rz(1.2787308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95555821) q[0];
sx q[0];
rz(-0.30534196) q[0];
sx q[0];
rz(0.22698639) q[0];
rz(1.3817894) q[1];
sx q[1];
rz(-2.558936) q[1];
sx q[1];
rz(1.7452128) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3811548) q[0];
sx q[0];
rz(-1.9312857) q[0];
sx q[0];
rz(-1.1973778) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7368083) q[2];
sx q[2];
rz(-1.976227) q[2];
sx q[2];
rz(-2.977598) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3815617) q[1];
sx q[1];
rz(-2.2027822) q[1];
sx q[1];
rz(-1.6595592) q[1];
x q[2];
rz(3.1110686) q[3];
sx q[3];
rz(-2.3168193) q[3];
sx q[3];
rz(-1.33547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9153626) q[2];
sx q[2];
rz(-1.9910944) q[2];
sx q[2];
rz(-0.076233141) q[2];
rz(1.7539615) q[3];
sx q[3];
rz(-2.1662655) q[3];
sx q[3];
rz(-1.4001728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48080322) q[0];
sx q[0];
rz(-1.5810409) q[0];
sx q[0];
rz(-1.8060818) q[0];
rz(0.90323365) q[1];
sx q[1];
rz(-1.8025554) q[1];
sx q[1];
rz(-2.9551771) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6064998) q[0];
sx q[0];
rz(-2.4267174) q[0];
sx q[0];
rz(-2.2660648) q[0];
rz(-2.8373986) q[2];
sx q[2];
rz(-0.59202164) q[2];
sx q[2];
rz(-1.6006058) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0707782) q[1];
sx q[1];
rz(-2.1970587) q[1];
sx q[1];
rz(0.29919821) q[1];
x q[2];
rz(2.3535054) q[3];
sx q[3];
rz(-2.4895146) q[3];
sx q[3];
rz(-1.8502082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.67941252) q[2];
sx q[2];
rz(-1.1932411) q[2];
sx q[2];
rz(1.818044) q[2];
rz(-1.1421674) q[3];
sx q[3];
rz(-0.78502941) q[3];
sx q[3];
rz(0.89046684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46397504) q[0];
sx q[0];
rz(-0.1828201) q[0];
sx q[0];
rz(0.63419813) q[0];
rz(2.0967261) q[1];
sx q[1];
rz(-2.2506782) q[1];
sx q[1];
rz(0.96010906) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38518181) q[0];
sx q[0];
rz(-2.7419081) q[0];
sx q[0];
rz(-1.5002285) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5627925) q[2];
sx q[2];
rz(-1.8718613) q[2];
sx q[2];
rz(2.0640399) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.83852488) q[1];
sx q[1];
rz(-0.52553287) q[1];
sx q[1];
rz(1.5938894) q[1];
x q[2];
rz(-2.3784901) q[3];
sx q[3];
rz(-1.0013442) q[3];
sx q[3];
rz(0.72009898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.94379696) q[2];
sx q[2];
rz(-0.65239492) q[2];
sx q[2];
rz(1.5956399) q[2];
rz(-1.724285) q[3];
sx q[3];
rz(-2.3167819) q[3];
sx q[3];
rz(1.5203169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5328131) q[0];
sx q[0];
rz(-2.9488035) q[0];
sx q[0];
rz(2.1667495) q[0];
rz(3.0335562) q[1];
sx q[1];
rz(-1.8861176) q[1];
sx q[1];
rz(-1.9727762) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9038531) q[0];
sx q[0];
rz(-2.0444336) q[0];
sx q[0];
rz(-2.6698378) q[0];
rz(-pi) q[1];
rz(0.76185267) q[2];
sx q[2];
rz(-1.5707865) q[2];
sx q[2];
rz(-2.5229483) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.16202422) q[1];
sx q[1];
rz(-1.8037533) q[1];
sx q[1];
rz(0.090090171) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1049012) q[3];
sx q[3];
rz(-0.49202575) q[3];
sx q[3];
rz(-1.2855315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.10890266) q[2];
sx q[2];
rz(-2.2638075) q[2];
sx q[2];
rz(0.6558134) q[2];
rz(-3.0374895) q[3];
sx q[3];
rz(-1.3452353) q[3];
sx q[3];
rz(-2.9534705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8769237) q[0];
sx q[0];
rz(-1.3035362) q[0];
sx q[0];
rz(1.4917829) q[0];
rz(-2.1022294) q[1];
sx q[1];
rz(-2.3014258) q[1];
sx q[1];
rz(-2.6731491) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.428498) q[0];
sx q[0];
rz(-1.575043) q[0];
sx q[0];
rz(-2.9539786) q[0];
x q[1];
rz(-0.55267398) q[2];
sx q[2];
rz(-2.6330593) q[2];
sx q[2];
rz(0.39297152) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.13889986) q[1];
sx q[1];
rz(-1.4093834) q[1];
sx q[1];
rz(2.3494865) q[1];
rz(-1.8519206) q[3];
sx q[3];
rz(-1.8418962) q[3];
sx q[3];
rz(-0.93096126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.09482065) q[2];
sx q[2];
rz(-1.580661) q[2];
sx q[2];
rz(2.2136733) q[2];
rz(1.3377442) q[3];
sx q[3];
rz(-0.51023054) q[3];
sx q[3];
rz(1.4891589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054166404) q[0];
sx q[0];
rz(-2.1091643) q[0];
sx q[0];
rz(-2.0637276) q[0];
rz(0.36733356) q[1];
sx q[1];
rz(-1.8108188) q[1];
sx q[1];
rz(-1.2841388) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3340942) q[0];
sx q[0];
rz(-2.0413412) q[0];
sx q[0];
rz(2.8583291) q[0];
x q[1];
rz(-0.91725332) q[2];
sx q[2];
rz(-1.4787357) q[2];
sx q[2];
rz(-1.1600509) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.974677) q[1];
sx q[1];
rz(-1.6568746) q[1];
sx q[1];
rz(0.15847107) q[1];
rz(1.359109) q[3];
sx q[3];
rz(-2.2943687) q[3];
sx q[3];
rz(-2.7016957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.0014570634) q[2];
sx q[2];
rz(-0.45150253) q[2];
sx q[2];
rz(2.7122811) q[2];
rz(-0.15520994) q[3];
sx q[3];
rz(-2.8838172) q[3];
sx q[3];
rz(1.3252873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24406381) q[0];
sx q[0];
rz(-1.5915992) q[0];
sx q[0];
rz(-1.5912548) q[0];
rz(2.2776729) q[1];
sx q[1];
rz(-0.37352957) q[1];
sx q[1];
rz(-1.4600798) q[1];
rz(0.39045329) q[2];
sx q[2];
rz(-2.5627372) q[2];
sx q[2];
rz(2.550203) q[2];
rz(1.1170837) q[3];
sx q[3];
rz(-2.4632005) q[3];
sx q[3];
rz(1.7076422) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
