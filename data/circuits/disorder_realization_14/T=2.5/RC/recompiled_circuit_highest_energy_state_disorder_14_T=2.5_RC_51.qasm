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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54017055) q[0];
sx q[0];
rz(-2.5937732) q[0];
sx q[0];
rz(-1.1146077) q[0];
rz(-1.363475) q[2];
sx q[2];
rz(-2.6229515) q[2];
sx q[2];
rz(2.867709) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0820056) q[1];
sx q[1];
rz(-1.4114037) q[1];
sx q[1];
rz(-2.6578147) q[1];
rz(-pi) q[2];
rz(-2.7685952) q[3];
sx q[3];
rz(-1.1524876) q[3];
sx q[3];
rz(-1.7278863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.5241549) q[2];
sx q[2];
rz(-1.9257156) q[2];
sx q[2];
rz(2.4784135) q[2];
rz(-3.0607306) q[3];
sx q[3];
rz(-2.9359449) q[3];
sx q[3];
rz(1.1801571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2422159) q[0];
sx q[0];
rz(-1.7689393) q[0];
sx q[0];
rz(-0.85897613) q[0];
rz(1.8513177) q[1];
sx q[1];
rz(-1.6764418) q[1];
sx q[1];
rz(-1.4069517) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0075571) q[0];
sx q[0];
rz(-1.7985744) q[0];
sx q[0];
rz(2.0640949) q[0];
x q[1];
rz(-1.2386049) q[2];
sx q[2];
rz(-1.6055577) q[2];
sx q[2];
rz(-1.8046276) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3000112) q[1];
sx q[1];
rz(-1.661013) q[1];
sx q[1];
rz(-2.0848227) q[1];
rz(-pi) q[2];
rz(0.57699109) q[3];
sx q[3];
rz(-1.9245879) q[3];
sx q[3];
rz(2.4959223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0591639) q[2];
sx q[2];
rz(-0.51247207) q[2];
sx q[2];
rz(1.3272237) q[2];
rz(-2.0969157) q[3];
sx q[3];
rz(-0.71377126) q[3];
sx q[3];
rz(0.41675848) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5663719) q[0];
sx q[0];
rz(-0.91082585) q[0];
sx q[0];
rz(-2.7161993) q[0];
rz(1.3771903) q[1];
sx q[1];
rz(-1.6433989) q[1];
sx q[1];
rz(1.707071) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1578428) q[0];
sx q[0];
rz(-2.2089777) q[0];
sx q[0];
rz(1.7179836) q[0];
rz(-pi) q[1];
rz(-0.70830958) q[2];
sx q[2];
rz(-1.596144) q[2];
sx q[2];
rz(0.011653221) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4340448) q[1];
sx q[1];
rz(-2.2128989) q[1];
sx q[1];
rz(-2.0944067) q[1];
x q[2];
rz(2.1892912) q[3];
sx q[3];
rz(-1.3188014) q[3];
sx q[3];
rz(-2.9661055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.44125685) q[2];
sx q[2];
rz(-1.8852899) q[2];
sx q[2];
rz(1.8505992) q[2];
rz(1.7839606) q[3];
sx q[3];
rz(-1.5057526) q[3];
sx q[3];
rz(-1.6443845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5212379) q[0];
sx q[0];
rz(-1.0076032) q[0];
sx q[0];
rz(2.3714016) q[0];
rz(-2.1417446) q[1];
sx q[1];
rz(-0.60037535) q[1];
sx q[1];
rz(-1.3410478) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65085852) q[0];
sx q[0];
rz(-0.58877173) q[0];
sx q[0];
rz(0.38932271) q[0];
x q[1];
rz(2.1140587) q[2];
sx q[2];
rz(-0.76342602) q[2];
sx q[2];
rz(-1.9486519) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8902258) q[1];
sx q[1];
rz(-0.40491762) q[1];
sx q[1];
rz(-0.7271073) q[1];
rz(-1.2854659) q[3];
sx q[3];
rz(-1.6661246) q[3];
sx q[3];
rz(0.25432977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8009214) q[2];
sx q[2];
rz(-1.069671) q[2];
sx q[2];
rz(0.7684024) q[2];
rz(0.10041222) q[3];
sx q[3];
rz(-1.6008335) q[3];
sx q[3];
rz(1.8628619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1860344) q[0];
sx q[0];
rz(-0.30534196) q[0];
sx q[0];
rz(0.22698639) q[0];
rz(-1.7598033) q[1];
sx q[1];
rz(-2.558936) q[1];
sx q[1];
rz(-1.3963799) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3811548) q[0];
sx q[0];
rz(-1.9312857) q[0];
sx q[0];
rz(-1.1973778) q[0];
x q[1];
rz(-0.41047217) q[2];
sx q[2];
rz(-1.7232401) q[2];
sx q[2];
rz(-1.6688011) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5310881) q[1];
sx q[1];
rz(-2.5042494) q[1];
sx q[1];
rz(3.0211041) q[1];
rz(0.030524039) q[3];
sx q[3];
rz(-2.3168193) q[3];
sx q[3];
rz(-1.8061226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9153626) q[2];
sx q[2];
rz(-1.9910944) q[2];
sx q[2];
rz(0.076233141) q[2];
rz(-1.7539615) q[3];
sx q[3];
rz(-2.1662655) q[3];
sx q[3];
rz(-1.7414198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.48080322) q[0];
sx q[0];
rz(-1.5605518) q[0];
sx q[0];
rz(1.8060818) q[0];
rz(2.238359) q[1];
sx q[1];
rz(-1.8025554) q[1];
sx q[1];
rz(-0.18641557) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53509287) q[0];
sx q[0];
rz(-0.71487521) q[0];
sx q[0];
rz(2.2660648) q[0];
rz(0.57045649) q[2];
sx q[2];
rz(-1.738731) q[2];
sx q[2];
rz(-2.8569375) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0708145) q[1];
sx q[1];
rz(-0.94453393) q[1];
sx q[1];
rz(-0.29919821) q[1];
x q[2];
rz(1.0746434) q[3];
sx q[3];
rz(-1.1285787) q[3];
sx q[3];
rz(0.38954716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.67941252) q[2];
sx q[2];
rz(-1.9483515) q[2];
sx q[2];
rz(-1.3235486) q[2];
rz(1.9994252) q[3];
sx q[3];
rz(-0.78502941) q[3];
sx q[3];
rz(0.89046684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6776176) q[0];
sx q[0];
rz(-0.1828201) q[0];
sx q[0];
rz(2.5073945) q[0];
rz(2.0967261) q[1];
sx q[1];
rz(-2.2506782) q[1];
sx q[1];
rz(-2.1814836) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38518181) q[0];
sx q[0];
rz(-0.39968458) q[0];
sx q[0];
rz(-1.5002285) q[0];
rz(0.025770806) q[2];
sx q[2];
rz(-2.8404245) q[2];
sx q[2];
rz(-1.1045375) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.81183103) q[1];
sx q[1];
rz(-1.0454181) q[1];
sx q[1];
rz(-3.1282022) q[1];
rz(-pi) q[2];
rz(2.3784901) q[3];
sx q[3];
rz(-1.0013442) q[3];
sx q[3];
rz(-0.72009898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1977957) q[2];
sx q[2];
rz(-0.65239492) q[2];
sx q[2];
rz(-1.5459527) q[2];
rz(1.724285) q[3];
sx q[3];
rz(-0.82481074) q[3];
sx q[3];
rz(1.5203169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6087795) q[0];
sx q[0];
rz(-0.19278917) q[0];
sx q[0];
rz(2.1667495) q[0];
rz(0.10803647) q[1];
sx q[1];
rz(-1.255475) q[1];
sx q[1];
rz(-1.9727762) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2377396) q[0];
sx q[0];
rz(-2.0444336) q[0];
sx q[0];
rz(0.4717548) q[0];
x q[1];
rz(-1.5708099) q[2];
sx q[2];
rz(-2.332649) q[2];
sx q[2];
rz(-2.1894313) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9795684) q[1];
sx q[1];
rz(-1.8037533) q[1];
sx q[1];
rz(0.090090171) q[1];
rz(-pi) q[2];
rz(-1.5511369) q[3];
sx q[3];
rz(-2.0624614) q[3];
sx q[3];
rz(1.2439072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.10890266) q[2];
sx q[2];
rz(-2.2638075) q[2];
sx q[2];
rz(-0.6558134) q[2];
rz(3.0374895) q[3];
sx q[3];
rz(-1.7963573) q[3];
sx q[3];
rz(-2.9534705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26466894) q[0];
sx q[0];
rz(-1.3035362) q[0];
sx q[0];
rz(1.4917829) q[0];
rz(-1.0393633) q[1];
sx q[1];
rz(-2.3014258) q[1];
sx q[1];
rz(2.6731491) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2830848) q[0];
sx q[0];
rz(-1.383184) q[0];
sx q[0];
rz(1.5664738) q[0];
rz(-pi) q[1];
rz(2.5889187) q[2];
sx q[2];
rz(-2.6330593) q[2];
sx q[2];
rz(0.39297152) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0026928) q[1];
sx q[1];
rz(-1.7322092) q[1];
sx q[1];
rz(0.79210613) q[1];
x q[2];
rz(1.2896721) q[3];
sx q[3];
rz(-1.8418962) q[3];
sx q[3];
rz(-0.93096126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.046772) q[2];
sx q[2];
rz(-1.580661) q[2];
sx q[2];
rz(-2.2136733) q[2];
rz(1.3377442) q[3];
sx q[3];
rz(-0.51023054) q[3];
sx q[3];
rz(-1.6524338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054166404) q[0];
sx q[0];
rz(-1.0324284) q[0];
sx q[0];
rz(1.077865) q[0];
rz(0.36733356) q[1];
sx q[1];
rz(-1.8108188) q[1];
sx q[1];
rz(1.8574538) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76332247) q[0];
sx q[0];
rz(-0.54369944) q[0];
sx q[0];
rz(-2.0732353) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7214891) q[2];
sx q[2];
rz(-0.65905276) q[2];
sx q[2];
rz(0.53021741) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.91032797) q[1];
sx q[1];
rz(-2.961425) q[1];
sx q[1];
rz(-2.641201) q[1];
x q[2];
rz(-2.4068042) q[3];
sx q[3];
rz(-1.7289203) q[3];
sx q[3];
rz(1.8693592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.0014570634) q[2];
sx q[2];
rz(-2.6900901) q[2];
sx q[2];
rz(2.7122811) q[2];
rz(2.9863827) q[3];
sx q[3];
rz(-0.25777543) q[3];
sx q[3];
rz(-1.3252873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24406381) q[0];
sx q[0];
rz(-1.5915992) q[0];
sx q[0];
rz(-1.5912548) q[0];
rz(-0.86391972) q[1];
sx q[1];
rz(-0.37352957) q[1];
sx q[1];
rz(-1.4600798) q[1];
rz(2.5979832) q[2];
sx q[2];
rz(-1.3610441) q[2];
sx q[2];
rz(1.3112031) q[2];
rz(-2.024509) q[3];
sx q[3];
rz(-2.4632005) q[3];
sx q[3];
rz(1.7076422) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
