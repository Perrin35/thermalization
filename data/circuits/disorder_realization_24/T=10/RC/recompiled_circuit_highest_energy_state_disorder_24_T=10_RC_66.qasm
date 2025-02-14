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
rz(3.0348294) q[0];
sx q[0];
rz(-1.3314629) q[0];
sx q[0];
rz(-1.7911628) q[0];
rz(0.30836937) q[1];
sx q[1];
rz(6.0344459) q[1];
sx q[1];
rz(11.587172) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5677264) q[0];
sx q[0];
rz(-1.4970333) q[0];
sx q[0];
rz(0.7598138) q[0];
x q[1];
rz(2.0151695) q[2];
sx q[2];
rz(-0.70894402) q[2];
sx q[2];
rz(-0.63731784) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.82047082) q[1];
sx q[1];
rz(-0.36013569) q[1];
sx q[1];
rz(0.19345691) q[1];
rz(-pi) q[2];
rz(-0.37110801) q[3];
sx q[3];
rz(-1.4083997) q[3];
sx q[3];
rz(1.1134256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2778492) q[2];
sx q[2];
rz(-1.2110445) q[2];
sx q[2];
rz(2.5085874) q[2];
rz(2.4999319) q[3];
sx q[3];
rz(-0.46052614) q[3];
sx q[3];
rz(-0.7707001) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68384701) q[0];
sx q[0];
rz(-2.1774543) q[0];
sx q[0];
rz(2.316851) q[0];
rz(-1.9174346) q[1];
sx q[1];
rz(-2.2861202) q[1];
sx q[1];
rz(2.8214084) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8267248) q[0];
sx q[0];
rz(-1.1629675) q[0];
sx q[0];
rz(0.5023707) q[0];
rz(-pi) q[1];
rz(1.2037781) q[2];
sx q[2];
rz(-0.53897714) q[2];
sx q[2];
rz(-2.0808329) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9588543) q[1];
sx q[1];
rz(-0.08402782) q[1];
sx q[1];
rz(-0.27167222) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8289964) q[3];
sx q[3];
rz(-0.22919433) q[3];
sx q[3];
rz(-0.75507873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0188633) q[2];
sx q[2];
rz(-1.9247232) q[2];
sx q[2];
rz(-1.8625721) q[2];
rz(3.0458798) q[3];
sx q[3];
rz(-2.0666104) q[3];
sx q[3];
rz(-1.3191282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9886446) q[0];
sx q[0];
rz(-0.64866346) q[0];
sx q[0];
rz(2.1107819) q[0];
rz(-2.5910494) q[1];
sx q[1];
rz(-0.85854733) q[1];
sx q[1];
rz(-1.8969089) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6129788) q[0];
sx q[0];
rz(-1.6242149) q[0];
sx q[0];
rz(-2.9821755) q[0];
rz(-pi) q[1];
rz(-0.95374505) q[2];
sx q[2];
rz(-1.707952) q[2];
sx q[2];
rz(1.1378261) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66164368) q[1];
sx q[1];
rz(-0.91001653) q[1];
sx q[1];
rz(0.26008545) q[1];
rz(-1.1641509) q[3];
sx q[3];
rz(-0.79504025) q[3];
sx q[3];
rz(-0.94889489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9024258) q[2];
sx q[2];
rz(-0.5639762) q[2];
sx q[2];
rz(-1.9695367) q[2];
rz(0.82550448) q[3];
sx q[3];
rz(-1.8768616) q[3];
sx q[3];
rz(-0.66506213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.091752) q[0];
sx q[0];
rz(-1.6850152) q[0];
sx q[0];
rz(-0.41393429) q[0];
rz(-0.44231689) q[1];
sx q[1];
rz(-2.1374173) q[1];
sx q[1];
rz(2.5962459) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53876153) q[0];
sx q[0];
rz(-1.6206546) q[0];
sx q[0];
rz(0.2754579) q[0];
x q[1];
rz(-2.4905861) q[2];
sx q[2];
rz(-0.88906139) q[2];
sx q[2];
rz(1.4780457) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2703823) q[1];
sx q[1];
rz(-2.5426513) q[1];
sx q[1];
rz(-2.0763055) q[1];
x q[2];
rz(-2.0796989) q[3];
sx q[3];
rz(-1.4413222) q[3];
sx q[3];
rz(-2.9751301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.44117323) q[2];
sx q[2];
rz(-1.8808695) q[2];
sx q[2];
rz(0.1612266) q[2];
rz(1.5876596) q[3];
sx q[3];
rz(-0.59993887) q[3];
sx q[3];
rz(-0.59066311) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.081838354) q[0];
sx q[0];
rz(-1.8738382) q[0];
sx q[0];
rz(-2.4181714) q[0];
rz(1.3149892) q[1];
sx q[1];
rz(-1.864121) q[1];
sx q[1];
rz(-1.0378729) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53053601) q[0];
sx q[0];
rz(-0.7270455) q[0];
sx q[0];
rz(0.1188936) q[0];
x q[1];
rz(-2.5332301) q[2];
sx q[2];
rz(-2.3181097) q[2];
sx q[2];
rz(1.057511) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.41748504) q[1];
sx q[1];
rz(-1.8383664) q[1];
sx q[1];
rz(-2.6883283) q[1];
rz(-0.54609046) q[3];
sx q[3];
rz(-1.6973064) q[3];
sx q[3];
rz(-1.5777335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3967241) q[2];
sx q[2];
rz(-2.0266271) q[2];
sx q[2];
rz(-2.569516) q[2];
rz(-2.5165596) q[3];
sx q[3];
rz(-1.0038989) q[3];
sx q[3];
rz(-1.9151789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52142757) q[0];
sx q[0];
rz(-3.1075952) q[0];
sx q[0];
rz(-0.52325621) q[0];
rz(-1.322586) q[1];
sx q[1];
rz(-1.6736284) q[1];
sx q[1];
rz(-0.62017131) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70411982) q[0];
sx q[0];
rz(-1.0123555) q[0];
sx q[0];
rz(2.0790786) q[0];
rz(1.6441395) q[2];
sx q[2];
rz(-2.0383325) q[2];
sx q[2];
rz(0.95638025) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.090724548) q[1];
sx q[1];
rz(-1.6584937) q[1];
sx q[1];
rz(2.1983653) q[1];
rz(-1.3955388) q[3];
sx q[3];
rz(-1.7850998) q[3];
sx q[3];
rz(2.8068723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.72914499) q[2];
sx q[2];
rz(-0.77249384) q[2];
sx q[2];
rz(-1.283851) q[2];
rz(0.9681975) q[3];
sx q[3];
rz(-1.7627534) q[3];
sx q[3];
rz(-1.5865954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9667483) q[0];
sx q[0];
rz(-2.0059678) q[0];
sx q[0];
rz(-1.211776) q[0];
rz(-2.1719596) q[1];
sx q[1];
rz(-1.8200834) q[1];
sx q[1];
rz(-1.8671573) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2454555) q[0];
sx q[0];
rz(-1.4198167) q[0];
sx q[0];
rz(2.2666988) q[0];
rz(-pi) q[1];
rz(-1.4757122) q[2];
sx q[2];
rz(-2.3543752) q[2];
sx q[2];
rz(-1.1457535) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.695622) q[1];
sx q[1];
rz(-1.2685163) q[1];
sx q[1];
rz(-0.92960371) q[1];
rz(-1.8936552) q[3];
sx q[3];
rz(-2.0104694) q[3];
sx q[3];
rz(-2.6722801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6756639) q[2];
sx q[2];
rz(-1.482654) q[2];
sx q[2];
rz(0.3234123) q[2];
rz(-0.0023500738) q[3];
sx q[3];
rz(-1.6833143) q[3];
sx q[3];
rz(-2.5200444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49331409) q[0];
sx q[0];
rz(-1.4819773) q[0];
sx q[0];
rz(-1.3508654) q[0];
rz(1.3510652) q[1];
sx q[1];
rz(-1.8565145) q[1];
sx q[1];
rz(-1.8538808) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0553594) q[0];
sx q[0];
rz(-1.4533193) q[0];
sx q[0];
rz(-1.8553977) q[0];
rz(-pi) q[1];
rz(-1.9016198) q[2];
sx q[2];
rz(-0.74218732) q[2];
sx q[2];
rz(-1.0617219) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1786036) q[1];
sx q[1];
rz(-1.4006536) q[1];
sx q[1];
rz(-0.32802204) q[1];
x q[2];
rz(-3.0898422) q[3];
sx q[3];
rz(-0.34933511) q[3];
sx q[3];
rz(-0.83788315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6172341) q[2];
sx q[2];
rz(-1.989216) q[2];
sx q[2];
rz(-2.1481245) q[2];
rz(2.1067545) q[3];
sx q[3];
rz(-2.5351758) q[3];
sx q[3];
rz(-2.9724227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33191037) q[0];
sx q[0];
rz(-1.9669635) q[0];
sx q[0];
rz(-1.6140953) q[0];
rz(2.2612259) q[1];
sx q[1];
rz(-2.7208734) q[1];
sx q[1];
rz(2.8299433) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6821996) q[0];
sx q[0];
rz(-1.2785204) q[0];
sx q[0];
rz(-2.9351013) q[0];
rz(-pi) q[1];
rz(-1.6402354) q[2];
sx q[2];
rz(-0.23190325) q[2];
sx q[2];
rz(-2.7551485) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5376725) q[1];
sx q[1];
rz(-2.2419856) q[1];
sx q[1];
rz(-1.0965986) q[1];
rz(2.8086189) q[3];
sx q[3];
rz(-2.1978488) q[3];
sx q[3];
rz(-2.0505333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46939048) q[2];
sx q[2];
rz(-0.92350525) q[2];
sx q[2];
rz(-1.7507318) q[2];
rz(-1.6145128) q[3];
sx q[3];
rz(-2.0938087) q[3];
sx q[3];
rz(-1.0745777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63956082) q[0];
sx q[0];
rz(-1.1447516) q[0];
sx q[0];
rz(0.61258739) q[0];
rz(-2.1514429) q[1];
sx q[1];
rz(-1.1724816) q[1];
sx q[1];
rz(1.6506857) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6263437) q[0];
sx q[0];
rz(-0.83413163) q[0];
sx q[0];
rz(2.5596928) q[0];
rz(2.886495) q[2];
sx q[2];
rz(-2.2146985) q[2];
sx q[2];
rz(1.7718499) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.13496357) q[1];
sx q[1];
rz(-1.4093168) q[1];
sx q[1];
rz(3.1352741) q[1];
rz(-pi) q[2];
rz(-1.321071) q[3];
sx q[3];
rz(-2.2533544) q[3];
sx q[3];
rz(-2.5988795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.426173) q[2];
sx q[2];
rz(-0.82948589) q[2];
sx q[2];
rz(3.1066011) q[2];
rz(1.4404826) q[3];
sx q[3];
rz(-0.8927497) q[3];
sx q[3];
rz(0.54097241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4424425) q[0];
sx q[0];
rz(-2.2461666) q[0];
sx q[0];
rz(2.9641892) q[0];
rz(1.9433446) q[1];
sx q[1];
rz(-1.6234963) q[1];
sx q[1];
rz(-2.1877098) q[1];
rz(2.4468471) q[2];
sx q[2];
rz(-0.80949819) q[2];
sx q[2];
rz(-0.037485952) q[2];
rz(-1.1097601) q[3];
sx q[3];
rz(-0.90519917) q[3];
sx q[3];
rz(-1.7310033) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
