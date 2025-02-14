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
rz(2.2344196) q[0];
sx q[0];
rz(-0.89651674) q[0];
sx q[0];
rz(-0.040891115) q[0];
rz(-1.2603124) q[1];
sx q[1];
rz(-2.6949096) q[1];
sx q[1];
rz(0.094951542) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1563383) q[0];
sx q[0];
rz(-1.3438936) q[0];
sx q[0];
rz(2.7309523) q[0];
x q[1];
rz(-0.55741252) q[2];
sx q[2];
rz(-1.3746309) q[2];
sx q[2];
rz(0.74166751) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.6247926) q[1];
sx q[1];
rz(-1.2007371) q[1];
sx q[1];
rz(-2.3462252) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0980174) q[3];
sx q[3];
rz(-1.1679139) q[3];
sx q[3];
rz(-1.948818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8591259) q[2];
sx q[2];
rz(-0.60255113) q[2];
sx q[2];
rz(-1.9762543) q[2];
rz(0.91853842) q[3];
sx q[3];
rz(-0.10401741) q[3];
sx q[3];
rz(-1.9281841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4024046) q[0];
sx q[0];
rz(-2.5027051) q[0];
sx q[0];
rz(-2.8366587) q[0];
rz(2.1091499) q[1];
sx q[1];
rz(-0.28034261) q[1];
sx q[1];
rz(0.84411821) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4979008) q[0];
sx q[0];
rz(-1.7232753) q[0];
sx q[0];
rz(1.8122079) q[0];
x q[1];
rz(1.3621037) q[2];
sx q[2];
rz(-1.450453) q[2];
sx q[2];
rz(2.142765) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1010437) q[1];
sx q[1];
rz(-1.1660468) q[1];
sx q[1];
rz(-1.8331265) q[1];
x q[2];
rz(0.02133298) q[3];
sx q[3];
rz(-0.6271242) q[3];
sx q[3];
rz(2.5950678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1673163) q[2];
sx q[2];
rz(-0.68334371) q[2];
sx q[2];
rz(0.58504504) q[2];
rz(0.85917464) q[3];
sx q[3];
rz(-1.1480568) q[3];
sx q[3];
rz(-1.1332716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0228731) q[0];
sx q[0];
rz(-2.3206503) q[0];
sx q[0];
rz(-0.098544772) q[0];
rz(0.56602829) q[1];
sx q[1];
rz(-0.84605828) q[1];
sx q[1];
rz(-0.50618323) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1126668) q[0];
sx q[0];
rz(-0.62917626) q[0];
sx q[0];
rz(-2.5553246) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82019851) q[2];
sx q[2];
rz(-1.9026105) q[2];
sx q[2];
rz(-2.2143176) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4892046) q[1];
sx q[1];
rz(-1.694177) q[1];
sx q[1];
rz(1.3493933) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1428301) q[3];
sx q[3];
rz(-1.9062098) q[3];
sx q[3];
rz(2.4728925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.026406) q[2];
sx q[2];
rz(-1.2135442) q[2];
sx q[2];
rz(3.0565267) q[2];
rz(2.3885942) q[3];
sx q[3];
rz(-2.4716061) q[3];
sx q[3];
rz(-1.6116713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6104777) q[0];
sx q[0];
rz(-1.5014638) q[0];
sx q[0];
rz(0.53989545) q[0];
rz(0.029622948) q[1];
sx q[1];
rz(-0.74631515) q[1];
sx q[1];
rz(-1.354904) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52068096) q[0];
sx q[0];
rz(-1.3518466) q[0];
sx q[0];
rz(-2.3205873) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7672655) q[2];
sx q[2];
rz(-2.3951963) q[2];
sx q[2];
rz(2.057892) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.190092) q[1];
sx q[1];
rz(-1.8279549) q[1];
sx q[1];
rz(-0.58488226) q[1];
rz(-0.79014312) q[3];
sx q[3];
rz(-1.5075659) q[3];
sx q[3];
rz(2.9823398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.70298355) q[2];
sx q[2];
rz(-2.1712124) q[2];
sx q[2];
rz(0.89540974) q[2];
rz(-1.4488975) q[3];
sx q[3];
rz(-1.9796895) q[3];
sx q[3];
rz(-0.40941516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92622906) q[0];
sx q[0];
rz(-2.138593) q[0];
sx q[0];
rz(-0.0005501752) q[0];
rz(-2.6161361) q[1];
sx q[1];
rz(-2.2976687) q[1];
sx q[1];
rz(2.1702683) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0063365) q[0];
sx q[0];
rz(-0.86733666) q[0];
sx q[0];
rz(1.9869861) q[0];
rz(-1.1575562) q[2];
sx q[2];
rz(-1.2451425) q[2];
sx q[2];
rz(-0.60631863) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.57572039) q[1];
sx q[1];
rz(-2.1196869) q[1];
sx q[1];
rz(0.68419211) q[1];
rz(-1.8841142) q[3];
sx q[3];
rz(-1.3660114) q[3];
sx q[3];
rz(1.1623032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9690669) q[2];
sx q[2];
rz(-1.0774287) q[2];
sx q[2];
rz(1.7871008) q[2];
rz(0.6012249) q[3];
sx q[3];
rz(-2.0982592) q[3];
sx q[3];
rz(-1.5662947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4300267) q[0];
sx q[0];
rz(-2.8373748) q[0];
sx q[0];
rz(2.8416204) q[0];
rz(-2.1638339) q[1];
sx q[1];
rz(-2.561196) q[1];
sx q[1];
rz(2.207644) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47563206) q[0];
sx q[0];
rz(-2.7802659) q[0];
sx q[0];
rz(-2.7368109) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0522473) q[2];
sx q[2];
rz(-2.1314959) q[2];
sx q[2];
rz(-2.581832) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0910513) q[1];
sx q[1];
rz(-0.81650298) q[1];
sx q[1];
rz(1.8920808) q[1];
x q[2];
rz(1.2190231) q[3];
sx q[3];
rz(-1.1922424) q[3];
sx q[3];
rz(1.5768676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0710435) q[2];
sx q[2];
rz(-2.1328378) q[2];
sx q[2];
rz(-2.0096931) q[2];
rz(-1.8848568) q[3];
sx q[3];
rz(-2.3102424) q[3];
sx q[3];
rz(1.7251714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8202332) q[0];
sx q[0];
rz(-1.2437404) q[0];
sx q[0];
rz(2.5309122) q[0];
rz(0.75675476) q[1];
sx q[1];
rz(-0.38833955) q[1];
sx q[1];
rz(1.6441708) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9349139) q[0];
sx q[0];
rz(-1.4467708) q[0];
sx q[0];
rz(-1.5453611) q[0];
rz(-pi) q[1];
rz(2.2708504) q[2];
sx q[2];
rz(-0.53360924) q[2];
sx q[2];
rz(2.6117532) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1042133) q[1];
sx q[1];
rz(-0.92475212) q[1];
sx q[1];
rz(-1.7039596) q[1];
x q[2];
rz(3.0277062) q[3];
sx q[3];
rz(-2.4166757) q[3];
sx q[3];
rz(2.8020086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9553767) q[2];
sx q[2];
rz(-2.4775041) q[2];
sx q[2];
rz(-0.60626283) q[2];
rz(-2.9210505) q[3];
sx q[3];
rz(-1.3371779) q[3];
sx q[3];
rz(-0.36673275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99899387) q[0];
sx q[0];
rz(-1.6479011) q[0];
sx q[0];
rz(0.9077453) q[0];
rz(1.5331462) q[1];
sx q[1];
rz(-0.95450675) q[1];
sx q[1];
rz(-0.018891637) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2174199) q[0];
sx q[0];
rz(-2.0402363) q[0];
sx q[0];
rz(1.4408235) q[0];
rz(-pi) q[1];
rz(-0.18757815) q[2];
sx q[2];
rz(-0.30454985) q[2];
sx q[2];
rz(-1.2414602) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0844473) q[1];
sx q[1];
rz(-1.6020163) q[1];
sx q[1];
rz(-2.679326) q[1];
rz(-pi) q[2];
rz(0.83062474) q[3];
sx q[3];
rz(-1.1234094) q[3];
sx q[3];
rz(-0.4448286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3733526) q[2];
sx q[2];
rz(-1.588593) q[2];
sx q[2];
rz(0.24660435) q[2];
rz(1.2982347) q[3];
sx q[3];
rz(-1.4539098) q[3];
sx q[3];
rz(1.8728144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75803718) q[0];
sx q[0];
rz(-1.0724496) q[0];
sx q[0];
rz(-0.86395907) q[0];
rz(-1.9268688) q[1];
sx q[1];
rz(-1.1178958) q[1];
sx q[1];
rz(2.5606959) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5340849) q[0];
sx q[0];
rz(-2.0568741) q[0];
sx q[0];
rz(-0.75774225) q[0];
rz(-pi) q[1];
rz(-0.038453416) q[2];
sx q[2];
rz(-2.3068301) q[2];
sx q[2];
rz(-2.9308386) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.89199663) q[1];
sx q[1];
rz(-2.432715) q[1];
sx q[1];
rz(-0.53541553) q[1];
rz(-pi) q[2];
rz(0.77928752) q[3];
sx q[3];
rz(-1.7699935) q[3];
sx q[3];
rz(1.3785386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8860127) q[2];
sx q[2];
rz(-2.6801127) q[2];
sx q[2];
rz(-2.9542921) q[2];
rz(2.1646132) q[3];
sx q[3];
rz(-1.6410442) q[3];
sx q[3];
rz(1.8627082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63000694) q[0];
sx q[0];
rz(-0.66660175) q[0];
sx q[0];
rz(-1.8774207) q[0];
rz(-0.97201792) q[1];
sx q[1];
rz(-0.7518026) q[1];
sx q[1];
rz(1.1690296) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.326871) q[0];
sx q[0];
rz(-1.4356614) q[0];
sx q[0];
rz(-1.9072471) q[0];
rz(-2.3832267) q[2];
sx q[2];
rz(-0.50687271) q[2];
sx q[2];
rz(3.1257983) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.059140511) q[1];
sx q[1];
rz(-0.881221) q[1];
sx q[1];
rz(0.13549681) q[1];
rz(-0.33785744) q[3];
sx q[3];
rz(-0.33379972) q[3];
sx q[3];
rz(-0.56188449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6675889) q[2];
sx q[2];
rz(-1.7550125) q[2];
sx q[2];
rz(-0.37707314) q[2];
rz(0.22370473) q[3];
sx q[3];
rz(-0.5880028) q[3];
sx q[3];
rz(1.912775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2558462) q[0];
sx q[0];
rz(-1.3930014) q[0];
sx q[0];
rz(-0.2572671) q[0];
rz(-2.5480351) q[1];
sx q[1];
rz(-2.3496353) q[1];
sx q[1];
rz(-1.4812462) q[1];
rz(2.3667546) q[2];
sx q[2];
rz(-0.88817468) q[2];
sx q[2];
rz(-1.9356288) q[2];
rz(1.0026686) q[3];
sx q[3];
rz(-0.46350652) q[3];
sx q[3];
rz(1.9476611) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
