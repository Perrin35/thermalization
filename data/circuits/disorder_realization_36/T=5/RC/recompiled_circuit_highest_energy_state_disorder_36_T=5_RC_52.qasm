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
rz(-0.90717301) q[0];
sx q[0];
rz(-2.2450759) q[0];
sx q[0];
rz(-3.1007015) q[0];
rz(-1.2603124) q[1];
sx q[1];
rz(-2.6949096) q[1];
sx q[1];
rz(0.094951542) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4584158) q[0];
sx q[0];
rz(-1.1712892) q[0];
sx q[0];
rz(1.3241121) q[0];
x q[1];
rz(0.35934885) q[2];
sx q[2];
rz(-0.58746019) q[2];
sx q[2];
rz(1.1321783) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5168001) q[1];
sx q[1];
rz(-1.9408556) q[1];
sx q[1];
rz(0.79536749) q[1];
rz(-pi) q[2];
rz(-2.6951249) q[3];
sx q[3];
rz(-2.0030336) q[3];
sx q[3];
rz(0.18011506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8591259) q[2];
sx q[2];
rz(-2.5390415) q[2];
sx q[2];
rz(1.1653384) q[2];
rz(2.2230542) q[3];
sx q[3];
rz(-0.10401741) q[3];
sx q[3];
rz(-1.2134086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4024046) q[0];
sx q[0];
rz(-2.5027051) q[0];
sx q[0];
rz(-0.30493394) q[0];
rz(-2.1091499) q[1];
sx q[1];
rz(-0.28034261) q[1];
sx q[1];
rz(-0.84411821) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88972487) q[0];
sx q[0];
rz(-1.8093523) q[0];
sx q[0];
rz(-0.15695928) q[0];
rz(-pi) q[1];
rz(1.0424625) q[2];
sx q[2];
rz(-2.9011263) q[2];
sx q[2];
rz(-0.056253091) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6355953) q[1];
sx q[1];
rz(-1.330101) q[1];
sx q[1];
rz(-0.41750112) q[1];
x q[2];
rz(0.62701608) q[3];
sx q[3];
rz(-1.5833143) q[3];
sx q[3];
rz(-1.0069969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1673163) q[2];
sx q[2];
rz(-0.68334371) q[2];
sx q[2];
rz(0.58504504) q[2];
rz(-0.85917464) q[3];
sx q[3];
rz(-1.9935358) q[3];
sx q[3];
rz(2.008321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1187196) q[0];
sx q[0];
rz(-0.82094231) q[0];
sx q[0];
rz(-0.098544772) q[0];
rz(0.56602829) q[1];
sx q[1];
rz(-2.2955344) q[1];
sx q[1];
rz(-2.6354094) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1126668) q[0];
sx q[0];
rz(-0.62917626) q[0];
sx q[0];
rz(-0.58626808) q[0];
rz(2.0385582) q[2];
sx q[2];
rz(-2.3342103) q[2];
sx q[2];
rz(0.30750662) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0876956) q[1];
sx q[1];
rz(-1.7904886) q[1];
sx q[1];
rz(3.0151574) q[1];
rz(-pi) q[2];
rz(2.1428895) q[3];
sx q[3];
rz(-0.65352189) q[3];
sx q[3];
rz(-1.7668949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.026406) q[2];
sx q[2];
rz(-1.2135442) q[2];
sx q[2];
rz(-0.085065993) q[2];
rz(-0.75299844) q[3];
sx q[3];
rz(-2.4716061) q[3];
sx q[3];
rz(-1.6116713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6104777) q[0];
sx q[0];
rz(-1.6401289) q[0];
sx q[0];
rz(2.6016972) q[0];
rz(3.1119697) q[1];
sx q[1];
rz(-2.3952775) q[1];
sx q[1];
rz(1.7866887) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8919075) q[0];
sx q[0];
rz(-2.2985795) q[0];
sx q[0];
rz(-2.84642) q[0];
x q[1];
rz(0.71075534) q[2];
sx q[2];
rz(-1.3198993) q[2];
sx q[2];
rz(-2.9352842) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.78614178) q[1];
sx q[1];
rz(-1.0075354) q[1];
sx q[1];
rz(-1.8763297) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0527086) q[3];
sx q[3];
rz(-2.349472) q[3];
sx q[3];
rz(1.3490335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.70298355) q[2];
sx q[2];
rz(-2.1712124) q[2];
sx q[2];
rz(-0.89540974) q[2];
rz(1.6926951) q[3];
sx q[3];
rz(-1.1619032) q[3];
sx q[3];
rz(-2.7321775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92622906) q[0];
sx q[0];
rz(-1.0029997) q[0];
sx q[0];
rz(-0.0005501752) q[0];
rz(-0.5254566) q[1];
sx q[1];
rz(-2.2976687) q[1];
sx q[1];
rz(-2.1702683) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0063365) q[0];
sx q[0];
rz(-0.86733666) q[0];
sx q[0];
rz(-1.9869861) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35324706) q[2];
sx q[2];
rz(-1.1804891) q[2];
sx q[2];
rz(2.3164904) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7156083) q[1];
sx q[1];
rz(-2.2930745) q[1];
sx q[1];
rz(-0.76894185) q[1];
rz(-pi) q[2];
rz(0.97782683) q[3];
sx q[3];
rz(-0.37245107) q[3];
sx q[3];
rz(2.1725871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1725258) q[2];
sx q[2];
rz(-2.064164) q[2];
sx q[2];
rz(1.7871008) q[2];
rz(-2.5403678) q[3];
sx q[3];
rz(-1.0433334) q[3];
sx q[3];
rz(1.5662947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4300267) q[0];
sx q[0];
rz(-0.30421782) q[0];
sx q[0];
rz(-2.8416204) q[0];
rz(0.9777588) q[1];
sx q[1];
rz(-2.561196) q[1];
sx q[1];
rz(-0.93394867) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4276036) q[0];
sx q[0];
rz(-1.4311218) q[0];
sx q[0];
rz(0.33433716) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63542346) q[2];
sx q[2];
rz(-0.72180702) q[2];
sx q[2];
rz(-1.8052764) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4377714) q[1];
sx q[1];
rz(-1.3385845) q[1];
sx q[1];
rz(2.3610568) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4277064) q[3];
sx q[3];
rz(-2.6306804) q[3];
sx q[3];
rz(0.78277222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0710435) q[2];
sx q[2];
rz(-1.0087548) q[2];
sx q[2];
rz(-2.0096931) q[2];
rz(1.8848568) q[3];
sx q[3];
rz(-0.83135024) q[3];
sx q[3];
rz(1.7251714) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8202332) q[0];
sx q[0];
rz(-1.8978523) q[0];
sx q[0];
rz(0.61068049) q[0];
rz(2.3848379) q[1];
sx q[1];
rz(-0.38833955) q[1];
sx q[1];
rz(-1.6441708) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9349139) q[0];
sx q[0];
rz(-1.6948218) q[0];
sx q[0];
rz(-1.5962315) q[0];
x q[1];
rz(2.7779105) q[2];
sx q[2];
rz(-1.1712345) q[2];
sx q[2];
rz(-1.3043208) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.25631902) q[1];
sx q[1];
rz(-2.4838964) q[1];
sx q[1];
rz(2.9672876) q[1];
rz(-pi) q[2];
rz(-1.4704711) q[3];
sx q[3];
rz(-0.85160321) q[3];
sx q[3];
rz(0.49121414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.18621592) q[2];
sx q[2];
rz(-0.66408855) q[2];
sx q[2];
rz(-2.5353298) q[2];
rz(-0.22054211) q[3];
sx q[3];
rz(-1.3371779) q[3];
sx q[3];
rz(0.36673275) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99899387) q[0];
sx q[0];
rz(-1.6479011) q[0];
sx q[0];
rz(-2.2338474) q[0];
rz(-1.6084464) q[1];
sx q[1];
rz(-2.1870859) q[1];
sx q[1];
rz(0.018891637) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9241727) q[0];
sx q[0];
rz(-1.1013563) q[0];
sx q[0];
rz(-1.7007692) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8420686) q[2];
sx q[2];
rz(-1.5148485) q[2];
sx q[2];
rz(2.6331226) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0844473) q[1];
sx q[1];
rz(-1.6020163) q[1];
sx q[1];
rz(-0.46226661) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1891892) q[3];
sx q[3];
rz(-0.84238543) q[3];
sx q[3];
rz(2.4579508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3733526) q[2];
sx q[2];
rz(-1.588593) q[2];
sx q[2];
rz(2.8949883) q[2];
rz(-1.2982347) q[3];
sx q[3];
rz(-1.6876829) q[3];
sx q[3];
rz(1.8728144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3835555) q[0];
sx q[0];
rz(-1.0724496) q[0];
sx q[0];
rz(-2.2776336) q[0];
rz(-1.2147238) q[1];
sx q[1];
rz(-2.0236969) q[1];
sx q[1];
rz(2.5606959) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6075077) q[0];
sx q[0];
rz(-2.0568741) q[0];
sx q[0];
rz(-2.3838504) q[0];
rz(2.3071981) q[2];
sx q[2];
rz(-1.5992924) q[2];
sx q[2];
rz(1.7557276) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.249596) q[1];
sx q[1];
rz(-0.70887762) q[1];
sx q[1];
rz(-2.6061771) q[1];
rz(-0.27972273) q[3];
sx q[3];
rz(-0.7991074) q[3];
sx q[3];
rz(-0.38995686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2555799) q[2];
sx q[2];
rz(-2.6801127) q[2];
sx q[2];
rz(-2.9542921) q[2];
rz(2.1646132) q[3];
sx q[3];
rz(-1.6410442) q[3];
sx q[3];
rz(-1.2788844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63000694) q[0];
sx q[0];
rz(-2.4749909) q[0];
sx q[0];
rz(-1.8774207) q[0];
rz(0.97201792) q[1];
sx q[1];
rz(-2.3897901) q[1];
sx q[1];
rz(1.1690296) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80315932) q[0];
sx q[0];
rz(-1.2375323) q[0];
sx q[0];
rz(0.14305556) q[0];
rz(1.9355785) q[2];
sx q[2];
rz(-1.9309461) q[2];
sx q[2];
rz(-0.84128887) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.15195512) q[1];
sx q[1];
rz(-2.440976) q[1];
sx q[1];
rz(1.7331657) q[1];
rz(-pi) q[2];
rz(2.8253978) q[3];
sx q[3];
rz(-1.4619816) q[3];
sx q[3];
rz(-1.8121882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6675889) q[2];
sx q[2];
rz(-1.7550125) q[2];
sx q[2];
rz(2.7645195) q[2];
rz(-2.9178879) q[3];
sx q[3];
rz(-0.5880028) q[3];
sx q[3];
rz(-1.2288176) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8857464) q[0];
sx q[0];
rz(-1.7485913) q[0];
sx q[0];
rz(2.8843256) q[0];
rz(-2.5480351) q[1];
sx q[1];
rz(-2.3496353) q[1];
sx q[1];
rz(-1.4812462) q[1];
rz(-0.72102265) q[2];
sx q[2];
rz(-0.99698721) q[2];
sx q[2];
rz(0.1884603) q[2];
rz(1.9695342) q[3];
sx q[3];
rz(-1.8137365) q[3];
sx q[3];
rz(0.89567281) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
