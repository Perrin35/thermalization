OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.18908137) q[0];
sx q[0];
rz(3.014325) q[0];
sx q[0];
rz(11.57796) q[0];
rz(2.610511) q[1];
sx q[1];
rz(-0.34871066) q[1];
sx q[1];
rz(-1.7928064) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3492337) q[0];
sx q[0];
rz(-1.4548737) q[0];
sx q[0];
rz(-0.16616343) q[0];
x q[1];
rz(-2.8785107) q[2];
sx q[2];
rz(-1.4882659) q[2];
sx q[2];
rz(-2.6334327) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8475436) q[1];
sx q[1];
rz(-0.8609035) q[1];
sx q[1];
rz(2.3178046) q[1];
x q[2];
rz(1.7705671) q[3];
sx q[3];
rz(-1.3817182) q[3];
sx q[3];
rz(-0.13669361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3794043) q[2];
sx q[2];
rz(-2.4193802) q[2];
sx q[2];
rz(-0.34525004) q[2];
rz(-2.9521862) q[3];
sx q[3];
rz(-0.49049401) q[3];
sx q[3];
rz(2.1729443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9703366) q[0];
sx q[0];
rz(-0.67983627) q[0];
sx q[0];
rz(-0.41980699) q[0];
rz(-0.35821113) q[1];
sx q[1];
rz(-2.5102291) q[1];
sx q[1];
rz(1.0158687) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.779218) q[0];
sx q[0];
rz(-1.172116) q[0];
sx q[0];
rz(0.53074145) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7989356) q[2];
sx q[2];
rz(-1.3698789) q[2];
sx q[2];
rz(1.9659496) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.8034755) q[1];
sx q[1];
rz(-2.7403767) q[1];
sx q[1];
rz(0.34715279) q[1];
rz(-2.696051) q[3];
sx q[3];
rz(-1.1098776) q[3];
sx q[3];
rz(-1.9259491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7887855) q[2];
sx q[2];
rz(-0.31987) q[2];
sx q[2];
rz(-1.1091728) q[2];
rz(-2.7099113) q[3];
sx q[3];
rz(-2.7717398) q[3];
sx q[3];
rz(-0.060401827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.1576841) q[0];
sx q[0];
rz(-1.2976054) q[0];
sx q[0];
rz(-2.2614959) q[0];
rz(-1.2616715) q[1];
sx q[1];
rz(-2.547956) q[1];
sx q[1];
rz(-0.18149158) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49849579) q[0];
sx q[0];
rz(-0.65319121) q[0];
sx q[0];
rz(2.3343587) q[0];
rz(0.65012765) q[2];
sx q[2];
rz(-2.1047154) q[2];
sx q[2];
rz(0.18661737) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.03141244) q[1];
sx q[1];
rz(-1.2978147) q[1];
sx q[1];
rz(0.024755342) q[1];
rz(-pi) q[2];
rz(-1.3550024) q[3];
sx q[3];
rz(-1.0991691) q[3];
sx q[3];
rz(2.140688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.024753831) q[2];
sx q[2];
rz(-1.1081868) q[2];
sx q[2];
rz(1.7987569) q[2];
rz(-1.3252307) q[3];
sx q[3];
rz(-2.467005) q[3];
sx q[3];
rz(-1.0935812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6670068) q[0];
sx q[0];
rz(-0.51303595) q[0];
sx q[0];
rz(-0.68853199) q[0];
rz(-2.9225598) q[1];
sx q[1];
rz(-0.34866798) q[1];
sx q[1];
rz(-2.9825488) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6200136) q[0];
sx q[0];
rz(-1.0431246) q[0];
sx q[0];
rz(2.2844727) q[0];
rz(-0.7775457) q[2];
sx q[2];
rz(-0.50707196) q[2];
sx q[2];
rz(0.61447243) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5608983) q[1];
sx q[1];
rz(-1.782522) q[1];
sx q[1];
rz(-3.1161518) q[1];
rz(-0.36985107) q[3];
sx q[3];
rz(-0.55579138) q[3];
sx q[3];
rz(-0.37214798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.29545414) q[2];
sx q[2];
rz(-1.4825772) q[2];
sx q[2];
rz(0.42567483) q[2];
rz(3.1221636) q[3];
sx q[3];
rz(-2.9255376) q[3];
sx q[3];
rz(2.4470636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.46538019) q[0];
sx q[0];
rz(-0.0066444962) q[0];
sx q[0];
rz(-0.12839578) q[0];
rz(-2.45576) q[1];
sx q[1];
rz(-0.99629712) q[1];
sx q[1];
rz(2.593186) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0869285) q[0];
sx q[0];
rz(-1.5318649) q[0];
sx q[0];
rz(-0.995308) q[0];
rz(-pi) q[1];
rz(-1.3850645) q[2];
sx q[2];
rz(-0.81132946) q[2];
sx q[2];
rz(-0.059543691) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4244528) q[1];
sx q[1];
rz(-1.5100749) q[1];
sx q[1];
rz(-1.1007924) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3324276) q[3];
sx q[3];
rz(-2.7276464) q[3];
sx q[3];
rz(0.64869374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2942865) q[2];
sx q[2];
rz(-2.8319478) q[2];
sx q[2];
rz(-1.4667286) q[2];
rz(1.0855801) q[3];
sx q[3];
rz(-1.3054566) q[3];
sx q[3];
rz(-0.87695688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0393031) q[0];
sx q[0];
rz(-1.9474494) q[0];
sx q[0];
rz(-2.0423245) q[0];
rz(0.33637235) q[1];
sx q[1];
rz(-0.66434324) q[1];
sx q[1];
rz(0.99463314) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1625527) q[0];
sx q[0];
rz(-2.131922) q[0];
sx q[0];
rz(-2.8812376) q[0];
rz(2.6380831) q[2];
sx q[2];
rz(-1.5739872) q[2];
sx q[2];
rz(-1.9838331) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5632538) q[1];
sx q[1];
rz(-1.564338) q[1];
sx q[1];
rz(1.2328813) q[1];
x q[2];
rz(1.1676222) q[3];
sx q[3];
rz(-0.91114985) q[3];
sx q[3];
rz(-0.50658222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2991335) q[2];
sx q[2];
rz(-0.9023388) q[2];
sx q[2];
rz(-0.4449521) q[2];
rz(-2.3343202) q[3];
sx q[3];
rz(-0.10891309) q[3];
sx q[3];
rz(2.3256433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.826236) q[0];
sx q[0];
rz(-2.5671791) q[0];
sx q[0];
rz(2.2454967) q[0];
rz(-2.232146) q[1];
sx q[1];
rz(-2.8912631) q[1];
sx q[1];
rz(3.0665841) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79338193) q[0];
sx q[0];
rz(-2.777034) q[0];
sx q[0];
rz(-1.7478463) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3415065) q[2];
sx q[2];
rz(-1.9817838) q[2];
sx q[2];
rz(-1.6192186) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8304886) q[1];
sx q[1];
rz(-1.5058335) q[1];
sx q[1];
rz(2.225609) q[1];
x q[2];
rz(2.1771031) q[3];
sx q[3];
rz(-0.74092591) q[3];
sx q[3];
rz(0.20674202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.16671495) q[2];
sx q[2];
rz(-2.0090943) q[2];
sx q[2];
rz(1.9141076) q[2];
rz(0.04315367) q[3];
sx q[3];
rz(-1.6601325) q[3];
sx q[3];
rz(-0.22668049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9880923) q[0];
sx q[0];
rz(-0.38689125) q[0];
sx q[0];
rz(0.41326997) q[0];
rz(-0.55832541) q[1];
sx q[1];
rz(-2.3007326) q[1];
sx q[1];
rz(-2.4024898) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16427134) q[0];
sx q[0];
rz(-0.56557206) q[0];
sx q[0];
rz(0.092202734) q[0];
rz(-pi) q[1];
rz(-1.1147538) q[2];
sx q[2];
rz(-1.5762435) q[2];
sx q[2];
rz(-1.3378439) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0308471) q[1];
sx q[1];
rz(-1.7304286) q[1];
sx q[1];
rz(-0.84407945) q[1];
rz(-1.7750752) q[3];
sx q[3];
rz(-0.26976997) q[3];
sx q[3];
rz(-1.1360053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83453137) q[2];
sx q[2];
rz(-0.31690159) q[2];
sx q[2];
rz(2.0397662) q[2];
rz(2.7448591) q[3];
sx q[3];
rz(-1.705403) q[3];
sx q[3];
rz(2.326899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48801625) q[0];
sx q[0];
rz(-0.62960136) q[0];
sx q[0];
rz(0.6189515) q[0];
rz(2.1221819) q[1];
sx q[1];
rz(-1.5827725) q[1];
sx q[1];
rz(2.6143262) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70452481) q[0];
sx q[0];
rz(-1.2906115) q[0];
sx q[0];
rz(1.8591451) q[0];
x q[1];
rz(2.451532) q[2];
sx q[2];
rz(-2.0482716) q[2];
sx q[2];
rz(2.8183187) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.077721715) q[1];
sx q[1];
rz(-1.9943024) q[1];
sx q[1];
rz(-1.0788171) q[1];
x q[2];
rz(-0.23993581) q[3];
sx q[3];
rz(-0.92471189) q[3];
sx q[3];
rz(2.5480888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8573389) q[2];
sx q[2];
rz(-1.3353835) q[2];
sx q[2];
rz(0.93150345) q[2];
rz(2.3305317) q[3];
sx q[3];
rz(-0.61412007) q[3];
sx q[3];
rz(-1.567599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60944027) q[0];
sx q[0];
rz(-2.092822) q[0];
sx q[0];
rz(0.47927454) q[0];
rz(-2.2562064) q[1];
sx q[1];
rz(-2.482174) q[1];
sx q[1];
rz(-0.17818174) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70885926) q[0];
sx q[0];
rz(-1.0646001) q[0];
sx q[0];
rz(-2.3082255) q[0];
rz(-pi) q[1];
rz(-0.41583305) q[2];
sx q[2];
rz(-0.73199474) q[2];
sx q[2];
rz(-2.7351565) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.61160117) q[1];
sx q[1];
rz(-2.0833263) q[1];
sx q[1];
rz(0.26228735) q[1];
x q[2];
rz(-0.43042572) q[3];
sx q[3];
rz(-1.6820551) q[3];
sx q[3];
rz(2.984897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6726154) q[2];
sx q[2];
rz(-0.72300935) q[2];
sx q[2];
rz(0.82328063) q[2];
rz(-3.0205884) q[3];
sx q[3];
rz(-0.76366097) q[3];
sx q[3];
rz(-2.1774489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.9085893) q[0];
sx q[0];
rz(-2.067726) q[0];
sx q[0];
rz(2.4698972) q[0];
rz(1.9227149) q[1];
sx q[1];
rz(-1.8119443) q[1];
sx q[1];
rz(-1.3655566) q[1];
rz(-3.1115816) q[2];
sx q[2];
rz(-1.9297615) q[2];
sx q[2];
rz(2.161138) q[2];
rz(-0.34006313) q[3];
sx q[3];
rz(-2.1652514) q[3];
sx q[3];
rz(1.3596331) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
