OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.31792274) q[0];
sx q[0];
rz(-1.5225141) q[0];
sx q[0];
rz(-0.063152753) q[0];
rz(-1.3178648) q[1];
sx q[1];
rz(-2.7422819) q[1];
sx q[1];
rz(0.77594405) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0948333) q[0];
sx q[0];
rz(-1.6164808) q[0];
sx q[0];
rz(1.5976357) q[0];
x q[1];
rz(0.93504436) q[2];
sx q[2];
rz(-1.9823977) q[2];
sx q[2];
rz(-1.7807478) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5580485) q[1];
sx q[1];
rz(-2.2698102) q[1];
sx q[1];
rz(-0.2322766) q[1];
rz(-pi) q[2];
rz(1.3315866) q[3];
sx q[3];
rz(-2.0612217) q[3];
sx q[3];
rz(-1.8338721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6905288) q[2];
sx q[2];
rz(-1.3961926) q[2];
sx q[2];
rz(-0.53831354) q[2];
rz(-2.8484143) q[3];
sx q[3];
rz(-1.5436951) q[3];
sx q[3];
rz(-1.258491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20741589) q[0];
sx q[0];
rz(-0.84686142) q[0];
sx q[0];
rz(2.9127981) q[0];
rz(3.0711807) q[1];
sx q[1];
rz(-1.8682624) q[1];
sx q[1];
rz(1.061903) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.966514) q[0];
sx q[0];
rz(-2.0850777) q[0];
sx q[0];
rz(-1.2660591) q[0];
x q[1];
rz(0.76490374) q[2];
sx q[2];
rz(-2.5203369) q[2];
sx q[2];
rz(2.3111631) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.50256601) q[1];
sx q[1];
rz(-0.4430534) q[1];
sx q[1];
rz(0.73213099) q[1];
rz(-pi) q[2];
rz(-0.44159378) q[3];
sx q[3];
rz(-0.46735763) q[3];
sx q[3];
rz(2.2668214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.215302) q[2];
sx q[2];
rz(-2.6965202) q[2];
sx q[2];
rz(-0.11623795) q[2];
rz(-0.91533533) q[3];
sx q[3];
rz(-1.4696308) q[3];
sx q[3];
rz(-2.5168929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4794469) q[0];
sx q[0];
rz(-1.5818469) q[0];
sx q[0];
rz(2.8060655) q[0];
rz(-1.7299995) q[1];
sx q[1];
rz(-2.2970707) q[1];
sx q[1];
rz(-1.1988877) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39349651) q[0];
sx q[0];
rz(-0.61433799) q[0];
sx q[0];
rz(-1.0636281) q[0];
x q[1];
rz(0.54036083) q[2];
sx q[2];
rz(-0.66695222) q[2];
sx q[2];
rz(-1.6467384) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.88628593) q[1];
sx q[1];
rz(-2.3706808) q[1];
sx q[1];
rz(0.16718276) q[1];
rz(0.24925225) q[3];
sx q[3];
rz(-1.0947168) q[3];
sx q[3];
rz(0.92633807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.18232839) q[2];
sx q[2];
rz(-2.3458643) q[2];
sx q[2];
rz(-1.1129334) q[2];
rz(0.5111323) q[3];
sx q[3];
rz(-1.7685578) q[3];
sx q[3];
rz(0.50890499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48753259) q[0];
sx q[0];
rz(-3.0656116) q[0];
sx q[0];
rz(2.7677166) q[0];
rz(-1.9591676) q[1];
sx q[1];
rz(-1.1610718) q[1];
sx q[1];
rz(-0.16071308) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041695853) q[0];
sx q[0];
rz(-1.4879491) q[0];
sx q[0];
rz(-3.003398) q[0];
rz(-pi) q[1];
rz(2.3473927) q[2];
sx q[2];
rz(-2.080285) q[2];
sx q[2];
rz(-1.8812665) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4121033) q[1];
sx q[1];
rz(-2.0827052) q[1];
sx q[1];
rz(0.44513925) q[1];
rz(-pi) q[2];
rz(2.8237958) q[3];
sx q[3];
rz(-1.2850015) q[3];
sx q[3];
rz(0.60493776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5339727) q[2];
sx q[2];
rz(-2.0325568) q[2];
sx q[2];
rz(-0.71471659) q[2];
rz(-0.25501069) q[3];
sx q[3];
rz(-1.3304354) q[3];
sx q[3];
rz(2.3431006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.052499972) q[0];
sx q[0];
rz(-2.8369501) q[0];
sx q[0];
rz(-2.3783045) q[0];
rz(2.8870562) q[1];
sx q[1];
rz(-1.7691879) q[1];
sx q[1];
rz(1.4483784) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0453706) q[0];
sx q[0];
rz(-1.9217446) q[0];
sx q[0];
rz(0.36349067) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4035475) q[2];
sx q[2];
rz(-2.3838701) q[2];
sx q[2];
rz(0.12239419) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4896488) q[1];
sx q[1];
rz(-1.3931511) q[1];
sx q[1];
rz(0.62477115) q[1];
rz(-2.4302825) q[3];
sx q[3];
rz(-1.2625202) q[3];
sx q[3];
rz(3.0570249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.86108565) q[2];
sx q[2];
rz(-2.0528767) q[2];
sx q[2];
rz(-2.1121934) q[2];
rz(-1.6086802) q[3];
sx q[3];
rz(-2.2423988) q[3];
sx q[3];
rz(-0.6338343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7480046) q[0];
sx q[0];
rz(-0.91359502) q[0];
sx q[0];
rz(-2.9314801) q[0];
rz(-2.4402319) q[1];
sx q[1];
rz(-1.7751834) q[1];
sx q[1];
rz(-2.0393541) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85066807) q[0];
sx q[0];
rz(-2.0966665) q[0];
sx q[0];
rz(0.53914244) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5731235) q[2];
sx q[2];
rz(-2.2330052) q[2];
sx q[2];
rz(2.4072449) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7581492) q[1];
sx q[1];
rz(-2.2510937) q[1];
sx q[1];
rz(-1.1636984) q[1];
rz(-pi) q[2];
rz(0.26138283) q[3];
sx q[3];
rz(-0.55703516) q[3];
sx q[3];
rz(1.5294531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.29409757) q[2];
sx q[2];
rz(-1.558344) q[2];
sx q[2];
rz(1.4964649) q[2];
rz(1.7485113) q[3];
sx q[3];
rz(-2.2025509) q[3];
sx q[3];
rz(-2.4699672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74712718) q[0];
sx q[0];
rz(-1.324993) q[0];
sx q[0];
rz(0.4766683) q[0];
rz(-1.7851625) q[1];
sx q[1];
rz(-1.2973659) q[1];
sx q[1];
rz(-0.93723255) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.40537) q[0];
sx q[0];
rz(-2.4511466) q[0];
sx q[0];
rz(-1.5140945) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4235882) q[2];
sx q[2];
rz(-1.5239626) q[2];
sx q[2];
rz(-0.61344922) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7644219) q[1];
sx q[1];
rz(-0.57767361) q[1];
sx q[1];
rz(-2.38297) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6680596) q[3];
sx q[3];
rz(-2.0174806) q[3];
sx q[3];
rz(0.47998715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6095907) q[2];
sx q[2];
rz(-1.3401745) q[2];
sx q[2];
rz(-1.0756005) q[2];
rz(2.7514451) q[3];
sx q[3];
rz(-1.8307999) q[3];
sx q[3];
rz(2.3384317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53668642) q[0];
sx q[0];
rz(-2.0669879) q[0];
sx q[0];
rz(-2.5965776) q[0];
rz(2.1489428) q[1];
sx q[1];
rz(-0.98572171) q[1];
sx q[1];
rz(1.4535905) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2611985) q[0];
sx q[0];
rz(-0.81982938) q[0];
sx q[0];
rz(0.96078787) q[0];
rz(1.7033061) q[2];
sx q[2];
rz(-1.514666) q[2];
sx q[2];
rz(-1.1289589) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7632227) q[1];
sx q[1];
rz(-1.6838131) q[1];
sx q[1];
rz(-2.3341353) q[1];
rz(1.5239632) q[3];
sx q[3];
rz(-1.5107656) q[3];
sx q[3];
rz(-1.6215274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.40889007) q[2];
sx q[2];
rz(-1.3971034) q[2];
sx q[2];
rz(-0.2962386) q[2];
rz(0.57684165) q[3];
sx q[3];
rz(-2.7448765) q[3];
sx q[3];
rz(1.860994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7883564) q[0];
sx q[0];
rz(-2.3529973) q[0];
sx q[0];
rz(2.9175135) q[0];
rz(-2.5375986) q[1];
sx q[1];
rz(-0.9915587) q[1];
sx q[1];
rz(-1.6557065) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1093134) q[0];
sx q[0];
rz(-1.8173738) q[0];
sx q[0];
rz(0.37301491) q[0];
x q[1];
rz(0.15023709) q[2];
sx q[2];
rz(-1.7662107) q[2];
sx q[2];
rz(0.11997151) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.60843492) q[1];
sx q[1];
rz(-2.4704993) q[1];
sx q[1];
rz(-0.41678269) q[1];
rz(1.8698244) q[3];
sx q[3];
rz(-1.097603) q[3];
sx q[3];
rz(-0.92834559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.47933444) q[2];
sx q[2];
rz(-2.0749638) q[2];
sx q[2];
rz(-2.7545641) q[2];
rz(-2.8582063) q[3];
sx q[3];
rz(-2.3895388) q[3];
sx q[3];
rz(-2.7391105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0249483) q[0];
sx q[0];
rz(-0.270917) q[0];
sx q[0];
rz(-1.1234294) q[0];
rz(-1.4943538) q[1];
sx q[1];
rz(-1.6554183) q[1];
sx q[1];
rz(-0.87596881) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7520053) q[0];
sx q[0];
rz(-1.8266062) q[0];
sx q[0];
rz(-2.385456) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89089762) q[2];
sx q[2];
rz(-2.4844137) q[2];
sx q[2];
rz(2.3466722) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0735262) q[1];
sx q[1];
rz(-1.7896393) q[1];
sx q[1];
rz(-0.88394758) q[1];
rz(-1.749529) q[3];
sx q[3];
rz(-1.1471841) q[3];
sx q[3];
rz(2.7152747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4103849) q[2];
sx q[2];
rz(-1.7816252) q[2];
sx q[2];
rz(0.13710055) q[2];
rz(1.3827263) q[3];
sx q[3];
rz(-2.2189249) q[3];
sx q[3];
rz(1.2285129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41505861) q[0];
sx q[0];
rz(-0.4701604) q[0];
sx q[0];
rz(-0.62248019) q[0];
rz(-0.38901916) q[1];
sx q[1];
rz(-0.35094378) q[1];
sx q[1];
rz(2.6019179) q[1];
rz(0.71826886) q[2];
sx q[2];
rz(-0.61189883) q[2];
sx q[2];
rz(1.9552678) q[2];
rz(2.9403654) q[3];
sx q[3];
rz(-1.22898) q[3];
sx q[3];
rz(-2.1635319) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
