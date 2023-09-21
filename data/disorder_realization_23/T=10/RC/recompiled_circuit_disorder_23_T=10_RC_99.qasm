OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2620579) q[0];
sx q[0];
rz(-1.7320002) q[0];
sx q[0];
rz(1.4341266) q[0];
rz(0.6342451) q[1];
sx q[1];
rz(-2.5399962) q[1];
sx q[1];
rz(2.7231725) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2375862) q[0];
sx q[0];
rz(-1.307784) q[0];
sx q[0];
rz(-2.0531274) q[0];
rz(-pi) q[1];
rz(-1.4346052) q[2];
sx q[2];
rz(-1.4601267) q[2];
sx q[2];
rz(-1.1510804) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.13222209) q[1];
sx q[1];
rz(-1.3271866) q[1];
sx q[1];
rz(1.1019215) q[1];
rz(-pi) q[2];
rz(-1.1316142) q[3];
sx q[3];
rz(-1.7712799) q[3];
sx q[3];
rz(0.85103121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.3216386) q[2];
sx q[2];
rz(-1.7724089) q[2];
sx q[2];
rz(2.3036172) q[2];
rz(-0.49301246) q[3];
sx q[3];
rz(-2.8686782) q[3];
sx q[3];
rz(0.078991927) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3943966) q[0];
sx q[0];
rz(-0.72421873) q[0];
sx q[0];
rz(-1.2778506) q[0];
rz(0.17678075) q[1];
sx q[1];
rz(-1.3143833) q[1];
sx q[1];
rz(2.7094254) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1457739) q[0];
sx q[0];
rz(-2.5403025) q[0];
sx q[0];
rz(2.6390618) q[0];
rz(1.1439267) q[2];
sx q[2];
rz(-1.9383213) q[2];
sx q[2];
rz(1.9232242) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4914815) q[1];
sx q[1];
rz(-1.1307798) q[1];
sx q[1];
rz(3.0122709) q[1];
rz(-pi) q[2];
rz(-3.1167332) q[3];
sx q[3];
rz(-2.1356574) q[3];
sx q[3];
rz(-0.75627518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5923578) q[2];
sx q[2];
rz(-1.8775512) q[2];
sx q[2];
rz(-2.6548927) q[2];
rz(-1.3782079) q[3];
sx q[3];
rz(-1.2599726) q[3];
sx q[3];
rz(-2.6087705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.040722672) q[0];
sx q[0];
rz(-0.7464872) q[0];
sx q[0];
rz(-0.41734636) q[0];
rz(1.4886645) q[1];
sx q[1];
rz(-2.5960943) q[1];
sx q[1];
rz(0.506385) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7931472) q[0];
sx q[0];
rz(-1.9217102) q[0];
sx q[0];
rz(0.1883513) q[0];
rz(-pi) q[1];
rz(2.8767013) q[2];
sx q[2];
rz(-0.75220097) q[2];
sx q[2];
rz(-1.0571935) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8295146) q[1];
sx q[1];
rz(-1.2631589) q[1];
sx q[1];
rz(-2.679146) q[1];
rz(1.1261602) q[3];
sx q[3];
rz(-2.7362842) q[3];
sx q[3];
rz(-2.0733548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4425519) q[2];
sx q[2];
rz(-2.6802345) q[2];
sx q[2];
rz(2.55012) q[2];
rz(0.58602035) q[3];
sx q[3];
rz(-1.2089217) q[3];
sx q[3];
rz(-1.4311786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9716924) q[0];
sx q[0];
rz(-0.62830347) q[0];
sx q[0];
rz(2.3024094) q[0];
rz(3.1160141) q[1];
sx q[1];
rz(-0.69568101) q[1];
sx q[1];
rz(1.5930088) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4913113) q[0];
sx q[0];
rz(-0.45931739) q[0];
sx q[0];
rz(-1.7772872) q[0];
x q[1];
rz(-0.80231248) q[2];
sx q[2];
rz(-2.1244086) q[2];
sx q[2];
rz(2.4137036) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7566484) q[1];
sx q[1];
rz(-1.952991) q[1];
sx q[1];
rz(-2.3734943) q[1];
rz(-pi) q[2];
rz(2.5769916) q[3];
sx q[3];
rz(-0.49800107) q[3];
sx q[3];
rz(-2.6548487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.30535355) q[2];
sx q[2];
rz(-0.88399115) q[2];
sx q[2];
rz(0.099686064) q[2];
rz(2.1827407) q[3];
sx q[3];
rz(-1.3189664) q[3];
sx q[3];
rz(-1.8166186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5754159) q[0];
sx q[0];
rz(-1.382099) q[0];
sx q[0];
rz(2.8856522) q[0];
rz(-0.4610962) q[1];
sx q[1];
rz(-1.0436811) q[1];
sx q[1];
rz(-2.3815313) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041199112) q[0];
sx q[0];
rz(-1.4677591) q[0];
sx q[0];
rz(0.11739199) q[0];
x q[1];
rz(-1.4270093) q[2];
sx q[2];
rz(-0.41848768) q[2];
sx q[2];
rz(2.8029122) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6872245) q[1];
sx q[1];
rz(-1.4441274) q[1];
sx q[1];
rz(1.6997937) q[1];
rz(-pi) q[2];
rz(-2.2171668) q[3];
sx q[3];
rz(-1.2741718) q[3];
sx q[3];
rz(0.38277205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5710859) q[2];
sx q[2];
rz(-1.7692302) q[2];
sx q[2];
rz(-0.6742397) q[2];
rz(2.9267866) q[3];
sx q[3];
rz(-0.45682296) q[3];
sx q[3];
rz(-3.1242483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.6102585) q[0];
sx q[0];
rz(-1.4704309) q[0];
sx q[0];
rz(1.1791139) q[0];
rz(2.9367661) q[1];
sx q[1];
rz(-2.3463459) q[1];
sx q[1];
rz(1.0669605) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2039295) q[0];
sx q[0];
rz(-1.5852889) q[0];
sx q[0];
rz(-0.020676215) q[0];
rz(2.3382171) q[2];
sx q[2];
rz(-0.57758812) q[2];
sx q[2];
rz(-0.20882777) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3031591) q[1];
sx q[1];
rz(-1.5232956) q[1];
sx q[1];
rz(0.82364239) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9818929) q[3];
sx q[3];
rz(-0.80280639) q[3];
sx q[3];
rz(1.6646977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.431488) q[2];
sx q[2];
rz(-1.2892712) q[2];
sx q[2];
rz(-2.5816494) q[2];
rz(0.7263178) q[3];
sx q[3];
rz(-2.8328219) q[3];
sx q[3];
rz(2.8360951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2840435) q[0];
sx q[0];
rz(-2.5950268) q[0];
sx q[0];
rz(1.7204826) q[0];
rz(-0.20206085) q[1];
sx q[1];
rz(-1.4338564) q[1];
sx q[1];
rz(-2.2834159) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35408033) q[0];
sx q[0];
rz(-1.4424099) q[0];
sx q[0];
rz(-1.7187198) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.467534) q[2];
sx q[2];
rz(-0.41245663) q[2];
sx q[2];
rz(-2.7445284) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9550025) q[1];
sx q[1];
rz(-1.5555256) q[1];
sx q[1];
rz(-0.028970684) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57111994) q[3];
sx q[3];
rz(-1.8652328) q[3];
sx q[3];
rz(-1.4049698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8537366) q[2];
sx q[2];
rz(-2.6617472) q[2];
sx q[2];
rz(-1.8161592) q[2];
rz(0.89007968) q[3];
sx q[3];
rz(-1.1471014) q[3];
sx q[3];
rz(-2.1530698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-1.4637852) q[0];
sx q[0];
rz(-2.5690434) q[0];
sx q[0];
rz(2.7668787) q[0];
rz(2.162714) q[1];
sx q[1];
rz(-2.4596877) q[1];
sx q[1];
rz(1.7920866) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6222854) q[0];
sx q[0];
rz(-1.1847727) q[0];
sx q[0];
rz(0.65667721) q[0];
rz(-pi) q[1];
rz(-0.026272341) q[2];
sx q[2];
rz(-2.4233344) q[2];
sx q[2];
rz(0.19933137) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.046557758) q[1];
sx q[1];
rz(-1.7525502) q[1];
sx q[1];
rz(0.05038105) q[1];
rz(-pi) q[2];
rz(-0.41534822) q[3];
sx q[3];
rz(-0.58598622) q[3];
sx q[3];
rz(-1.8893482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4902041) q[2];
sx q[2];
rz(-2.3888402) q[2];
sx q[2];
rz(0.46869579) q[2];
rz(1.1941341) q[3];
sx q[3];
rz(-1.8374551) q[3];
sx q[3];
rz(1.7780001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63012183) q[0];
sx q[0];
rz(-2.2352495) q[0];
sx q[0];
rz(-1.8796896) q[0];
rz(2.966554) q[1];
sx q[1];
rz(-1.9997528) q[1];
sx q[1];
rz(1.5375686) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3106829) q[0];
sx q[0];
rz(-1.3498107) q[0];
sx q[0];
rz(-1.882878) q[0];
x q[1];
rz(1.4666918) q[2];
sx q[2];
rz(-1.1541919) q[2];
sx q[2];
rz(-0.57550752) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.21737145) q[1];
sx q[1];
rz(-1.4573759) q[1];
sx q[1];
rz(-0.70593112) q[1];
rz(-pi) q[2];
rz(-0.36112862) q[3];
sx q[3];
rz(-2.0054521) q[3];
sx q[3];
rz(1.1772732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.039915446) q[2];
sx q[2];
rz(-2.1890409) q[2];
sx q[2];
rz(-2.3802479) q[2];
rz(-2.2411761) q[3];
sx q[3];
rz(-0.59949985) q[3];
sx q[3];
rz(0.049023978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.802357) q[0];
sx q[0];
rz(-2.673322) q[0];
sx q[0];
rz(-0.21690579) q[0];
rz(-2.5096109) q[1];
sx q[1];
rz(-1.6501553) q[1];
sx q[1];
rz(2.1868618) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58669421) q[0];
sx q[0];
rz(-2.2838755) q[0];
sx q[0];
rz(-1.9999534) q[0];
rz(-pi) q[1];
x q[1];
rz(0.63072272) q[2];
sx q[2];
rz(-2.4927757) q[2];
sx q[2];
rz(-2.0231252) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6001544) q[1];
sx q[1];
rz(-0.46240515) q[1];
sx q[1];
rz(2.9701783) q[1];
x q[2];
rz(-2.1438164) q[3];
sx q[3];
rz(-0.7191092) q[3];
sx q[3];
rz(2.5620808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1251936) q[2];
sx q[2];
rz(-1.9026326) q[2];
sx q[2];
rz(-2.1968502) q[2];
rz(2.7567806) q[3];
sx q[3];
rz(-1.1184357) q[3];
sx q[3];
rz(2.184536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64086296) q[0];
sx q[0];
rz(-0.50518112) q[0];
sx q[0];
rz(1.5541979) q[0];
rz(0.8846994) q[1];
sx q[1];
rz(-0.90507602) q[1];
sx q[1];
rz(-0.25837635) q[1];
rz(0.89818556) q[2];
sx q[2];
rz(-2.2834416) q[2];
sx q[2];
rz(-1.7590547) q[2];
rz(-0.20529071) q[3];
sx q[3];
rz(-1.1355573) q[3];
sx q[3];
rz(-0.60187403) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
