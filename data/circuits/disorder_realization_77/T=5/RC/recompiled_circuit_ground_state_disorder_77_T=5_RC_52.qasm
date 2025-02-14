OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.90280688) q[0];
sx q[0];
rz(2.0030608) q[0];
sx q[0];
rz(8.3264405) q[0];
rz(1.1605473) q[1];
sx q[1];
rz(-1.1359954) q[1];
sx q[1];
rz(2.5392037) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1404323) q[0];
sx q[0];
rz(-0.70297697) q[0];
sx q[0];
rz(0.70859895) q[0];
rz(-pi) q[1];
rz(-2.4279934) q[2];
sx q[2];
rz(-1.6125792) q[2];
sx q[2];
rz(0.080619911) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5307926) q[1];
sx q[1];
rz(-1.4565174) q[1];
sx q[1];
rz(2.0487993) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6528879) q[3];
sx q[3];
rz(-2.7118976) q[3];
sx q[3];
rz(-1.8703509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41523734) q[2];
sx q[2];
rz(-1.0865612) q[2];
sx q[2];
rz(-1.497867) q[2];
rz(-0.1114397) q[3];
sx q[3];
rz(-1.2846416) q[3];
sx q[3];
rz(-2.1845412) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9621256) q[0];
sx q[0];
rz(-2.2158556) q[0];
sx q[0];
rz(-0.16394462) q[0];
rz(-0.97490772) q[1];
sx q[1];
rz(-2.4539852) q[1];
sx q[1];
rz(-1.499739) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1346171) q[0];
sx q[0];
rz(-2.3762759) q[0];
sx q[0];
rz(1.2045443) q[0];
x q[1];
rz(2.2664505) q[2];
sx q[2];
rz(-2.2457457) q[2];
sx q[2];
rz(-0.067719134) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.837599) q[1];
sx q[1];
rz(-1.493206) q[1];
sx q[1];
rz(-2.4291888) q[1];
x q[2];
rz(-2.9256066) q[3];
sx q[3];
rz(-1.2519998) q[3];
sx q[3];
rz(-2.2270798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.46453309) q[2];
sx q[2];
rz(-0.50731069) q[2];
sx q[2];
rz(2.9110009) q[2];
rz(-2.4584037) q[3];
sx q[3];
rz(-0.69470996) q[3];
sx q[3];
rz(0.75146875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36000073) q[0];
sx q[0];
rz(-1.3104023) q[0];
sx q[0];
rz(2.7202284) q[0];
rz(-1.6257809) q[1];
sx q[1];
rz(-2.6457364) q[1];
sx q[1];
rz(2.1741507) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7712647) q[0];
sx q[0];
rz(-1.2441085) q[0];
sx q[0];
rz(0.26854923) q[0];
rz(-pi) q[1];
rz(2.967796) q[2];
sx q[2];
rz(-2.0325543) q[2];
sx q[2];
rz(-2.6231036) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.9989717) q[1];
sx q[1];
rz(-1.843286) q[1];
sx q[1];
rz(0.15010826) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5021867) q[3];
sx q[3];
rz(-2.5192776) q[3];
sx q[3];
rz(0.50014225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4030054) q[2];
sx q[2];
rz(-1.0120729) q[2];
sx q[2];
rz(2.9834413) q[2];
rz(-0.75913298) q[3];
sx q[3];
rz(-1.682155) q[3];
sx q[3];
rz(-2.7706743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-0.89449969) q[0];
sx q[0];
rz(-2.5573754) q[0];
sx q[0];
rz(-1.1868813) q[0];
rz(3.0184556) q[1];
sx q[1];
rz(-1.5703399) q[1];
sx q[1];
rz(-1.8298967) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0054811) q[0];
sx q[0];
rz(-1.4887047) q[0];
sx q[0];
rz(-1.422922) q[0];
rz(-pi) q[1];
rz(-3.0391215) q[2];
sx q[2];
rz(-1.6471286) q[2];
sx q[2];
rz(-1.9988521) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0801436) q[1];
sx q[1];
rz(-1.176136) q[1];
sx q[1];
rz(0.45545642) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5667241) q[3];
sx q[3];
rz(-1.5713672) q[3];
sx q[3];
rz(-0.79889311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9436403) q[2];
sx q[2];
rz(-1.3200878) q[2];
sx q[2];
rz(2.2894335) q[2];
rz(-0.85396829) q[3];
sx q[3];
rz(-1.9210509) q[3];
sx q[3];
rz(1.0796237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3581486) q[0];
sx q[0];
rz(-1.3842979) q[0];
sx q[0];
rz(-0.070601687) q[0];
rz(-1.0380896) q[1];
sx q[1];
rz(-1.9990653) q[1];
sx q[1];
rz(-2.3322754) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3883078) q[0];
sx q[0];
rz(-1.4986254) q[0];
sx q[0];
rz(0.57467527) q[0];
rz(0.90427601) q[2];
sx q[2];
rz(-2.1993786) q[2];
sx q[2];
rz(1.5360127) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0910803) q[1];
sx q[1];
rz(-1.3111171) q[1];
sx q[1];
rz(2.5804094) q[1];
x q[2];
rz(-2.9633459) q[3];
sx q[3];
rz(-2.2615454) q[3];
sx q[3];
rz(2.8028298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.223083) q[2];
sx q[2];
rz(-0.56879908) q[2];
sx q[2];
rz(-1.9522379) q[2];
rz(-3.0009624) q[3];
sx q[3];
rz(-2.6719002) q[3];
sx q[3];
rz(-2.752059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73914948) q[0];
sx q[0];
rz(-0.91366714) q[0];
sx q[0];
rz(-1.5341349) q[0];
rz(0.88273478) q[1];
sx q[1];
rz(-2.2969756) q[1];
sx q[1];
rz(-2.5260063) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7360216) q[0];
sx q[0];
rz(-1.6628436) q[0];
sx q[0];
rz(-2.7476601) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65290218) q[2];
sx q[2];
rz(-0.98207475) q[2];
sx q[2];
rz(-1.9984286) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5639685) q[1];
sx q[1];
rz(-1.7610487) q[1];
sx q[1];
rz(-3.0012194) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5869523) q[3];
sx q[3];
rz(-2.5211341) q[3];
sx q[3];
rz(2.4718493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.1659871) q[2];
sx q[2];
rz(-0.96419445) q[2];
sx q[2];
rz(-1.1791505) q[2];
rz(-0.21373448) q[3];
sx q[3];
rz(-1.7515747) q[3];
sx q[3];
rz(-2.8960622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8971276) q[0];
sx q[0];
rz(-1.1470969) q[0];
sx q[0];
rz(0.044064673) q[0];
rz(0.38912082) q[1];
sx q[1];
rz(-2.6521284) q[1];
sx q[1];
rz(0.74180952) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.378063) q[0];
sx q[0];
rz(-2.6776095) q[0];
sx q[0];
rz(0.65391175) q[0];
rz(0.82382186) q[2];
sx q[2];
rz(-1.3082817) q[2];
sx q[2];
rz(1.4854625) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9981813) q[1];
sx q[1];
rz(-0.3340958) q[1];
sx q[1];
rz(-0.27004529) q[1];
rz(-pi) q[2];
rz(2.4398285) q[3];
sx q[3];
rz(-1.7265111) q[3];
sx q[3];
rz(-3.0879741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6645708) q[2];
sx q[2];
rz(-1.7844113) q[2];
sx q[2];
rz(0.50194293) q[2];
rz(1.4255514) q[3];
sx q[3];
rz(-2.612096) q[3];
sx q[3];
rz(0.79290032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75941706) q[0];
sx q[0];
rz(-1.4074396) q[0];
sx q[0];
rz(-0.39696804) q[0];
rz(-1.601903) q[1];
sx q[1];
rz(-0.27286467) q[1];
sx q[1];
rz(2.4698965) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3701733) q[0];
sx q[0];
rz(-1.1259997) q[0];
sx q[0];
rz(3.0014787) q[0];
rz(-pi) q[1];
rz(2.7717088) q[2];
sx q[2];
rz(-1.9871759) q[2];
sx q[2];
rz(0.86341349) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.54597832) q[1];
sx q[1];
rz(-0.71450662) q[1];
sx q[1];
rz(-2.6564068) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3109452) q[3];
sx q[3];
rz(-2.0214571) q[3];
sx q[3];
rz(-0.78307952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4063065) q[2];
sx q[2];
rz(-0.51484171) q[2];
sx q[2];
rz(1.7740645) q[2];
rz(-0.96274084) q[3];
sx q[3];
rz(-1.5509501) q[3];
sx q[3];
rz(2.4875557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3108567) q[0];
sx q[0];
rz(-1.5565358) q[0];
sx q[0];
rz(0.32064015) q[0];
rz(-2.2036208) q[1];
sx q[1];
rz(-1.2057949) q[1];
sx q[1];
rz(1.7373614) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7119634) q[0];
sx q[0];
rz(-0.14664635) q[0];
sx q[0];
rz(0.78369998) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.088555468) q[2];
sx q[2];
rz(-2.8635396) q[2];
sx q[2];
rz(-2.0012358) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.93130805) q[1];
sx q[1];
rz(-3.0979994) q[1];
sx q[1];
rz(1.0924073) q[1];
rz(-pi) q[2];
x q[2];
rz(0.25038428) q[3];
sx q[3];
rz(-0.77810771) q[3];
sx q[3];
rz(0.66533662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2001026) q[2];
sx q[2];
rz(-0.76277554) q[2];
sx q[2];
rz(0.74083677) q[2];
rz(-2.0790993) q[3];
sx q[3];
rz(-2.0334358) q[3];
sx q[3];
rz(1.197193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.9366539) q[0];
sx q[0];
rz(-2.7118201) q[0];
sx q[0];
rz(2.3858261) q[0];
rz(1.3142122) q[1];
sx q[1];
rz(-1.5145489) q[1];
sx q[1];
rz(-2.4155713) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8963658) q[0];
sx q[0];
rz(-0.78316044) q[0];
sx q[0];
rz(-3.0928715) q[0];
x q[1];
rz(-2.2560202) q[2];
sx q[2];
rz(-1.4902334) q[2];
sx q[2];
rz(-2.0256113) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0418813) q[1];
sx q[1];
rz(-1.9181632) q[1];
sx q[1];
rz(-2.9798286) q[1];
rz(-pi) q[2];
rz(1.8733379) q[3];
sx q[3];
rz(-2.6740251) q[3];
sx q[3];
rz(-2.9085161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3296335) q[2];
sx q[2];
rz(-2.6128431) q[2];
sx q[2];
rz(-0.94919666) q[2];
rz(-0.1720998) q[3];
sx q[3];
rz(-2.1642919) q[3];
sx q[3];
rz(-2.7026091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(3.1127472) q[0];
sx q[0];
rz(-1.5486568) q[0];
sx q[0];
rz(-1.5966709) q[0];
rz(-1.159809) q[1];
sx q[1];
rz(-1.3216959) q[1];
sx q[1];
rz(-1.6285223) q[1];
rz(2.8904758) q[2];
sx q[2];
rz(-0.56240964) q[2];
sx q[2];
rz(3.0887847) q[2];
rz(-3.0712745) q[3];
sx q[3];
rz(-2.0587772) q[3];
sx q[3];
rz(-1.8785431) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
