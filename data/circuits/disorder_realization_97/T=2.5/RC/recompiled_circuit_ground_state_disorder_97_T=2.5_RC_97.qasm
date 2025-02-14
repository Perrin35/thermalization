OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.41807732) q[0];
sx q[0];
rz(-0.014208566) q[0];
sx q[0];
rz(-0.040810458) q[0];
rz(3.1006676) q[1];
sx q[1];
rz(-2.0857781) q[1];
sx q[1];
rz(-2.5850886) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2262629) q[0];
sx q[0];
rz(-1.2253083) q[0];
sx q[0];
rz(0.19476632) q[0];
rz(-pi) q[1];
x q[1];
rz(2.56968) q[2];
sx q[2];
rz(-1.3777857) q[2];
sx q[2];
rz(-0.62096769) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8260658) q[1];
sx q[1];
rz(-0.64441753) q[1];
sx q[1];
rz(0.63912005) q[1];
rz(-pi) q[2];
rz(1.3338148) q[3];
sx q[3];
rz(-1.4795884) q[3];
sx q[3];
rz(2.3438675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1540264) q[2];
sx q[2];
rz(-1.7511) q[2];
sx q[2];
rz(-1.7517368) q[2];
rz(-3.0912345) q[3];
sx q[3];
rz(-1.4202838) q[3];
sx q[3];
rz(-1.095298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9339226) q[0];
sx q[0];
rz(-1.4140797) q[0];
sx q[0];
rz(0.06047824) q[0];
rz(0.98840797) q[1];
sx q[1];
rz(-1.6840839) q[1];
sx q[1];
rz(-2.9284533) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52046257) q[0];
sx q[0];
rz(-0.91875568) q[0];
sx q[0];
rz(2.0921642) q[0];
rz(-pi) q[1];
rz(-1.408281) q[2];
sx q[2];
rz(-1.6846416) q[2];
sx q[2];
rz(1.3448037) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4276398) q[1];
sx q[1];
rz(-0.25498495) q[1];
sx q[1];
rz(-2.990635) q[1];
rz(2.0673236) q[3];
sx q[3];
rz(-1.4098216) q[3];
sx q[3];
rz(-2.9825008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9947784) q[2];
sx q[2];
rz(-1.5922567) q[2];
sx q[2];
rz(3.1352622) q[2];
rz(1.7019042) q[3];
sx q[3];
rz(-2.3531745) q[3];
sx q[3];
rz(2.7009713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45563662) q[0];
sx q[0];
rz(-0.39304471) q[0];
sx q[0];
rz(-1.8259557) q[0];
rz(1.3872967) q[1];
sx q[1];
rz(-0.60391128) q[1];
sx q[1];
rz(-3.0953298) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9994557) q[0];
sx q[0];
rz(-1.5408514) q[0];
sx q[0];
rz(-1.5718476) q[0];
x q[1];
rz(-1.1836902) q[2];
sx q[2];
rz(-2.0069242) q[2];
sx q[2];
rz(1.3189463) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.97291482) q[1];
sx q[1];
rz(-0.70446223) q[1];
sx q[1];
rz(2.5903914) q[1];
x q[2];
rz(-0.1667695) q[3];
sx q[3];
rz(-0.96264364) q[3];
sx q[3];
rz(-1.2603354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5837273) q[2];
sx q[2];
rz(-1.1745619) q[2];
sx q[2];
rz(-2.2306856) q[2];
rz(1.9117428) q[3];
sx q[3];
rz(-2.3807447) q[3];
sx q[3];
rz(2.7228444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8987027) q[0];
sx q[0];
rz(-1.3985343) q[0];
sx q[0];
rz(-0.42309764) q[0];
rz(-0.93370071) q[1];
sx q[1];
rz(-1.9405148) q[1];
sx q[1];
rz(-2.7074714) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1541464) q[0];
sx q[0];
rz(-0.29329663) q[0];
sx q[0];
rz(3.0213839) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66017326) q[2];
sx q[2];
rz(-1.6962564) q[2];
sx q[2];
rz(-1.2646641) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4695376) q[1];
sx q[1];
rz(-2.5651882) q[1];
sx q[1];
rz(0.23643929) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7755327) q[3];
sx q[3];
rz(-1.4582485) q[3];
sx q[3];
rz(1.1074668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9551129) q[2];
sx q[2];
rz(-3.0035512) q[2];
sx q[2];
rz(-1.0008) q[2];
rz(0.75438386) q[3];
sx q[3];
rz(-2.3013134) q[3];
sx q[3];
rz(1.3365356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56483018) q[0];
sx q[0];
rz(-2.5046528) q[0];
sx q[0];
rz(1.7621967) q[0];
rz(-2.5881536) q[1];
sx q[1];
rz(-1.484) q[1];
sx q[1];
rz(0.71054202) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5003649) q[0];
sx q[0];
rz(-2.3661748) q[0];
sx q[0];
rz(1.7983318) q[0];
rz(-pi) q[1];
rz(1.5132873) q[2];
sx q[2];
rz(-1.3293666) q[2];
sx q[2];
rz(-1.4032569) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.48747739) q[1];
sx q[1];
rz(-0.23329188) q[1];
sx q[1];
rz(-1.7196991) q[1];
x q[2];
rz(-0.54039012) q[3];
sx q[3];
rz(-1.7149263) q[3];
sx q[3];
rz(-1.5976417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9589299) q[2];
sx q[2];
rz(-1.6561597) q[2];
sx q[2];
rz(0.043969285) q[2];
rz(-0.26642695) q[3];
sx q[3];
rz(-0.1259585) q[3];
sx q[3];
rz(0.21336475) q[3];
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
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5037395) q[0];
sx q[0];
rz(-1.5218691) q[0];
sx q[0];
rz(0.36852401) q[0];
rz(0.35399327) q[1];
sx q[1];
rz(-2.0123672) q[1];
sx q[1];
rz(2.0781719) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89404384) q[0];
sx q[0];
rz(-3.0008033) q[0];
sx q[0];
rz(-1.6611499) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11234653) q[2];
sx q[2];
rz(-1.8469166) q[2];
sx q[2];
rz(2.744855) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8700347) q[1];
sx q[1];
rz(-1.9816625) q[1];
sx q[1];
rz(0.4510057) q[1];
rz(-0.20078378) q[3];
sx q[3];
rz(-2.3827083) q[3];
sx q[3];
rz(-0.52606612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4949558) q[2];
sx q[2];
rz(-1.8821303) q[2];
sx q[2];
rz(-1.4259526) q[2];
rz(-0.87092733) q[3];
sx q[3];
rz(-2.0723074) q[3];
sx q[3];
rz(2.9072185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17212412) q[0];
sx q[0];
rz(-0.5640465) q[0];
sx q[0];
rz(-3.1308351) q[0];
rz(2.5116008) q[1];
sx q[1];
rz(-1.0456345) q[1];
sx q[1];
rz(-0.90525544) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0765823) q[0];
sx q[0];
rz(-2.6040051) q[0];
sx q[0];
rz(-1.697886) q[0];
x q[1];
rz(-2.0834604) q[2];
sx q[2];
rz(-1.8638637) q[2];
sx q[2];
rz(-0.91171564) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.71996321) q[1];
sx q[1];
rz(-1.0698228) q[1];
sx q[1];
rz(1.73354) q[1];
rz(-pi) q[2];
rz(0.42845528) q[3];
sx q[3];
rz(-0.96142229) q[3];
sx q[3];
rz(-0.46655015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.68911997) q[2];
sx q[2];
rz(-1.8931188) q[2];
sx q[2];
rz(-0.11736891) q[2];
rz(-0.72702414) q[3];
sx q[3];
rz(-0.89594642) q[3];
sx q[3];
rz(-1.1400384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.82379502) q[0];
sx q[0];
rz(-1.057484) q[0];
sx q[0];
rz(2.8463038) q[0];
rz(1.3914039) q[1];
sx q[1];
rz(-2.4877805) q[1];
sx q[1];
rz(-0.48694912) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3759252) q[0];
sx q[0];
rz(-1.4698979) q[0];
sx q[0];
rz(-0.19815066) q[0];
x q[1];
rz(-1.0560965) q[2];
sx q[2];
rz(-1.2916512) q[2];
sx q[2];
rz(-2.5893958) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0436957) q[1];
sx q[1];
rz(-2.0123031) q[1];
sx q[1];
rz(2.312078) q[1];
x q[2];
rz(1.7510909) q[3];
sx q[3];
rz(-1.3077985) q[3];
sx q[3];
rz(0.20151787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.22750172) q[2];
sx q[2];
rz(-0.5946244) q[2];
sx q[2];
rz(-1.5957069) q[2];
rz(2.5325736) q[3];
sx q[3];
rz(-1.1142497) q[3];
sx q[3];
rz(2.473089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85536999) q[0];
sx q[0];
rz(-0.90317059) q[0];
sx q[0];
rz(1.632593) q[0];
rz(1.9980994) q[1];
sx q[1];
rz(-1.1610463) q[1];
sx q[1];
rz(-3.031409) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0116029) q[0];
sx q[0];
rz(-2.5429259) q[0];
sx q[0];
rz(2.9574139) q[0];
rz(0.25709935) q[2];
sx q[2];
rz(-2.0201224) q[2];
sx q[2];
rz(-1.7420235) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5852803) q[1];
sx q[1];
rz(-1.1729128) q[1];
sx q[1];
rz(2.3817987) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1391644) q[3];
sx q[3];
rz(-0.49730992) q[3];
sx q[3];
rz(2.8893472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1630359) q[2];
sx q[2];
rz(-2.0608993) q[2];
sx q[2];
rz(2.2829096) q[2];
rz(-1.6363755) q[3];
sx q[3];
rz(-1.7965763) q[3];
sx q[3];
rz(-1.4428562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34864146) q[0];
sx q[0];
rz(-2.0703147) q[0];
sx q[0];
rz(0.68565482) q[0];
rz(-2.3220956) q[1];
sx q[1];
rz(-1.6623431) q[1];
sx q[1];
rz(2.562838) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98812719) q[0];
sx q[0];
rz(-2.1442226) q[0];
sx q[0];
rz(0.98116409) q[0];
rz(2.1376325) q[2];
sx q[2];
rz(-1.7805779) q[2];
sx q[2];
rz(-0.45406859) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.600665) q[1];
sx q[1];
rz(-1.3131716) q[1];
sx q[1];
rz(-0.39389807) q[1];
rz(-pi) q[2];
rz(1.6195504) q[3];
sx q[3];
rz(-1.6798868) q[3];
sx q[3];
rz(-1.297612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3026352) q[2];
sx q[2];
rz(-0.18582782) q[2];
sx q[2];
rz(0.1798943) q[2];
rz(2.5707865) q[3];
sx q[3];
rz(-2.4195318) q[3];
sx q[3];
rz(2.4898237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30052341) q[0];
sx q[0];
rz(-0.64995926) q[0];
sx q[0];
rz(1.7267701) q[0];
rz(0.48814804) q[1];
sx q[1];
rz(-2.585325) q[1];
sx q[1];
rz(-1.5811031) q[1];
rz(-0.93964259) q[2];
sx q[2];
rz(-1.8710536) q[2];
sx q[2];
rz(-1.6761129) q[2];
rz(0.53849738) q[3];
sx q[3];
rz(-1.8168728) q[3];
sx q[3];
rz(-1.3088165) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
