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
rz(-1.1385318) q[0];
sx q[0];
rz(1.0983374) q[0];
rz(1.1605473) q[1];
sx q[1];
rz(-1.1359954) q[1];
sx q[1];
rz(2.5392037) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9928474) q[0];
sx q[0];
rz(-2.0050328) q[0];
sx q[0];
rz(-2.56987) q[0];
rz(-1.5155529) q[2];
sx q[2];
rz(-2.2836402) q[2];
sx q[2];
rz(-1.6875658) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8849395) q[1];
sx q[1];
rz(-2.6511484) q[1];
sx q[1];
rz(-1.3262733) q[1];
x q[2];
rz(-1.9992152) q[3];
sx q[3];
rz(-1.6049634) q[3];
sx q[3];
rz(-0.37421253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.41523734) q[2];
sx q[2];
rz(-1.0865612) q[2];
sx q[2];
rz(1.497867) q[2];
rz(3.030153) q[3];
sx q[3];
rz(-1.8569511) q[3];
sx q[3];
rz(2.1845412) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1794671) q[0];
sx q[0];
rz(-0.92573708) q[0];
sx q[0];
rz(0.16394462) q[0];
rz(-2.1666849) q[1];
sx q[1];
rz(-2.4539852) q[1];
sx q[1];
rz(1.499739) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3079118) q[0];
sx q[0];
rz(-1.8215067) q[0];
sx q[0];
rz(-2.3019019) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80625466) q[2];
sx q[2];
rz(-1.0467751) q[2];
sx q[2];
rz(-1.1576705) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1999535) q[1];
sx q[1];
rz(-2.2806045) q[1];
sx q[1];
rz(1.468424) q[1];
x q[2];
rz(0.21598609) q[3];
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
rz(0.683189) q[3];
sx q[3];
rz(-2.4468827) q[3];
sx q[3];
rz(2.3901239) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7815919) q[0];
sx q[0];
rz(-1.8311904) q[0];
sx q[0];
rz(-0.42136425) q[0];
rz(-1.5158117) q[1];
sx q[1];
rz(-0.49585626) q[1];
sx q[1];
rz(-0.96744195) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.370328) q[0];
sx q[0];
rz(-1.8974842) q[0];
sx q[0];
rz(2.8730434) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2363628) q[2];
sx q[2];
rz(-2.6504271) q[2];
sx q[2];
rz(-2.9984891) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6306259) q[1];
sx q[1];
rz(-0.31019638) q[1];
sx q[1];
rz(-2.0621745) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2382965) q[3];
sx q[3];
rz(-2.1070814) q[3];
sx q[3];
rz(3.0474049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4030054) q[2];
sx q[2];
rz(-2.1295197) q[2];
sx q[2];
rz(-2.9834413) q[2];
rz(-0.75913298) q[3];
sx q[3];
rz(-1.4594376) q[3];
sx q[3];
rz(2.7706743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.247093) q[0];
sx q[0];
rz(-0.58421725) q[0];
sx q[0];
rz(1.9547113) q[0];
rz(-0.12313708) q[1];
sx q[1];
rz(-1.5703399) q[1];
sx q[1];
rz(-1.8298967) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7191221) q[0];
sx q[0];
rz(-1.718169) q[0];
sx q[0];
rz(-3.0585994) q[0];
x q[1];
rz(-1.6475296) q[2];
sx q[2];
rz(-1.6729681) q[2];
sx q[2];
rz(-0.42021423) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8183916) q[1];
sx q[1];
rz(-1.9889327) q[1];
sx q[1];
rz(-2.0050603) q[1];
rz(-pi) q[2];
rz(-1.4315286) q[3];
sx q[3];
rz(-0.0041120681) q[3];
sx q[3];
rz(-0.63263661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9436403) q[2];
sx q[2];
rz(-1.8215048) q[2];
sx q[2];
rz(-2.2894335) q[2];
rz(-2.2876244) q[3];
sx q[3];
rz(-1.9210509) q[3];
sx q[3];
rz(2.061969) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78344408) q[0];
sx q[0];
rz(-1.3842979) q[0];
sx q[0];
rz(-0.070601687) q[0];
rz(2.103503) q[1];
sx q[1];
rz(-1.9990653) q[1];
sx q[1];
rz(-2.3322754) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7532849) q[0];
sx q[0];
rz(-1.6429672) q[0];
sx q[0];
rz(-0.57467527) q[0];
rz(-pi) q[1];
rz(2.2373166) q[2];
sx q[2];
rz(-2.1993786) q[2];
sx q[2];
rz(1.6055799) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0496724) q[1];
sx q[1];
rz(-0.61245239) q[1];
sx q[1];
rz(-2.6785707) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3594987) q[3];
sx q[3];
rz(-2.431892) q[3];
sx q[3];
rz(3.078408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.91850963) q[2];
sx q[2];
rz(-0.56879908) q[2];
sx q[2];
rz(-1.1893547) q[2];
rz(-3.0009624) q[3];
sx q[3];
rz(-0.46969241) q[3];
sx q[3];
rz(2.752059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73914948) q[0];
sx q[0];
rz(-0.91366714) q[0];
sx q[0];
rz(1.5341349) q[0];
rz(-2.2588579) q[1];
sx q[1];
rz(-0.84461707) q[1];
sx q[1];
rz(-0.61558634) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12703757) q[0];
sx q[0];
rz(-1.1786228) q[0];
sx q[0];
rz(-1.4711625) q[0];
rz(-pi) q[1];
rz(2.2698055) q[2];
sx q[2];
rz(-2.1004371) q[2];
sx q[2];
rz(-0.82923928) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.96645725) q[1];
sx q[1];
rz(-1.7086204) q[1];
sx q[1];
rz(-1.7628927) q[1];
rz(-pi) q[2];
rz(-0.011544051) q[3];
sx q[3];
rz(-0.95043105) q[3];
sx q[3];
rz(0.64988713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9756056) q[2];
sx q[2];
rz(-2.1773982) q[2];
sx q[2];
rz(1.1791505) q[2];
rz(0.21373448) q[3];
sx q[3];
rz(-1.390018) q[3];
sx q[3];
rz(-2.8960622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8971276) q[0];
sx q[0];
rz(-1.9944958) q[0];
sx q[0];
rz(-0.044064673) q[0];
rz(-0.38912082) q[1];
sx q[1];
rz(-0.48946425) q[1];
sx q[1];
rz(0.74180952) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.378063) q[0];
sx q[0];
rz(-0.46398315) q[0];
sx q[0];
rz(0.65391175) q[0];
x q[1];
rz(1.1941694) q[2];
sx q[2];
rz(-0.78327006) q[2];
sx q[2];
rz(2.7830091) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.85838137) q[1];
sx q[1];
rz(-1.2492531) q[1];
sx q[1];
rz(1.47846) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70176418) q[3];
sx q[3];
rz(-1.7265111) q[3];
sx q[3];
rz(-0.053618535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.47702181) q[2];
sx q[2];
rz(-1.3571813) q[2];
sx q[2];
rz(-0.50194293) q[2];
rz(-1.4255514) q[3];
sx q[3];
rz(-0.52949667) q[3];
sx q[3];
rz(-2.3486923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75941706) q[0];
sx q[0];
rz(-1.4074396) q[0];
sx q[0];
rz(2.7446246) q[0];
rz(1.5396897) q[1];
sx q[1];
rz(-2.868728) q[1];
sx q[1];
rz(-2.4698965) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8803589) q[0];
sx q[0];
rz(-1.6971998) q[0];
sx q[0];
rz(-1.1221627) q[0];
rz(-pi) q[1];
rz(-2.7717088) q[2];
sx q[2];
rz(-1.9871759) q[2];
sx q[2];
rz(2.2781792) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.063350141) q[1];
sx q[1];
rz(-0.95253175) q[1];
sx q[1];
rz(1.1863913) q[1];
x q[2];
rz(-2.193161) q[3];
sx q[3];
rz(-0.84377224) q[3];
sx q[3];
rz(-1.2326944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4063065) q[2];
sx q[2];
rz(-0.51484171) q[2];
sx q[2];
rz(-1.3675281) q[2];
rz(-2.1788518) q[3];
sx q[3];
rz(-1.5509501) q[3];
sx q[3];
rz(0.65403691) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3108567) q[0];
sx q[0];
rz(-1.5565358) q[0];
sx q[0];
rz(-0.32064015) q[0];
rz(-0.93797183) q[1];
sx q[1];
rz(-1.9357977) q[1];
sx q[1];
rz(1.7373614) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42962924) q[0];
sx q[0];
rz(-2.9949463) q[0];
sx q[0];
rz(-2.3578927) q[0];
x q[1];
rz(-1.5455568) q[2];
sx q[2];
rz(-1.2938616) q[2];
sx q[2];
rz(-1.0482839) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2102846) q[1];
sx q[1];
rz(-0.043593229) q[1];
sx q[1];
rz(1.0924073) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8102989) q[3];
sx q[3];
rz(-2.31862) q[3];
sx q[3];
rz(2.8209958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9414901) q[2];
sx q[2];
rz(-0.76277554) q[2];
sx q[2];
rz(-2.4007559) q[2];
rz(2.0790993) q[3];
sx q[3];
rz(-2.0334358) q[3];
sx q[3];
rz(1.9443996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9366539) q[0];
sx q[0];
rz(-0.42977253) q[0];
sx q[0];
rz(-2.3858261) q[0];
rz(1.3142122) q[1];
sx q[1];
rz(-1.5145489) q[1];
sx q[1];
rz(0.72602138) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8963658) q[0];
sx q[0];
rz(-2.3584322) q[0];
sx q[0];
rz(0.048721192) q[0];
x q[1];
rz(0.88557247) q[2];
sx q[2];
rz(-1.6513593) q[2];
sx q[2];
rz(2.0256113) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0418813) q[1];
sx q[1];
rz(-1.9181632) q[1];
sx q[1];
rz(-0.16176407) q[1];
x q[2];
rz(-1.1216702) q[3];
sx q[3];
rz(-1.4360997) q[3];
sx q[3];
rz(-2.0756034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8119592) q[2];
sx q[2];
rz(-2.6128431) q[2];
sx q[2];
rz(0.94919666) q[2];
rz(0.1720998) q[3];
sx q[3];
rz(-2.1642919) q[3];
sx q[3];
rz(-0.43898359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1127472) q[0];
sx q[0];
rz(-1.5486568) q[0];
sx q[0];
rz(-1.5966709) q[0];
rz(-1.9817837) q[1];
sx q[1];
rz(-1.8198967) q[1];
sx q[1];
rz(1.5130704) q[1];
rz(0.25111689) q[2];
sx q[2];
rz(-2.579183) q[2];
sx q[2];
rz(-0.05280799) q[2];
rz(1.439194) q[3];
sx q[3];
rz(-2.6489762) q[3];
sx q[3];
rz(1.4121644) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
