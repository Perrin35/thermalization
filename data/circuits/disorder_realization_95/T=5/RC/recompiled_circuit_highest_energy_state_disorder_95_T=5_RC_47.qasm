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
rz(-1.7669825) q[0];
sx q[0];
rz(-0.97281015) q[0];
sx q[0];
rz(1.9109803) q[0];
rz(1.5988916) q[1];
sx q[1];
rz(-2.7886432) q[1];
sx q[1];
rz(-0.43500873) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9821138) q[0];
sx q[0];
rz(-1.4712988) q[0];
sx q[0];
rz(-1.2364819) q[0];
rz(-pi) q[1];
rz(1.4392008) q[2];
sx q[2];
rz(-1.7965266) q[2];
sx q[2];
rz(2.7011958) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5953345) q[1];
sx q[1];
rz(-1.8500684) q[1];
sx q[1];
rz(-0.78650773) q[1];
rz(-1.6023956) q[3];
sx q[3];
rz(-1.5571619) q[3];
sx q[3];
rz(-2.0054833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1970485) q[2];
sx q[2];
rz(-0.41434449) q[2];
sx q[2];
rz(-2.22331) q[2];
rz(-0.36327547) q[3];
sx q[3];
rz(-1.8922292) q[3];
sx q[3];
rz(-1.3242807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0927703) q[0];
sx q[0];
rz(-0.73960441) q[0];
sx q[0];
rz(-2.0059465) q[0];
rz(0.28226918) q[1];
sx q[1];
rz(-1.5156563) q[1];
sx q[1];
rz(-2.938882) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1644703) q[0];
sx q[0];
rz(-2.8852069) q[0];
sx q[0];
rz(0.77986832) q[0];
rz(-pi) q[1];
x q[1];
rz(1.545867) q[2];
sx q[2];
rz(-2.6597616) q[2];
sx q[2];
rz(-2.3584751) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9276638) q[1];
sx q[1];
rz(-2.0239554) q[1];
sx q[1];
rz(-1.434668) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.0052350076) q[3];
sx q[3];
rz(-1.5810895) q[3];
sx q[3];
rz(2.0106842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.89381605) q[2];
sx q[2];
rz(-2.717369) q[2];
sx q[2];
rz(-0.068537863) q[2];
rz(2.2872772) q[3];
sx q[3];
rz(-1.2645384) q[3];
sx q[3];
rz(-0.0059277047) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0900367) q[0];
sx q[0];
rz(-2.5904901) q[0];
sx q[0];
rz(2.0447482) q[0];
rz(-0.54496533) q[1];
sx q[1];
rz(-2.4577591) q[1];
sx q[1];
rz(-0.1942689) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.504963) q[0];
sx q[0];
rz(-1.2218814) q[0];
sx q[0];
rz(1.1303085) q[0];
x q[1];
rz(1.0610007) q[2];
sx q[2];
rz(-2.5246127) q[2];
sx q[2];
rz(1.4088907) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4020445) q[1];
sx q[1];
rz(-0.6303936) q[1];
sx q[1];
rz(2.5044548) q[1];
rz(-pi) q[2];
rz(1.9241289) q[3];
sx q[3];
rz(-1.4446745) q[3];
sx q[3];
rz(-2.452891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.37476173) q[2];
sx q[2];
rz(-1.7685879) q[2];
sx q[2];
rz(-2.4135446) q[2];
rz(-0.22045615) q[3];
sx q[3];
rz(-2.9900592) q[3];
sx q[3];
rz(-2.5098586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24882889) q[0];
sx q[0];
rz(-1.4735824) q[0];
sx q[0];
rz(-2.072075) q[0];
rz(-0.88691521) q[1];
sx q[1];
rz(-0.506217) q[1];
sx q[1];
rz(1.3217529) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89267545) q[0];
sx q[0];
rz(-1.6098611) q[0];
sx q[0];
rz(-0.43401735) q[0];
rz(0.063912674) q[2];
sx q[2];
rz(-2.9021429) q[2];
sx q[2];
rz(-2.4452345) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.48044887) q[1];
sx q[1];
rz(-2.6972507) q[1];
sx q[1];
rz(-1.737887) q[1];
rz(-2.6245313) q[3];
sx q[3];
rz(-1.0365465) q[3];
sx q[3];
rz(1.8815966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.47361031) q[2];
sx q[2];
rz(-2.0986291) q[2];
sx q[2];
rz(1.854151) q[2];
rz(2.1227664) q[3];
sx q[3];
rz(-2.5835218) q[3];
sx q[3];
rz(2.4655925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2126615) q[0];
sx q[0];
rz(-0.23290817) q[0];
sx q[0];
rz(2.2610597) q[0];
rz(3.0048043) q[1];
sx q[1];
rz(-0.22577481) q[1];
sx q[1];
rz(0.61019623) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5870415) q[0];
sx q[0];
rz(-2.4992538) q[0];
sx q[0];
rz(1.5498398) q[0];
rz(0.32322804) q[2];
sx q[2];
rz(-0.45219496) q[2];
sx q[2];
rz(-1.5791073) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5049664) q[1];
sx q[1];
rz(-1.7646878) q[1];
sx q[1];
rz(-2.8411179) q[1];
rz(-pi) q[2];
rz(2.8585096) q[3];
sx q[3];
rz(-0.45520458) q[3];
sx q[3];
rz(-2.9415123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.60787624) q[2];
sx q[2];
rz(-2.4166962) q[2];
sx q[2];
rz(1.2302715) q[2];
rz(2.2599334) q[3];
sx q[3];
rz(-1.6917546) q[3];
sx q[3];
rz(-1.3151883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6864258) q[0];
sx q[0];
rz(-1.5436341) q[0];
sx q[0];
rz(-1.0536739) q[0];
rz(1.3764489) q[1];
sx q[1];
rz(-1.1076628) q[1];
sx q[1];
rz(2.460316) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2864425) q[0];
sx q[0];
rz(-1.5871829) q[0];
sx q[0];
rz(0.70006242) q[0];
rz(-0.7448911) q[2];
sx q[2];
rz(-1.3110263) q[2];
sx q[2];
rz(-1.3448754) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.14398814) q[1];
sx q[1];
rz(-1.0790842) q[1];
sx q[1];
rz(-2.8365447) q[1];
rz(-2.2627955) q[3];
sx q[3];
rz(-1.5419772) q[3];
sx q[3];
rz(-1.6520713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.072558746) q[2];
sx q[2];
rz(-2.7636038) q[2];
sx q[2];
rz(0.47522137) q[2];
rz(-1.2230988) q[3];
sx q[3];
rz(-2.1971072) q[3];
sx q[3];
rz(-0.16330115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-1.908919) q[0];
sx q[0];
rz(-0.95218807) q[0];
sx q[0];
rz(-2.9222144) q[0];
rz(1.0064005) q[1];
sx q[1];
rz(-1.718113) q[1];
sx q[1];
rz(-1.2054319) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69898283) q[0];
sx q[0];
rz(-0.97215334) q[0];
sx q[0];
rz(3.1151616) q[0];
rz(-pi) q[1];
rz(2.7458887) q[2];
sx q[2];
rz(-1.323441) q[2];
sx q[2];
rz(-2.5620714) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0675065) q[1];
sx q[1];
rz(-1.5879022) q[1];
sx q[1];
rz(-2.9422804) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1262283) q[3];
sx q[3];
rz(-1.6939112) q[3];
sx q[3];
rz(0.11492226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9203303) q[2];
sx q[2];
rz(-2.7483676) q[2];
sx q[2];
rz(0.66366759) q[2];
rz(-2.4237733) q[3];
sx q[3];
rz(-0.65687537) q[3];
sx q[3];
rz(-0.17620152) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8127036) q[0];
sx q[0];
rz(-2.1825574) q[0];
sx q[0];
rz(1.2411728) q[0];
rz(-2.7136956) q[1];
sx q[1];
rz(-1.5792081) q[1];
sx q[1];
rz(1.7180299) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5964029) q[0];
sx q[0];
rz(-1.3529142) q[0];
sx q[0];
rz(-1.4942188) q[0];
rz(-pi) q[1];
rz(-0.85278089) q[2];
sx q[2];
rz(-1.4508392) q[2];
sx q[2];
rz(0.44242417) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9729303) q[1];
sx q[1];
rz(-0.902938) q[1];
sx q[1];
rz(0.027110312) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0377722) q[3];
sx q[3];
rz(-0.33664931) q[3];
sx q[3];
rz(2.3262084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1443892) q[2];
sx q[2];
rz(-2.453697) q[2];
sx q[2];
rz(0.29590657) q[2];
rz(1.0812673) q[3];
sx q[3];
rz(-1.5238949) q[3];
sx q[3];
rz(-2.9937939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5645771) q[0];
sx q[0];
rz(-0.58203375) q[0];
sx q[0];
rz(-0.70593315) q[0];
rz(0.18711807) q[1];
sx q[1];
rz(-2.3046604) q[1];
sx q[1];
rz(-0.44609889) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68413823) q[0];
sx q[0];
rz(-1.7230526) q[0];
sx q[0];
rz(1.5726552) q[0];
rz(2.1438333) q[2];
sx q[2];
rz(-0.54513844) q[2];
sx q[2];
rz(1.0447431) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2526568) q[1];
sx q[1];
rz(-1.2351079) q[1];
sx q[1];
rz(-0.97519213) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8817552) q[3];
sx q[3];
rz(-2.5492955) q[3];
sx q[3];
rz(-0.01334503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.17911653) q[2];
sx q[2];
rz(-1.0255145) q[2];
sx q[2];
rz(1.4013438) q[2];
rz(-2.5380747) q[3];
sx q[3];
rz(-2.9602435) q[3];
sx q[3];
rz(0.13874273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.03183455) q[0];
sx q[0];
rz(-1.7367481) q[0];
sx q[0];
rz(3.0314714) q[0];
rz(1.3212063) q[1];
sx q[1];
rz(-0.92480129) q[1];
sx q[1];
rz(-0.22470156) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6075144) q[0];
sx q[0];
rz(-0.67994962) q[0];
sx q[0];
rz(2.2282766) q[0];
rz(-2.0761833) q[2];
sx q[2];
rz(-0.32534625) q[2];
sx q[2];
rz(2.4869652) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4555295) q[1];
sx q[1];
rz(-1.7357329) q[1];
sx q[1];
rz(2.4071724) q[1];
rz(-1.3984001) q[3];
sx q[3];
rz(-1.5443003) q[3];
sx q[3];
rz(-0.66988984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8297742) q[2];
sx q[2];
rz(-2.8303787) q[2];
sx q[2];
rz(-2.6149926) q[2];
rz(0.38463587) q[3];
sx q[3];
rz(-1.8943818) q[3];
sx q[3];
rz(-0.23469901) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2100691) q[0];
sx q[0];
rz(-2.7061404) q[0];
sx q[0];
rz(2.7192116) q[0];
rz(-0.79413636) q[1];
sx q[1];
rz(-1.7105449) q[1];
sx q[1];
rz(1.7078043) q[1];
rz(-0.88243816) q[2];
sx q[2];
rz(-2.0954162) q[2];
sx q[2];
rz(1.4388234) q[2];
rz(-0.41153367) q[3];
sx q[3];
rz(-1.2306662) q[3];
sx q[3];
rz(-1.6231404) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
