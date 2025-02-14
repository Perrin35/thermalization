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
rz(-2.0432552) q[0];
rz(1.1605473) q[1];
sx q[1];
rz(2.0055973) q[1];
sx q[1];
rz(10.027167) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1404323) q[0];
sx q[0];
rz(-2.4386157) q[0];
sx q[0];
rz(-0.70859895) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4279934) q[2];
sx q[2];
rz(-1.6125792) q[2];
sx q[2];
rz(0.080619911) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.25665313) q[1];
sx q[1];
rz(-0.49044427) q[1];
sx q[1];
rz(1.8153193) q[1];
x q[2];
rz(-1.9992152) q[3];
sx q[3];
rz(-1.6049634) q[3];
sx q[3];
rz(2.7673801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7263553) q[2];
sx q[2];
rz(-1.0865612) q[2];
sx q[2];
rz(-1.497867) q[2];
rz(-0.1114397) q[3];
sx q[3];
rz(-1.8569511) q[3];
sx q[3];
rz(2.1845412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1794671) q[0];
sx q[0];
rz(-0.92573708) q[0];
sx q[0];
rz(-0.16394462) q[0];
rz(0.97490772) q[1];
sx q[1];
rz(-2.4539852) q[1];
sx q[1];
rz(-1.6418537) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83368084) q[0];
sx q[0];
rz(-1.320086) q[0];
sx q[0];
rz(-2.3019019) q[0];
rz(-2.4663839) q[2];
sx q[2];
rz(-2.2134502) q[2];
sx q[2];
rz(-0.86057885) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.837599) q[1];
sx q[1];
rz(-1.493206) q[1];
sx q[1];
rz(-0.71240387) q[1];
x q[2];
rz(0.99489958) q[3];
sx q[3];
rz(-0.38299503) q[3];
sx q[3];
rz(0.3037616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6770596) q[2];
sx q[2];
rz(-2.634282) q[2];
sx q[2];
rz(2.9110009) q[2];
rz(-2.4584037) q[3];
sx q[3];
rz(-2.4468827) q[3];
sx q[3];
rz(2.3901239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7815919) q[0];
sx q[0];
rz(-1.3104023) q[0];
sx q[0];
rz(-2.7202284) q[0];
rz(1.5158117) q[1];
sx q[1];
rz(-0.49585626) q[1];
sx q[1];
rz(0.96744195) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4801299) q[0];
sx q[0];
rz(-0.41981831) q[0];
sx q[0];
rz(0.90645193) q[0];
x q[1];
rz(0.17379667) q[2];
sx q[2];
rz(-2.0325543) q[2];
sx q[2];
rz(2.6231036) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.53114226) q[1];
sx q[1];
rz(-1.4262661) q[1];
sx q[1];
rz(1.2953613) q[1];
x q[2];
rz(-1.9032962) q[3];
sx q[3];
rz(-2.1070814) q[3];
sx q[3];
rz(-3.0474049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4030054) q[2];
sx q[2];
rz(-2.1295197) q[2];
sx q[2];
rz(-2.9834413) q[2];
rz(-2.3824597) q[3];
sx q[3];
rz(-1.682155) q[3];
sx q[3];
rz(2.7706743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.89449969) q[0];
sx q[0];
rz(-0.58421725) q[0];
sx q[0];
rz(1.9547113) q[0];
rz(3.0184556) q[1];
sx q[1];
rz(-1.5712527) q[1];
sx q[1];
rz(-1.311696) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9378915) q[0];
sx q[0];
rz(-0.16898705) q[0];
sx q[0];
rz(1.0615055) q[0];
rz(1.4940631) q[2];
sx q[2];
rz(-1.4686246) q[2];
sx q[2];
rz(0.42021423) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.32320106) q[1];
sx q[1];
rz(-1.15266) q[1];
sx q[1];
rz(1.1365324) q[1];
rz(-1.5748686) q[3];
sx q[3];
rz(-1.5713672) q[3];
sx q[3];
rz(-2.3426995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.19795236) q[2];
sx q[2];
rz(-1.8215048) q[2];
sx q[2];
rz(-0.8521592) q[2];
rz(-2.2876244) q[3];
sx q[3];
rz(-1.2205418) q[3];
sx q[3];
rz(-2.061969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3581486) q[0];
sx q[0];
rz(-1.3842979) q[0];
sx q[0];
rz(3.070991) q[0];
rz(-2.103503) q[1];
sx q[1];
rz(-1.9990653) q[1];
sx q[1];
rz(-0.80931726) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22915086) q[0];
sx q[0];
rz(-0.99780592) q[0];
sx q[0];
rz(-1.4848764) q[0];
rz(2.3951934) q[2];
sx q[2];
rz(-2.0944907) q[2];
sx q[2];
rz(0.46800287) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.091920225) q[1];
sx q[1];
rz(-0.61245239) q[1];
sx q[1];
rz(2.6785707) q[1];
x q[2];
rz(-2.9633459) q[3];
sx q[3];
rz(-2.2615454) q[3];
sx q[3];
rz(-0.33876283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.91850963) q[2];
sx q[2];
rz(-0.56879908) q[2];
sx q[2];
rz(-1.1893547) q[2];
rz(-3.0009624) q[3];
sx q[3];
rz(-0.46969241) q[3];
sx q[3];
rz(-0.3895337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73914948) q[0];
sx q[0];
rz(-0.91366714) q[0];
sx q[0];
rz(1.5341349) q[0];
rz(2.2588579) q[1];
sx q[1];
rz(-2.2969756) q[1];
sx q[1];
rz(-0.61558634) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7587335) q[0];
sx q[0];
rz(-0.40399562) q[0];
sx q[0];
rz(2.9055779) q[0];
rz(0.87178715) q[2];
sx q[2];
rz(-1.0411556) q[2];
sx q[2];
rz(2.3123534) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.57762417) q[1];
sx q[1];
rz(-1.3805439) q[1];
sx q[1];
rz(-0.14037324) q[1];
rz(1.5869523) q[3];
sx q[3];
rz(-0.62045853) q[3];
sx q[3];
rz(-2.4718493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9756056) q[2];
sx q[2];
rz(-0.96419445) q[2];
sx q[2];
rz(-1.1791505) q[2];
rz(-2.9278582) q[3];
sx q[3];
rz(-1.390018) q[3];
sx q[3];
rz(-2.8960622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24446503) q[0];
sx q[0];
rz(-1.1470969) q[0];
sx q[0];
rz(-3.097528) q[0];
rz(-0.38912082) q[1];
sx q[1];
rz(-0.48946425) q[1];
sx q[1];
rz(-2.3997831) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7334866) q[0];
sx q[0];
rz(-1.2950962) q[0];
sx q[0];
rz(2.7635126) q[0];
x q[1];
rz(-1.1941694) q[2];
sx q[2];
rz(-2.3583226) q[2];
sx q[2];
rz(2.7830091) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9981813) q[1];
sx q[1];
rz(-0.3340958) q[1];
sx q[1];
rz(0.27004529) q[1];
rz(2.4398285) q[3];
sx q[3];
rz(-1.7265111) q[3];
sx q[3];
rz(-3.0879741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47702181) q[2];
sx q[2];
rz(-1.3571813) q[2];
sx q[2];
rz(2.6396497) q[2];
rz(-1.7160412) q[3];
sx q[3];
rz(-0.52949667) q[3];
sx q[3];
rz(2.3486923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3821756) q[0];
sx q[0];
rz(-1.734153) q[0];
sx q[0];
rz(-0.39696804) q[0];
rz(-1.601903) q[1];
sx q[1];
rz(-2.868728) q[1];
sx q[1];
rz(-2.4698965) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8803589) q[0];
sx q[0];
rz(-1.6971998) q[0];
sx q[0];
rz(2.0194299) q[0];
rz(-pi) q[1];
rz(-1.1279068) q[2];
sx q[2];
rz(-1.2338362) q[2];
sx q[2];
rz(0.86293399) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4038329) q[1];
sx q[1];
rz(-1.8813526) q[1];
sx q[1];
rz(-0.65447372) q[1];
x q[2];
rz(0.94843169) q[3];
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
rz(0.73528618) q[2];
sx q[2];
rz(-2.6267509) q[2];
sx q[2];
rz(1.3675281) q[2];
rz(0.96274084) q[3];
sx q[3];
rz(-1.5906426) q[3];
sx q[3];
rz(2.4875557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.830736) q[0];
sx q[0];
rz(-1.5850569) q[0];
sx q[0];
rz(-0.32064015) q[0];
rz(0.93797183) q[1];
sx q[1];
rz(-1.9357977) q[1];
sx q[1];
rz(-1.7373614) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9228678) q[0];
sx q[0];
rz(-1.6744807) q[0];
sx q[0];
rz(-1.4669048) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27701849) q[2];
sx q[2];
rz(-1.5465186) q[2];
sx q[2];
rz(-2.6259823) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.45253048) q[1];
sx q[1];
rz(-1.5320996) q[1];
sx q[1];
rz(-0.020078151) q[1];
rz(2.8912084) q[3];
sx q[3];
rz(-2.3634849) q[3];
sx q[3];
rz(-2.476256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2001026) q[2];
sx q[2];
rz(-2.3788171) q[2];
sx q[2];
rz(-0.74083677) q[2];
rz(-2.0790993) q[3];
sx q[3];
rz(-1.1081568) q[3];
sx q[3];
rz(-1.197193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9366539) q[0];
sx q[0];
rz(-2.7118201) q[0];
sx q[0];
rz(-2.3858261) q[0];
rz(1.3142122) q[1];
sx q[1];
rz(-1.6270437) q[1];
sx q[1];
rz(-0.72602138) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8963658) q[0];
sx q[0];
rz(-2.3584322) q[0];
sx q[0];
rz(-0.048721192) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2560202) q[2];
sx q[2];
rz(-1.6513593) q[2];
sx q[2];
rz(-1.1159814) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5948702) q[1];
sx q[1];
rz(-2.7597962) q[1];
sx q[1];
rz(1.1522271) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8733379) q[3];
sx q[3];
rz(-2.6740251) q[3];
sx q[3];
rz(0.23307652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3296335) q[2];
sx q[2];
rz(-2.6128431) q[2];
sx q[2];
rz(2.192396) q[2];
rz(-0.1720998) q[3];
sx q[3];
rz(-0.97730079) q[3];
sx q[3];
rz(2.7026091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1127472) q[0];
sx q[0];
rz(-1.5929359) q[0];
sx q[0];
rz(1.5449217) q[0];
rz(-1.9817837) q[1];
sx q[1];
rz(-1.8198967) q[1];
sx q[1];
rz(1.5130704) q[1];
rz(2.5934577) q[2];
sx q[2];
rz(-1.4379063) q[2];
sx q[2];
rz(-1.4098991) q[2];
rz(-1.439194) q[3];
sx q[3];
rz(-0.49261649) q[3];
sx q[3];
rz(-1.7294283) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
