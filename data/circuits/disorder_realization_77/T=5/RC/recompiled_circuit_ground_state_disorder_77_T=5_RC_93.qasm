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
rz(-1.1359954) q[1];
sx q[1];
rz(2.5392037) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9928474) q[0];
sx q[0];
rz(-1.1365599) q[0];
sx q[0];
rz(2.56987) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7135993) q[2];
sx q[2];
rz(-1.5290135) q[2];
sx q[2];
rz(0.080619911) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5307926) q[1];
sx q[1];
rz(-1.4565174) q[1];
sx q[1];
rz(-2.0487993) q[1];
rz(-pi) q[2];
rz(-1.9992152) q[3];
sx q[3];
rz(-1.5366293) q[3];
sx q[3];
rz(-2.7673801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.41523734) q[2];
sx q[2];
rz(-1.0865612) q[2];
sx q[2];
rz(1.6437257) q[2];
rz(-3.030153) q[3];
sx q[3];
rz(-1.2846416) q[3];
sx q[3];
rz(2.1845412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-1.1794671) q[0];
sx q[0];
rz(-0.92573708) q[0];
sx q[0];
rz(2.977648) q[0];
rz(-2.1666849) q[1];
sx q[1];
rz(-2.4539852) q[1];
sx q[1];
rz(-1.6418537) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83368084) q[0];
sx q[0];
rz(-1.8215067) q[0];
sx q[0];
rz(2.3019019) q[0];
rz(-pi) q[1];
rz(0.87514211) q[2];
sx q[2];
rz(-2.2457457) q[2];
sx q[2];
rz(0.067719134) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7852868) q[1];
sx q[1];
rz(-0.71588022) q[1];
sx q[1];
rz(-0.11838491) q[1];
x q[2];
rz(-1.2449366) q[3];
sx q[3];
rz(-1.3658524) q[3];
sx q[3];
rz(0.72494331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.46453309) q[2];
sx q[2];
rz(-0.50731069) q[2];
sx q[2];
rz(2.9110009) q[2];
rz(-0.683189) q[3];
sx q[3];
rz(-0.69470996) q[3];
sx q[3];
rz(2.3901239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36000073) q[0];
sx q[0];
rz(-1.8311904) q[0];
sx q[0];
rz(-2.7202284) q[0];
rz(-1.5158117) q[1];
sx q[1];
rz(-2.6457364) q[1];
sx q[1];
rz(0.96744195) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4801299) q[0];
sx q[0];
rz(-2.7217743) q[0];
sx q[0];
rz(2.2351407) q[0];
rz(-pi) q[1];
rz(1.9052299) q[2];
sx q[2];
rz(-2.6504271) q[2];
sx q[2];
rz(0.14310357) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6306259) q[1];
sx q[1];
rz(-0.31019638) q[1];
sx q[1];
rz(-2.0621745) q[1];
rz(-1.2382965) q[3];
sx q[3];
rz(-2.1070814) q[3];
sx q[3];
rz(-0.09418776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4030054) q[2];
sx q[2];
rz(-1.0120729) q[2];
sx q[2];
rz(0.15815132) q[2];
rz(-0.75913298) q[3];
sx q[3];
rz(-1.682155) q[3];
sx q[3];
rz(0.37091836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89449969) q[0];
sx q[0];
rz(-2.5573754) q[0];
sx q[0];
rz(1.1868813) q[0];
rz(-0.12313708) q[1];
sx q[1];
rz(-1.5703399) q[1];
sx q[1];
rz(1.311696) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4224706) q[0];
sx q[0];
rz(-1.4234237) q[0];
sx q[0];
rz(-3.0585994) q[0];
x q[1];
rz(-2.4995828) q[2];
sx q[2];
rz(-0.12769709) q[2];
sx q[2];
rz(-2.0754432) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8183916) q[1];
sx q[1];
rz(-1.15266) q[1];
sx q[1];
rz(1.1365324) q[1];
rz(-1.7100641) q[3];
sx q[3];
rz(-0.0041120681) q[3];
sx q[3];
rz(-2.508956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.19795236) q[2];
sx q[2];
rz(-1.8215048) q[2];
sx q[2];
rz(0.8521592) q[2];
rz(-0.85396829) q[3];
sx q[3];
rz(-1.2205418) q[3];
sx q[3];
rz(2.061969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3581486) q[0];
sx q[0];
rz(-1.3842979) q[0];
sx q[0];
rz(3.070991) q[0];
rz(2.103503) q[1];
sx q[1];
rz(-1.9990653) q[1];
sx q[1];
rz(0.80931726) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3883078) q[0];
sx q[0];
rz(-1.6429672) q[0];
sx q[0];
rz(-0.57467527) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4368183) q[2];
sx q[2];
rz(-2.2597729) q[2];
sx q[2];
rz(2.5344684) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0505123) q[1];
sx q[1];
rz(-1.3111171) q[1];
sx q[1];
rz(-0.56118324) q[1];
x q[2];
rz(-1.7820939) q[3];
sx q[3];
rz(-0.70970067) q[3];
sx q[3];
rz(-3.078408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.91850963) q[2];
sx q[2];
rz(-2.5727936) q[2];
sx q[2];
rz(1.9522379) q[2];
rz(3.0009624) q[3];
sx q[3];
rz(-2.6719002) q[3];
sx q[3];
rz(-0.3895337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73914948) q[0];
sx q[0];
rz(-2.2279255) q[0];
sx q[0];
rz(1.5341349) q[0];
rz(-2.2588579) q[1];
sx q[1];
rz(-2.2969756) q[1];
sx q[1];
rz(0.61558634) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7587335) q[0];
sx q[0];
rz(-0.40399562) q[0];
sx q[0];
rz(0.23601471) q[0];
rz(-pi) q[1];
rz(0.87178715) q[2];
sx q[2];
rz(-1.0411556) q[2];
sx q[2];
rz(2.3123534) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5639685) q[1];
sx q[1];
rz(-1.7610487) q[1];
sx q[1];
rz(3.0012194) q[1];
x q[2];
rz(-3.1300486) q[3];
sx q[3];
rz(-0.95043105) q[3];
sx q[3];
rz(-0.64988713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9756056) q[2];
sx q[2];
rz(-2.1773982) q[2];
sx q[2];
rz(1.1791505) q[2];
rz(2.9278582) q[3];
sx q[3];
rz(-1.7515747) q[3];
sx q[3];
rz(-2.8960622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24446503) q[0];
sx q[0];
rz(-1.9944958) q[0];
sx q[0];
rz(-0.044064673) q[0];
rz(2.7524718) q[1];
sx q[1];
rz(-0.48946425) q[1];
sx q[1];
rz(-2.3997831) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76352966) q[0];
sx q[0];
rz(-2.6776095) q[0];
sx q[0];
rz(2.4876809) q[0];
x q[1];
rz(0.35105437) q[2];
sx q[2];
rz(-2.2864955) q[2];
sx q[2];
rz(-0.15049105) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.68315893) q[1];
sx q[1];
rz(-1.6583879) q[1];
sx q[1];
rz(-2.8187672) q[1];
rz(0.70176418) q[3];
sx q[3];
rz(-1.4150815) q[3];
sx q[3];
rz(-3.0879741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6645708) q[2];
sx q[2];
rz(-1.3571813) q[2];
sx q[2];
rz(0.50194293) q[2];
rz(-1.7160412) q[3];
sx q[3];
rz(-2.612096) q[3];
sx q[3];
rz(-2.3486923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75941706) q[0];
sx q[0];
rz(-1.734153) q[0];
sx q[0];
rz(0.39696804) q[0];
rz(-1.601903) q[1];
sx q[1];
rz(-0.27286467) q[1];
sx q[1];
rz(2.4698965) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3701733) q[0];
sx q[0];
rz(-2.015593) q[0];
sx q[0];
rz(3.0014787) q[0];
rz(-pi) q[1];
rz(2.256088) q[2];
sx q[2];
rz(-0.54965082) q[2];
sx q[2];
rz(-0.099121475) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.063350141) q[1];
sx q[1];
rz(-0.95253175) q[1];
sx q[1];
rz(1.9552014) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3109452) q[3];
sx q[3];
rz(-1.1201356) q[3];
sx q[3];
rz(-0.78307952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.73528618) q[2];
sx q[2];
rz(-2.6267509) q[2];
sx q[2];
rz(1.7740645) q[2];
rz(2.1788518) q[3];
sx q[3];
rz(-1.5906426) q[3];
sx q[3];
rz(-2.4875557) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.830736) q[0];
sx q[0];
rz(-1.5565358) q[0];
sx q[0];
rz(-0.32064015) q[0];
rz(-2.2036208) q[1];
sx q[1];
rz(-1.9357977) q[1];
sx q[1];
rz(-1.7373614) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2187248) q[0];
sx q[0];
rz(-1.6744807) q[0];
sx q[0];
rz(-1.4669048) q[0];
rz(-pi) q[1];
rz(-3.0530372) q[2];
sx q[2];
rz(-2.8635396) q[2];
sx q[2];
rz(-1.1403569) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.117489) q[1];
sx q[1];
rz(-1.5908594) q[1];
sx q[1];
rz(-1.5320918) q[1];
x q[2];
rz(2.3793166) q[3];
sx q[3];
rz(-1.395985) q[3];
sx q[3];
rz(-2.0559514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2001026) q[2];
sx q[2];
rz(-0.76277554) q[2];
sx q[2];
rz(0.74083677) q[2];
rz(-1.0624933) q[3];
sx q[3];
rz(-2.0334358) q[3];
sx q[3];
rz(-1.197193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9366539) q[0];
sx q[0];
rz(-2.7118201) q[0];
sx q[0];
rz(-0.75576654) q[0];
rz(1.8273805) q[1];
sx q[1];
rz(-1.6270437) q[1];
sx q[1];
rz(-2.4155713) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36011115) q[0];
sx q[0];
rz(-1.5364293) q[0];
sx q[0];
rz(0.78256677) q[0];
rz(-pi) q[1];
x q[1];
rz(0.88557247) q[2];
sx q[2];
rz(-1.6513593) q[2];
sx q[2];
rz(-1.1159814) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.0997114) q[1];
sx q[1];
rz(-1.2234294) q[1];
sx q[1];
rz(0.16176407) q[1];
rz(1.8733379) q[3];
sx q[3];
rz(-2.6740251) q[3];
sx q[3];
rz(-2.9085161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3296335) q[2];
sx q[2];
rz(-2.6128431) q[2];
sx q[2];
rz(-2.192396) q[2];
rz(2.9694929) q[3];
sx q[3];
rz(-2.1642919) q[3];
sx q[3];
rz(-2.7026091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
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
rz(-2.5934577) q[2];
sx q[2];
rz(-1.7036864) q[2];
sx q[2];
rz(1.7316935) q[2];
rz(-1.7023986) q[3];
sx q[3];
rz(-2.6489762) q[3];
sx q[3];
rz(1.4121644) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
