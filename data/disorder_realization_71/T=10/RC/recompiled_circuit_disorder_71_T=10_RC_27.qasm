OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.527737) q[0];
sx q[0];
rz(-1.4976488) q[0];
sx q[0];
rz(0.82984501) q[0];
rz(3.9217477) q[1];
sx q[1];
rz(5.2182066) q[1];
sx q[1];
rz(10.301104) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012955879) q[0];
sx q[0];
rz(-2.6814333) q[0];
sx q[0];
rz(-0.53928661) q[0];
x q[1];
rz(-2.949602) q[2];
sx q[2];
rz(-2.1477094) q[2];
sx q[2];
rz(0.38919762) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6688924) q[1];
sx q[1];
rz(-1.9007705) q[1];
sx q[1];
rz(-0.015923576) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7238486) q[3];
sx q[3];
rz(-1.9733841) q[3];
sx q[3];
rz(-0.11868417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.17065389) q[2];
sx q[2];
rz(-1.2761513) q[2];
sx q[2];
rz(2.412964) q[2];
rz(2.6206) q[3];
sx q[3];
rz(-0.96121585) q[3];
sx q[3];
rz(-2.9339824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.306863) q[0];
sx q[0];
rz(-1.1704209) q[0];
sx q[0];
rz(1.1215425) q[0];
rz(0.25575486) q[1];
sx q[1];
rz(-1.47822) q[1];
sx q[1];
rz(2.2671525) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16337285) q[0];
sx q[0];
rz(-1.997943) q[0];
sx q[0];
rz(0.33421974) q[0];
rz(-0.99526694) q[2];
sx q[2];
rz(-1.4412291) q[2];
sx q[2];
rz(-2.1339983) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2494617) q[1];
sx q[1];
rz(-0.79626894) q[1];
sx q[1];
rz(-0.65278058) q[1];
rz(0.053697649) q[3];
sx q[3];
rz(-1.3018381) q[3];
sx q[3];
rz(-2.4503436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4008537) q[2];
sx q[2];
rz(-2.6541371) q[2];
sx q[2];
rz(2.7056616) q[2];
rz(0.68108264) q[3];
sx q[3];
rz(-2.3705132) q[3];
sx q[3];
rz(0.40288231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23713672) q[0];
sx q[0];
rz(-2.853892) q[0];
sx q[0];
rz(2.0667734) q[0];
rz(2.3020321) q[1];
sx q[1];
rz(-2.3222175) q[1];
sx q[1];
rz(0.39594617) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.031035) q[0];
sx q[0];
rz(-1.9651411) q[0];
sx q[0];
rz(2.1179384) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.01404889) q[2];
sx q[2];
rz(-1.8648246) q[2];
sx q[2];
rz(-2.369898) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.695895) q[1];
sx q[1];
rz(-0.88765111) q[1];
sx q[1];
rz(-1.5831069) q[1];
rz(-3.0523236) q[3];
sx q[3];
rz(-2.8850728) q[3];
sx q[3];
rz(1.741011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5376771) q[2];
sx q[2];
rz(-1.8751514) q[2];
sx q[2];
rz(-0.23920693) q[2];
rz(0.075332969) q[3];
sx q[3];
rz(-1.1139261) q[3];
sx q[3];
rz(-1.3031561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4199715) q[0];
sx q[0];
rz(-2.2890838) q[0];
sx q[0];
rz(-0.81992942) q[0];
rz(0.48768249) q[1];
sx q[1];
rz(-0.90355021) q[1];
sx q[1];
rz(2.908169) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6452091) q[0];
sx q[0];
rz(-1.6312946) q[0];
sx q[0];
rz(-1.4738826) q[0];
rz(1.3350305) q[2];
sx q[2];
rz(-2.1996017) q[2];
sx q[2];
rz(0.51970383) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9277716) q[1];
sx q[1];
rz(-1.908761) q[1];
sx q[1];
rz(-2.4073699) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38846429) q[3];
sx q[3];
rz(-0.82287517) q[3];
sx q[3];
rz(0.32271656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0507811) q[2];
sx q[2];
rz(-1.2819141) q[2];
sx q[2];
rz(2.3941669) q[2];
rz(-2.9181972) q[3];
sx q[3];
rz(-2.5441393) q[3];
sx q[3];
rz(-2.8994765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4500047) q[0];
sx q[0];
rz(-1.7084028) q[0];
sx q[0];
rz(-0.064237021) q[0];
rz(2.1977987) q[1];
sx q[1];
rz(-2.4143024) q[1];
sx q[1];
rz(-0.76104004) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8262186) q[0];
sx q[0];
rz(-1.5407019) q[0];
sx q[0];
rz(1.8316359) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3741578) q[2];
sx q[2];
rz(-0.27786294) q[2];
sx q[2];
rz(2.1674736) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3984822) q[1];
sx q[1];
rz(-2.8201436) q[1];
sx q[1];
rz(0.076996315) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.28273545) q[3];
sx q[3];
rz(-1.3021819) q[3];
sx q[3];
rz(-0.80073592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.80660194) q[2];
sx q[2];
rz(-2.0709753) q[2];
sx q[2];
rz(-0.09207329) q[2];
rz(2.4798685) q[3];
sx q[3];
rz(-2.3501553) q[3];
sx q[3];
rz(-1.7592643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46835607) q[0];
sx q[0];
rz(-1.3156923) q[0];
sx q[0];
rz(-3.1150505) q[0];
rz(-2.2684855) q[1];
sx q[1];
rz(-1.1353506) q[1];
sx q[1];
rz(-2.81566) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8317141) q[0];
sx q[0];
rz(-1.9448115) q[0];
sx q[0];
rz(-0.563234) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2665777) q[2];
sx q[2];
rz(-1.9982669) q[2];
sx q[2];
rz(1.8311335) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.95820108) q[1];
sx q[1];
rz(-0.32918731) q[1];
sx q[1];
rz(-2.7126461) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16485729) q[3];
sx q[3];
rz(-0.10247173) q[3];
sx q[3];
rz(-2.7285679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.59763336) q[2];
sx q[2];
rz(-1.8240857) q[2];
sx q[2];
rz(-0.65417543) q[2];
rz(1.7116961) q[3];
sx q[3];
rz(-1.9734029) q[3];
sx q[3];
rz(2.6005319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8687246) q[0];
sx q[0];
rz(-1.4725715) q[0];
sx q[0];
rz(0.72189271) q[0];
rz(-1.4121274) q[1];
sx q[1];
rz(-2.3528603) q[1];
sx q[1];
rz(0.11925764) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3463466) q[0];
sx q[0];
rz(-2.3756785) q[0];
sx q[0];
rz(0.38608293) q[0];
rz(-pi) q[1];
rz(-3.0316333) q[2];
sx q[2];
rz(-1.7320776) q[2];
sx q[2];
rz(1.9732628) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.81469401) q[1];
sx q[1];
rz(-2.6242995) q[1];
sx q[1];
rz(-1.3002404) q[1];
rz(-pi) q[2];
rz(-1.9414385) q[3];
sx q[3];
rz(-0.89768411) q[3];
sx q[3];
rz(2.9692269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6922336) q[2];
sx q[2];
rz(-1.3204152) q[2];
sx q[2];
rz(-1.3593486) q[2];
rz(2.3826777) q[3];
sx q[3];
rz(-0.24154285) q[3];
sx q[3];
rz(-2.6045077) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3100202) q[0];
sx q[0];
rz(-2.4825403) q[0];
sx q[0];
rz(2.0781562) q[0];
rz(2.8670782) q[1];
sx q[1];
rz(-1.9332705) q[1];
sx q[1];
rz(-2.2559821) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0793593) q[0];
sx q[0];
rz(-1.0842807) q[0];
sx q[0];
rz(-1.3079206) q[0];
rz(1.5947072) q[2];
sx q[2];
rz(-0.57966053) q[2];
sx q[2];
rz(-1.5469345) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.209621) q[1];
sx q[1];
rz(-2.0721657) q[1];
sx q[1];
rz(-0.45100905) q[1];
rz(-pi) q[2];
x q[2];
rz(0.76836821) q[3];
sx q[3];
rz(-1.0230912) q[3];
sx q[3];
rz(2.2636569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.72835913) q[2];
sx q[2];
rz(-0.76247549) q[2];
sx q[2];
rz(1.1317066) q[2];
rz(-2.0570095) q[3];
sx q[3];
rz(-2.0621433) q[3];
sx q[3];
rz(-1.926698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8383012) q[0];
sx q[0];
rz(-1.6475995) q[0];
sx q[0];
rz(2.0595179) q[0];
rz(1.8661631) q[1];
sx q[1];
rz(-1.0042896) q[1];
sx q[1];
rz(-2.0057604) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58763114) q[0];
sx q[0];
rz(-0.56092867) q[0];
sx q[0];
rz(1.2558623) q[0];
rz(2.7068107) q[2];
sx q[2];
rz(-1.5506622) q[2];
sx q[2];
rz(0.15205631) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.51582) q[1];
sx q[1];
rz(-1.9966218) q[1];
sx q[1];
rz(-1.6243402) q[1];
x q[2];
rz(-0.42580749) q[3];
sx q[3];
rz(-0.96447456) q[3];
sx q[3];
rz(0.013289277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.11848005) q[2];
sx q[2];
rz(-1.2597522) q[2];
sx q[2];
rz(-0.87289587) q[2];
rz(2.2980799) q[3];
sx q[3];
rz(-2.8819363) q[3];
sx q[3];
rz(0.18994722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2492367) q[0];
sx q[0];
rz(-1.7974239) q[0];
sx q[0];
rz(-1.2783485) q[0];
rz(2.1168013) q[1];
sx q[1];
rz(-2.0147851) q[1];
sx q[1];
rz(-1.9445673) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68162936) q[0];
sx q[0];
rz(-2.3869793) q[0];
sx q[0];
rz(-0.90453903) q[0];
rz(-1.3329266) q[2];
sx q[2];
rz(-1.4537721) q[2];
sx q[2];
rz(0.34840096) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4492606) q[1];
sx q[1];
rz(-2.3093866) q[1];
sx q[1];
rz(0.87528054) q[1];
x q[2];
rz(2.5246546) q[3];
sx q[3];
rz(-0.60041282) q[3];
sx q[3];
rz(2.2850349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3060351) q[2];
sx q[2];
rz(-2.4287537) q[2];
sx q[2];
rz(-0.79997921) q[2];
rz(-1.9647313) q[3];
sx q[3];
rz(-1.7088944) q[3];
sx q[3];
rz(2.1879788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4326614) q[0];
sx q[0];
rz(-0.14871696) q[0];
sx q[0];
rz(0.8014252) q[0];
rz(-2.6196383) q[1];
sx q[1];
rz(-0.83871651) q[1];
sx q[1];
rz(-2.9768859) q[1];
rz(-0.75688731) q[2];
sx q[2];
rz(-0.49513985) q[2];
sx q[2];
rz(-2.0078299) q[2];
rz(2.342631) q[3];
sx q[3];
rz(-0.81294717) q[3];
sx q[3];
rz(-3.0019928) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];