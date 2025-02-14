OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2312343) q[0];
sx q[0];
rz(-0.8690106) q[0];
sx q[0];
rz(1.0847217) q[0];
rz(-1.1551069) q[1];
sx q[1];
rz(-0.81973633) q[1];
sx q[1];
rz(-0.91135946) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7408352) q[0];
sx q[0];
rz(-0.82159737) q[0];
sx q[0];
rz(1.2375968) q[0];
rz(-1.4248104) q[2];
sx q[2];
rz(-0.8643736) q[2];
sx q[2];
rz(1.3024769) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5973917) q[1];
sx q[1];
rz(-2.4185838) q[1];
sx q[1];
rz(0.71031481) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0769071) q[3];
sx q[3];
rz(-2.0035335) q[3];
sx q[3];
rz(0.81165867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.34065166) q[2];
sx q[2];
rz(-2.030535) q[2];
sx q[2];
rz(-3.0621373) q[2];
rz(0.60845145) q[3];
sx q[3];
rz(-1.918101) q[3];
sx q[3];
rz(-0.89103812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4873753) q[0];
sx q[0];
rz(-0.86367622) q[0];
sx q[0];
rz(-1.3551711) q[0];
rz(-1.0379418) q[1];
sx q[1];
rz(-0.67960056) q[1];
sx q[1];
rz(-2.3341446) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3231182) q[0];
sx q[0];
rz(-1.2809296) q[0];
sx q[0];
rz(0.94375837) q[0];
rz(2.1464426) q[2];
sx q[2];
rz(-1.2479221) q[2];
sx q[2];
rz(0.49597464) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3088919) q[1];
sx q[1];
rz(-1.8101793) q[1];
sx q[1];
rz(1.6390087) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2427166) q[3];
sx q[3];
rz(-0.17016093) q[3];
sx q[3];
rz(2.8191322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.21343931) q[2];
sx q[2];
rz(-1.5779053) q[2];
sx q[2];
rz(1.6820924) q[2];
rz(-2.6835119) q[3];
sx q[3];
rz(-2.5638678) q[3];
sx q[3];
rz(-0.95988449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15993519) q[0];
sx q[0];
rz(-1.0502879) q[0];
sx q[0];
rz(-0.67498573) q[0];
rz(0.58810294) q[1];
sx q[1];
rz(-0.85398483) q[1];
sx q[1];
rz(-1.810422) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0494306) q[0];
sx q[0];
rz(-1.9665008) q[0];
sx q[0];
rz(0.53073287) q[0];
rz(-pi) q[1];
rz(2.0841062) q[2];
sx q[2];
rz(-0.9710487) q[2];
sx q[2];
rz(2.5059003) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2053255) q[1];
sx q[1];
rz(-0.85204879) q[1];
sx q[1];
rz(-2.8014552) q[1];
rz(-pi) q[2];
rz(2.939091) q[3];
sx q[3];
rz(-1.080386) q[3];
sx q[3];
rz(-2.5449076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5171234) q[2];
sx q[2];
rz(-2.7915967) q[2];
sx q[2];
rz(-0.2207174) q[2];
rz(0.80495009) q[3];
sx q[3];
rz(-1.5419818) q[3];
sx q[3];
rz(1.8055003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6014366) q[0];
sx q[0];
rz(-0.24208459) q[0];
sx q[0];
rz(-2.5228187) q[0];
rz(-0.082911804) q[1];
sx q[1];
rz(-0.34268788) q[1];
sx q[1];
rz(1.5203016) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5280209) q[0];
sx q[0];
rz(-0.53410406) q[0];
sx q[0];
rz(-1.7604802) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5624406) q[2];
sx q[2];
rz(-1.3716099) q[2];
sx q[2];
rz(2.3483089) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6449086) q[1];
sx q[1];
rz(-1.0606979) q[1];
sx q[1];
rz(-0.094385191) q[1];
rz(-pi) q[2];
rz(-1.3906562) q[3];
sx q[3];
rz(-1.2015752) q[3];
sx q[3];
rz(-0.1843017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.80359047) q[2];
sx q[2];
rz(-3.0323196) q[2];
sx q[2];
rz(-2.9962311) q[2];
rz(0.27211443) q[3];
sx q[3];
rz(-1.4067255) q[3];
sx q[3];
rz(-1.1468148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1266601) q[0];
sx q[0];
rz(-1.8155875) q[0];
sx q[0];
rz(3.0354011) q[0];
rz(2.1821678) q[1];
sx q[1];
rz(-2.5966849) q[1];
sx q[1];
rz(2.9471961) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4721699) q[0];
sx q[0];
rz(-1.484777) q[0];
sx q[0];
rz(-2.4079313) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4559559) q[2];
sx q[2];
rz(-1.5756338) q[2];
sx q[2];
rz(-2.9034535) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0982617) q[1];
sx q[1];
rz(-2.1594467) q[1];
sx q[1];
rz(1.5215988) q[1];
x q[2];
rz(0.077094519) q[3];
sx q[3];
rz(-1.1655679) q[3];
sx q[3];
rz(-0.73540686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.83547366) q[2];
sx q[2];
rz(-1.4543616) q[2];
sx q[2];
rz(-0.56841889) q[2];
rz(3.0889555) q[3];
sx q[3];
rz(-1.5024374) q[3];
sx q[3];
rz(-0.13681017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1133872) q[0];
sx q[0];
rz(-0.55332342) q[0];
sx q[0];
rz(-2.9521039) q[0];
rz(1.0264617) q[1];
sx q[1];
rz(-0.84954134) q[1];
sx q[1];
rz(2.386327) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3657235) q[0];
sx q[0];
rz(-2.2138174) q[0];
sx q[0];
rz(-2.4588206) q[0];
rz(2.6195141) q[2];
sx q[2];
rz(-1.0713801) q[2];
sx q[2];
rz(-1.5273818) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.64769563) q[1];
sx q[1];
rz(-1.4129479) q[1];
sx q[1];
rz(1.3998019) q[1];
rz(-pi) q[2];
rz(-2.7868458) q[3];
sx q[3];
rz(-1.915853) q[3];
sx q[3];
rz(1.7999032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3834164) q[2];
sx q[2];
rz(-1.5551609) q[2];
sx q[2];
rz(1.7002534) q[2];
rz(1.3759184) q[3];
sx q[3];
rz(-0.91840363) q[3];
sx q[3];
rz(2.4413696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7021084) q[0];
sx q[0];
rz(-1.0478042) q[0];
sx q[0];
rz(1.9973607) q[0];
rz(2.8817835) q[1];
sx q[1];
rz(-0.83994284) q[1];
sx q[1];
rz(0.98794404) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0555094) q[0];
sx q[0];
rz(-1.4713773) q[0];
sx q[0];
rz(1.5539597) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6778474) q[2];
sx q[2];
rz(-0.79163359) q[2];
sx q[2];
rz(0.82198696) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0007947) q[1];
sx q[1];
rz(-0.70604815) q[1];
sx q[1];
rz(-2.7435859) q[1];
x q[2];
rz(-2.8389205) q[3];
sx q[3];
rz(-2.1244123) q[3];
sx q[3];
rz(0.68676567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40954956) q[2];
sx q[2];
rz(-1.6403551) q[2];
sx q[2];
rz(0.6130971) q[2];
rz(2.4260855) q[3];
sx q[3];
rz(-0.77176538) q[3];
sx q[3];
rz(-2.7806921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.9382984) q[0];
sx q[0];
rz(-2.2977915) q[0];
sx q[0];
rz(-0.26813689) q[0];
rz(0.2306436) q[1];
sx q[1];
rz(-1.5985039) q[1];
sx q[1];
rz(0.014009744) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1582571) q[0];
sx q[0];
rz(-1.3620485) q[0];
sx q[0];
rz(0.76738417) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1300541) q[2];
sx q[2];
rz(-1.7157555) q[2];
sx q[2];
rz(-0.10565378) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.18302984) q[1];
sx q[1];
rz(-1.5261478) q[1];
sx q[1];
rz(-0.53739287) q[1];
rz(-pi) q[2];
rz(1.713578) q[3];
sx q[3];
rz(-1.6001978) q[3];
sx q[3];
rz(0.93522302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1676499) q[2];
sx q[2];
rz(-1.3385945) q[2];
sx q[2];
rz(-0.30926427) q[2];
rz(-0.15657982) q[3];
sx q[3];
rz(-1.4045818) q[3];
sx q[3];
rz(1.0369161) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2579047) q[0];
sx q[0];
rz(-0.64249277) q[0];
sx q[0];
rz(-0.25074348) q[0];
rz(0.58654395) q[1];
sx q[1];
rz(-1.5637014) q[1];
sx q[1];
rz(2.3768545) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4737807) q[0];
sx q[0];
rz(-2.7062253) q[0];
sx q[0];
rz(2.3171168) q[0];
x q[1];
rz(-0.7711507) q[2];
sx q[2];
rz(-1.7131107) q[2];
sx q[2];
rz(-2.3754295) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.918461) q[1];
sx q[1];
rz(-1.4479965) q[1];
sx q[1];
rz(-1.5801013) q[1];
rz(-pi) q[2];
rz(2.7627691) q[3];
sx q[3];
rz(-1.6215542) q[3];
sx q[3];
rz(1.7876884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5137198) q[2];
sx q[2];
rz(-2.3508115) q[2];
sx q[2];
rz(-2.9676843) q[2];
rz(2.3993313) q[3];
sx q[3];
rz(-2.3895013) q[3];
sx q[3];
rz(-3.098587) q[3];
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
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0088418) q[0];
sx q[0];
rz(-0.77021563) q[0];
sx q[0];
rz(0.26790628) q[0];
rz(-1.8065037) q[1];
sx q[1];
rz(-0.91064149) q[1];
sx q[1];
rz(1.3722027) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1233032) q[0];
sx q[0];
rz(-0.51890131) q[0];
sx q[0];
rz(-2.0730134) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8453061) q[2];
sx q[2];
rz(-1.6673267) q[2];
sx q[2];
rz(1.5636843) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7160373) q[1];
sx q[1];
rz(-2.4074984) q[1];
sx q[1];
rz(-0.74931637) q[1];
x q[2];
rz(-0.80268152) q[3];
sx q[3];
rz(-0.44375989) q[3];
sx q[3];
rz(-2.8783523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.043896) q[2];
sx q[2];
rz(-1.0363657) q[2];
sx q[2];
rz(1.3807266) q[2];
rz(-1.1287639) q[3];
sx q[3];
rz(-0.94996101) q[3];
sx q[3];
rz(-1.5560163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8783405) q[0];
sx q[0];
rz(-1.5338407) q[0];
sx q[0];
rz(-1.1080678) q[0];
rz(0.68645984) q[1];
sx q[1];
rz(-1.7484799) q[1];
sx q[1];
rz(1.9304986) q[1];
rz(-2.6471967) q[2];
sx q[2];
rz(-1.313268) q[2];
sx q[2];
rz(3.1368844) q[2];
rz(-2.108528) q[3];
sx q[3];
rz(-1.929639) q[3];
sx q[3];
rz(0.7076984) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
