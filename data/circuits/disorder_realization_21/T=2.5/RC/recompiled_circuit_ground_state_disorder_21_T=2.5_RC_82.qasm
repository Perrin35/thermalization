OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2709687) q[0];
sx q[0];
rz(-0.55611098) q[0];
sx q[0];
rz(2.1882353) q[0];
rz(0.019729992) q[1];
sx q[1];
rz(1.035773) q[1];
sx q[1];
rz(8.4761578) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9069288) q[0];
sx q[0];
rz(-1.0415823) q[0];
sx q[0];
rz(-2.0509023) q[0];
rz(1.5736012) q[2];
sx q[2];
rz(-1.8082779) q[2];
sx q[2];
rz(-1.5495007) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1780121) q[1];
sx q[1];
rz(-0.90998703) q[1];
sx q[1];
rz(-1.0170487) q[1];
rz(-pi) q[2];
rz(0.3462195) q[3];
sx q[3];
rz(-2.0967212) q[3];
sx q[3];
rz(0.093737515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3775776) q[2];
sx q[2];
rz(-3.0674051) q[2];
sx q[2];
rz(0.78262502) q[2];
rz(3.022656) q[3];
sx q[3];
rz(-2.1075893) q[3];
sx q[3];
rz(2.9940166) q[3];
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
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19621944) q[0];
sx q[0];
rz(-2.0202899) q[0];
sx q[0];
rz(-2.8837606) q[0];
rz(3.0505772) q[1];
sx q[1];
rz(-1.0992522) q[1];
sx q[1];
rz(1.4965422) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89423075) q[0];
sx q[0];
rz(-1.5024439) q[0];
sx q[0];
rz(2.1240881) q[0];
rz(-pi) q[1];
rz(1.2256669) q[2];
sx q[2];
rz(-1.218443) q[2];
sx q[2];
rz(0.27252588) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2014034) q[1];
sx q[1];
rz(-1.1588105) q[1];
sx q[1];
rz(-3.1088016) q[1];
rz(2.4739059) q[3];
sx q[3];
rz(-1.2419309) q[3];
sx q[3];
rz(-1.2038976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1220793) q[2];
sx q[2];
rz(-1.2747108) q[2];
sx q[2];
rz(2.7369734) q[2];
rz(0.98958611) q[3];
sx q[3];
rz(-1.8414958) q[3];
sx q[3];
rz(2.5185481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0290381) q[0];
sx q[0];
rz(-0.14415388) q[0];
sx q[0];
rz(-0.83830225) q[0];
rz(-0.84838947) q[1];
sx q[1];
rz(-1.9219857) q[1];
sx q[1];
rz(-0.99260509) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86322376) q[0];
sx q[0];
rz(-1.7107757) q[0];
sx q[0];
rz(2.4158258) q[0];
x q[1];
rz(-1.6279334) q[2];
sx q[2];
rz(-1.2369725) q[2];
sx q[2];
rz(1.305507) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5287244) q[1];
sx q[1];
rz(-0.49002417) q[1];
sx q[1];
rz(-2.571066) q[1];
x q[2];
rz(-1.6818524) q[3];
sx q[3];
rz(-1.8594311) q[3];
sx q[3];
rz(1.9185324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.8695716) q[2];
sx q[2];
rz(-1.1676936) q[2];
sx q[2];
rz(2.0167548) q[2];
rz(0.62075067) q[3];
sx q[3];
rz(-2.1706332) q[3];
sx q[3];
rz(-3.0264405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89897412) q[0];
sx q[0];
rz(-1.9816575) q[0];
sx q[0];
rz(2.9677891) q[0];
rz(2.4558892) q[1];
sx q[1];
rz(-1.655429) q[1];
sx q[1];
rz(0.83820835) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3371171) q[0];
sx q[0];
rz(-0.67401471) q[0];
sx q[0];
rz(-2.3524257) q[0];
rz(-pi) q[1];
rz(-2.4160552) q[2];
sx q[2];
rz(-0.71070403) q[2];
sx q[2];
rz(-2.2827471) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.9058075) q[1];
sx q[1];
rz(-2.2407994) q[1];
sx q[1];
rz(-1.3895821) q[1];
x q[2];
rz(0.46282262) q[3];
sx q[3];
rz(-1.5372542) q[3];
sx q[3];
rz(2.7963146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.54245) q[2];
sx q[2];
rz(-1.6833545) q[2];
sx q[2];
rz(-0.35169265) q[2];
rz(-1.1951949) q[3];
sx q[3];
rz(-1.1870793) q[3];
sx q[3];
rz(-3.1174507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46566063) q[0];
sx q[0];
rz(-2.6213578) q[0];
sx q[0];
rz(-0.35292536) q[0];
rz(-2.832761) q[1];
sx q[1];
rz(-2.1052723) q[1];
sx q[1];
rz(0.028506361) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4052947) q[0];
sx q[0];
rz(-0.60966821) q[0];
sx q[0];
rz(1.2627312) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79645313) q[2];
sx q[2];
rz(-2.6203794) q[2];
sx q[2];
rz(-3.1375743) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9526279) q[1];
sx q[1];
rz(-2.3213534) q[1];
sx q[1];
rz(2.8431975) q[1];
x q[2];
rz(2.0621544) q[3];
sx q[3];
rz(-0.48499987) q[3];
sx q[3];
rz(0.62937832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.21294022) q[2];
sx q[2];
rz(-3.0781367) q[2];
sx q[2];
rz(0.97079903) q[2];
rz(1.0468696) q[3];
sx q[3];
rz(-1.6887083) q[3];
sx q[3];
rz(2.651732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30188072) q[0];
sx q[0];
rz(-1.3588926) q[0];
sx q[0];
rz(-0.090959892) q[0];
rz(1.92314) q[1];
sx q[1];
rz(-0.33088845) q[1];
sx q[1];
rz(-2.5023696) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1341742) q[0];
sx q[0];
rz(-1.245541) q[0];
sx q[0];
rz(-1.4283309) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84750643) q[2];
sx q[2];
rz(-2.9581262) q[2];
sx q[2];
rz(3.1364721) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7057695) q[1];
sx q[1];
rz(-2.2652049) q[1];
sx q[1];
rz(-0.69470508) q[1];
x q[2];
rz(1.8973783) q[3];
sx q[3];
rz(-1.5726798) q[3];
sx q[3];
rz(2.7060946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.062139221) q[2];
sx q[2];
rz(-2.1663351) q[2];
sx q[2];
rz(-2.6677168) q[2];
rz(-1.8114932) q[3];
sx q[3];
rz(-0.88089839) q[3];
sx q[3];
rz(2.7661095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5652931) q[0];
sx q[0];
rz(-1.0813035) q[0];
sx q[0];
rz(-2.9190049) q[0];
rz(1.0635771) q[1];
sx q[1];
rz(-1.56366) q[1];
sx q[1];
rz(-1.9482013) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28119606) q[0];
sx q[0];
rz(-1.687107) q[0];
sx q[0];
rz(2.3915245) q[0];
rz(1.7980099) q[2];
sx q[2];
rz(-1.4526723) q[2];
sx q[2];
rz(0.83641499) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.64456731) q[1];
sx q[1];
rz(-0.46176061) q[1];
sx q[1];
rz(-0.30739947) q[1];
rz(2.1126595) q[3];
sx q[3];
rz(-0.30908424) q[3];
sx q[3];
rz(2.2268471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0003164) q[2];
sx q[2];
rz(-2.2681984) q[2];
sx q[2];
rz(2.5013962) q[2];
rz(0.021942465) q[3];
sx q[3];
rz(-1.4098189) q[3];
sx q[3];
rz(-2.791361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6496277) q[0];
sx q[0];
rz(-1.7438629) q[0];
sx q[0];
rz(-0.52038991) q[0];
rz(-2.261816) q[1];
sx q[1];
rz(-1.9089411) q[1];
sx q[1];
rz(-2.2487776) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8917903) q[0];
sx q[0];
rz(-0.098488657) q[0];
sx q[0];
rz(-2.6633224) q[0];
x q[1];
rz(-0.8142578) q[2];
sx q[2];
rz(-2.149048) q[2];
sx q[2];
rz(-1.1922497) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0860111) q[1];
sx q[1];
rz(-2.308929) q[1];
sx q[1];
rz(0.38565454) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51629169) q[3];
sx q[3];
rz(-1.9441351) q[3];
sx q[3];
rz(-0.0084127154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8799379) q[2];
sx q[2];
rz(-2.0377906) q[2];
sx q[2];
rz(-3.0547764) q[2];
rz(-0.9451198) q[3];
sx q[3];
rz(-2.7961531) q[3];
sx q[3];
rz(1.1300348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.9574808) q[0];
sx q[0];
rz(-0.090395398) q[0];
sx q[0];
rz(0.064067319) q[0];
rz(2.0059026) q[1];
sx q[1];
rz(-1.7013197) q[1];
sx q[1];
rz(-0.51220977) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6653888) q[0];
sx q[0];
rz(-1.8955766) q[0];
sx q[0];
rz(-2.7074714) q[0];
x q[1];
rz(-0.86821235) q[2];
sx q[2];
rz(-1.3033452) q[2];
sx q[2];
rz(-0.78154678) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.87742245) q[1];
sx q[1];
rz(-1.8296731) q[1];
sx q[1];
rz(-1.6593462) q[1];
rz(-pi) q[2];
rz(-1.5548315) q[3];
sx q[3];
rz(-1.1619688) q[3];
sx q[3];
rz(-0.4086993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8168489) q[2];
sx q[2];
rz(-0.73255676) q[2];
sx q[2];
rz(2.0762439) q[2];
rz(-1.6522853) q[3];
sx q[3];
rz(-2.1842712) q[3];
sx q[3];
rz(0.092441946) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5548993) q[0];
sx q[0];
rz(-2.7296992) q[0];
sx q[0];
rz(0.33083415) q[0];
rz(1.9031485) q[1];
sx q[1];
rz(-2.5336783) q[1];
sx q[1];
rz(-2.8668561) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.559812) q[0];
sx q[0];
rz(-1.1413478) q[0];
sx q[0];
rz(0.37508022) q[0];
rz(0.27336911) q[2];
sx q[2];
rz(-1.6022575) q[2];
sx q[2];
rz(2.4431369) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6152321) q[1];
sx q[1];
rz(-0.52605275) q[1];
sx q[1];
rz(1.2895209) q[1];
rz(-0.072321691) q[3];
sx q[3];
rz(-2.157271) q[3];
sx q[3];
rz(-0.14432913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1642509) q[2];
sx q[2];
rz(-2.3873316) q[2];
sx q[2];
rz(1.784262) q[2];
rz(2.6409798) q[3];
sx q[3];
rz(-1.5784135) q[3];
sx q[3];
rz(-1.64465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5951344) q[0];
sx q[0];
rz(-1.559579) q[0];
sx q[0];
rz(2.5478242) q[0];
rz(-3.0733227) q[1];
sx q[1];
rz(-1.9207813) q[1];
sx q[1];
rz(-3.1178738) q[1];
rz(2.7958128) q[2];
sx q[2];
rz(-1.7263392) q[2];
sx q[2];
rz(1.6172258) q[2];
rz(-2.9636737) q[3];
sx q[3];
rz(-1.3861227) q[3];
sx q[3];
rz(1.9021195) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
