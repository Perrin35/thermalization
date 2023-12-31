OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8863949) q[0];
sx q[0];
rz(-1.2502517) q[0];
sx q[0];
rz(-1.3347081) q[0];
rz(2.788738) q[1];
sx q[1];
rz(3.3021441) q[1];
sx q[1];
rz(8.4488206) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1057518) q[0];
sx q[0];
rz(-2.5076206) q[0];
sx q[0];
rz(-2.7266154) q[0];
rz(-1.0078148) q[2];
sx q[2];
rz(-1.8554167) q[2];
sx q[2];
rz(-2.9993338) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0426892) q[1];
sx q[1];
rz(-0.96682036) q[1];
sx q[1];
rz(-2.9788115) q[1];
x q[2];
rz(-1.9777771) q[3];
sx q[3];
rz(-0.80899901) q[3];
sx q[3];
rz(0.33363261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7261937) q[2];
sx q[2];
rz(-1.4582863) q[2];
sx q[2];
rz(0.78380084) q[2];
rz(-0.33256724) q[3];
sx q[3];
rz(-2.9794897) q[3];
sx q[3];
rz(-1.8012841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48288229) q[0];
sx q[0];
rz(-1.4504526) q[0];
sx q[0];
rz(-2.9630307) q[0];
rz(-1.3372955) q[1];
sx q[1];
rz(-2.6006915) q[1];
sx q[1];
rz(-0.006342412) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77942383) q[0];
sx q[0];
rz(-2.9075025) q[0];
sx q[0];
rz(-0.97658821) q[0];
x q[1];
rz(-2.6132934) q[2];
sx q[2];
rz(-2.106278) q[2];
sx q[2];
rz(0.9165802) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.22390631) q[1];
sx q[1];
rz(-0.31653857) q[1];
sx q[1];
rz(2.9427337) q[1];
rz(-pi) q[2];
rz(-2.7553495) q[3];
sx q[3];
rz(-2.5442903) q[3];
sx q[3];
rz(2.951705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1938842) q[2];
sx q[2];
rz(-2.7145553) q[2];
sx q[2];
rz(-2.2634899) q[2];
rz(-0.86205035) q[3];
sx q[3];
rz(-2.2191007) q[3];
sx q[3];
rz(3.1356964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55643117) q[0];
sx q[0];
rz(-2.1142024) q[0];
sx q[0];
rz(0.01097824) q[0];
rz(0.36704656) q[1];
sx q[1];
rz(-1.234602) q[1];
sx q[1];
rz(-3.045851) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.028713) q[0];
sx q[0];
rz(-0.23447795) q[0];
sx q[0];
rz(-2.2857091) q[0];
rz(-2.4143366) q[2];
sx q[2];
rz(-0.63823344) q[2];
sx q[2];
rz(2.9122695) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0248191) q[1];
sx q[1];
rz(-2.1110592) q[1];
sx q[1];
rz(2.5953672) q[1];
x q[2];
rz(-1.5355574) q[3];
sx q[3];
rz(-2.1217151) q[3];
sx q[3];
rz(-0.86722022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0694971) q[2];
sx q[2];
rz(-2.3196689) q[2];
sx q[2];
rz(0.83479184) q[2];
rz(2.9299724) q[3];
sx q[3];
rz(-1.9112588) q[3];
sx q[3];
rz(2.766818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19462207) q[0];
sx q[0];
rz(-1.2397543) q[0];
sx q[0];
rz(-0.34657493) q[0];
rz(0.52571458) q[1];
sx q[1];
rz(-2.321967) q[1];
sx q[1];
rz(1.0338763) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3165986) q[0];
sx q[0];
rz(-2.1567417) q[0];
sx q[0];
rz(-0.98866776) q[0];
rz(2.1448574) q[2];
sx q[2];
rz(-2.0248932) q[2];
sx q[2];
rz(0.10276375) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.911072) q[1];
sx q[1];
rz(-1.6358747) q[1];
sx q[1];
rz(1.7622403) q[1];
x q[2];
rz(1.2437808) q[3];
sx q[3];
rz(-1.3033086) q[3];
sx q[3];
rz(0.063484065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.45450777) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(-1.3467849) q[2];
rz(0.421031) q[3];
sx q[3];
rz(-2.1249168) q[3];
sx q[3];
rz(0.69798654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7088257) q[0];
sx q[0];
rz(-2.3481752) q[0];
sx q[0];
rz(0.41473266) q[0];
rz(-1.3955836) q[1];
sx q[1];
rz(-2.4826629) q[1];
sx q[1];
rz(2.5674852) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5363279) q[0];
sx q[0];
rz(-0.61075532) q[0];
sx q[0];
rz(-1.3461793) q[0];
rz(-0.71530576) q[2];
sx q[2];
rz(-1.5878521) q[2];
sx q[2];
rz(-2.7280083) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.50358665) q[1];
sx q[1];
rz(-1.9203016) q[1];
sx q[1];
rz(-1.867678) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8168713) q[3];
sx q[3];
rz(-0.68021357) q[3];
sx q[3];
rz(2.084793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.85270143) q[2];
sx q[2];
rz(-0.39704278) q[2];
sx q[2];
rz(-1.627702) q[2];
rz(0.04143516) q[3];
sx q[3];
rz(-1.8824717) q[3];
sx q[3];
rz(-3.0630625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0742652) q[0];
sx q[0];
rz(-1.7794309) q[0];
sx q[0];
rz(-2.4434027) q[0];
rz(2.9746338) q[1];
sx q[1];
rz(-1.0792462) q[1];
sx q[1];
rz(0.73227698) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4740144) q[0];
sx q[0];
rz(-1.7828373) q[0];
sx q[0];
rz(-1.408512) q[0];
rz(0.90768355) q[2];
sx q[2];
rz(-0.92539061) q[2];
sx q[2];
rz(-1.4484608) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.60019833) q[1];
sx q[1];
rz(-1.7957893) q[1];
sx q[1];
rz(-0.66585983) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8090217) q[3];
sx q[3];
rz(-0.59623527) q[3];
sx q[3];
rz(-0.81070825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3884864) q[2];
sx q[2];
rz(-0.33375084) q[2];
sx q[2];
rz(1.1614655) q[2];
rz(0.30900624) q[3];
sx q[3];
rz(-1.2523264) q[3];
sx q[3];
rz(1.9741612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5724065) q[0];
sx q[0];
rz(-2.321406) q[0];
sx q[0];
rz(-0.15643315) q[0];
rz(2.6898443) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(-2.3715473) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5891588) q[0];
sx q[0];
rz(-1.7434412) q[0];
sx q[0];
rz(2.5347559) q[0];
x q[1];
rz(-2.7355843) q[2];
sx q[2];
rz(-1.0205262) q[2];
sx q[2];
rz(1.5232435) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.60939497) q[1];
sx q[1];
rz(-0.97071338) q[1];
sx q[1];
rz(2.1983912) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6386912) q[3];
sx q[3];
rz(-1.8637878) q[3];
sx q[3];
rz(-2.6754232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.6440755) q[2];
sx q[2];
rz(-0.13921177) q[2];
sx q[2];
rz(-1.0151803) q[2];
rz(-1.4366359) q[3];
sx q[3];
rz(-2.5865343) q[3];
sx q[3];
rz(-2.5456583) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62676936) q[0];
sx q[0];
rz(-2.5829782) q[0];
sx q[0];
rz(-2.8998937) q[0];
rz(0.73879755) q[1];
sx q[1];
rz(-2.6530478) q[1];
sx q[1];
rz(-2.7752005) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1358571) q[0];
sx q[0];
rz(-2.2861087) q[0];
sx q[0];
rz(2.4321796) q[0];
rz(1.1364469) q[2];
sx q[2];
rz(-0.98698101) q[2];
sx q[2];
rz(-0.39603147) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.39216343) q[1];
sx q[1];
rz(-2.2656419) q[1];
sx q[1];
rz(1.3543345) q[1];
rz(-pi) q[2];
rz(1.0189692) q[3];
sx q[3];
rz(-1.6120211) q[3];
sx q[3];
rz(3.105203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0176795) q[2];
sx q[2];
rz(-1.9423449) q[2];
sx q[2];
rz(2.8477342) q[2];
rz(-0.014523225) q[3];
sx q[3];
rz(-1.9745275) q[3];
sx q[3];
rz(0.6706388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83207399) q[0];
sx q[0];
rz(-0.80219769) q[0];
sx q[0];
rz(-0.88395399) q[0];
rz(2.4813095) q[1];
sx q[1];
rz(-0.98799223) q[1];
sx q[1];
rz(0.79137897) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26830772) q[0];
sx q[0];
rz(-2.377016) q[0];
sx q[0];
rz(-1.1029878) q[0];
x q[1];
rz(2.0558526) q[2];
sx q[2];
rz(-1.1527449) q[2];
sx q[2];
rz(2.5784091) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1014175) q[1];
sx q[1];
rz(-1.1244332) q[1];
sx q[1];
rz(1.1091713) q[1];
x q[2];
rz(1.8459686) q[3];
sx q[3];
rz(-2.5857946) q[3];
sx q[3];
rz(0.2273699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0749977) q[2];
sx q[2];
rz(-2.834088) q[2];
sx q[2];
rz(-1.0212612) q[2];
rz(2.7630473) q[3];
sx q[3];
rz(-0.56423855) q[3];
sx q[3];
rz(-0.59797257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1148949) q[0];
sx q[0];
rz(-0.4037936) q[0];
sx q[0];
rz(0.60780203) q[0];
rz(-2.9027477) q[1];
sx q[1];
rz(-0.74964476) q[1];
sx q[1];
rz(1.7609319) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5138807) q[0];
sx q[0];
rz(-0.6427592) q[0];
sx q[0];
rz(1.3520665) q[0];
rz(-pi) q[1];
x q[1];
rz(0.025823921) q[2];
sx q[2];
rz(-1.1354453) q[2];
sx q[2];
rz(0.045407427) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7271991) q[1];
sx q[1];
rz(-0.52800035) q[1];
sx q[1];
rz(0.90331932) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4019743) q[3];
sx q[3];
rz(-1.7572548) q[3];
sx q[3];
rz(-2.8846915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7878788) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(0.33622462) q[2];
rz(-2.0119038) q[3];
sx q[3];
rz(-1.7166398) q[3];
sx q[3];
rz(1.1415793) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0040141) q[0];
sx q[0];
rz(-1.7332358) q[0];
sx q[0];
rz(-1.9144203) q[0];
rz(2.5972988) q[1];
sx q[1];
rz(-1.2398564) q[1];
sx q[1];
rz(1.6287631) q[1];
rz(-0.33967321) q[2];
sx q[2];
rz(-2.2876231) q[2];
sx q[2];
rz(-3.0291578) q[2];
rz(1.2119157) q[3];
sx q[3];
rz(-0.59960312) q[3];
sx q[3];
rz(-3.0740769) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
