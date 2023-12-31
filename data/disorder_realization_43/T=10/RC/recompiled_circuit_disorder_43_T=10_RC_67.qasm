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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1057518) q[0];
sx q[0];
rz(-2.5076206) q[0];
sx q[0];
rz(0.41497725) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33306723) q[2];
sx q[2];
rz(-1.0330079) q[2];
sx q[2];
rz(-1.2531467) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9589899) q[1];
sx q[1];
rz(-2.5187153) q[1];
sx q[1];
rz(-1.8014924) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39335143) q[3];
sx q[3];
rz(-2.297612) q[3];
sx q[3];
rz(2.2497183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41539899) q[2];
sx q[2];
rz(-1.6833064) q[2];
sx q[2];
rz(0.78380084) q[2];
rz(-2.8090254) q[3];
sx q[3];
rz(-0.16210292) q[3];
sx q[3];
rz(-1.8012841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48288229) q[0];
sx q[0];
rz(-1.4504526) q[0];
sx q[0];
rz(2.9630307) q[0];
rz(1.8042971) q[1];
sx q[1];
rz(-0.54090118) q[1];
sx q[1];
rz(-3.1352502) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20992499) q[0];
sx q[0];
rz(-1.7010265) q[0];
sx q[0];
rz(-1.3757214) q[0];
rz(0.86654051) q[2];
sx q[2];
rz(-2.4079977) q[2];
sx q[2];
rz(-1.3726335) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1576707) q[1];
sx q[1];
rz(-1.6323286) q[1];
sx q[1];
rz(0.31069777) q[1];
rz(-pi) q[2];
rz(2.7553495) q[3];
sx q[3];
rz(-2.5442903) q[3];
sx q[3];
rz(-2.951705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1938842) q[2];
sx q[2];
rz(-0.42703736) q[2];
sx q[2];
rz(0.87810278) q[2];
rz(-2.2795423) q[3];
sx q[3];
rz(-0.92249191) q[3];
sx q[3];
rz(-0.0058962065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55643117) q[0];
sx q[0];
rz(-1.0273902) q[0];
sx q[0];
rz(3.1306144) q[0];
rz(-0.36704656) q[1];
sx q[1];
rz(-1.9069907) q[1];
sx q[1];
rz(-3.045851) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.028713) q[0];
sx q[0];
rz(-0.23447795) q[0];
sx q[0];
rz(2.2857091) q[0];
rz(-pi) q[1];
rz(2.6355866) q[2];
sx q[2];
rz(-1.9780469) q[2];
sx q[2];
rz(-2.4207052) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1167736) q[1];
sx q[1];
rz(-2.1110592) q[1];
sx q[1];
rz(2.5953672) q[1];
x q[2];
rz(3.0843094) q[3];
sx q[3];
rz(-0.55192845) q[3];
sx q[3];
rz(0.79997593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0694971) q[2];
sx q[2];
rz(-2.3196689) q[2];
sx q[2];
rz(2.3068008) q[2];
rz(-2.9299724) q[3];
sx q[3];
rz(-1.2303338) q[3];
sx q[3];
rz(2.766818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19462207) q[0];
sx q[0];
rz(-1.2397543) q[0];
sx q[0];
rz(-2.7950177) q[0];
rz(2.6158781) q[1];
sx q[1];
rz(-0.81962568) q[1];
sx q[1];
rz(-2.1077164) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.046712) q[0];
sx q[0];
rz(-1.0949507) q[0];
sx q[0];
rz(-0.67142077) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1448574) q[2];
sx q[2];
rz(-2.0248932) q[2];
sx q[2];
rz(3.0388289) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.66400601) q[1];
sx q[1];
rz(-2.9395182) q[1];
sx q[1];
rz(-1.2408153) q[1];
x q[2];
rz(0.86446188) q[3];
sx q[3];
rz(-0.41941386) q[3];
sx q[3];
rz(0.97233397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.45450777) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(1.3467849) q[2];
rz(0.421031) q[3];
sx q[3];
rz(-2.1249168) q[3];
sx q[3];
rz(-2.4436061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.43276697) q[0];
sx q[0];
rz(-2.3481752) q[0];
sx q[0];
rz(0.41473266) q[0];
rz(1.3955836) q[1];
sx q[1];
rz(-2.4826629) q[1];
sx q[1];
rz(-2.5674852) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78050437) q[0];
sx q[0];
rz(-1.6988806) q[0];
sx q[0];
rz(-0.9719406) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71530576) q[2];
sx q[2];
rz(-1.5537405) q[2];
sx q[2];
rz(0.41358435) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.50358665) q[1];
sx q[1];
rz(-1.9203016) q[1];
sx q[1];
rz(-1.867678) q[1];
x q[2];
rz(-1.8168713) q[3];
sx q[3];
rz(-0.68021357) q[3];
sx q[3];
rz(1.0567997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.85270143) q[2];
sx q[2];
rz(-0.39704278) q[2];
sx q[2];
rz(1.627702) q[2];
rz(-0.04143516) q[3];
sx q[3];
rz(-1.8824717) q[3];
sx q[3];
rz(-0.07853011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0742652) q[0];
sx q[0];
rz(-1.7794309) q[0];
sx q[0];
rz(2.4434027) q[0];
rz(2.9746338) q[1];
sx q[1];
rz(-1.0792462) q[1];
sx q[1];
rz(-2.4093157) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66757827) q[0];
sx q[0];
rz(-1.3587553) q[0];
sx q[0];
rz(1.7330806) q[0];
rz(-pi) q[1];
rz(-0.68533021) q[2];
sx q[2];
rz(-0.88980674) q[2];
sx q[2];
rz(0.53390098) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9975035) q[1];
sx q[1];
rz(-0.92458506) q[1];
sx q[1];
rz(-1.287582) q[1];
rz(-pi) q[2];
rz(-2.5712588) q[3];
sx q[3];
rz(-1.7551646) q[3];
sx q[3];
rz(1.0384699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3884864) q[2];
sx q[2];
rz(-2.8078418) q[2];
sx q[2];
rz(1.9801271) q[2];
rz(0.30900624) q[3];
sx q[3];
rz(-1.2523264) q[3];
sx q[3];
rz(1.9741612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56918615) q[0];
sx q[0];
rz(-2.321406) q[0];
sx q[0];
rz(0.15643315) q[0];
rz(-0.45174831) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(-2.3715473) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2609445) q[0];
sx q[0];
rz(-0.62793193) q[0];
sx q[0];
rz(-0.29675608) q[0];
x q[1];
rz(2.7355843) q[2];
sx q[2];
rz(-2.1210665) q[2];
sx q[2];
rz(-1.6183491) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8416482) q[1];
sx q[1];
rz(-2.3023459) q[1];
sx q[1];
rz(0.70920918) q[1];
rz(-1.2392063) q[3];
sx q[3];
rz(-1.0911897) q[3];
sx q[3];
rz(-0.94707205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.6440755) q[2];
sx q[2];
rz(-3.0023809) q[2];
sx q[2];
rz(2.1264123) q[2];
rz(1.7049568) q[3];
sx q[3];
rz(-2.5865343) q[3];
sx q[3];
rz(-2.5456583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62676936) q[0];
sx q[0];
rz(-2.5829782) q[0];
sx q[0];
rz(2.8998937) q[0];
rz(2.4027951) q[1];
sx q[1];
rz(-2.6530478) q[1];
sx q[1];
rz(-0.36639211) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.217427) q[0];
sx q[0];
rz(-2.1806742) q[0];
sx q[0];
rz(-0.92745552) q[0];
rz(-pi) q[1];
rz(1.1364469) q[2];
sx q[2];
rz(-0.98698101) q[2];
sx q[2];
rz(2.7455612) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.061325039) q[1];
sx q[1];
rz(-2.419201) q[1];
sx q[1];
rz(-0.25218833) q[1];
x q[2];
rz(1.4922769) q[3];
sx q[3];
rz(-2.588387) q[3];
sx q[3];
rz(-1.5403403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0176795) q[2];
sx q[2];
rz(-1.1992477) q[2];
sx q[2];
rz(-0.29385847) q[2];
rz(-0.014523225) q[3];
sx q[3];
rz(-1.9745275) q[3];
sx q[3];
rz(0.6706388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83207399) q[0];
sx q[0];
rz(-2.339395) q[0];
sx q[0];
rz(0.88395399) q[0];
rz(-2.4813095) q[1];
sx q[1];
rz(-2.1536004) q[1];
sx q[1];
rz(0.79137897) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7991379) q[0];
sx q[0];
rz(-0.90478169) q[0];
sx q[0];
rz(2.7333583) q[0];
x q[1];
rz(0.4653761) q[2];
sx q[2];
rz(-2.010979) q[2];
sx q[2];
rz(0.79681764) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2453354) q[1];
sx q[1];
rz(-0.63071139) q[1];
sx q[1];
rz(0.74951042) q[1];
x q[2];
rz(1.8459686) q[3];
sx q[3];
rz(-2.5857946) q[3];
sx q[3];
rz(-2.9142227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0749977) q[2];
sx q[2];
rz(-2.834088) q[2];
sx q[2];
rz(-1.0212612) q[2];
rz(-0.3785454) q[3];
sx q[3];
rz(-0.56423855) q[3];
sx q[3];
rz(2.5436201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1148949) q[0];
sx q[0];
rz(-2.7377991) q[0];
sx q[0];
rz(-0.60780203) q[0];
rz(-2.9027477) q[1];
sx q[1];
rz(-2.3919479) q[1];
sx q[1];
rz(1.3806608) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35683435) q[0];
sx q[0];
rz(-0.94576242) q[0];
sx q[0];
rz(-2.9805095) q[0];
rz(2.0062749) q[2];
sx q[2];
rz(-1.594211) q[2];
sx q[2];
rz(1.6053111) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7271991) q[1];
sx q[1];
rz(-0.52800035) q[1];
sx q[1];
rz(0.90331932) q[1];
x q[2];
rz(-2.4139666) q[3];
sx q[3];
rz(-2.8907223) q[3];
sx q[3];
rz(-0.48654702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7878788) q[2];
sx q[2];
rz(-1.6106662) q[2];
sx q[2];
rz(2.805368) q[2];
rz(1.1296889) q[3];
sx q[3];
rz(-1.7166398) q[3];
sx q[3];
rz(1.1415793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0040141) q[0];
sx q[0];
rz(-1.7332358) q[0];
sx q[0];
rz(-1.9144203) q[0];
rz(-0.54429383) q[1];
sx q[1];
rz(-1.2398564) q[1];
sx q[1];
rz(1.6287631) q[1];
rz(-1.2056247) q[2];
sx q[2];
rz(-2.361459) q[2];
sx q[2];
rz(-0.38103719) q[2];
rz(-0.23562283) q[3];
sx q[3];
rz(-2.1274673) q[3];
sx q[3];
rz(-2.6475788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
