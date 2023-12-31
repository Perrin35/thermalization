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
rz(1.8068846) q[0];
rz(-0.35285464) q[1];
sx q[1];
rz(-0.1605514) q[1];
sx q[1];
rz(0.97595739) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53544331) q[0];
sx q[0];
rz(-2.1436999) q[0];
sx q[0];
rz(-1.2826305) q[0];
rz(2.8085254) q[2];
sx q[2];
rz(-1.0330079) q[2];
sx q[2];
rz(1.2531467) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.18260278) q[1];
sx q[1];
rz(-0.62287736) q[1];
sx q[1];
rz(1.3401003) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7482412) q[3];
sx q[3];
rz(-2.297612) q[3];
sx q[3];
rz(-0.8918744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.41539899) q[2];
sx q[2];
rz(-1.4582863) q[2];
sx q[2];
rz(0.78380084) q[2];
rz(0.33256724) q[3];
sx q[3];
rz(-0.16210292) q[3];
sx q[3];
rz(1.3403085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48288229) q[0];
sx q[0];
rz(-1.4504526) q[0];
sx q[0];
rz(-0.17856199) q[0];
rz(-1.8042971) q[1];
sx q[1];
rz(-0.54090118) q[1];
sx q[1];
rz(3.1352502) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9316677) q[0];
sx q[0];
rz(-1.4405662) q[0];
sx q[0];
rz(-1.3757214) q[0];
x q[1];
rz(-0.52829929) q[2];
sx q[2];
rz(-2.106278) q[2];
sx q[2];
rz(2.2250125) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.22390631) q[1];
sx q[1];
rz(-2.8250541) q[1];
sx q[1];
rz(2.9427337) q[1];
rz(-pi) q[2];
rz(2.7553495) q[3];
sx q[3];
rz(-2.5442903) q[3];
sx q[3];
rz(-2.951705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1938842) q[2];
sx q[2];
rz(-2.7145553) q[2];
sx q[2];
rz(-2.2634899) q[2];
rz(-2.2795423) q[3];
sx q[3];
rz(-0.92249191) q[3];
sx q[3];
rz(-0.0058962065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55643117) q[0];
sx q[0];
rz(-1.0273902) q[0];
sx q[0];
rz(0.01097824) q[0];
rz(-2.7745461) q[1];
sx q[1];
rz(-1.9069907) q[1];
sx q[1];
rz(3.045851) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84155267) q[0];
sx q[0];
rz(-1.3944355) q[0];
sx q[0];
rz(-0.15533133) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50600608) q[2];
sx q[2];
rz(-1.9780469) q[2];
sx q[2];
rz(-0.72088748) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9853471) q[1];
sx q[1];
rz(-0.74838446) q[1];
sx q[1];
rz(-0.85703874) q[1];
rz(1.5355574) q[3];
sx q[3];
rz(-1.0198776) q[3];
sx q[3];
rz(-0.86722022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0720955) q[2];
sx q[2];
rz(-0.82192373) q[2];
sx q[2];
rz(-0.83479184) q[2];
rz(2.9299724) q[3];
sx q[3];
rz(-1.2303338) q[3];
sx q[3];
rz(0.37477469) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19462207) q[0];
sx q[0];
rz(-1.2397543) q[0];
sx q[0];
rz(0.34657493) q[0];
rz(0.52571458) q[1];
sx q[1];
rz(-2.321967) q[1];
sx q[1];
rz(-2.1077164) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82499408) q[0];
sx q[0];
rz(-0.98485095) q[0];
sx q[0];
rz(-0.98866776) q[0];
x q[1];
rz(2.1448574) q[2];
sx q[2];
rz(-2.0248932) q[2];
sx q[2];
rz(0.10276375) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8139207) q[1];
sx q[1];
rz(-1.76183) q[1];
sx q[1];
rz(-3.0753067) q[1];
rz(-0.86446188) q[3];
sx q[3];
rz(-2.7221788) q[3];
sx q[3];
rz(-2.1692587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6870849) q[2];
sx q[2];
rz(-1.7093753) q[2];
sx q[2];
rz(-1.7948077) q[2];
rz(-2.7205617) q[3];
sx q[3];
rz(-2.1249168) q[3];
sx q[3];
rz(-2.4436061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43276697) q[0];
sx q[0];
rz(-2.3481752) q[0];
sx q[0];
rz(-0.41473266) q[0];
rz(-1.3955836) q[1];
sx q[1];
rz(-0.65892977) q[1];
sx q[1];
rz(-2.5674852) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5363279) q[0];
sx q[0];
rz(-2.5308373) q[0];
sx q[0];
rz(-1.7954134) q[0];
rz(1.5933883) q[2];
sx q[2];
rz(-2.2859757) q[2];
sx q[2];
rz(-1.1720282) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9088604) q[1];
sx q[1];
rz(-2.6870011) q[1];
sx q[1];
rz(2.4652387) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2361034) q[3];
sx q[3];
rz(-1.4169766) q[3];
sx q[3];
rz(0.32115768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2888912) q[2];
sx q[2];
rz(-2.7445499) q[2];
sx q[2];
rz(-1.5138907) q[2];
rz(-3.1001575) q[3];
sx q[3];
rz(-1.2591209) q[3];
sx q[3];
rz(3.0630625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0742652) q[0];
sx q[0];
rz(-1.7794309) q[0];
sx q[0];
rz(0.69818991) q[0];
rz(-0.16695887) q[1];
sx q[1];
rz(-1.0792462) q[1];
sx q[1];
rz(-2.4093157) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4740144) q[0];
sx q[0];
rz(-1.3587553) q[0];
sx q[0];
rz(1.408512) q[0];
rz(-2.3789669) q[2];
sx q[2];
rz(-1.0566933) q[2];
sx q[2];
rz(-2.5800173) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9975035) q[1];
sx q[1];
rz(-0.92458506) q[1];
sx q[1];
rz(-1.8540107) q[1];
rz(-pi) q[2];
rz(-2.8090217) q[3];
sx q[3];
rz(-2.5453574) q[3];
sx q[3];
rz(2.3308844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3884864) q[2];
sx q[2];
rz(-0.33375084) q[2];
sx q[2];
rz(-1.1614655) q[2];
rz(-2.8325864) q[3];
sx q[3];
rz(-1.2523264) q[3];
sx q[3];
rz(-1.1674315) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5724065) q[0];
sx q[0];
rz(-0.82018667) q[0];
sx q[0];
rz(-2.9851595) q[0];
rz(2.6898443) q[1];
sx q[1];
rz(-2.2765171) q[1];
sx q[1];
rz(0.77004534) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2419287) q[0];
sx q[0];
rz(-0.97424346) q[0];
sx q[0];
rz(-1.3616189) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1595575) q[2];
sx q[2];
rz(-1.9141478) q[2];
sx q[2];
rz(0.26867586) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7913831) q[1];
sx q[1];
rz(-2.0767127) q[1];
sx q[1];
rz(-0.70178589) q[1];
rz(-0.55925925) q[3];
sx q[3];
rz(-2.5659605) q[3];
sx q[3];
rz(1.5881133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4975171) q[2];
sx q[2];
rz(-0.13921177) q[2];
sx q[2];
rz(1.0151803) q[2];
rz(-1.7049568) q[3];
sx q[3];
rz(-2.5865343) q[3];
sx q[3];
rz(2.5456583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5148233) q[0];
sx q[0];
rz(-2.5829782) q[0];
sx q[0];
rz(2.8998937) q[0];
rz(0.73879755) q[1];
sx q[1];
rz(-2.6530478) q[1];
sx q[1];
rz(0.36639211) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1358571) q[0];
sx q[0];
rz(-2.2861087) q[0];
sx q[0];
rz(-2.4321796) q[0];
rz(-pi) q[1];
rz(2.0051458) q[2];
sx q[2];
rz(-2.1546116) q[2];
sx q[2];
rz(2.7455612) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.39216343) q[1];
sx q[1];
rz(-0.87595075) q[1];
sx q[1];
rz(-1.3543345) q[1];
rz(1.4922769) q[3];
sx q[3];
rz(-0.55320569) q[3];
sx q[3];
rz(-1.6012524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0176795) q[2];
sx q[2];
rz(-1.1992477) q[2];
sx q[2];
rz(2.8477342) q[2];
rz(-3.1270694) q[3];
sx q[3];
rz(-1.1670651) q[3];
sx q[3];
rz(0.6706388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83207399) q[0];
sx q[0];
rz(-2.339395) q[0];
sx q[0];
rz(-2.2576387) q[0];
rz(-0.66028315) q[1];
sx q[1];
rz(-0.98799223) q[1];
sx q[1];
rz(-2.3502137) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6521097) q[0];
sx q[0];
rz(-1.8882505) q[0];
sx q[0];
rz(0.86274685) q[0];
x q[1];
rz(-2.0558526) q[2];
sx q[2];
rz(-1.9888478) q[2];
sx q[2];
rz(2.5784091) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1014175) q[1];
sx q[1];
rz(-1.1244332) q[1];
sx q[1];
rz(1.1091713) q[1];
rz(-pi) q[2];
x q[2];
rz(1.295624) q[3];
sx q[3];
rz(-2.5857946) q[3];
sx q[3];
rz(2.9142227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0749977) q[2];
sx q[2];
rz(-0.30750465) q[2];
sx q[2];
rz(-1.0212612) q[2];
rz(0.3785454) q[3];
sx q[3];
rz(-0.56423855) q[3];
sx q[3];
rz(-2.5436201) q[3];
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
rz(-3.1148949) q[0];
sx q[0];
rz(-2.7377991) q[0];
sx q[0];
rz(2.5337906) q[0];
rz(-2.9027477) q[1];
sx q[1];
rz(-2.3919479) q[1];
sx q[1];
rz(1.3806608) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5138807) q[0];
sx q[0];
rz(-2.4988334) q[0];
sx q[0];
rz(1.7895262) q[0];
x q[1];
rz(1.1353178) q[2];
sx q[2];
rz(-1.594211) q[2];
sx q[2];
rz(1.5362816) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.32523649) q[1];
sx q[1];
rz(-1.9776077) q[1];
sx q[1];
rz(-0.34646323) q[1];
x q[2];
rz(0.72762604) q[3];
sx q[3];
rz(-0.25087038) q[3];
sx q[3];
rz(0.48654702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7878788) q[2];
sx q[2];
rz(-1.5309265) q[2];
sx q[2];
rz(0.33622462) q[2];
rz(2.0119038) q[3];
sx q[3];
rz(-1.4249529) q[3];
sx q[3];
rz(1.1415793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
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
rz(-0.82472807) q[2];
sx q[2];
rz(-1.8246973) q[2];
sx q[2];
rz(1.4551103) q[2];
rz(2.9059698) q[3];
sx q[3];
rz(-2.1274673) q[3];
sx q[3];
rz(-2.6475788) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
